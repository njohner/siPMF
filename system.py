"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains the :class:`System` object as well as a function to Load a :class:`System`
from a file.
"""
import os,subprocess,logging
import numpy as npy
import matplotlib.pyplot as plt
import itertools
import pickle
from window import Window
from other import *

__all__=('LoadSystem','System')

def LoadSystem(filename):
  """
  Loads a :class:`System` object from a file (using pickle)

  :param filename: The path to the file. A ".pkl" extension will be added automatically to the filename.
  :type filename:  :class:`str`
  """
  f=open(filename+".pkl","r")
  system=pickle.load(f)
  f.close()
  return system

class System():
  """
  The :class:`System` is the top class in the hierarchical description of a PMF calculation
  by umbrella sampling (:class:`System` contains :class:`Window` comprising several :class:`Phase`, which each
  have one :class:`Job`). The corresponding simulations will be carried out in a corresponding hierarchy of directories.
  The :class:`System` contains all the information about the current status
  of the simulations and is therefore the only object needed to restart a SiPMF calculation.
  """
  def __repr__(self):
    return "System({0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12})".format(self.basedir,self.cv_list,self.init_input_fname,self.run_input_fname,self.init_job_fname,self.run_job_fname,self.data_filename,self.init_nstep,self.run_nstep,self.n_data,self.max_E1,self.max_E2,self.temperature)

  def __init__(self,basedir,cv_list,init_input_fname,run_input_fname,init_job_fname,run_job_fname,data_filename,init_nstep,run_nstep,n_data,max_E1,max_E2,temperature,check_fnames=[],target_cv_vals=[],adapt_spring_constants=True,adapt_window_centers=True):
    """
    :param basedir: The root directory in which the PMF calculation will be performed. Windows and
     phases will correspond to subdirectories of *basedir*.
    :param cv_list: List of the collective variables (CV) used in the umbrella sampling
    :param init_input_fname: *basedir/init_input_fname* is the MD input file used for initializing new windows.
    :param run_input_fname: *basedir/run_input_fname* is the MD input file used to sample a window.
    :param init_job_fname: *basedir/init_job_fname* is the job file used to submit initialization runs to the cluster
    :param run_job_fname: *basedir/run_job_fname* is the job file used to submit sampling runs to the cluster
    :param data_filename: is the name of the file containing the values of the CVs. One such file should be generated
     for each sampling phase of every window by the corresponding job and be found in the corresponding directory.
    :param init_nstep: The number of steps for initialization phases
    :param run_nstep: The number of steps for sampling phases
    :param n_data: The number of data wanted for each window.
    :param max_E1: Lower boundary of free energy threshold to decide whether to extend the simulation to neighboring windows or not.
    :param max_E2: Upper boundary of free energy threshold to decide whether to extend the simulation to neighboring windows or not.
    :param temperature: The temperature at which WHAM is performed.
    :param check_fnames: Filenames that file be checked to exist to determine whether a phase has finished properly 
    :param target_cv_vals: Target values of the CVs. Once the system has reached these values it will only use max_E1 as energy threshold to generate new windows.
    :param adapt_spring_constants: If spring constants should be automatically adapted for each window.
    :param adapt_window_centers: If window centers should be automatically adapted for each window.

    :type basedir: :class:`str`
    :type cv_list: :class:`list` (:class:`~other.CollectiveVariable`)
    :type init_input_fname: :class:`str`
    :type run_input_fname: :class:`str`
    :type init_job_fname: :class:`str`
    :type run_job_fname: :class:`str`
    :type data_filename: :class:`str`
    :type init_nstep: :class:`int`
    :type run_nstep: :class:`int`
    :type n_data: :class:`int`
    :type max_E1: :class:`float`
    :type max_E2: :class:`float`
    :type temperature: :class:`float`
    :type check_fnames: :class:`list` (:class:`str`)
    :type target_cv_vals: :class:`list` (:class:`tuple` (:class:`float` ) )
    """
    self.basedir=basedir
    self.pmf_dir=os.path.join(basedir,"PMF")
    subprocess.call(["mkdir",self.pmf_dir])
    self.simu_dir=os.path.join(basedir,"windows")
    subprocess.call(["mkdir",self.simu_dir])
    self.hist_dir=os.path.join(basedir,"Histogram")
    subprocess.call(["mkdir",self.hist_dir])
    self.temperature=temperature
    self.check_fnames=check_fnames
    self.max_E1=max_E1
    self.max_E2=max_E2
    self.max_E_plot=1.2*max_E2
    self.cv_list=cv_list
    self.init_input_fname=init_input_fname
    self.run_input_fname=run_input_fname
    self.init_job_fname=init_job_fname
    self.run_job_fname=run_job_fname
    self.run_nstep=run_nstep
    self.init_nstep=init_nstep
    self.n_data=n_data
    self.dimensionality=len(cv_list)
    self.windows=[]
    self.unfinished_jobs=[]
    self.updated_windows=[]
    self.data_filename=data_filename
    self.path_to_pmf_input=os.path.join(self.pmf_dir,"pmf_input.txt")
    self.path_to_pmf_output=os.path.join(self.pmf_dir,"pmf.txt")
    self.pmf=None
    self.target_cv_vals=target_cv_vals
    self.reached_target=False
    self.adapt_spring_constants=adapt_spring_constants
    self.adapt_window_centers=adapt_window_centers
    
  def Save(self,filename):
    """
    Save the :class:`System` to a file using *pickle*.

    :param filename: The :class:`System` will be saved to *basedir/filename.pkl*
    :type filename: :class:`str`
    """
    f=open(os.path.join(self.basedir,filename+".pkl"),"w")
    pickle.dump(self,f)
    f.close()

  def Initialize(self,cv_values,spring_constants,init_restartdir):
    """
    Generate the first window.

    :param cv_values: The values to which the collective variables should be restrained for this window.
    :param spring_constants: The spring constants used to restrain the collective variables
    :param init_restartdir: *basedir/init_restartdir* should point to the directory containing the files 
     from which the simulation of the first phase of this window should be restarted. 
    
    :type cv_values: :class:`list` (:class:`float`)
    :type spring_constants: :class:`list` (:class:`float`)
    :type init_restartdir: :class:`str`
    """
    self.windows.append(Window(self,cv_values,spring_constants))
    self.init_restartdir=init_restartdir
    self.updated_windows.append(self.windows[-1])

  def UpdateUnfinishedJobList(self,environment):
    """
    Check whether running jobs are still in the queue and update 
    the list of unfinished jobs and updated windows.

    :param environment:  The environment used to check the job status
    :type environment: :class:`~environment.Environment`
    """
    n_crashed=0
    to_remove=[]
    for job in self.unfinished_jobs:
      job.UpdateStatus(environment)
      if job.queue_status=="finished":
        to_remove.append(job)
    for job in to_remove:
      self.unfinished_jobs.remove(job)
      self.updated_windows.append(job.phase.window)
      if not job.success:
        n_crashed+=1
        job.phase.window.last_phase_n_crashed+=1
        subprocess.call(["mv",job.phase.outdir,job.phase.outdir+"_back"+str(job.phase.window.last_phase_n_crashed)])
        job.phase.window.phases.remove(job.phase)
        if len(job.phase.window.phases)==0:job.phase.window.is_new=True
      else:
        job.phase.window.last_phase_n_crashed=0
      if job.phase.window.last_phase_n_crashed>=2:
        logging.error("{0} crashed twice, please verify why and correct error before restarting".format(job.phase))
        raise RuntimeError("{0} crashed twice, please verify why and correct error before restarting".format(job.phase))
    return len(self.updated_windows),n_crashed

  def SubmitNewJobs(self,environment):
    """
    Submit the next series of jobs. This does not create new windows, only go through the
    existing windows and submit the next job (initialization or run) for that window if necessary
    (if not enough data has been collected yet for that window).

    :param environment:  The environment used to submit the jobs
    :type environment: :class:`~environment.Environment`
    """
    for window in self.updated_windows:
      window.UpdateDataCount()
    to_remove=[]
    n_new_jobs=0
    for window in self.updated_windows:
      if window.n_data>=self.n_data:
        to_remove.append(window)
      else:
        window.SubmitNextPhase(environment)
        to_remove.append(window)
        n_new_jobs+=1
    for window in to_remove:
      self.updated_windows.remove(window)
    return n_new_jobs

  def GetPathToInitInputFile(self):
    """
    Get the path to the MD input file used to generate new windows (initialization phase)
    """
    return os.path.join(self.basedir,self.init_input_fname)

  def GetPathToRunInputFile(self):
    """
    Get the path to the MD input file used to sample windows (run phase)
    """
    return os.path.join(self.basedir,self.run_input_fname)

  def GetPathToInitJobFile(self):
    """
    Get the path to the job submission file for initialization phases
    """
    return os.path.join(self.basedir,self.init_job_fname)

  def GetPathToRunJobFile(self):
    """
    Get the path to the job submission file for run phases
    """
    return os.path.join(self.basedir,self.run_job_fname)

  def AddWindow(self,cv_values,spring_constants,shifts=None,parent_window=None):
    """
    Add a new window to the system

    :param cv_values: Values of the collective variables for this window
    :param spring_constants: Spring constants for the constraining potentials
    :param parent_window: Window from which the new window should be created

    :type cv_values: :class:`list` (:class:`float`)
    :type spring_constants: :class:`list` (:class:`float`)
    :type parent_window: :class:`~window.Window`
    """
    self.windows.append(Window(self,cv_values,spring_constants,shifts,parent_window))
    self.updated_windows.append(self.windows[-1])

  def FindWindow(self,cv_values,spring_constants=None):
    """
    Find window with given values of the cvs

    :param cv_values: Values of the collective variables for the window
    :type cv_values: :class:`list` (:class:`float`)
    """
    for w in self.windows:
      if tuple(w.cv_values)==tuple(cv_values):
        if not spring_constants:return w
        elif tuple(w.spring_constants)==tuple(spring_constants):return w
    else:
      return

  def FindWindows(self,cv_values,spring_constants=None):
    """
    Find windows with given values of the cvs

    :param cv_values: Values of the collective variables for the window
    :type cv_values: :class:`list` (:class:`float`)
    """
    w_list=[]
    for w in self.windows:
      if tuple(w.cv_values)==tuple(cv_values):
        if not spring_constants:w_list.append(w)
        elif tuple(w.spring_constants)==tuple(spring_constants):w_list.append(w)
    return w_list

  def UpdateDataFiles(self,n_skip=0,n_tot=-1,new_only=True):
    """
    Updates the datafiles of all windows. For each *window*, it takes the data 
    from all the run phases and writes it into a datafile, skipping the first n_skip data points
    and adding a maximum of n_tot data points for each window.
    *n_tot=-1* means there is no maximal number of data points.
    
    :param n_skip: The number of data points to skip.
    :param n_tot: The total number of data points used to calculate the PMF.
    :param new_only: Only add the data from phases that have not yet been added to the data files.
    :type n_skip: :class:`int`
    :type n_tot: :class:`int`
    :type new_only: :class:`bool`
    """
    for window in self.windows:
      window.UpdateDataFile(n_skip,n_tot,new_only)

  def CalculatePMF(self,environment):
    """
    Calculate the PMF. The function generates the meta file listing all the data files and
    their corresponding CV values and spring constants and then calls WHAM.

    :param environment: The environment used to call WHAM
    :type environment: :class:`~environment.Environment`
    """

    pmf_cmd=[environment.wham_executable]
    if len(self.cv_list)==1:
      cv=self.cv_list[0]
      if not cv.periodicity:pmf_cmd.extend([cv.wham_min_value,cv.wham_max_value,cv.wham_num_bins])
      else:pmf_cmd.extend(["P"+str(cv.periodicity),cv.wham_min_value,cv.wham_max_value,cv.wham_num_bins])
    elif len(self.cv_list)==2:
      for x,cv in zip(["Px=","Py="],self.cv_list):
        if not cv.periodicity:p=x+str(0)
        else:p=x+str(cv.periodicity)
        pmf_cmd.extend([p,cv.wham_min_value,cv.wham_max_value,cv.wham_num_bins])
    else:
      logging.error("Can only work with one or two collective variables")
      raise ValueError("Can only work with one or two collective variables")
    pmf_cmd.append(0.01)#tolerance
    pmf_cmd.append(self.temperature)
    num_pads=max([cv.num_pads for cv in self.cv_list])
    pmf_cmd.append(num_pads)#Number of pads for periodic variables
    f=open(self.path_to_pmf_input,"w")
    for window in self.windows:
      if not os.path.isfile(window.path_to_datafile):continue
      l=[window.path_to_datafile]
      l.extend([cvv+cvs for cvv,cvs in zip(window.cv_values,window.cv_shifts)])
      l.extend(window.spring_constants)
      f.write(" ".join([str(el) for el in l])+"\n")
    f.close()
    pmf_cmd.append(self.path_to_pmf_input)
    pmf_cmd.append(self.path_to_pmf_output)
    if len(self.cv_list)==2:pmf_cmd.append(1)#Use automatic masking
    pmf_cmd=[str(el) for el in pmf_cmd]
    print "".join(subprocess.list2cmdline(pmf_cmd))
    t0=time.time()
    subprocess.call(pmf_cmd)
    logging.info("WHAM finished in {0}s".format(time.time()-t0))

  def ReadPMFFile(self):
    """
    Read the PMF file generated by the *CalculatePMF* function.
    """
    f=open(self.path_to_pmf_output,"r")
    nd=self.dimensionality+1
    pmf=[[] for i in range(nd)]
    for l in f:
      if l.startswith("#"):continue
      if l=="\n":continue
      s=l.split()
      if len(s)<nd:
        print l
        continue
      for i in range(nd):pmf[i].append(float(s[i]))
    f.close()
    pmf=npy.array(pmf)
    if self.dimensionality==1:
      self.pmf=PMF(pmf[0,:],pmf[1,:],self.cv_list,1.25*self.max_E_plot)
    else:
      self.pmf=PMF(pmf[:-1,:].transpose(),pmf[-1,:],self.cv_list,1.25*self.max_E_plot)

  def UpdatePMF(self,environment,n_skip=0,n_tot=-1,new_only=True,fname_extension=""):
    """
    Calculates the PMF (*CalculatePMF*), then reads the ouput PMF file (*ReadPMFFile*)
    and plots the new PMF. Finally, using the PMF, it assigns a free energy value to each window.

    :param environment: The environment used to call WHAM
    :type environment: :class:`~environment.Environment`
    """
    logging.info("Updating PMF")
    self.UpdateDataFiles(n_skip,n_tot,new_only)
    self.CalculatePMF(environment)
    self.ReadPMFFile()
    self.PlotPMF(fname_extension)
    #Windows get assigned the minimal free energy
    steps=[npy.arange(-cv.step_size/2.,cv.step_size/2.,cv.bin_size) for cv in self.cv_list]
    delta_cv_list=list(itertools.product(*steps))
    for window in self.windows:
      El=[float(self.pmf.interpolator.__call__(tuple(npy.array(window.cv_values)-npy.array(delta_cv)))) for delta_cv in delta_cv_list]
      window.free_energy=min([el for el in El if not npy.isnan(el)])
      window.curvatures=self.pmf.GetCurvatures(npy.array(window.cv_values),npy.array([cv.step_size/2. for cv in self.cv_list]))
  
  def ShiftWindowFreeEnergies(self,min_val=0):
    """
    Shift the free energies of the windows such that the window with the lowest free
    energy has a free energy of min_val.
    
    :param min_val: the free energy assigned to the window with lowest free energy
    :type min_val: :class:`float`
    """
    fes=[w.free_energy for w in self.windows]
    shift=min_val-min(fes)
    for w in self.windows:
      w.free_energy+=shift
    return shift

  def PlotPMF(self,fname_extension=""):
    """
    Plot the PMF.
    """
    if self.pmf:
      filename="pmf_{0}{1}".format(len(self.windows),fname_extension)
      self.pmf.Plot(self.pmf_dir,filename,max_E=self.max_E_plot,windows=self.windows)
    else:
      logging.info("PMF has to be initialized before it can be plotted.")


  def GenerateNewWindows(self,environment):
    """
    Generates new windows expected to have low free energy from windows with low free energy themselves.
    The function first updates the PMF (*UpdatePMF*) and then goes through all the windows
    in the system to find potential new windows to generate. Specifically it will start with a maximal
    free energy threshold of *max_E1*, go through all the windows and for those with free energy below
    the current threshold, it will look at all potential neighbor windows. It estimates the free
    of each potential new window and if it is below the free energy threshold, it adds that window to the system.
    If no new window was generated, it increases the threshold in 5 steps up to *max_E2*, each time repeating
    the procedure described above. If a new window can be generated from several
    windows, the one with lowest free energy will be used as parent.

    :param environment: The environment used to call WHAM
    :type environment: :class:`~environment.Environment`
    """
    self.UpdatePMF(environment)
    self.PlotHistogram()
    fe_shift=self.ShiftWindowFreeEnergies()
    current_windows=[tuple(w.cv_values) for w in self.windows]
    steps=[[-cv.step_size,0,cv.step_size] for cv in self.cv_list]
    delta_cv_list=list(itertools.product(*steps))
    delta_cv_list.remove(tuple([0 for cv in self.cv_list]))
    if not self.reached_target:fe_step=(self.max_E2-self.max_E1)/5.
    else:fe_step=1.0
    for max_free_energy in npy.arange(self.max_E1,self.max_E2+fe_step/2.,fe_step):
      new_windows={}
      for window in self.windows:
        if window.free_energy>max_free_energy:continue
        for step in delta_cv_list:
          new_cv_vals=[el1+el2 for el1,el2 in zip(window.cv_values,step)]
          for i,(cv_val,cv) in enumerate(zip(new_cv_vals,self.cv_list)):
            if not cv.periodicity:continue
            if cv_val>cv.max_value:new_cv_vals[i]=cv_val-cv.periodicity
            elif cv_val<=cv.min_value:new_cv_vals[i]=cv_val+cv.periodicity
          new_cv_vals=tuple(new_cv_vals)
          if new_cv_vals in current_windows:continue
          if any([cv_val>cv.max_value for cv_val,cv in zip(new_cv_vals,self.cv_list)]):continue
          if any([cv_val<cv.min_value for cv_val,cv in zip(new_cv_vals,self.cv_list)]):continue
          if not new_cv_vals in new_windows:
            new_windows[new_cv_vals]={}
            new_windows[new_cv_vals]["parent"]=window
          else:
            if window.free_energy<new_windows[new_cv_vals]["parent"].free_energy:
              new_windows[new_cv_vals]["parent"]=window
      #Estimate what the free energy in that window could be
      steps2=[npy.arange(-cv.step_size/2.,cv.step_size/2.,cv.bin_size) for cv in self.cv_list]
      delta_cv_list=list(itertools.product(*steps))
      for cv_values in new_windows:
        El=[float(self.pmf.interpolator.__call__(tuple(npy.array(cv_values)-npy.array(delta_cv)))) for delta_cv in delta_cv_list]
        try:new_windows[cv_values]["free_energy"]=min([el for el in El if not npy.isnan(el)])
        except:new_windows[cv_values]["free_energy"]=npy.nan
      #Add the new windows
      n_new_windows=0
      for cv_values in new_windows:
        if new_windows[cv_values]["free_energy"]+fe_shift<max_free_energy:
          parent=new_windows[cv_values]["parent"]
          if self.adapt_spring_constants:
            curvatures=parent.curvatures
            spring_constants=[]
            for i in range(self.dimensionality):
              spring_constant=min(max(self.cv_list[i].min_spring_constant,abs(curvatures[i])),self.cv_list[i].max_spring_constant) 
              spring_constants.append(int(100*spring_constant/self.cv_list[i].min_spring_constant)*self.cv_list[i].min_spring_constant/100)
            logging.info("cv values for new window {0}, curvatures {1}, spring constants {2}".format(cv_values,curvatures,spring_constants))
          else:
            spring_constants=[self.cv_list[i].min_spring_constant for i in range(self.dimensionality)]
          if self.adapt_window_centers:
            p=parent.cv_values
            shifts=[]
            for i in range(self.dimensionality):
              step_size=self.cv_list[i].step_size
              step=npy.zeros(self.dimensionality)
              step[i]=cv_values[i]-p[i]
              if self.cv_list[i].periodicity:
                if step[i]>self.cv_list[i].periodicity/2.0:
                  step[i]-=self.cv_list[i].periodicity
                elif step[i]<-self.cv_list[i].periodicity/2.0:
                  step[i]+=self.cv_list[i].periodicity
              F=2.0*(self.pmf.GetValue(p+step/2.0)-self.pmf.GetValue(p))/step_size
              shift=npy.sign(step[i])*F/(2.0*spring_constants[i])
              shift=npy.sign(shift)*min(self.cv_list[i].max_shift,abs(shift))
              shifts.append(shift)
            logging.info("cv values for new window {0}, spring constants {1}, shifts {2}".format(cv_values,spring_constants,shifts))
          else:
            shifts=None
          self.AddWindow(cv_values,spring_constants,shifts,parent)
          n_new_windows+=1
          if cv_values in self.target_cv_vals:
            self.reached_target=True
            self.max_E2=self.max_E1
            logging.info("Reached target CV value: cv={0}. Setting max_E2=max_E1".format(cv_values))
      if n_new_windows>=1:return n_new_windows,max_free_energy
    return n_new_windows,max_free_energy

  def CalculateWindowsHistConvergence(self,environment,n_skip_list,n_tot_list,update_data_files=True,pool_windows=False):
    logging.info("Calculating histogram convergence for each window")
    if update_data_files:self.UpdateDataFiles(0,-1,False)#Make sure all the data is in the data files.
    windows_convergence_list=[]
    if not pool_windows:windows_list=[[w] for w in self.windows]
    else:
      cv_values=[]
      windows_list=[]
      for w in self.windows:
        if tuple(w.cv_values) in cv_values:continue
        else:
          cv_values.append(tuple(w.cv_values))
          windows_list.append(self.FindWindows(w.cv_values))
    print "working with {0} window pools containing {1} windows each".format(len(windows_list),npy.average([len(el) for el in windows_list]))
    for wl in windows_list:
      data_list=[]
      for window in wl:
        data=[[] for i in range(self.dimensionality)]
        if not os.path.isfile(window.path_to_datafile):continue
        f=open(window.path_to_datafile,"r")
        for l in f:
          s=l.split()
          for i in range(self.dimensionality):data[i].append(float(s[i+1]))
        f.close()
        data_list.append(data)
      hist_range=[(cv.min_value,cv.max_value) for cv in self.cv_list]
      bins=[self.cv_list[0].num_bins,self.cv_list[1].num_bins]
      hist_list=[]
      for n_skip,n_tot in zip(n_skip_list,n_tot_list):
        d=[[] for i in range(self.dimensionality)]
        for data in data_list:
          for i,el in enumerate(data):d[i].extend(el[n_skip:n_skip+n_tot])
        #d=[el[n_skip:n_skip+n_tot] for el in data]
        nd=float(len(d[0]))
        hist_list.append(npy.histogramdd(d,range=hist_range,bins=bins)[0]/nd)
      ref_hist=hist_list[-1]
      windows_convergence_list.append([0.5*npy.sum(npy.abs(h-ref_hist)) for h in hist_list[:-1]])
    convergence=[npy.average([c[i] for c in windows_convergence_list]) for i in range(len(windows_convergence_list[0]))]
    convergence_std=[npy.std([c[i] for c in windows_convergence_list]) for i in range(len(windows_convergence_list[0]))]
    return windows_convergence_list,convergence,convergence_std

  def CalculatePMFConvergence(self,environment,n_skip_list,n_tot_list,max_E):
    """
    This function generates a set of PMFs from subsets of the whole data,
    defined by *n_skip_list* and *n_tot_list* and calculates their convergence.
    The last PMF of the set is used as reference. Only points of the PMF with
    an energy below *max_E* are used in the calculation.

    :param n_skip_list: List of the number of data points to skip.
    :param n_tot_list: List of the total number of data points used to calculate the PMF.
    :param max_E: maximal energy of considered points.
    :type n_skip_list: :class:`list` (:class:`int`)
    :type n_tot_list: :class:`list` (:class:`int`)
    :type max_E: :class:`float`
    """
    #Calculate the PMFs
    pmf_list=[]
    for n_skip,n_data in zip(n_skip_list,n_tot_list):
      self.UpdatePMF(environment,n_skip,n_data,new_only=False,fname_extension="_skip{0}_tot{1}".format(n_skip,n_data))
      pmf_list.append(npy.array(self.pmf.values))
    #Find common mask
    ref_pmf=pmf_list[-1]
    m=ref_pmf<max_E
    print npy.sum(m)
    for pmf in pmf_list[:-1]:
      m=m*(pmf<max_E)
      print npy.sum(m)
    #Now we normalize all pmfs
    for pmf in pmf_list:
      pmf-=npy.average(pmf[m])
    #Now we calculate the convergence
    ref_pmf=pmf_list[-1][m]
    res_list=[]
    for pmf in pmf_list[:-1]:
      d=pmf[m]-ref_pmf
      res_list.append(npy.sqrt(npy.sum(d*d)/d.size))
    return res_list

  def PlotHistogram(self,fname_extension=""):
    """
    Plots the histogram of the accumulated data.
    """
    data=[[] for i in range(self.dimensionality)]
    for window in self.windows:
      if not os.path.isfile(window.path_to_datafile):continue
      f=open(window.path_to_datafile,"r")
      for l in f:
        s=l.split()
        for i in range(self.dimensionality):data[i].append(float(s[i+1]))
      f.close()
    filename="histogram_{0}{1}".format(len(self.windows),fname_extension)
    hist_range=[(cv.min_value,cv.max_value) for cv in self.cv_list]
    plt.figure()
    if self.dimensionality==2:
      bins=[self.cv_list[0].num_bins,self.cv_list[1].num_bins]
      hist_range=[(self.cv_list[0].min_value,self.cv_list[0].max_value),(self.cv_list[1].min_value,self.cv_list[1].max_value)]
      plt.hist2d(data[0],data[1],range=hist_range,bins=bins)
      if self.cv_list[0].units:plt.xlabel("{0} [{1}]".format(self.cv_list[0].name,self.cv_list[0].units))
      else:plt.xlabel("{0}".format(self.cv_list[0].name))
      if self.cv_list[1].units:plt.ylabel("{0} [{1}]".format(self.cv_list[1].name,self.cv_list[1].units))
      else:plt.ylabel("{0}".format(self.cv_list[1].name))
      plt.colorbar()
    elif self.dimensionality==1:
      bins=self.cv_list[0].num_bins
      plt.hist(data[0],bins=bins,range=hist_range[0])
      if self.cv_list[0].units:plt.xlabel("{0} [{1}]".format(self.cv_list[0].name,self.cv_list[0].units))
      else:plt.xlabel("{0}".format(self.cv_list[0].name))
      plt.ylabel("Count")
    plt.savefig(os.path.join(self.hist_dir,filename))
    plt.close()

  def SetNewBaseDir(self,basedir):
    self.basedir=basedir
    self.pmf_dir=os.path.join(basedir,"PMF")
    self.simu_dir=os.path.join(basedir,"windows")
    self.hist_dir=os.path.join(basedir,"Histogram")
    self.path_to_pmf_input=os.path.join(self.pmf_dir,"pmf_input.txt")
    self.path_to_pmf_output=os.path.join(self.pmf_dir,"pmf.txt")
    for window in self.windows:
      window.subdir=os.path.join(self.simu_dir,window.name)
      window.path_to_datafile=os.path.join(window.subdir,"data.txt")
      for phase in window.phases:
        phase.outdir=os.path.join(window.subdir,phase.name)
        phase.path_to_datafile=os.path.join(phase.outdir,self.data_filename)
        #We do not update jobs as this is not useful
    return


