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
    return "System({0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11})".format(self.basedir,self.cv_list,self.init_input_fname,self.run_input_fname,self.init_job_fname,self.run_job_fname,self.data_filename,self.init_nstep,self.run_nstep,self.n_data,self.max_E,self.temperature)

  def __init__(self,basedir,cv_list,init_input_fname,run_input_fname,init_job_fname,run_job_fname,data_filename,init_nstep,run_nstep,n_data,max_E,temperature,check_fnames):
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
    :param max_E: The free energy threshold to decide whether to extend the simulation to neighboring windows or not.
    :param temperature: The temperature at which WHAM is performed.

    :type basedir: :class:`str`
    :type cv_list: :class:`list`(:class:`CollectiveVariable`)
    :type init_input_fname: :class:`str`
    :type run_input_fname: :class:`str`
    :type init_job_fname: :class:`str`
    :type run_job_fname: :class:`str`
    :type data_filename: :class:`str`
    :type init_nstep: :class:`int`
    :type run_nstep: :class:`int`
    :type n_data: :class:`int`
    :type max_E: :class:`float`
    :type temperature: :class:`float`
    """
    self.basedir=basedir
    self.pmf_dir=os.path.join(basedir,"PMF")
    subprocess.call(["mkdir",self.pmf_dir])
    self.simu_dir=os.path.join(basedir,"windows")
    subprocess.call(["mkdir",self.simu_dir])
    self.temperature=temperature
    self.check_fnames=check_fnames
    self.max_E=max_E
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
    
    :type cv_values: :class:`list`(:class:`float`)
    :type spring_constants: :class:`list`(:class:`float`)
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
    :type environemnt: :class:`Environment`
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
    :type environemnt: :class:`Environment`
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

  def AddWindow(self,cv_values,spring_constants,parent_window=None):
    """
    Add a new window to the system

    :param cv_values: Values of the collective variables for this window
    :param spring_constants: Spring constants for the constraining potentials
    :param parent_window: Window from which the new window should be created

    :type cv_values: :class:`list`(:class:`float`)
    :type spring_constants: :class:`list`(:class:`float`)
    :type parent_window: :class:`Window`
    """
    self.windows.append(Window(self,cv_values,spring_constants,parent_window))
    self.updated_windows.append(self.windows[-1])

  def CalculatePMF(self,environment):
    """
    Calculate the PMF. The function generates the meta file listing all the data files and
    their corresponding CV values and spring constants and then calls WHAM.

    :param environment: The environment used to call WHAM
    :type environment: :class:`Environment`
    """
    pmf_cmd=[environment.wham_executable]
    if len(self.cv_list)==1:
      cv=cv_list[0]
      if not cv.periodicity:pmf_cmd.extend([cv.min_value,cv.max_value,cv.num_bins])
      else:pmf_cmd.extend(["P"+str(cv.periodicity),cv.min_value,cv.max_value,cv.num_bins])
    elif len(self.cv_list)==2:
      for x,cv in zip(["Px=","Py="],self.cv_list):
        if not cv.periodicity:p=x+str(0)
        else:p=x+str(cv.periodicity)
        pmf_cmd.extend([p,cv.min_value,cv.max_value,cv.num_bins])
    else:
      logging.error("Can only work with one or two collective variables")
      raise ValueError("Can only work with one or two collective variables")
    pmf_cmd.append(0.01)#tolerance
    pmf_cmd.append(self.temperature)
    pmf_cmd.append(0)#Number of pads for periodic variables
    f=open(self.path_to_pmf_input,"w")
    for window in self.windows:
      for phase in window.phases:
        if phase.type=="initialization":continue
        l=[phase.path_to_datafile]
        l.extend(phase.window.cv_values)
        l.extend(phase.window.spring_constants)
        f.write(" ".join([str(el) for el in l])+"\n")
    f.close()
    pmf_cmd.append(self.path_to_pmf_input)
    pmf_cmd.append(self.path_to_pmf_output)
    if len(self.cv_list)==2:pmf_cmd.append(1)#Use automatic masking
    pmf_cmd=[str(el) for el in pmf_cmd]
    print "".join(subprocess.list2cmdline(pmf_cmd))
    subprocess.call(pmf_cmd)

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
    self.pmf=PMF(pmf[:-1,:].transpose(),pmf[-1,:],self.cv_list)

  def UpdatePMF(self,environment):
    """
    Calculates the PMF (*CalculatePMF*), then reads the ouput PMF file (*ReadPMFFile*)
    and plots the new PMF. Finally, using the PMF, it assigns a free energy value to each window.

    :param environment: The environment used to call WHAM
    :type environment: :class:`Environment`
    """
    self.CalculatePMF(environment)
    self.ReadPMFFile()
    self.PlotPMF()
    #Windows get assigned the minimal free energy
    steps=[npy.arange(-cv.step_size/2.,cv.step_size/2.,cv.bin_size) for cv in self.cv_list]
    delta_cv_list=list(itertools.product(*steps))
    for window in self.windows:
      El=[float(self.pmf.interpolator.__call__(tuple(npy.array(window.cv_values)-npy.array(delta_cv)))) for delta_cv in delta_cv_list]
      window.free_energy=min(El)

  def PlotPMF(self):
    """
    Plot the PMF.
    """
    if self.pmf:
      filename="pmf_{0}".format(len(self.windows))
      self.pmf.Plot(self.pmf_dir,filename,max_E=1.5*self.max_E,windows=self.windows)
    else:
      logging.info("PMF has to be initialized before it can be plotted.")


  def GenerateNewWindows(self,environment):
    """
    Generates new windows around windows with free energy lower than *max_E* of the system.
    The function first updates the PMF (*UpdatePMF*) and then goes through all the windows
    in the system to decide where to generate new ones. If a new window can be generated from several
    windows, the one with lowest free energy will be used as parent.

    :param environment: The environment used to call WHAM
    :type environment: :class:`Environment`
    """
    self.UpdatePMF(environment)
    current_windows=[tuple(w.cv_values) for w in self.windows]
    steps=[[-cv.step_size,0,cv.step_size] for cv in self.cv_list]
    delta_cv_list=list(itertools.product(*steps))
    delta_cv_list.remove((0,0))
    new_windows={}
    for window in self.windows:
      if window.free_energy>self.max_E:continue
      for step in delta_cv_list:
        new_cv_vals=[el1+el2 for el1,el2 in zip(window.cv_values,step)]
        for i,(cv_val,cv) in enumerate(zip(new_cv_vals,self.cv_list)):
          if not cv.periodicity:continue
          if cv_val>cv.max_value:new_cv_vals[i]=cv_val-cv.periodicity
          elif cv_val<cv.min_value:new_cv_vals[i]=cv_val+cv.periodicity
        new_cv_vals=tuple(new_cv_vals)
        if new_cv_vals in current_windows:continue
        if any([cv_val>cv.max_value for cv_val,cv in zip(new_cv_vals,self.cv_list)]):continue
        if any([cv_val<cv.min_value for cv_val,cv in zip(new_cv_vals,self.cv_list)]):continue
        if not new_cv_vals in new_windows:new_windows.update({new_cv_vals:window})
        else:
          if window.free_energy<new_windows[new_cv_vals].free_energy:
            new_windows[new_cv_vals]=window
    for cv_values in new_windows:
      self.AddWindow(cv_values,new_windows[cv_values].spring_constants,new_windows[cv_values])
    return len(new_windows)
