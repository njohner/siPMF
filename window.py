"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This file contains the :class:`Window` object which represents a simulation window.
"""
import os,subprocess,logging
import numpy as npy
import scipy.interpolate
import matplotlib.pyplot as plt
import itertools
from phase import Phase

class Window():
  """
  This class is second in the hierarchical structure used in SiPMF. The :class:`System` contains 
  several (usually many) class:`Window` which correspond to the :class:`System` restrained to a small
  region of the CV space using a quadratic restraining potential. The restraining potential has to
  be of the form **1/2 spring_constant(x-cv_value)**. Notice the 1/2 in the definition used here, whereas some
  MD codes do not have this factor. **Carefully check that the potential defined here matches
  the potential defined in the MD code you are using!**
  """
  def __repr__(self):
    return "Window({0},{1},{2},{3})".format(self.system,self.cv_values,self.spring_constants,self.parent)

  def __str__(self):
    if self.parent:
      return "{0} initialized form {1}".format(self.name,self.parent.name)
    else:return "{0}".format(self.name)

  def __init__(self,system,cv_values,spring_constants,parent=None):
    """
    :param system: The system to which the window belongs
    :param cv_values: The centers of the quadratic potentials restraining the CVs.
    :param spring_constants: The spring constants of the potentials restraining the CVs
    :param parent: The parent window
    :type system: :class:`~system.System`
    :type cv_values: :class:`list` (:class:`float`)
    :type spring_constants: :class:`list` (:class:`float`)
    :type parent: :class:`~window.Window`
    """
    self.cv_names=[cv.name for cv in system.cv_list]
    self.cv_values=cv_values
    self.spring_constants=spring_constants
    self.is_new=True
    self.n_data=0
    self.n_run_phases=0
    self.phases=[]
    self.system=system
    self.name="_".join(["".join([cvn,str(cvv)]) for cvn,cvv in zip(self.cv_names,self.cv_values)])
    self.subdir=os.path.join(system.simu_dir,self.name)
    self.parent=parent
    logging.info("New window: {0}".format(self))
    if os.path.isdir(self.subdir):
      logging.error("Directory already exists, program stops to avoid overwriting {0}.".format(self.subdir))
      raise IOError("Directory already exists, program stops to avoid overwriting {0}.".format(self.subdir))
    r=subprocess.call(["mkdir",self.subdir])
    if r!=0:
      logging.error("Problem creating output directory {0}.".format(self.subdir))
      raise IOError("Problem creating output directory {0}.".format(self.subdir))
    self.last_phase_n_crashed=0
    self.path_to_datafile=os.path.join(self.subdir,"data.txt")
    self.datafile_n_data_tot=0
    self.datafile_n_data_skipped=0

  def SubmitNextPhase(self,environment):
    """
    Automatically creates the appropriate next :class:`Phase` and corresponding :class:`Job` and
    submits it to the cluster. What the appropriate next phase is, is determined as follows:

    - If the window does not have a parent and does not contain any phase yet, the new phase will be 
    a run phase using as restart the *System.init_restartdir*
    - If the window has a parent phase does not contain any phase yet, the new phase will be 
    an initialization phase using as restart the last phase of the parent phase.
    - If the window already contains one or several phases, the new phase will be 
    an run phase using as restart the last phase of this window.

    :param environment: The environment used to submit the job to the cluster.
    :type environment: :class:`~environment.Environment`
    """
    if self.is_new:
      if self.parent:
        phase_name="initialization"
        phase_type="initialization"
        restart_phase=self.parent.phases[-1]
        self.phases.append(Phase(self,phase_name,phase_type,restart_phase))
      else:
        phase_name="phase1"
        phase_type="run"
        self.phases.append(Phase(self,phase_name,phase_type))
      self.is_new=False
    else:  
      phase_name="phase"+str(self.n_run_phases+1)
      phase_type="run"
      self.phases.append(Phase(self,phase_name,phase_type,self.phases[-1]))
    next_phase=self.phases[-1]
    next_phase.Initialize()
    next_phase.job.Submit(environment)
          
  def UpdateDataCount(self):
    """
    Update the total number of data accumulated for this window (sum over the data in each phase of the window).
    """
    self.n_data=0
    for phase in self.phases:
      phase.UpdateDataCount()
      self.n_data+=phase.GetDataCount()

  def UpdateDataFile(self,n_skip=0,n_tot=-1,new_only=True):
    """
    Updates the window's datafile. It takes the data from all the run phases and
    writes it into its own datafile, skipping the first n_skip data points
    and adding a maximum of n_tot data points.
    *n_tot=-1* means there is no maximal number of data points.

    :param n_skip: The number of data points to skip.
    :param n_tot: The total number of data points used to calculate the PMF.
    :param new_only: Only add the data from phases that have not yet been added to the data files.
    :type n_skip: :class:`int`
    :type n_tot: :class:`int`
    :type new_only: :class:`bool`
    """
    if not new_only:
      self.datafile_n_data_tot=0
      self.datafile_n_data_skipped=0
      if os.path.isfile(self.path_to_datafile):
        subprocess.call(["rm",self.path_to_datafile])
    for phase in self.phases:
      phase.AddDataToWindow(n_skip,n_tot,new_only)

  def ReadDataFile(self):
    """
    Reads the datafile of the window and returns a tuple with
    the list of times as the first element. The second element in
    the tuple is a list containing one list of values for each CV.
    """
    f=open(self.path_to_datafile,"r")
    ll=f.readlines()
    f.close()
    t=[]
    cvs=[[] for i in range(self.system.dimensionality)]
    for l in ll:
      s=l.split()
      t.append(float(s[0]))
      for i in range(self.system.dimensionality):cvs[i].append(float(s[i+1]))
    return (t,cvs)

