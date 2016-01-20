import os,subprocess,logging
import numpy as npy
import scipy.interpolate
import matplotlib.pyplot as plt
import itertools
from phase import Phase

class Window():
  def __repr__(self):
    return "Window({0},{1},{2},{3})".format(self.system,self.cv_values,self.spring_constants,self.parent)

  def __str__(self):
    if self.parent:
      return "{0} initialized form {1}".format(self.name,self.parent.name)
    else:return "{0}".format(self.name)

  def __init__(self,system,cv_values,spring_constants,parent=None):
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

  def SubmitNextPhase(self,environment):
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
    self.n_data=0
    for phase in self.phases:
      phase.UpdateDataCount()
      self.n_data+=phase.GetDataCount()
  
