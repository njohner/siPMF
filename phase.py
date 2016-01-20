import os,subprocess
from job import Job

class Phase():
  def __repr__(self):
    return "Phase({0},{1},{2},{3})".format(self.window,self.name,self.type,self.parent_phase)
  
  def __str__(self):
    cvs=" , ".join([cv.name+"="+str(cv_val) for cv,cv_val in zip(self.window.system.cv_list,self.window.cv_values)])
    if self.parent_phase:
      cvs2=" , ".join([cv.name+"="+str(cv_val) for cv,cv_val in zip(self.parent_phase.window.system.cv_list,self.parent_phase.window.cv_values)])
      return "Phase {0}: This is a {1} phase with {2} and initialized form {3}".format(self.name,self.type,cvs,cvs2)
    else:return "Phase {0}: This is a {1} phase with {2}".format(self.name,self.type,cvs)
  
  def __init__(self,window,phase_name,phase_type,parent_phase=None):
    self.window=window
    self.name=phase_name
    self.type=phase_type
    self.n_data=0
    self.parent_phase=parent_phase
    if self.parent_phase:
      self.restartdir=parent_phase.outdir
      #self.restartname=parent_phase.outname
    else:
      self.restartdir=self.window.system.init_restartdir
      #self.restartname=self.window.system.init_restartname
    self.outdir=os.path.join(self.window.subdir,self.name)
    self.path_to_datafile=os.path.join(self.outdir,self.window.system.data_filename)
    #self.outname=self.name

  def Initialize(self):
    print "New",self
    r=subprocess.call(["mkdir",self.outdir])
    if r!=0:raise IOError("Problem creating an output directory.")
    self.job=Job(self)

  def UpdateDataCount(self):
    if self.type=="initialization":self.n_data=0
    else:
      if not os.path.isfile(self.path_to_datafile):self.n_data=0
      else:
        f=open(self.path_to_datafile,"r")
        self.n_data=len([1 for line in f if not line.startswith("#")])
        f.close()

  def GetDataCount(self):
    return self.n_data
