import time

from environment import Environment
from system import System
from window import Window
from phase import Phase
from job import Job
from other import PMF,CollectiveVariable

class SiPMF():
  def __repr__(self):
    return "SiPMF({0},{1})".format(self.system,self.environment)

  def __init__(self,system,environment):
    self.system=system
    self.environment=environment

  def Run(self,max_time,max_jobs,sleep_length):
    njobs=0
    t0=time.time()
    continue_flag=True
    submit_flag=True
    if len(self.system.windows)==0:
      print "System does not contain any Window."
      print "Make sure to initialize the system before running it."
      return
    c=0
    while continue_flag:
      c+=1
      continue_flag=False
      if njobs>=max_jobs or time.time()-t0>=max_time:submit_flag=False
      n_updated_windows=self.system.UpdateUnfinishedJobList(self.environment)
      n_running_jobs=len(self.system.unfinished_jobs)
      #If some windows finised during last sleep (or if there are new windows)
      #We submit new jobs
      print c,continue_flag,submit_flag,n_updated_windows,n_running_jobs
      if n_updated_windows>0 and submit_flag:
        nj=self.system.SubmitNewJobs(self.environment)
        njobs+=nj
        n_running_jobs=len(self.system.unfinished_jobs)
        self.system.Save("siPMF_state")
      #If there are no running jobs, this means all current windows are finished
      #So we generate new windows
      if n_running_jobs==0:
        n_new_windows=self.system.GenerateNewWindows(self.environment)
        if n_new_windows!=0 and submit_flag:
          nj=self.system.SubmitNewJobs(self.environment)
          njobs+=nj
          n_running_jobs=len(self.system.unfinished_jobs)
          self.system.Save("siPMF_state")
        elif n_new_windows==0:
          print "No new window was generated."
      #Now we update the flags
      if n_running_jobs>0:continue_flag=True
      if njobs<max_jobs and time.time()-t0<max_time:submit_flag=True
      if continue_flag:time.sleep(sleep_length)
    #Make sure the PMF is up to date before saving and stopping
    self.system.UpdatePMF(self.environment)
    self.system.Save("siPMF_state")


