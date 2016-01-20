import time,logging,os

from environment import Environment
from system import System
from window import Window
from phase import Phase
from job import Job
from other import PMF,CollectiveVariable

class SiPMF():
  def __repr__(self):
    return "SiPMF({0},{1})".format(self.system,self.environment)

  def __init__(self,system,environment,log_filename="siPMF.log"):
    self.system=system
    self.environment=environment
    logging.basicConfig(filename=os.path.join(self.system.basedir,"siPMF.log"),level=logging.INFO,format='%(asctime)s: %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M %p')

  def Run(self,max_time,max_jobs,sleep_length):
    njobs=0
    n_running_jobs=0
    n_finished_jobs=0
    t0=time.time()
    continue_flag=True
    submit_flag=True
    if len(self.system.windows)==0:
      print "System does not contain any Window."
      print "Make sure to initialize the system before running it."
      return
    c=0
    logging.info("Starting the calculation with max_time={0}s,max_jobs={1},sleep_length={2}s".format(max_time,max_jobs,sleep_length))
    while continue_flag:
      c+=1
      save_flag=False
      continue_flag=False
      if njobs>=max_jobs or time.time()-t0>=max_time:submit_flag=False
      n_updated_windows=self.system.UpdateUnfinishedJobList(self.environment)
      n_finished_jobs=n_running_jobs-len(self.system.unfinished_jobs)
      if n_finished_jobs>0:logging.info("{0} jobs finished".format(n_finished_jobs))
      n_running_jobs=len(self.system.unfinished_jobs)
      #If some windows finised during last sleep (or if there are new windows)
      #We submit new jobs
      if n_updated_windows>0 and submit_flag:
        nj=self.system.SubmitNewJobs(self.environment)
        njobs+=nj
        n_running_jobs=len(self.system.unfinished_jobs)
        if nj>0:
          save_flag=True
          logging.info("Submitted {0} new jobs, {1} jobs are running or in the queue".format(nj,n_running_jobs))
      #If there are no running jobs, this means all current windows are finished
      #So we generate new windows
      if n_running_jobs==0:
        logging.info("No more jobs in the queue, checking whether to generate new windows")
        n_new_windows=self.system.GenerateNewWindows(self.environment)
        if n_new_windows!=0:
          save_flag=True
          logging.info("Generate {0} new windows. Total of {1} windows".format(n_new_windows,len(self.system.windows)))
        if n_new_windows!=0 and submit_flag:
          nj=self.system.SubmitNewJobs(self.environment)
          njobs+=nj
          n_running_jobs=len(self.system.unfinished_jobs)
          logging.info("Submitted {0} new jobs, {1} jobs are running or in the queue".format(nj,n_running_jobs))
        elif n_new_windows==0:
          logging.info("No new windows were generated")
      #Now we update the flags
      if n_running_jobs>0:continue_flag=True
      if njobs<max_jobs and time.time()-t0<max_time:submit_flag=True
      else:logging.info("Reached maximum number of jobs or time, no new jobs will be submitted.")
      if save_flag:
        self.system.Save("siPMF_state")
        logging.info("Saving the system.")
      if continue_flag:time.sleep(sleep_length)
    #Make sure the PMF is up to date before saving and stopping
    self.system.UpdatePMF(self.environment)
    self.system.Save("siPMF_state")
    logging.info("Stopping.")

