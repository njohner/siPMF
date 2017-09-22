"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains the :class:`siPMF` class
"""
import time
import logging


class SiPMF():
  """
  The :class:`siPMF` class defines the process that will supervise the whole calculation of
  the free energy landscape. It will survey the status of the jobs, submit new jobs, create
  new windows and so on.
  """

  def __repr__(self):
    return "SiPMF({0},{1})".format(self.system, self.environment)

  def __init__(self, system, environment):
    """
    :param system: The system that will be studied
    :param environment: The environment
    :type system: :class:`~system.System`
    :type environment: :class:`~environment.Environment`
    """
    self.system = system
    self.environment = environment
    # We add all windows to the list of updated windows, to make sure
    # they will get checked to see whether a new phase should be generated
    # or not for each window. This allows to add sampling by changing system.n_data
    to_skip = []
    for job in self.system.unfinished_jobs:
      to_skip.append(job.phase.window)
    for w in self.system.windows:
      if w in to_skip:
        continue
      if w in self.system.updated_windows:
        continue
      self.system.updated_windows.append(w)

  def Run(self, max_time, max_jobs, sleep_length, generate_new_windows=True):
    """
    Run the process to explore the free energy landscape. The process is an infinite loop in which
    it will sleep for some time, then when it wakes up it checks the status of the jobs in the queue.
    If some jobs have finished, it will submit new phases for the corresponding windows if necessary
    (i.e not enough data accumulated yet). When all the windows have been properly sampled, and no more
    jobs are running, it will create new windows and start sampling those.
    The process either stops when it cannot create any new windows or when it has run for more than a preset
    time or if it has submitted more jobs than a preset maximal number of jobs. For any of these stopping
    conditions, the process will keep on running until there are no jobs left in the queue, then it will
    update the PMF once more, save its state and stop.
    When running it saves its state everytime it has changed. It also recalculates and plots the PMF
    every time it generates new windows.

    :param max_time: Maximal time (in seconds) the process will run for
    :param max_jobs: Maximal number of jobs the process will submit
    :param sleep_length: time (in seconds) the process will sleep between two cycles.
    :type max_time: :class:`int`
    :type max_jobs: :class:`int`
    :type sleep_length: :class:`int`
    """
    njobs = 0
    n_running_jobs = 0
    n_finished_jobs = 0
    t0 = time.time()
    continue_flag = True
    submit_flag = True
    if len(self.system.windows) == 0:
      print "System does not contain any Window."
      print "Make sure to initialize the system before running it."
      return
    c = 0
    logging.info("Starting the calculation with max_time={0}s,max_jobs={1},sleep_length={2}s".format(
        max_time, max_jobs, sleep_length))
    while continue_flag:
      c += 1
      save_flag = False
      continue_flag = False
      if njobs >= max_jobs or time.time() - t0 >= max_time:
        submit_flag = False
      n_updated_windows, n_crashed_jobs = self.system.UpdateUnfinishedJobList(
          self.environment)
      n_finished_jobs = n_running_jobs - len(self.system.unfinished_jobs)
      if n_finished_jobs > 0:
        logging.info("{0} jobs finished among which {1} crashed".format(
            n_finished_jobs, n_crashed_jobs))
      n_running_jobs = len(self.system.unfinished_jobs)
      # If some windows finised during last sleep (or if there are new windows)
      # We submit new jobs
      if n_updated_windows > 0 and submit_flag:
        nj = self.system.SubmitNewJobs(self.environment)
        njobs += nj
        n_running_jobs = len(self.system.unfinished_jobs)
        if nj > 0:
          save_flag = True
          logging.info("Submitted {0} new jobs, {1} jobs are running or in the queue".format(
              nj, n_running_jobs))
      # If there are no running jobs, this means all current windows are finished
      # So we generate new windows
      if n_running_jobs == 0:
        if generate_new_windows:
          logging.info(
              "No more jobs in the queue, checking whether to generate new windows")
          n_new_windows, fe_threshold = self.system.GenerateNewWindows(
              self.environment)
        else:
          logging.info(
              "No more jobs in the queue and no new windows will be generated (generate_new_windows=False)")
          n_new_windows = 0
        if n_new_windows != 0:
          save_flag = True
          logging.info("Generate {0} new windows with free energy threshold={1}. Total of {2} windows".format(
              n_new_windows, fe_threshold, len(self.system.windows)))
        if n_new_windows != 0 and submit_flag:
          nj = self.system.SubmitNewJobs(self.environment)
          njobs += nj
          n_running_jobs = len(self.system.unfinished_jobs)
          logging.info("Submitted {0} new jobs, {1} jobs are running or in the queue".format(
              nj, n_running_jobs))
        elif n_new_windows == 0:
          logging.info("No new windows were generated")
      # Now we update the flags
      if n_running_jobs > 0:
        continue_flag = True
      if njobs < max_jobs and time.time() - t0 < max_time:
        submit_flag = True
      else:
        logging.info(
            "Reached maximum number of jobs or time, no new jobs will be submitted.")
      if save_flag:
        self.system.Save("siPMF_state")
        logging.info("Saving the system.")
      if continue_flag:
        time.sleep(sleep_length)
    # Make sure the PMF is up to date before saving and stopping
    self.system.UpdatePMF(self.environment)
    self.system.Save("siPMF_state")
    logging.info("Stopping.")
