"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This file contains the :class:`Job` object which represents a job on a cluster.
"""
import os


class Job():
  """
  This class is the bottom of the hierarchical structure used in SiPMF. There is one :class:`Job`
  object for every :class:`Phase` of every :class:`Window`. This object is used to submit a job
  to the cluster and check its status
  """

  def __repr__(self):
    return "Job({0})".format(self.phase)

  def __init__(self, phase):
    """
    :param phase: The corresponding simulation phase.
    :type phase:  :class:`~phase.Phase`
    """
    self.phase = phase
    self.queue_status = "To submit"
    self.success = None
    self.GenerateInputFile()
    self.GenerateJobFile()

  def GenerateJobFile(self):
    """
    Generate the Job submission file. This function reads the template job file (either the *initialization job file*
    or the *run job file*) and replaces several fields by their corresponding values, notably
    {BASEDIR}, {RESTARTDIR}, {OUTPUTDIR}, {WINDOW}, {PHASE} and {INPUTFILE}. For the initialization job submission file we also replace fields for the parent window, namely
    {PARENT_WINDOW} and {PARENT_PHASE}.
    """
    to_replace = self.GetJobReplacementDict()
    if self.phase.type == "initialization":
      f = open(self.phase.window.system.GetPathToInitJobFile(), "r")
      to_replace.update(self.GetInitJobReplacementDict())
      self.path_to_job_file = os.path.join(
          self.phase.outdir, self.phase.window.system.init_job_fname)
    elif self.phase.type == "run":
      f = open(self.phase.window.system.GetPathToRunJobFile(), "r")
      to_replace.update(self.GetRunJobReplacementDict())
      self.path_to_job_file = os.path.join(
          self.phase.outdir, self.phase.window.system.run_job_fname)
    f_data = f.readlines()
    f.close()
    f_data_new = []
    for line in f_data:
      for key in to_replace:
        line = line.replace(key, str(to_replace[key]))
      f_data_new.append(line)
    outf = open(self.path_to_job_file, "w")
    outf.write("".join(f_data_new))
    outf.close()

  def GetJobReplacementDict(self):
    """
    This function returns a dictionary containing the fields that will be replaced both in the run and in the
    initialization Job submission files.
    """
    to_replace = {"{BASEDIR}": self.phase.window.system.basedir,
                  "{RESTARTDIR}": self.phase.restartdir,
                  "{OUTPUTDIR}": self.phase.outdir,
                  "{WINDOW}": self.phase.window.name,
                  "{PHASE}": self.phase.name,
                  "{INPUTFILE}": self.path_to_input_file}
    return to_replace

  def GetInitJobReplacementDict(self):
    """
    This function returns a dictionary containing the fields that will be replaced only in the
    initialization Job submission files.
    """
    to_replace = {"{PARENT_WINDOW}": self.phase.window.parent.name,
                  "{PARENT_PHASE}": self.phase.parent_phase.name}
    return to_replace

  def GetRunJobReplacementDict(self):
    """
    This function returns a dictionary containing the fields that will be replaced only in the
    run Job submission files.
    """
    to_replace = {}
    return to_replace

  def GenerateInputFile(self):
    """
    Generate the MD input file. This function reads the template input file (either the *initialization input file*
    or the *run input file*) and replaces several fields by their corresponding values, notably:
    {BASEDIR},{RESTARTDIR},{OUTPUTDIR} as well as the values and spring constants of the CVs. Specifically
    for each CV a field containing its name ({cvname}) is replaced by its value and {cvname_K} is replaced 
    by the spring constant. For the initialization we also replace fields for the collective variables
    from the parent phase {PARENT_cvname} by the value.
    Finally in the initialization input, {INIT_NSTEP} will be replaced by its value (system.init_nstep)
    and for a run phase {RUN_NSTEP} is replaced by system.run_nstep.
    """
    to_replace = self.GetInputReplacementDict()
    if self.phase.type == "initialization":
      f = open(self.phase.window.system.GetPathToInitInputFile(), "r")
      to_replace.update(self.GetInitInputReplacementDict())
      self.path_to_input_file = os.path.join(
          self.phase.outdir, self.phase.window.system.init_input_fname)
    elif self.phase.type == "run":
      f = open(self.phase.window.system.GetPathToRunInputFile(), "r")
      to_replace.update(self.GetRunInputReplacementDict())
      self.path_to_input_file = os.path.join(
          self.phase.outdir, self.phase.window.system.run_input_fname)
    f_data = f.readlines()
    f.close()
    f_data_new = []
    for line in f_data:
      for key in to_replace:
        line = line.replace(key, str(to_replace[key]))
      f_data_new.append(line)
    outf = open(self.path_to_input_file, "w")
    outf.write("".join(f_data_new))
    outf.close()

  def GetInputReplacementDict(self):
    """
    This function returns a dictionary containing the fields that will be replaced both in the run and in the
    initialization MD input files.
    """
    to_replace = {"{BASEDIR}": self.phase.window.system.basedir,
                  "{RESTARTDIR}": self.phase.restartdir,
                  "{OUTPUTDIR}": self.phase.outdir,
                  "{TEMPERATURE}": self.phase.window.system.temperature}
    for cvn, cvv, cvs in zip(self.phase.window.cv_names, self.phase.window.cv_values, self.phase.window.cv_shifts):
      to_replace["{" + cvn + "}"] = cvv + cvs
    for cvn, cvk in zip(self.phase.window.cv_names, self.phase.window.spring_constants):
      to_replace["{" + cvn + "_K}"] = cvk
    return to_replace

  def GetInitInputReplacementDict(self):
    """
    This function returns a dictionary containing the fields that will be replaced only in the
    initialization MD input files.
    """
    to_replace = {"{INIT_NSTEP}": self.phase.window.system.init_nstep}
    for pcvn, pcvv, pcvs, cvv, cvs, cv in zip(self.phase.window.parent.cv_names, self.phase.window.parent.cv_values, self.phase.window.parent.cv_shifts, self.phase.window.cv_values, self.phase.window.cv_shifts, self.phase.window.system.cv_list):
      if not cv.periodicity:
        to_replace["{PARENT_" + pcvn + "}"] = pcvv + pcvs
      else:
        # Careful with periodic CVs, not to pull accross the whole range when crossing a boundary
        if abs(cvv - pcvv) < cv.periodicity / 2.0:
          to_replace["{PARENT_" + pcvn + "}"] = pcvv + pcvs
        elif cvv - pcvv < -cv.periodicity / 2.0:
          to_replace["{PARENT_" + pcvn + "}"] = pcvv + pcvs - cv.periodicity
        elif cvv - pcvv > cv.periodicity / 2.0:
          to_replace["{PARENT_" + pcvn + "}"] = pcvv + pcvs + cv.periodicity
    for cvn, cvk in zip(self.phase.window.parent.cv_names, self.phase.window.parent.spring_constants):
      to_replace["{PARENT_" + cvn + "_K}"] = cvk
    return to_replace

  def GetRunInputReplacementDict(self):
    """
    This function returns a dictionary containing the fields that will be replaced only in the
    run MD input files.
    """
    to_replace = {"{RUN_NSTEP}": self.phase.window.system.run_nstep}
    return to_replace

  def Submit(self, environment):
    """
    Submit the job to the cluster

    :param environment: The environment used to submit the job
    :type environment: :class:`~environment.Environment`
    """
    self.jid = environment.qsub(self.phase.outdir, self.path_to_job_file)
    self.status = "submitted"
    self.phase.window.system.unfinished_jobs.append(self)

  def UpdateStatus(self, environment):
    """
    Check whether the job is still in the queue

    :param environment: The environment used to check the job status
    :type environment: :class:`~environment.Environment`
    """
    if self.queue_status != "finished":
      self.queue_status = environment.qstat(self.jid)
      if self.queue_status == "finished":
        self.success = True
        for fname in self.phase.window.system.check_fnames:
          if not os.path.isfile(os.path.join(self.phase.outdir, fname)):
            self.success = False
            break
        if self.success and self.phase.type == "run":
          self.phase.window.n_run_phases += 1
