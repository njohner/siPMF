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

  def __init__(self,phase):
    """
    :param phase: The corresponding simulation phase.
    :type phase:  :class:`Phase`
    """
    self.phase=phase
    self.queue_status="To submit"
    self.success=None
    self.GenerateInputFile()
    if self.phase.type=="initialization":self.path_to_job_file=os.path.join(self.phase.window.system.basedir,self.phase.window.system.init_job_fname)
    elif self.phase.type=="run":self.path_to_job_file=os.path.join(self.phase.window.system.basedir,self.phase.window.system.run_job_fname)
  
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
    to_replace=self.GetInputReplacementDict()
    if self.phase.type=="initialization":
      f=open(self.phase.window.system.GetPathToInitInputFile(),"r")
      to_replace.update(self.GetInitInputReplacementDict())
      self.path_to_job_file=os.path.join(self.phase.outdir,self.phase.window.system.init_input_fname)
    elif self.phase.type=="run":
      f=open(self.phase.window.system.GetPathToRunInputFile(),"r")
      to_replace.update(self.GetRunInputReplacementDict())
      self.path_to_job_file=os.path.join(self.phase.outdir,self.phase.window.system.run_input_fname)
    f_data=f.readlines()
    f.close()
    f_data_new=[]
    for line in f_data:
      for key in to_replace:line=line.replace(key,str(to_replace[key]))
      f_data_new.append(line)
    outf=open(self.path_to_job_file,"w")
    outf.write("".join(f_data_new))
    outf.close()

  def GetInputReplacementDict(self):
    """
    This function returns a dictionary containing the fields that will be replaced both in the run and in the
    initialization MD input files.
    """
    to_replace={"{BASEDIR}":self.phase.window.system.basedir,"{RESTARTDIR}":self.phase.restartdir}#,"_RESTARTNAME":self.phase.restartname}
    to_replace.update({"{OUTPUTDIR}":self.phase.outdir})#,"_OUTPUTNAME":self.outname})
    for cvn,cvv in zip(self.phase.window.cv_names,self.phase.window.cv_values):to_replace["{"+cvn+"}"]=cvv
    for cvn,cvk in zip(self.phase.window.cv_names,self.phase.window.spring_constants):to_replace["{"+cvn+"_K}"]=cvk
    return to_replace

  def GetInitInputReplacementDict(self):
    """
    This function returns a dictionary containing the fields that will be replaced only in the
    initialization MD input files.
    """
    to_replace={"{INIT_NSTEP}":self.phase.window.system.init_nstep}
    for cvn,cvv in zip(self.phase.window.parent.cv_names,self.phase.window.parent.cv_values):to_replace["{PARENT_"+cvn+"}"]=cvv
    for cvn,cvk in zip(self.phase.window.parent.cv_names,self.phase.window.parent.spring_constants):to_replace["{PARENT_"+cvn+"_K}"]=cvk
    return to_replace

  def GetRunInputReplacementDict(self):
    """
    This function returns a dictionary containing the fields that will be replaced only in the
    run MD input files.
    """
    to_replace={"{RUN_NSTEP}":self.phase.window.system.run_nstep}
    return to_replace

  def Submit(self,environment):
    """
    Submit the job to the cluster

    :param environment: The environment used to submit the job
    :type environment: :class:`Environment`
    """
    self.jid=environment.qsub(self.phase.outdir,self.path_to_job_file)
    self.status="submitted"
    self.phase.window.system.unfinished_jobs.append(self)

  def UpdateStatus(self,environment):
    """
    Check whether the job is still in the queue

    :param environment: The environment used to check the job status
    :type environment: :class:`Environment`
    """
    if self.queue_status!="finished":
      self.queue_status=environment.qstat(self.jid)
      if self.queue_status=="finished":
        self.success=True
        for fname in self.phase.window.system.check_fnames:
          if not os.path.isfile(os.path.join(self.phase.outdir,fname)):
            self.success=False
            break
        if self.success and self.phase.type=="run":
          self.phase.window.n_run_phases+=1



