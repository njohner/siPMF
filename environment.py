"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This file contains the :class:`Environment` object which represents the computational environment.
"""
import os,subprocess

class Environment():
  """
  This class represents the computational environment in which the software is run.
  It defines the functions used to communicate with the queuing system and the path to
  the WHAM executable.
  """
  def __init__(self,qsub_command,jid_pos,qstat_command,jid_flag,wham_executable):
    """
    :param qsub_command: Command used to submit a job to the queuing system. On SGE this should be "qsub"
    :param jid_pos: Position of the job ID in the string returned by the *qsub_command*
    :param qstat_command: Command used to check the status of a job. On SGE this should be "qstat"
    :param jid_flag: Flag that should be added to the *qstat_command* to check the status of a job with 
     a specific job ID. On SGE this should be "-j"
    :param wham_executable: Path to the wham executable (either the 1D wham or 2D wham, depending on the number of CVs in the system)
    :type qsub_command: :class:`str`
    :type jid_pos: :class:`int`
    :type qstat_command: :class:`str`
    :type jid_flag: :class:`str`
    :type wham_executable: :class:`str`
    """
    self.qsub=self.DefineQsub(qsub_command,jid_pos)
    self.qstat=self.DefineQstat(qstat_command,jid_flag)
    self.wham_executable=wham_executable
  
  def DefineQsub(self,qsub_command,jid_pos):
    def qsub(run_directory,path_to_job_file):
      #os.chdir(run_directory)
      out=subprocess.check_output([qsub_command,path_to_job_file],cwd=run_directory)
      jid=out.split()[jid_pos]
      return jid
    return qsub

  def DefineQstat(self,qstat_command,jid_flag):
    def qstat(jid):
      try:
        out=subprocess.check_output([qstat_command,jid_flag,jid])
      except:
        return "finished"
      return "in queue"
    return qstat

