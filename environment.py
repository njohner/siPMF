import os,subprocess

class Environment():
  def __init__(self,qsub_command,jid_pos,qstat_command,jid_flag,wham_executable):
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

