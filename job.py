import os

class Job():
  def __repr__(self):
    return "Job({0})".format(self.phase)

  def __init__(self,phase):
    self.phase=phase
    self.status="To submit"
    self.GenerateInputFile()
    if self.phase.type=="initialization":self.path_to_job_file=os.path.join(self.phase.window.system.basedir,self.phase.window.system.init_job_fname)
    elif self.phase.type=="run":self.path_to_job_file=os.path.join(self.phase.window.system.basedir,self.phase.window.system.run_job_fname)
  
  def GenerateInputFile(self):
    if self.phase.type=="initialization":
      f=open(self.phase.window.system.GetPathToInitInputFile(),"r")
      replace=self.ReplaceInInitInput
      self.path_to_job_file=os.path.join(self.phase.outdir,self.phase.window.system.init_input_fname)
    elif self.phase.type=="run":
      f=open(self.phase.window.system.GetPathToRunInputFile(),"r")
      replace=self.ReplaceInRunInput
      self.path_to_job_file=os.path.join(self.phase.outdir,self.phase.window.system.run_input_fname)
    f_data=f.readlines()
    f.close()
    f_data=replace(f_data)
    outf=open(self.path_to_job_file,"w")
    outf.write("".join(f_data))
    outf.close()

  def ReplaceInInput(self,f_data):
    to_replace={"{BASEDIR}":self.phase.window.system.basedir,"{RESTARTDIR}":self.phase.restartdir}#,"_RESTARTNAME":self.phase.restartname}
    to_replace.update({"{OUTPUTDIR}":self.phase.outdir})#,"_OUTPUTNAME":self.outname})
    for cvn,cvv in zip(self.phase.window.cv_names,self.phase.window.cv_values):to_replace["{"+cvn+"}"]=cvv
    for cvn,cvk in zip(self.phase.window.cv_names,self.phase.window.spring_constants):to_replace["{"+cvn+"_K}"]=cvk
    f_data_new=[]
    for line in f_data:
      for key in to_replace:line=line.replace(key,str(to_replace[key]))
      f_data_new.append(line)
    return f_data_new

  def ReplaceInInitInput(self,f_data):
    to_replace={"{INIT_NSTEP}":self.phase.window.system.init_nstep}
    f_data=self.ReplaceInInput(f_data)
    for cvn,cvv in zip(self.phase.window.parent.cv_names,self.phase.window.parent.cv_values):to_replace["{PARENT_"+cvn+"}"]=cvv
    for cvn,cvk in zip(self.phase.window.parent.cv_names,self.phase.window.parent.spring_constants):to_replace["{PARENT_"+cvn+"_K}"]=cvk
    f_data_new=[]
    for line in f_data:
      for key in to_replace:line=line.replace(key,str(to_replace[key]))
      f_data_new.append(line)
    return f_data_new

  def ReplaceInRunInput(self,f_data):
    to_replace={"{RUN_NSTEP}":self.phase.window.system.run_nstep}
    f_data=self.ReplaceInInput(f_data)
    f_data_new=[]
    for line in f_data:
      for key in to_replace:line=line.replace(key,str(to_replace[key]))
      f_data_new.append(line)
    return f_data_new

  def Submit(self,environment):
    self.jid=environment.qsub(self.phase.outdir,self.path_to_job_file)
    self.status="submitted"
    self.phase.window.system.unfinished_jobs.append(self)

  def UpdateStatus(self,environment):
    if self.status!="finished":
      self.status=environment.qstat(self.jid)
      if self.status=="finished" and self.phase.type=="run":
        self.phase.window.n_run_phases+=1


