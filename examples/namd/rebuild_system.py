import sys,os,logging
basedir=os.getcwd()
sys.path.append(os.path.join(basedir,"../../.."))
from siPMF import *

#Setup the log file. This has to be done before the other objects are created.
logging.basicConfig(filename=os.path.join(basedir,"siPMF.log"),level=logging.INFO,format='%(asctime)s: %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M %p')

#We define the environment:
qsub_command="qsub"
jid_pos=2
qstat_command="qstat"
jid_flag="-j"
wham_executable="path/to/wham"
environment=Environment(qsub_command,jid_pos,qstat_command,jid_flag,wham_executable)

#We define the collective variables
name="PHI"
min_value=-180
max_value=180
step_size=15
num_bins=72
periodicity=360
spring_constant=0.05
cv1=CollectiveVariable(name,min_value,max_value,step_size,num_bins,spring_constant,periodicity=periodicity)
name="PSI"
cv2=CollectiveVariable(name,min_value,max_value,step_size,num_bins,spring_constant,periodicity=periodicity)

#Now we define the system
cv_list=[cv1,cv2]
init_input_fname="initialize.in"
run_input_fname="run.in"
init_job_fname="submit_init.lsf"
run_job_fname="submit_run.lsf"
data_filename="alanine.colvars.traj"
init_nstep=1000
run_nstep=5000
n_data=run_nstep/15
max_E1=1.0
max_E2=6.0
temperature=303
check_fnames=["alanine.restart.coor","alanine.restart.vel","alanine.restart.xsc"]
target_cv_vals=[(105,130)]
alanine_system=System(basedir,cv_list,init_input_fname,run_input_fname,init_job_fname,run_job_fname,data_filename,init_nstep,run_nstep,n_data,max_E1,max_E2,temperature,check_fnames,target_cv_vals)
alanine_system.init_restartdir=os.path.join(basedir,"../data/equilibration")
from system import RebuildWindowsAndPhasesFromDirectoryTree
RebuildWindowsAndPhasesFromDirectoryTree(alanine_system)
alanine_system.UpdateDataFiles(new_only=False)
alanine_system.UpdateDataCounts()


