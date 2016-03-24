import sys,os,logging
basedir=os.getcwd()
sys.path.append(os.path.join(basedir,"../.."))
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
cv1=CollectiveVariable(name,min_value,max_value,step_size,num_bins,periodicity)
name="PSI"
cv2=CollectiveVariable(name,min_value,max_value,step_size,num_bins,periodicity)

#Now we define the system
cv_list=[cv1,cv2]
init_input_fname="initialize.in"
run_input_fname="run.in"
init_job_fname="namd_init.lsf"
run_job_fname="namd_run.lsf"
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

#Initialize the system, i.e. generate the starting window.
cv_values=[90,90]
spring_constants=[0.05,0.05]
init_restartdir=os.path.join(basedir,"../data/equilibration")
alanine_system.Initialize(cv_values,spring_constants,init_restartdir)

#Now we create the SiPMF object and start the calculation
app=SiPMF(alanine_system,environment)
max_time=86400 #Run for maximum one day
max_jobs=1000 #Run a maximum of 1000 jobs
sleep_length=60 #wait 60 between checks of the status of the jobs
app.Run(max_time,max_jobs,sleep_length)

