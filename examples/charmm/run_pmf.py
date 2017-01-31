import matplotlib
matplotlib.use("Agg")
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
wham_executable="/scicore/home/schwede/johner/wham/wham/wham"
environment=Environment(qsub_command,jid_pos,qstat_command,jid_flag,wham_executable)

#We define the collective variables
name="Z3"
min_value=-15
max_value=15
step_size=1.0
num_bins=310
periodicity=None
min_k=15.
max_k=15.
max_shift=0.
cv1=CollectiveVariable(name,min_value,max_value,step_size,num_bins,min_k,max_k,max_shift,periodicity,units="A")

#Now we define the system
cv_list=[cv1]
init_input_fname="initialize.in"
run_input_fname="run.in"
init_job_fname="submit_init.lsf"
run_job_fname="submit_run.lsf"
data_filename="pot_pmf.dat"
init_nstep=1000
run_nstep=5000
n_data=run_nstep/15
max_E1=90.0
max_E2=100.0
temperature=303
check_fnames=["pot_pmf.crd","pot_pmf.rst"]
alanine_system=System(basedir,cv_list,init_input_fname,run_input_fname,init_job_fname,run_job_fname,data_filename,init_nstep,run_nstep,n_data,max_E1,max_E2,temperature,check_fnames)

#Initialize the system, i.e. generate the starting window.
spring_constants=[15.0]
init_restartdir=os.path.join(basedir,"equilibration")
cv_values=[0.0]
alanine_system.Initialize(cv_values,spring_constants,init_restartdir)

#Now we create the SiPMF object and start the calculation
app=SiPMF(alanine_system,environment)
max_time=86400 #Run for maximum one day
max_jobs=1000 #Run a maximum of 1000 jobs
sleep_length=60 #wait 60 between checks of the status of the jobs
app.Run(max_time,max_jobs,sleep_length)

