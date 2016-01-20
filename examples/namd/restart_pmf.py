import sys,os
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

#We load the system
from system import LoadSystem
alanine_system=LoadSystem("siPMF_state")

#Make sure the PMF is up to date
alanine_system.UpdatePMF(environment)

#Now we create the SiPMF object and restart the calculation
app=SiPMF(alanine_system,environment)
max_time=86400 #Run for maximum one day
max_jobs=1000 #Run a maximum of 1000 jobs
sleep_length=60 #wait 60 between checks of the status of the jobs
app.Run(max_time,max_jobs,sleep_length)

