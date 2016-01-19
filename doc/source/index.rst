.. SiPMF documentation master file, created by
   sphinx-quickstart on Fri Jan 15 15:38:56 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

******************************************************************************************
A simple implementation of the self-learning adaptative umbrella sampling method 
******************************************************************************************

Introduction
==============
Molecular dynamics (MD) simulations offer the possibility to study biological mechanisms in much more details than any experimental technique. Nevertheless obtaining relevant information from MD simulations is challenging and often unfeasible with simple simulations of the system of interest. One reason for this is that most processes of interest happen on timescales similar or longer than the simulation time within reach with current HPC clusters, leading to insufficient data and therefore results that do not have any statistical significance. 
Many strategies to accelerate sampling or to calculate the free energy landscape governing the process of interest have therefore been developped. Enhanced sampling methods are typically used for exploratory work when there is very little knowledge of the conformational changes associated with the studied process, whereas free-energy calculations are much more efficient but presuppose a good prior understanding of the process. Indeed these approaches require the definition of collective-variables (CV, also called reaction coordinates) which describe the process under study and can be used to calculate the potential of mean force (PMF, or free energy landscape) governing the biological process.
Most modern methods to calculate PMFs (such as metadynamics or Adaptative Biasing Force) are adaptative, in that they create a representation of the PMF on the fly and use it to bias the ongoing simulation to improve the sampling of the free-energy surface. Their advantage are typically that the simulations are easy to setup, requiring little input apart from the definition of the CVsm and that they concentrate the sampling on relevant regions (regions of low free-energy). Their disadvantage is that they typically require to cross the free energy surface several times to converge. An older method, called umbrella sampling, does not suffer the same problem as the PMF is calculated by combining the data obtained from many simulations restrained to small portions of the free-energy surface (called windows), which avoids the need to diffuse on the free-energy surface. This method is also embarassingly parallelizable as the simulations of different windows are independent. The downside is that the starting configurations for each window have to be generated, which can be tedious, and that the windows typically cover the full free-energy landscape instead of only the regions of interest. The self-learning adaptative umbrella sampling method combines the best of both worlds by automatizing the generation of new windows and restricting it to regions of low-free energy by using a partial PMF obtained from the already simulated windows to decide in which directions to expand the exploration of the free energy surface.  A complete description of the method can be found in ref. [1]_, but a quick summary, which should suffice to understand how the method works and its implementation, is provided below.

Method
============
self-learning adaptative umbrella sampling
--------------------------------------------
The self-learning adaptative umbrella sampling method is based on umbrella sampling, which uses a restraining potential applied to the collective variables to focus the sampling to a particular region of the free energy landscape, called a window. The data obtained from many such windows, which have to partially overlap, is combined to calculate the PMF, using the wheighted-histogram analysis method (WHAM). Using many narrow windows was shown to be the optimal approach for calculating PMFs. The downside of the classical umbrella sampling is that the generation of the windows can be time consuming and that simply setting up windows to explore the full free energy landscape leads to simulation time wasted in regions of very high free energy which are usually irrelevant for the understanding of the studied system.
The adaptative umbrella sampling method alleviates both these issues by automatizing the generation of windows and limiting the exploration of the free energy landscape to regions of low free energies. This is achieved with a simple procedure which can be summarized as follows:
1. Simulate the system in a given window
2. Calculate the PMF from the simulated windows using WHAM (partial PMF).
3. Construct new windows around the ones corresponding to regions of low free-energies.
4. Repeat (simulate the system in the newly constructed windows, recalculate the PMF, generate new windows, etc.)

Convergence speedup using a linear exploration scheme
-------------------------------------------------------
One important problem with free-energy calculations is the existence of slow degrees of freedom, orthogonal to the CVs studied. Typically for any free-energy calculation to yield meaningfull results, you have to move slowly enough in the CV space that all other degrees always remain at equilibrium. In umbrella sampling, this means that when a new window is generated (i.e. the center of the constraining potentials is moved), some time will be needed to allow the system to equilibrate before data is collected to calculate the PMF. This equilibration time will be longer is the system is further from equilibrium. Hence our implementation only proceeds to create new windows once the sampling of the window from which they are created (called the *parent window*) is finished. This should ensure that the *parent window* is well equilibrated so that the child window should not be too much out of equilibrium. Moreover the child window is created by slowly moving the constraining potential form its position in the parent window to its new position. This is called the *initialization phase*, whereas the following simulations performed for an initialized window, which have fixed constraining potentials are called *run phases*.
Clearly it is also advised to use small windows, i.e. high constraining potentials and small steps in the CV space when creating a new window.

Implementation
===================

Dependencies
--------------
The method is implemented entirely in **python 2.7** but relies on several standard scientific python packages, notably numpy, scipy, matplotlib and pickle. It also relies on an external software for the WHAM calculations, using the **implementation of WHAM by Alan Grossfield** [2]_.

Applicability
---------------
The implementation of WHAM used [2]_, effectively limits the use of the current tool to **1 and 2-dimensional systems**.
The implementation should be general enough to allow support for different MD packages such as **NAMD** (tested), **CHARMM** (untested) and **Gromacs** (untested). 
The implementation is meant for use on HPC clusters, with simulations being submitted to a queuing system, and should be adaptable to different systems. For now it has been tested in **SGE** and will shortly be tested on **SLURM**. 


Code architecture and classes
-------------------------------
As described in the **Method** section, umbrella sampling uses many simulation *Windows*, each corresponding to a simulation restrained to a small region of the CV space. The data from these windows can then be combined to determine the PMF. The generation and sampling of each *Window* is accomplished with successive simulations, called *Phases*. Of course, each *Phase* uses as initial conditions the final state of another simulation (typically the positions, velocities and periodic boundaries), we call this the *Parent Phase*. 

New *Windows* are generated from neighboring *Windows*, so that the change in restraining potential is small and the system remains as close to equilibrium as possible. The *Window* from which a new *Window* is generated is called its *Parent Window* and the simulation used to generate the new *Window* (or *Child Window*) is called the *Initialization Phase*. During this *Initialization Phase* the position of the constraining potential is moved from its position in the *Parent Window* to its position in the *Child Window*. So when generating a new *Window*, the *parent phase* for the *initialization phase* will be the last *phase* of the *parent window*.

This Hierarchical structure is translated both in the python objects and in the data structure (directory tree) generated by the software. In the following we describe these classes and give the names of some of their attributes and methods in brackets. The object at the top of the hierarchy is the :doc:`System <system>` class, which contains all the information about the ongoing calculation, notably the path to the directory where the calculations are done (*basedir*), the list of collective variables (*cv_list*), filenames of the different template files (MD input files, job submision files) and importantly it contains a list of all the *Windows* in the *System*. The :doc:`Window <window>` class is the next object in the hierarchy and represents a simulation window. It notably contains the path to corresponding subdirectory (*subdir*) which is *basedir/Window.name*. Among its other attributes there are the values of the CVs (center of the restraining potentials) and the spring constants used for that *Window* as well as a reference to its *parent window* and to the object above it in the hierarchy, the *System*. Of course a :doc:`Window <window>` also contains a list of the *phases* that were run for that *Window*. So the :doc:`Phase <phase>` class is the next object of the hierarchy and represents a simulation phase in a **Window**. Again it contains the path to the output directory (*outdir*) which is *subdir/Phase.name*, whereto the corresponding simulation output should be written (notably the *datafile* and the files needed to restart the next simulation). A *Phase* also has a reference to its *Window* (*window*) and to its *parent phase* (*parent_phase*). Finally, each *Phase* corresponds to a :doc:`Job <job>`, which represents a job on the HPC cluster. The *Job* class has methods to submit itself on the cluster (*Submit*) and check whether the job is still on the cluster or not (*UpdateStatus*). It also has methods to generate the MD input files needed to run the simulation (*GenerateInputFile*).

Apart from this hierarchy, there are 3 more classes in *SiPMF*. First the `CollectiveVariable <collective_variable>` class, representing a CV. It has a name (*name*), which is used as a field that will be replaced in the MD inputs by the value of the CV for a particular window. The it has a range (*min_value* and *max_value*) and a period (*periodicity*), number of bins used in the PMF calculation (*num_bins*) and the size of the steps used for that CV when generating a new window (*step_size*). The second object is the :doc:`PMF <pmf>` class, which represents the free energy landscape. The *PMF* is a list of points in the CV space (*points*) and the corresponding free energy (*values*). It also has an *interpolator* which is used to return the free energy of any point in the CV space (*GetValue*). The last class is the :doc:`Environment <environment>` class represents the environment of the cluster and defines communication with the queuing system and WHAM. Its main methods are *qsub* and *qstat* in reference to the corresponding SGE commands for submitting a job and checking the status of a job. As an attribute it has the path to the WHAM executable (*wham_executable*).

Each class and their methods are described in the code using docstrings, which can be accessed in python with the *help()* command or found in the corresponding documentation pages below:

- The :doc:`SiPMF <si_pmf>` class is the object used to control the calculation. It defines the process that will run in the background to periodically check the simulations, decide whether to create new windows or submit new jobs. 
- The :doc:`Environment <environment>` class represents the environment of the cluster and defines communication with the queuing system and WHAM.
- The :doc:`System <system>` class contains all the information about the simulated system
- The :doc:`Window <window>` class represents a simulation window
- The :doc:`Phase <phase>` class represents a simulation phase in a **Window**
- The :doc:`Job <job>` class represents a job on the HPC cluster.
- The :doc:`PMF <pmf>` class represents the PMF
- The :doc:`CollectiveVariable <collective_variable>` class represents a collective variable.


Usage
==================

Description
--------------
For each window and for each phase, a specific MD input file has to be generated, making sure that the potentials constraining the CVs are centered at the right values and that the simulation starts with the correct positions and velocities. This is achieved by replacing fields in a template input file provided by the user with the appropriate values. These fields are delimited with curly braces in the input file, for example *{RESTARTDIR}* will be replaced by the path pointing to the directory of the parent phase.


Installation
=====================

Be aware that these modules depend on several standard scientific python packages, notably numpy,
scipy, matplotlib and pickle. It also relies on an external software for the WHAM calculations, using the implementation of WHAM by Alan Grossfield [2]_.





Bibliography
---------------------------------------
.. [1] Wojtas-Niziurski, Wojciech, Yilin Meng, Benoı̂t Roux, and Simon Bernèche. “Self-Learning Adaptive Umbrella Sampling Method for the Determination of Free Energy Landscapes in Multiple Dimensions.” Journal of Chemical Theory and Computation 9, no. 4 (2013): 1885–95.

.. [2] Grossfield, Alan, "WHAM: the weighted histogram analysis method", Version 2.0.9

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

