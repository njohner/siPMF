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
------------------------------------------------------
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

Usage
--------------



Classes
------------
Each class and their methods are described in the code using docstrings, which can be accessed in python with the *help()* command or found in the documentation pages for each object below:

- The :doc:`SiPMF <si_pmf>` class is the object used to control the calculation. It defines the process that will run in the background to periodically check the simulations, decide whether to create new windows or submit new jobs. 
- The :doc:`Environment <environment>` class represents the environment of the cluster and defines communication with the queuing system and WHAM.
- The :doc:`System <system>` class contains all the information about the simulated system
- The :doc:`Window <window>` class represents a simulation window
- The :doc:`Phase <phase>` class represents a simulation phase in a **Window**
- The :doc:`Job <job>` class represents a job on the HPC cluster.
- The :doc:`PMF <pmf>` class represents the PMF
- The :doc:`CollectiveVariable <collective_variable>` class represents a collective variable.


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

