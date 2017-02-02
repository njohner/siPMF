This is an implementation of the self-learning adaptative umbrella sampling method (ref. [1]),
used to calculate free energy landscapes from molecular dynamics (MD) simulations.

This implementation of self-learning adaptative umbrella sampling is meant for calculating a PMF on an HPC cluster. 
It will automatically generate simulation Windows and Phases (see the Implementation section for details) 
and submit the corresponding jobs to the queuing system of the cluster. 
It will periodically check the status of the submitted jobs, calculate the current PMF and 
decide whether to generate new windows and/or submit new jobs.

The tool is implemented in Python 2.7 and works with NAMD and CHARMM. It should also work
with other MD simulation packages such as GROMACS, although this has not been tested yet.

Documentation
----------------------
The complete documentation is included in docstrings in each python file. It can be found online at http://njohner.github.io/siPMF/ 
or it can be assembled into html files using sphinx from the source files. To build the documentation, simply go to the *doc* folder and :

```shell
make html
```
This should generate a set of html files in *doc/build/html*. To view the documentation, open the *index.html* file.

Examples on how tu use this tool with NAMD and CHARMM can be found in the *example* folder.
 

References
----------------------

[1] Wojtas-Niziurski, Wojciech, Yilin Meng, Benoı̂t Roux, and Simon Bernèche. "Self-Learning Adaptive Umbrella Sampling Method for the Determination of Free Energy Landscapes in Multiple Dimensions." Journal of Chemical Theory and Computation 9, no. 4 (2013): 1885–95.


