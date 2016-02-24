"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This file contains the :class:`Phase` object which represents a simulation phase.
"""
import os,subprocess
from job import Job
import logging

class Phase():
  """
  This class is at the bottom of the hierarchical structure used in SiPMF. Every class:`Window`
  can have several :class:`Phase` corresponding to successive simulations.
  """
  def __repr__(self):
    return "Phase({0},{1},{2},{3})".format(self.window,self.name,self.type,self.parent_phase)
  
  def __str__(self):
    cvs=" , ".join([cv.name+"="+str(cv_val) for cv,cv_val in zip(self.window.system.cv_list,self.window.cv_values)])
    if self.parent_phase:
      cvs2=" , ".join([cv.name+"="+str(cv_val) for cv,cv_val in zip(self.parent_phase.window.system.cv_list,self.parent_phase.window.cv_values)])
      return "{0} with {1} and initialized form {2}".format(self.name,cvs,cvs2)
    else:return "{0} with {1}".format(self.name,cvs)
  
  def __init__(self,window,phase_name,phase_type,parent_phase=None):
    """
    :param window:  The window to which the phase belongs
    :param phase_name:  The name of the phase
    :param phase_type:  Either *Initialization* or *run* phase.
    :param parent_phase: The phase from which this one is restarted.
    :type window: :class:`~window.Window`
    :type phase_name: :class:`str`
    :type phase_type: :class:`str`
    :type parent_phase: :class:`~phase.Phase`
    """
    self.window=window
    self.name=phase_name
    self.type=phase_type
    self.n_data=0
    self.parent_phase=parent_phase
    if self.parent_phase:
      self.restartdir=parent_phase.outdir
      #self.restartname=parent_phase.outname
    else:
      self.restartdir=self.window.system.init_restartdir
      #self.restartname=self.window.system.init_restartname
    self.outdir=os.path.join(self.window.subdir,self.name)
    self.path_to_datafile=os.path.join(self.outdir,self.window.system.data_filename)
    #self.outname=self.name

  def Initialize(self):
    """
    Initialize the phase by creating its output directory, and setting up the simulation :class:`Job`, i.e.
    preparing the MD input file.
    """
    logging.info("New phase: {0}".format(self))
    if os.path.isdir(self.outdir):
      logging.error("Directory already exists, program stops to avoid overwriting {0}.".format(self.outdir))
      raise IOError("Directory already exists, program stops to avoid overwriting {0}.".format(self.outdir))
    r=subprocess.call(["mkdir",self.outdir])
    if r!=0:
      logging.error("Problem creating output directory {0}.".format(self.outdir))
      raise IOError("Problem creating output directory {0}.".format(self.outdir))
    self.job=Job(self)

  def UpdateDataCount(self):
    """
    Count how much data has been accumulated in this phase (number of lines in its *datafile*).
    """
    if self.type=="initialization":self.n_data=0
    else:
      if not os.path.isfile(self.path_to_datafile):self.n_data=0
      else:
        f=open(self.path_to_datafile,"r")
        self.n_data=len([1 for line in f if not line.startswith("#")])
        f.close()

  def GetDataCount(self):
    """
    Get the number of data points accumulated for this phase (number of lines in its *datafile*).
    """
    return self.n_data

  def AddDataToWindow(self,n_skip=0,n_tot=-1,new_only=True):
    """
    Add the data of this phase to the data file of the window.
    *n_tot=-1* means there is no maximal number of data points added.

    :param n_skip: The number of data points to skip.
    :param n_tot: The total number of data points used to calculate the PMF.
    :param new_only: Only add the data from phases that have not yet been added to the data files.
    :type n_skip: :class:`int`
    :type n_tot: :class:`int`
    :type new_only: :class:`bool`
    """
    if self.type=="initialization":
      return 
    elif not os.path.isfile(self.path_to_datafile):
      return
    elif new_only and self.data_added_to_window==True:
      return
    f_out=open(self.window.path_to_datafile,"a")
    f_in=open(self.path_to_datafile,"r")
    for line in f_in:
      if line.startswith("#"):continue
      elif self.window.datafile_n_data_skipped<n_skip:
        self.window.datafile_n_data_skipped+=1
      elif n_tot>0 and self.window.datafile_n_data_tot>=n_tot:break
      else:
        f_out.write(line)
        self.window.datafile_n_data_tot+=1
    f_in.close()
    f_out.close()
    self.data_added_to_window=True
    return



