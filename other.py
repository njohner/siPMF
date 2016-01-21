"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This file contains the :class:`CollectiveVariable` and the :class:`PMF` class.
"""
import os
import scipy.interpolate
from scipy.interpolate import NearestNDInterpolator
import matplotlib.pyplot as plt
import numpy as npy

class CollectiveVariable():
  """
  This class represents a collective variable (CV or reaction coordinate).
  The name of the CV will give rise to the corresponding fields that are replaced
  in the MD input files. The minimal and maximal values determine the range in which
  the windows can be generated. The number of bins are used for the calculation of the PMF.
  The periodicity of the collective variables is used both for generating new windows and for the
  calculation of the PMF. The step size defines the distance between two neighboring windows.
  Finally the units are just used in the labels of the plots of the PMF.
  """
  def __repr__(self):
    return "CollectiveVariable({0},{1},{2},{3},{4},{5})".format(self.name,self.min_value,self.max_value,self.step_size,self.num_bins,self.periodicity)

  def __init__(self,cv_name,min_value,max_value,step_size,num_bins,periodicity=None,units=""):
    """
    :param cv_name: The name of the CV
    :param min_value: The minimal value of the CV
    :param max_value: The maximal value of the CV
    :param step_size: The size of the step for this CV
    :param num_bins: number of bins
    :param periodicity: periodicity of the CV
    :param units: the units of the CV
    :type cv_name: :class:`str`
    :type min_value: :class:`float`
    :type max_value: :class:`float`
    :type step_size: :class:`float`
    :type num_bins: :class:`int`
    :type periodicity: :class:`float`
    :type units: :class:`str`
    """
    self.name=cv_name
    self.min_value=min_value
    self.max_value=max_value
    self.num_bins=num_bins
    self.step_size=step_size
    self.bin_size=(self.max_value-self.min_value)/float(num_bins)
    self.periodicity=periodicity
    self.units=units
  
class PMF():
  def __repr__(self):
    return "PMF({0},{1})".format(self.points,self.values)

  def __init__(self,points,values,cv_list):
    """
    :param points: values of the CVs at which the PMF is found in values.
    :param values: Free energy corresponding to CV values in points
    :param cv_list: List of the CVs 
    :type points: :class:`numpy.array`
    :type values: :class:`numpy.array`
    :type cv_list: :class:`list`(:class:`CollectiveVariable`)
    """
    self.points=points
    self.values=values
    self.cv_list=cv_list
    self.dimensionality=len(self.cv_list)
    self.interpolator=scipy.interpolate.interpnd.LinearNDInterpolator(self.points,self.values)
    #self.interpolator=NearestNDInterpolator(self.points,self.values)
  def GetValue(self,point):
    """
    Free energy at a certain position on the free energy surface.

    :param point: values of the CV for which we want the free energy.
    """
    return self.interpolator(point)

  def Plot(self,outputdir,filename,n_levels=20,max_E=None,windows=None,energy_units=""):
    """
    Plot the PMF.
    """
    if self.dimensionality==2:
      X=self.points[:,0]
      Y=self.points[:,1]
      Z=self.values
      if max_E:Z=npy.array([min(el,max_E) for el in Z])
      nb_x=self.cv_list[0].num_bins
      nb_y=self.cv_list[1].num_bins
      X=X.reshape([nb_x,nb_y])
      Y=Y.reshape([nb_x,nb_y])
      Z=Z.reshape([nb_x,nb_y])
      plt.figure()
      plt.contourf(X,Y,Z,n_levels)
      plt.colorbar(label=energy_units)
      if self.cv_list[0].units:plt.xlabel("{0} [{1}]".format(self.cv_list[0].name,self.cv_list[0].units))
      else:plt.xlabel("{0}".format(self.cv_list[0].name))
      if self.cv_list[1].units:plt.ylabel("{0} [{1}]".format(self.cv_list[1].name,self.cv_list[1].units))
      else:plt.ylabel("{0}".format(self.cv_list[1].name))
      if windows:
        p1=self.cv_list[0].periodicity
        p2=self.cv_list[1].periodicity
        if not p1:p1=0
        if not p2:p2=0
        for w in windows:
          if w.parent:
            p=w.parent
            x=p.cv_values[0]
            y=p.cv_values[1]
            dx=w.cv_values[0]-x
            dy=w.cv_values[1]-y
            if p1 and dx>p1/2.:dx=dx-p1
            if p1 and dx<-p1/2.:dx=dx+p1
            if p2 and dy>p2/2.:dy=dy-p2
            if p2 and dy<-p2/2.:dy=dy+p2
            plt.arrow(x,y,dx,dy,length_includes_head=True,head_length=3,head_width=3)
      plt.savefig(os.path.join(outputdir,filename+".png"))
      plt.close()
    elif self.dimensionality==1:
      plt.figure()
      plt.plot(pmf.points,pmf.values)
      if max_E:plt.ylim([0,max_E])
      if self.cv_list[0].units:plt.xlabel("{0} [{1}]".format(self.cv_list[0].name,self.cv_list[0].units))
      else:plt.xlabel("{0}".format(self.cv_list[0].name))
      plt.ylabel("Free Energy")
      plt.savefig(os.path.join(outputdir,filename+".png"))
      plt.close()
