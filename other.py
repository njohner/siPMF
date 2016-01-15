import scipy.interpolate
from scipy.interpolate import NearestNDInterpolator
import matplotlib.pyplot as plt

class CollectiveVariable():
  def __repr__(self):
    return "CollectiveVariable({0},{1},{2},{3},{4},{5})".format(self.name,self.min_value,self.max_value,self.step_size,self.num_bins,self.periodicity)

  def __init__(self,cv_name,min_value,max_value,step_size,num_bins,periodicity=None):
    self.name=cv_name
    self.min_value=min_value
    self.max_value=max_value
    self.num_bins=num_bins
    self.step_size=step_size
    self.bin_size=(self.max_value-self.min_value)/float(num_bins)
    self.periodicity=periodicity
  
class PMF():
  def __repr__(self):
    return "PMF({0},{1})".format(self.points,self.values)

  def __init__(self,points,values):
    self.points=points
    self.values=values
    self.interpolator=scipy.interpolate.interpnd.LinearNDInterpolator(self.points,self.values)
    #self.interpolator=NearestNDInterpolator(self.points,self.values)
  def GetValue(self,point):
    return self.interpolator(point)


