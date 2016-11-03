"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains the analysis function for trajectories genereated with a siPMF run.
"""
import os,subprocess,logging
import numpy as npy
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from ost import *

__all__=('BuildTrajectory')

def BuildTrajectory(system,windows,eh,traj_filename,stride=1):
  """
  Makes a trajectory from the simulations of a list of windows

  :param system: The system for which to generate a trajectory
  :param windows: The list of windows from which to generate the trajectory
  :param eh: The entity
  :param traj_filename: The name of the trajectory files

  :type system:  :class:`~system.System`
  :type windows:  :class:`list` (:class:`~window.Window`)
  :type eh: :class:`mol.EntityHandle`
  :type traj_filename: :class:`str`
  """
  t_out=mol.CreateCoordGroup(eh.atoms)
  for w in windows:
    p=w.phases[-1]
    t=io.LoadCHARMMTraj(eh,os.path.join(p.outdir,traj_filename),stride=stride)
    print os.path.join(p.outdir,traj_filename),t.GetFrameCount()
    vl=t.GetFramePositions(t.GetFrameCount()-1)
    t_out.AddFrame(vl)
  return t_out


def AnalyzeWindows(system,eh,analysis_function,function_arguments,traj_filename,first_phase=0,first_frame=0,stride=1,windows=None):
  if not windows:windows=system.windows
  res_dict={}
  for w in windows:
    res=[]
    for p in w.phases[first_phase:]:
      t=io.LoadCHARMMTraj(eh,os.path.join(p.outdir,traj_filename),stride=stride)
      res.extend(analysis_function(t,*function_arguments)[first_frame:])
    res_dict[w]=res
  return res_dict

def PlotWindowAverage(system,res_dict,title=None,outname=None):
  dimensionality=system.dimensionality
  if dimensionality>2:
    print "can only plot averages for systems with dimensionality 1 or 2"
    return
  data=[]
  for w in res_dict:
    m=npy.average(res_dict[w])
    if npy.isnan(m):continue
    res=list(w.cv_values)
    res.append(m)
    data.append(res)
  data=npy.array(data)
  f=plt.figure()
  if dimensionality==1:
    plt.plot(data[:,0],data[:,1],'o')
    plt.xlabel("{0} [{1}]".format(system.cv_list[0].name,system.cv_list[0].units))
    plt.ylabel("Average")
  elif dimensionality==2:
    xl=[]
    for cv in system.cv_list:
      xl.append(npy.linspace(cv.min_value,cv.max_value,cv.num_bins))
    xi,yi=npy.meshgrid(*xl)
    zi=griddata(data[:,0],data[:,1],data[:,2], xi, yi)
    plt.contourf(xi, yi, zi)
    plt.colorbar()
    plt.xlabel("{0} [{1}]".format(system.cv_list[0].name,system.cv_list[0].units))
    plt.ylabel("{0} [{1}]".format(system.cv_list[1].name,system.cv_list[1].units))
    cv=system.cv_list[0]
    plt.xlim(cv.wham_min_value,cv.wham_max_value)
    cv=system.cv_list[1]
    plt.ylim(cv.wham_min_value,cv.wham_max_value)
  if title:plt.title(title)
  if outname:plt.savefig(outname)
  return f


