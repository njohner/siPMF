"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains the analysis function for trajectories genereated with a siPMF run.
"""
import os,subprocess,logging
import numpy as npy
import matplotlib.pyplot as plt
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

