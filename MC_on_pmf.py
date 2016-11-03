"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains the analysis function for trajectories genereated with a siPMF run.
"""
import os,subprocess,logging
import numpy as npy
import matplotlib.pyplot as plt
#from matplotlib.mlab import griddata
from ost import *
import random
import itertools
from scipy.interpolate import griddata

__all__=('GenerateTrajectory','PlotTrajOnPMF')

kB=0.0019872041#kcal/mol

#def CalculateDiffusionConstant(system):
class Node():
  def __repr__(self):
    return "Node({0},{1},{2})".format(self.cv_values,self.energy,self.diffusion_constants)

  def __init__(self,cv_values,energy,diffusion_constants):
    self.cv_values=cv_values
    self.bonds=[]
    self.rates=[]
    self.energy=energy
    self.diffusion_constants=diffusion_constants
  
  def AddBond(self,bond,rate):
    self.bonds.append(bond)
    self.rates.append(rate)

def MakeGrid(system,num_bins=None):
  if not num_bins:num_bins=[cv.num_bins for cv in system.cv_list]
  xl=[]
  for cv,nb in zip(system.cv_list,num_bins):
    xl.append(npy.linspace(cv.min_value,cv.max_value,nb))
  step_sizes=npy.array([x[1]-x[0] for x in xl])
  grid=list(itertools.product(*xl))
  return grid,step_sizes

class Grid():
  def __init__(self,grid,values):
    self.grid=grid
    self.values=values
  def GetValue(point):
    return self.values(self.grid.index(point))

def PMFOnGrid(pmf,grid):
  return Grid(grid,griddata(pmf.points,pmf.values,grid))

def DiffusionConstantsOnGrid(system,grid):
  cv_vals=[w.cv_values for w in system.windows]
  diff_grids=[]
  for i in range(system.dimensionality):
    diff_coeffs=npy.array([w.diffusion_constants[i] for w in system.windows])
    diff_points=npy.array([w.cv_values for w in system.windows])
    g1=griddata(diff_points,diff_coeffs,grid)
    g2=griddata(diff_points,diff_coeffs,grid,method='nearest')
    g1[npy.where(npy.isnan(g1))]=g2[npy.where(npy.isnan(g1))]
    diff_grids.append(Grid(grid,g1))
  return diff_grids

class Bond():
  def __repr__(self):
    return "Bond({0},{1},{2},{3},{4},{5})".format(self.first,self.second,self.length,self.D,self.kf,self.kb)
  def __init__(self,node1,node2,length,D,kf,kb):
    self.first=node1
    self.second=node2
    self.length=length
    self.D=D
    self.kf=kf
    self.kb=kb
  def Other(self,node):
    if self.first==node:return self.second
    elif self.second==node:return self.first
    else:raise(IOError("could not find node {0} in bond {1}".format(node,self)))

class Network():
  def __init__(self):
    self.nodes=[]
    self.bonds=[]
  
  def FindNode(self,cv_values):
    for n in self.nodes:
      if tuple(n.cv_values)==tuple(cv_values):return n
    return None

  def FindClosestNode(self,cv_values):
    n=self.nodes[0]
    cv_values=npy.array(cv_values)
    v=(npy.array(n.cv_values)-cv_values)
    d=npy.dot(v,v)
    for ni in self.nodes:
      v=(npy.array(ni.cv_values)-cv_values)
      di=npy.dot(v,v)
      if di<d:
        n=ni
        d=di
    return n

  def InitializeFromSystem(self,system,num_bins=None,max_E=None):
    if not num_bins:num_bins=[cv.num_bins for cv in system.cv_list]
    if not max_E:max_E=system.max_E_plot
    grid,step_sizes=MakeGrid(system,num_bins)
    pmf_grid=PMFOnGrid(system.pmf,grid)
    D_grids=DiffusionConstantsOnGrid(system,grid)
    for i in range(len(grid)):
      if pmf_grid.values[i]<max_E:self.AddNode(grid[i],pmf_grid.values[i],[D.values[i] for D in D_grids])
    for i,node1 in enumerate(self.nodes):
      for node2 in self.nodes[i+1:]:
        v=npy.abs(npy.array(node1.cv_values)-npy.array(node2.cv_values))
        if all(v<=1.1*step_sizes) and npy.sum(v<=0.1*step_sizes)==1:self.LinkNodes(node1,node2,system.temperature)

  def BrownianDynamics(self,init_pos,nsteps=1000,final_nodes=[]):
    from numpy.random import choice
    n0=self.FindNode(init_pos)
    if not n0:n0=self.FindClosestNode(init_pos)
    n=n0
    nl=[n0]
    tl=[0.0]
    for i in range(nsteps):
      lambd=npy.sum(n.rates)
      holding_time=random.expovariate(lambd)
      tl.append(tl[-1]+holding_time)
      bi=choice(len(n.bonds),p=n.rates/lambd)
      n=n.bonds[bi].Other(n)
      nl.append(n)
      if n in final_nodes:break
    return  tl,nl
  """
  def InitializeFromPMF(self,pmf,step_sizes,max_E):
    step_sizes=npy.array(step_sizes)
    l=[]
    for cv,step in zip(pmf.cv_list,step_sizes):
      l.append(npy.arange(cv.min_value,cv.max_value+step/2.,step))
    grid_iter=itertools.product(*l)
    for el in grid_iter:
      if pmf.GetValue(el)<max_E:self.AddNode(el)
    for i,node1 in enumerate(self.nodes):
      for node2 in self.nodes[i:]:
        v=npy.abs(npy.array(node1.cv_values)-npy.array(node2.cv_values))
        if all(v<step_sizes):self.LinkNodes(node1,node2)
  """
  def AddNode(self,cv_values,energy,diffusion_constants):
    self.nodes.append(Node(cv_values,energy,diffusion_constants))
  
  def LinkNodes(self,node1,node2,temperature):
    v=npy.array(node1.cv_values)-npy.array(node2.cv_values)
    d=npy.dot(v,v)**0.5
    c=npy.abs(v)/(2.*d)
    D=npy.dot(c,node1.diffusion_constants)+npy.dot(c,node2.diffusion_constants)
    kf=D/(d*d)*npy.exp(-(node2.energy-node1.energy)/(kB*temperature))
    kb=D/(d*d)*npy.exp(-(node1.energy-node2.energy)/(kB*temperature))
    bond=Bond(node1,node2,d,D,kf,kb)
    self.bonds.append(bond)
    node1.AddBond(bond,kf)
    node2.AddBond(bond,kb)


def GenerateStateNetworkFromPMF(system,step_sizes,max_E=None):
  if not len(step_sizes)==system.dimensionality:
    raise(IOError("step_sizes should contain one step size for each dimensionality of the system"))
  if not max_E:max_E=system.max_E_plot
  network=Network()
  network.InitializeFromPMF(system.pmf,step_sizes,max_E)


def _brownian_dynamics(div_pot,x0=0.0,T=300.0,D=1.0,time_step=0.01,nsteps=1000):
  """
  div_pot should be in units of kbT/distance
  """
  import random
  mu=0.0
  sigma=1.0
  xl=[]
  tl=[]
  x=x0
  for i in range(nsteps):
    #v=-div_pot(x)/gamma+(2*D)**0.5*random.gauss(mu,sigma)/(time_step**0.5)
    v=-div_pot(x)*D+(2*D)**0.5*random.gauss(mu,sigma)/(time_step**0.5)
    x=x+v*time_step
    xl.append(x)
    tl.append(i*time_step)
  return npy.array(tl),npy.array(xl)


def _test_laplace_transform():
  a=10
  w=10
  t=npy.arange(0,10,0.001)
  y=npy.exp(-a*t)*npy.cos(w*t)
  x,yl=LaplaceTransform(t,y)
  xl=(x+a)/((x+a)**2.0+w**2.0)
  plt.figure()
  plt.plot(x,yl)
  plt.plot(x[::10],xl[::10],'x')
  plt.savefig("test.png")
  plt.close()

def autocorrelation(x):
  x=npy.array(x)
  sm=npy.mean(x*x)
  #ac=[npy.corrcoef(x[i:],x[:-i])[0,1] for i in range(1,len(x)/2)]
  ac=[1.0]
  ac.extend([npy.mean(x[i:]*x[:-i])/sm for i in range(1,int(0.9*len(x)))])
  return npy.array(ac)

def TrapezoidalIntegration(x,y):
  x=npy.array(x)
  y=npy.array(y)
  y=(y[:-1]+y[1:])/2.
  dx=x[1:]-x[:-1]
  return npy.sum(dx*y)

def LaplaceTransform(x,y,s_min=0,s_max=5,s_step=0.01):
  x=npy.array(x)
  y=npy.array(y)
  s=npy.arange(s_min,s_max,s_step)
  yl=[]
  for si in s:
    z=npy.exp(-si*x)*y
    yl.append(TrapezoidalIntegration(x,z))
  return (s,npy.array(yl))


def GenerateTrajectory(system,init_cv,step_sizes,temperature=300,max_nstep=1000):
  current_cv=npy.array(init_cv)
  current_E=system.pmf.GetValue(current_cv)
  traj=[current_cv]
  kb=1.98*0.001 #in kCal/K/mol
  for i in range(max_nstep):
    step=npy.array([(random.random()-0.5)*ss for ss in step_sizes])
    E=system.pmf.GetValue(current_cv+step)
    if random.random()<npy.exp(-(E-current_E)/(kb*temperature)):
      current_cv=current_cv+step
      current_E=system.pmf.GetValue(current_cv)
      traj.append(current_cv)
  return traj

def PlotTrajOnPMF(system,traj,outputdir,filename,skip=100,n_levels=None,max_E=None):
    energy_units="kcal/mol"
    if not max_E:max_E=system.max_E_plot
    if not n_levels:n_levels=int(max_E)
    if system.pmf.dimensionality==2:
      X=system.pmf.points[:,0]
      Y=system.pmf.points[:,1]
      Z=system.pmf.values
      if max_E:Z=npy.array([min(el,max_E) for el in Z])
      num_pads=max([system.pmf.cv_list[0].num_pads,system.pmf.cv_list[1].num_pads])
      nb_x=system.pmf.cv_list[0].wham_num_bins+2*num_pads
      nb_y=system.pmf.cv_list[1].wham_num_bins+2*num_pads
      X=X.reshape([nb_x,nb_y])
      Y=Y.reshape([nb_x,nb_y])
      Z=Z.reshape([nb_x,nb_y])
      plt.figure()
      plt.contourf(X,Y,Z,n_levels)
      plt.colorbar(label=energy_units)
      plt.contour(X,Y,Z,n_levels,colors="k")
      if system.pmf.cv_list[0].units:plt.xlabel("{0} [{1}]".format(system.pmf.cv_list[0].name,system.pmf.cv_list[0].units))
      else:plt.xlabel("{0}".format(system.pmf.cv_list[0].name))
      if system.pmf.cv_list[1].units:plt.ylabel("{0} [{1}]".format(system.pmf.cv_list[1].name,system.pmf.cv_list[1].units))
      else:plt.ylabel("{0}".format(system.pmf.cv_list[1].name))
      plt.plot([el.cv_values[0] for el in traj[::skip]],[el.cv_values[1] for el in traj[::skip]],color="m")
      plt.savefig(os.path.join(outputdir,filename+".png"))
      plt.close()
    elif system.pmf.dimensionality==1:
      plt.figure()
      plt.plot(system.pmf.points,system.pmf.values)
      if max_E:plt.ylim([0,max_E])
      if system.pmf.cv_list[0].units:plt.xlabel("{0} [{1}]".format(system.pmf.cv_list[0].name,system.pmf.cv_list[0].units))
      else:plt.xlabel("{0}".format(system.pmf.cv_list[0].name))
      plt.ylabel("Free Energy")
      plt.savefig(os.path.join(outputdir,filename+".png"))
      plt.close()


class DiffusionConstants():
  def __init__(self,points,values):
    self.points=points
    self.values=values

"""
tau=1.0
xl=[]
for i in range(10000):
  xl.append(random.expovariate(tau))
print tau,npy.mean(xl)
y,x,h=plt.hist(xl,bins=100)
plt.close()
y/=y[0]
x=0.5*(x[1:]+x[:-1])
x2=npy.arange(0,x[-1],0.1)
plt.figure()
plt.plot(x,y)
plt.plot(x2,npy.exp(-x2*tau),linestyle='--',color='r')
plt.show()
"""

