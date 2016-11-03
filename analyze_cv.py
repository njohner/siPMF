"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This file contains the functions to analyse the collective variables in the
different :class:`Window` of the :class:`System`
"""
import os,subprocess
import numpy as npy
import matplotlib.pyplot as plt

def _OutputDir(system):
  outdir=os.path.join(system.basedir,"cv_analysis")
  if not os.path.isdir(outdir):subprocess.call(["mkdir",outdir])
  return outdir

def PlotCollectiveVariables(system,update_datafiles=True,markersize=1.0,linewidth=1.0):
  outdir=_OutputDir(system)
  if update_datafiles:system.UpdateDataFiles(new_only=False)
  for w in system.windows:
    t,cvs=w.ReadDataFile()
    for i,cv in enumerate(cvs):
      plt.figure()
      plt.plot(t,cv,'.',markersize=markersize)
      plt.hlines(w.cv_values[i],t[0],t[-1],linestyle="--",color="r",linewidth=linewidth,zorder=4)
      plt.xlim([t[0],t[-1]])
      plt.title("{0} for window {1}".format(system.cv_list[i].name,w.name))
      plt.xlabel("Time")
      plt.ylabel(system.cv_list[i].name+" "+system.cv_list[i].units)
      plt.savefig(os.path.join(outdir,"{0}_{1}.png".format(w.name,system.cv_list[i].name)))
      plt.close()
  return

def PlotAutocorrelations(system,shifts=[],update_datafiles=True,all_windows=True):
  outdir=_OutputDir(system)
  if update_datafiles:system.UpdateDataFiles(new_only=False)
  for w in system.windows:
    if (not all_windows) and hasattr(w,"auto_correlation_times"):continue
    t,cvs=w.ReadDataFile()
    ndata=len(t)
    if not shifts:
      dts=[]
      for i in range(ndata):
        dt=2**i
        if dt<ndata:dts.append(dt)
        else:break
    else:
      dts=[el for el in shifts if el<ndata]
    auto_corr_times=[]
    for i,cv in enumerate(cvs):
      cv=npy.array(cv)
      cl=[]
      for dt in dts:cl.append(npy.corrcoef(cv[:-dt],cv[dt:])[0,1])
      bools=npy.array(cl)<0.1
      auto_corr_time=dts[-1]
      for j in range(len(dts)-1,-1,-1):
        if not bools[j]:break
        auto_corr_time=dts[j]
      plt.figure()
      plt.plot(dts,cl,'-o')
      plt.xlim([1,dts[-1]+1])
      plt.xscale("log")
      plt.xlabel("Time shift")
      plt.ylabel("Autocorrelation")   
      plt.title("Autocorr for {0} for window {1}:{2}".format(system.cv_list[i].name,w.name,auto_corr_time))
      plt.savefig(os.path.join(outdir,"{0}_{1}_autocorr.png".format(w.name,system.cv_list[i].name)))
      plt.close()
      auto_corr_times.append(auto_corr_time)
    w.auto_correlation_times=auto_corr_times
  cv1=[]
  cv2=[]
  auto_corr=[]
  for w in system.windows:
    cv1.append(w.cv_values[0])
    cv2.append(w.cv_values[1])
    auto_corr.append(npy.log10(npy.max(w.auto_correlation_times)))
  plt.figure()
  plt.tricontourf(cv1,cv2,auto_corr, 20)
  plt.xlabel(system.cv_list[0].name+" "+system.cv_list[0].units)
  plt.ylabel(system.cv_list[1].name+" "+system.cv_list[1].units)
  cbar=plt.colorbar()
  cbar.set_label("log10 decorrelation time", rotation=270)
  plt.title("Decorrelation times")
  plt.savefig(os.path.join(outdir,"decorrelation_times.png"))
  plt.close()
  return

