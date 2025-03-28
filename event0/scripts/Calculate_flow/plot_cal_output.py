#!/usr/bin/env python3

from numpy import *
from sys import argv, exit
from os import path, mkdir
import h5py
from scipy import interpolate
import shutil
import pylab as pl
import numpy as np

aa = np.loadtxt("QAQBall")

MAMB_0p0_3p0 = aa[:,0]
MAMB_0p3_3p0 = aa[:,1]
MAMB_0p5_3p0 = aa[:,2]
meanpTch = aa[:,3]
meanpTch_nevent = aa[:,4]

QAQB1_0p0_3p0 = aa[:,5]
QAQB1_0p3_3p0 = aa[:,6]
QAQB1_0p5_3p0 = aa[:,7]
meanpTpi = aa[:,8]
meanpTpi_nevent = aa[:,9]

QAQB2_0p0_3p0 = aa[:,10]
QAQB2_0p3_3p0 = aa[:,11]
QAQB2_0p5_3p0 = aa[:,12]
meanpTka = aa[:,13]
meanpTka_nevent = aa[:,14]

QAQB3_0p0_3p0 = aa[:,15]
QAQB3_0p3_3p0 = aa[:,16]
QAQB3_0p5_3p0 = aa[:,17]
meanpTpr = aa[:,18]
meanpTpr_nevent = aa[:,19]


Nchbins = np.array([5, 15,25,35,45,55,65,75,85,95,105])

#pl.yscale('log')
Nchall = np.loadtxt("Nchall")
Nchbins = Nchall[:,1]/Nchall[:,0]
exec("pl.figure({})".format(0))
data = np.loadtxt("v2_033.txt")
Nchdata = data[:,0]
data = np.loadtxt("v2_003.txt")
pl.plot(Nchdata, data[:,1],'ro',label='CMS data')
pl.errorbar(Nchdata, data[:,1],data[:,2], capsize=0, ls='none', color='r', elinewidth=2)
pl.plot(Nchbins, (QAQB2_0p0_3p0/MAMB_0p0_3p0)**0.5,'r',label='$0.0 < p^{j}_{T} < 3.0 GeV/c$')
data = np.loadtxt("v2_033.txt")
pl.plot(Nchdata, data[:,1],'bh',label='CMS data')
pl.errorbar(Nchdata, data[:,1],data[:,2], capsize=0, ls='none', color='b', elinewidth=2)
pl.plot(Nchbins, (QAQB2_0p3_3p0/MAMB_0p3_3p0)**0.5,'b',label='$0.3 < p^{j}_{T} < 3.0 GeV/c$')
data = np.loadtxt("v2_053.txt")
pl.plot(Nchdata, data[:,1],'g*',label='CMS data')
pl.errorbar(Nchdata, data[:,1],data[:,2], capsize=0, ls='none', color='g', elinewidth=2)
pl.plot(Nchbins, (QAQB2_0p5_3p0/MAMB_0p5_3p0)**0.5,'g',label='$0.5 < p^{j}_{T} < 3.0 GeV/c$')
pl.legend( loc='upper right' , fontsize=12)
#pl.xlim(0,0.25)
#pl.ylim(0.0,0.32)
#pl.title('Incoherent, e+U, '.format(plotname[ip]), fontsize=15)# give plot a title
pl.xlabel('$Nch^{j}$', fontsize=15)# make axis labels
pl.ylabel('$v_{2}\{2, |\Delta\eta|>2\}$', fontsize=12)
pl.xticks(fontsize=15)
pl.yticks(fontsize=13)
#pl.show()# show the plot on the screen
pl.savefig('v2.png', dpi=120)
    
exec("pl.figure({})".format(30))
pl.plot(Nchbins, meanpTpi/meanpTpi_nevent,'r',label='$\pi$')
pl.plot(Nchbins, meanpTka/meanpTka_nevent,'b',label='$K$')
pl.plot(Nchbins, meanpTpr/meanpTpr_nevent,'g',label='$P$')
pl.legend( loc='upper right' , fontsize=12)
#pl.xlim(0,0.25)
#pl.ylim(0.0,0.32)
#pl.title('Incoherent, e+U, '.format(plotname[ip]), fontsize=15)# give plot a title
pl.xlabel('$Nch^{j}$', fontsize=15)# make axis labels
pl.ylabel('$<p_{T}>$', fontsize=12)
pl.xticks(fontsize=15)
pl.yticks(fontsize=13)
#pl.show()# show the plot on the screen
pl.savefig('mean_pT.png', dpi=120)

