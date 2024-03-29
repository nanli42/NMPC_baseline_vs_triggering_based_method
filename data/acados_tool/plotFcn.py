#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

# author: Daniel Kloeser

from acados_tool.readDataFcn import getTrack
from acados_tool.time2spatial import transformProj2Orig,transformOrig2Proj
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

def plotTrackProj(simX, opt=1, T_opt=None, vel=True, whole=True, show=True, traj_color='b-', traj_name=""):
    if (opt==1):
        filename = 'LMS_Track.txt'
    if (opt==2):
        filename = 'bigger_track.txt'

    # scale
    gain = 1

    # load track
    s=simX[:,0]
    n=simX[:,1]
    alpha=simX[:,2]
    v=simX[:,3]
    distance=0.17 * gain
    # transform data
    [x, y, _, _] = transformProj2Orig(s, n, alpha, v, filename)

    # plot racetrack map
    #Setup plot
    if (show):
        if (opt==1):
            plt.figure(figsize=(10, 10), dpi=80)
            plt.ylim(bottom=-1.75*gain,top=0.35*gain)
            plt.xlim(left=-1.1*gain,right=1.6*gain)
        if (opt==2):
            plt.figure(figsize=(10, 10), dpi=80)
            plt.ylim(bottom=-2.0*gain,top=2.0*gain)
            plt.xlim(left=-1.5*gain,right=2.0*gain)
        plt.ylabel('y[m]')
        plt.xlabel('x[m]')

        # Plot center line
        [Sref,Xref,Yref,Psiref,_] = getTrack(filename)
        plt.plot(Xref,Yref,'--',color='k')

        # Draw Trackboundaries
        Xboundleft=Xref-distance*np.sin(Psiref)
        Yboundleft=Yref+distance*np.cos(Psiref)
        Xboundright=Xref+distance*np.sin(Psiref)
        Yboundright=Yref-distance*np.cos(Psiref)
        plt.plot(Xboundleft,Yboundleft,color='k',linewidth=1)
        plt.plot(Xboundright,Yboundright,color='k',linewidth=1)

        ax = plt.gca()
        ax.set_aspect('equal', 'box')

        # Put markers for s values
        if (opt==1):
            xi=np.zeros(9)
            yi=np.zeros(9)
            xi1=np.zeros(9)
            yi1=np.zeros(9)
            xi2=np.zeros(9)
            yi2=np.zeros(9)
            for ii in range(0,8*gain,gain):
                i = int(ii / gain)
                try:
                    k = list(Sref).index(ii + min(abs(Sref - ii)))
                except:
                    k = list(Sref).index(ii - min(abs(Sref - ii)))
                [_,nrefi,_,_]=transformOrig2Proj(Xref[k],Yref[k],Psiref[k],0)
                [xi[i],yi[i],_,_]=transformProj2Orig(Sref[k],nrefi+0.26*gain,0,0)
                plt.text(xi[i], yi[i], '{}m'.format(ii), fontsize=8,horizontalalignment='center',verticalalignment='center')
                [xi1[i],yi1[i],_,_]=transformProj2Orig(Sref[k],nrefi+0.12*gain,0,0)
                [xi2[i],yi2[i],_,_]=transformProj2Orig(Sref[k],nrefi+0.15*gain,0,0)
                plt.plot([xi1[i],xi2[i]],[yi1[i],yi2[i]],color='black')
        if (opt==2):
            N = 18
            xi=np.zeros(N)
            yi=np.zeros(N)
            xi1=np.zeros(N)
            yi1=np.zeros(N)
            xi2=np.zeros(N)
            yi2=np.zeros(N)
            for ii in range(0,N*gain,gain):
                i = int(ii / gain)
                try:
                    k = list(Sref).index(ii + min(abs(Sref - ii)))
                except:
                    k = list(Sref).index(ii - min(abs(Sref - ii)))
                [_,nrefi,_,_]=transformOrig2Proj(Xref[k],Yref[k],Psiref[k],0,filename)
                [xi[i],yi[i],_,_]=transformProj2Orig(Sref[k],nrefi+0.14*gain,0,0,filename)
                plt.text(xi[i], yi[i], '{}m'.format(ii), fontsize=8,horizontalalignment='center',verticalalignment='center')
                [xi1[i],yi1[i],_,_]=transformProj2Orig(Sref[k],nrefi+0.09*gain,0,0,filename)
                [xi2[i],yi2[i],_,_]=transformProj2Orig(Sref[k],nrefi+0.11*gain,0,0,filename)
                plt.plot([xi1[i],xi2[i]],[yi1[i],yi2[i]],color='black')

    # Draw trajectory
    if (whole):
        plt.plot(x, y, traj_color, label=traj_name)

    if (vel):
        # Draw driven trajectory
        heatmap = plt.scatter(x,y, c=v, cmap=cm.rainbow, edgecolor='none', marker='o')
        cbar = plt.colorbar(heatmap, fraction=0.035)
        cbar.set_label("velocity in [m/s]")

def plotRes(simX, str):
    # plot results
    plt.figure()
    t = np.linspace(0.0, simX.shape[0], simX.shape[0])
    plt.plot(t, simX[:,:])
    plt.ylabel('x')
    plt.xlabel('t')
    plt.legend(str)
    plt.grid(True)
