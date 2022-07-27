#!/usr/bin/env python
import sys, os
#from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Times new Roman'
plt.rcParams['font.size'] = 16
#plt.rcParams["font.weight"] = 60
#plt.rcParams["axes.labelweight"] = 60
plt.rcParams["axes.labelsize"] = 16
#plt.rcParams["axes.ticksize"] = 16
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.width'] = 2
import numpy as np
from scipy.optimize import curve_fit, leastsq
from copy import deepcopy
import pickle as pkl


# Band with errors
dataplus = 'pkl_files/data_0.1.pkl'
datafile = 'pkl_files/data_0.0.pkl'
dataminus = 'pkl_files/data_-0.1.pkl'
alpha=1

if __name__ == "__main__":
    data_in = pkl.load(open(dataplus,'rb'),encoding='latin1')
    data_plus={}
    pots_plus,phs=[],[]
    nsteps=2
    for dat in data_in['rate_map']:
        pot_plus,ph=np.around(dat[0][0],3),np.around(dat[0][1],3)
        if pot_plus not in data_plus:
            data_plus[pot_plus] = {}
        data_plus[pot_plus][ph] = dat[1]
        if pot_plus not in pots_plus:
            pots_plus.append(pot_plus)
        if ph not in phs:
            phs.append(ph)

    #data_in = pkl.load(open('data_py3.pkl','rb'),encoding='latin1')
    data_in = pkl.load(open(dataminus,'rb'),encoding='latin1')
    data_minus={}
    pots_minus,phs=[],[]
    nsteps=2
    for dat in data_in['rate_map']:
        pot_minus,ph=np.around(dat[0][0],3),np.around(dat[0][1],3)
        if pot_minus not in data_minus:
            data_minus[pot_minus] = {}
        data_minus[pot_minus][ph] = dat[1]
        if pot_minus not in pots_minus:
            pots_minus.append(pot_minus)
        if ph not in phs:
            phs.append(ph)

    #data_in = pkl.load(open('data_py3.pkl','rb'),encoding='latin1')
    data_in = pkl.load(open(datafile,'rb'),encoding='latin1')
    data={}
    pots,phs=[],[]
    nsteps=2
    for dat in data_in['rate_map']:
        pot,ph=np.around(dat[0][0],3),np.around(dat[0][1],3)
        if pot not in data:
            data[pot] = {}
        data[pot][ph] = dat[1]
        if pot not in pots:
            pots.append(pot)
        if ph not in phs:
            phs.append(ph)

    print(data.keys())
    print(data_plus.keys())
    print(data_minus.keys())
    print(len(data.keys()))
    print(len(data_plus.keys()))
    print(len(data_minus.keys()))

    X=np.array(sorted(pots))
    X_plus=np.array(sorted(pots_plus))
    X_minus=np.array(sorted(pots_minus))
    scale='RHE'
    RHE_pots=X+0.4
    Y=np.array(sorted(phs))

    fig,ax=plt.subplots(nsteps,2,sharex=True)

    for col in range(nsteps):
     for istep in range(2):
        XY=np.ones((len(X),len(Y)))*0.5
        rate=np.ones((len(X),len(Y)))*0.5
        for ix,x in enumerate(X):
            for iy,y in enumerate(Y):
                            try:
                                #XY[ix][iy]=data[x][y][istep]/np.sum(data[x][y][:3])
                                if col == 1:
                                    XY[ix][iy]=data[x][y][istep]/np.sum(data[x][y][:nsteps])
                                else:
                                    rate[ix][iy]=data[x][y][istep]#/np.sum(data[x][y][:nsteps])
                            except Exception as e:
                                XY[ix][iy]=np.nan#data[x][y][istep]#/np.sum(data[x][y][:3])
                                rate[ix][iy]=1e-20#np.nan#data[x][y][istep]#/np.sum(data[x][y][:nsteps])
                                print('Here')
                                print(repr(e))
                                print(x,y)
                                pass

        XY_plus=np.ones((len(X_plus),len(Y)))*0.5
        rate_plus=np.ones((len(X_plus),len(Y)))*0.5
        for ix,x in enumerate(X_plus):
            for iy,y in enumerate(Y):
                            try:
                                #XY_plus[ix][iy]=data_plus[x][y][istep]/np.sum(data_plus[x][y][:3])
                                if col == 1:
                                    XY_plus[ix][iy]=data_plus[x][y][istep]/np.sum(data_plus[x][y][:nsteps])
                                else:
                                    rate_plus[ix][iy]=data_plus[x][y][istep]#/np.sum(data_plus[x][y][:nsteps])
                            except Exception as e:
                                XY_plus[ix][iy]=np.nan#data_plus[x][y][istep]#/np.sum(data_plus[x][y][:3])
                                rate_plus[ix][iy]=1e-20#np.nan#data_plus[x][y][istep]#/np.sum(data_plus[x][y][:nsteps])
                                print('Here plus')
                                print(repr(e))
                                print(x,y)
                                pass

        XY_minus=np.ones((len(X_minus),len(Y)))*0.5
        rate_minus=np.ones((len(X_minus),len(Y)))*0.5
        for ix,x in enumerate(X_minus):
            for iy,y in enumerate(Y):
                            try:
                                #XY_minus[ix][iy]=data_minus[x][y][istep]/np.sum(data_minus[x][y][:3])
                                if col == 1:
                                    XY_minus[ix][iy]=data_minus[x][y][istep]/np.sum(data_minus[x][y][:nsteps])
                                    #print(XY_minus)
                                else:
                                    rate_minus[ix][iy]=data_minus[x][y][istep]#/np.sum(data_minus[x][y][:nsteps])
                            except Exception as e:
                                XY_minus[ix][iy]=np.nan#data_minus[x][y][istep]#/np.sum(data_minus[x][y][:3])
                                rate_minus[ix][iy]=1e-20#np.nan#data_minus[x][y][istep]#/np.sum(data_minus[x][y][:nsteps])
                                print('Here minus')
                                print(repr(e))
                                print(x,y)
                                pass

        if col == 0:
          if len(phs) == 1:
              #print(XY)
              print(f'selectivity {len(XY)}')
              #print(len(XY))
              print(f'Rate plus {len(XY_plus)}')
              #print(XY_plus[0])
              print(f'Rate minus {len(XY_minus)}')
              #print(XY_minus)
              b = ax[istep][0].plot(sorted(RHE_pots),rate,'-k',linewidth=2)
              flat_rate_plus = [item for sublist in rate_plus for item in sublist]
              flat_rate_minus = [item for sublist in rate_minus for item in sublist]
              ax[istep][0].fill_between(sorted(RHE_pots),flat_rate_plus,flat_rate_minus,alpha=alpha)#,color)
              ax[istep][0].set_yscale('log')
              ax[istep][0].set_ylim([1e-4,1e4])
          else:
           b = ax[istep][0].imshow(rate.T,
                interpolation='bicubic',
                cmap=cm.jet,
                   origin='lower', extent=[X.min(), X.max(), Y.min(), Y.max()],norm=LogNorm(),#,
                    vmin=1e-15,
                    vmax=1e5,#)
                    aspect='auto')#, vmin=-abs(alldata[:,2]).max())

        else:
          if len(phs) == 1:
              flat_XY_plus = [item for sublist in XY_plus for item in sublist]
              flat_XY_minus = [item for sublist in XY_minus for item in sublist]
              #print(len(flat_XY_plus))
              #das
              a = ax[istep][1].plot(sorted(RHE_pots),XY,'-k',linewidth=2)
              cut=-1
              ax[istep][1].fill_between(sorted(RHE_pots[:cut]),flat_XY_plus[:cut],flat_XY_minus[:cut],alpha=alpha)#,color)
          else:
           a = ax[istep][1].imshow(XY.T,
                interpolation='bicubic',
                cmap=cm.RdYlGn,
                   origin='lower', extent=[X.min(), X.max(), Y.min(), Y.max()],#norm=LogNorm(),#,
                    vmin=0,
                    vmax=1,
                    aspect='auto')#, vmin=-abs(alldata[:,2]).max())

        if istep == nsteps-1:
#            fig.colorbar(a,ax=ax[0][1],fraction=0.146, pad=.04,label='Selectivity')
            for thisax in ax[istep]:
                thisax.set_xlabel('U$_{\mathrm{RHE}}$ [V]')
                #ax[istep][0].set_xlabel('U$_{SHE}$ [V]')
        #else:
            #for thisax in ax[istep]:
            #    thisax.set_xticks([])
        ax[istep][0].set_ylabel('TOF / s$^{-1}$')
        ax[istep][1].set_ylim([-0.1,1.1])
        ax[istep][1].set_ylabel('Selectivity')
        ax[istep][1].yaxis.set_ticks_position('right')
        ax[istep][1].yaxis.set_label_position('right')
        #eli:

#            fig.colorbar(a,ax=ax[:,0],fraction=0.046, pad=0.04,label='Turnover frequency [s$^{-1}$]')

#     fig.subplots_adjust(right=0.8)
#     cbar_ax = fig.add_axes([0.45, 0.15, 0.05, 0.7])
    if len(phs) > 1:
     axins1 = inset_axes(ax[0][0],
                    width="100%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    loc='lower left',
                     bbox_to_anchor=(0, 1.3, 1, 1),
                  bbox_transform=ax[0][0].transAxes)
    #axins1.xaxis.set_ticks_position('top')
    #axins1.yaxis.set_ticks_position('left')
     fig.colorbar(b, cax=axins1,orientation='horizontal',fraction=0.07,anchor=(1.0,1.0))#cbar_ax)

     axins2 = inset_axes(ax[0][1],
                    width="100%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    loc='lower left',
                    bbox_to_anchor=(0, 1.3, 1, 1),
                    bbox_transform=ax[0][1].transAxes)
     fig.colorbar(a, cax=axins2,orientation='horizontal',fraction=0.07,anchor=(1.0,1.0))#cbar_ax)

     ax[0][0].annotate('Turnover frequency [s$^{-1}$]',(-0.8,14),color='k',bbox=dict(facecolor='w', edgecolor=None))
     ax[0][1].annotate('Selectivity',(-1.0,14),color='k',bbox=dict(facecolor='w', edgecolor=None))
     fig.subplots_adjust(wspace=0.2,hspace=0.01)
    else:
     ax[1][1].text(-0.725,0.94,'HCOO$^-$',color='k',bbox=dict(facecolor='w', edgecolor='w'),ha='center',va='bottom',fontsize=20)
     ax[0][1].text(-0.725,0.94,'CO',color='k',bbox=dict(facecolor='w', edgecolor='w'),ha='center',va='bottom',fontsize=20)
     for row in ax:
         for col in row:
             col.set_xlim([-0.7,-0.2])
     fig.subplots_adjust(wspace=0.1,hspace=0.1)

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.1,hspace=0.1)
    plt.savefig('Fig6b.pdf')
    plt.show()
