
# -*- coding: utf-8 -*-
"""
updated stylesheet
@author: mscholz
"""
import matplotlib as mpl
import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.collections import LineCollection
from scipy.ndimage.filters import gaussian_filter1d

################################################
#
# custom color palette
################################################
# coolors palette
apple = '#08a045'
green = '#8EA604'
rust = '#b02e0c'
navy = '#05324d'
fire = '#d00000'
sunflower = '#ffba08'
steel = '#5F7C8D'
gray = '#878787'
# prediction paper
# shades of red, dark to light
R0, R1, R2 = '#651119ff', '#b0202eff', '#d15144ff'
Rs = [R0, R1, R2]
# shades of blue
B0, B1, B2 = '#2e2f48ff', '#2b497aff', '#647a9eff'
Bs = [B0, B1, B2]
# shades of viridis
V0, V1, V2, V3, V4 = '#403f85ff', '#006e90ff', '#03cea4ff', '#c3de24ff', '#f1e524ff'
Vs = [V0, V1, V2, V3, V4]
# line plot shades
L0, L1, L2, L3 = ['#1a5477ff', '#0d8d9bff', '#ce5c00ff', '#f0a202ff']
Ls = [L0, L1, L2, L3]
# neutrals
N0, N1, N2 = '#383936ff', '#8b8b8bff', '#d1d1d1ff'
Ns = [N0, N1, N2]
# make a transition cmap - can be extended for any colors
transientcmap = mpl.colors.ListedColormap([mpl.colors.to_rgb(B1), mpl.colors.to_rgb(R1)], name='transient', N=None)

################################################
#
# define axis stuff
#
################################################
#mpl.rc('font', **{'sans-serif' : 'FiraSans','family' : 'sans-serif'})
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Deja Vu'#'Fira Sans'
#mpl.rcParams['font.weight'] = 'regular'
#mpl.rcParams['figure.titleweight'] = 'medium'
mpl.rc('text.latex', preamble='\usepackage{sfmath}')
mpl.rcParams['image.cmap'] = 'viridis'

axescolor = 'k'
mpl.rcParams["axes.edgecolor"]=axescolor
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.spines.top"] = False
# text
mpl.rcParams["text.color"]='k'
mpl.rcParams["ytick.color"]=axescolor
mpl.rcParams["xtick.color"]=axescolor
mpl.rcParams["axes.labelcolor"]='k'
mpl.rcParams["savefig.format"] ='pdf'
# change legend properties
mpl.rcParams["legend.frameon"]=False
mpl.rcParams["legend.labelspacing"]=0.25
mpl.rcParams["legend.labelspacing"]=0.25
#mpl.rcParams['text.usetex'] =True
mpl.rcParams["font.size"] = 12
mpl.rcParams["axes.labelsize"]=  18
mpl.rcParams["xtick.labelsize"]=  18
mpl.rcParams["ytick.labelsize"]=  18
mpl.rcParams["axes.labelpad"] = 0

# suddenly this isn't imported from stylesheet anymore...
mpl.rcParams["axes.labelsize"] = 14
mpl.rcParams["xtick.labelsize"] = 14
mpl.rcParams["ytick.labelsize"] = 14
mpl.rcParams["font.size"] = 12
fs = mpl.rcParams["font.size"]
#=============================================================================#
#                           moving axes
#=============================================================================#

def moveAxes(ax, action, step ):
    if action =='left':
        pos = ax.get_position().get_points()
        pos[:,0] -=step
        
    if action =='right':
        pos = ax.get_position().get_points()
        pos[:,0] +=step
        
    if action =='down':
        pos = ax.get_position().get_points()
        pos[:,1] -=step
    if action =='up':
        pos = ax.get_position().get_points()
        pos[:,1] +=step
    if action =='scale':
        pos = ax.get_position().get_points()
        pos[1,:] +=step/2.
        pos[0,:] -=step/2.
    if action =='scaley':
        pos = ax.get_position().get_points()
        pos[1,1] +=step/2.
        pos[0,1] -=step/2.
    if action =='scalex':
        pos = ax.get_position().get_points()
        pos[1,0] +=step/2.
        pos[0,0] -=step/2.
        
    posNew = mpl.transforms.Bbox(pos)
    ax.set_position(posNew)

#=============================================================================#
#                           clean away spines
#=============================================================================#
def cleanAxes(ax, where='all'):
    '''remove plot spines, ticks, and labels. Either removes both, left or bottom axes.'''
    if where=='all':
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_yticks([])
        ax.set_xticks([])
    elif where=='x':
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
    elif where=='y':
        ax.spines['left'].set_visible(False)
        ax.set_yticks([])
    elif where=='ticks':
        ax.set_yticks([])
        ax.set_xticks([])
    else:
        print 'Command not found. Use "x" or "y" or "all"'
    
#=============================================================================#
#                           align two plots
#=============================================================================# 
def alignAxes(ax1, ax2, where='x'):
    """move axes such that the x or y corners align. Reference is ax1, ax2 gets moved."""
    if where =='xspan':
        x0 = ax1.get_position().get_points()[0][0]
        x1 = ax1.get_position().get_points()[1][0]
        pos = ax2.get_position().get_points()
        pos[0][0] = x0
        pos[1][0] = x1
        ax2.set_position(mpl.transforms.Bbox(pos))  
    if where =='yspan':
        y0 = ax1.get_position().get_points()[0][1]
        y1 = ax1.get_position().get_points()[1][1]
        pos = ax2.get_position().get_points()
        pos[0][1] = y0
        pos[1][1] = y1
        ax2.set_position(mpl.transforms.Bbox(pos))  
    if where =='x':
        x0 = ax1.get_position().get_points()[0][0]
        pos = ax2.get_position().get_points()
        diffx = pos[0][0]-x0
        pos[0][0] = x0
        pos[1][0] -= diffx
        ax2.set_position(mpl.transforms.Bbox(pos))  
    if where =='y':
        y0 = ax1.get_position().get_points()[0][1]
        y1 = ax1.get_position().get_points()[1][1]
        pos = ax2.get_position().get_points()
        diffy = pos[0][1]-y0
        pos[0][1] = y0
        pos[1][1] -= diffy
        ax2.set_position(mpl.transforms.Bbox(pos))
        
    else:
        print 'specify alignment, either enter "x" or "y" or "xspan/yspan."'

#=============================================================================#
#                          create inset axis
#=============================================================================#
def createInset(fig, ax, size, loc = (0,0)):
    """creates an inset with a size in percent of in an existing axis."""
    pos = ax.get_position().get_points()
    x0 = ax.get_position().get_points()[0][0]
    y0 = ax.get_position().get_points()[0][1]
    w = pos[1,0] - pos[0,0]
    h =  pos[1,1] - pos[0,1]
    rect = x0 + loc[0]*w, y0+loc[1]*h, size[0]*w, size[1]*h
    
    axInset = fig.add_axes(rect)    
    return axInset

def annotateInset(ax0, ax1, pos0, pos1):
    """two datapoints pos1, pos2 in data coordinates in their respective axes with a line."""
    trans = ax1.transData + ax0.transData.inverted()
    x0, y0 = pos0
    x1, y1 = trans.transform(pos1)
    ax0.plot([x0,x1], [y0, y1], color = 'k', lw = 1, linestyle=':')

#=============================================================================#
#                           colored text legend
#=============================================================================#
def txtLegend(ax, txt, color, loc, fs = 12, hz = 'left', vt = 'center', data = False, **kwargs):
    """Make a colored text instead of a classical legend.
    Location is in axes coordinates"""
    if data:
        ax.text(loc[0] , loc[1], txt, color = color,  horizontalalignment=hz,\
             verticalalignment=vt, fontsize = fs, **kwargs)
    else:
        return ax.text(loc[0] , loc[1], txt, color = color,  horizontalalignment=hz,\
             verticalalignment=vt, transform=ax.transAxes, fontsize = fs, **kwargs)

#=============================================================================#
#                          scale bar
#=============================================================================#
def scalebar(ax, size, txt, x0 = 0.1, color = 'k', vertical = True, lw = 2, pad = 0.05):
    """add a vertical or horizontal scale bar. Size is a tuple in data coordinates.
    """
    axis_to_data = ax.transAxes + ax.transData.inverted()
    if vertical:
        x = axis_to_data.transform((x0, 0))[0]
        ymid = np.min(size) + np.diff(size)*0.5
        ax.vlines(x, size[0], size[1], color = color, linewidth = lw)
        xpad = axis_to_data.transform((x0+pad, 0))[0]
        ax.text(xpad,  ymid, txt, color = color, verticalalignment='center', horizontalalignment='left')
    else:
        y = axis_to_data.transform((0, x0))[1]
        xmid = np.min(size) + np.diff(size)*0.5
        ax.hlines(y, size[0], size[1], color = color, linewidth = lw)
        ypad = axis_to_data.transform((0,x0+pad))[1]
        ax.text(xmid,  ypad, txt, color = color, verticalalignment='center', horizontalalignment='center')
        


#=============================================================================#
#                          shaded error bars
#=============================================================================#
def shadedError(ax, xdata, samples, axis = 0, color = 'r', alpha = 0.5, sem = False):
    """Given an array of shape = (Nsamples, nPoints ) where Nsamples is repeated samples of each xData point,
    it will plot the mean and standard deviation or if sem = True the standrd error of the mean.
    """
    mean = np.mean(samples, axis = axis)
    std = np.std(samples, axis = axis)
    N = samples.shape[axis]
    ax.plot(xdata, mean, color = color, zorder = 5)
    if sem:
        sem = std/np.sqrt(N)
        ax.fill_between(xdata, mean-sem, mean+sem, color = color, alpha = alpha)
    else:
        ax.fill_between(xdata, mean-std, mean+std, color = color, alpha = alpha)
        
        

#=============================================================================#
#                          scatter error bars
#=============================================================================#
def scatterError(ax, xdata, samples, axis = 0, color = 'r', alpha = 0.5, sem = False, marker = 'o'):
    """Given an array of shape = (Nsamples, nPoints ) where Nsamples is repeated samples of each xData point,
    it will plot the mean and standard deviation or if sem = True the standrd error of the mean.
    """
    mean = np.mean(samples, axis = axis)
    std = np.std(samples, axis = axis)
    N = samples.shape[axis]
    ax.scatter(xdata, mean, color = color, zorder = 5, marker = marker)
    if sem:
        sem = std/np.sqrt(N)
        ax.errorbar(xdata, mean, yerr=sem, color = color, alpha = alpha, linestyle = 'none')
    else:
        ax.errorbar(xdata, mean, yerr= +std, color = color, alpha = alpha, linestyle = 'none')
    

#=============================================================================#
#                          heatmap with correct extent
#=============================================================================#
def plotHeatmap(T, Y,  ax, vmin=-2, vmax=2):
    """nice looking heatmap for neural dynamics.
	T is the time axis in units of sec, Y is the 2d-array of (Nneurons, Ntimepoints), ax is the axis in which it should be plotted."""
    cax1 = ax.imshow(Y, aspect='auto', interpolation='none', origin='lower',extent=[T[0],T[-1],len(Y),0],vmax=vmax, vmin=vmin)
    ax.set_xlabel('Time (s)')
    ax.set_yticks(np.arange(0, len(Y),25))
    ax.set_ylabel("Neuron")
    return cax1
