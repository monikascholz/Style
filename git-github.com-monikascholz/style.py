# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 10:50:27 2015
style file for matplotlib.
@author: monika
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

#=============================================================================#
#                           Define UC colors
#=============================================================================#
UCyellow = ['#FFA319','#FFB547','#CC8214']
UCorange = ['#C16622','#D49464','#874718']
UCred    = ['#8F3931','#B1746F','#642822']
UCgreen  = ['#8A9045','#ADB17D','#616530','#58593F','#8A8B79','#3E3E23']
UCblue   = ['#155F83','#5B8FA8','#0F425C']
UCviolet = ['#350E20','#725663']
UCgray   = ['#767676','#D6D6CE']

UCmain   = '#800000'

#=============================================================================#
#                         Define style of plots
#=============================================================================#
sns.set_style("ticks")
sns.set_context("notebook", font_scale=2, rc={"lines.linewidth": 2, })
axescolor = UCgray[0]
mpl.rcParams["axes.edgecolor"]=axescolor
mpl.rcParams["text.color"]=axescolor
mpl.rcParams["ytick.color"]=axescolor
mpl.rcParams["xtick.color"]=axescolor
mpl.rcParams["axes.labelcolor"]=axescolor
mpl.rcParams["savefig.format"] ='svg'
mpl.rcParams['text.usetex'] =True
mpl.rc('font', **{'sans-serif' : 'FiraSans','family' : 'sans-serif'})
mpl.rc('text.latex', preamble='\usepackage{sfmath}')


#=============================================================================#
#                         Some plotting helper functions
#=============================================================================#
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 11:44:57 2015
Custom plotting libary
@author: monika
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

mpl.rc('font', **{'sans-serif' : 'Arial','family' : 'sans-serif'})
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['ytick.major.pad'] = 5
mpl.rcParams['xtick.major.pad'] = 5
plt.rcParams['pdf.fonttype'] = 42
plt.rc('font', **{'size':18})
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Arial'
mpl.rcParams['mathtext.it'] = 'Arial:italic'
mpl.rcParams['mathtext.bf'] = 'Arial:bold'
fs = 18
#=============================================================================#
#                           Define UC colors
#=============================================================================#
UCyellow = ['#FFA319','#FFB547','#CC8214']
UCorange = ['#C16622','#D49464','#874718']
UCred    = ['#8F3931','#B1746F','#642822']
UCgreen  = ['#8A9045','#ADB17D','#616530','#58593F','#8A8B79','#3E3E23']
UCblue   = ['#155F83','#5B8FA8','#0F425C']
UCviolet = ['#350E20','#725663']
UCgray   = ['#767676','#D6D6CE']
UCmain   = '#800000'

#=============================================================================#
#                           make custom UC colormaps
#=============================================================================#
def make_cmap_from_color(color, dark = False):
    """takes a hex color code and converts to a linear 
    colormap to white or black(dark=True)."""
   
    r,g,b = mpl.colors.hex2color(color) 
    if dark:
        c0 = 0.0
    else:
        c0 = 1.0
    cdict = {'red':   ((0.0,  c0, c0),(1.0,  r, r)),
             'green': ((0.0, c0, c0),(1.0, g, g)),
             'blue':  ((0.0,  c0, c0),(1.0,  b, b))}
    return mpl.colors.LinearSegmentedColormap('customcmap',cdict,256)
#=============================================================================#
#                           Histograms 1D and 2D
#=============================================================================#
def hist2d(x,y,nBinsX,nBinsY,rngX=None,rngY=None):
    if rngX == None and rngY == None:
        h2d, xp, yp = np.histogram2d(y,x,bins=(nBinsY,nBinsX), normed = True)
    else:
        h2d, xp, yp = np.histogram2d(y,x,bins=(nBinsY,nBinsX),range=[rngY,rngX], normed = True)
    extent = [yp[0],yp[-1],xp[0],xp[-1]]
    return h2d, extent
    
def histogram(data, bins, normed=False):
    """simple numpy hist wrapper."""
    hist, bins = np.histogram(data, bins, normed=normed)
    x = bins[:-1]+0.5*(bins[1]-bins[0])
    return x, hist

#=============================================================================#
#               Plotting functions - box, bar and violinplots.
#=============================================================================#
def boxplot(data, labels, colors):
    """custom box plots.
            Args:
               data (array, (n x m)): n categories, m data points.
               labels array (n): x-tick label for each box.
               colors array (n): facecolor of box.
              
            Returns:
                matplotlib boxplot object.
    """
    bp = plt.boxplot(data,whis = [5, 95],\
        patch_artist=True, bootstrap=None, usermedians=None, conf_intervals=None, meanline=False, \
        showmeans=False, showcaps=True, showbox=True, showfliers=True, boxprops={'facecolor':'red', 'lw':1.5}, labels=labels, \
        flierprops={'color':"k", 'lw':1.5}, medianprops={'color':"k", 'lw':1.5}, meanprops=None, capprops={'color':"k", 'lw':1.5},\
        whiskerprops={'color':"k", 'lw':1.5, 'linestyle':'-'}, manage_xticks=True, hold=None)
    for ind, box in enumerate(bp['boxes']):
        box.set(color=colors[ind], linewidth=1.5, alpha=1,edgecolor="None")
    return bp
    
#=============================================================================#
#                           Plotting functions
#=============================================================================#
def barplot(data, labels, colors, rotate=0, sem = False, locations = None, alpha=1, orientation = 'vertical'):
    """custom bar plots.
            Args:
               data (array, (n x m)): n categories, m data points.
               labels array (n): x-tick label for each box.
               colors array (n): facecolor of box.
               rotate (bool): rotate axis labels, keywords: 'vertical' or float for angle
               sem (bool): plot standard error of mean as error (default: Stdev)
               locations array(n): locations of bar plots. 
              
            Returns:
                matplotlib barplot object.
    """
    
    width = 0.5
    if locations==None:
        left = np.linspace(width,width*len(data)*1.5, len(data))
    else:
        left = locations
    if len(data.shape)<2: # deal with uneven lengths
        
        height = [np.nanmean(d) for d in data]
        yerr  = [np.nanstd(d) for d in data]#np.nanstd(data, axis= 1)/np.sqrt(data.shape[1])
        if sem:
            yerr  = [np.nanstd(d)/np.sqrt(len(d)) for d in data]#np.nanstd(data, axis= 1)/np.sqrt(data.shape[1])
    else:
        #left = np.linspace(width,width*data.shape[0]*1.1, data.shape[0])
        height = np.nanmean(data, axis= 1)
        yerr  = np.nanstd(data, axis= 1)#/np.sqrt(data.shape[1])
        if sem:
            yerr /=  np.sqrt(data.shape[1])
    if orientation =='vertical':
        bp = plt.bar(left, height,width, color=colors, alpha=alpha, ecolor='k', linewidth=0, align='center')
        plt.errorbar(left,height, yerr=yerr,linestyle='None', color='k', linewidth = 2)
        plt.xticks(left, labels, rotation=rotate)
        
    else:
        bottom = left
        left = 0
        bp = plt.barh(bottom, height,width,color=colors, alpha=alpha, ecolor='k', linewidth=0, align='center')
        plt.errorbar(height,bottom, xerr=yerr,linestyle='None', color='k', linewidth = 2)
        plt.yticks(bottom, labels, rotation=rotate)
   
    
    return bp    
    
def violinplot(data, labels, colors, positions, points = True, maincolor= 'k'):
    """custom violin plots.
            Args:
               data (array, (m x n)): n categories, m data points.
               labels array (n): x-tick label for each box.
               colors array (n): facecolor of box.
               positions array(n): location of violins.
              
            Returns:
                matplotlib boxplot object.
    """
    w = 0.8#np.diff(positions)[0]*0.8
    
    bp = plt.violinplot(data, positions= positions, vert=True, widths=w, \
    showmeans=False, showextrema=False, showmedians=True, points=100, bw_method=0.25, hold=None)
    plt.boxplot(data,notch=1,positions=positions,vert=1)
    for ind, box in enumerate(bp['bodies']):
        box.set(color=colors[ind], edgecolor=colors[ind],linewidth=1.5, alpha=1)
    bp['cmedians'].set(color= maincolor,linewidth=1.5)
    #bp['cmaxes'].set(color= maincolor,linewidth=1.5)
    #bp['cbars'].set(color= maincolor,linewidth=1.5)
    #bp['cmins'].set(color= maincolor,linewidth=1.5)
    
    plt.xticks(positions,labels, rotation =20)
    
    #### hackish!
    if points:
        for index, d in enumerate(data):
            if len(d)==len(positions) or len(data) != len(positions):
                jitter = positions + (1-2*np.random.rand(len(d)))*w/4.
                plt.scatter(jitter, d, color = maincolor, s = 8)
            else:
                jitter = positions[index] + (1-2*np.random.rand(len(d)))*w/4.
                plt.scatter(jitter, d, color = maincolor, s = 8)
        
        
    return bp

def line_plot(xdata,ydata, color, linestyle='-', alpha = 1):
    ax = plt.plot(xdata,ydata, color = color, lw=2, linestyle=linestyle, alpha = alpha)
    return ax


def histogram_plot(data,bins, color, linestyle='-'):
    b, h = histogram(data, bins, normed=True)
    ax = plt.step(b,h,color=color, lw=2, linestyle=linestyle, where='post')
    return ax, h, b


def scatter_plot(x,y, color):
    ax = plt.plot(x,y,linestyle='none',marker= '.',mfc = 'none',mec = color,\
    markersize=10, mew=2)
    return ax


def joint_histogram(x,y,fig,outer_grid,bins_on, bins_off,color,alpha, xlabel, ylabel):
    """joit histogram and marginals"""
    xmin ,xmax = bins_on[0],bins_on[-1]
    ymin, ymax = bins_off[0],bins_off[-1]
    #xmin, ymin = np.min(x),np.min(y)
    #xmax, ymax = np.max(x),np.max(y)
    #xmax, xmin = tuple(np.array([xmax, xmin]) + 0.15*(xmax - xmin)*np.array([1, -1]))
    #ymax, ymin = tuple(np.array([ymax, ymin]) + 0.15*(ymax - ymin)*np.array([1, -1]))
    gs = mpl.gridspec.GridSpecFromSubplotSpec(2, 2,
            subplot_spec=outer_grid, wspace=-0.1, hspace=-0.1,width_ratios=[3, 1], height_ratios = [1, 4])    
    #gs = mpl.gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios = [1, 4])
    #Create scatter plot
    #ax = plt.subplot(gs[1, 0])
    ax = plt.Subplot(fig, gs[1,0])
    cmap = make_cmap_from_color(color, dark = False)
    sns.kdeplot(x,y, bw=1,gridsize=(ymax-ymin), cut=ymax, clip=(ymin,ymax), cmap=cmap,ax =ax, shade=True)
    #cax = ax.scatter(x, y, color=color, alpha=alpha)
    
    
    #Turn off all axes
    ax.set_xticks([0,5])
    ax.set_yticks([0,5])
    for sp in ax.spines.values():
        sp.set_visible(False)
    
    #Create Y-marginal (right)
    #axr = plt.subplot(gs[1, 1], sharey=ax, xticks=[], yticks=[],frameon = False, xlim=(0, 1), ylim = (ymin, ymax))
    axr = plt.Subplot(fig, gs[1,1], sharey=ax, xticks=[], yticks=[],frameon = False, ylim = (ymin, ymax+1))
    
    axr.hist(y,bins_off, color = color, orientation = 'horizontal', normed = True)
    axr.set_ylabel(ylabel, rotation=-90, va='top')
    axr.set_xlim([0,2])
    #Create X-marginal (top)
    #axt = plt.subplot(gs[0,0], sharex=ax,frameon = False, xticks = [], yticks = [], xlim = (xmin, xmax), ylim=(0, 1))
    axt = plt.Subplot(fig, gs[0,0], sharex=ax, xticks=[], yticks=[],frameon = False, xlim=(xmin, xmax+1))
    axt.set_xlabel(xlabel)     
    axt.hist(x, bins_on,color = color, normed = True)
    axt.set_ylim([0,4])
#    axt.axis('off')
#    axr.axis('off')
    
    fig.add_subplot(ax)
    fig.add_subplot(axr)
    fig.add_subplot(axt)
    
    
    




