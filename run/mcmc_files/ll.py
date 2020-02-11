import matplotlib.pyplot as plt
#from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
#from matplotlib import ticker
from matplotlib.colors import LogNorm
plt.rc('text', usetex=False)

# you will likely wish to edit these parameters
npars=8
nbins=40
mcmcfilename='mcmctrace.csv'
limitsfilename='mcmc_names.dat'
outputfilename='mcmc.pdf'
#you might wish to edit these
myfontsize=61.0/(1.0+sqrt(npars))
marginpadding=1.2*myfontsize/18.0
marginpadding=0.55*marginpadding   #smaller to account for turning off ticklabels
figsep=0.3

#hopefully, everything below will just work
squaresize=1.0
totalxsize=2.0*marginpadding+npars*squaresize+(npars-0.5)*figsep
totalysize=totalxsize
plt.figure(figsize=(totalxsize,totalysize))
fig = plt.figure(1)
marginpadding=marginpadding/totalxsize
squaresize=squaresize/totalxsize
figsep=figsep/totalxsize
counts1d = np.ndarray((nbins,),float)
xarray = np.ndarray((nbins,),float)

mcmcdata = np.loadtxt(mcmcfilename,delimiter=',',skiprows=1,comments='#')
mcmcdata = mcmcdata.transpose()
npts = mcmcdata[0].size
print('npts=', npts)
#limitsdata = np.loadtxt(limitsfilename,dtype=np.str)
#limitsdata=limitsdata.inverse()

limitsdata = [np.array(map(str, line.split())) for line in open(limitsfilename)]

for ipar in range (0,npars):
  nstrings = limitsdata[ipar].size
  print ('ipar=',ipar,' nstrings=',nstrings)
  for istring in range(5,nstrings):
    limitsdata[ipar][4] = str(limitsdata[ipar][4])+' '+str(limitsdata[ipar][istring])
  #print 'label=', limitsdata[ipar][4]

for ipanel in range (0,npars):
  for jpanel in range (0,npars):
    xmin=float(limitsdata[ipanel][2])
    xmax=float(limitsdata[ipanel][3])
    ymin=float(limitsdata[jpanel][2])
    ymax=float(limitsdata[jpanel][3])
    ax = fig.add_axes([marginpadding+ipanel*(squaresize+figsep),marginpadding+(npars-jpanel-1)*(squaresize+figsep),squaresize,squaresize])
    ax.tick_params(axis='y',which='both',labelleft='off',labelsize=myfontsize)
    ax.set_yticks(np.arange(ymin,ymax,0.5*(ymax-ymin)), minor=False)
    ax.tick_params(axis='x',which='both',labelbottom='off',labelsize=myfontsize)
    ax.set_xticks(np.arange(xmin,xmax,0.5*(xmax-xmin)), minor=False)
    if ipanel == 0:
      ax.tick_params(axis='y',which='both',left='on',labelleft='on',labelright='off', labelsize=myfontsize)
      ax.set_yticks([ymin,0.5*(ymin+ymax),ymax], minor=False)
      ax.set_yticklabels([])
      #ax.set_yticklabels([ymin,0.5*(ymin+ymax),ymax], minor=False, family='serif',rotation=90)
      ax.set_ylabel(limitsdata[jpanel][4],size=1.5*myfontsize, family='serif')
      
    if ipanel == npars-1:
      ax.yaxis.set_label_position('right')
      ax.tick_params(axis='y',which='both',left='off',right='on',labelleft='off',labelright='on', labelsize=myfontsize)
      ax.set_yticks([ymin,0.5*(ymin+ymax),ymax], minor=False)
      ax.set_yticklabels([])
      #ax.set_yticklabels([ymin,0.5*(ymin+ymax),ymax], minor=False, family='serif', rotation=-90)
      ax.set_ylabel(limitsdata[jpanel][4],size=1.5*myfontsize, family='serif', rotation=-90, va='bottom')
    
    if jpanel == 0:
      ax.xaxis.set_label_position('top')
      ax.tick_params(axis='x',which='both',bottom='off',top='on',labelbottom='off', labeltop='on',labelsize=myfontsize)
      ax.set_xticks([xmin,0.5*(xmax+xmin),xmax], minor=False)
      ax.set_xticklabels([])
      #ax.set_xticklabels([xmin,0.5*(xmax+xmin),xmax], minor=False, family='serif', size=myfontsize)
      ax.set_xlabel(limitsdata[ipanel][4],size=1.5*myfontsize, family='serif')
       
    if jpanel == npars-1:
      ax.tick_params(axis='x',which='both',bottom='on',labelbottom='on',labelsize=myfontsize)
      ax.set_xticks([xmin,0.5*(xmax+xmin),xmax], minor=False)
      ax.set_xticklabels([])
      #ax.set_xticklabels([xmin,0.5*(xmax+xmin),xmax], minor=False, family='serif', size=myfontsize)
      ax.set_xlabel(limitsdata[ipanel][4],size=1.5*myfontsize, family='serif')
      
    if ipanel != jpanel:
      counts, xbins, ybins, Image =  ax.hist2d(mcmcdata[ipanel],mcmcdata[jpanel], bins=nbins, norm=LogNorm(), cmap=plt.cm.YlGnBu)
      counts = counts.transpose()
      maxcounts = np.nanmax(counts)
      altmaxcounts = 0
      for i in range(0,nbins):
        for j in range(0,nbins):
          if counts[i][j] > altmaxcounts:
            altmaxcounts = counts[i][j]
            imax=i
            jmax=j
      
      naverage=0
      countsbar=0
      for i in range(imax-1,imax+2):
        for j in range(jmax-1,jmax+2):
          print ('test: i=',i,' j=',j)
          if i!=imax and j!=jmax and i<nbins and i>=0 and j<nbins and j>=0:
            naverage+=1
            countsbar+=counts[i][j]
      
      countsbar=countsbar/naverage

      print ('ipanel=', ipanel, ' jpanel=', jpanel, ' maxcounts=', maxcounts, ' altmaxcounts=', altmaxcounts, ' countsbar=', countsbar)
      level3=countsbar*exp(-0.5)
      level2=countsbar*exp(-2.0)
      level1=countsbar*exp(-4.5)

   
      ax.contour(counts,extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()], linewidths=3,levels=(level1,level2,level3),colors=('darkred','darkgreen','darkblue'),figsize=(10.0+marginpadding+0.02,10.0+marginpadding+0.02))
      plt.xlim(xmin-0.01*(xmax-xmin),xmax+0.01*(xmax-xmin))
      plt.ylim(ymin-0.01*(ymax-ymin),ymax+0.01*(ymax-ymin))
      
    if ipanel==jpanel:
      for ibin in range(0,nbins):
        maxcount=0
        counts1d[ibin]=0;
        xarray[ibin] = xmin+(ibin+0.5)*(xmax-xmin)/float(nbins)
        
      for i in range(0,npts):
        x=mcmcdata[ipanel][i]
        ibin=np.floor(0.99999*nbins*(x-xmin)/(xmax-xmin))
        if ibin >= nbins or ibin < 0:
          print ('ibin=', ibin, '???', 'x = ', x, ' xmax=',x)
        else:
          counts1d[int(ibin)] = counts1d[int(ibin)]+1
          if counts1d[int(ibin)] > maxcount:
            maxcount = counts1d[int(ibin)]

      counts1d = counts1d/float(maxcount)      
      ax.set_yticks([0,0.5,1.0], minor=False)
      ax.tick_params(axis='y',which='both',labelleft='off',labelright='off',labelsize=myfontsize)
      plt.ylim(0,1.05)
      plt.xlim(xmin-0.01*(xmax-xmin),xmax+0.01*(xmax-xmin))
      plt.plot(xarray,counts1d,linestyle='-',linewidth=3,color='r')
    
#  plt.show()
plt.savefig(outputfilename,format='pdf')
os.system("open -a Preview mcmc.pdf")


quit()
