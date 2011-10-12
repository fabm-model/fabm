from pylab import *
import numpy as np
from datetime import datetime

a=np.loadtxt('output.dat',skiprows=1, dtype={ \
          'names':('date','timetime','swr','temp','salt','n','p','z','d'), \
          'formats':('S10','S8','f4','f4','f4','f4','f4','f4','f4')})

dat=asarray([ list(aa)[2:] for aa in a ])
dt=asarray([ num2date(datestr2num(aa[0]+' '+aa[1])) for aa in a ])

plottime=date2num(dt)-date2num(datetime(2002,1,1,0,0,0))

plot(plottime,dat[:,0],'k-')
xlabel('days')
ylabel(u'swr [W/m\u00b2]')
savefig('swr.png')
