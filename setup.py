# modified rdsetupsimulf from DWF.
 
# This program is the part of simultaran that reads in the data and does the
# Fourier analysis.  It is separated so that simultaran can run on katmai or
# other machines that do not have sac.  Writes out the data ready for use in
# simultaran  

import numpy
import pysac
import math
import scipy
import obspy
#import subroutines

# positive fourier transform
def frt(up, fr, nall, delt):
    w = numpy.empty(8192)
    theta = twopi * fr * delt
    c = numpy.cos(theta) * 2

    w[0] = up[nall-1]
    w[1] = c* w[0] + up[nall-2]
    for i in range(2, nall):
        w[i] = c * w[i-1] - w[i-2] + up[nall - i + 1]

    zre = (w[nall-1] - w[nall-3] + up[0]) * delt / 2.0
    zim = w[nall-2] * numpy.sin(theta) * delt 
    prz = phase(zre, zim) # phase
    arz = numpy.sqrt(zre*zre+zim*zim) # amplitude
    return arz, prz


def phase(x, y):
    if x < 0:
        phi = numpy.arctan(y/x) + pi
    if x == 0 and y < 0:
        phi = 1.5*pi
    if x == 0 and y == 0:
        phi = 0
    if x == 0 and y > 0:
        phi = 0.5*pi
    if x > 0 and y < 0:
        phi = numpy.arctan(y/x) + twopi
    if x > 0 and (y == 0 or y > 0):
        phi = numpy.arctan(y/x)
    return phi/twopi-int(phi/twopi)

real4 = numpy.float32
real8 = numpy.float64
integer4 = numpy.int32
double_precision = numpy.float64

maxnfreq, maxnsta, maxpts, maxnodes, maxevnts = 2, 400, 8192, 1500, 400
pi = real4(3.1415928 )
convdeg = pi/real4(180)
circ = real4(6371)*convdeg
twopi = pi*real4(2)
cattnfac = real4(0.25e-3)

staph, staamp = numpy.empty((maxevnts,maxnsta,maxnfreq), dtype=real4), numpy.empty((maxevnts,maxnsta,maxnfreq), dtype=real4)
freq = numpy.empty(maxnfreq, dtype=real4)
stadist, staazi = numpy.empty((maxevnts,maxnsta), dtype=real4), numpy.empty((maxevnts,maxnsta), dtype=real4)
bazi, stadelt = numpy.empty((maxevnts,maxnsta), dtype=real4), numpy.empty((maxevnts,maxnsta), dtype=real4)
stacor, geomsprd = numpy.empty(maxnsta, dtype=real4), numpy.empty(maxnsta, dtype=real4)

beg = numpy.empty((maxevnts,maxnsta), dtype=real4)
delt = numpy.empty(maxnsta, dtype=real4)
stalat = numpy.empty((maxevnts,maxnsta), dtype=real4)
stalon = numpy.empty((maxevnts,maxnsta), dtype=real4)
nsta = numpy.empty(maxevnts, dtype=integer4)
nfreq = integer4()
nstapts = numpy.empty(maxnsta, dtype=integer4)
istanum = numpy.empty(maxnsta, dtype=integer4)
nevents = integer4()
fn = numpy.empty((maxevnts,maxnsta), dtype=object)

# read list of files and frequencies to be analyzed and file to output results
def parse_eqlist(file):
    global nevents, nfreq
    eqlist = open(file, 'r')
    nevents = integer4(eqlist.readline())
    for iev in range(nevents):
        line = eqlist.readline().split(' ')
        nsta[iev] = integer4(line[0])
        for ista in range(nsta[iev]):
            fn[iev][ista] = eqlist.readline().splitlines()[0]

    nfreq = integer4(eqlist.readline())
    for i in range(nfreq):
        freq[i] = real4(eqlist.readline())

def parse_staampcore(file):
    istacor, tempcor = numpy.empty(maxnsta, dtype=integer4), numpy.empty(maxnsta, dtype=real4)
    f3 = open(file, 'r')
    nstacor = integer4(f3.readline())
    for i in range(nstacor):
        line = f3.readline().split()
        istacor[i] = integer4(line[0])
        tempcor[i] = real4(line[1])
    for i in range(nstacor):
        stacor[istacor[i]] = tempcor[i]
    f3.close()

def process_data():
    global nevents, nfreq
    for iev in range(nevents):
        for ista in range(nsta[iev]):
            sac = pysac.SACTrace.read(fn[iev][ista])
            data = sac.data
            nstapts[ista] = len(data)
            beg[iev][ista] = sac.b
            stadist[iev][ista] = sac.dist
            staazi[iev][ista] = sac.az
            bazi[iev][ista] = sac.baz
            stadelt[iev][ista] = sac.gcarc
            stalat[iev][ista] = sac.stla
            stalon[iev][ista] = sac.stlo
            delt[ista] = sac.delta

            geomsprd[ista] = numpy.sqrt(numpy.sin(stadelt[iev][ista] * convdeg))
            data = numpy.roll(data, 2) # need to figure this out

            for ifreq in range(nfreq):
                staamp[iev][ista][ifreq], staph[iev][ista][ifreq] = frt(data[:maxpts], freq[ifreq], maxpts, delt[ista])
                attneffect = math.exp(cattnfac*(stadist[iev,ista]-stadist[iev,0]))
                #staamp[iev][ista][ifreq] = staamp[iev][ista][ifreq]*attneffect*(geomsprd[ista]/stacor[int(istanum[ista])])
                staamp[iev][ista][ifreq] = staamp[iev][ista][ifreq]*attneffect*(geomsprd[ista]/1)

def write(file):
    global nevents, nfreq
    f1 = open(file, 'a')
    for iev in range(nevents):
        f1.write("\t" + str(iev+1) + "\n")
        for ista in range(int(nsta[iev])):
            f1.write("  " + str('{:.4f}'.format(round(beg[iev][ista], 4))) + "\n")
            f1.write("   " + str('{:9.3f}'.format(round(stadist[iev][ista], 3))) + "\t" + str('{:9.5f}'.format(round(staazi[iev][ista], 5))) + "\t" + str('{:9.5f}'.format(round(bazi[iev][ista], 5))) + "\t" + str('{:9.6f}'.format(round(stadelt[iev][ista], 6))) + "\t" + str('{:9.5f}'.format(round(stalat[iev][ista], 6))) + "\t" + str('{:9.5f}'.format(round(stalon[iev][ista], 5))) + "\n")
            for ifreq in range(nfreq):
                f1.write("  " + str(round(staamp[iev][ista][ifreq], 10)) + "  " + str(round(staph[iev][ista][ifreq], 10)) + "\n")
    f1.write(' ')
    f1.close

def main():
    """Run all functions"""
    # change for user input
    parse_eqlist('data/eqlist.25.txt')
    parse_staampcore('staampcore.ave')
    process_data()
    write('phamp.182setup')

if __name__=="__main__":
     main()
