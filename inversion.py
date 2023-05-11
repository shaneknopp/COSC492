import numpy
import sys
numpy.set_printoptions(threshold=2000000)    
numpy.set_printoptions(precision=15)

real4 = numpy.float32
real = numpy.float64
real8 = numpy.float64
integer4 = numpy.int32
double_precision = numpy.float64

vals = numpy.empty((6,6), dtype=real4, order='F')
vals[0,0] = 0.96992457
vals[0,1] = 0.24238300
vals[0,2] = 7.70767033E-02
vals[0,3] = -0.14315228
vals[0,4] = -0.12523445
vals[0,5] = 0.10469636

vals[1,0] = 0.14555798
vals[1,1] = 0.16776943
vals[1,2] = 0.20685530
vals[1,3] = -0.66514605
vals[1,4] = 1.1042516
vals[1,5] = 1.0296481

vals[2,0] = 0.20436168
vals[2,1] = 0.19450855
vals[2,2] = 0.29393995
vals[2,3] = -8.56770948E-03
vals[2,4] = 0.14588107
vals[2,5] = 0.91137785

# vals[2,0] = -0.14851128
# vals[2,1] = -0.18504372
# vals[2,2] = 4.90637980E-02
# vals[2,3] = -0.91343498
# vals[2,4] = -1.7378793
# vals[2,5] = 0.46495962

# vals[3,0] = 0.28000993
# vals[3,1] = -0.29721546
# vals[3,2] = 6.37818128E-04
# vals[3,3] = 0.19280085
# vals[3,4] = -0.13585463
# vals[3,5] = 0.68889910

# vals[4,0] = 5.4167457
# vals[4,1] = 5.4537840
# vals[4,2] = -0.49230266
# vals[4,3] = -0.45710397
# vals[4,4] = -0.23354781
# vals[4,5] = 0.21970394

# vals[5,0] = 0.78464282
# vals[5,1] = 0.15630294
# vals[5,2] = 5.41635007E-02
# vals[5,3] = 0.66943502
# vals[5,4] = 1.28932670E-02
# vals[5,5] = -8.18159133E-02

def gohead(slat, slon, delta, azim):
    dtor = real4(3.1415928)/real4(180)
    slt = slat*dtor
    dlta = delta*dtor
    azm = azim*dtor
    flat = numpy.arcsin(numpy.cos(slt)*numpy.sin(dlta)*numpy.cos(azm)+numpy.sin(slt)*numpy.cos(dlta))
    flon = numpy.arctan2(numpy.sin(dlta)*numpy.sin(azm),
                         numpy.cos(slt)*numpy.cos(dlta)-numpy.sin(slt)*numpy.sin(dlta)*numpy.cos(azm))
    flat = flat/dtor
    flon = slon + flon/dtor
    if(flon > real4(360)):
        flon = flon - real4(360)
    if(flon < real4(-360)):
        flon = flon + real4(360)
    return flat, flon

def disthead(slat, slon, flat, flon, flag):
    dtor = real4(3.1415928)/real4(180)
    slt = slat*dtor
    sln = slon*dtor
    flt = flat*dtor
    fln = flon*dtor
    delta = numpy.arccos(numpy.sin(slt)*numpy.sin(flt)+numpy.cos(slt)*numpy.cos(flt)*numpy.cos(fln-sln))
    azim = numpy.arctan2(numpy.sin(fln-sln)*numpy.cos(flt),
                         numpy.sin(flt)*numpy.cos(slt) - numpy.cos(fln-sln)*numpy.cos(flt)*numpy.sin(slt))
    #print('delta1', float(delta))
    #print('azim1', float(azim))
    delta = delta/dtor
    azim = azim/dtor
    return delta, azim

yhi = real4()
yb = real4(1.0E+16)
def anneal(p,y,mp,np,ndim,pb,ftol,iter,temptr):
    #print('1')
    flag = False
    global yhi, yb

    psum = numpy.empty(30000)
    tt=-temptr
    ytry = real4(2)
    ysave = real4()
    ynhi = real4()
    ylo = real4(1)
    #print('3')
    for n in range(ndim): # go to 1
        sum = real4(0)
        for m in range(ndim+1):
            sum = sum+p[m,n]
        psum[n] = sum
    while(1):
        #print('2')
        if(ytry >= ysave and ytry >= ynhi and ytry > ylo and flag):
            #print('3')
            for n in range(ndim): # go to 1
                sum = real4(0)
                for m in range(ndim+1):
                    sum = sum+p[m,n]
                psum[n] = sum
        flag = False
        ilo=integer4(0) # go to 2
        inhi=integer4(0)
        r1 = ran1()
        r2 = ran1()
        #print('start')
        #print(r1)
        #print(r2)
        #print(tt)
        #print(numpy.log(r1))
        #print(numpy.log(r1) * tt)
        #print(numpy.log(r2))
        #print(numpy.log(r2) * tt)
        #print(y[0])
        #print(y[1])
        ihi=integer4(1)
        ylo=y[0]+tt*numpy.log(r1)
        #ynhi=ylo
        yhi=y[1]+tt*numpy.log(r2)
        #print(ylo)
        #print(yhi)
        #print('end')

        if(ylo > yhi):
            #print('4')
            ihi = integer4(0)
            inhi = integer4(1)
            ilo = integer4(1)
            ynhi = yhi
            yhi = ylo
            ylo = ynhi
        
        for i in range(2, ndim+1):
            yt = y[i]+tt*numpy.log(ran1())
            if(yt <= ylo):
                #print('5')
                ilo=i
                ylo=yt
            if(yt > yhi):
                #print('6')
                inhi = ihi
                ynhi = yhi
                ihi = i
                yhi = yt
            elif(yt > ynhi):
                #print('6.1')
                inhi = i
                ynhi = yt
        
        rtol=real4(2)*numpy.abs(yhi-ylo)/(numpy.abs(yhi)+numpy.abs(ylo))

        if(rtol<ftol or iter < integer4(0)):
            #print('7')
            swap=y[0]
            y[0] = y[ilo]
            y[ilo] = swap
            for n in range(ndim):
                swap=p[0,n]
                p[0,n] = p[ilo,n]
                p[ilo,n] = swap
            return
        iter = iter-integer4(2)
        ytry = amotsa(p,y,psum,mp,np,ndim,pb,ihi,real4(-1), tt)
        if(ytry <= ylo):
            #print('8')
            ytry=amotsa(p,y,psum,mp,np,ndim,pb,ihi,real4(2), tt)
        elif(ytry >= ynhi):
            #print('9')
            ysave=yhi
            ytry=amotsa(p,y,psum,mp,np,ndim,pb,ihi,real4(0.5), tt)
            if(ytry >= ysave):
                #print('10')
                for i in range(ndim+1):
                    #print('11')
                    if(i!=ilo):
                        #print('12')
                        for j in range(ndim):
                            psum[j]=real4(0.5)*(p[i,j]+p[ilo,j])
                            p[i,j]=psum[j]
                        #print('psum',psum)
                        y[i] = misfit(psum)
                iter=iter-ndim
                flag = True
        else:
            iter=iter+integer4(1)
            #print('13')

def amotsa(p,y,psum,mp,np,ndim,pb,ihi,fac, tt):
    #print('14')
    #print(p)
    #print(y)
    #print(pb)
    global yhi, yb
    ptry = numpy.empty(6000, dtype=real4)
    ptry.fill(0)
    fac1 = (real4(1)-fac)/real4(ndim)
    fac2 = fac1-fac
    for j in range(ndim):
        ptry[j] = psum[j]*fac1-p[ihi,j]*fac2
        #print(ptry[j])
    ytry=misfit(ptry)
    #print(ptry)
    #print(ytry)
    #print(yb)
    if(ytry<=yb):
        #print('15')
        for j in range(ndim):
            pb[j]=ptry[j]
        yb=ytry
    yflu=ytry-tt*numpy.log(ran1())
    #print(yflu)
    #print(yhi)
    if(yflu<yhi):
        #print('16')
        y[ihi] = ytry
        yhi=yflu
        for j in range(ndim):
            psum[j] = psum[j] - p[ihi,j] + ptry[j]
            p[ihi,j] = ptry[j]
    #print('17')
    return yflu


idum = integer4()
iy = integer4()
iv2 = numpy.empty(32, dtype=integer4)
iv2.fill(0)
def ran1():
    global iy
    global idum
    IA = integer4(16807)
    IM = integer4(2147483647)
    AM = real4(1)/real4(IM)
    IQ = integer4(127773)
    IR = integer4(2836)
    NTAB = integer4(32)
    NDIV = integer4(integer4(1)+(IM-integer4(1))//NTAB)
    EPS=real4(1.2e-7)
    RNMX=real4(1)-EPS
    if(idum < 0 or iy == 0):
        idum = numpy.maximum(-idum,integer4(1))
        for i in range(NTAB+7, -1, -1): # changed index
            k=int(idum/IQ)
            idum=IA*(idum-k*IQ)-IR*k
            if(idum < 0):
                idum = idum + IM
            if(i < NTAB):
                iv2[i] = idum
        iy=iv2[0]
    k=integer4(idum/IQ)
    idum=IA*(idum-k*IQ)-IR*k
    if(idum < 0):
        idum = idum+IM
    j = integer4(integer4(1)+iy/NDIV)-integer4(1)
    iy=iv2[j]
    iv2[j] = idum
    return min(AM*iy, RNMX)

def dludcmp(a, n, np, indx):
    d = double_precision(1)
    sum = double_precision()
    dum = double_precision()
    vv = numpy.empty(10000, dtype=double_precision, order='F')

    for i in range(n):
        aamax = double_precision()
        for j in range(n):
            if(numpy.abs(a[i,j]) > aamax): # finds max absolute value in np x np section of a in each row
                aamax = numpy.abs(a[i,j])
        if(aamax == real8(0)): 
            print('singular matrix in ludcmp')
            exit()
        vv[i] = real8(1)/aamax # nornmalization

    for j in range(n):
        for i in range(j): # no -1
            sum = a[i,j]
            for k in range(i): # no -1
                sum = sum - a[i,k]*a[k,j]
            a[i,j] = sum
        aamax = double_precision()
        for i in range(j,n): # should be good
            sum = a[i,j]
            for k in range(j): # no -1
                sum = sum - a[i,k]*a[k,j]
            a[i,j] = sum
            dum = vv[i]*numpy.abs(sum)
            if(dum >= aamax):
                imax = i
                aamax = dum
        if(j != imax):
            for k in range(n):
                dum = a[imax, k]
                a[imax, k] = a[j,k]
                a[j,k] = dum
            d=-d
            vv[imax] = vv[j]
        indx[j] = imax
        if(a[j,j] == real8(0)):
            a[j,j] = real8(1.0E-20)
        if(j!=n-1): # n should be -1
            dum = real8(1)/a[j,j]
            for i in range(j+1, n):
                a[i,j] = a[i,j]*dum
    return d



def dlubksb(a, n, np, indx, b):
    ii = -1
    sum = double_precision()
    for i in range(n):
        ll = indx[i]
        sum = b[ll]
        b[ll] = b[i]
        if(ii!=-1):
            for j in range(ii, i): # no -1
                sum=sum-a[i,j]*b[j]
        elif(sum != 0):
            ii=i
        b[i] = sum
    for i in range(n-1,-1,-1):
        sum=b[i]
        for j in range(i+1,n):
            sum=sum-a[i,j]*b[j]
        b[i] = sum/a[i,i]
    return

maxnfreq=integer4(2)
maxnsta=integer4(400)
maxpts=integer4(8192)
nparam = integer4(10000)
maxnobs = integer4(30000)
maxnodes = integer4(1500)
maxevnts = integer4(400)
maxnxints = integer4(501)
maxnyints = integer4(501)
ndeg = integer4(81)

dph = numpy.empty((maxevnts,maxnsta), order='F', dtype=real4)
staph = numpy.empty((maxevnts,maxnsta,maxnfreq), order='F', dtype=real4)
staamp = numpy.empty((maxevnts,maxnsta,maxnfreq), order='F', dtype=real4)
freq = numpy.empty(maxnfreq, dtype=real4)
g = numpy.empty((maxnobs,nparam), order='F', dtype=real4)
stddevdata = numpy.empty(maxevnts, dtype=real4)
stddevdata.fill(0.0)
change = numpy.empty(nparam, dtype=double_precision)
gtdcmm = numpy.empty(nparam, dtype=double_precision)
gtg = numpy.empty((nparam,nparam), order='F', dtype=double_precision)
gtginv = numpy.empty((nparam,nparam), order='F', dtype=double_precision)
ddd = double_precision()
savegtg = numpy.empty((nparam,nparam), order='F', dtype=double_precision)
chmax = real8(1)
gtd = numpy.empty(nparam, dtype=real4)
stddev = numpy.empty(nparam, dtype=real4)
covinv = numpy.empty(nparam, dtype=real4)
stadist = numpy.empty((maxevnts,maxnsta), order='F', dtype=real4)
staazi = numpy.empty((maxevnts,maxnsta), order='F', dtype=real4)
ysta = numpy.empty((maxevnts,maxnsta), order='F', dtype=real4)
xsta = numpy.empty((maxevnts,maxnsta), order='F', dtype=real4)
streal = numpy.empty((maxevnts,maxnsta,maxnfreq), order='F', dtype=real4)
stimag = numpy.empty((maxevnts,maxnsta,maxnfreq), order='F', dtype=real4)
rloc = numpy.empty((maxevnts,maxnsta), order='F', dtype=real4)
azloc = numpy.empty((maxevnts,maxnsta), order='F', dtype=real4)
attnfac = numpy.empty(maxnfreq, dtype=real4)
xmin = numpy.empty((maxevnts,ndeg), order='F', dtype=real4)
amprms = numpy.empty((maxevnts,maxnfreq), order='F', dtype=real4)

bazi = numpy.empty((maxevnts,maxnsta), order='F', dtype=real4)
stadelt = numpy.empty((maxevnts,maxnsta), order='F', dtype=real4)
d = numpy.empty(maxnobs, dtype=real4)
beg = numpy.empty((maxevnts,maxnsta), order='F', dtype=real4)
origmod = numpy.empty(nparam, dtype=real4)
crrntmod = numpy.empty(nparam, dtype=real4)
boxlat = numpy.empty(4, dtype=real4)
boxlon = numpy.empty(4, dtype=real4)
applat = real4(1)
applon = real4(1)
nodelat = numpy.empty(maxnodes, dtype=real4)
nodelon = numpy.empty(maxnodes, dtype=real4)
nodevel = numpy.empty(maxnodes, dtype=real4)
nodecos2 = numpy.empty(maxnodes, dtype=real4)
nodesin2 = numpy.empty(maxnodes, dtype=real4)
xbox = numpy.empty((maxevnts, 4), order='F', dtype=real4)
ybox = numpy.empty((maxevnts, 4), order='F', dtype=real4)
xnode = numpy.empty((maxevnts, maxnodes), order='F', dtype=real4)
ynode = numpy.empty((maxevnts, maxnodes), order='F', dtype=real4)
stalat = numpy.empty((maxevnts, maxnsta), order='F', dtype=real4)
stalon = numpy.empty((maxevnts, maxnsta), order='F', dtype=real4)
tazim = real4(1)
delta = real4(1)
tazimref = real4(1)
bazimnd = numpy.empty(maxnodes, dtype=real4)

wgtnode = numpy.empty((maxnsta, maxnodes, ndeg), order='F', dtype=real4)
wgtnode1 = numpy.empty((maxnsta, maxnodes,ndeg), order='F', dtype=real4)
ampwgtnode1 = numpy.empty((maxnsta, maxnodes, ndeg), order='F', dtype=real4)

sensitivity = numpy.empty((maxnxints, maxnyints), order='F', dtype=real4)
ampsens = numpy.empty((maxnxints, maxnyints), order='F', dtype=real4)
dtime = numpy.empty(maxnsta, dtype=real4)
avslow = numpy.empty(maxnsta, dtype=real4)
dtime1 = numpy.empty(maxnsta, dtype=real4)
dtime2 = numpy.empty(maxnsta, dtype=real4)
avslow1 = numpy.empty(maxnsta, dtype=real4)
avslow2 = numpy.empty(maxnsta, dtype=real4)
appvel = numpy.empty(maxnodes, dtype=real4)
vage = numpy.empty(maxnodes, dtype=real4)
cos2node = numpy.empty((maxnodes, ndeg), order='F', dtype=real4)
sin2node = numpy.empty((maxnodes, ndeg), order='F', dtype=real4)
startamp1 = numpy.empty(maxevnts, dtype=real4)
startamp2 = numpy.empty(maxevnts, dtype=real4)
stazim1 = numpy.empty(maxevnts, dtype=real4)
stazim2 = numpy.empty(maxevnts, dtype=real4)
stphase1 = numpy.empty(maxevnts, dtype=real4)
stphase2 = numpy.empty(maxevnts, dtype=real4)
pv = numpy.empty(6, dtype=real4)
p = numpy.empty((7, 6), order='F', dtype=real4)
smsft = numpy.empty(7, dtype=real4)
pb = numpy.empty(6, dtype=real4)
annlscl = numpy.empty(6, dtype=real4)
ppb = numpy.empty(6, dtype=real4)
pppb = numpy.empty(6, dtype=real4)
minstd =real4(1)
rmsphase = numpy.empty(maxevnts, dtype=real4)
rmsamp = numpy.empty(maxevnts, dtype=real4)
sortrms = numpy.empty(maxevnts, dtype=real4)
rmsdata = numpy.empty(maxevnts, dtype=real4)
damp1per = numpy.empty((maxevnts, maxnsta), order='F', dtype=real4)
damp2per = numpy.empty((maxevnts, maxnsta), order='F', dtype=real4)
dphase1 = numpy.empty((maxevnts, maxnsta), order='F', dtype=real4)
dphase2 = numpy.empty((maxevnts, maxnsta), order='F', dtype=real4)
phase1 = numpy.empty((maxevnts, maxnsta), order='F', dtype=real4)
phase2 = numpy.empty((maxevnts, maxnsta), order='F', dtype=real4)
phase = numpy.empty((maxnsta, ndeg), order='F', dtype=real4)
dphase = numpy.empty((maxnsta, ndeg), order='F', dtype=real4)
dampper = numpy.empty((maxnsta, ndeg), order='F', dtype=real4)
unifvel =real4(1)
ampmult = numpy.empty(maxnsta, dtype=real4)

dxnode = real4()
dynode = real4()
nxkern = integer4()
xbegkern = real4()
dxkern = real4()
nykern = integer4()
ybegkern = real4()
dykern = real4()

nsta = numpy.empty(maxevnts, dtype=integer4)
iref = numpy.empty(maxevnts, dtype=integer4)
idnode = numpy.empty(maxnodes, dtype=integer4)
nfreq = integer4(1)
nstapts = numpy.empty(maxnsta, dtype=integer4)
indx = numpy.empty(nparam, dtype=integer4)
nstacor = integer4(1)
nobs = integer4(1)
istanum = numpy.empty((maxevnts, maxnsta), order='F', dtype=integer4)
istacor = numpy.empty(maxnsta, dtype=integer4)
nnodes = integer4(1)
nevents = integer4(1)
idnum = numpy.empty(maxevnts, dtype=integer4)
nevntsta = numpy.empty(maxnsta, dtype=integer4)
istavar = numpy.empty(maxnsta, dtype=integer4)
blank = numpy.empty
iblank1 = numpy.empty
iblank2 = numpy.empty
iblank3 = numpy.empty

foutput0 = ''
fsummary0 = ''
fvariance0 = ''
foutput = ''
fn = numpy.empty((maxevnts,maxnsta), dtype=object, order='F')
fsummary = ''
dirpath = ''
evtdatatime = ''
fnname = ''
finvrsnodes = ''
fftinput = ''
fvariance = ''
fmaxavamp = ''
fvelarea  = ''
ftemp = ''
fvelin = ''
fvelout0 = ''
fvelout = ''
fvelnodes  = ''
sensfn = ''
dummy = ''

findid = numpy.empty
pi = real4(3.1415928)
convdeg = pi/real4(180)
circ = real4(6371)*convdeg
twopi = pi*real4(2)
pb.fill(0)


def misfit(p):
    phase1 = numpy.empty(maxnsta, dtype=real4)
    phase2 = numpy.empty(maxnsta, dtype=real4)
    dtime1 = numpy.empty(maxnsta, dtype=real4)
    dtime2 = numpy.empty(maxnsta, dtype=real4)
    startamp1 = p[0]
    startamp2 = p[1]
    stazim1 = p[2]
    stazim2 = p[3]
    stphase1 = p[4]
    stphase2 = p[5]

    misfit = real4()

    if(stazim1 < real4(-40)*convdeg):
        stazim1 = real4(-40)*convdeg
    if(stazim1 > real4(40)*convdeg):
        stazim1 = real4(40)*convdeg
    if(stazim2 < real4(-40)*convdeg):
        stazim2 = real4(-40)*convdeg
    if(stazim2 > real4(40)*convdeg):
        stazim2 = real4(40)*convdeg


    ideg1 = integer4(stazim1/convdeg+(ndeg-real4(1))/real4(2))
    ideg2 = integer4(stazim2/convdeg+(ndeg-real4(1))/real4(2))

    for ista in range(nsta[iev]):
        xstatemp = xsta[iev,ista]*numpy.cos(stazim1) + ysta[iev,ista]*numpy.sin(stazim1)   
        ystatemp = -xsta[iev,ista]*numpy.sin(stazim1) + ysta[iev,ista]*numpy.cos(stazim1)
        dtime1[ista] = (xstatemp-xmin[iev,ideg1])*(real4(1)/unifvel)
        #print('xstatemp', xstatemp)
        #print('xmin', xmin[iev,ideg1])
        #print('dtime', dtime1[ista])

    for ista in range(nsta[iev]):
        xstatemp = xsta[iev,ista]*numpy.cos(stazim2) + ysta[iev,ista]*numpy.sin(stazim2)    
        ystatemp = -xsta[iev,ista]*numpy.sin(stazim2) + ysta[iev,ista]*numpy.cos(stazim2)
        dtime2[ista] = (xstatemp-xmin[iev,ideg2])*(real4(1)/unifvel)

    for ista in range(nsta[iev]):
        phase1[ista] = dtime1[ista]*freq[0]
        phase2[ista] = dtime2[ista]*freq[0]

    for ista in range(nsta[iev]):
        prefase1 = ((phase1[ista]+dphase[ista,ideg1]) - (phase1[iref[iev]]+dphase[iref[iev],ideg1])) + stphase1
        prefase2 = ((phase2[ista]+dphase[ista,ideg2]) - (phase2[iref[iev]]+dphase[iref[iev],ideg2])) + stphase2

        staamp1 = startamp1*(real4(1)+dampper[ista,ideg1])
        staamp2 = startamp2*(real4(1)+dampper[ista,ideg2])

        cosph1 = numpy.cos(prefase1*twopi)
        cosph2 = numpy.cos(prefase2*twopi)
        sinph1 = numpy.sin(prefase1*twopi)
        sinph2 = numpy.sin(prefase2*twopi)

        prereal = staamp1*cosph1+staamp2*cosph2
        preimag = real4(-1)*(staamp1*sinph1+ staamp2*sinph2)
        prereal = prereal*ampmult[istavar[istanum[iev,ista]-1]-1]
        preimag = preimag*ampmult[istavar[istanum[iev,ista]-1]-1]

        kreal = ista + naddat
        kimag = kreal + nsta[iev]
        d[kreal] = (streal[iev,ista,ifreq] - prereal)/stddevdata[iev]
        d[kimag] = (stimag[iev,ista,ifreq] - preimag)/stddevdata[iev]
        misfit = misfit + d[kreal]**real4(2) + d[kimag]**real4(2)
    return misfit

eqlist = open('data/eqlist.25.txt', 'r')

nevents = integer4(eqlist.readline())
nobs = integer4()
for iev in range(nevents):
    line = eqlist.readline().split(' ')
    nsta[iev] = integer4(line[0])
    idnum[iev] = integer4(line[1])
    nobs = nobs + integer4(2) * nsta[iev]
    for i in range(nsta[iev]):
        temp = eqlist.readline().splitlines()
        fn[iev][i] = temp[0] # looks suspicious in fortran, but good here

nfreq = integer4(eqlist.readline().replace("\n",""))
tempfreq = real4(eqlist.readline())
#for i in range(int(nfreq)-1):
#    freq[i] = tempfreq
freq[0] = tempfreq
foutput0 = eqlist.readline().replace("\n","")
fsummary0 = eqlist.readline().replace("\n","")
finvrsnodes = eqlist.readline().replace("\n","")
fftinput = eqlist.readline().replace("\n","")
fvariance0 = eqlist.readline().replace("\n","")
dummy = eqlist.readline().replace("\n","")
ftemp = eqlist.readline().replace("\n","")
unifvel = eqlist.readline().replace("\n","")
unifvel = real4(unifvel)
line = eqlist.readline().split(' ')
iterlimit = real4(line[0]) # 0 in fortran, but correct here
wlambda = real4(line[1])
dampvel = real4(line[2])
dampaniso = real4(line[3].replace("\n",""))
fvelnodes = eqlist.readline().replace("\n","")
sensfn = eqlist.readline().replace("\n","")
fvelout0 = eqlist.readline().replace("\n","")
fvelout = fvelout0
#some weird stuff fortran is doing with filenames that I dont think i need

f10 = open(foutput0, 'a')
f11 = open(fsummary0, 'a')
#f12 = open(fftinput, 'a')
f12 = open('phamp.182new2', 'r')
f13 = open(ftemp, 'a')
#f15 = open(finvrsnodes, 'r')
f15 = open('nodes_0.5extended', 'r')
f16 = open(fvariance0, 'a')
f30 = open(fvelout, 'a')
#f66 = open(sensfn, 'r')
f66 = open('sens182s100km.dat', 'r')

for ista in range(maxnsta):
    nevntsta[ista] = 0

#get station numbers from path
for iev in range(nevents):
    for ista in range(nsta[iev]):
        istanum[iev,ista] = integer4(fn[iev,ista][46:49]) # hard coded
        nevntsta[istanum[iev,ista]-1] = nevntsta[istanum[iev,ista]-1] + 1  # note no element in [0] because stations start at 1

# station identification
jstacnt = integer4()
for ista in range(maxnsta):
    if(nevntsta[ista] > 0):
        jstacnt = jstacnt+integer4(1)
        istavar[ista] = jstacnt

for ixkern in range(maxnxints):
    for iykern in range(maxnyints):
        sensitivity[ixkern,iykern] = 0
        ampsens[ixkern,iykern] = 0

line = f66.readline().split(' ')
nxkern = integer4(line[0])
xbegkern = integer4(line[1])
dxkern = integer4(line[2])
line = f66.readline().split(' ')
nykern = integer4(line[0])
ybegkern = integer4(line[1])
dykern = integer4(line[2])

for ixkern in range(nxkern):
    for iykern in range(nykern):
        line = f66.readline().split(' ')
        sensitivity[ixkern,iykern] = line[2]
        ampsens[ixkern,iykern] = line[3]

dummy = f15.readline()
nnodes = integer4(f15.readline())
for i in range(nnodes):
    line = f15.readline().split(' ')
    nodelat[i] = real4(line[0])
    nodelon[i] = real4(line[1])

for i in range(4):
    line = f15.readline().split(' ')
    boxlat[i] = line[0]
    boxlon[i] = line[1]

ncol = f15.readline().split(' ')
ncol = integer4(ncol[0])
line = f15.readline().split(' ')
dxnode = real4(line[0])
dynode = real4(line[1].replace("\n",""))

f60 = open('areanodestartvel.182', 'r')
dummy = f60.readline()
for i in range(nnodes):
    line = f60.readline().split(' ')
    xx = line[0]
    nodevel[i] = real4(line[1])

iarea = real4(1)
for i in range(nnodes):
    idnode[i] = real4(1)

f10.write(str(foutput0) + '\n')
f11.write(str(foutput0) + '\n')
for iev in range(nevents):
    iv = integer4(f12.readline())
    for ista in range(nsta[iev]):
        beg[iev,ista] = real4(f12.readline())
        line = f12.readline().split()
        stadist[iev,ista] = real4(line[0])
        staazi[iev,ista] = real4(line[1])
        bazi[iev,ista] = real4(line[2])
        stadelt[iev,ista] = real4(line[3])
        stalat[iev,ista] = real4(line[4])
        stalon[iev,ista] = real4(line[5])
        line = f12.readline().split()
        for ifreq in range(nfreq):
            staamp[iev,ista,ifreq] = line[0]
            staph[iev,ista,ifreq] = line[1]

f10.write(str(nfreq) + " nfreq")
for ifreq in range(nfreq):
    f10.write(str(freq[ifreq]) + '\n')
    f11.write(str(freq[ifreq]) + '\n')    
    print(freq[ifreq])
    for iev in range(nevents):
        stddevdata[iev] = real4(0.2)

    npnoamp = integer4(6)*nevents + nnodes + integer4(2)*integer4(iarea)
    np = npnoamp + jstacnt
    kj = integer4(nnodes/ncol)
    i6 = integer4(6)*nevents

    annlscl[0] = .25
    annlscl[1] = .25
    annlscl[2] = real4(15)*convdeg
    annlscl[3] = real4(15)*convdeg
    annlscl[4] = .2
    annlscl[5] = .2
    
    for iev in range(nevents):
        ip = (iev)*integer4(6)
        covinv[0+ip] = real4(1)/(real4(.4)**real4(2))
        covinv[1+ip] = real4(1)/(real4(.4)**real4(2))
        covinv[2+ip] = real4(1)/((real4(10)*convdeg)**real4(2))
        covinv[3+ip] = real4(1)/((real4(10)*convdeg)**real4(2))
        covinv[4+ip] = real4(1)/(real4(.25)**real4(2))
        covinv[5+ip] = real4(1)/(real4(.25)**real4(2))

        startamp1[iev] = .5
        startamp2[iev] = .5
        stazim1[iev] = real4(-7)*convdeg
        stazim2[iev] = real4(7)*convdeg
        stphase1[iev] = 0
        stphase2[iev] = 0

        origmod[0+ip] = startamp1[iev]
        origmod[1+ip] = startamp2[iev]
        origmod[2+ip] = stazim1[iev]
        origmod[3+ip] = stazim2[iev]
        origmod[4+ip] = stphase1[iev]
        origmod[5+ip] = stphase2[iev]

        crrntmod[0+ip] = startamp1[iev]
        crrntmod[1+ip] = startamp2[iev]
        crrntmod[2+ip] = stazim1[iev]
        crrntmod[3+ip] = stazim2[iev]
        crrntmod[4+ip] = stphase1[iev]
        crrntmod[5+ip] = stphase2[iev]

    for ii in range(nnodes):
        ip = i6 + ii
        origmod[ip] = nodevel[ii]
        covinv[ip] = real4(1)/(dampvel**real4(2))
        crrntmod[ip] = origmod[ip]

    for i in range(integer4(iarea)):
        ipp = i6 + nnodes + i
        ippp = ipp + integer4(iarea)
        origmod[ipp] = 0
        origmod[ippp] = 0
        covinv[ipp] = real4(1)/(dampaniso**real4(2))
        covinv[ippp] = real4(1)/(dampaniso**real4(2))
        crrntmod[ipp] = origmod[ipp]
        crrntmod[ippp] = origmod[ippp]

    for ii in range(jstacnt):
        ampmult[ii] = 1
        ip = npnoamp+ii
        origmod[ip] = ampmult[ii] # fortran has 1 in array
        crrntmod[ip] = origmod[ip]
        covinv[ip] = real4(1)/(real4(.3)**real4(2))

    varfac2 = real4(10)
    for itp in range(1): # madness
        ityp = nnodes*(itp)
        # right end
        for ijk in range(i6,i6+kj):
            jk = ijk+ityp
            covinv[jk] = covinv[jk]/varfac2
        #left end
        for ijk in range(i6+nnodes-kj,i6+nnodes):
            jk = ijk+ityp
            covinv[jk] = covinv[jk]/varfac2
        #top
        for ijk in range(i6+kj,i6+nnodes-1*kj,kj):
            jk = ijk+ityp
            covinv[jk] = covinv[jk]/varfac2
        #bottom
        for ijk in range(i6+2*kj-1,i6+nnodes-kj+1,kj):
            jk = ijk+ityp
            covinv[jk] = covinv[jk]/varfac2


    for iev in range(nevents):
        amplarge = real4(0)
        amprms[iev,ifreq] = 0
        iref[iev] = 1
        stasmall = real4(100000)
        for ista in range(nsta[iev]):
            amprms[iev,ifreq] = amprms[iev,ifreq] + staamp[iev,ista,ifreq]**real4(2)
            if(staamp[iev,ista,ifreq] > amplarge):
                amplarge = staamp[iev,ista,ifreq]
                iref[iev] = ista

        amprms[iev,ifreq] = numpy.sqrt(amprms[iev,ifreq]/nsta[iev])

        xsta[iev, iref[iev]] = 0
        ysta[iev, iref[iev]] = 0
        rloc[iev, iref[iev]] = 0
        azloc[iev, iref[iev]] = 0

        for ista in range(nsta[iev]):
            if(ista != iref[iev]):
                xsta[iev,ista] = stadist[iev,ista] - stadist[iev,iref[iev]]
                azidiff = staazi[iev,iref[iev]] - staazi[iev,ista]
                if(azidiff > 180):
                    azidiff = azidiff - real4(360)
                if(azidiff < -180):
                    azidiff = azidiff + real4(360)
                ysta[iev,ista] = circ*numpy.sin(stadelt[iev,ista]*convdeg)*azidiff
                rloc[iev,ista] = numpy.sqrt(xsta[iev,ista]*xsta[iev,ista]+ysta[iev,ista]*ysta[iev,ista])
                azloc[iev,ista] = numpy.arctan2(ysta[iev,ista],xsta[iev,ista])

        applat, applon = gohead(stalat[iev,iref[iev]],stalon[iev,iref[iev]],stadelt[iev,iref[iev]],bazi[iev,iref[iev]])
        delta, tazimref = disthead(applat,applon,stalat[iev,iref[iev]],stalon[iev,iref[iev]], 0)

        appcirc = stadist[iev,iref[iev]]/stadelt[iev,iref[iev]]
        for inode in range(nnodes):
            delta, tazim = disthead(applat,applon,nodelat[inode],nodelon[inode], 0)
            delta, bazimnd[inode] = disthead(nodelat[inode],nodelon[inode],applat,applon, 0)
            xnode[iev,inode] = appcirc*delta - stadist[iev,iref[iev]]
            azidiff = tazimref - tazim
            if(azidiff > 180):
                azidiff = azidiff-real4(360)
            if(azidiff < -180):
                azidiff = azidiff+real4(360)
            ynode[iev,inode] = appcirc*numpy.sin(delta*convdeg)*azidiff

        for ibox in range(4):
            delta, tazim = disthead(applat,applon,boxlat[ibox],boxlon[ibox], 1)
            #print('delta' , float(delta))
            #print('lat', float(boxlat[ibox]))
            #print('lon', float(boxlon[ibox]))
            xbox[iev,ibox] = appcirc*delta - stadist[iev,iref[iev]]
            #print('xbox', xbox[iev, ibox])
            azidiff = tazimref - tazim
            if(azidiff > 180):
                azidiff = azidiff-real4(360)
            if(azidiff < -180):
                azidiff = azidiff+real4(360)
            ybox[iev,ibox] = appcirc*numpy.sin(delta*convdeg)*azidiff

        for ista in range(nsta[iev]):
            dph[iev,ista] = staph[iev,ista,ifreq]-staph[iev,iref[iev], ifreq] + freq[ifreq]*(beg[iev,ista]-beg[iev,iref[iev]])
            streal[iev,ista,ifreq]=staamp[iev,ista,ifreq]*numpy.cos(dph[iev,ista]*twopi)/amprms[iev,ifreq]
            stimag[iev,ista,ifreq]=-staamp[iev,ista,ifreq]*numpy.sin(dph[iev,ista]*twopi)/amprms[iev,ifreq]

    for iev in range(nevents):
        for ideg in range(ndeg):
            stazim1[iev] = (real4(ideg) - (real4(ndeg)-real4(1))/real4(2))*convdeg
            iflag = 0
            templen0 = xbox[iev,0]*numpy.cos(stazim1[iev]) + ybox[iev,0]*numpy.sin(stazim1[iev])
            for ibox in range(4):
                xmintemp = xbox[iev,ibox]*numpy.cos(stazim1[iev]) + ybox[iev,ibox]*numpy.sin(stazim1[iev])
                #print('xmintemp', xmintemp)
                #print('templen0', templen0)
                if(xmintemp < templen0):
                    templen0 = xmintemp
                    iflag = ibox
            xmin[iev,ideg] = xbox[iev,iflag]*numpy.cos(stazim1[iev]) + ybox[iev,iflag]*numpy.sin(stazim1[iev])

    iter = integer4(1)
    icnt = integer4(1)
    nobs = nobs #+ integer4(1)

    dstacnt = jstacnt # keep an eye out here

    sumampcor = real4(0)
    while(True):
        print(iter, icnt)
        for ii in range(jstacnt):
            sumampcor = sumampcor + ampmult[ii]
        d[nobs] = (dstacnt - sumampcor)/real4(1.0e-4)
        
        g.fill(0.0)

        for ii in range(jstacnt):
            icol = npnoamp+ii
            g[nobs,icol] = real4(1.0)/real4(1.0e-4)

        naddat = integer4(0)

        for iev in range(nevents):
            print(iev)

            wgtnode1.fill(0)
            ampwgtnode1.fill(0)
            wgtnode.fill(0)
            dphase.fill(0)
            dampper.fill(0)

            for ideg in range(ndeg):
                stazim1[iev] = real4((ideg) - (ndeg-real4(1))/real4(2))*convdeg
                cos2node_scalar = numpy.cos(real4(2)*(convdeg*bazi[iev,iref[iev]] - stazim1[iev]))
                sin2node_scalar = numpy.sin(real4(2)*(convdeg*bazi[iev,iref[iev]] - stazim1[iev]))
                cos2node[:, ideg] = cos2node_scalar
                sin2node[:, ideg] = sin2node_scalar

                jjj = i6+nnodes
                jjjj = jjj+1
                appvel[:nnodes] = crrntmod[i6:jjj] + cos2node_scalar*crrntmod[jjj] + sin2node_scalar*crrntmod[jjjj]

                for ista in range(nsta[iev]):
                    xstatemp = xsta[iev,ista]*numpy.cos(stazim1[iev]) + ysta[iev,ista]*numpy.sin(stazim1[iev])
                    ystatemp = -xsta[iev,ista]*numpy.sin(stazim1[iev]) + ysta[iev,ista]*numpy.cos(stazim1[iev])
                    senssum1 = 0

                    for ii in range(nnodes):
                        xnodetemp = xnode[iev,ii]*numpy.cos(stazim1[iev]) + ynode[iev,ii]*numpy.sin(stazim1[iev])
                        ynodetemp = -xnode[iev,ii]*numpy.sin(stazim1[iev]) + ynode[iev,ii]*numpy.cos(stazim1[iev])
                        xstanode = xnodetemp - xstatemp
                        ystanode = ynodetemp - ystatemp
                        if(xnodetemp >= xmin[iev,ideg]):
                            ixindex = integer4(xstanode/dxkern) + integer4((nxkern+integer4(1))/integer4(2)) - integer4(1) 
                            iyindex = integer4(ystanode/dykern) + integer4((nykern+integer4(1))/integer4(2)) - integer4(1)
                            if(not (ixindex < 0 or ixindex > nxkern or iyindex < 0 or iyindex > nxkern)):
                                wgtnode1[ista,ii,ideg] = sensitivity[ixindex,iyindex]*(dxnode*dynode)/(dxkern*dykern)
                                ampwgtnode1[ista,ii,ideg] = ampsens[ixindex,iyindex]*(dxnode*dynode)/(dxkern*dykern)

                    dphase[ista,ideg] = (real4(1)/twopi)/unifvel*numpy.sum(wgtnode1[ista,:,ideg]*(appvel - unifvel))
                    dampper[ista,ideg] = numpy.sum(ampwgtnode1[ista,:,ideg]*(appvel-unifvel))/unifvel

            sys.stderr.write('before anneal')
            ip = (iev)*integer4(6) #removed iev-1
            #start simulated annealing setup
            nrestart = integer4(3)
            bestfit = real4(1.0E+16)
            irestart = integer4(0)
            #set up original simplex of points to start simulated annealing process
            #starting with current model and one point perturbed in direction of 
            #each variable
            for ismplx in range(7):
                iprtrb = ismplx - integer4(1)
                for ismp in range(6):
                    p[ismplx,ismp] = crrntmod[ip+ismp]
                    if(ismp == iprtrb):
                        p[ismplx,ismp] = crrntmod[ip+ismp]+annlscl[ismp]
                    pv[ismp] = p[ismplx,ismp]
                
                #calculate misfit at each vertex
                #print('pv', pv)
                smsft[ismplx] = misfit(pv)
                #print('misfit', smsft[ismplx])

            # set starting temperature, tolerance, cooling factor, loops, and iteration counter
            flag = False
            temptr = real4(100)
            tempinit = temptr
            cfac = real4(0.5)
            nloop = integer4(14)
            niter = integer4(100)
            ftol = real4(0.0001)
            loopann = integer4(0)
            while(True):
                if(irestart <= nrestart and not loopann <= nloop and flag):
                    temptr = real4(100) # 498
                    tempinit = temptr
                    cfac = real4(0.5)
                    nloop = integer4(14)
                    niter = integer4(100)
                    ftol = real4(0.0001)
                    loopann = integer4(0)

                iterann = niter # 499
                #print('begin anneal')
                #print(p)
                #print(smsft)
                #print(pb)
                #print('begin anneal')
                anneal(p,smsft,integer4(7),integer4(6),integer4(6),pb,ftol,iterann,temptr)
                loopann = loopann + integer4(1)
                if(loopann <= nloop):
                    temptr = temptr*cfac
                    continue

                irestart = irestart + integer4(1)

        # use three other starting models to make sure that best model is not missed
        # due to difficulty of finding absolute minimum
                if(irestart <= nrestart):
                    flag = True
                    if((irestart==1) or (irestart==2)):
                        pppb[0] = 0.7
                        pppb[1] = 0.3
                        pppb[4] = 0.0
                        pppb[5] = 0.0
                        pppb[2] = real4()*convdeg
                        pppb[3] = real4(10)*convdeg
                    if(irestart == 1):
                        pppb[3] = real4(10)*convdeg
                    if(irestart == 3):
                        pppb[0] = 0.5
                        pppb[1] = 0.5
                        pppb[4] = 0.0
                        pppb[5] = 0.0
                        pppb[2] = real4(7)*convdeg
                        pppb[3] = real4(-7)*convdeg
                    for ismplx in range(7):
                        iprtrb = ismplx-integer4(1)
                        for ismp in range(6):
                            p[ismplx,ismp] = pppb[ismp]
                            if(ismp == iprtrb):
                                p[ismplx, ismp] = pppb[ismp] + annlscl[ismp]
                            pv[ismp] = p[ismplx,ismp]

                        # calculate misfit at each vertex
                        smsft[ismplx] = misfit(pv)
                else:
                    break

            #print('end anneal')
            #print(p)
            #print(smsft)
            #print(pb)
            #print('end anneal')
            #sys.stderr.write('after anneal')
            startamp1[iev] = pb[0]
            startamp2[iev] = pb[1]
            stazim1[iev] = pb[2]
            stazim2[iev] = pb[3]
            stphase1[iev] = pb[4]
            stphase2[iev] = pb[5]
            #startamp1[iev] = vals[iev, 0]
            #startamp2[iev] = vals[iev, 1]
            #stazim1[iev] = vals[iev, 2]
            #stazim2[iev] = vals[iev, 3]
            #stphase1[iev] = vals[iev, 4]
            #stphase2[iev] = vals[iev, 5]

            if(stazim1[iev] < -38*convdeg):
                stazim1[iev] = real4(-38)*convdeg
            if(stazim1[iev] > 38*convdeg):
                stazim1[iev] = real4(38)*convdeg
            if(stazim2[iev] < -38*convdeg):
                stazim2[iev] = real4(-38)*convdeg
            if(stazim2[iev] > 38*convdeg):
                stazim2[iev] = real4(38)*convdeg

            ideg1 = integer4(stazim1[iev]/convdeg+(ndeg-real4(1))/real4(2))
            ideg2 = integer4(stazim2[iev]/convdeg+(ndeg-real4(1))/real4(2))

            for ista in range(nsta[iev]):
                xstatemp = xsta[iev,ista]*numpy.cos(stazim1[iev])+ysta[iev,ista]*numpy.sin(stazim1[iev])
                ystatemp = -xsta[iev,ista]*numpy.sin(stazim1[iev])+ysta[iev,ista]*numpy.cos(stazim1[iev])
                dtime1[ista] = (xstatemp-xmin[iev,ideg1])*(real4(1)/unifvel)
            
            #calculate for second plane wave

            for ista in range(nsta[iev]):
                xstatemp = xsta[iev,ista]*numpy.cos(stazim2[iev])+ysta[iev,ista]*numpy.sin(stazim2[iev])
                ystatemp = -xsta[iev,ista]*numpy.sin(stazim2[iev])+ysta[iev,ista]*numpy.cos(stazim2[iev])
                dtime2[ista] = (xstatemp-xmin[iev,ideg2])*(real4(1)/unifvel)

            # find minimum xbox
            for ista in range(nsta[iev]):
                phase1[iev,ista] = dtime1[ista]*freq[0]
                phase2[iev,ista] = dtime2[ista]*freq[0]
                dphase1[iev,ista] = dphase[ista,ideg1]
                dphase2[iev,ista] = dphase[ista,ideg2]
                damp1per[iev,ista] = dampper[ista,ideg1]
                damp2per[iev,ista] = dampper[ista,ideg2]


            #calulate the change of ph1,ph2, damp1,damp2 with respect to azimuth at the reference station
            #number of intervals = nints, interval length = actxint ista

            aziminc = 2*convdeg
            stazim1[iev] = stazim1[iev] + aziminc
            stazim2[iev] = stazim2[iev] + aziminc

            ideg1 = integer4(stazim1[iev]/convdeg+(ndeg-1)/2)
            ideg2 = integer4(stazim2[iev]/convdeg+(ndeg-1)/2)

            xstatemp = xsta[iev,iref[iev]]*numpy.cos(stazim1[iev]) + ysta[iev,iref[iev]]*numpy.sin(stazim1[iev])        
            ystatemp = -xsta[iev,iref[iev]]*numpy.sin(stazim1[iev]) + ysta[iev,iref[iev]]*numpy.cos(stazim1[iev])  

            # normalize weights to be equivalent to distance
            dtime1ref = (xstatemp-xmin[iev,ideg1])*(real4(1.0)/unifvel)

            # calculate for sencond plane wave ist
            # initialize weights
            xstatemp = xsta[iev,iref[iev]]*numpy.cos(stazim2[iev]) + ysta[iev,iref[iev]]*numpy.sin(stazim2[iev])               
            ystatemp = -xsta[iev,iref[iev]]*numpy.sin(stazim2[iev]) + ysta[iev,iref[iev]]*numpy.cos(stazim2[iev])     

            # normalize weights to be equivalent to distance
            dtime2ref = (xstatemp-xmin[iev,ideg2])*(real4(1.0)/unifvel)

            # initialize weights
            # calculate for first plane wave
            phase1ref = dtime1ref*freq[0]
            phase2ref = dtime2ref*freq[0]
            dphase1ref = dphase[iref[iev],ideg1]
            dphase2ref = dphase[iref[iev],ideg2]
            damp1perref = dampper[iref[iev],ideg1]
            damp2perref = dampper[iref[iev],ideg2]
            stazim1[iev] = stazim1[iev] - aziminc
            stazim2[iev] = stazim2[iev] - aziminc


            # end of for reference station
            for ista in range(nsta[iev]):
                lol = istavar[istanum[iev,ista]-1]-1
                
                # number of intervals = nints, interval length = actxint
                aziminc = 2.*convdeg
                stazim1[iev] = stazim1[iev] + aziminc
                stazim2[iev] = stazim2[iev] + aziminc

                ideg1 = integer4(stazim1[iev]/convdeg+(ndeg-1)/2)
                ideg2 = integer4(stazim2[iev]/convdeg+(ndeg-1)/2)
                xstatemp = xsta[iev,ista]*numpy.cos(stazim1[iev]) + ysta[iev,ista]*numpy.sin(stazim1[iev])            
                ystatemp = -xsta[iev,ista]*numpy.sin(stazim1[iev]) + ysta[iev,ista]*numpy.cos(stazim1[iev])     

                # normalize weights to be equivalent to distance
                dtime1temp = (xstatemp-xmin[iev,ideg1])*(real4(1.0)/unifvel)
    
                # calculate for sencond plane wave
                # initialize weights
                xstatemp = xsta[iev,ista]*numpy.cos(stazim2[iev]) + ysta[iev,ista]*numpy.sin(stazim2[iev])               
                ystatemp = -xsta[iev,ista]*numpy.sin(stazim2[iev]) + ysta[iev,ista]*numpy.cos(stazim2[iev])     

                # normalize weights to be equivalent to distance
                dtime2temp = (xstatemp-xmin[iev,ideg2])*(real4(1.0)/unifvel)

                # initialize weights
                # calculate for first plane wave
                phase1temp = dtime1temp*freq[0]
                phase2temp = dtime2temp*freq[0]
                dphase1temp = dphase[ista,ideg1]
                dphase2temp = dphase[ista,ideg2]
                damp1pertemp = dampper[ista,ideg1]
                damp2pertemp = dampper[ista,ideg2]
                stazim1[iev] = stazim1[iev] - aziminc
                stazim2[iev] = stazim2[iev] - aziminc

                # normalize weights to be equivalent to distance
                parph1azim = ((phase1temp+dphase1temp)-(phase1[iev,ista]+dphase1[iev,ista]))/aziminc - ((phase1ref+dphase1ref) -(phase1[iev,iref[iev]]+dphase1[iev,iref[iev]]))/aziminc
                parph2azim = ((phase2temp+dphase2temp)-(phase2[iev,ista]+dphase2[iev,ista]))/aziminc - ((phase2ref+dphase2ref) -(phase2[iev,iref[iev]]+dphase2[iev,iref[iev]]))/aziminc
                pardamp1azim = (damp1pertemp - damp1per[iev,ista])/aziminc
                pardamp2azim = (damp2pertemp - damp2per[iev,ista])/aziminc

                ideg1 = integer4(stazim1[iev]/convdeg+real4((ndeg-1)/2))
                ideg2 = integer4(stazim2[iev]/convdeg+real4((ndeg-1)/2))

                prefase1 = ((phase1[iev,ista]+dphase1[iev,ista]) - (phase1[iev,iref[iev]]+dphase1[iev,iref[iev]])) + stphase1[iev]
                prefase2 = ((phase2[iev,ista]+dphase2[iev,ista]) - (phase2[iev,iref[iev]]+dphase2[iev,iref[iev]])) + stphase2[iev]

                staamp1 = startamp1[iev]*(real4(1)+damp1per[iev,ista])
                staamp2 = startamp2[iev]*(real4(1)+damp2per[iev,ista])

                cosph1 = numpy.cos(prefase1*twopi)
                cosph2 = numpy.cos(prefase2*twopi)
                sinph1 = numpy.sin(prefase1*twopi)
                sinph2 = numpy.sin(prefase2*twopi)

                prereal = staamp1*cosph1+staamp2*cosph2  
                preimag = -real4(1.0)*(staamp1*sinph1+ staamp2*sinph2)
                prereal = prereal*ampmult[lol]
                preimag = preimag*ampmult[lol]

                # data vect or and partial derivatives listed event by event with all
                # real data for first event followed by imaginary data, then onto next event
                # d contains misfit to starting model
                kreal = ista + naddat
                kimag = kreal + nsta[iev]
                d[kreal] = (streal[iev,ista,ifreq] - prereal)/stddevdata[iev]
                d[kimag] = (stimag[iev,ista,ifreq] - preimag)/stddevdata[iev]

                g[kreal,npnoamp+lol] = prereal/ampmult[lol]/stddevdata[iev]
                g[kimag,npnoamp+lol] = preimag/ampmult[lol]/stddevdata[iev]

                ampadj = ampmult[lol]
                atte1 = real4(1)
                atte2 = real4(1)


                for ii in range(nnodes):
                    parph1v = (real4(1.0)/twopi)*wgtnode1[ista,ii,ideg1]/unifvel - (real4(1.0)/twopi)*wgtnode1[iref[iev],ii,ideg1]/unifvel
                    parph2v = (real4(1.0)/twopi)*wgtnode1[ista,ii,ideg2]/unifvel - (real4(1.0)/twopi)*wgtnode1[iref[iev],ii,ideg2]/unifvel
                    parph1cs = (real4(1.0)/twopi)*cos2node[ii,ideg1]*wgtnode1[ista,ii,ideg1]/unifvel - (real4(1.0)/twopi)*cos2node[ii,ideg1]*wgtnode1[iref[iev],ii,ideg1]/unifvel
                    parph2cs = (real4(1.0)/twopi)*cos2node[ii,ideg2]*wgtnode1[ista,ii,ideg2]/unifvel - (real4(1.0)/twopi)*cos2node[ii,ideg2]*wgtnode1[iref[iev],ii,ideg2]/unifvel
                    parph1sn = (real4(1.0)/twopi)*sin2node[ii,ideg1]*wgtnode1[ista,ii,ideg1]/unifvel - (real4(1.0)/twopi)*sin2node[ii,ideg1]*wgtnode1[iref[iev],ii,ideg1]/unifvel
                    parph2sn = (real4(1.0)/twopi)*sin2node[ii,ideg2]*wgtnode1[ista,ii,ideg2]/unifvel - (real4(1.0)/twopi)*sin2node[ii,ideg2]*wgtnode1[iref[iev],ii,ideg2]/unifvel
                    paramp1v = startamp1[iev]*ampwgtnode1[ista,ii,ideg1]/unifvel 
                    paramp2v = startamp2[iev]*ampwgtnode1[ista,ii,ideg2]/unifvel
                    paramp1cs  = startamp1[iev]*cos2node[ii,ideg1]*ampwgtnode1[ista,ii,ideg1]/unifvel
                    paramp2cs  = startamp2[iev]*cos2node[ii,ideg2]*ampwgtnode1[ista,ii,ideg2]/unifvel
                    paramp1sn  = startamp1[iev]*sin2node[ii,ideg1]*ampwgtnode1[ista,ii,ideg1]/unifvel
                    paramp2sn  = startamp2[iev]*sin2node[ii,ideg2]*ampwgtnode1[ista,ii,ideg2]/unifvel


                    # partial derivatives with respect to velocity, & cos2theta &sin2theta

                    jjjarea=i6 + nnodes + idnode[ii] - 1
                    jjjjarea=jjjarea+integer4(iarea)

                    g[kreal,i6+ii]=(-startamp1[iev]*(real4(1.0)+damp1per[iev,ista])*twopi*sinph1*parph1v*atte1-startamp2[iev]*(real4(1.0)+damp2per[iev,ista])*twopi*sinph2*parph2v*atte2+ paramp1v*cosph1*atte1 + paramp2v*cosph2*atte2)  *ampadj/stddevdata[iev]
                    g[kimag,i6+ii]=(-startamp1[iev]*(real4(1.0)+damp1per[iev,ista])*twopi*cosph1*parph1v*atte1-startamp2[iev]*(real4(1.0)+damp2per[iev,ista])*twopi*cosph2*parph2v*atte2- (paramp1v*sinph1*atte1  +  paramp2v*sinph2*atte2)) *ampadj/stddevdata[iev]
                    g[kreal,jjjarea]=(-startamp1[iev]*(real4(1.0)+damp1per[iev,ista])*twopi*sinph1*parph1cs*atte1-startamp2[iev]*(real4(1.0)+damp2per[iev,ista])*twopi*sinph2*parph2cs*atte2+ paramp1cs*cosph1*atte1 + paramp2cs*cosph2*atte2) *ampadj/stddevdata[iev]+g[kreal,jjjarea]
                    g[kimag,jjjarea]=(-startamp1[iev]*(real4(1.0)+damp1per[iev,ista])*twopi*cosph1*parph1cs*atte1-startamp2[iev]*(real4(1.0)+damp2per[iev,ista]) *twopi*cosph2*parph2cs*atte2 - (paramp1cs*sinph1*atte1  +  paramp2cs*sinph2*atte2))*ampadj/stddevdata[iev]+g[kimag,jjjarea]
                    g[kreal,jjjjarea]=(-startamp1[iev]*(real4(1.0)+damp1per[iev,ista])*twopi*sinph1*parph1sn*atte1-startamp2[iev]*(real4(1.0)+damp2per[iev,ista])*twopi*sinph2*parph2sn*atte2+ paramp1sn*cosph1*atte1 + paramp2sn*cosph2*atte2)*ampadj/stddevdata[iev]+g[kreal,jjjjarea]
                    g[kimag,jjjjarea]=(-startamp1[iev]*(real4(1.0)+damp1per[iev,ista])*twopi*cosph1*parph1sn*atte1 -startamp2[iev]*(real4(1.0)+damp2per[iev,ista])*twopi*cosph2*parph2sn*atte2- (paramp1sn*sinph1*atte1  +  paramp2sn*sinph2*atte2))*ampadj/stddevdata[iev]+g[kimag,jjjjarea]

                #partial derivatives in order are with respect to amplitudes,
                #azimuths, starting phases and  slowness
                ip = (iev)*6 # remove -1

                g[kreal,0+ip] = (real4(1.0)+damp1per[iev,ista])*cosph1*atte1*ampadj/stddevdata[iev]
                g[kreal,1+ip] = (real4(1.0)+damp2per[iev,ista])*cosph2*atte2*ampadj/stddevdata[iev]
                g[kreal,2+ip] = (-startamp1[iev]*(real4(1.0)+damp1per[iev,ista])*sinph1*parph1azim*twopi+startamp1[iev]*pardamp1azim*cosph1)*atte1*ampadj/stddevdata[iev]
                g[kreal,3+ip] = (-startamp2[iev]*(real4(1.0)+damp2per[iev,ista])*sinph2*parph2azim*twopi+ startamp2[iev]*pardamp2azim*cosph2)*atte2*ampadj/stddevdata[iev]
                g[kreal,4+ip] = -startamp1[iev]*(real4(1.0)+damp1per[iev,ista])*sinph1 *twopi*atte1*ampadj/stddevdata[iev]
                g[kreal,5+ip] = -startamp2[iev]*(real4(1.0)+damp2per[iev,ista])*sinph2*twopi*atte2 *ampadj/stddevdata[iev]
                g[kimag,0+ip] = -(real4(1.0)+damp1per[iev,ista])*sinph1*atte1*ampadj/stddevdata[iev]
                g[kimag,1+ip] = -(real4(1.0)+damp2per[iev,ista])*sinph2*atte2*ampadj/stddevdata[iev]
                g[kimag,2+ip] = -(startamp1[iev]*(real4(1.0)+damp1per[iev,ista])*cosph1*parph1azim*twopi+startamp1[iev]*pardamp1azim*sinph1)*atte1*ampadj/stddevdata[iev]
                g[kimag,3+ip] = -(startamp2[iev]*(real4(1.0)+damp2per[iev,ista])*cosph2*parph2azim*twopi+startamp2[iev]*pardamp2azim*sinph2)*atte2*ampadj/stddevdata[iev]
                g[kimag,4+ip] = -startamp1[iev]*(real4(1.0)+damp1per[iev,ista])*cosph1*twopi*atte1*ampadj/stddevdata[iev]
                g[kimag,5+ip] = -startamp2[iev]*(real4(1.0)+damp2per[iev,ista])*cosph2*twopi*atte2*ampadj/stddevdata[iev]
            naddat = naddat + 2*nsta[iev]
            
        # save g and write to compare with Fortran
        # calculate gtg and gtd
        # sys.stderr.write('befote gtg')
        for j in range(np):
            gtd[j] = 0.0
            for i in range(nobs+1):
                gtd[j] = gtd[j] + g[i,j]*d[i]
                #print(gtd[j])

            # add to gtd Tarantola term penalizing misfit to original starting model but skip for wave parameters
            if(j < i6):
                gtdcmm[j] = gtd[j]
                #print(gtdcmm[j])
            if(j >= i6):
                gtdcmm[j] = gtd[j] - covinv[j]*(crrntmod[j]-origmod[j])
                #print(gtdcmm[j])
            for jj in range(0, j+1):
                gtg[jj,j] = 0.0
                for i in range(nobs+1):
                    gtg[jj,j] = gtg[jj,j]+g[i,jj]*g[i,j]

                gtg[j,jj] = gtg[jj,j]
                savegtg[j,jj] = gtg[j,jj]
                savegtg[jj,j] = gtg[jj,j]
            gtg[j,j] = gtg[j,j] + covinv[j]
        #Invert gtg.  gtg will be destroyed.  
        #Not the most efficient approach because doesn't take advantage of 
        #symmetry of gtg.  Use LU decomposition from Press et al.


        for i in range(np):
            gtginv[i,i] = 1.0

        #numpy.save('gtgsave2.npy', savegtg)
        #sys.stderr.write('before dlud')
        #numpy.save('gtg.npy', gtg)
        #gtg = numpy.load('gtgsave.npy')
        ddd = dludcmp(gtg,np,nparam,indx) # indx pass by value
        #numpy.save('gtgmodified.npy', gtg)
        #numpy.save('indx.npy', indx)
        #numpy.save('gtdcmm.npy', gtdcmm)
        #numpy.save('d.npy', d)
        #numpy.save('crrntmod.npy', crrntmod)
        #gtg = numpy.load('gtgmodified.npy')
        #indx = numpy.load('indx.npy')
        #gtdcmm = numpy.load('gtdcmm.npy')
        #d = numpy.load('d.npy')
        #crrntmod = numpy.load('crrntmod.npy')

        #sys.stderr.write('before dlub')
        for j in range(np):
            dlubksb(gtg,np,nparam,indx,gtginv[:,j])
            
        #numpy.save('gtginv.npy', gtginv)
        #gtginv = numpy.load('gtginv.npy')
        #savegtg = numpy.load('gtgsave2.npy')

        for i in range(np):
            change[i] = 0.0
            for j in range(np):
                change[i] = change[i] + gtdcmm[j] * gtginv[i,j]
        
        # Find rank (sum of diagonals of resolution matrix), i.e., number of
        # pieces of information or number of independent model parameters
        # rank1 is contribution from source wave terms, rank2 from velocity
        # variables 
        rank = real4(0)
        rank1 = real4(0.0)
        rank2 = real4(0.0)
        for i in range(np):
            resparm = real4(0.0)
            for j in range(np):
                resparm = resparm + gtginv[i,j]*savegtg[j,i]
            if(i <= i6):
                rank1 = rank1 + resparm
            else:
                rank2 = rank2 + resparm
        
        rank = rank1 + rank2

        # Update current model

        naddat = integer4(0)
        for iev in range(nevents):
            ip = (iev)*6
            startamp1[iev] = startamp1[iev] + change[0+ip]
            startamp2[iev] = startamp2[iev] + change[1+ip]
            stazim1[iev] = stazim1[iev] + change[2+ip]
            stazim2[iev] = stazim2[iev] + change[3+ip]

            if(stazim1[iev] < real4(-40)*convdeg):
                stazim1[iev] = real4(-40)*convdeg
            if(stazim1[iev] > real4(40)*convdeg):
                stazim1[iev] = real4(40)*convdeg
            if(stazim2[iev] < real4(-40)*convdeg):
                stazim2[iev] = real4(-40)*convdeg
            if(stazim2[iev] > real4(40)*convdeg):
                stazim2[iev] = real4(40)*convdeg

            stphase1[iev] = stphase1[iev] + change[4+ip]
            stphase2[iev] = stphase2[iev] + change[5+ip]
            crrntmod[0+ip] = startamp1[iev]
            crrntmod[1+ip] = startamp2[iev]
            crrntmod[2+ip] = stazim1[iev]
            crrntmod[3+ip] = stazim2[iev]
            crrntmod[4+ip] = stphase1[iev]
            crrntmod[5+ip] = stphase2[iev]

            sumsq = real4(0.0)
            sumsqph = real4(0.0)
            
            for ista in range(nsta[iev]):
                dresid1=d[ista+naddat]*stddevdata[iev]
                dresid2=d[ista+naddat+nsta[iev]]*stddevdata[iev]
                predamp = amprms[iev,ifreq]*numpy.sqrt((streal[iev,ista,ifreq]-dresid1)**real4(2) + (stimag[iev,ista,ifreq]-dresid2)**real4(2))
                prefase = numpy.arctan2(-(stimag[iev,ista,ifreq]-dresid2),(streal[iev,ista,ifreq]-dresid1))/twopi
                if(prefase - dph[iev,ista] > real4(0.5)):
                    prefase = prefase-real4(1.0)
                if(prefase - dph[iev,ista] < -real4(0.5)):
                    prefase = prefase+real4(1.0)

                sumsq = sumsq+(d[ista+naddat]**real4(2)+d[ista+naddat+nsta[iev]]**real4(2)) *stddevdata[iev]**real4(2)
                sumsqph = sumsqph + (dph[iev,ista]-prefase)**real4(2)

            rmsph = numpy.sqrt(sumsqph/nsta[iev])/freq[ifreq]
            sigma2 = numpy.sqrt(sumsq/(2*nsta[iev]-real4(8)))
            naddat = naddat + integer4(2)*nsta[iev]

            f10.write('iev ' + str(idnum[iev]) + '  sigma2 ' + str(sigma2) + ' rmsph ' + str(rmsph) + '\n')
            f10.write(str(stazim1[iev]/convdeg) + "\t" + str(startamp1[iev]) + "\t" + str(stphase1[iev]) + '\n')
            f10.write(str(stazim2[iev]/convdeg) + "\t" + str(startamp2[iev]) + "\t" + str(stphase2[iev]) + '\n')

        for ii in range(nnodes):
            crrntmod[ii+i6] = crrntmod[ii+i6] + change[ii+i6]

        for ii in range(integer4(iarea)):
            iii = ii+i6 + nnodes
            iiii = iii+integer4(iarea)
            crrntmod[iii] = crrntmod[iii] + change[iii]
            crrntmod[iiii] = crrntmod[iiii] + change[iiii]
        
        for ii in range(jstacnt):
            ip = ii + npnoamp
            ampmult[ii] = ampmult[ii] + change[ip]
            crrntmod[ip] = ampmult[ii]

        chmax = real4(0.0)
        for ii in range(nnodes):
            ip = i6 + ii
            chmax = numpy.maximum(chmax, numpy.abs(change[ip]))
        icnt = icnt + 1

        print(icnt, chmax)
        if ((chmax > real4(0.0005)) and (icnt <= iterlimit)):
            continue

        if(iter == 1):
            icnttot = icnt
            iter = real4(2)
            minstd = real4(.03)
            naddat = integer4(0)
            for iev in range(nevents):
                sumsq = real4(0.0)
                for ista in range(2*nsta[iev]):
                    sumsq = sumsq + d[ista+naddat]**real4(2)
                sigma2 = sumsq/(real4(2)*nsta[iev]-real4(8))
                stddevdata[iev] = stddevdata[iev]*numpy.sqrt(sigma2)
                if(stddevdata[iev] < minstd):
                    stddevdata[iev] = minstd
                naddat = naddat + 2*nsta[iev]
            icnt = 1
            continue
        icnttot = icnt
        icnttot = icnt
        break

    sumsq = real4(0.0)
    nddat = integer4(0)
    totsumsqph = real4(0.0)
    totsumsqamp = real4(0.0)

    for iev in range(nevents):
        f10.write('event ' + '\n')
        f10.write(str(idnum[iev]) + '\n')
        sumsqph = real4(0.0)
        sumsqamp = real4(0.0)
        sumsqtemp = real4(0.0)
        for ista in range(nsta[iev]):
            dresid1=d[ista+naddat]*stddevdata[iev]
            dresid2=d[ista+naddat+nsta[iev]]*stddevdata[iev]
            predamp = amprms[iev,ifreq]*(numpy.sqrt((streal[iev,ista,ifreq]-dresid1)**2+(stimag[iev,ista,ifreq]-dresid2**2)))
            prefase = numpy.arctan2(-(stimag[iev,ista,ifreq]-dresid2), (streal[iev,ista,ifreq]-dresid1))/twopi

            if(prefase - dph[iev,ista] > real4(0.5)):
                prefase = prefase-real4(1.0)
            if(prefase - dph[iev,ista] < real4(-0.5)):
                prefase = prefase+real4(1.0)
            f10.write(str(ista) + "\t" + str(streal[iev,ista,ifreq]) + "\t" + str(dresid1) + "\t" + str(stimag[iev,ista,ifreq]) + "\t" + str(dresid2) + "\t" + str(staamp[iev,ista,ifreq]) + "\t" + str(predamp) + "\t" + str(dph[iev,ista]) + "\t" + str(prefase) + '\n')

            sumsq = sumsq+(d[ista+naddat]**real4(2)+d[ista+naddat+nsta[iev]]**real4(2)*stddevdata[iev]**real4(2))
            sumsqtemp = sumsqtemp+(d[ista+naddat]**real4(2)+d[ista+naddat+nsta[iev]]**real4(2)*stddevdata[iev]**real4(2))

            sumsqph = sumsqph + (dph[iev,ista]-prefase)**real4(2)
            sumsqamp = sumsqamp + ((predamp-staamp[iev,ista,ifreq])/amprms[iev,ifreq])**real4(2)
        
        totsumsqph = totsumsqph + sumsqph/freq[ifreq]**real4(2)
        totsumsqamp = totsumsqamp + sumsqamp
        rmsphase[iev] = numpy.sqrt(sumsqph/nsta[iev])/freq[ifreq]
        rmsamp[iev] = numpy.sqrt(sumsqamp/nsta[iev])
        rmsdata[iev] = numpy.sqrt(sumsqtemp/(real4(2)*nsta[iev]))
        naddat = naddat + integer4(2)*nsta[iev]

    rmsph = numpy.sqrt(real4(2)*totsumsqph/nobs)
    rmsamplitude = numpy.sqrt(real4(2)*totsumsqamp/nobs)
    sigma2 = numpy.sqrt(sumsq/nobs)

    for j in range(np):
        stddev[j] = numpy.sqrt(gtginv[j, j])

    f10.write('nobs ' + str(nobs) + 'rank ' + str(rank) + ' rank from vel params ' + str(rank2) + '\n')
    f10.write(str(icnttot) + "\t" + str(icnt) + ' iterations ' + str(sigma2) + ' unnormalized rms misfit ' + str(rmsph) + ' rmsphase misfit, s ' + 'r ms amp misfit ' + str(rmsamplitude) + '\n')

    f11.write('nobs ' + str(nobs) + 'rank ' + str(rank) + ' rank from vel params ' + str(rank2) + '\n')
    f11.write(str(icnttot) + "\t" + str(icnt) + ' iterations ' + str(sigma2) + ' unnormalized rms misfit ' + str(rmsph) + ' rmsphase misfit, s ' + ' rms amp misfit ' + str(rmsamplitude) + '\n')

    imed1 = (nevents+integer4(1))//integer4(2)
    imed2 = (nevents+integer4(2))//integer4(2)
    for iev in range(nevents):
        sortrms[iev] = rmsphase[iev]
    
    sortrms = numpy.sort(sortrms)

    rmsmedian = sortrms[imed1]+sortrms[imed2]/real4(2)
    f10.write('median event misfit ' + str(rmsmedian) + '\n')
    f11.write('median event misfit ' + str(rmsmedian) + '\n')

    for iev in range(nevents):
        ip = iev*6
        f10.write('event ' + "\t" + str(idnum[iev]) + "\t" + str(rmsdata[iev]) +  "\t" + ' data std dev ' +  str(rmsphase[iev]) + ' rms phase misfit s' + ' rms amp misfit ' + str(rmsamp[iev]) + '\n')
        f11.write('event ' + "\t" + str(idnum[iev]) + "\t" + str(rmsdata[iev]) +  "\t" + ' data std dev ' + str(rmsphase[iev]) + ' rms phase misfit s' + ' rms amp misfit ' + str(rmsamp[iev]) + '\n')

        wvaz1 = stazim1[iev]/convdeg
        wvaz2 = stazim2[iev]/convdeg
        stdwvaz1 = stddev[ip+2]/convdeg
        stdwvaz2 = stddev[ip+3]/convdeg
        f10.write(str(wvaz1) + "\t" + str(startamp1[iev]) + "\t" + str(stphase1[iev]) + '\n')
        f10.write(str(wvaz2) + "\t" + str(startamp2[iev]) + "\t" + str(stphase2[iev]) + '\n')
        f11.write(str(wvaz1) + "\t" + str(startamp1[iev]) + "\t" + str(stphase1[iev]) + '\n')
        f11.write(str(wvaz2) + "\t" + str(startamp2[iev]) + "\t" + str(stphase2[iev]) + '\n')

    f16.write(str(nnodes) + '\n')
    for ii in range(nnodes):
        for jj in range(nnodes):
            f16.write(str(gtginv[ii+i6, jj+i6]) + '\n')
    
    for ii in range(nnodes):
        f10.write(str(ii) + "\t" + str(crrntmod[ii+i6]) + "\t" + str(stddev[ii+i6]) + '\n')
        f11.write(str(ii) + "\t" + str(crrntmod[ii+i6]) + "\t" + str(stddev[ii+i6]) + '\n')
        f30.write(str(ii) + "\t" + str(crrntmod[ii+i6]) + "\t" + str(stddev[ii+i6]) + '\n')

f30.close()
f26 = open('resmatrix.dat', 'a')
for i in range(i6, i6+nnodes):
    for j in range(i6, i6+nnodes):
        resolution = 0.0
        for k in range(np):
            resolution = resolution+gtg[i,k]*savegtg[k,j]
        f26.write(str(resolution) + '\n')

f26.close()

for ii in range(integer4(iarea)):
    iii + i6 + nnodes + ii
    iiii = iii+integer4(iarea)
    f10.write(str(crrntmod[iii]) + "\t" +  str(stddev[iii]) +  "\t" +  str(crrntmod[iiii]) + "\t" + str(stddev[iiii]) + '\n')
    f11.write(str(crrntmod[iii]) +  "\t" +  str(stddev[iii]) +  "\t" +  str(crrntmod[iiii]) + "\t" + str(stddev[iiii]) + '\n')

print('fifth')
f10.close()
f11.close()
f12.close()
f13.close()
f15.close()
f16.close()
