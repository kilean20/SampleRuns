import numpy as np

npt = 64
nturn = 1
nturnPlot = 512

NL_nu = 0.3
NL_L  = 1.8
NL_c  = 0.01
NL_t  = 0.4
alfx = np.tan(np.pi*NL_nu)
betx = NL_L/np.sin(2.0*np.pi*NL_nu)
k = 2*alfx/betx


def MLI2norm(data_in,bet0=1.0,sign=1):
    data=data_in.copy()
    data[:,5] = -data[:,5]*bet0
    data[:,1] = (data[:,0]*alfx*sign/np.sqrt(betx) + data[:,1]/(1+data[:,5])*np.sqrt(betx))/NL_c
    data[:,3] = (data[:,2]*alfx*sign/np.sqrt(betx) + data[:,3]/(1+data[:,5])*np.sqrt(betx))/NL_c
    data[:,0] = data[:,0]/(np.sqrt(betx)*NL_c)
    data[:,2] = data[:,2]/(np.sqrt(betx)*NL_c)
    return data
    
def norm2MLI(data_in,bet0=1.0,sign=1):
    data=data_in.copy()
    data[:,1] = (-data[:,0]*alfx*sign + data[:,1])*NL_c/np.sqrt(betx)*(1+data[:,5])
    data[:,3] = (-data[:,2]*alfx*sign + data[:,3])*NL_c/np.sqrt(betx)*(1+data[:,5])
    data[:,0] = data[:,0]*np.sqrt(betx)*NL_c
    data[:,2] = data[:,2]*np.sqrt(betx)*NL_c
    data[:,5] = -data[:,5]/bet0
    return data

def getTBT(npt,nturn,fname='rays.out'):
    TBT = np.loadtxt(fname)
    TBT = TBT[:npt*nturn,:]
    out = np.zeros([npt,nturn,7])
    for i in range(nturn):
      out[:,i,:] = TBT[i*npt:(i+1)*npt,:].reshape([npt,7])
      out[:,i,:] = MLI2norm(out[:,i,:],bet0=0.999942)
    return out

TBT0 = np.loadtxt('./track1/rays.in')
TBT0 = MLI2norm(TBT0,bet0=1.0,sign=-1)

TBT1 = getTBT(npt,nturn,'./track1/rays.out')
TBT2 = getTBT(npt,nturn,'./track2/rays.out')
TBT3 = getTBT(npt,nturn,'./track3/rays.out')

print(TBT0[0,:])
print(TBT1[0,0,:6])
print(TBT2[0,0,:6])
print(TBT3[0,0,:6])
