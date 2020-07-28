import pImpactR as impact
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy as copy
import pickle
int = np.vectorize(int)


order = 3
nCore_y = 2
nCore_z = 2
nturn = 100

emitGeomRMS = 3.3e-6
sigma_K = 2.0e-3

NL_t = 0.4
NL_L = 1.8
NL_c = 0.01
NL_nu = 0.3
NL_beta = 0.5*NL_L/np.tan(np.pi*NL_nu)


##########################################
## prepare pipe info
##########################################
pipe_info = np.loadtxt('../../pipeinfo.in')
index_NLcenter = np.argmin(pipe_info[:,1])

L_tot = pipe_info[-1,0]
s_NLmid = pipe_info[index_NLcenter,0]

pipe_info[:,0] = pipe_info[:,0] - s_NLmid
pipe_info[:index_NLcenter,0] = pipe_info[:index_NLcenter,0] + L_tot
pipe_info = np.concatenate((pipe_info[index_NLcenter:,:],pipe_info[:index_NLcenter+1,:]),axis=0)

pipe_info[-1,0] = L_tot

np.savetxt('pipeinfo.in',pipe_info)


##########################################
## read Impact input file
##########################################
beam,lattice = impact.readInputFile('../test_iota_v8_4_SextOff_NLon.in')
beam.nCore_y = nCore_y
beam.nCore_z = nCore_z
ke = beam.kinetic_energy
mass = beam.mass
freq = beam.frequency


##########################################
## prepare particles
##########################################
ap1_x,ap1_y = 4e-3, 5.5e-3
nx = 50
testP = 0.0
sig = 0.015
k=0
q_m = beam.multi_charge.q_m[0]
pTest = []
for x in np.linspace(-ap1_x*0.95,ap1_x*0.95,nx):
    for y in np.linspace(-ap1_y*0.95,ap1_y*0.95,nx):
        circ = np.sqrt((x/ap1_x)**2 + (y/ap1_y)**2)
        arcR = np.sqrt(((x-1.3*ap1_x)/ap1_x)**2 + (y/(1.1*ap1_y))**2)
        arcL = np.sqrt(((x+1.3*ap1_x)/ap1_x)**2 + (y/(1.1*ap1_y))**2)
        if circ < 0.95 and y>0:
            if x> 0 and 0.6 < arcR and (circ > 0.75 or  arcR < 0.8):
                pTest.append([x,0.0,y,0.0,0.0,0.0,q_m,0.0,k])
                k=k+1
            elif x<0 and 0.6 < arcL and (circ > 0.75 or  arcL < 0.8):
                pTest.append([x,0.0,y,0.0,0.0,0.0,q_m,0.0,k])
                k=k+1
pTest0 = np.array(pTest)
pTest1 = copy(pTest0)
pTest2 = copy(pTest0)
npt = len(pTest0)

pTest0[:, 5] = -2*sigma_K*ke
pTest1[:,-1] = pTest1[:,-1] + npt
pTest2[:,5] =  2*sigma_K*ke
pTest2[:,-1] = pTest2[:,-1] + 2*npt

pTest = np.concatenate((pTest0,pTest1,pTest2))

impact.writeParticleData(pTest,ke,mass,freq)
beam.n_particles = 3*npt

##########################################
## prepare lattice
##########################################
cleanLat = impact.clearLattice(lattice)
for item in cleanLat:
    if item.type == 'RFkick':
        item.vmax = 0.0
    if 'length' in item:
        item.n_sckick = int(np.ceil(item.length*50))
        item.n_map = 1

for iNL,elem in enumerate(cleanLat):
    if 'nonlinear' in elem.type:
        break

NL0 = impact.getElem('nonlinear_insert_sliced')

NL0.length = 0.9
NL0.start_position = 0.0
NL0.total_length = 1.8

NL0.tune_advance = 0.3
NL0.strength_t = 0.4
NL0.transverse_scale_c = 0.01

NL0.n_map = 45
NL0.n_sckick = 1

NL1 = copy(NL0)
NL1.start_position = 0.9

cleanLat =  [NL1] + cleanLat[iNL+1:] + cleanLat[:iNL] + [NL0]

# add QFF
cleanLat = impact.addHardEdgeQuad(cleanLat)



##########################################
## identity sextupoles
##########################################

Sext=[]
flagNewSext = False
j=-1
for i,elem in enumerate(cleanLat):
    if elem.type in ['quad','dipole']:
        flagNewSext = False
    if elem.type == 'multipole_thin':
        if not flagNewSext:
            flagNewSext = True
            j = j+1
            Sext.append(elem)
            elem.sext_family = j
        else:
            cleanLat[i] = Sext[j]
            
nSext = len(Sext)

##########################################
## prepare I/O
##########################################

loop = impact.getElem('loop')

writeIn = impact.getElem('write_raw_ptcl')
writeIn.turn = 1
writeIn.file_id = 100000
writeIn.format_id = 2

writetmp = impact.getElem('write_raw_ptcl')
writetmp.format_id = 2
writetmp.file_id = 100000 + nturn

writeOut = impact.getElem('write_raw_ptcl')
writeOut.format_id = 2
writeOut.turn = nturn
writeOut.file_id = 100000 - nturn

latticeF = [loop, impact.getElem('pipeinfo'), writeIn] + cleanLat + [writetmp]



##########################################
## getDA data function
##########################################
from time import sleep

def getDA_data(SextStr,sec=1.0):
    beam.distribution.distribution_type = 'ReadFile'
    
    # ==== run foward (from NL mid point) ====   
    for i,item in enumerate(Sext):
        item.KL_sext = SextStr[i]
        
    loop.turns = 10
    writetmp.turn = 10
    impact.writeInputFile(beam, latticeF)
    impact.run(beam,order=order)
    sleep(sec)
    nptSurvived = impact.readLostAt(-1)
    if nptSurvived < int(3*npt*0.5):
        return 10, nptSurvived
    
    beam.distribution.distribution_type = 'ReadFile_binary'
    beam.distribution.file_id = writetmp.file_id
    impact.writeInputFile(beam, latticeF)
    nSplit = int(nturn/10)
    for i in range(nSplit-1):
        impact.run(beam,order=order)
        sleep(sec)
        nptSurvived = impact.readLostAt(-1)
        if nptSurvived < int(3*npt*0.5):
            return 10*(i+2), nptSurvived
    if nSplit*10 != nturn:
        loop.turns = nturn - nSplit*10
        writetmp.turn = nturn - nSplit*10
        impact.writeInputFile(beam, latticeF)
        impact.run(beam,order=order)
        nptSurvived = impact.readLostAt(-1)
        if nptSurvived < int(3*npt*0.5):
            return nturn, nptSurvived
        

    # ==== run backward  ====
    loop.turns = nturn  
    beam.distribution.distribution_type = 'ReadFile_binary'
    beam.distribution.file_id = writetmp.file_id
    
    cleanLat_backward = impact.getInverseLattice(cleanLat)
    latticeB = [loop] + cleanLat_backward + [writeOut]    
      
    impact.writeInputFile(beam, latticeB)
    sleep(2*sec)
    impact.run(beam,order=order)
    sleep(2*sec)
    
    pDataOut = impact.readParticleData(writeOut.file_id, ke,mass,freq, writeOut.format_id)
    
    iSurvived = np.in1d(pTest[:,-1],pDataOut[:,-1])
    pDataIn = pTest[iSurvived,:]
    
    diff =  np.sqrt(np.abs(
                       (pDataIn[:,0] - pDataOut[:,0])**2/(NL_beta*NL_c*NL_c) + 
                       (pDataIn[:,2] - pDataOut[:,2])**2/(NL_beta*NL_c*NL_c) +
                       (pDataIn[:,1] - pDataOut[:,1])**2*NL_beta/(NL_c*NL_c) + 
                       (pDataIn[:,3] - pDataOut[:,3])**2*NL_beta/(NL_c*NL_c) 
                    ))
    ind = int(pDataIn[:,-1])
    ind0 = ind<  npt
    ind1 = (  npt<=ind) * (ind<2*npt)
    ind2 = (2*npt<=ind) * (ind<3*npt)
    return ind[ind0],ind[ind1],ind[ind2], diff[ind0], diff[ind1], diff[ind2]





##########################################
## collect data function
##########################################

def collect_data(init,sec=1.0,
                 data_regular = {'inputs':[],'outputs':[]},
                 data_chaotic = {'inputs':[],'outputs':[]}):
    for x in init:
        y = getDA_data(x,sec)
        if len(y) == 6:
            data_regular['inputs'].append(copy(x))
            data_regular['outputs'].append(copy(y))
        else:
            data_chaotic['inputs'].append(copy(x))
            data_chaotic['outputs'].append(copy(y))
    
    data = {'regular':data_regular,
            'chaotic':data_chaotic}
            
    return data