import differential_evolution as diffevol
import numpy as np
import os
import shutil

with open("strHead") as f:
  strHead = f.read()
  
with open("strTail") as f:
  strTail = f.read()


sextName = ['sa1rv = ', 'sc1rv = ', 'sc2rv = ', 'sd1rv = ', 'se1rv = ','se2rv = ',
            'se2lv = ','se1lv = ', 'sd1lv = ', 'sc2lv = ', 'sc1lv = ', 'sa1lv = ',
            'sadd1v = ','sadd2v = ','sadd3v = ','sadd4v = ','sadd5v = ','sadd6v = ']

def getStrSext(data):
    f='\n'
    for i in range(len(sextName)):
      f = f + sextName[i]+str(data[i])+' \n'
    return f+'\n'

sextData = np.random.randn(18)
strSext = getStrSext(sextData)

with open("mli.in","w") as f:
    f.write(strHead+strSext+strTail)


ref3 = np.zeros(18)
ref3[2]=0.727273845
ref3[15]=0.727273845

def objFunc(arg): 

    
    target = pm.opt.id_generator()  # generage random directory name
    while os.path.exists(target):  
        target = pm.opt.id_generator()   
    os.chdir(target) # cd to the randome directory and
    
    strSext = getStrSext(arg)
    with open("mli.in","w") as f:
    f.write(strHead+strSext+strTail)
    os.system('srun -n 1 mli.x')

    G = np.zeros(233) #np.zeros(683-450)
    k=0
    with open("mli.out") as f:
        for i, line in enumerate(f):
            if i >= 450 and i < 683:
                G[k]=float(line[23:-1])
                k=k+1
    # 3rd order
    obj1 = np.sum((G[:18]-ref3)**2)
    # 4rd order 
    obj2 = (twissX.betx - arg[1])**4 + (5.0*twissX.alfx + 5.0*arg[2])**4 +\
           (twissY.bety - arg[1])**4 + (5.0*twissY.alfy - 5.0*arg[2])**4
    os.chdir('..')
    shutil.rmtree(target)
    return obj1 + 10*obj2




data = np.zeros(683-450)

k=0
with open("mli.out") as f:
    for i, line in enumerate(f):
        if i >= 450 and i < 683:
            data[k]=float(line[23:-1])
            k=k+1


diffevol.getNP()
#os.system('time srun -n 1 mli.x')
