#===================================================== 
# ele2ez.py
#
# Read Elegant-generated file with lattice  
# and format it in an 'easy'  way suitable for generating an IMPACTz input file. 
#
# Elegant file has to be prepared using the comand:
# >sddsprintout $1.par $1.elmnts.txt -col=ElementName  -col=ElementParameter -col=ParameterValue -col=ElementType
# 
# The main output file (from this script) is "beamline.txt", with data structure:
#   Name, type, s, length, K0L, e1, e2,  K1L, k2L, tilt, HGAP, FINT, rf Voltage (V), rf Phase, rf Freq (Hz), element no. 
#   Note: rf Phase is in deg; crest is at phase=0 (Litrack convention).
 
# Other output files are 'sbendList.txt' and 'cavList.txt', 
# containing lists of all the bend magnets and cavities;
# 'purgedEleFile.txt' is  an intermediate file for debugging.

# 2-24-2014 MV 
# 3-24-2014 Include element MATR
# 11-17-2014 MV: Change "lcls2sch.elmnts.txt" to "lcls2scs.elmnts.txt"
# 11-24-2015 CM: Include elements RFDF and MONI
#=====================================================   
import os, sys, string, pdb, math
 
#----------------------------
# Element types in the Elegant file recognized by this script and no. of parameters (to be expanded ...)
notypEl  = 6 +3 +1 +2 # dimention of array  typEl
typEl = ['DRIF', 'CSRDRIFT', 'QUAD', 'CSRCSBEND', 'SBEND', 'RFCW', 'RCOL', 'ECOL', 'ROTATE', 'MATR', 'RFDF', 'MONI']
noElPar = [2,       21 ,       28,       63,         1 ,     31,     6,      7,       1,       2,      26,     11] 
 

#----------------------------
# Read in elegant output file

fileele = open('NGLS2.lat','r')
eleline = range(0,90000)
ielefile = 0 
while 1:
   line1 = fileele.readline()
   if not line1: break
   eleline[ielefile] = string.split(line1)
   print '  ',ielefile, eleline[ielefile]
   ielefile = ielefile + 1 


#----------------------------
# Purge by retaining only the physical elements; store result in eleline2
#eleline2 = range(0,ielefile)
#i = 4 # Skip  header to avoid problems with empty records 
#ix = 0
#while (i <= ielefile-1): 
#   j = 0 
#   while (j <= notypEl - 1): 
#      if eleline[i][3] == typEl[j]: 
#         eleline2[ix] = eleline[i]
#         ix = ix + 1
#      j = j + 1
#   i = i + 1


eleline2 = range(0,ielefile)
i = 4 # Skip  header to avoid problems with empty records 
ix = 0
while (i <= ielefile-1): 
   j = 0 
   irecognized = 0
   while (j <= notypEl - 1): 
      if eleline[i][3] == typEl[j]: 
         eleline2[ix] = eleline[i]
         ix = ix + 1
         irecognized = 1  
      j = j + 1
   if  (irecognized==0) and  (eleline[i][1]=='L')  and ( math.fabs(float(eleline[i][2])) > 0.): # Stop if we find unrecognized el. type w/ L>0
    print ' UNRECOGNIZED FINITE-LENGTH ELEMENT TYPE : ', eleline[i][3]
    print '                            ELEMENT NAME : ', eleline[i][0]
    print '                          ELEMENT LENGTH : ', eleline[i][2]
    print '                                       i = ', i
    print ' ... exiting'
    sys.exit(0)     
   i = i + 1



#----------------------------
# Write out purged beamline 
purgedEleFile1 = open('purgedEleFile.txt', 'w')
outx = range(0,4)
i=0   
while (i <= ix-1): 
   outx[0]= eleline2[i][0]+ ' ' # element name
   outx[1]= eleline2[i][1]+ ' ' # name of parameter
   outx[2]= eleline2[i][2]+ ' ' # value of parameter
   outx[3]= eleline2[i][3]+ ' ' # element type
   outxx = string.join(outx) + '\n'
   purgedEleFile1.write(outxx)
   i = i + 1
    
#----------------------------
# Reorganize beam line parameters as follows:
#   Element name, Element type, s,  length, K0L, e1, e2,  K1L, k2L, tilt, HGAP, FINT, rf Voltage, rf Phase, rf Freq, el. no.   
#        0            1         2     3      4   5    6    7    8     9    10     11       12        13       14       15

i = 0   
j = 0
mybeamline = range(0,10000)
icountEl = 0  # counr phys element
icountSbend = 0
icountCav = 0
sbendList = range(0,1000)
cavList =  range(0,10000)
sbendFile = open('sbendList.txt', 'w') 
cavFile = open('cavList.txt', 'w') 
beamlineFile = open('beamline.txt', 'w') 

cavities = range(0,10000)
spath = 0 
while (i <= ix-1): 
   mybeamlineX = ['0 ','0 ',  '0 ','0 ',  '0 ','0 ',  '0 ','0 ',  '0 ','0 ',  '0 ','0 ',  '0 ','0 ',  '0 ']
   
   if  (eleline2[i][3]=='DRIF'):
        spath = spath + float(eleline2[i][2]) 
        mybeamlineX[0] = eleline2[i][0]+ ' '  # Name 
        mybeamlineX[1] = eleline2[i][3]+ ' '  # Type
        mybeamlineX[2] = str(spath)           # s
        mybeamlineX[3] = eleline2[i][2]+ ' '  # Length      
        icountEl = icountEl + 1                      
        beamlineFile.write(string.join(mybeamlineX) + ' '+str(icountEl) + ' \n')
        i = i + (noElPar[0] -1) 
        
        
   elif  (eleline2[i][3]=='RCOL'):              # Treat this a drift
        spath = spath + float(eleline2[i][2])    
        mybeamlineX[0] = eleline2[i][0]+ ' '  # Name 
        mybeamlineX[1] = 'DRIF'+ ' '          # Type
        mybeamlineX[2] = str(spath)           # s        
        mybeamlineX[3] = eleline2[i][2]+ ' '  # Length      
        icountEl = icountEl + 1                      
        beamlineFile.write(string.join(mybeamlineX) + ' '+str(icountEl) + ' \n') 
        i = i + (noElPar[6] -1)     
        
   elif  (eleline2[i][3]=='ECOL'):            # Treat this a drift
        spath = spath + float(eleline2[i][2])    
        mybeamlineX[0] = eleline2[i][0]+ ' '  # Name 
        mybeamlineX[1] = 'DRIF'+ ' '          # Type
        mybeamlineX[2] = str(spath)           # s
        mybeamlineX[3] = eleline2[i][2]+ ' '  # Length      
        icountEl = icountEl + 1                      
        beamlineFile.write(string.join(mybeamlineX) + ' '+str(icountEl) + ' \n') 
        i = i + (noElPar[7] -1)                

   elif  (eleline2[i][3]=='ROTATE'):          
        spath = spath + float(eleline2[i][2])       
        mybeamlineX[0] = eleline2[i][0]+ ' '  # Name 
        mybeamlineX[1] = eleline2[i][3]+ ' '  # Type
        mybeamlineX[2] = str(spath)           # s
        mybeamlineX[3] = ' 0.0'+ ' '          # Length     
        mybeamlineX[9] = eleline2[i][2]+ ' '  # Tilt angle             
        icountEl = icountEl + 1                      
        beamlineFile.write(string.join(mybeamlineX) + ' '+str(icountEl) + ' \n') 
        i = i + (noElPar[8] -1)      

   elif (eleline2[i][3]=='CSRDRIFT'):
        spath = spath + float(eleline2[i][2])    
        mybeamlineX[0] = eleline2[i][0]+ ' '  # Name 
        mybeamlineX[1] = eleline2[i][3]+ ' '  # Type
        mybeamlineX[2] = str(spath)           # s
        mybeamlineX[3] = eleline2[i][2]+ ' '  # Length  
        icountEl = icountEl + 1        
        beamlineFile.write(string.join(mybeamlineX) + ' '+str(icountEl) + ' \n') 
        i = i + (noElPar[1] -1) 
                                
   elif (eleline2[i][3]=='QUAD'):
        spath = spath + float(eleline2[i][2])    
        qlength = float(eleline2[i][2])
        qk1 = float(eleline2[i+1][2])
        k1 = qlength*qk1
        mybeamlineX[0] = eleline2[i][0]+ ' '  # Name 
        mybeamlineX[1] = eleline2[i][3]+ ' '  # Type
        mybeamlineX[2] = str(spath)           # s
        mybeamlineX[3] = eleline2[i][2]+ ' '  # Length  
        mybeamlineX[7] = str(k1) + ' '        # K1 
        icountEl = icountEl + 1            
        beamlineFile.write(string.join(mybeamlineX) + ' '+str(icountEl) + ' \n')        
        i = i + (noElPar[2] -1)  

   elif (eleline2[i][3]=='CSRCSBEND'):
        spath = spath + float(eleline2[i][2])    
        k0L = float(eleline2[i+1][2])
        mybeamlineX[0] = eleline2[i][0]+ ' '      # Name 
        mybeamlineX[1] = eleline2[i][3]+ ' '      # Type
        mybeamlineX[2] = str(spath)               # s
        mybeamlineX[3] = eleline2[i][2]+ ' '      # Length  
        mybeamlineX[4] = str(k0L) + ' '           # K0L (=angle)
        mybeamlineX[5] = eleline2[i+10][2]+ ' '   # E1  
        mybeamlineX[6] = eleline2[i+11][2]+ ' '   # E2          
        mybeamlineX[9] = eleline2[i+12][2]+ ' '   # TILT        
        mybeamlineX[10] = eleline2[i+15][2]+ ' '  # HGAP  (half gap)
        mybeamlineX[11] = eleline2[i+16][2]+ ' '  # FINT  (normalized field integral)     
        icountEl = icountEl + 1             
        beamlineFile.write(string.join(mybeamlineX) + ' '+str(icountEl) + ' \n')        
        sbendFile.write(string.join(mybeamlineX) + '\n')
        icountSbend = icountSbend + 1
        i = i + (noElPar[3] -1)  

   elif (eleline2[i][3]=='RFCW'):
        spath = spath + float(eleline2[i][2])    
        eleRfPhase = float(eleline2[i+3][2])
        litrackPhase = eleRfPhase -90.e0 
        mybeamlineX[0] = eleline2[i][0]+ ' '        # Name 
        mybeamlineX[1] = eleline2[i][3]+ ' '        # Type
        mybeamlineX[2] = str(spath)                 # s
        mybeamlineX[3] = eleline2[i][2]+ ' '        # Length  
        mybeamlineX[12] = eleline2[i+2][2]+ ' '     # Voltage (MV)
        mybeamlineX[13] = str(litrackPhase)+ ' '    # Phase (deg)
        mybeamlineX[14] = eleline2[i+4][2]+ ' '     # Frequency (Hz)
        icountEl = icountEl + 1        
        cavFile.write(str(icountCav+1)+' '+string.join(mybeamlineX)  + '\n')
        beamlineFile.write(string.join(mybeamlineX) + ' '+str(icountEl) + ' \n') 
        icountCav = icountCav + 1
        i = i + (noElPar[5] -1)  

   elif  (eleline2[i][3]=='MATR'):            
        spath = spath + float(eleline2[i][2])    
        mybeamlineX[0] = eleline2[i][0]+ ' '  # Name 
        mybeamlineX[1] = eleline2[i][3]+ ' '  # Type
        mybeamlineX[2] = str(spath)           # s        
        mybeamlineX[3] = eleline2[i][2]+ ' '  # Length      
        icountEl = icountEl + 1                      
        beamlineFile.write(string.join(mybeamlineX) + ' '+str(icountEl) + ' \n') 
        i = i + (noElPar[9] -1)             

   elif (eleline2[i][3]=='RFDF'):
        spath = spath + float(eleline2[i][2])
        eleRfPhase = float(eleline2[i+1][2])
        litrackPhase = eleRfPhase -90.e0
        mybeamlineX[0] = eleline2[i][0]+ ' '  # Name
        mybeamlineX[1] = eleline2[i][3]+ ' '  # Type
        mybeamlineX[2] = str(spath)           # s
        mybeamlineX[3] = eleline2[i][2]+ ' '  # Length
        mybeamlineX[12] = eleline2[i+4][2]+ ' '     # Voltage (MV)
        mybeamlineX[13] = str(litrackPhase)+ ' '    # Phase (deg)
        mybeamlineX[14] = eleline2[i+3][2]+ ' '     # Frequency (Hz)
        icountEl = icountEl + 1
        beamlineFile.write(string.join(mybeamlineX) + ' '+str(icountEl) + ' \n')
        i = i + (noElPar[10] -1)

   elif  (eleline2[i][3]=='MONI'):
        spath = spath + float(eleline2[i][2])
        mybeamlineX[0] = eleline2[i][0]+ ' '  # Name
        mybeamlineX[1] = eleline2[i][3]+ ' '  # Type
        mybeamlineX[2] = str(spath)           # s
        mybeamlineX[3] = eleline2[i][2]+ ' '  # Length
        icountEl = icountEl + 1
        beamlineFile.write(string.join(mybeamlineX) + ' '+str(icountEl) + ' \n')
        i = i + (noElPar[11] -1)
   
   i = i + 1 

cavFile.close()
sbendFile.close()
beamlineFile.close()

#----------------------------
# Find voltages contributed by  each linac section
# The user  should open the file cavFile.txt and identfy the 
# LAST cavity for each linac section L1, HL, L2, L3 and 
# define the vector beginLS
noSections=4
beginLS = range(0,noSections-1)
beginLS = [16-1, 32-1,128-1,288-1]
volL = range(0,noSections)

cavFile = open('cavList.txt', 'r') 

i=0
vL = 0
dE = 0
while (i <= beginLS[0]):  # L1
        myline = cavFile.readline()
        if not myline: break
        myline2 = string.split(myline)
#        print ' myline2[13]=', myline2
#        print ' myline2[13]=', float(myline2[13])
        vcavStr = myline2[13]
        phscavStr = myline2[14]
        vcav = float(vcavStr)
        phscav = float(phscavStr)  # phase in deg (Litrack convention)
        
        vL = vL + vcav
        dE = dE + vcav*math.cos(phscav*3.141592653589793/180)
        i = i + 1 

volL[0] = vL        
print ' L1 dE (MeV) =', dE*1.e-6
print ' L1  V  (MV) =', volL[0]*1.e-6


vL = 0
dE = 0
while (i  <=  beginLS[1] ):  # HL
        myline = cavFile.readline()
        if not myline: break
        myline2 = string.split(myline)
#        print ' myline[12]=', myline2
#        print ' myline[12]=', float(myline2[12])
        vcavStr = myline2[13]
        phscavStr = myline2[14]
        vcav = float(vcavStr)
        phscav = float(phscavStr)  # phase in deg (Litrack convention)
        
        vL = vL + vcav
        dE = dE + vcav*math.cos(phscav*3.141592653589793/180)
        i = i + 1 
volL[1] = vL        
print ' HL dE (MeV) =', dE*1.e-6
print ' HL  V  (MV) =', volL[1]*1.e-6


vL = 0
dE = 0
while (i  <=  beginLS[2] ):  # L2
        myline = cavFile.readline()
        if not myline: break
        myline2 = string.split(myline)
        vcavStr = myline2[13]
        phscavStr = myline2[14]
        vcav = float(vcavStr)
        phscav = float(phscavStr)  # phase in deg (Litrack convention)
        
        vL = vL + vcav
        dE = dE + vcav*math.cos(phscav*3.141592653589793/180)
        i = i + 1 
volL[2] = vL        
print ' L2 dE (MeV) =', dE*1.e-6
print ' L2  V  (MV) =', volL[2]*1.e-6

vL = 0
dE = 0
while (i  <=  beginLS[3] ):  # L3
        myline = cavFile.readline()
        if not myline: break
        myline2 = string.split(myline)
        vcavStr = myline2[13]
        phscavStr = myline2[14]
        vcav = float(vcavStr)
        phscav = float(phscavStr)  # phase in deg (Litrack convention)
        
        vL = vL + vcav
        dE = dE + vcav*math.cos(phscav*3.141592653589793/180)
        i = i + 1 
volL[3] = vL        
print ' L3 dE (MeV) =', dE*1.e-6
print ' L3  V  (MV) =', volL[3]*1.e-6
