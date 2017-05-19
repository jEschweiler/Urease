import numpy as np
from ctypes import cdll, c_int, c_float, byref 
import os
class impact(object):
        
    ## load IMPACT library
    # @param library path 
    def __init__(self, libfile="/home/joeesch/anaconda3/lib/python3.5/libimpact.so"):

        try:
            self.libs=cdll.LoadLibrary(libfile)        
            self.libs.pa2tjm.restype = c_float

        except:
            raise Exception("loading library %s failed!"%libfile)

        #declare output variables
        self.ccs=c_float()
        self.sem=c_float()
        self.niter= c_int()
        

    ## compute CCS using the PA method as implemented in IMPACT library.
    # @param points xyz coordinates of atoms, Angstrom (nx3 numpy array)
    # @param radii van der Waals radii associated to every point (numpy array with n elements)
    # @param a power-law factor for calibration with TJM
    # @param b power-law exponent for calibration with TJM
    # @retval ccs_tjm CCS
    # @retva errl sem standardor
    # @retval niter number of iterations
    def get_ccs(self, points, radii, a=0.842611, b=1.051280):
        
        #create ctypes for intput data
        unravel=np.ravel(points)
        cpoints=(c_float*len(unravel))(*unravel)
        cradii=(c_float*len(radii))(*radii)
        natoms=(c_int)(len(radii))
        
        #call library, and rescale obtained value using exponential law
        self.libs.ccs_from_atoms_defaults(natoms, byref(cpoints), byref(cradii), byref(self.ccs), byref(self.sem), byref(self.niter))
        PA = self.ccs
        #ccs_tjm=self.libs.pa2tjm(c_float(a), c_float(b), self.ccs)
        
        return PA.value#self.sem.value, #self.niter.value
    



im = impact()


files=os.listdir(".")
#fnames=[x for x in files if x.endswith('_coordfile_new.npy')] #Get the files
#refc = np.load("05_2qi9_refcoords.npy")
fnames =["coordfile_new.npy"]
good_models = []
i = 0
g = 0
l = 0
we = 0

ABC = np.array([0,1,2,3,8,9,10,11,16,17,18,19])
ABCMDFG1 = np.array([0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19])
ABCMD = np.array([0,1,2,3,4,5,8,9,10,11,16,17,18,19])
ABCMDFG2 = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])

onegood = 0
twogood = 0
threegood = 0
for f in fnames:
	good_models = []
	#refc = np.load("refcoords.npy") 
	#rc = np.delete(refc,1,0)
	#rc1 = np.delete(rc, 1,0)
	#refc = rc1
	#print(np.shape(refc))
	radii = np.load("radii.npy")
	#print(radii[ABCMD])
	#radii1 = [radii[0],radii[1],radii[2],radii[3],radii[4],radii[5],radii[6],radii[7],radii[8],radii[9],radii[14],radii[15],radii[16]]
	#radii2 = [radii[0],radii[1],radii[2],radii[3],radii[4],radii[5],radii[6],radii[7],radii[8],radii[9],radii[10],radii[11],radii[12],radii[13],radii[14],radii[15],radii[16]]
	#rr = np.delete(radii, 1, 0)
	#rr1 = np.delete(rr, 1, 0)
	#radii = rr1
	#refccs = 2400
	ccslist1 = []
	ccslist2 = []
	ccslist3 = []
	eccslist1 = []
	eccslist2 = []
	eccslist3 = []
	arr = np.load(f)

	goodmodels3 = []
	goodmodels2 = []
	goodmodels1 = []
	goodmodels0 = []
	#print(np.shape(arr))
	Alist = []
	AClist = []
	print(np.shape(arr))
	length = np.shape(arr)[0]
	i = 0
	for struct in arr:
		if i%100 == 0:
			print(i/length*100)
		i += 1
		ccsABCMDFG3 = (im.get_ccs(struct,radii))
		ABCMDFG3er = ((23417 - ccsABCMDFG3)/23417)
		if  abs(ABCMDFG3er) < 0.03:
			goodmodels3.append(struct)
			ccsABCMDFG2 = im.get_ccs(struct[ABCMDFG2], radii[ABCMDFG2])
			ABCMDFG2er = ((19646 - ccsABCMDFG2)/19646)
			if abs(ABCMDFG2er) < 0.03:
				goodmodels2.append(struct)	
				#ccsABCMDFG = im.get_ccs(struct[ABCMDFG1], radii[ABCMDFG1])
				#ABCMDFGer = ((15899 - ccsABCMDFG)/15899)
				#if abs(ABCMDFGer)  < 0.05:
				#goodmodels1.append(struct)	
				ccsABCMD = im.get_ccs(struct[ABCMD], radii[ABCMD])
				ABCMDer = ((12879 - ccsABCMD)/12879)
				if abs(ABCMDer) < 0.03:
					goodmodels0.append(struct)
					
						
					
					
				
			
		#if abs(ABCMDFGer) < 0.05 and abs(ABCMDFG2er) < 0.05 and abs(ABCMDFG3er) < 0.05: good_models.append(struct); we +=1;		


		#if abs(er) < 0.03: good_models.append(struct); we+= 1
		#if er > 0: g += 1
		#if er < 0: l +=1
		#ccslist1.append(ccsABCMD)
		#ccslist2.append(ccsABCMDFG2)
		#ccslist3.append(ccsABCMDFG3)

		#eccslist1.append(ABCMD)#er)
		#eccslist2.append(ABCMDFG2er)
		#eccslist3.append(ABCMDFG3er)
		#Alist.append(Aer)
		#AClist.append(ACer)	
	#ccslist = np.array(ccslist)
	#good_models0 = np.array(good_models)
	#print(i)
	print(np.shape(goodmodels0))
	#print(np.shape(goodmodels1))
	print(np.shape(goodmodels2))
	print(np.shape(goodmodels3))
	#print(np.mean(eccslist2), np.std(Alist))
	#print(refccs)
	#print(radii)
	#print(im.get_ccs(refc, radii))




	#np.save(f[:10]+"_ccslist.npy", ccslist)
	np.save("goodmodels0.npy", goodmodels0)
	#np.save("goodmodels1.npy", goodmodels1)
	np.save("goodmodels2.npy", goodmodels2)
	np.save("goodmodels3.npy", goodmodels3)
	#for li in [ccslist, ccslistAAA,ccslistCCC,ccslistAC3]:
	#	print(np.mean(li), np.std(li))
		
