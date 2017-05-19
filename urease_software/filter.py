import numpy as np
from ctypes import cdll, c_int, c_float, byref 
import os
import scipy.spatial.distance as dis
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
    
import itertools
def pairwise(iterable):
	a,b = itertools.tee(iterable)
	next(b,None)
	return list(zip(a,b))

def isconnected(pd):
	pd = dis.squareform(pd)
	if np.any(np.all(pd > -0, axis = 0)):
		return True
	else: return False
colors = [(1,.4,.4),(.4,.4,1),(.4,1,.4),(1,.4,1),(.4,1,1),(1,.7,.4),(1,.4,.7)] 	
def writepym(coords, r, fname):
	w = open(fname, 'w')    
	w.write('from pymol.cgo import *'+ '\n')
	w.write('from pymol import cmd'+ '\n')
	w.write('from pymol.vfont import plain' + '\n' + 'data={}' + '\n' + "curdata=[]" + '\n') 
	for index, c in enumerate(coords):
		w.write("k='Protein" + str(index) +  " geometry'" +'\n'+ "if not k in data.keys():" +'\n'+"   data[k]=[]"+'\n'+'curdata=['+'\n'+'COLOR,' + str(colors[index][0])+","+str(colors[index][1])+","+ str(colors[index][2])+"," + '\n' + 'SPHERE,'+ str(c[0])+ ','+ str(c[1])+',' + str(c[2])+','+ str(r[index]) +'\n')
		w.write("]"+"\n"+"k='Protein" + str(index) +  " geometry'" + '\n' + "if k in data.keys():" + "\n" + "   data[k]= data[k]+curdata"+'\n'+"else:" +'\n' +"   data[k]= curdata"+"\n")
	w.write("for k in data.keys():" + "\n" + "   cmd.load_cgo(data[k], k, 1)" +"\n"+ "data= {}")
	w.close()

im = impact()


files=os.listdir(".")
fnames=[x for x in files if x.endswith('coordfile.npy')] #Get the files
fnames = ['coordfile.npy']
i = 0
for f in fnames:
	ccslist = []*10000
	arr = np.load(f)
	radii = np.load("radii.npy")
	newcoords = []
	#print(radii)
	#radii = pairwise(radii)
	#radii = dis.squareform(radii)
	radii = list(itertools.combinations(radii, 2))
	r = []
	for i in radii:
		r1 = i[0] + i[1]
		#print(r1)
		r.append(r1)
	p = 0	
	n = 0
	for i in arr:
		pd = dis.pdist(i)
		#pd = dis.squareform(pd)
		pd = (pd-r)/r
		#print(pd)
		
		D = i[5]
		D.shape = (1,3)
		F = i[6]
		F.shape = (1,3)
		G = i[7]
		G.shape = (1,3)
		DF = dis.cdist(D,F, "euclidean")
		FG = dis.cdist(F,G,"euclidean")
		#if DF > 44: DFc += 1
		#if FG > 41: FGc += 1

		if np.any(pd <-0.7) or isconnected(pd):
			#writepym(i, np.load("radii.npy"), str(n)+".pym")
			n += 1
		elif DF < 44 and FG < 41:
			newcoords.append(i)
			p += 1
	print(p)

	np.save("coordfile_new.npy", newcoords)

