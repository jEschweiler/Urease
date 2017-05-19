import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math
import scipy.spatial.distance as dist

colors = [(1,.4,.4),(.4,.4,1),(.4,1,.4),(1,.4,1),(.4,1,1),(1,.7,.4),(1,.4,.7)]
def writepym(i,coords,radii):
	pymfilename= i + ".pym"
	pymfile=open(pymfilename, "w")
	pymfile.write('from pymol.cgo import *'+ '\n')
	pymfile.write('from pymol import cmd'+ '\n')
	pymfile.write('from pymol.vfont import plain' + '\n' + 'data={}' + '\n' + "curdata=[]" + '\n')
	#print(x for x in enumerate(coords))
	for item in enumerate(coords):
		#print(colors[item[0]][0])
		#print(item)
		pymfile.write("k='Protein" + str(item[0]) +  " geometry'" +'\n'+ "if not k in data.keys():" +'\n'+"   data[k]=[]"+'\n'+'curdata=['+'\n'+'COLOR,' + str(colors[item[0]%7][0])+","+str(colors[item[0]%7][1])+","+ str(colors[item[0]%7][2])+"," + '\n' + 'SPHERE,'+ str(item[1][0])+ ','+ str(item[1][1])+',' + str(item[1][2])+','+ str(radii[item[0]]) +'\n')
		pymfile.write("]"+"\n"+"k='Protein" + str(item[0]) +  " geometry'" + '\n' + "if k in data.keys():" + "\n" + "   data[k]= data[k]+curdata"+'\n'+"else:" +'\n' +"   data[k]= curdata"+"\n")

	pymfile.write("for k in data.keys():" + "\n" + "   cmd.load_cgo(data[k], k, 1)" +"\n"+ "data= {}")
	pymfile.close()
	#print("wrote ", pymfilename)


def get_ccs(coords,radii):
			return impact.PA(im.get_ccs(coords, radii))
			
def align(coords):
	coords= coords -  np.mean(coords, axis = 0)

def kabsch(coord, ref, app):
	C = np.dot(np.transpose(coord), ref)
	V, S, W = np.linalg.svd(C)
	#print("VSW", V,S,W)
	d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
	if d:
		S[-1] = -S[-1]
		V[:,-1] = -V[:,-1]

	# Create Rotation matrix U
	U = np.dot(V, W)

	# Rotate coord
	kcoord = np.dot(app, U)
	
	return kcoord

def rmsd(coord, ref):
	sd = (coord - ref)**2
	ssd = np.mean(sd)
	rmsd = np.sqrt(ssd)
	return rmsd
	
def angle(n1, n2,n3):
	d1 = np.linalg.norm(n1-n2)
	d2 = np.linalg.norm(n1-n3)
	d3 = np.linalg.norm(n2-n3)
	AB = n2 - n1
	BC = n3 - n2
	angle = math.acos(np.dot(AB,BC)/(d1*d3))
	return math.degrees(angle)

def translatei(ecoords):
	t = len(ecoords)

errors = []

coords = np.load("goodmodels0.npy")
#index = np.random.choice(len(coords), 1000, replace = False)
#print(index)
#coords = coords[index]
#np.save("coords_random.npy", coords)
radii = np.load("radii.npy")
CCC = np.array([2,3,10,11,18,19])
AC3 = np.array([0,2,3,8,10,11,16,18,19])

##actually MD
MDFG = np.array([4,5,12,13,20,21])
print(np.shape(coords))
rmsdlist = []
radii = (np.load("radii.npy"))
l = 0
P = 0
length = len(coords)
z = 0
for i in coords:
	if z%100 == 0:
		print(z/length*100)
	z += 1
	iCi = i[CCC]
	Ci = iCi - np.mean(iCi, axis = 0)
	im = i - np.mean(iCi, axis = 0)
	#ACi = np.array([i[0],i[2],i[5],i[7],i[10],i[12]])
	
	rmsdlist1 = []
	#refd = dist.pdist(i, "euclidean")
	for j in coords:
		jCj = j[CCC]
		Cj = jCj - np.mean(jCj, axis = 0)
		jm = j - np.mean(jCj, axis = 0)
		k = kabsch(Cj, Ci, jm)
		kAC = k[AC3]
		#kMDFG = k[MDFG]
		r = rmsd(kAC, im[AC3])
		#r2 = rmsd(kMDFG, im[MDFG])
		
		#r = rmsd(k,i)
		#print(r, r2)
		#d = dist.pdist(j, "euclidean")
		#print(np.shape(d)[0])
		#rmsdd = (2/(np.shape(d)[0]**2 - np.shape(d)[0])*np.sum((d - refd)**2))**.5
		rmsdlist1.append(r)
	rmsdlist.append(rmsdlist1)
	
print(np.shape(rmsdlist))
#print(radii[np.array([2,9,16])])
#print(radii[np.array([0,2,5,7,10,12])])
np.save("rmsdmat.npy",rmsdlist)

#plt.show()


			
