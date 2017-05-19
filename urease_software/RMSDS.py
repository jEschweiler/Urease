import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math
import scipy.spatial.distance as dist

colors =[(0,.6,.6),(1,0,.5),(1,1,.2),(1,1,.2),(.6,1,1),(.8,0,.8),(0,.9,0),(0,.6,.6),(1,0,.5),(1,1,.2),(1,1,.2),(.6,1,1),(.8,0,.8),(0,.9,0),(0,.6,.6),(1,0,.5),(1,1,.2),(1,1,.2),(.6,1,1),(.8,0,.8),(0,.9,0)]	
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
#amb0 = np.load("goodmodels_amb_0.npy")
#amb0 = np.mean(amb0, axis = 0)
#amb1 = np.load("goodmodels_amb_1.npy")
#amb1 = np.mean(amb1, axis = 0)
ref = np.load("refcoords3.npy")
refT = np.load("refcoordsT.npy")
#ref = refT

#coords = np.load("goodmodels0.npy")
radii = np.load("radii.npy")
#CCCmodel = np.array([2,3,10,11,18,19])
CCC = np.array([2,3,9,10,16,17])
#AC3 = np.array([0,2,3,8,10,11,16,18,19])
AC3 = np.array([0,2,3,7,9,10,14,16,17])
radii = np.load("radii.npy")

#amb0str = np.delete(amb0, np.array([4,12,20]), axis = 0)
#amb1str = np.delete(amb1, np.array([4,12,20]), axis = 0)
radiistr = np.delete(radii, np.array([4,12,20]), axis = 0)



#amb0centered = amb0str - np.mean(amb0str[CCCref], axis = 0)
#amb0CCCcentered = amb0str[CCCref] - np.mean(amb0str[CCCref], axis = 0)

#amb1centered = amb1str - np.mean(amb1str[CCCref], axis = 0)
#amb1CCCcentered = amb1str[CCCref] - np.mean(amb1str[CCCref], axis = 0)

centered_ref = ref - np.mean(ref[CCC], axis = 0)
#centered_ref = ref - np.mean(ref[AC3], axis = 0)
#refCCCcentered = ref[CCCref] - np.mean(ref[CCCref], axis = 0)

#k0 = kabsch(amb0CCCcentered, refCCCcentered, amb0centered)
#k1 = kabsch(amb1CCCcentered, refCCCcentered, amb1centered)

#k0centered = k0 - np.mean(k0, axis = 0)
#k1centered = k1 - np.mean(k1, axis = 0)

#rmsd0 = rmsd(k0, refcentered)
#r#msd1 = rmsd(k1, refcentered)
#import scipy.spatial.distance as dist
#dm0 = dist.pdist(k0)
#dm1 = dist.pdist(k1)
#dmR = dist.pdist(refcentered)

#drmsd0 = rmsd(dm0,dmR)
#drmsd1 = rmsd(dm1,dmR)

#coords = np.concatenate((k1, refcentered), axis = 0)
#radii = np.concatenate((radiistr, radiistr))
#writepym("kabsch1", coords, radii)
#print(rmsd0, rmsd1)
#print(drmsd0, drmsd1)

import matplotlib.pyplot as plt


clusters = ["cluster_0.npy", "cluster_1.npy", "cluster_2.npy","cluster_3.npy"]

for cluster in clusters:
	clust = np.load(cluster)
	cluster_rmsds = []
	i = 0
	for structure in clust:
		stripped_structure = np.delete(structure, np.array([4,12,20]), axis = 0)
		centered_structure = stripped_structure - np.mean(stripped_structure[CCC], axis = 0)
		#print(np.mean(centered_structure[CCC], axis = 0),np.mean(centered_ref[CCC], axis = 0))
		k = kabsch(centered_structure[CCC], centered_ref[CCC], centered_structure)
		#print(np.mean(k[AC3], axis = 0), np.mean(centered_ref[AC3], axis = 0))
		#k1 = k - np.mean(k[AC3], axis = 0)
		#ref1 = centered_ref - np.mean(centered_ref[AC3], axis = 0)
		rms = rmsd(k, centered_ref)
		cluster_rmsds.append(rms)
		if i == 0:
			coords = np.concatenate((k, centered_ref), axis = 0)
			radii = np.concatenate((radiistr, radiistr))
			writepym("kabsch1", coords, radii)
		i += 1
	plt.clf()
	plt.hist(cluster_rmsds)
	plt.show()


	
		
		
		
		



##actually MD
MDFG = np.array([4,5,12,13,20,21])
#print(np.shape(coords))
rmsdlist = []
radii = (np.load("radii.npy"))
l = 0
P = 0
#for i in coords:
#	iCi = i[CCC]
#	Ci = iCi - np.mean(iCi, axis = 0)
#	im = i - np.mean(iCi, axis = 0)
#	#ACi = np.array([i[0],i[2],i[5],i[7],i[10],i[12]])
#	
#	rmsdlist1 = []
#	#refd = dist.pdist(i, "euclidean")
#	for j in coords:
#		jCj = j[CCC]
#		Cj = jCj - np.mean(jCj, axis = 0)
#		jm = j - np.mean(jCj, axis = 0)
#		k = kabsch(Cj, Ci, jm)
#		kAC = k[AC3]
#		kMDFG = k[MDFG]
#		r = rmsd(kAC, im[AC3])
#		r2 = rmsd(kMDFG, im[MDFG])
#		if l == 0:
#			if r < 1: writepym(str(1),k,radii); l +=1
#		if P == 0:
#			if r > 13: writepym(str(2),k,radii); P += 1
#		
#		#r = rmsd(k,i)
#		print(r, r2)
#		#d = dist.pdist(j, "euclidean")
#		#print(np.shape(d)[0])
#		#rmsdd = (2/(np.shape(d)[0]**2 - np.shape(d)[0])*np.sum((d - refd)**2))**.5
#		rmsdlist1.append(r)
#	rmsdlist.append(rmsdlist1)
#	
#print(np.shape(rmsdlist))
##print(radii[np.array([2,9,16])])
##print(radii[np.array([0,2,5,7,10,12])])
#np.save("rmsdmat.npy",rmsdlist)

##plt.show()


			
