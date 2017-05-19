import numpy as np
import scipy.spatial.distance as dist
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.lines as mplines
import scipy.cluster.hierarchy as clust
import os

def kabsch(coord, ref,app):
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
#colors = [(1,.4,.4),(.4,.4,1),(.4,1,.4),(1,.4,1),(.4,1,1),(1,.7,.4),(1,.4,.7)]	

colors = [(0,.6,.6),(1,0,.5),(1,1,.2),(1,1,.2),(.8,.4,0),(.6,1,1),(.8,0,.8),(0,.9,0),(0,.6,.6),(1,0,.5),(1,1,.2),(1,1,.2),(.8,.4,0),(.6,1,1),(.8,0,.8),(0,.9,0),(0,.6,.6),(1,0,.5),(1,1,.2),(1,1,.2),(.8,.4,0),(.6,1,1),(.8,0,.8),(0,.9,0)]	
def writepym(i,coords,radii):
	pymfilename= i + ".pym"
	pymfile=open(pymfilename, "w")
	pymfile.write('from pymol.cgo import *'+ '\n')
	pymfile.write('from pymol import cmd'+ '\n')
	pymfile.write('from pymol.vfont import plain' + '\n' + 'data={}' + '\n' + "curdata=[]" + '\n')
	#print(x for x in enumerate(coords))
	for item in enumerate(coords):
		#print(colors[item[0]][0],colors[item[0]][1], colors[item[0]][2])
		#print(colors[item[0]][0])
		#print(item)
		pymfile.write("k='Protein" + str(item[0]) +  " geometry'" +'\n'+ "if not k in data.keys():" +'\n'+"   data[k]=[]"+'\n'+'curdata=['+'\n'+'COLOR,' + str(colors[item[0]%8][0])+","+str(colors[item[0]%8][1])+","+ str(colors[item[0]%8][2])+"," + '\n' + 'SPHERE,'+ str(item[1][0])+ ','+ str(item[1][1])+',' + str(item[1][2])+','+ str(radii[item[0]]) +'\n')
		pymfile.write("]"+"\n"+"k='Protein" + str(item[0]) +  " geometry'" + '\n' + "if k in data.keys():" + "\n" + "   data[k]= data[k]+curdata"+'\n'+"else:" +'\n' +"   data[k]= curdata"+"\n")

	pymfile.write("for k in data.keys():" + "\n" + "   cmd.load_cgo(data[k], k, 1)" +"\n"+ "data= {}")
	pymfile.close()


files=os.listdir(".")
#refs=[x for x in files if x.endswith('1k90_refcoords.npy')]
np.set_printoptions(threshold=1000000)
pdf = PdfPages("corrected.pdf")

# Read the pairwise distance matrix (discard row and column labels).
#fname = "corrected-res.csv"
distmat = np.load("rmsdmat.npy")

# Calculate the mean of the pairwise similarities.
ii = np.tril_indices(distmat.shape[0], -1)
pwise = distmat[ii]
mdist = np.mean(pwise)
print(mdist)
#print(pwise)

# Generate a historgram of the pairwise similarities.
plt.clf()
plt.hist(pwise, 20, color='lightblue')
plt.xlabel("Similarity")#, size=17)
plt.ylabel("Frequency")#, size=17)
pdf.savefig()


# Do the clustering
h = clust.average(distmat)

# Plot the dendrogram
plt.figure(figsize=(16,10))
plt.figure(linewidth=100.0)
plt.clf()
ax = plt.axes()
for pos in 'right','bottom','top':
    ax.spines[pos].set_color('none')
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('outward', 10))
x=clust.dendrogram(h)

#plt.getp(x)
pdf.savefig()

pdf.close()

#ll = clust.leaves_list(h)
#print(len(ll))
tree = clust.to_tree(h)
#print(tree)
#ctree = clust.cut_tree(h, height = 150)
#print(np.shape(ctree))
ctree = clust.cut_tree(h, n_clusters = 2)
leaves = clust.leaves_list(h)
#print(np.shape(ctree))
ctree = np.reshape(ctree, len(leaves))
#print(np.shape(leaves))
#print(np.shape(ctree))
#print(np.vstack((leaves,ctree)))


files=os.listdir(".")
files=[x for x in files if x.startswith('tetramer_model_')]
print(len(files))

n_clusters = np.max(ctree) + 1
#print(n_clusters)
clusters = [[] for i in range(n_clusters)]

CCC = np.array([2,3,10,11,18,19])
AC3 = np.array([0,2,3,8,10,11,16,18,19])
#MDFG = np.array([4,5,6,7,12,13,14,15,20,21,22,23])
##actually MD
MDFG = np.array([4,5,12,13,20,21])



for i, leaf in enumerate(leaves):
	cluster = ctree[i]
	structure = np.load("goodmodels0.npy")[i]
#	print(len(clusters))
#	print(cluster)
	clusters[cluster].append(structure)


rmsdlist = []
coordlist = []
for clust in clusters:
	l = len(clust)
	av = round(l / 2, -1)
	av = int(av)
	crmsdlist = []
	alignedcoordlist = []
	for o,st in enumerate(clust):
		strmsdlist = []
		stCst = st[CCC]
		stC = stCst - np.mean(stCst, axis = 0)
		st3 = st - np.mean(st, axis = 0)

		#ik = i[np.array([2,7,12])]
		#ikm = ik - np.mean(ik, axis = 0)
		#im = i - np.mean(i, axis = 0)
		#print(i)




		for st2 in clust:
			st2Cst = st2[CCC]
			st2C = st2Cst - np.mean(st2Cst, axis = 0)
			st23 = st2 - np.mean(st2Cst, axis = 0)


			k = kabsch(st2C, stC, st23)
			k = k - np.mean(k, axis =0)
			#r2 = rmsd(k[np.array([3,4,8,9,13,14])], st3[np.array([3,4,8,9,13,14])])
			r = rmsd(k, st3)
			#print(r, r2)
			#r = rmsd(st, k)
			strmsdlist.append(r)
			if o == av:
				alignedcoordlist.append(k)
				#print(r)

			
			#jm = j - np.mean(j, axis = 0)
			#jk = j[np.array([2,7,12])]
			#jkm = jk - np.mean(jk, axis = 0)


			#k = kabsch(jkm, ikm, jm)
			#k = k - np.mean(k, axis =0)
			#r = rmsd(k[np.array([3,4,8,9,13,14])], im[np.array([3,4,8,9,13,14])])
			#r2 = rmsd(k[np.array([2,7,12])], im[np.array([2,7,12])])
			#print(i)
			#print(r, r2)
			#rmsdlist1.append(r)








		crmsdlist.append(strmsdlist)
	#print(alignedcoordlist)	
	rmsdlist.append(crmsdlist)
	coordlist.append(alignedcoordlist)
radii = np.load("radii.npy")
clustcoords = []
for i,item in enumerate(coordlist):
	print(np.shape(item))
	mean = np.mean(item, axis = 0)
	med = round(len(item)/2)
	writepym("cluster_mean_"+str(i), mean, radii)
	#writepym("cluster_med_"+str(i), item[med],radii)
	#print(item))
	np.save("cluster_"+str(i)+".npy", item)
	#print("std ", np.std(item, axis = 0))
	clustcoords.append(mean)
	
np.save("clust_av_coordsn.npy",clustcoords)
	
m = []		
for cl in rmsdlist:
	mean = np.mean(cl)
	m.append(mean)
	print(mean)

print(np.mean(m))
	
