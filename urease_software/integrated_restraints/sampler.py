import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.pmi
import sys
import os.path
import IMP.display
import numpy as np
import math
import matplotlib.pyplot as plt
import time

colors = [(1,.4,.4),(.4,.4,1),(.4,1,.4),(1,.4,1),(.4,1,1),(1,.7,.4),(1,.4,.7)]

class MassSpecSystem:
	def __init__(self,model,fname):
		self.model = model
		self.fname = fname
		self.chains = []
		self.ccs = []
		self.radii = []
		self.ccs_radii = []
		self.distances = []
		self.raw_connect = []
		self.node_labels = []
		self.node_structure = []
		self.composite = []
		self.distance_restraints = []
		self.connectivity_restraint = []
		self.nparticles = 0
		self.nrestraints = 0
		self.size = 0
		self.ds = []
		self.idx = []
		self.rigid_bodies = []
		self.ptypes = []
		self.restraints = []
		self.max_score = 100
		self.distance_force_constant = 10 
		self.connectivity_force_constant = 10
		self.refcoords = []
		self.sympairs = [(0,1),(0,3),(1,3)]
		self.symres = []
	
	def convert_node_names_node_indices(self, child):
		child_index = []
		for i in range(len(child)):
			child_index.append(self.chains.index(child[i]))
		return child_index
		
	def get_node(self, node):
		tmp = node.split("=")[1].split(" ")
		tmp = list(filter(None,tmp))	
		if len(tmp) is not 1:
			#print(list(tmp)[0])
			child1 = tmp[0].replace(","," ").strip().replace("["," ").replace("]"," ").strip().split()
			child1 = self.convert_node_names_node_indices(child1)       
			parent = tmp[1].strip()    
		if len(tmp) is 1:
			#print(tmp[0])
			child1 = tmp[0].replace(","," ").strip().replace("["," ").replace("]"," ").strip().split()
			child1 = self.convert_node_names_node_indices(child1)       
			parent = None
		return [child1, parent]
		
	def parse_raw_connectivity(self):
		node_labels = []
		node_structure = []
		for i in range(len(self.raw_connect)):
			node_labels.append(self.raw_connect[i].split("=")[0].strip())
			node_structure.append(self.get_node(self.raw_connect[i]))
		if (len(node_labels) == len(node_structure)):
			self.node_labels = node_labels
			self.node_structure = node_structure
		else:
			#print("error while reading tree")
			sys.exit(0)
			
	def read_restraints(self): 
		refcoords = []
		if not os.path.exists(self.fname):
			#print("%s does not exist"%(self.fname))
			sys.exit(0)
		else:
			f = open(self.fname, 'r')
			lines = f.readlines()
			counter = 0 
			for line in lines:
				if not line.startswith('##CCS'):
					if not line.startswith('##DISTANCE'):
						if not line.startswith('##CONNECT'):
							if not line.startswith('##REF'):
								if (counter is 0):
									self.ccs_radii.append(line)
								if (counter is 1):
									self.distances.append(line)
								if (counter is 2):
									self.raw_connect.append(line)
								if (counter is 3):
									refcoords.append(line)  
							else:
								counter = 3         
						else:
							counter = 2
					else:
						counter = 1
			f.close() 
                              
       # loop over
			for i in range(len(self.ccs_radii)):
				if i is 0:
					self.composite.append(self.ccs_radii[i].split(",")[0].strip())
					self.composite.append(self.ccs_radii[i].split(",")[1].strip())
					self.composite.append(self.ccs_radii[i].split(",")[2].strip())
				else:
					chain = self.ccs_radii[i].replace('"','').split(",")[0].strip()
					if (len(chain) is 1):
						self.chains.append(chain)
						self.ccs.append(self.ccs_radii[i].split(",")[1].strip())
						self.radii.append(self.ccs_radii[i].split(",")[2].strip())  

			self.nparticles = len(self.chains)
			self.nrestraints = len(self.distances)     
			self.size = float(self.composite[2])
			#print("SELFSIZE",self.size)

       # get connectivity information from restraint file
			self.parse_raw_connectivity()
       #print(self.refcoords)
			for item in refcoords:
				item = item.translate([None, "'[]()\n"])
				item = item.split(":")[1]
				item = item.replace("[","").replace("]","")
				F = np.fromstring(item,sep = " ")
				self.refcoords.append(F)

			
	def setup_system(self):

		self.bb = IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(0,0,0),(IMP.algebra.Vector3D(self.size*2, self.size*2, self.size*2)))
		#print(self.bb)
		self.ps = [IMP.Particle(self.model) for x in range(self.nparticles)]
		self.idx = self.model.get_particle_indexes()
		self.rs = [IMP.pmi.Resolution.setup_particle(self.model, x, 300) for x in self.idx]
		for i in range(self.nparticles):
			init_coors = IMP.algebra.get_random_vector_in(self.bb)
			#print(init_coors)
			self.ds.append(IMP.core.XYZR.setup_particle(self.model, self.idx[i], IMP.algebra.Sphere3D(IMP.algebra.Vector3D(init_coors[0], init_coors[1], init_coors[2]), float(self.radii[i]))))
			#print(self.ds)  
			self.rigid_bodies.append(IMP.core.RigidBody.setup_particle(self.model, self.ds[i],IMP.algebra.ReferenceFrame3D()))
			#print("149")
			self.rigid_bodies[i].set_coordinates(IMP.algebra.Vector3D(init_coors[0], init_coors[1], init_coors[2]))
			#print("151")
			self.rigid_bodies[i].set_coordinates_are_optimized(True)
			#print("153")
			IMP.atom.Mass.setup_particle(self.model, self.idx[i], 1.0)
			#print("155")

	def setup_symmetry_restraint(self):
		for i,sp in enumerate([self.ds[8],self.ds[16]]):
			IMP.core.Reference.setup_particle(sp,self.ds[0])
			tr = IMP.algebra.Transformation3D(
			IMP.algebra.get_rotation_about_axis(IMP.algebra.get_basis_vector_3d(2),np.pi*(2/3) * (i + 1)),IMP.algebra.Vector3D(0, 0, 0))
			sm = IMP.core.TransformationSymmetry(tr)
			c = IMP.core.SingletonConstraint(sm, None,sp)
			self.symres.append(c)
		for j, k in enumerate([self.ds[9],self.ds[17]]):
			IMP.core.Reference.setup_particle(k, self.ds[1])
			tra = IMP.algebra.Transformation3D(
			IMP.algebra.get_rotation_about_axis(IMP.algebra.get_basis_vector_3d(2),np.pi*(2/3) * (j + 1)),IMP.algebra.Vector3D(0, 0, 0))
			sma = IMP.core.TransformationSymmetry(tra)
			d = IMP.core.SingletonConstraint(sma,None,k)
			self.symres.append(d)
			#print("SELF.SYMRES ", self.symres)
		for j, k in enumerate([self.ds[10],self.ds[18]]):
			IMP.core.Reference.setup_particle(k, self.ds[2])
			tra = IMP.algebra.Transformation3D(
			IMP.algebra.get_rotation_about_axis(IMP.algebra.get_basis_vector_3d(2),np.pi*(2/3) * (j + 1)),IMP.algebra.Vector3D(0, 0, 0))
			sma = IMP.core.TransformationSymmetry(tra)
			d = IMP.core.SingletonConstraint(sma,None,k)
			self.symres.append(d)
		for i,sp in enumerate([self.ds[11],self.ds[19]]):
			IMP.core.Reference.setup_particle(sp,self.ds[3])
			tr = IMP.algebra.Transformation3D(
			IMP.algebra.get_rotation_about_axis(IMP.algebra.get_basis_vector_3d(2),np.pi*(2/3) * (i + 1)),IMP.algebra.Vector3D(0, 0, 0))
			sm = IMP.core.TransformationSymmetry(tr)
			c = IMP.core.SingletonConstraint(sm, None,sp)
			self.symres.append(c)
		for j, k in enumerate([self.ds[12],self.ds[20]]):
			IMP.core.Reference.setup_particle(k, self.ds[4])
			tra = IMP.algebra.Transformation3D(
			IMP.algebra.get_rotation_about_axis(IMP.algebra.get_basis_vector_3d(2),np.pi*(2/3) * (j + 1)),IMP.algebra.Vector3D(0, 0, 0))
			sma = IMP.core.TransformationSymmetry(tra)
			d = IMP.core.SingletonConstraint(sma,None,k)
			self.symres.append(d)
			#print("SELF.SYMRES ", self.symres)
		for j, k in enumerate([self.ds[13],self.ds[21]]):
			IMP.core.Reference.setup_particle(k, self.ds[5])
			tra = IMP.algebra.Transformation3D(
			IMP.algebra.get_rotation_about_axis(IMP.algebra.get_basis_vector_3d(2),np.pi*(2/3) * (j + 1)),IMP.algebra.Vector3D(0, 0, 0))
			sma = IMP.core.TransformationSymmetry(tra)
			d = IMP.core.SingletonConstraint(sma,None,k)
			self.symres.append(d)
			#print("SELF.SYMRES ", self.symres)
		for j, k in enumerate([self.ds[14],self.ds[22]]):
			IMP.core.Reference.setup_particle(k, self.ds[6])
			tra = IMP.algebra.Transformation3D(
			IMP.algebra.get_rotation_about_axis(IMP.algebra.get_basis_vector_3d(2),np.pi*(2/3) * (j + 1)),IMP.algebra.Vector3D(0, 0, 0))
			sma = IMP.core.TransformationSymmetry(tra)
			d = IMP.core.SingletonConstraint(sma,None,k)
			self.symres.append(d)
			#print("SELF.SYMRES ", self.symres)
		for j, k in enumerate([self.ds[15],self.ds[23]]):
			IMP.core.Reference.setup_particle(k, self.ds[7])
			tra = IMP.algebra.Transformation3D(
			IMP.algebra.get_rotation_about_axis(IMP.algebra.get_basis_vector_3d(2),np.pi*(2/3) * (j + 1)),IMP.algebra.Vector3D(0, 0, 0))
			sma = IMP.core.TransformationSymmetry(tra)
			d = IMP.core.SingletonConstraint(sma,None,k)
			self.symres.append(d)
			#print("SELF.SYMRES ", self.symres)
#		for j, k in enumerate([self.ds[12],self.ds[19]]):
#			IMP.core.Reference.setup_particle(k, self.ds[5])
#			tra = IMP.algebra.Transformation3D(
#			IMP.algebra.get_rotation_about_axis(IMP.algebra.get_basis_vector_3d(2),np.pi*(2/3) * (j + 1)),IMP.algebra.Vector3D(0, 0, 0))
#			sma = IMP.core.TransformationSymmetry(tra)
#			d = IMP.core.SingletonConstraint(sma,None,k)
#			self.symres.append(d)
#			print("SELF.SYMRES ", self.symres)
#		for j, k in enumerate([self.ds[13],self.ds[20]]):
#			IMP.core.Reference.setup_particle(k, self.ds[6])
#			tra = IMP.algebra.Transformation3D(
#			IMP.algebra.get_rotation_about_axis(IMP.algebra.get_basis_vector_3d(2),np.pi*(2/3) * (j + 1)),IMP.algebra.Vector3D(0, 0, 0))
#			sma = IMP.core.TransformationSymmetry(tra)
#			d = IMP.core.SingletonConstraint(sma,None,k)
#			self.symres.append(d)
#			print("SELF.SYMRES ", self.symres)
			
	def setup_distance_restraints(self, d_res = 0):
		tmp_distance_restraints = self.distances[d_res]
		tmp_distance_restraints = tmp_distance_restraints.translate([None, "'[]()\n"])
		tmp_distance_restraints = tmp_distance_restraints.replace("\n","").split(";")
		tmp_distance_restraints = list(filter(None, tmp_distance_restraints))
		n_distance_restraints = len(tmp_distance_restraints)
		
		for i in range(n_distance_restraints):
			distance_restraint = tmp_distance_restraints[i].replace(",","").replace("[","").replace("]","").replace("'",'').split()
			distance_restraint = list(filter(None, distance_restraint))
			#dref.append(",".join(distance_restraint))
			if distance_restraint:
				pi = int(self.chains.index(distance_restraint[0]))
				pj = int(self.chains.index(distance_restraint[1]))
				rij = float(distance_restraint[2].replace("]",""))
				rij = rij - float(self.radii[pi]) - float(self.radii[pj])
				#print(rij)
				self.distance_restraints.append(IMP.atom.create_distance_restraint(self.rigid_bodies[pi], self.rigid_bodies[pj], rij, self.distance_force_constant))			

	def setup_connectivity_restraints(self):  
     # Create MS connectivity restraint
		#hw = IMP.core.HarmonicWell((-18.25,27.5), self.connectivity_force_constant)
		
		#hw = IMP.core.HarmonicWell((16,58), 100)# range for ALL protein complexes
		hw = IMP.core.HarmonicWell((20,46), 100)
		ss = IMP.core.DistancePairScore(hw)
		self.connectivity_restraint = IMP.core.MSConnectivityRestraint(self.model, ss)
		self.connectivity_restraint.set_maximum_score(0)

     # Connectivity taken from the composite information -- set chain ID -- read in from the restraint file
     #print("COMPOSITE ", self.composite)
		for i in range(len(self.composite[0])):
			index_from_chain = self.chains.index(self.composite[0][i])
			self.ptypes.append(self.connectivity_restraint.add_type([self.ds[index_from_chain]]))
       #print(self.ds[index_from_chain].get_radius())

     # Connectivity taken from the composite information -- set chain ID -- read in from the restraint file
		for i in range(len(self.node_labels)):
			node = self.node_structure[i]
			node_label = self.node_structure[i][0]
			node_parent = self.node_structure[i][1]
			if node_parent is None:
				self.connectivity_restraint.add_composite(node_label)
			else:
				self.connectivity_restraint.add_composite(node_label, self.node_labels.index(node_parent)) 


	def collect_restraints(self):       
		restraints = []
		for i in self.distance_restraints:
			restraints.append(i)
		restraints.append(self.connectivity_restraint)
		self.restraints = restraints

	def setup_restraints(self, d_res = 0):
		#print("SETUP RESTRAINTS") 
		self.distance_restraints = [] 
		self.connectivity_restraint = []      
		self.setup_distance_restraints(d_res = d_res)
		self.setup_symmetry_restraint()
		self.setup_connectivity_restraints()
		self.collect_restraints()
		#print(self.connectivity_restraint.get_connected_pairs())
		#print(self.connectivity_restraint.get_pair_score())


class MassSpecDynamics:
   #'MC sampling options class'
	def __init__(self, system, scoring_function, initial_temperature = 1000, final_temperature = 100, mc_cool_cycles = 500, mc_cool_steps = 5000, mc_cycles = 1000, mc_steps = 1000, optimization_cycles = 10):
		self.system = system
		self.scoring_function = scoring_function
		self.initial_temperature = initial_temperature
		self.final_temperature = final_temperature
		self.mc_cool_cycles = mc_cool_cycles
		self.mc_cool_steps = mc_cool_steps
		self.mc_cycles = mc_cycles
		self.mc_steps = mc_steps
		self.optimization_cycles = optimization_cycles
		self.movers = []
		self.optimizer = IMP.core.MonteCarlo(self.system.model)
		self.print_annealing = True

	def get_coordinates_xyz(self, header, output):
     #""" writes out coordinates in xyz format
     #"""         
		output.write("%s\n" % len(self.system.ds))
		output.write("%s\n" % header)
		for index, particle in enumerate(self.system.ds):
			outi = "C"+str(index)+" "+str(particle.get_coordinates()[0])+" "+str(particle.get_coordinates()[1])+" "+str(particle.get_coordinates()[2])
			output.write("%s\n" % outi)
         
         
	def writepym(self, fname):
		w = open(fname, 'w')    
		w.write('from pymol.cgo import *'+ '\n')
		w.write('from pymol import cmd'+ '\n')
		w.write('from pymol.vfont import plain' + '\n' + 'data={}' + '\n' + "curdata=[]" + '\n') 
		for index, particle in enumerate(self.system.ds):
			w.write("k='Protein" + str(index) +  " geometry'" +'\n'+ "if not k in data.keys():" +'\n'+"   data[k]=[]"+'\n'+'curdata=['+'\n'+'COLOR,' + str(colors[index][0])+","+str(colors[index][1])+","+ str(colors[index][2])+"," + '\n' + 'SPHERE,'+ str(particle.get_coordinates()[0])+ ','+ str(particle.get_coordinates()[1])+',' + str(particle.get_coordinates()[2])+','+ str(particle.get_radius()) +'\n')
			w.write("]"+"\n"+"k='Protein" + str(index) +  " geometry'" + '\n' + "if k in data.keys():" + "\n" + "   data[k]= data[k]+curdata"+'\n'+"else:" +'\n' +"   data[k]= curdata"+"\n")
		w.write("for k in data.keys():" + "\n" + "   cmd.load_cgo(data[k], k, 1)" +"\n"+ "data= {}")
		w.close()
		
	def initialize_MC(self):
     #""" initialize MC optimizer
     #"""         
     # Initialize Monte Carlo sampler     
		self.optimizer.set_return_best(True)
		self.optimizer.set_score_threshold(self.system.max_score*2)
		self.optimizer.set_scoring_function(self.scoring_function)
		self.movers = []
		self.print_annealing = False
 
		# Accumulate movers that are need the Monte Carlo sampler    
		for rbd in self.system.rigid_bodies:
			self.movers.append(IMP.core.RigidBodyMover(rbd, 1, 2))
  
		# Add movers to the Monte Carlo sampler
		self.optimizer.add_movers(self.movers)


			
		
	def run_MC(self, i,name):
		y = 0
     #""" optimizies scoring function using Monte Carlo sampling
     #"""         
     
     # Setup MC optimizer
		
		self.initialize_MC()     
		self.optimizer.set_return_best(False)
		coords = []
		scores = []
		
		for mc in range(self.optimization_cycles): ## RANDOMIZE COORDINATES
			for particle in self.system.ds:
				init_coors = IMP.algebra.get_random_vector_in(self.system.bb)
				particle.set_coordinates(init_coors)
				
			if y != 0:
				for i in [0,1,2,3,5,8,9,10,11,13,16,17,18,19,21]:
				#self.system.ds[i].set_coordinates_are_optimized(False)
					self.optimizer.add_mover(self.movers[i])
			cycle = "annealing"
			T = self.initial_temperature
			print("INIT TEMP", T)
			start = time.time() 
			score = []  
			temp = []             	 
			#for cc in range(self.mc_cool_cycles*7):
			print(self.optimizer.get_movers())
			for cc in range(self.mc_cool_cycles*10):	
			#for cc in range(500):   ## DO THE ANNEALING MC
				T = 0.999*T
				self.optimizer.set_kt(T)
				self.optimizer.optimize(self.mc_cool_steps)
				score.append(self.scoring_function.evaluate(False))
				temp.append(T)
			#plt.clf()
			#plt.plot(temp, score)
			#plt.ylim([0,100000])
			#plt.show()
			print("FINAL TEMP,", T)
			T = 400
			stop = time.time()
			print("annealing time", stop - start)
			print(y)
			
			
			
			#for particle in self.system.ds[np.array([0,1,2,3])]:
			for i in [0,1,2,3,5,8,9,10,11,13,16,17,18,19,21]:
				#self.system.ds[i].set_coordinates_are_optimized(False)
				self.optimizer.remove_mover(self.movers[i])
				#init_coors = IMP.algebra.get_random_vector_in(self.system.bb)
				#particle.set_coordinates_are_optimized(False)




			print(self.optimizer.get_movers())	
			#for c in range(self.mc_cycles): ## DO THE CONSTANT TEMP MC
			start = time.time()
			score = []
			step = []
			for c in range(1000): ## DO THE CONSTANT TEMP MC
				self.optimizer.set_kt(T)
				#self.optimizer.optimize(self.mc_steps)
				self.optimizer.optimize(1000000)
				if c%100 == 0:   ## WRITE LIST TO AN ARRAY EVERY 100TH CYCLE
					score.append(self.scoring_function.evaluate(False))
					step.append(c)
					y+=1
					list = []
					for i,particle in enumerate(self.system.ds):
						x = particle.get_coordinates()
						list.append(x)
						x = np.array(list)
					if type(coords) is np.ndarray:
						coords = coords.tolist()
					if type(scores) is np.ndarray:
						scores = scores.tolist()
					coords.append(list)
					#print(list)
					scores.append(self.scoring_function.evaluate(False))

			stop = time.time()
			print("MC Time", stop - start)


		coords = np.array(coords)
		scores = np.array(scores)
		#np.save(name.replace(".pdb1_restraints.txt","")+ "_" + str("%02d" % i) + "_" +"coordfile.npy",coords)
		np.save("coordfile.npy",coords) ## ALL THE STRUCTURES
		np.save(name.replace(".pdb1_restraints.txt","")+ "_" + str("%02d" % i) + "_" +"scorefile.npy",scores)
