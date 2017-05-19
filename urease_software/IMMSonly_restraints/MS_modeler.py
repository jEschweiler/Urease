""" 
 Authors: Aaron T. Frank and Joseph Eschweiler
 Monte Carlo simulation script for testing the use of MS data
 Optimization is carried out using a simulated annealing Monte Carlo procedure followed by constant temperature Monte Carlo
"""
import warnings
from itertools import repeat
from optparse import OptionParser
import warnings, sys
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.pmi
import sys
import sampler
import numpy as np
import random
import os

def main():
	# parse command line
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)    
	# MassSpecSystem
	parser.add_option("--distance_force_constant",dest="distance_force_constant",type="float",default=100.0,help="force constant used in distance restraints [default: %default]")    
	parser.add_option("--connectivity_mean",dest="connectivity_mean",type="float",default=1.0,help="the mean value for the connectivity restraints [default: %default]")
	parser.add_option("--max_score",dest="max_score",type="float",default=100.0,help="max score used for connectivity restraints  [default: %default]")
	# MassSpecDynamics
	parser.add_option("--optimization_cycles",dest="optimization_cycles",type="int",default=10,help="number of optimization cycles [default: %default]")
	parser.add_option("--connectivity_force_constant",dest="connectivity_force_constant",type="float",default=0.1,help="force constant used in connectivity restraints [default: %default]")
	parser.add_option("--initial_temperature",dest="initial_temperature",type="int",default=1000,help="initial temperature for simulated annealing MC [default: %default]")
	parser.add_option("--final_temperature",dest="final_temperature",type="float",default=100,help="temperature for constant temperature MC [default: %default]")
	parser.add_option("--mc_cool_cycles",dest="mc_cool_cycles",type="int",default=500,help="number of simulated annealing cooling cycles [default: %default]")
	parser.add_option("--mc_cool_steps",dest="mc_cool_steps",type="int",default=1000,help="number of mc steps per cooling cycle [default: %default]")
	parser.add_option("--mc_cycles",dest="mc_cycles",type="int",default=1000,help="number of constant T mc cycles [default: %default]")    
	parser.add_option("--mc_steps",dest="mc_steps",type="int",default=1000,help="number of constant T mc steps [default: %default]")    
	parser.add_option("--print_annealing",dest="print_annealing",metavar="BOOL",default=False,action="store_true",help="should coordinates be printed out during annealing? [default: %default]")
	parser.add_option("--verbose",dest="verbose",metavar="BOOL",default=False,action="store_true",help="verbose output [default: %default]")
	(options_parse, args) = parser.parse_args()
	nr = []
	f= 0
	if len(args) != 1:
		parser.error("incorrect number of arguments")
	else:        
		# set log level
		IMP.set_deprecation_exceptions(False)
		files=os.listdir(".")
		files=[x for x in files if x.endswith('_restraints.txt')]
		for fname in files:
        # Initialize from the command-line
        #fname=args[0]
        # Create instance of an IMP model object
			m = IMP.Model()
			# Create instance of MassSpecSystem
			s = sampler.MassSpecSystem(m, fname)       
			s.max_score = options_parse.max_score
			s.distance_force_constant = options_parse.distance_force_constant
			s.connectivity_force_constant = options_parse.connectivity_force_constant             
			s.read_restraints()
			print("LINE 115") 
			s.setup_system()
			print("LINE 115") 

			for sr in s.symres:
				m.add_score_state(sr)
			ccslist = []
			for item in s.ccs_radii:
				ccs = item.split(",")[1]
				ccslist.append(ccs)
			radii = np.array(s.radii, dtype = "float")
			refcoords = np.array(s.refcoords)
			refccs = np.array(ccslist,dtype = "float")
			np.save(fname.replace(".pdb1_restraints.txt","") + "_" +"radii.npy", radii)
			np.save(fname.replace(".pdb1_restraints.txt","") + "_" +"refcoords.npy", refcoords)
			np.save(fname.replace(".pdb1_restraints.txt","") + "_" +"refccs.npy", refccs)
			dref = open(fname.replace(".pdb1_restraints.txt","") + "_" + "drestraints.txt", "w")
			dres = []
			#print(radii,refcoords,refccs)
			# loop over distance restraints and run MC optimization
			#for i in range(1):
			#if len(s.distances) > 10:
			#	#print(len(s.distances))
			#	step = int((len(s.distances)-2)/float(7))
			#	print("STEP",step)
			#	#rnd = random.randint(1,len(s.distances)%10)
			#	rnd = 1
			#	ind = [0]
			#	for i in range(8):
			#		#print(i)
			#		ind.append(rnd + i*step)
			#	#print(ind)
			#	dtmp = [s.distances[x] for x in ind]
			#	dtmp.append(s.distances[-1])
			#	s.distances = dtmp

			for dis in enumerate(s.distances):
				dref.write(str(dis[0])+" : "+ "%s\n" % dis[1])
			dref.close()
			for i in range(len(s.distances)): 

				mm = IMP.Model()
				ss = sampler.MassSpecSystem(mm, fname) 
				print("LINE 115")       
				ss.max_score = options_parse.max_score
				ss.distance_force_constant = options_parse.distance_force_constant
				ss.connectivity_force_constant = options_parse.connectivity_force_constant             
				ss.read_restraints()

				ss.setup_system()

				for sr in ss.symres:
					m.add_score_state(sr)



				print(s.fname + "RESTRAINT SET" + str(i))
				#dref.write("%s\n" % s.distances[i])
				nres = len(s.distance_restraints)
				nr.append(nres)      

				ss.setup_restraints(d_res = i)
				print(ss.restraints)

				# Create instance of IMP RestraintsScoringFunction
				sf = IMP.core.RestraintsScoringFunction(ss.restraints, "scoring function")
				#print(s.restraints)
				# Create instance of MassSpecDynamics
				o = sampler.MassSpecDynamics(ss, sf)        
				o.initial_temperature = options_parse.initial_temperature
				o.final_temperature = options_parse.final_temperature
				o.mc_cool_cycles = options_parse.mc_cool_cycles
				o.mc_cool_steps = options_parse.mc_cool_steps
				o.mc_cycles = options_parse.mc_cycles
				o.mc_steps = options_parse.mc_steps
				o.optimization_cycles = options_parse.optimization_cycles
				o.print_annealing = options_parse.print_annealing
				
						# Run optimization
						#if i is 0:
							#o.print_xyzr_css_score_header()
				coords = o.run_MC(i,fname)
				f +=1

			dref.close()
	np.save("nr.npy",nr)
if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
