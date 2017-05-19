# Urease

IMMS_modeler

Authors: Joseph D. Eschweiler and Aaron T. Frank (joeesch@umich.edu, atfrank@umich.edu) 
 Monte Carlo search for ensembles of structures satisfying IM-MS restraints
Optimization is carried out using a simulated annealing Monte Carlo procedure followed by constant temperature Monte Carlo
Note: This script only generates ensembles of structures; it does not do any further analysis. 
A general workflow for full analysis would be
1.	IMMS_modeler.py to generate an ensemble
2.	filter.py to filter out unphysical models
3.	CCS_calc.py to calculate CCS and filter models based on CCS
4.	get_rmsdmat.py to generate a RMSD matrix for all of the filtered models
5.	cluster.py to perform hierarchical clustering of the filtered models
6.	density.py to visualize clusters (this script will require changes to the source code in most cases)
Inputs: Restraint File
Outputs: coordinate file (numpy array), scores File (numpy array)
Note: I usually only use the scores file to test for convergence – All the scores should be about the same, rendering them pretty meaningless. If you have multimodal populations of scores, it means that some models aren’t reaching a satisfactory minimum. You can change the scoring function or just discard these models
Dependencies: 
Integrated Modeling Platform V2.6 or higher
Numpy
(Must have sampler.py in the directory, or alternatively hard coded into the script)
Usage: 
python IMMS_modeler.py “restraint_file_name”  --options
Options:
--distance_force_constant, type="float", default=100.0,  [force constant used in distance restraints]
--max_score, type="float", default=100.0, [max score allowed for connectivity restraints]
--optimization_cycles, type="int", default=10, [number of optimization cycles]
--connectivity_force_constant, type="float", default=100, [force constant used in connectivity restraints]
--initial_temperature, type="int", default=1000, [initial temperature for simulated annealing MC]
--final_temperature, type="float", default=100, [temperature for constant temperature MC]
--mc_cool_cycles, type="int", default=500, [number of simulated annealing cooling cycles]
--mc_cool_step, type = “int”, default=1000, [number of mc steps per cooling cycle]
--mc_cycles", type="int", default=10000, [number of constant T mc cycles]
--mc_steps, type="int",default=10000, [number of constant T mc steps]
--print_annealing , metavar="BOOL",default=False, action="store_true", [should coordinates be printed out during annealing? ]
--verbose, metavar="BOOL", default=False,action="store_true",[verbose output ?])


Choosing Parameters: 
Each Optimization cycle will output mc_cycles / 100 strctures. The default is 100 structures per optimization cycle. 
Each optimization cycle starts with new randomized coordinates, and may find new local minima to explore. Change the optimization_cycles to generate more sets of 100 related structures. 
The total amount of structures output will be optimization_cycles * mc_cycles / 100. Default = 1000
Default temperature and mc_step parameters should be fine for most applications. For extremely large systems >12 subunits, see the modified code in the forthcoming Urease documents. 
