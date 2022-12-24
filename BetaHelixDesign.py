### DESCRIPTION ###
# Python script used for remodeling and design of 
#  tau beta-helix protein scaffold

### Import Rosetta/PyRosetta
import pyrosetta
import pyrosetta.toolbox

__author__ = "Aina Adekunle"
__copyright__ = "Copyright 2022, Aina Adekunle"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Aina Adekunle"
__email__ = "aina@phas.ubc.ca"
__status__ = "Production"


### Initialize PyRosetta
pyrosetta.init("-ignore_unrecognized_res 1 -ex1 -ex2aro -detect_disulf 0")


### REMODELING STAGE ###

### Clean the AphaFold model pdb file
pyrosetta.toolbox.cleanATOM("bhp_model.pdb")

### Import the structure into rosetta space
start_pose = pyrosetta.pose_from_pdb("bhp_model.clean.pdb")
pose = start_pose.clone()

#seq=pose.sequence()
#print(seq)


### Create the constraint file
cst_text ='''CST::BEGIN
  TEMPLATE::   ATOM_MAP: 1 atom_type: CAbb
  TEMPLATE::   ATOM_MAP: 1 is_backbone
  TEMPLATE::   ATOM_MAP: 1 residue1: ARNDCQEGHILKMFSTWYVP
  TEMPLATE::   ATOM_MAP: 2 atom_type: CAbb
  TEMPLATE::   ATOM_MAP: 2 is_backbone
  TEMPLATE::   ATOM_MAP: 2 residue1: ARNDCQEGHILKMFSTWYVP

  CONSTRAINT:: distanceAB:    3.50   0.20 100.00  0
CST::END'''

cstfile = "cstfile"
with open(cstfile, "w") as f:
    f.write(cst_text)

#!cat {cstfile}


# create a blueprint file
bpfile = "./bpfile"
#rm_residues = list(range(73,93))
with open(bpfile, "w") as f:
    for i in range(92):
        if i == 0:
            f.write(f'{i+1} {seq[i]} . CST1A\n')

        elif (i >= 73) and (i <= 75):
            f.write(f'{i+1} {seq[i]} L ALLAA\n')

        elif (i >= 76) and (i <= 88):
            f.write(f'{i+1} {seq[i]} H ALLAA\n')

        elif (i >= 89) and (i <= 90):
            f.write(f'{i+1} {seq[i]} H ALLAA\n')

        elif i == 91:
            f.write(f'{i+1} {seq[i]} H ALLAA CST1B\n')

        else:
            f.write(f'{i+1} {seq[i]} .\n')

        
#!cat {bpfile}

### Create a Remodel mover and set options 
rm = pyrosetta.rosetta.protocols.forge.remodel.RemodelMover()

pyrosetta.rosetta.basic.options.set_boolean_option('remodel:design:find_neighbors', True)
pyrosetta.rosetta.basic.options.set_boolean_option('remodel:quick_and_dirty', True)
pyrosetta.rosetta.basic.options.set_file_option('remodel:blueprint', bpfile)
pyrosetta.rosetta.basic.options.set_file_option('enzdes:cstfile', cstfile)

### Register and apply Remodel options
rm.register_options()
rm.apply(pose)

#pose.sequence()

### Save the remodeled structure
pose.dump_pdb("./bhp_remodel.pdb")



### DESIGN STAGE ###

### Clean the remodeled pdb file
pyrosetta.toolbox.cleanATOM("bhp_remodel.pdb")

### Import the structure into rosetta space
start_pose = pyrosetta.pose_from_pdb("bhp_remodel.clean.pdb")
pose = start_pose.clone()

### Set up score function
scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

#print(start_pose.sequence())
#start_pose.dump_pdb("start_pose_bhp.pdb")


### Fast relax the initial strucure
relax = pyrosetta.rosetta.protocols.relax.FastRelax()
relax.set_scorefxn(scorefxn)
relax.apply(start_pose)

start_pose.dump_pdb('bhp_remodel_relaxed.pdb')


### Make a list of which residues to be designed
design_res = [17, 18, 19, 36, 37, 38, 55, 56, 57, 74, 75, 76, 
	      77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92]


### Create a residue file for design specification in chain "A".
resfile = "resfile"
with open(resfile, "w") as f:
    f.write("NATAA\n")
    f.write("start\n")
    for i in design_res:
        f.write("{0} {1} ALLAA\n".format(i, "A"))
        
#!cat {resfile}


### Create a task factory. The task factory accepts all the task operations
tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
# Include standard operations
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
# Include the resfile
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.ReadResfile(resfile))

# Create PackerTask and view it
# packer_task = tf.create_task_and_apply_taskoperations(start_pose)
# print(packer_task)


### Set up a MoveMapFactory to specify which torsions are free to minimize
mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
mmf.all_bb(setting=True)
mmf.all_bondangles(setting=True)
mmf.all_bondlengths(setting=True)
mmf.all_chi(setting=True)
mmf.all_jumps(setting=True)
mmf.set_cartesian(setting=True)


### Set up FastRelax and FastDesign protocols
fr = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1)
fr.cartesian(True)
fr.set_task_factory(tf)
fr.set_movemap_factory(mmf)
fr.min_type("lbfgs_armijo_nonmonotone") 

### Run FastDesign to generate N=100 decoys
job = pyrosetta.PyJobDistributor("decoy", 100, scorefxn, compress=False)
# The above constructs a job distibutor that will create 100 decoys
#  named decoy_1.pdb to filename_N.pdb and a score file, decoy.fasc.
# The PyJobDistributor will not overwrite a file already in existence.
# When initialized, the next available output file is started as an 
#  in-progress file.

job.native_pose = start_pose
# If a native pose is provided, a column of RMSDs will be included in the
#  score file.

while not job.job_complete:
    pose.assign(start_pose)
    fr.apply(pose)
    #%time fr.apply(pose)
    #print(pose.sequence())
    job.output_decoy(pose)

