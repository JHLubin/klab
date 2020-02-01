""" Executed from loop_swap_protocol folder """
import loop_align_updates as la
import pickle
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
    OperateOnResidueSubset, RestrictAbsentCanonicalAASRLT
from pyrosetta.rosetta.protocols.grafting import CCDEndsGraftMover
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

# Load data
with open('protease_db_3.pkl', 'rb') as o:
    db_dict = pickle.load(o)

# Set poses
query_pose = pose_from_pdb('tev.pdb')
subject_name = '1lvm'
subject_pose = pose_from_pdb('aligned_pdbs/{}.pdb'.format(subject_name))
lvm_dat = db_dict['1lvm']

# All anine-swap the start pose
sf = get_fa_scorefxn()
tf = TaskFactory()
restriction = RestrictAbsentCanonicalAASRLT()
restriction.aas_to_keep('A')
selection = ResidueIndexSelector('1-{}'.format(subject_pose.total_residue()))
tf.push_back(OperateOnResidueSubset(restriction, selection))
pt = tf.create_task_and_apply_taskoperations(subject_pose)
prm = PackRotamersMover(sf, pt)
all_A_pose = Pose(subject_pose)
prm.apply(all_A_pose)
all_A_pose.dump_pdb('self_swaps_1lvm_test/1lvm_all_A.pdb')

# Put in loops
for loop_name, loop in lvm_dat.loop_maps.items():
    print(loop_name)

    # Setting up poses
    swapped_pose = Pose(query_pose)
    loop_pose = Pose(all_A_pose, 
        loop.N_outside_overlap_residue.subject_pose_number, 
        loop.C_outside_overlap_residue.subject_pose_number)
    loop_pose.dump_pdb('self_swaps_1lvm_test/loop_{}_loop_only.pdb'.format(loop_name))

    # Setting up CCDEndsGraftMover
    ccdgm = CCDEndsGraftMover()
    ccdgm.set_insert_region(loop.N_splice_residue.query_pose_number,
        loop.C_splice_residue.query_pose_number)
    ccdgm.set_piece(loop_pose, 
        loop.N_overlap_size, loop.C_overlap_size)

    # Applying mover and scoring and dumping the pose
    ccdgm.apply(swapped_pose)
    sf(swapped_pose)
    pp.dump_pdb('self_swaps_1lvm_test/loop_{}_tev_insert.pdb'.format(loop_name))
