from pyrosetta import *
from pyrosetta.rosetta.protocols.grafting import CCDEndsGraftMover
with open('changeable_loops_v2.csv', 'r') as r:
looplist = r.readlines()

init()
tev = pose_from_pdb('tev.pdb')
sf = get_fa_scorefxn()

failed_loops = []
for loop in looplist[1:]:
    try:
        # Read in database output
        loop_specs = loop.split(',')
        site = int(loop_specs[0])
        source = loop_specs[2]
        query_start = int(loop_specs[7])
        query_end = int(loop_specs[8])
        subject_start = int(loop_specs[9]) - 2 # -2 is Correcting a bug
        subject_end = int(loop_specs[10]) - 2 # -2 is Correcting a bug
        n_overlap = int(loop_specs[11])
        c_overlap = int(loop_specs[12])

        # Load pose, make subpose including only the loop
        pose = pose_from_pdb('aligned_pdbs/{}.pdb'.format(source))
        # loop_subpose = Pose(pose, subject_start - n_overlap, subject_end + c_overlap)
        loop_subpose = Pose(pose, subject_start, subject_end) # Not using overlaps

        # Setting up CCDEndsGraftMover
        ccdgm = CCDEndsGraftMover()
        ccdgm.set_insert_region(query_start, query_end + 1)
        # ccdgm.set_piece(loop_subpose, n_overlap, c_overlap)
        ccdgm.set_piece(loop_subpose, 0, 0) # Not using overlaps

        # Applying mover and scoring and dumping the pose
        pp = Pose(tev)
        ccdgm.apply(pp)
        loop_name = 'exchanged_loops/loop{}_{}.pdb'.format(site, source)
        sf(pose)
        pp.dump_pdb(loop_name)
    except:
        failed_loops.append(loop)

print(*failed_loops)