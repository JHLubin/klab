#!/usr/bin/env python

from pyrosetta import * ; init( extra_options='-ignore_zero_occupancy false' )# -extra_res /home/wah49/database/params_files/BIC_files/BIC.params' )
import math, csv, sys, optparse, os, ast, shutil
from pyrosetta.toolbox import *
from transform import *
from SyPRIS_pkgs import *
#import rosetta.core.init
#import rosetta.core.conformation.symmetry
#import rosetta.core.pose.symmetry
#import rosetta.core.scoring.symmetry
#import rosetta.protocols.simple_moves.symmetry

def expand_res_list_to_all_chains(res_list, chain_num, max_res):
    inc = int(max_res)/int(chain_num)
    new_list = []
    for chain in xrange(chain_num):
        for res in res_list:
            new_list.append(int(res)+(chain*inc))
    return new_list

def energy_eval(pose, wt, catalytic_res, design_shell, repack_shell, name):
    catalytic_res = expand_res_list_to_all_chains(catalytic_res, pose.num_chains()-1, len(pose.sequence().rstrip('X')))
    design_shell = expand_res_list_to_all_chains(design_shell, pose.num_chains()-1, len(pose.sequence().rstrip('X')))
    repack_shell = expand_res_list_to_all_chains(repack_shell, pose.num_chains()-1, len(pose.sequence().rstrip('X')))
    
    total_energy = 0.
    shell_energy = 0.
    catalytic_res_energy = 0.
    cat_res_energy_minus_fa_elec_fa_dun_fa_internal_sol = 0.
    delta_total = 0.
    delta_shell = 0.
    delta_design_shell_fa_elec = 0.
    holes = 0.
    atom_pair_cst = pose.energies().total_energies()[rosetta.core.scoring.atom_pair_constraint]
    angle_cst = pose.energies().total_energies()[rosetta.core.scoring.angle_constraint]
    dihedral_cst = pose.energies().total_energies()[rosetta.core.scoring.dihedral_constraint]

    for res in xrange(1, len(pose.sequence().rstrip('X'))+1):
        pose_res_energy = pose.energies().residue_total_energies(res)[rosetta.core.scoring.total_score]
        print 'res: %s  energy: %s' % (res, pose_res_energy)
        wt_res_energy = wt.energies().residue_total_energies(res)[rosetta.core.scoring.total_score]
        total_energy += pose_res_energy
        delta_total += (pose_res_energy - wt_res_energy)
        if res in design_shell + repack_shell:
            shell_energy += pose_res_energy
            delta_shell += (pose_res_energy - wt_res_energy)
            if res in design_shell:
                pose_fa_elec = pose.energies().residue_total_energies(res)[rosetta.core.scoring.fa_elec]
                wt_fa_elec = wt.energies().residue_total_energies(res)[rosetta.core.scoring.fa_elec]
                delta_design_shell_fa_elec += (pose_fa_elec - wt_fa_elec)
        if res in catalytic_res:
            pose_fa_elec = pose.energies().residue_total_energies(res)[rosetta.core.scoring.fa_elec]
            pose_fa_dun = pose.energies().residue_total_energies(res)[rosetta.core.scoring.fa_dun]
            pose_fa_intra_sol_xover4 = pose.energies().residue_total_energies(res)[rosetta.core.scoring.fa_intra_sol_xover4]
            catalytic_res_energy += pose_res_energy
            cat_res_energy_minus_fa_elec_fa_dun_fa_internal_sol += (pose_res_energy - sum([pose_fa_elec, pose_fa_dun, pose_fa_intra_sol_xover4]))

    output_line = [name, catalytic_res_energy, cat_res_energy_minus_fa_elec_fa_dun_fa_internal_sol,\
                   delta_total, delta_shell, delta_design_shell_fa_elec, atom_pair_cst, angle_cst, dihedral_cst]
    output_line = [str(x) for x in output_line]
    outline = ','.join(output_line) + '\n'
    return outline

def symmetrize_pose(pose, symm):
    pose_symm_data = rosetta.core.conformation.symmetry.SymmData(pose.size(), pose.num_jump())
    pose_symm_data.read_symmetry_data_from_file(symm)
    rosetta.core.pose.symmetry.make_symmetric_pose(pose, pose_symm_data)

def enzyme_shell_selector(pose, probe_file, scaffold_file, shell_distance, exclude_res=[]):
    syp_obj = SyprisCofactor(probe_file, scaffold=scaffold_file)
    shell_pdb_res = syp_obj.find_residues_around_ligand(shell_distance)
    if shell_pdb_res is None:
        return []
    shell_pose_res = []
    for res_num in shell_pdb_res:
        pose_res_num = pose.pdb_info().pdb2pose('A', int(res_num))
        shell_pose_res.append(str(int(pose_res_num)))
    output_res_list = [x for x in shell_pose_res if x not in exclude_res]
    output_res_list.sort()
    return output_res_list


def make_cst_func(xyz):
    if len(xyz) == 2:
        dist = get_mag_difference(xyz[0], xyz[1])
        func = rosetta.core.scoring.constraints.BoundFunc(dist-0.25, dist+0.25, 0.5, "dis")
    elif len(xyz) == 3: 
        ang = get_angle(xyz[0], xyz[1], xyz[2])
        func = rosetta.core.scoring.func.HarmonicFunc(ang, math.radians(2.5)) #5.0 angle std
    elif len(xyz) == 4: 
        dihedral = get_dihedral(xyz[0], xyz[1], xyz[2], xyz[3])
        func = rosetta.core.scoring.func.CircularHarmonicFunc(dihedral, math.radians(5.0))
    return func


def parse_probe_for_csts(probe, CMP_chain_closest_to_A):
    '''parse_probe_for_csts will separate the probe into chains and then locate pairs of atoms
    across chains which have small distances. It will then determine the 6-dofs needed for enzdes
    style constraints.'''
    #final output dictionary
    master_constraint_functions = {}
    #backbone atom_types
    bb = ['N','CA','C','O']
    #chain_dict is a pairing of atom_names and xyzs
    all_xyz, all_chain_names, all_mag_vals = pdb_to_xyz(probe)
    chain_dict = dict(zip(all_chain_names, all_xyz))
    #separate the probe into symmetric chains
    probe_by_chain = symmetric_chain_separate(probe, True)
    #need to define chain1 as the chain that most closely matches chainA in the pdb
    #need the index of the CMP_chain_closest_to_A
    for i,x in enumerate(probe_by_chain):
        if x[0][21] == CMP_chain_closest_to_A:
            chain_index = i
    residue_blocks_of_main_chain = block_pdb_by_res(probe_by_chain[chain_index])
    for block_index, residue_block in enumerate(residue_blocks_of_main_chain):
        chain1_xyz, chain1_names, chain1_mags = pdb_to_xyz(residue_block)
        residue_pair_id = chain1_names[0].split('_')[0]
        #chain1_xyz, chain1_names, chain1_mags = pdb_to_xyz(probe_by_chain[chain_index])
        inner_chain_functions = {}
        
        for chain2 in probe_by_chain[1:]:
            chain2_xyz, chain2_names, chain2_mags = pdb_to_xyz(block_pdb_by_res(chain2)[block_index])
            inner_mags = []
            chain1_mag_keys = {}
            chain2_mag_keys = {}
            for i in xrange(len(chain1_xyz)):
                if chain1_names[i].split('_')[3] in bb:
                    continue
                symmetric_atom_mag = get_mag_difference(chain1_xyz[i], chain2_xyz[i])
                if symmetric_atom_mag > 0.3:
                    chain1_mag_keys[str(symmetric_atom_mag)]  = chain1_names[i]
                    chain2_mag_keys[str(symmetric_atom_mag)]  = chain2_names[i]
                    inner_mags.append(symmetric_atom_mag)
            inner_mags.sort()

            chain1_atoms = []
            chain2_atoms = []
            for mag in inner_mags[:3]:
                
                chain1_atoms.append(chain1_mag_keys[str(mag)])
                chain2_atoms.append(chain2_mag_keys[str(mag)])
            
            #from here on we can pull the distance, angle, and dihedral values using the chain#_atoms lists which are a list of 3
            #atoms that are considered to be the closest to the center of the CMP and not within the backbone set
            # inner_chain_functions = { chain1_chain# : { atom#_atom# : func . . .} . . .}
            chain_pair_name = '_'.join([chain1_names[0].split('_')[1], chain2_names[0].split('_')[1]])
            
            #these combos are the 6-dofs distance, 2x angle, 3x dihedral indices.
            combinations = [ (0,0), (1,0,0), (0,0,1), (2,1,0,0), (1,0,0,1), (0,0,1,2) ] 
            cst_functions = {}
            atom_name_combo_convert = {}
            for combo in combinations:
                #break up the combination into chain1 and chain2
                chain1_inds = combo[:combo.index(0)+1]
                chain2_inds = combo[combo.index(0)+1:]
                #find the full names for the given combo: resnum_chain_restype_atomtype
                full_names = [chain1_atoms[x] for x in chain1_inds] + [chain2_atoms[x] for x in chain2_inds]
                #pull the atomtypes from the full name and string them together for the combo
                atom_names = '_'.join([x.split('_')[3] for x in full_names])
                #find the xyz_list for each full_name
                xyz_list = [chain_dict[x] for x in full_names]
                #for each atom-name combo, add the corresponding function
                cst_functions[atom_names] = [make_cst_func(xyz_list), combo]
            
            inner_chain_functions[chain_pair_name] = cst_functions
        master_constraint_functions[residue_pair_id] = inner_chain_functions

    # inner_chain_functions = { chain1_chain# : { atom#_atom# : func . . .} . . .}
    print master_constraint_functions
    return master_constraint_functions


def apply_constraints_from_probe(pose, probe, coord_cst=False):
    if coord_cst == True:
        #write this in later
        return
    probe = read_file(probe)
    # generate the constraint functions for inner-chain ncAA
    functions = parse_probe_for_csts(probe, 'A')
    
    for residue, chain_pair_functions in functions.iteritems():
        for chain_pair, atom_set_functions in chain_pair_functions.iteritems():
            pose_res1 = pose.pdb_info().pdb2pose(chain_pair.split('_')[0], int(residue))
            pose_res2 = pose.pdb_info().pdb2pose(chain_pair.split('_')[1], int(residue))
            #for atom_set, func_combo in atom_set_functions.iteritems():
                #function, combo = func_combo
            for atom_set_string, (function, combo) in atom_set_functions.iteritems():
                print atom_set_string
                atom_set = atom_set_string.split('_')
                print atom_set
                #chain1_atom_names = [atom_set[i] for i,x in enumerate(combo[:combo.index(0)+1])]
                chain1_atom_names = atom_set[:combo.index(0)+1]
                print chain1_atom_names
                #create atom object for each atom name in set for first chain
                chain1_atoms = [rosetta.core.id.AtomID(pose.residue(pose_res1).atom_index(name), pose_res1) for name in chain1_atom_names]
                for atom in chain1_atoms:
                    print atom
                #chain2_atom_names = [atom_set[i] for i,x in enumerate(combo[combo.index(0)+1:])]
                chain2_atom_names = atom_set[combo.index(0)+1:]
                print chain2_atom_names
                #create atom object for each atom name in set for second chain
                chain2_atoms = [rosetta.core.id.AtomID(pose.residue(pose_res2).atom_index(name), pose_res2) for name in chain2_atom_names]
                for atom in chain2_atoms:
                    print atom
                #combine the atom objects for easy assignment
                all_atoms = chain1_atoms + chain2_atoms
                print all_atoms
                #determine type of constraint based on number of atom names
                if len(all_atoms) == 2:
                    constraint = rosetta.core.scoring.constraints.AtomPairConstraint(all_atoms[0], all_atoms[1], function)
                elif len(all_atoms) == 3:
                    constraint = rosetta.core.scoring.constraints.AngleConstraint(all_atoms[0], all_atoms[1], all_atoms[2], function)
                elif len(all_atoms) == 4:
                    constraint = rosetta.core.scoring.constraints.DihedralConstraint(all_atoms[0], all_atoms[1], all_atoms[2], all_atoms[3], function)
                #apply the constraint to the pose
                pose.add_constraint(constraint)
                print '\n~~~~~~~~ %s_%s_%s ~~~~~~~~\n' % (residue, chain_pair, atom_set_string)
                print pose.constraint_set()


def mutate_res(pose_input, rotamer_vals, residue, resType, ncAA):
    print 'res:', residue
    print 'restype:', convert_resname(resType)
    if convert_resname(resType) == None:
        ncAA = rosetta.core.conformation.ResidueFactory.create_residue( ncAA.name_map(resType) )
        pose_input.replace_residue(int(residue), ncAA, True)
    else:
        pose_input = mutate_residue(pose_input, int(residue), '%s' % str(resType[0])) #BPY

    for num, chi in enumerate(rotamer_vals):
        if chi == 'X':
            break
        print 'chi type:', chi
        print 'chi:', num+1
        print 'resind:', int(residue)
        if chi < 0.0:
            chi += 360.
        pose_input.set_chi(num+1, int(residue), float(chi))

    return pose_input

def fast_relax(pose, scorefxn, design_shell, repack_shell, repeats=1):
    design_shell = [int(x) for x in design_shell]
    repack_shell = [int(x) for x in repack_shell]
    #task_factory = rosetta.core.pack.task.TaskFactory()
    task_pack = standard_packer_task(pose)
    #task_pack = rosetta.core.pack.task.TaskFactory.create_packer_task(pose)
    #task_pack.restrict_to_repacking()
    #task_pack.temporarily_fix_everything()
    move_map = MoveMap()
    move_map.set_chi(False) #set to fixed rotamer
    move_map.set_bb(False) #set to fixed backbone
    #for res in repack_shell + design_shell:
    #sequence will return the single letter code + the virtual residues added during symmetry
    #we want to strip off the trailing virtual residues (denoted with single letter 'X')
    for res in xrange(1, len(pose.sequence().rstrip('X'))+1):
        if res in repack_shell:
            #task_pack.temporarily_set_pack_residue(res, True) #allow residue to repack
            #task_pack.temporarily_set_pack_residue(res, False) #allow residue to repack
            #task_pack.pack_residue(res) #allow residue to repack
            task_pack.residue_task(res).restrict_to_repacking()
            move_map.set_chi(res, True) #sample rotamers
        elif res in design_shell:
            print 'design_shell'
            print res
            #task_pack.design_residue(res) #allow designable residues to be designable
            move_map.set_bb(res, True) #allow designable residues to also sample backbone position
            move_map.set_chi(res, True) #sample rotamers
        else:
            #task_pack.temporarily_set_pack_residue(res, False) #allow residue to repack
            task_pack.residue_task(res).prevent_repacking()
    print task_pack
    #task_factory.modify_task(pose, task_pack)
    rosetta.core.pose.symmetry.make_symmetric_movemap(pose, move_map)
    pack_mover = rosetta.protocols.simple_moves.symmetry.SymPackRotamersMover(scorefxn, task_pack)
    min_mover = rosetta.protocols.simple_moves.symmetry.SymMinMover()
    min_mover.movemap(move_map)
    min_mover.score_function(scorefxn)
    #fast_relax_mover = rosetta.protocols.relax.FastRelax(scorefxn, standard_repeats=1)
    #fast_relax_mover.set_movemap(move_map)
    #fast_relax_mover.show()
    #sys.exit()
    #fast_relax_mover.set_enable_design(False)
    #fast_relax_mover.set_task_factory(task_pack)
    monte_carlo_mover = MonteCarlo(pose, scorefxn, 0.8)
    original_fa_rep_weight = scorefxn.get_weight(rosetta.core.scoring.fa_rep)
    #for i in xrange(1):#4):
    for i in xrange(repeats):
        for fa_rep_weight_multiplier in [0.02, 0.25, 0.55, 1.]:
            fa_rep_weight = original_fa_rep_weight * fa_rep_weight_multiplier
            scorefxn.set_weight(rosetta.core.scoring.fa_rep, fa_rep_weight)
            pack_mover.apply(pose)
            min_mover.apply(pose)
            #min
        
        #fast_relax_mover.apply(pose)
        monte_carlo_mover.boltzmann(pose)

    scorefxn(pose)
    scorefxn.show(pose)
    return
    


def post_SyPRIS_prep(pdb_file_name, scaffold_file, all_red_rot_vals, rotamer_restypes, rotamer_res_inds, symm_file, cof_res_dict, orig_pose, ncAA, nstruct, MC_repeats):
    #need rotamer_vals, pre_resn, rotamer_res all on same index
    scorefxn = rosetta.core.scoring.symmetry.SymmetricScoreFunction()
    scorefxn.add_weights_from_file('/home/wah49/Rosetta/main/database/scoring/weights/beta_nov15_cst.wts')
    non_symm_score = get_fa_scorefxn()
    non_symm_score.set_weight(rosetta.core.scoring.atom_pair_constraint, 1.0)
    non_symm_score.set_weight(rosetta.core.scoring.angle_constraint, 1.0)
    non_symm_score.set_weight(rosetta.core.scoring.dihedral_constraint, 1.0)
    print scaffold_file
    print pdb_file_name
    pose_from_file(orig_pose, scaffold_file) #('chainA_temp.pdb')
    wild_type = Pose()
    wild_type.assign(orig_pose)
    energy_lines_per_decoy = []
    for pose_num in xrange(1, nstruct+1):
        pose_input = Pose()
        pose_input.assign(orig_pose)
        output_name = '%s_design_%i.pdb' % (pdb_file_name[:-4], pose_num)
        #os.remove('chainA_temp.pdb')
        res_conversion = {}
        pose_focus_string = []
        for index, pre_resn in enumerate(rotamer_res_inds): 
            residue = pose_input.pdb_info().pdb2pose('A', int(pre_resn))
            print 'rosetta res', residue
            print 'chis to apply', all_red_rot_vals[index]
            print 'restype', rotamer_restypes[index]
            #all_red_rot_vals[index] is the chis with which to apply to new res
            #residue = the pose # from the res number given in rotamer_res_inds
            #rotamer_restypes[index] is the residue type we would like to mutate too
            pose_input = mutate_res(pose_input, all_red_rot_vals[index], residue, rotamer_restypes[index], ncAA)
            #res_conversion will have the chelant model res index as a key and the new pose resi as value
            res_conversion[cof_res_dict[str(pre_resn)]] = str(residue)
            pose_focus_string.append(str(residue))
        
        #find design_shell which will be 5A larger than the clash distance
        design_shell = enzyme_shell_selector(pose_input, './new_bis_his_matches/c2_matches/4cqb/' + pdb_file_name, scaffold_file, 6.5, pose_focus_string)
        #find repack_shell which will be 10A larger than the clash distance
        repack_shell = enzyme_shell_selector(pose_input, './new_bis_his_matches/c2_matches/4cqb/' + pdb_file_name, scaffold_file, 10.5, pose_focus_string + design_shell)
        print 'new_des_shell:\n', design_shell
        print 'new_rep_shell:\n', repack_shell
        #task_pack = standard_packer_task(pose_input)
        symmetrize_pose(pose_input, symm_file)
        if pose_num == 1:
            symmetrize_pose(wild_type, symm_file)
        apply_constraints_from_probe(pose_input, './new_bis_his_matches/c2_matches/4cqb/' + pdb_file_name)
        print pose_input.constraint_set()
        #sys.exit()
        #print pose_input.constraint_set()
        scorefxn(pose_input)
        if pose_num == 1:
            scorefxn(wild_type)
        scorefxn.show(pose_input)
        scorefxn.show(wild_type)
        fast_relax(pose_input, scorefxn, design_shell, repack_shell, MC_repeats)
        #evaluate pose_decoy energy and output a table
        #scorefxn(pose_input)
        energy_lines_per_decoy.append(energy_eval(pose_input, wild_type, pose_focus_string, design_shell, repack_shell, output_name))
        #last_res = pose_input.pdb_info().pdb2pose('A', int(last_pdb_res))
        pose_input.dump_pdb(output_name)
        #add_hetatms('%s_chainA.pdb' % pdb_file_name[:-4], hetatm_lines) 
        #create_coordCST(pdb_file_name, res_conversion, hetatm_lines, last_res)
        #make_flags_file(pdb_file_name, symm_file, res_conversion)
    print energy_lines_per_decoy
    with open('scores.sc', 'a') as fH:
        fH.writelines(energy_lines_per_decoy)

    '''
    Make the working directories update depending on the files
    Remove any location-based hard-code
    #I want to add nstruct to the input. Make copys of the pose before fast relax and do it to each
    #I want to change monteCarlo to 4 before production
    #I want to make a table of energy and rmsd 
       #table should include:
            #-delta total
            #-delta repack + design shell
            #-residue energy total of forced mutations + ligand (subtract the fa_dun and fa_elec and maybe internal_sol)
            #-delta fa_elec of design shell
            -holes?
            #-all contraint energies (atom_pair, angle, dihedral)
    Ideally, I would like to make this post script an option in the main SyPRIS module
    Any useful functions that could be generalized in a pyrosetta_PKGs should be converted and written there.
    #'''
    
    return 

def main(argv):

    parser = optparse.OptionParser(usage="\n\nTake SyPRIS output and mutate matched residues.")

    parser.add_option('--pdb-path', dest = 'pdb_path',
        help = 'The path to the pdbs listed in the outcsv \n \
                Example: /home/user/my_pdbs/')

    parser.add_option('-o', dest = 'final_out_path',
        help = 'The path to where the output INPUT files should be stored. \n \
                Example: /home/user/my_pdbs/')

    parser.add_option('--rot-lib', dest = 'rot_lib',
        default = '',
        help="Rotamer library dictionary file: ( default = '' )")

    parser.add_option('--params-file', dest = 'params_file',
        default = '',
        help="Specify path to params file when mutating to ncAA.")

    parser.add_option('--nstruct', type=int, dest = 'nstruct',
        default = 1,
        help="Specify number of unqiue decoys to design.")

    parser.add_option('--monte-carlo-repeats', type=int, dest = 'monte_carlo_repeats',
        default = 1,
        help="Specify number of times a single is passed through iterative design and subject to MonteCarlo.")

    
    (options,args) = parser.parse_args()
   
    ncAA_opt = Pose()
    try:
        ncAA = generate_nonstandard_residue_set( ncAA_opt, ['%s' % options.params_file]  )
    except RuntimeError:
        ncAA = []
    #read off the pdb, and send over
    if options.rot_lib:
        with open(options.rot_lib, 'r') as f:
            s = f.read()
            my_dict = ast.literal_eval(s)
    

    pdb_file_line = options.pdb_path
    pdb_file_line = pdb_file_line.strip('\r\n')
    split_file_line = pdb_file_line.split(',')

    full_scaffold = split_file_line[3]
    pdb_file_name = full_scaffold[:-4] + '_' + '_'.join(split_file_line[1:3]) + '_score_' + split_file_line[-1] + '.pdb'
    chis_to_change = split_file_line[4:-1]

    #pdb_obj = Transform(read_file('./new_bis_his_matches/c2_matches/' + pdb_file_name))
    pdb_file_name_split = pdb_file_name.split('_')
    #7 for C_2_up/down
    #4 if the C_2_up wasn't there
    pdb_split_res_types = pdb_file_name_split[ pdb_file_name_split.index('relaxed')+7:\
                                               pdb_file_name_split.index('score') ][::2]
    pdb_split_match_types = pdb_file_name_split[ pdb_file_name_split.index('relaxed')+7:\
                                                 pdb_file_name_split.index('score') ][1::2]
    print '\n'
    print pdb_split_res_types
    print pdb_split_match_types
    print '\n'
    rotamers = []
    rotamer_restypes = []
    rotamer_res_inds = []
    cof_res_dict = {}
    for ind, res_type in enumerate(pdb_split_res_types):
        residue_in_scaff = ''
        chelant_resi = ''
        rot_ind = 0
        for i, chara in enumerate(pdb_split_match_types[ind]):
            try:
                a = int(chara)
                chelant_resi += chara
            except ValueError:
                residue_in_scaff += chara
                if len(residue_in_scaff) == 3:
                    rot_ind = (i+1)
                    break
        print residue_in_scaff
        print pdb_split_match_types[ind]
        pre_rotamer = pdb_split_match_types[ind][rot_ind:]
    
        print pre_rotamer
        rotamer = ''
        bad_indices = []
        for ind, char in enumerate(pre_rotamer):
            print char
            if char in ['a','b']:
                bad_indices.append(ind +1)
            else:
                if ind not in bad_indices:
                    rotamer += char

        rotamers.append(rotamer)
        rotamer_restypes.append(residue_in_scaff)
        rotamer_res_inds.append(res_type[:-3])
        cof_res_dict[res_type[:-3]] = chelant_resi 
    
    if not rotamers:
        sys.exit()
    print 'rotamers:', rotamers
    print 'restypes:', rotamer_restypes
    print 'res indices:', rotamer_res_inds
    print 'cof_res_dict:\n', cof_res_dict
    pdb_obj = Transform(read_file('./new_bis_his_matches/c2_matches/4cqb/' + pdb_file_name))
    cof_residue_blocks = block_pdb_by_res(read_file('./new_bis_his_matches/c2_matches/4cqb/' + pdb_file_name))
    cof_residue_blocks = [ x for x in cof_residue_blocks if x[0][17:20].strip().upper() not in rotamer_restypes ]
    hetatm_lines = unblock_pdb_by_res(cof_residue_blocks)
    hetatm_lines = [ x for x in hetatm_lines if x[21] == 'A' ] 
    print 'hetatm_lines', hetatm_lines

    all_red_rot_vals = []
    for indicy, rotamer in enumerate(rotamers):
        chelant_resi_key = rotamer_res_inds[indicy]
        rotamer_vals = get_first_three_chis_of_resnum(copy_transform_object(pdb_obj), [int(cof_res_dict[chelant_resi_key])])
        all_red_rot_vals.append(rotamer_vals[:])
        
    all_red_rot_vals[0] = chis_to_change[:]
    print 'all chis:', all_red_rot_vals
    #print chis_to_change
    #sys.exit()
    scaff_pruned_name = '_'.join(pdb_file_name_split[:(pdb_file_name_split.index('centered')+4)])
    database_path = '/home/wah49/database/pdbs/SCAFF/c2_pdbs/postRelax/'
    scaffold_file = database_path + scaff_pruned_name + '_INPUT.pdb'
    symm_file_name = database_path + scaff_pruned_name + '.symm'
    #symm_file_name = '_'.join(pdb_file_name_split[:(pdb_file_name_split.index('standard')-1)]) + '.symm'

    scaff_obj = Transform(read_file(scaffold_file))
    last_pdb_res = scaff_obj.get_xyz_info(-1)[-1].split('_')[0]
    
    post_SyPRIS_prep(pdb_file_name, scaffold_file, all_red_rot_vals, rotamer_restypes, rotamer_res_inds,\
                     symm_file_name, cof_res_dict, ncAA_opt, ncAA, options.nstruct, options.monte_carlo_repeats)

if __name__ == '__main__':
    #pyrosetta.init( extra_options='-ignore_zero_occupancy false')# -extra_res /home/wah49/database/params_files/BIC_files/BIC.params' )
    #res_set2 = generate_nonstandard_residue_set( Vector1( ['/home/wah49/params_files/CoO4.params'] ) )
    main(sys.argv[1:])
