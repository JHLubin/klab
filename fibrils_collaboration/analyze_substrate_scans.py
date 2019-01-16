import argparse
from glob import glob
from operator import itemgetter
from os.path import basename, join
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import ScoreType as st
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector, NeighborhoodResidueSelector
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from statistics import stdev

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("directory", type=str, 
		help="Read fasc files from what directory?")
	parser.add_argument("-z", "--zipped", action="store_true", default=False,
		help="Add .gz suffix to decoys?")
	parser.add_argument("-kf", "--keep_folder", action="store_true", 
		default=False, help="Keep the folder name of the decoy? (By default, \
		only the basename will appear.)")
	parser.add_argument("-c", "--constraints", type=str, default='htra1_protease.cst',
		help="Pick a constraints file for the enzyme")
	args = parser.parse_args()
	return args


def collect_fasc_files(folder):
	"""
	For a given folder, collects a list of all .fasc files
	Excludes any combined_reports files produced by condense_fasc.py
	"""
	# Search for .fasc files
	base_name = basename(folder.rstrip('/'))
	folder_search = join(folder, "*.fasc")
	fasc_files = glob(folder_search)
	fasc_files.sort()

	# Remove aggregate reports from condense_fasc.py
	for n, i in enumerate(fasc_files):
		if 'combined_reports' in i:
			fasc_files.pop(n)

	return sorted(fasc_files)


def text_to_numbers(in_list):
	""" 
	Takes a list that includes numbers and text all stored as strings.
	Returns the same list with the numbers converted to floats.
	"""
	new_list = []
	for i in in_list:
		try:
			new_list.append(float(i))
		except ValueError:
			new_list.append(i)

	return new_list


def fasc_to_dict_list(fasc):
	"""
	Reading in a .fasc file, converts each line to a dict, with score terms 
	as keys and number values. A list of these dicts is collected, and sorted
	based on total score.
	"""
	# Read in fasc file
	with open(fasc, 'r') as read:
		f_lines = read.readlines()
		f_lines.pop(0) # First line is not useful

	# Collecting all lines in the fasc file
	lines_data = []
	for l in f_lines:
		# Converting line from a string to lists
		text_2_list = l.split()
		line_headings = [i[:-1] for i in text_2_list[0::2]] # remove colons
		line_data = text_to_numbers(text_2_list[1::2]) # convert numbers

		# Combining lists into a dict
		line_dict = {}
		for i in range(len(line_headings)):
			line_dict[line_headings[i]] = line_data[i]

		# Adding dict to list
		lines_data.append(line_dict)

	# Sort list and return
	return sorted(lines_data, key=itemgetter('total_score'))


def add_unconstrained_scores(score_dict):
	"""
	If a scores dict includes constraints, calculate the total constraint score 
	and the total unconstrained score, and add these to the dict
	"""
	cst_terms = ['atom_pair_constraint', 'coordinate_constraint', 
				'angle_constraint', 'dihedral_constraint']
	# Calculate constraint score, if constraints are present
	cst_score = 0
	for ct in cst_terms:
		try:
			cst_score += score_dict[ct]
		except:
			pass

	# Calculate unconstrained total score
	uc_tot = score_dict['total_score'] - cst_score

	# Add scores into dict
	score_dict['cst_total'] = cst_score
	score_dict['unconstrained_total_score'] = uc_tot
	return


def remove_outliers(dict_list):
	"""
	Checks through a list of scoreline dicts and knocks out any with outlying
	constraint scores. This is done by finding the minimum constraint score 
	and the standard deviation in constraint scores, and removing any lines 
	that have constraint scores more than two standard deviations higher than 
	the minimum from the list.
	"""
	constraint_scores = [i['cst_total'] for i in dict_list]
	min_cst = min(constraint_scores)
	cst_sd = stdev(constraint_scores)
	threshold = min_cst + (2 * cst_sd)

	for n, i in enumerate(dict_list):
		if i['cst_total'] > threshold:
			dict_list.pop(n)

	return dict_list


def filter_by_constraints(dict_list):
	"""
	Adds total constraint and total unconstrained scores to score lines, 
	removes any with outlying high constraint scores, and sorts by 
	unconstrained total scores.
	"""
	# Add in extra scoring terms
	for i in dict_list:
		add_unconstrained_scores(i)

	# Eliminating outliers
	cleaned_list = remove_outliers(dict_list)

	# Sort list and return
	return sorted(cleaned_list, key=itemgetter('unconstrained_total_score'))


def selector_to_list(pose, selector):
	""" Converts a selector output vector to a list of selected residues """
	selection_vector = selector.apply(pose)
	selection_list = []
	for i in range(len(selection_vector)): 
		if selection_vector[i+1]==1:
			selection_list.append(i+1)

	return selection_list 


def fix_file(pdb):
	bn = basename(pdb)
	last_res = bn[19]

	aa = {'A':'ALA', 'C':'CYS','D':'ASP','E':'GLU', 'F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}

	with open(pdb, 'r') as r:
		lines = r.readlines()

	lines.pop(6)
	lines.pop(6)

	for l in lines[2:6]:
		lines[lines.index(l)]=l.replace('VAL',aa[last_res])

	with open(pdb,'w') as w:
		w.writelines(lines)


def get_pose_discrimination_scores(pdb):
	fix_file(pdb)
	num = basename(pdb)[11:14]
	pose = pose_from_pdb(pdb)

	# Enzdes constraints
	cstm = AddOrRemoveMatchCsts()
	cstm.set_cst_action(ADD_NEW)
	cstm.apply(pose)

	sf = create_score_function('ref2015_cst')
	sf(pose)

	pep_res = selector_to_list(pose, ChainSelector('B'))

	neigh = NeighborhoodResidueSelector()
	neigh.set_distance(8)
	neigh.set_focus_selector(ChainSelector('B'))
	neigh.set_include_focus_in_subset(False)

	prot_res = selector_to_list(pose, neigh)

	pep_score = 0
	for i in pep_res:
		pep_score += pose.energies().residue_total_energies(i)[st.total_score]

	prot_score = 0
	for i in prot_res:
		prot_score += pose.energies().residue_total_energies(i)[st.total_score]

	cst_score = 0
	for i in [st.atom_pair_constraint, st.coordinate_constraint, 
			st.angle_constraint, st.dihedral_constraint]:
		cst_score += pose.energies().total_energies()[i]

	discriminator = 1 * prot_score + 1 * pep_score + 3.5 * cst_score

	print(prot_score, pep_score, cst_score, discriminator)

	return [num, prot_score, pep_score, cst_score, discriminator]


def matching_seqs(seq_1, seq_2, frame_size, threshold):
	for i in range(len(seq_1))[:(frame_size-1)]:
		aframe = seq_1[i:i+frame_size]
		for j in range(len(tau_seq))[:(frame_size-1)]:
		    bframe = seq_2[j:j+frame_size]
		    match = 0
		    for k in range(frame_size):
		        if aframe[k] == bframe[k]:
		            match += 1
		    if match>=threshold:
		        print(i, aframe, j, bframe)


def main(args):
	# Get list of fasc files
	fasc_files = collect_fasc_files(args.directory)

	best_decoys = []
	for f in fasc_files:
		# Convert file to dict list, clean up outliers
		f_converted = fasc_to_dict_list(f)
		f_cleaned = filter_by_constraints(f_converted)

		# Collect single best model
		best_decoys.append(f_cleaned[0])

	#for n, i in enumerate(best_decoys):
	#	with open('dummy.txt', 'a') as w:
	#		if n == 0:
	#			w.write('\t'.join(list(i.keys()))+'\n')
	#		w.write('\t'.join([str(i) for i in list(i.values())])+'\n')

	# Get residue scores
	pdbs_to_check = ['scan_simulations/new_cat_scan_asyn/htra1_scan_005_MDVFM_6.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_006_DVFMK_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_007_VFMKG_7.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_008_FMKGL_7.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_009_MKGLS_1.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_010_KGLSK_1.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_011_GLSKA_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_012_LSKAK_1.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_013_SKAKE_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_014_KAKEG_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_015_AKEGV_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_016_KEGVV_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_017_EGVVA_1.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_018_GVVAA_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_019_VVAAA_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_020_VAAAE_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_021_AAAEK_1.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_022_AAEKT_3.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_023_AEKTK_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_024_EKTKQ_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_025_KTKQG_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_026_TKQGV_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_027_KQGVA_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_028_QGVAE_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_029_GVAEA_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_030_VAEAA_9.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_031_AEAAG_9.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_032_EAAGK_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_033_AAGKT_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_034_AGKTK_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_035_GKTKE_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_036_KTKEG_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_037_TKEGV_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_038_KEGVL_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_039_EGVLY_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_040_GVLYV_6.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_041_VLYVG_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_042_LYVGS_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_043_YVGSK_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_044_VGSKT_1.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_045_GSKTK_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_046_SKTKE_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_047_KTKEG_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_048_TKEGV_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_049_KEGVV_1.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_050_EGVVH_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_051_GVVHG_6.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_052_VVHGV_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_053_VHGVA_3.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_054_HGVAT_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_055_GVATV_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_056_VATVA_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_057_ATVAE_7.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_058_TVAEK_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_059_VAEKT_7.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_060_AEKTK_3.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_061_EKTKE_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_062_KTKEQ_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_063_TKEQV_6.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_064_KEQVT_3.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_065_EQVTN_9.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_066_QVTNV_6.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_067_VTNVG_7.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_068_TNVGG_3.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_069_NVGGA_9.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_070_VGGAV_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_071_GGAVV_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_072_GAVVT_1.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_073_AVVTG_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_074_VVTGV_6.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_075_VTGVT_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_076_TGVTA_6.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_077_GVTAV_6.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_078_VTAVA_7.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_079_TAVAQ_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_080_AVAQK_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_081_VAQKT_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_082_AQKTV_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_083_QKTVE_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_084_KTVEG_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_085_TVEGA_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_086_VEGAG_1.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_087_EGAGS_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_088_GAGSI_3.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_089_AGSIA_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_090_GSIAA_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_091_SIAAA_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_092_IAAAT_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_093_AAATG_9.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_094_AATGF_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_095_ATGFV_3.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_096_TGFVK_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_097_GFVKK_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_098_FVKKD_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_099_VKKDQ_9.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_100_KKDQL_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_101_KDQLG_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_102_DQLGK_3.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_103_QLGKN_9.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_104_LGKNE_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_105_GKNEE_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_106_KNEEG_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_107_NEEGA_7.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_108_EEGAP_6.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_109_EGAPQ_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_110_GAPQE_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_111_APQEG_1.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_112_PQEGI_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_113_QEGIL_3.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_114_EGILE_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_115_GILED_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_116_ILEDM_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_117_LEDMP_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_118_EDMPV_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_119_DMPVD_1.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_120_MPVDP_8.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_121_PVDPD_7.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_122_VDPDN_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_123_DPDNE_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_124_PDNEA_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_125_DNEAY_7.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_126_NEAYE_3.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_127_EAYEM_9.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_128_AYEMP_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_129_YEMPS_0.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_130_EMPSE_5.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_131_MPSEE_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_132_PSEEG_4.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_133_SEEGY_6.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_134_EEGYQ_1.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_135_EGYQD_2.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_136_GYQDY_9.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_137_YQDYE_6.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_138_QDYEP_9.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_139_DYEPE_9.pdb', 'scan_simulations/new_cat_scan_asyn/htra1_scan_140_YEPEA_4.pdb']
	for i in pdbs_to_check:
		o = get_pose_discrimination_scores(i)

		with open('dummy.txt', 'a') as w:
			w.write('\t'.join([str(i) for i in o])+'\n')

if __name__ == '__main__':
	args = parse_args()

	opts = '-cst_fa_weight 1.0 -run:preserve_header -enzdes::cstfile {}'
	init(opts.format(args.constraints))

	main(args)