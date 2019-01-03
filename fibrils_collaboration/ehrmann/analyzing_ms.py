#!/usr/bin/python
from numpy import count_nonzero
tau_seq = 'MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL'

"""
with open('fibril_wt_set-d-selected-columns-relevant-2N4R-records.txt', 'r') as r:
ms_raw = [line.strip().split('\t') for line in r.readlines()]

header=ms_raw.pop(0)

# Cuts treated as downstream, i.e. a cut at residue 2 breaks the bond between residues 2 and 3, not 1 and 2
cut_counts = [[0 for j in range(440)] for i in range(4)]
cut_intensities = [[0 for j in range(440)] for i in range(4)]
for line in ms_raw:
	if not 'N/A' in line[0:2]:
		cleave_sites = [int(line[0])-1, int(line[1])]
		# cleavage is upstream of peptide start, downstream of peptide end
		expt_count = [float(count_nonzero([int(i) for i in line[ 7:10]])) / 3, float(count_nonzero([int(i) for i in line[10:13]])) / 3, float(count_nonzero([int(i) for i in line[13:16]])) / 3, float(count_nonzero([int(i) for i in line[16:19]])) / 3]
		expt_intensity = [sum([int(i) for i in line[ 7:10]]), sum([int(i) for i in line[10:13]]), sum([int(i) for i in line[13:16]]), sum([int(i) for i in line[16:19]])]
		for i in range(4):
			for site in cleave_sites:
				if site != 441:
					cut_counts[i][site - 1] += expt_count[i]
					cut_intensities[i][site - 1] += expt_intensity[i]


with open('intensities_by_site.txt', 'w') as o:
o.write('site\t' + '\t'.join(['10min', '60min', '180min', 'overnight'] * 2) + '\n')
for n, i in enumerate(range(1,440)):
ints = [str(i)]
ints += [str(cut_counts[x][n]) for x in range(4)]
ints += [str(cut_intensities[x][n]) for x in range(4)]
o.write('\t'.join(ints) + '\n')
"""


with open('Export_MS_data.txt', 'r') as r:
	txt=r.readlines()
with open('ehrmann_tau_ms_excerpts.txt','w') as w:
	w.write(txt[0])
	for l in txt[2::2]:
		if l.split('\t')[0] in tau_seq:
			w.write(txt[1])
			w.write(l)

def pick_expt_columns(first_expt_col, num_expts, num_cond, reps):
	"""
	Generates lists of the columns within the mass spec data that correspond 
	to each experiment and condition

	Requires list to be ordered by UMSAP

	first_expt_col is the first column with experimental data
	num_expts is the number of different experimental sets
		ex: WT, WT+SA, WT+SA_dPDZ
	num_cond is the number of conditions for each experimental set
		ex: 4 if samples were taken at 10min, 1 hour, 3 hours, and overnight
	reps is the number of repetitions for each point, ex: 3
	"""
	num_cols_per_expt = num_cond * reps
	last_expt_col = first_expt_col + num_expts * num_cols_per_expt
	e_cols = []

	for i in range(first_expt_col, last_expt_col, num_cols_per_expt):
		c = [range(j, j + reps) for j in range(i, i + num_cols_per_expt, reps)]
		e_cols.append(c)

	return e_cols


class expt_ms_set:
	def __init__(self, ms_set, expt_cols, prot_len):
		self.ms_set = ms_set
		self.expt_cols = expt_cols
		self.prot_len = prot_len
		self.cut_counts = [[0 for j in range(prot_len - 1)] for i in range(len(self.expt_cols))]
		self.cut_intensities = [[0 for j in range(prot_len - 1)] for i in range(len(self.expt_cols))]

		self.populate_counts()

	def populate_counts(self):
		"""
		Go through each line (with each line corresponding to an identified 
		protein fragment) of the data file and do the folowing:
		1: Identify the cleavage sites demonstrated by that fragment. These 
		are the residue just before the fragment and the last residue in the 
		fragment, unless either is a terminal residue. These numbers are in 
		the first column of each line.
		2: In the columns specific to this data set (given at construction),
		collect two values across all repetitions for each sample condition.
		The values are the average number of reps in which any cleavage 
		activity was identified, and the total intensity of the identified 
		fragment.
		3: Add the collected values to the appropriate sites in a list of all
		sites in the target protein. Thus multiple fragments with an end 
		corresponding to a cut site will produce values for that site which 
		are the sum of all counts/sums for each fragment.
		"""
		for line in self.ms_set:
			# cleavage is upstream of peptide start, downstream of peptide end
			cleave_sites = [int(line[0])-1, int(line[1])]

			frag_counts = []
			frag_intensities = []
			for c in self.expt_cols:
				frag_counts.append(float(count_nonzero([int(line[i]) for i in c])) / len(c))
				frag_intensities.append(sum([int(line[i]) for i in c]))

			for i in range(len(self.expt_cols)):
				for site in cleave_sites:
					if site not in [1, self.prot_len]:
						self.cut_counts[i][site - 1] += frag_counts[i]
						self.cut_intensities[i][site - 1] += frag_intensities[i]

my_cols = pick_expt_columns(7, 6, 4, 3)
the_sets = []
for x in my_cols:
	the_sets.append(expt_ms_set(ms_raw, x, 441))

lines_out = []
for i in range(440):
	line = []
	for j in range(6):
		for k in range(4):
			line.append(the_sets[j].cut_counts[k][i])
		for k in range(4):
			line.append(the_sets[j].cut_intensities[k][i])
	lines_out.append(line)

with open('intensities_by_site.txt', 'w') as o:
	for l in lines_out:
		line_text = '\t'.join([str(i) for i in l]) + '\n'
		o.write(line_text)
