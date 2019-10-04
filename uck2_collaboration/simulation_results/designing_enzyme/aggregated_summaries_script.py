from glob import glob

ala = glob('all_ala/*.csv') 
wts = glob('wt_start/*.csv')

in_headers = []

for i in ala + wts:
	with open(i, 'r') as r:
		in_headers.append(r.readline().strip())

master_header = []

for n, h in enumerate(in_headers):
	this_head = h.split(',')
	if n == 0:
		master_header = this_head
		continue
	for x, e in enumerate(this_head):
		if e not in master_header:
			prev_elem = master_header.index(this_head[x-1])
			master_header.insert(prev_elem + 1, e)

master_header.insert(0,'group')
master_out = []
master_out.append(master_header)

for f in ala + wts:
	with open(f, 'r') as r:
		summary_lines = r.readlines()
		for y, l in enumerate(summary_lines):
			input_line = l.strip().split(',')
			if y == 0:
				head_gaps = []
				for x, e in enumerate(master_header):
					if e not in input_line:
						head_gaps.append(x)
				continue
			out_line = input_line[:]
			for hg in head_gaps:
				out_line.insert(hg, '')
			if f in ala:
				out_line.insert(0, 'all_ala')
			if f in wts:
				out_line.insert(0, 'wt_start')
			master_out.append(out_line)

with open('full_summary.csv', 'w') as w:
	for line in master_out:
		w.write(','.join([str(i) for i in line]))

mut_list = []

for n, i in enumerate(master_header):
	if n < 11:
		continue
	for line in master_out[1:]:
		if line[n + 1] != '':
			mut = i + line[n + 1]
			if mut not in mut_list:
				mut_list.append(mut)

master_list = []
for line in master_out[1:]:
	rep_line = line[:12]
	for n in range(12,len(line)):
		if line[n] in ['', '\n']:
			rep_line.append('')
		else:
			mut = master_header[n - 1]+ line[n]
			rep_line.append(mut)
	master_list.append(rep_line)

mut_freq = []

for i in mut_list:
	count = 0
	for line in master_list:
		if i in line:
			count += 1
	mut_freq.append(count)

mut_rates = [round(float(i)/len(master_out[1:]),3) for i in mut_freq]

for i in zip(mut_list, mut_rates):
	print(str(i).replace('(', '').replace(')', '').replace("'","").replace(',','\t'))

lig_list = {'5amc', '5cpc', '5dzc', '5ec', 'c35amc', 'cycloc', 'm5c', 'wtc', '5amu', '5cpu', '5dzu', '5eu', 'c35amu', 'cyclou', 'm5u', 'wtu'}
mut_for_ligand = []

for mut in mut_list:
	yes_ligs = []
	for lig in lig_list:
		for line in master_list:
			if mut in line:
				if any(lig in i for i in line):
					yes_ligs.append(lig)
					break
	mut_for_ligand.append(yes_ligs)

for i in zip(mut_list, mut_rates, mut_for_ligand):
	print(i)