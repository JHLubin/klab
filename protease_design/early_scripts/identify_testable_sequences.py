#!/usr/bin/python
p6_list = {
	'D': ['R138K', 'V173I', 'T175K', 'T175R', 'T175Q', 'D183Y'],
	'E': ['R138K', 'V173I', 'D183R', 'D183F', 'D183L'],
	'F': ['T175L', 'T175V', 'T175I', 'T175K'],
	'H': ['V173L', 'C174A', 'T175L', 'T175V'],
	'L': ['V173M', 'V173L', 'C174A', 'T175L', 'T175I', 'T175K']
	'M': ['R138M', 'C174A', 'T175I', 'T175V'],
	'P': ['T175I'],
	'Q': ['T175E', 'T175I'],
	'R': ['R138E', 'R138V', 'R138I', 'V173I', 'T175E'],
	'T': ['T175I'],
	'V': ['V173I','C174A'],
	'W': ['R138I', 'V173I', 'T175I', 'T175V']}
p5_list = {
	'A': ['C174Y', 'C174H'], # Depending on p3 charge
	'C': ['C174Y', 'C174H'], # Depending on p3 charge
	'E': ['C174A'],}
p4_list = {}
p3_list = {}
p2_list = {}

with open('pep_seq_lists/all_HCV_uncleaved.txt', 'r') as r:
	list_uncleaved = [i.strip() for i in r.readlines()]

with open('pep_seq_lists/all_HCV_cleaved.txt', 'r') as s:
	list_cleaved = [i.strip() for i in s.readlines()]

uncleaved = []
untested = []

for p6 in p6_list:
	for p5 in p5_list:
		for p4 in p4_list:
			for p3 in p3_list:
				for p2 in p2_list:
					seq = p6 + p5 + p4 + p3 + p2 + 'C'
					if seq in list_uncleaved:
						uncleaved.append(seq)
					elif seq not in list_cleaved:
						untested.append(seq)

print "uncleaved:"
for seq in uncleaved:
	print '\t', seq

print "untested:"
for seq in untested:
	print '\t', seq

