cd ml_training/
with open('uniq_cleaved_uncleaved_structseq_binary.csv', 'rb') as r:
	seqs = r.readlines()

with open('WT_nextseq_middle_structseq_binary.csv', 'rb') as r:
	seqs += r.readlines()

aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

[i.decode("utf-8")  for i in seqs[0].split(b',')]

ident_seqs = ['DDEQR', 'DDSQD', 'DDSQR', 'DITVP', 'DNPVP', 'DQQDD', 'DQSED', 'DQSER', 'EDHQD', 'EDQQD', 'EDSQD', 'EIDVP', 'EITVP', 'EKQED', 'EKSED', 'ENPVP', 'ENSQD', 'ENSVP', 'ENTVP', 'EQDVP', 'EQHED', 'EQNEK', 'EQQED', 'EQSED', 'EQSEN', 'ERHED', 'ERQED', 'ERSED', 'EWHED', 'KIEVP', 'KKEED', 'KNEQR', 'KQEED', 'KREED', 'NNSED', 'QDEQH', 'QHSSD', 'QQEEE', 'QQEEH', 'QQEEK', 'QQEER', 'QQESH', 'QQSED', 'QRSED', 'RDEQH', 'RDEQK', 'RHEED', 'RHESD', 'RKDEE', 'RKEED', 'RQEED', 'RQEEQ', 'RQEER', 'RQESD', 'RRDEH', 'RREED', 'RREEK', 'RREEQ', 'RSEQD', 'RWEED', 'SQQED', 'TIDIP', 'TIDVP', 'TIEVP', 'TQEVP']

other_seqs = ['DEMEE', 'DVDAR', 'FEDFQ', 'LEEFF', 'LEEYQ']

for l in seqs:
	sqcol = [int(i)  for i in l.split(b',')[:100]]
	s = ''
	for i in range(20,101,20):
		frame = sqcol[i-20:i]
		s+=aas[frame.index(1)]
	if s in ident_seqs:
		print(s, l)
# No hits

l2r = []
for isq in ident_seqs+other_seqs:
	ssum = 0
	for n, i in enumerate(isq[::-1]):
		ssum += aas.index(i)*(20**n)
	l2r.append(ssum)

l2r.sort()
relevants = []
with open('structure_features.csv', 'rb') as r:
	for i, line in enumerate(r):
		if i-1 in l2r:
			relevants.append(line)

# header: sequence,prot,pept,cst,amber
sqchecks = [i.decode("utf-8").split(',') for i in relevants]

for s in sqchecks:
	seq = s[0]
	disc = 1 * float(s[1]) + 1 * float(s[2]) + 3.5 * float(s[3]) + 0.5 * float(s[4])
	print(seq, disc)

stat_max_min = {'CLEAVED':[-80,'',-80,''], 'MIDDLE':[-80,'',-80,''], 'UNCLEAVED':[-80,'',-80,'']}
for l in seqs:
	sqcol = [int(i)  for i in l.split(b',')[:100]]+[float(i) for i in l.split(b',')[100:104]] + [l.split(b',')[104].decode("utf-8").replace('\n','')]
	s = ''
	for i in range(20,101,20):
		frame = sqcol[i-20:i]
		s+=aas[frame.index(1)]
	discrim = 1 * float(sqcol[100]) + 1 * float(sqcol[101]) + 3.5 * float(sqcol[102]) + 0.5 * float(sqcol[103])
	stat = sqcol[104]
	if discrim > stat_max_min[stat][0]:
		stat_max_min[stat][0] = discrim
		stat_max_min[stat][1] = s
	if discrim < stat_max_min[stat][2]:
		stat_max_min[stat][2] = discrim
		stat_max_min[stat][3] = s

""" 
{'CLEAVED': [-62.075405, 'RAAKR', -101.57879450000001, 'DEMEE'],
 'MIDDLE': [-62.0405435, 'KRTRG', -91.84807095000001, 'DPGMP'],
 'UNCLEAVED': [-59.59537900000001, 'KKNKR', -95.89509395, 'EEPCH']}
"""

"""
DDEQR -91.008927
DDSQD -91.383711
DDSQR -86.80095700000001
DEMEE -101.57879450000001
DITVP -86.24992599999999
DNPVP -88.1186605
DQQDD -90.09225049999999
DQSED -95.28553249999999
DQSER -88.252106
DVDAR -82.36064999999999
EDHQD -89.4285845
EDQQD -92.1335105
EDSQD -93.20112
EIDVP -90.91604099999999
EITVP -86.39840500000001
EKQED -90.93129149999999
EKSED -90.4899645
ENPVP -89.39027
ENSQD -88.921534
ENSVP -88.165986
ENTVP -87.90341399999998
EQDVP -90.6553105
EQHED -95.0163185
EQNEK -89.60940550000001
EQQED -95.30623449999999
EQSED -95.78073499999999
EQSEN -92.36306400000001
ERHED -87.11646400000001
ERQED -89.59599450000002
ERSED -91.43300550000001
EWHED -91.7741385
FEDFQ -83.3509535
KIEVP -78.048999
KKEED -82.582509
KNEQR -75.987473
KQEED -87.888544
KREED -84.94144049999998
LEEFF -86.2917415
LEEYQ -87.871025
NNSED -87.71768349999999
QDEQH -87.23800250000001
QHSSD -81.65148049999999
QQEEE -93.9059205
QQEEH -90.8511925
QQEEK -87.16291949999999
QQEER -88.52771500000001
QQESH -82.2175145
QQSED -88.626356
QRSED -86.18238349999999
RDEQH -80.83677750000001
RDEQK -79.5090615
RHEED -87.4543375
RHESD -79.58954650000001
RKDEE -79.33683049999999
RKEED -82.06831749999999
RQEED -88.656296
RQEEQ -85.10617950000001
RQEER -82.5373455
RQESD -81.111642
RRDEH -82.91263314999999
RREED -83.66036700000001
RREEK -78.36381600000001
RREEQ -83.98775065000002
RSEQD -79.160188
RWEED -86.734922
SQQED -87.75695200000001
TIDIP -82.065458
TIDVP -83.5380085
TIEVP -84.44197049999998
TQEVP -83.958235
"""