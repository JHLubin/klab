from pyrosetta import *
from glob import glob
from pyrosetta.rosetta.core.simple_metrics.metrics import SequenceMetric
from pyrosetta.rosetta.core.select.residue_selector import \
    AndResidueSelector, ChainSelector, InterGroupInterfaceByVectorSelector,\
        NotResidueSelector, OrResidueSelector, ResidueIndexSelector
init('-mute all')
pdbs = glob('*designed*.pdb')
pdbs.sort()

orig_seq = a.calculate(pose)
for pdb in pdbs:
    p = pose_from_pdb(pdb)
    p_seq = a.calculate(p)
    muts = []
    for n, i in enumerate(orig_seq):
        if p_seq[n] != i:
            muts.append(i+p.pdb_info().pose2pdb(n+1).split()[0]+p_seq[n])
    print(pdb, ' '.join(muts))

pose = pose_from_pdb('../start_proteases/HCV.pdb')
orig_seq = a.calculate(pose)
chb = ChainSelector('B')
b = SequenceMetric(chb)
pep_seq = b.calculate(pose)
pdbs = glob('*designed*.pdb')
pdbs.sort()
for pdb in pdbs:
    p = pose_from_pdb(pdb)
    print(pdb, b.calculate(p)[1:6])

with open('uniq_cleaved_uncleaved_structseq_binary.csv', 'rb') as r:
    seqs = r.readlines()
with open('WT_nextseq_middle_structseq_binary.csv', 'rb') as r:
    seqs += r.readlines()
aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'V', 'W']
substring = [int(i)  for i in seqs[0].split(b',')[:100]]

ident_seqs = ['DDEQR', 'DDSQD', 'DDSQR', 'DITVP', 'DNPVP', 'DQQDD', 'DQSED', 'DQSER', 'EDHQD', 'EDQQD', 'EDSQD', 'EIDVP', 'EITVP', 'EKQED', 'EKSED', 'ENPVP', 'ENSQD', 'ENSVP', 'ENTVP', 'EQDVP', 'EQHED', 'EQNEK', 'EQQED', 'EQSED', 'EQSEN', 'ERHED', 'ERQED', 'ERSED', 'EWHED', 'KIEVP', 'KKEED', 'KNEQR', 'KQEED', 'KREED', 'NNSED', 'QDEQH', 'QHSSD', 'QQEEE', 'QQEEH', 'QQEEK', 'QQEER', 'QQESH', 'QQSED', 'QRSED', 'RDEQH', 'RDEQK', 'RHEED', 'RHESD', 'RKDEE', 'RKEED', 'RQEED', 'RQEEQ', 'RQEER', 'RQESD', 'RRDEH', 'RREED', 'RREEK', 'RREEQ', 'RSEQD', 'RWEED', 'SQQED', 'TIDIP', 'TIDVP', 'TIEVP', 'TQEVP']
other_seqs = ['DEMEE', 'DVDAR', 'FEDFQ', 'LEEFF', 'LEEYQ']
for l in seqs:
    sqcol = [int(i)  for i in l.split(b',')[:100]]
    s = ''
    for i in range(20,101,20):
        frame = sqcol[i-20:i]
        s+=aas[frame.index(1)]
    if s in ident_seqs:
        print(l)

aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
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

sqchecks = [i.decode("utf-8").split(',') for i in relevants]
[x[0] + float(x[j]) for x in sqchecks[:5] j in [1:3]]
for s in sqchecks:
    seq = s[0]
    disc = 1 * float(s[1]) + 1 * float(s[2]) + 3.5 * float(s[3]) + 0.5 * float(s[4])
    print(seq, disc)

meh=['0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-72.9722,-0.62781,0.423458, 3.081625,UNCLEAVED', '0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-76.1073,0.05851,0.430077, 7.954625,UNCLEAVED', '0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-74.6781,-3.86377,0.159481, -1.278125,UNCLEAVED', '0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-74.2416,-0.37881,0.408195, -2.81975,UNCLEAVED', '0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-72.6746,-3.19253,0.378113, -0.856499999999998,UNCLEAVED']
meh = [i.split(',')[-5:] for i in meh]
meh = [i[-1]+i[:-1] for i in meh]
meh = [[i[-1]]+i[:-1] for i in meh]
for s in meh:
    seq = s[0]
    disc = 1 * float(s[1]) + 1 * float(s[2]) + 3.5 * float(s[3]) + 0.5 * float(s[4])
    print(seq, disc)
for l in seqs:
    sqcol = [int(i)  for i in l.split(b',')[:100]]
    s = ''
    for i in range(20,101,20):
        frame = sqcol[i-20:i]
        s+=aas[frame.index(1)]
    if s in ident_seqs:
        print(l)
seqs
'0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,-75.2536,-1.82782,0.423983, -3.00625,CLEAVED\n'.split(',')[100]
stat_max_min = {'CLEAVED':[-80,'',-80,''], 'UNCLEAVED':[-80,'',-80,''], 'MIDDLE':[-80,'',-80,'']}
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
stat_max_min
for l in seqs:
    sqcol = [int(i)  for i in l.split(b',')[:100]]
    s = ''
    for i in range(20,101,20):
        frame = sqcol[i-20:i]
        s+=aas[frame.index(1)]
    if s in ident_seqs+other_seqs:
        print(l)
for l in seqs:
    sqcol = [int(i)  for i in l.split(b',')[:100]]
    s = ''
    for i in range(20,101,20):
        frame = sqcol[i-20:i]
        s+=aas[frame.index(1)]
    if s in ident_seqs+other_seqs:
        print(s, l)
len(seqs)