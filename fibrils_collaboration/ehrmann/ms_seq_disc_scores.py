from numpy import mean
from math import exp

# input aa_dist data from Excel as str_set
l_set=[i.split() for i in str_set]
d_set=[]

for l in l_set:
    d={}
    for i in range(1,len(l),2):
        d[l[i]]=float(l[i+1])
    d_set.append(d)

outs=[]
for i in d_set:
    l = []
    for j in sorted(d_set[0].keys()):
        l.append(i[j])
    outs.append(l)

# output list to put back into Excel
outs_str = ['\t'.join([str(i) for i in l]) for l in outs]
for i in outs_str:
    print i

scoring = [[mean([outs[i][j],outs[i+12][j]]) for j in range(20)] for i in range(12)]
inds = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
a_syn_seq = 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'

for i in range(5,len(a_syn_seq)-6):
frame = a_syn_seq[i-5:i+7]
score = 0
for n, r in enumerate(frame):
score += scoring[n][inds.index(r)]
print i, frame, score

for i in range(5,len(a_syn_seq)-6):
    frame = a_syn_seq[i-5:i+7]
    score = 1
    for n, r in enumerate(frame):
        score *= scoring[n][inds.index(r)]
    print i, frame, score

"""
trim_scores = [[] for i in range(12)]
for i in range(12):
for j in range(20):
if scoring[i][j]>0.1:
trim_scores[i].append(scoring[i][j])

btot = 0
for a in range(len(trim_scores[0])):
for b in range(len(trim_scores[1])):
for c in range(len(trim_scores[2])):
for d in range(len(trim_scores[3])):
for e in range(len(trim_scores[4])):
for f in range(len(trim_scores[5])):
for g in range(len(trim_scores[6])):
for h in range(len(trim_scores[7])):
for i in range(len(trim_scores[8])):
for j in range(len(trim_scores[9])):
for k in range(len(trim_scores[10])):
for l in range(len(trim_scores[11])):
tot = trim_scores[0][a] + trim_scores[1][b] + trim_scores[2][c] + trim_scores[3][d] + trim_scores[4][e] + trim_scores[5][f] + trim_scores[6][g] + trim_scores[7][h] + trim_scores[8][i] + trim_scores[9][j] + trim_scores[10][k] + trim_scores[11][l] 
btot+=exp(-tot)
"""

bw = [[exp(-i) for i in ss] for ss in scoring]
bol_score = [[i/sum(b) for i in b] for b in bw]

for i in range(5,len(a_syn_seq)-6):
    frame = a_syn_seq[i-5:i+7]
    score = 1
    for n, r in enumerate(frame):
        score *= bol_score[n][inds.index(r)]
    print i, frame, score