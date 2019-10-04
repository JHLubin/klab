hide all
show cartoon, all
select triad, res 220+250+328
select hole, res 325-328
select pep, chain B
#select mutable, res 203+204+205+218+221+281+283+284+285+304+307+308+309+316+323+324+325+326+327+329+330+343+344+345+346+347+348+351+352+353
select mutable, res 203+285+304+307+308+309+316+323+324+325+326+327+343+344+345+346+347+348+351+353
#select mutable, res 185+186+187+188+190+199+200+201+202+203+204+205+217+218+221+225+283+284+285+286+287+288+302+307+308+309+316+323+324+325+326+327+329+343+344+345+346+347+348+351+352+353

show sticks, triad and not (name N+C+O)
show lines, pep and not (name N+C+O)
show lines, mutable and not (name N+C+O)
# show sticks, hole and (name N+C+CA+O)
# hide cartoon, hole and not triad
#show lines, pep and res 8
hide everything, e. H

#color cyan, hole
color cyan
color magenta, triad
color yellow, pep
color purpleblue, mutable
color green, cat_relaxed
util.cnc

label (cat_relaxed and mutable and n. CA), "(%s, %s)" % (resi, resn)

deselect
set seq_view, 1



