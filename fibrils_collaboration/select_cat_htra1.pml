hide all
show cartoon, all
select triad, res 220+250+328
select hole, res 325-328
select pep, chain B
select mutable, res 203-205+218+221+283+285+307-307+316+323-327+329+343-348+351-353

show sticks, triad and not (name N+C+O)
# show sticks, hole and (name N+C+CA+O)
# hide cartoon, hole and not triad
show lines, pep and not (name N+C+O)
show lines, pep and res 8
show lines, mutable and not (name N+C+O)
hide everything, e. H

color cyan, hole
color magenta, triad
color yellow, pep
color purpleblue, mutable
util.cnc

deselect
set seq_view, 1



