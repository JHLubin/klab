hide all
show cartoon, all
select triad, res 220+250+328
select hole, res 325-328
select pep, chain B
show sticks, triad and not (name N+C+O)
show sticks, hole and (name N+C+CA+O)
hide cartoon, hole and not triad
show lines, pep and not (name N+C+O)
show lines, pep and res 8
hide everything, e. H
color cyan, hole
color magenta, triad
color yellow, pep
util.cnc
deselect