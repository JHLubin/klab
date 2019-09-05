select catalytic, res 33+62+169+174
#select mutable, res 61+65+88+89+114
#select mutable, res 29+30+63+65+68+72+73+76+77+79-84+112+117+136-138+166+169+170+173+175+176+181+184+185+188+189+192
select mutable, res 29+63+65+68+72+73+81-83+112+117+137+138+166+169+170+173+175+176+184+185+188+189

select ligand, chain B

show sticks, catalytic and not (name N+C+O)
show lines, mutable and not (name N+C+O)
hide everything, e. H

color cyan
color purpleblue, mutable
color magenta, catalytic
color yellow, ligand
util.cnc()