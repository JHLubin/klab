#!/usr/bin/perl

open (F,"$ARGV[0]"); #list of pdb names
$alignto = $ARGV[1]; #align to this PDB
open (G, ">$ARGV[2]"); #output pml file

print G "
#
# example PyMOL script, save as run with pymol -cq align.pml
# needs cealign plugin in pymol and numpy installed
# add the following two lines to .pymolrc 
# run /path/to/cealign/directory/qkabsch.py  
# run /path/to/cealign/directory/cealign.py  

from pymol import cmd

";

while(<F>){
chomp($_);
print G "cmd.load(\"$_.pdb\")\n";
}
print G "alignto(\"$alignto")\n";

seek(F,0,0);

while(<F>){
chomp($_);
print G "cmd.save(\"$_.pdb\", \"$_\")\n";
}
