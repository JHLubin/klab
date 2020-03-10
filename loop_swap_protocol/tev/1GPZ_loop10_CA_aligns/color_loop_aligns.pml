from pymol import cmd

#load ..\tev.pdb
#loadall *.pdb

color cyan
color green, tev

models = cmd.get_names('objects')
models.pop(0)
color = 'magenta'

commands = []
for m in models: commands.append('({} and res {}-{}+{}-{})'.format(m, str(int(m[14:17]) + 317), str(int(m[18:21]) + 317), str(int(m[23:26]) + 317), str(int(m[27:30]) + 317))) 

for c in commands:cmd.color('%s' % color,c) 

cmd.create('best','aligned_loop_N303-303_C321-321')
cmd.color('%s' % 'white','best')
cmd.color('%s' % 'orange','best'+' and res 620+638') 

util.cnc

hide cartoon, not tev
hide cartoon, tev and res 139-157
show sticks, tev and res 138-158 and n. N+CA+C
show lines, not tev and n. N+CA+C
hide lines, best
show sticks, best and n. N+CA+C