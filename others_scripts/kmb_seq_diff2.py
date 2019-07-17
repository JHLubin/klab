# rlc color_b.py version 6.0

import colorsys,sys
from pymol import cmd

# main function called from within PyMOL
def seq_diff(obj1, obj2, doColor=False, doColor2=False):
    if cmd.count_atoms(obj1+' and name CA') == 0:
        print('%s is empty' % obj1)
        return
    if cmd.count_atoms(obj2+' and name CA') == 0:
        print('%s is empty' % obj2)
        return

    if cmd.count_atoms(obj1+' and name CA') != cmd.count_atoms(obj2+' and name CA'):
        print('Selections have different number of residues.')
        return

    stored.resn1, stored.resn2 = [], []
    stored.resi1, stored.resi2 = [], []
    stored.chain1, stored.chain2 = [], []

    cmd.iterate(obj1 + ' and name CA', 'stored.resn1.append(resn)')
    cmd.iterate(obj2 + ' and name CA', 'stored.resn2.append(resn)')

    cmd.iterate(obj1 + ' and name CA', 'stored.resi1.append(resi)')
    cmd.iterate(obj2 + ' and name CA', 'stored.resi2.append(resi)')

    cmd.iterate(obj1 + ' and name CA', 'stored.chain1.append(chain)')
    cmd.iterate(obj2 + ' and name CA', 'stored.chain2.append(chain)')

    #print len(stored.resn1)

    mutant_selection = []
    reg_selection = []

    for n1, n2, i1, i2, c1, c2 in zip(stored.resn1, stored.resn2, stored.resi1, stored.resi2, stored.chain1, stored.chain2):
        if c1 == '':
            c1 = '""'
        if c2 == '':
            c2 = '""'
        #print 'Resi1: '+i1
        #print 'Resi2: '+i2
        if n1 != n2:
            #print 'Mutation at resi: '+i2
            mutant_selection.append('%s and resi %s and chain %s' % (obj2, i2, c2))
            reg_selection.append('%s and resi %s and chain %s' % (obj1, i1, c1))

    if mutant_selection == '':
        print('No mutations found')
        #raise CmdException

    total_selection = ''
    total_color = ''
    for m in mutant_selection:
        total_selection += '('+m+') or '
        total_color += '('+m+' and name C*) or '

    total_selection = total_selection[:-4]
    total_color = total_color[:-4]

    cmd.select("%s_diff" % (obj2), total_selection )

    if doColor:
        cmd.color('%s' % doColor,total_color)


    second_selection = ''
    second_color = ''
    for n in reg_selection:
        second_selection += '('+n+') or '
        second_color += '('+n+' and name C*) or '

    second_selection = second_selection[:-4]
    second_color = second_color[:-4]

    cmd.select("%s_diff" % (obj1), second_selection )

    if doColor2:
        cmd.color('%s' % doColor2,second_color)
        
    #cmd.deselect()

cmd.extend("seq_diff",seq_diff)

cmd.auto_arg[0]['seq_diff'] = cmd.auto_arg[0]['align']
cmd.auto_arg[1]['seq_diff'] = cmd.auto_arg[1]['align']

# vi: ts=1:sw=4:smarttab:expandtab
