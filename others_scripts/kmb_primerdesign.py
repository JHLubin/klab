#!/usr/bin/python
import os, sys
sys.path.append('/home/kmb413/biopython')
sys.path.append('/home/kmb413/mechanize-0.2.5')

from Bio import SeqIO
from Bio.Seq import Seq

from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from rosetta import *
init('-ignore_unrecognized_res')

import re
import mechanize

import argparse

codon_usage_table = { 'A':'GCG', 'C':'TGC', 'D':'GAT', 'E':'GAA', 'F':'TTT', 'G':'GGC', 'H':'CAT', 'I':'ATT', 'K':'AAA', 'L':'CTG', 'M':'ATG', 'N':'AAC', 'P':'CCG', 'Q':'CAG', 'R':'CGT', 'S':'AGC', 'T':'ACC', 'V':'GTG', 'W':'TGG', 'Y':'TAT' }

def get_aa_seq(pdb_file):
    pose = pose_from_pdb(pdb_file)
    chains = pose.split_by_chain()
    chain_sequences = []
    for chain in chains:
        aa_seq = chain.sequence()
        chain_sequences.append(aa_seq)
    return chain_sequences

def align_aa_seq(temp_seq, aa_seq):
    alignments = pairwise2.align.globalxx(temp_seq, aa_seq)
    return alignments

def find_gaps(aligned_seq):
    return [i for i, ltr in enumerate(aligned_seq) if ltr == '-' ]

def parse_primers(primer_html):
    p1 = re.sub('<br />', '', primer_html)
    p2 = p1.split(';')
    p3 = [x for x in p2 if x != '' and '*' not in x]
    #print(p3)
    parsed_primer = { 'forward' : p3[1], 'reverse' : p3[2], 'gc_content' : p3[3], 'location' : p3[4], 'melt_temp' : p3[5], 'mismatched_bases' : p3[6], 'length' : p3[7], 'mutation_type' : p3[8], 'length_5p' : p3[9], 'forward_mass' : p3[10], 'length_3p' : p3[11], 'reverse_mass' : p3[12][:p3[12].index('<p />')] }

    return parsed_primer

def parse_primer_params( primer_params_file ):
    
    with open(primer_params_file) as pfile:
        primer_params_lines = pfile.readlines()

    primer_param_dict = {}
    for primer_line in primer_params_lines:
        if primer_line[0] != '#':
            split_primer_line = primer_line.split('=')
            key = "".join(split_primer_line[0].split())
            value = "".join(split_primer_line[1].split())
            if key == 'ends_in_GC' or key == 'mut_at_center':
                if value == 'True' or value == 'T' or value == 'true' or value == 't' or value == '1':
                    value = True
                else:
                    value = False
            primer_param_dict[key] = value

    return primer_param_dict

def get_primers(template_sequence, mutation_code, primer_params):
    br = mechanize.Browser()
    url = "http://www.bioinformatics.org/primerx/cgi-bin/DNA_1.cgi"
    response = br.open(url)

    br.form = list(br.forms())[0]
    protocol_control = br.form.find_control("protocol")
    br[protocol_control.name] = ["QuikChange"]
    
    seq_control = br.form.find_control("sequence_textarea")
    br[seq_control.name] = template_sequence

    code_control = br.form.find_control("code")
    br[code_control.name] = mutation_code

    response = br.submit()
    br.form = list(br.forms())[0]
    br.set_all_readonly(False)

    min_tm_control = br.form.find_control("min_Tm")
    max_tm_control = br.form.find_control("max_Tm")
    min_gc_control = br.form.find_control("min_GC")
    max_gc_control = br.form.find_control("max_GC")
    min_len_control = br.form.find_control("min_length")
    max_len_control = br.form.find_control("max_length")
    min_5pflank_control = br.form.find_control("min_5p_flank")
    max_5pflank_control = br.form.find_control("max_5p_flank")
    min_3pflank_control = br.form.find_control("min_3p_flank")
    max_3pflank_control = br.form.find_control("max_3p_flank")
    ends_in_GC_control = br.form.find_control("ends_in_GC")
    mut_at_center_control = br.form.find_control("mut_at_center")

    if primer_params.has_key('min_Tm'):
        br[min_tm_control.name] = primer_params['min_Tm']
    if primer_params.has_key('max_Tm'):
        br[max_tm_control.name] = primer_params['max_Tm']
    if primer_params.has_key('min_GC'):
        br[min_gc_control.name] = primer_params['min_GC']
    if primer_params.has_key('max_GC'):
        br[max_gc_control.name] = primer_params['max_GC']
    if primer_params.has_key('min_length'):
        br[min_len_control.name] = primer_params['min_length']
    if primer_params.has_key('max_length'):
        br[max_len_control.name] = primer_params['max_length']
    if primer_params.has_key('min_5p_flank'):
        br[min_5pflank_control.name] = primer_params['min_5p_flank']
    if primer_params.has_key('max_5p_flank'):
        br[max_5pflank_control.name] = primer_params['max_5p_flank']
    if primer_params.has_key('min_3p_flank'):
        br[min_3pflank_control.name] = primer_params['min_3p_flank']
    if primer_params.has_key('max_3p_flank'):
        br[max_3pflank_control.name] = primer_params['max_3p_flank']
    if primer_params.has_key('ends_in_GC'):
        br.find_control("ends_in_GC").items[0].selected=primer_params['ends_in_GC']  ##True/False
    if primer_params.has_key('mut_at_center'):
        br.find_control("mut_at_center").items[0].selected=primer_params['mut_at_center'] ##True/False

    print('Primer Parameters:')
    print("min_Tm = {mntm}\tmax_Tm = {mxtm}\tmin_length = {mnlen}\t\tmax_length = {mxlen}".format(\
    mntm=min_tm_control.value, mxtm=max_tm_control.value, mnlen=min_len_control.value, mxlen=max_len_control.value))
    print("min_gc = {mngc}\tmax_gc = {mxgc}\tends_in_GC = {endgc}\tmut_at_center = {mtcen}".format(\
    mngc=min_gc_control.value, mxgc=max_gc_control.value, endgc=ends_in_GC_control.value, mtcen=mut_at_center_control.value))
    print("min_5p_flank = {mn5p}\tmax_5p_flank = {mx5p}\tmin_3p_flank = {mn3p}\tmax_3p_flank = {mx3p}".format(\
    mn5p=min_5pflank_control.value, mx5p=max_5pflank_control.value, mn3p=min_3pflank_control.value, mx3p=max_3pflank_control.value))

    response = br.submit()

    html_string = response.get_data()
    nospace = re.sub('&nbsp', '', html_string)
    primer_chunks = nospace.split('Primer pair')
    count = 0
    
    parsed_primers_list = []
    for primer in primer_chunks[1:]:
        count += 1
        parsed_primers = parse_primers(primer)
        parsed_primers_list.append(parsed_primers)

    return parsed_primers_list

def write_primers_to_outfile(outfile, list_of_primers):
    primer_count = 0
    with open(outfile,'w') as ofile:
        ofile.write("Primer Pair No.\tForward (5'->3')\tReverse (5'->3')\t\
Length (bp)\tTm (C)\tGC Content (%)\tMismatched Bases\t5p Length (bp)\t\
3p Length (bp)\tForward MW (Da)\tReverse MW (Da)\tLocation\tMutation Type\n")

        for primer in list_of_primers:
            primer_count += 1
            ofile.write("{count}\t{fwd}\
            \t{rev}\t{ln}\t{mt}\t{gc}\t\
            {mb}\t{l5p}\t{l3p}\t{fmass}\t\
            {rmass}\t{loc}\t{mtype}\n".format(count=primer_count, \
            fwd=primer['forward'].split(': ')[1].split()[1],\
            rev=primer['reverse'].split(': ')[1].split()[1],\
            ln=primer['length'].split(': ')[1].split()[0],\
            mt=primer['melt_temp'].split(': ')[1][:-2],\
            gc=primer['gc_content'].split(': ')[1][:-1],\
            mb=primer['mismatched_bases'].split(': ')[1],\
            l5p=primer['length_5p'].split(': ')[1].split()[0],\
            l3p=primer['length_3p'].split(': ')[1].split()[0],\
            fmass=primer['forward_mass'].split(': ')[1].split()[0],\
            rmass=primer['reverse_mass'].split(': ')[1].split()[0],\
            loc=primer['location'].split(': ')[1],\
            mtype=primer['mutation_type'].split(': ')[1] ))

################################################################
################################################################
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--template_fasta", action="store")
    parser.add_argument("-o", "--outfile", action="store")
    parser.add_argument("-l", "--design_list", action="store")
    parser.add_argument("-d", "--design", action="store")
    parser.add_argument("-p", "--primer_params", action="store")
    parser.add_argument("-b", "--prepend_num", action="store")
    args = parser.parse_args()

    ## Parse template FASTA seq
    fasta_file = args.template_fasta
    if not fasta_file:
        sys.exit("\033[91mERROR: No FASTA file input.\033[0m")

    handle = open(fasta_file)
    for template_dna in SeqIO.parse(handle, "fasta"):
        print(template_dna.id)
        #TODO Save chains in list?
    handle.close()
    print("===== Template DNA Sequence:=====\n" + str(template_dna.seq))

    ## Don't translate the prepended part
    prepended_bps = args.prepend_num
    if not prepended_bps:
        prepended_bps = 0
    else:
        prepended_bps = int(prepended_bps)

    prepended_seq = str(template_dna.seq[:prepended_bps])
    print("\n\tPrepended Sequence: " + prepended_seq)
    template_protein_dna = template_dna.seq[prepended_bps:]
    print("\n\tTemplate Protein DNA: " + template_protein_dna)

    ## Translate to Protein Seq
    template_aa = template_protein_dna.translate()
    print("\n=====Template Amino Acid Sequence:=====\n" + str(template_aa))

    ## Gather Designs 
    ## TODO Ignore commented-out designs
    list_of_designs = []
    if args.design:
        list_of_designs.append(args.design)
    if args.design_list:
        with open(args.design_list) as dfile:
            design_lines = dfile.readlines()
        for dline in design_lines:
            list_of_designs.append(dline.rstrip())

    ## Parse user input
    for design in list_of_designs:
        ##### Parse PDB Designs #####
        pdb_designs = {}
        if '.pdb' in design:
            list_of_designs.pop(list_of_designs.index(design))
            design_chains_seqs = get_aa_seq(design)
            for design_aa_seq in design_chains_seqs: ###ONLY CHAIN A FOR NOW
                alignments = align_aa_seq(template_aa, design_aa_seq)
                #print(alignments[0])
                gap_indices = find_gaps( alignments[0][0] )
                #print(gap_indices)
                
                mut_ranges = {}
                for gap_index in gap_indices:
                    #print("GapIndex = "+str(gap_index))
                    adjacent = False
                    for key in mut_ranges.keys():
                        if gap_index == key + mut_ranges[key]:
                            mut_ranges[key] += 1
                            adjacent = True
                    if not adjacent:
                        #print("Not adjacent")
                        mut_ranges[gap_index] = 1
                #print mut_ranges
                #print(sorted(mut_ranges.keys()))

                dashes = 0
                for mut in sorted(mut_ranges.keys()):
                    dashes += mut_ranges[mut]
                    list_of_designs.append("{native}{index}{design}".format(\
                    native=alignments[0][0][mut-mut_ranges[mut]:mut],\
                    index=mut-dashes+1,\
                    design=alignments[0][1][ mut:mut+mut_ranges[mut]]\
                    ))
        ##### /Parse PDB Designs #####

    print("\nDesigns:")
    print(", ".join(list_of_designs))

    if len(list_of_designs) == 0:
        sys.exit("\033[93mWARNING: No designs found.\033[0m")

    for design in list_of_designs:
        ### Find DNA seq index of mutations in template
        print('\n===== Finding DNA Range for {} ====\n'.format(design))
        match = re.match(r"([a-z]+)([0-9]+)([a-z]+)", design, re.I)
        if match:
            mutation = match.groups()
        else:
            sys.exit("\033[91mERROR: Mutation is not in {Char(s)}{Int}{Char(s)} format.\033[0m")
        #print(mutation)

        template_aa_index_begin = int(mutation[1])-1
        template_aa_index_end = template_aa_index_begin+len(mutation[0])

        if template_aa[ template_aa_index_begin : template_aa_index_end ] != mutation[0]:
            sys.exit("\033[91mERROR: Amino Acid #{index} is {aa}, not {user_aa}.\033[0m\n...{t1}\033[91m{t}\033[0m{t2}...\n".format(\
            index=template_aa_index_begin+1, aa=template_aa[ template_aa_index_begin : template_aa_index_end ], user_aa=mutation[0], \
            t1=template_aa[ template_aa_index_begin - 10 : template_aa_index_begin ], t=template_aa[ template_aa_index_begin: template_aa_index_end ], \
            t2=template_aa[ template_aa_index_end+1:template_aa_index_end + 11] ))
        
        mutation_design_codons = ''
        mut_index = 0
        while mut_index < len(mutation[2]):
            
            if mutation[2][mut_index].islower():
                mutation_design_codon = mutation[2][mut_index:mut_index+3]
                mut_index += 3
            elif mutation[2][mut_index] not in codon_usage_table.keys():
                sys.exit("\033[91mERROR: {} is not a valid amino acid.\n\033[0m".format(mutation[2][mut_index]))
            else:
                mutation_design_codon = codon_usage_table[mutation[2][mut_index]]
                mut_index += 1
            
            mutation_design_codons += mutation_design_codon

        template_dna_index = template_aa_index_begin*3

        ## Get template seq from -60 of first mutation to +60 of last mutation
        left_addition = 60
        right_addition = 60
        prepend = ''
        if template_dna_index < left_addition:
            if prepended_bps < left_addition - template_dna_index:
                prepend = prepended_seq
            else:
                prepend_start = -1*(left_addition-template_dna_index)
                prepend = prepended_seq[ prepend_start: ]
            left_addition = template_dna_index

        if len(template_dna) < right_addition:
            right_addition = len(template_dna) - right_addition
        
        mutation_dna_range_left = prepend + str(template_protein_dna)[template_dna_index - left_addition : template_dna_index]
        mutation_dna_range_target = str(template_protein_dna)[template_dna_index : template_dna_index+(3*len(mutation[0]))]
        mutation_dna_range_right = str(template_protein_dna)[template_dna_index+(3*len(mutation[0])) : template_dna_index + 3 + right_addition]

        template_range_pretty = mutation_dna_range_left + '\033[91m' + mutation_dna_range_target + '\033[0m' + mutation_dna_range_right
        template_range = mutation_dna_range_left + mutation_dna_range_target + mutation_dna_range_right
        
        mutation_range_pretty = mutation_dna_range_left + '\033[91m' + mutation_design_codons + '\033[0m' + mutation_dna_range_right
        mutation_range = mutation_dna_range_left + mutation_design_codons + mutation_dna_range_right
        
        mutation_code = str(template_protein_dna[template_dna_index : template_dna_index+(3*len(mutation[0]))]) + str(left_addition+1+len(prepend)) + mutation_design_codons
        
        print("     TEMPLATE\t" + template_range_pretty)
        print("     MUTATION\t" + mutation_range_pretty)
        print("MUTATION CODE\t" + mutation_code + '\n')

        primer_params = {}
        primer_params_file = args.primer_params
        if primer_params_file:
            primer_params = parse_primer_params( primer_params_file )

        list_of_primers = get_primers(template_range, mutation_code, primer_params)

        if len(list_of_primers) != 0:
            outfile = args.outfile
            if not outfile:
                outfile = "primers_{fas}_{des}.out".format(fas=re.sub('.fasta','',fasta_file), des=design)
            
            plural = ''
            if len(list_of_primers) > 1:
                plural = 's'
            
            print("\n\033[92m{l} primer pair{s} found.".format(l=len(list_of_primers), s=plural ) )
            print("Writing primers to outfile:  {}\033[0m".format(outfile))
            write_primers_to_outfile(outfile, list_of_primers)
        else:
            print("\n\033[91mNo primers found.\n\033[93mPossible reasons:\n\t1. The target site for mutation in the DNA sequence you provided was too near one of the ends. Try submitting a longer template sequence with base pairs added to both sides.\n\t2. The DNA sequence you provided had a GC content that was particularly very high or very low. Try adjusting the minimum and maximum melting temperature and %GC content accordingly.\n\t3. The mutation you entered affected too many base pairs. This leads to a high mismatch, which consequently lowers the melting temperature. Try entering a more minor mutation, or lowering the minimum melting temperature.\n\t4. Your constraints were too stringent. Try making them a little more flexible.\033[0m")


if __name__ == "__main__":
    main()
    
