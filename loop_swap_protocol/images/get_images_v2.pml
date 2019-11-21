set orthoscopic, on
set ray_trace_color, black
set ray_trace_mode, 1
load tev.pdb
select start_loop, tev and res 142-153
select start_loop_in, tev and res 143-152
select pep, tev and chain C
color gray60, tev and chain A
color yellow, tev and chain C
color cyan, start_loop
set cartoon_transparency, 0.4, tev and not (pep or start_loop)
create orig_loop, start_loop
set cartoon_transparency, 0.7, orig_loop
load aligned_pdbs/1KY9.pdb; create 1KY9_10,  1KY9 and res 201-212; delete 1KY9; pair_fit 1KY9_10///201+212/C+CA+N,  tev///142+153/C+CA+N;  load aligned_pdbs/1GPZ.pdb; create 1GPZ_10,  1GPZ and res 623-639; delete 1GPZ; pair_fit 1GPZ_10///623+639/C+CA+N,  tev///142+153/C+CA+N;  load aligned_pdbs/5EDM.pdb; create 5EDM_10,  5EDM and res 494-513; delete 5EDM; pair_fit 5EDM_10///494+513/C+CA+N,  tev///142+153/C+CA+N; 
color green, not tev or orig_loop
cmd.set_view((0.8107369542121887, 0.5591747164726257, 0.1731332391500473, 0.5358902812004089, -0.8280136585235596, 0.16482844948768616, 0.23550602793693542, -0.040846485644578934, -0.970989465713501, -0.006182016339153051, 0.00023031607270240784, -19.096799850463867, 53.47949981689453, 65.19629669189453, 28.549774169921875, -146.4417266845703, 185.23509216308594, 60.0))
# Original loop
hide everything, not tev
ray 1000, 1000
# 1KY9
show cartoon, orig_loop
hide everything, start_loop_in
show cartoon, 1KY9_10
ray 1000, 1000
# 1GPZ
hide everything, 1KY9_10
show cartoon, 1GPZ_10
ray 1000, 1000
# 5EDM
hide everything, 1GPZ_10
show cartoon, 5EDM_10
ray 1000, 1000

