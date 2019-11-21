load tev.pdb
load aligned_pdbs/2M9P.pdb; create 2M9P_10,  2M9P and res 187-198; delete 2M9P; pair_fit 2M9P_10///187+198/C+CA+N,  tev///142+153/C+CA+N
load aligned_pdbs/1EQ9.pdb; create 1EQ9_10,  1EQ9 and res 183-198; delete 1EQ9; pair_fit 1EQ9_10///183+198/C+CA+N,  tev///142+154/C+CA+N; 
select start_loop, tev and res 142-153
select start_loop_in, tev and res 143-152
color gray60, tev and chain A
color yellow, tev and chain C
color cyan, start_loop
color green, not tev
set orthoscopic, on
cmd.set_view((0.6336492896080017, 0.6316966414451599, 0.4465373456478119, 0.606505811214447, -0.7639744877815247, 0.22010622918605804, 0.4801729619503021, 0.13136525452136993, -0.867253839969635, -0.006182016339153051, 0.00023031607270240784, -18.032573699951172, 53.47949981689453, 65.19629669189453, 28.549774169921875, -6.318345546722412, 42.9832763671875, 60.0))
set ray_trace_color, black
set ray_trace_mode, 1
# Original loop
hide everything, 2M9P_10 or 1EQ9_10
ray
# 2M9P
show cartoon, 2M9P_10
hide everything, start_loop_in
ray
# 1EQ9
hide everything, 2M9P_10
show cartoon, 1EQ9_10
ray

