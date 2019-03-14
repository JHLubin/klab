select relaxed, *relaxed*
color green, relaxed
select designed, *designed*
color cyan, designed
select protease, chain A and designed
select peptide, chain B
select cat_res, res 72+96+154
select des_res, res 58+70+138+147+150-153+169-177+183
select pep_recognition, res 198-203 and designed
select shell_res, res 51-57+59+69+71+73-75+121+123-124+134-135+137+139+146+148-149+155-156+158+168+178-182 and designed
color yellow, peptide and designed
color magenta, cat_res
color purpleblue, des_res and designed
util.cnc

hide
show cartoon
show sticks, cat_res 
hide sticks, name C+N+O
show lines, peptide des_res
hide lines, name C+N+O
set line_width, 3
hide everything, elem H
set seq_view, 1
deselect