precision 3
method ml
set szr 0
set sv 0
set st0 0.2
set d 0
set p 0
depends v frequency
format subjectID which_task word_nonword RESPONSE frequency correct TIME
load *English_final.dat
log parms_E_hf_lf_final.txt