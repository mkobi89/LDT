precision 3
method ml
set szr 0
set sv 0
set st0 0.2
set d 0
set p 0
depends v condition
format subjectID which_task word_nonword RESPONSE condition correct TIME
load *Switch_final.dat
log parms_SW_condition_final.txt