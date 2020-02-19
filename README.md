# COMSOL
Garfield project for IBF simulation in a Micromegas detector

Don't forget to change in the source code the unit length (in ComponentComsol.cc) --> it has to be in um and not m, replace 100.0 with 1.e-4



To draw the gain:
1) ./gain
2) root -l -q Convolute.C
3) root -l -q Analyse.C

To compute the IBF:
1) ./ibf
2) root -l -q Analyse.C
