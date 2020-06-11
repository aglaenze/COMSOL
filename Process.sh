#!/bin/bash

hvList={
for var in "$@"
do
	echo $var
	hvList=$hvList$var,
done
hvList=${hvList%?}}
#echo $hvList

gasName="Ar-iC4H10"
modelNum=1


#root -l -q "AddSignalTrees.C($modelNum, \"$gasName\", $hvList)"
#root -l -q "Convolute.C($modelNum, \"$gasName\", $hvList)"
#root -l -q "Analyse.C($modelNum, \"$gasName\", $hvList)"


# boucle for uniquement pour modeles 1, 16 et 17
for ((k=0;k<2;k++)); do
hv1=$((340+$k*20))
hv2=$((540+$k*20))
hvList="{$hv1, $hv2}"
#echo $hvList
#
#root -l -q "AddSignalTrees.C($modelNum, \"$gasName\", $hvList)"
#root -l -q "Convolute.C($modelNum, \"$gasName\", $hvList)"
root -l -q "Analyse.C($modelNum, \"$gasName\", $hvList)"
done

