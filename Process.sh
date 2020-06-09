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
modelNum=8


root -l -q "AddSignalTrees.C($modelNum, \"$gasName\", $hvList)"
root -l -q "Convolute.C($modelNum, \"$gasName\", $hvList)"
root -l -q "Analyse.C($modelNum, \"$gasName\", $hvList)"

