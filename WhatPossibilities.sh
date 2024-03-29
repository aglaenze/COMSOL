#!/bin/bash

maxModel=0

for ((k=1;k<100;k++)); do
#dataFolder="COMSOL_data/model$k"
dataFolder="/Volumes/Transcend/COMSOL_data/model$k"
if [ -d $dataFolder ]
then
maxModel=$k
#else
#break
fi
done


# function that will display electic field configuration in a given rep
display()
{
#repertoire="$COMSOL/COMSOL_data/model$num"
repertoire="/Volumes/Transcend/COMSOL_data/model$num"
if [ -d $repertoire ]
then
echo ""
echo "Possible electric field configurations for MODEL $num"
cd $repertoire
#y=$(eval "ls ewfield*.txt | wc -l")
#echo "$y possible electric field configuration"
for entry in `ls ewfield*.txt`; do
    #echo $entry
    echo "$entry" | tr 'ewfield' ' ' | tr '-' ' ' | tr '.txt' ' '
done
#else
#echo "repository COMSOL_data/model$num does not exist (yet)"
#echo ""
#exit
fi
}

# Here it starts
echo ""
# Checks what model number(s) we're interested in
if [ -z $1 ]
then
        echo "Possible electric field configurations for all models"
        echo "If you want to know for one specific model, write ./WhatPossibilities.sh \$modelNum"
else
if [ $1 -lt 1 ] || [ $maxModel -lt $1 ]
then
echo "Model number has to be an integer between 1 and $maxModel"
echo ""
exit
fi
fi

# And then call the display function
if [ -z $1 ]
then
for ((k=1;k<=$maxModel;k++)); do
num=$k
display
done
else
num=$1
display
fi


echo ""











