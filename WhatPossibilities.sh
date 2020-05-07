#!/bin/bash

# function that will display electic field configuration in na given rep
display()
{
repertoire="$COMSOL/COMSOL_data/model$num"
if [ ! -d $repertoire ]
then
echo "repository COMSOL_data/model$num does not exist (yet)"
#echo ""
#exit
else
cd $repertoire
#y=$(eval "ls ewfield*.txt | wc -l")
#echo "$y possible electric field configuration"
for entry in `ls ewfield*.txt`; do
    #echo $entry
    echo "$entry" | tr 'ewfield' ' ' | tr '-' ' ' | tr '.txt' ' '
done
fi
}

# Here it starts
echo ""
# Checks what model number(s) we're interested in
if [ -z $1 ]
then
        echo "Possible electric field configurations for all models"
        echo "(If you want to know for one specific model, write ./WhatPossibilities.sh \$modelNum)"
else
if [ 1 -le $1 ] && [ $1 -le 14 ]
then
echo "Possible electric field configurations for model $1"
else
echo "Model number has to be an integer between 1 and 14"
echo ""
exit
fi
fi

# And then call the display function
if [ -z $1 ]
then
for ((k=1;k<=14;k++)); do
echo ""
echo "Model $k"
num=$k
display
done
else
num=$1
display $1
fi


echo ""











