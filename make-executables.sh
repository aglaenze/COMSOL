#!/bin/bash
source $GARFIELD_HOME/install/share/Garfield/setupGarfield.sh

if ! [ -d $COMSOL/build ]
then
mkdir $COMSOL/build
fi

cd $COMSOL/build
cmake $COMSOL
wait
make

mv avalanche $COMSOL
mv plotField $COMSOL
mv spectrumFe55 $COMSOL
mv signal $COMSOL

cd $COMSOL
