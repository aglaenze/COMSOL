#!/bin/bash

cd $COMSOL/build
cmake $COMSOL
wait
make

mv avalanche $COMSOL
mv plotField $COMSOL
mv spectrumFe55 $COMSOL
mv signal $COMSOL
