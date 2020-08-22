#!/bin/bash

rm job.sub

export LD_LIBRARY_PATH=${GARFIELD_HOME}/Install/lib:$LD_LIBRARY_PATH

numberOfJobs=5
phraseIn='rootFiles-Ar-iC4H10-model10-fe-signal-noibf-350-430-530-830-950-0.root'
phraseOut='rootFiles-Ar-iC4H10-model10-fe-signal-noibf-350-430-530-830-950-$(ProcId).root'
hv='350 430 530 830 950'
touch job.sub

echo 'executable      = signal' >> job.sub
echo 'arguments       = '$hv 0 >> job.sub
echo 'output          = output/ex.$(ClusterId).$(ProcId).out' >> job.sub
echo 'error           = error/ex.$(ClusterId).$(ProcId).err' >> job.sub
echo 'log             = log/ex.$(ClusterId).log' >> job.sub
#echo 'transfer_output_files   =' "$phrase" >> job.sub
echo 'transfer_output_remaps  = "'$phraseIn'"="'$phraseOut'" ' >> job.sub
echo 'queue' $numberOfJobs >> job.sub

make
condor_submit job.sub
