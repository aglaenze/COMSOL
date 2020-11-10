#!/bin/bash

numberOfJobs=400

# variables

modelNum=21
hv='390 560 900 1240 1340'
gasName='Ar-iC4H10'     # Ar-iC4H10 or Ne or Ar-CO2
nEvents=10            # number of events to simulate
computeIBF=1
useFeSource=0
testMode=0

## end of variables


delete() {
for entry in $filesToDelete
do
if test -f $entry;
then
#echo Deleting $entry
rm $entry
fi
done
}

filesToDelete="Include/*.d Include/*.pcm Include/*.so *.d *.so *.pcm job.sub log/* error/* output/* *.txt"
delete

touch input.txt
#echo Creating new input.txt
echo '# variables' >> input.txt
echo 'modelNum ='  $modelNum >> input.txt
echo 'gasName =' $gasName'     # Ar-iC4H10 or Ne or Ar-CO2'  >> input.txt
echo 'nEvents =' $nEvents'            # number of events to simulate'  >> input.txt
echo 'computeIBF =' $computeIBF '         # if false, it will only compute the number of amplification electrons in the avalanche (in signal.C)'  >> input.txt
echo 'useFeSource =' $useFeSource'        # in signal.C: in false, will only simulate one ionisation in the drift region / if true, it will simulate a photon that converts into electrons in the drift region'  >> input.txt
echo 'testMode =' $testMode'		# in signal.C: tests the code on a reduced number of events'  >> input.txt
echo '# to draw the avalanche'  >> input.txt
echo 'plotDrift2D = 0'  >> input.txt
echo 'plotDrift3D = 1'  >> input.txt


### Create the string hv with - between numbers
for entry in $hv; do
    hvString=$hvString$entry-
done
hvString="${hvString%?}"	# delete last -
#echo $hvString

### Create filename depending on the variables
if [ $useFeSource == 1 ]
then
filename=fe-
fi
filename=${filename}signal-
if [ $computeIBF == 0 ]
then
filename=${filename}noibf-
fi
filename=${filename}$hvString

phraseIn=rootFiles-$gasName-model$modelNum-$filename-0.root
phraseOut=rootFiles/$gasName/model$modelNum/$filename'-$(ProcId).root'
dataFolder='COMSOL_data/model'$modelNum
inputMeshFile=$dataFolder/mesh.mphtxt
inputMatFile=$dataFolder/dielectrics.dat
inputFieldFile=$dataFolder/ewfield-$hvString.txt

# make sure output directory exists
createFolder ()
{
if ! test -d $folder;
then
echo Creating new folder $folder
mkdir $folder
fi
}
folder=rootFiles
createFolder
folder=rootFiles/$gasName
createFolder
folder=rootFiles/$gasName/model$modelNum
createFolder

echo
echo Computing signal root files for model $modelNum, gas = $gasName and HV = $hv
echo

if [ $testMode == 0 ]
then

touch job.sub
#echo Creating new job.sub
echo 'executable        = signal' >> job.sub
echo 'arguments = '$hv 0 >> job.sub
echo 'output            = output/ex.$(ClusterId).$(ProcId).out' >> job.sub
#echo 'input          = input.txt' $inputMeshFile $inputMatFile $inputFieldFile >> job.sub
echo 'input             = '$dataFolder/ >> job.sub
echo 'transfer_input_files = input.txt' >> job.sub
echo 'error             = error/ex.$(ClusterId).$(ProcId).err' >> job.sub
echo 'log               = log/ex.$(ClusterId).$(ProcId).log' >> job.sub
echo 'getenv            = true' >> job.sub
#echo 'environment      = LD_LIBRARY_PATH=${GARFIELD_HOME}/Install/lib:$LD_LIBRARY_PATH' >> job.sub
#echo 'transfer_output_files   =' "$phrase" >> job.sub
echo '+MaxRuntime       = 460000' >> job.sub	# that is, more than 5 days
echo 'request_cpus      = 2' >> job.sub
#echo '++JobFlavour = "tomorrow"' >> job.sub # tomorrow for one day or testmatch for 3 days, nextweek for one week, workday for 8 hours, longlunch for 2 hours
echo 'transfer_output_remaps  = "'$phraseIn'='$phraseOut'" ' >> job.sub
echo 'queue' $numberOfJobs >> job.sub

make
condor_submit job.sub

else

make
./signal $hv

fi


filesToDelete="Include/*.d Include/*.pcm Include/*.so *.d *.so *.pcm"
delete

