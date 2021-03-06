# COMSOL
Garfield project for ion back flow simulation (IBF) in a Micromegas detector

Don't forget to change in the source code the unit length (in ComponentComsol.cc) --> it has to be in um and not m, replace ```unit = 100.0``` with 1.e-4



## How to use

1) Create a file input.txt, that contains the variables: gas name, model number, ... (see further an example of input.txt)  
And then
```./signal HV1 HV2 HV3 ... + saveNum```
3) ```root -l -q "AddSignalTrees.C($modelNum, \"$gasName\", $hvList)"``` to add the signal root files
4) ```root -l -q "Convolute.C($modelNum, \"$gasName\", $hvList)"``` to convolute the gain obtained in signal root files with the spectrum of Fe (obtained with ./spectrumFe55, see Description of macros)
5) ```root -l -q "Analyse.C($modelNum, \"$gasName\", $hvList)"```

Alternatively, steps 2 to 4 can be written in a Process.sh executable that would look like this:  
```
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

root -l -q "AddSignalTrees.C($modelNum, \"$gasName\", $hvList)"
root -l -q "Convolute.C($modelNum, \"$gasName\", $hvList)"
root -l -q "Analyse.C($modelNum, \"$gasName\", $hvList)"
```  

And then: ```./Process.sh HV1 HV2 HV3 ... ```  

Don't forget to change gas name and model number inside Process.sh.


## Description of macros in there

- WhatPossibilities.sh: says which electric field configuration exist for one given model, or for all models, with the format: HV1 HV2 HV3 ... So that it can be copy-pasted in ```./signal HV1 HV2 HV3 ... savenum```
To use: ```./WhatPossibilities.sh $modelNum ```(for one model only)
or just ```./WhatPossibilities.sh``` for all models

### For simulations in Garfield

- ComputeSignal.sh: automatically starts ./signal for the different electric configurations of a model, by chunks of 4

- avalanche: to draw one avalanche signal, started by one electron in the drift region  
To use: First change input.txt with right variables, then ```make; ./avalanche HV1 HV2 HV3 ...```

- plotField : to draw electric field + voltage in the detector: broad view + zoom view on the amplification region  
To use: First change input.txt with right variables, then ```make; ./plotField HV1 HV2 HV3 ...```

- spectrumFe55: simulates nEvents photons that convert into electrons in the detector, then stores the number of primary electrons produced in a TH1 in a rootfile
[update] also stores in a TH1 the proportion of photons that converted, as a function of the thickness of the detector  
To use: ```make; ./spectrumFe55```

- signal: simulates the signal generated by avalanches started by one electron, and stores the number of secondary electrons + ibf ratio + current signals (and their integrals) in TTrees . There is also a possibility to compute signal generated by an incoming photon that converts in the drift region and produces a few hundreds of electrons, by saying say useFeSource=1 in input.txt.
To use: First change input.txt with right variables, then ```make; ./signal HV1 HV2 HV3 ... saveNum```


### For the analysis of the simulations: this uses the output of signal.C

- AddSignalTrees.C: merges the root trees generated by signal.C

- Convolute.C: convolutes the number of secondary electrons (obtained by counting and by extracting the number of charges) with the number of primaries. Output = histograms of amplified Fe spectra

- Analyse.C: complete analysis of the simulation: draws the gain, IBF, detector configuration

- ResVsIbf: draws resolution against IBF.

- Spectrum.C: macro to compare 1) spectra convoluted with a Fe source offline, after signal.C simulation and 2) simulation of a Fe source


For model1 only:
- GetGain.C: takes the output files of Convolute.C to get the gain of the spectra
- GetIbf.C: take the output files of signal or feSignal to read the currents on the pad electrode and the drift electrode and compute the ibf

### For curiosity, not related to the detector model

- PhotonConversionRatio.C: to draw the proportion of photons that converted in the detector as a function of the thickness of the detector, and compare it to the theoretical calculation using the cross section (not found! --> find the cross section of photon at 5.9keV in Ar to finish this)


### In the folder Include:

- Initiate.C: functions to initate gas + electric field for Garfield++
- Utils.C contains a few useful functions
- Data.C contains data for ZZBOT (simple Micromegas) to compare with simulations
- Geometry.C: functions related to the geometry of detectors
- Functions.C: functions for the analysis (fit functions etc)
- Transparency.C: functions to draw the position of ions where they stop and compute the transparency of a given electrode



## In the folder unused (you should NOT use it):

To compute gain and IBF by counting the electrons and ions:

To draw the gain:  
1) ``` ./gain```
2) ```root -l -q Convolute.C```
3) ```root -l -q Analyse.C```

To compute the IBF:
1) ```./ibf```
2) ```root -l -q Analyse.C```

- drawGeometry : draws a geometry of solid boxes (shape of the detector) in 3D, kind of useless

- feSignal: simulates the Fe source with a given rate of events, then the photons convert into primary electrons, each of which then starts an avalanche [takes too long]

- GetNumberOfCharges.C: takes in input the signal root files. It looks at the current = f(time); each time the current goes above a threshold, it'a a new event. For each event, it integrates the current above the thrshold to extract the number of charges (should correspond to the gain)

- Analyse.C : uses the root files previously generated with gain.C and ibf.C to analyse the data
--> fit gain and ibf
--> draw gain and ibf curve
--> compare with real data

- GetIbf.C: to get IBF by looking at 1D distribution of current --> gives an underestimated IBF



# Appendix

### Example of input.txt

```
# variables
modelNum = 10
gasName = Ar-iC4H10     # Ar-iC4H10 or Ne or Ar-CO2
nEvents = 100	        # number of events to simulate
computeIBF = 0          # if false, it will only compute the number of amplification electrons in the avalanche (in signal.C)
useFeSource = 1		# in signal.C: in false, will only simulate one ionisation in the drift region / if true, it will simulate a photon that converts into electrons in the drift region
# to draw the avalanche
plotDrift2D = 0
plotDrift3D = 1
```

### Executable ComputeSignal.sh to locally parallelize jobs

```
#!/bin/bash

task()
{
echo "Vmesh = $1, Vdrift = $2 and Index = $3"
echo ./signal $1 $2 $3
./signal $1 $2 $3
}

# run jobs in parallel, by packets of 5
N=5
for ((k=0;k<=6;k++)); do
V1=$((340+$k*20))
V2=$((540+$k*20))
for ((i=1;i<=5;i++)); do
# chunks of 5
((j=j%N)); ((j++==0)) && wait
task "$V1" "$V2" "$i"&
#echo "$V $i"
done
done
```  
