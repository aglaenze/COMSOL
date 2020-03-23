# COMSOL
Garfield project for IBF simulation in a Micromegas detector

Don't forget to change in the source code the unit length (in ComponentComsol.cc) --> it has to be in um and not m, replace 100.0 with 1.e-4

- parameters.C: contains geometry characterization of each model

- plotField : to draw electric field + voltage in the detector: broad view + zoom view on the amplification region

- avalanche: to draw one avalanche signal, started by one electron in the drift region

- spectrumFe55: simulates nEvents photons that convert into electrons in the detector, then stores the number of primary electrons produced in a TH1 in a rootfile

- signal: simulates the signal generated by avalanches started by one electron (x nEvents), and stores the number of secondary electrons + current signals in root TTrees

- feSignal: simulates the Fe source with a given rate of events, then the photons convert into primary electrons, each of which then starts an avalanche

- GetNumberOfCharges.C: takes in input the signal root files, to integrate the current and extract the number of charges that generated a signalm (should correspond to the gain)

- Convolute.C: convolutes the number of secondary electrons (obtained by counting and by extracting the number of charges) with the number of primaries. Output = histograms of amplified Fe spectra

- Analyse.C : uses the previously generated root files to analyse the data
--> fit gain and ibf
--> draw gain and ibf curve

- GetIbf.C: take the output files of signal or feSignal to read the currents on the pad electrode and the drift electrode and compute the ibf





Note:
_Utils.C and parameters.C have to be updated / at least checked when using another model than 1


In obsolete, to compute gain and IBF by counting the electrons and ions :

To draw the gain:
1) ./gain
2) root -l -q Convolute.C
3) root -l -q Analyse.C

To compute the IBF:
1) ./ibf
2) root -l -q Analyse.C

- drawGeometry : draws a geometry of solid boxes (shape of the detector) in 3D, kind of useless
