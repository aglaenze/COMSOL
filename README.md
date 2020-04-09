# COMSOL
Garfield project for IBF simulation in a Micromegas detector

For more details look at https://docs.google.com/presentation/d/1y05BzpoZSenbcIpKmfioo0CFtZaqDmReeGPT6FDoHwY/edit?usp=sharing


Don't forget to change in the source code the unit length (in ComponentComsol.cc) --> it has to be in um and not m, replace 100.0 with 1.e-4

- ComputeSignal.sh: automatically starts ./signal for the different electric configurations of a model, by chunks of 4

- parameters.C: contains geometry characterization of each model

- plotField : to draw electric field + voltage in the detector: broad view + zoom view on the amplification region

- avalanche: to draw one avalanche signal, started by one electron in the drift region

- spectrumFe55: simulates nEvents photons that convert into electrons in the detector, then stores the number of primary electrons produced in a TH1 in a rootfile
[update] stores in a TH1 the proportion of photons that converted, as a function of the thickness of the detector

- signal: simulates the signal generated by avalanches started by one electron, and stores the number of secondary electrons + ibf ratio + current signals (and their integrals) in TTrees

- Convolute.C: convolutes the number of secondary electrons (obtained by counting and by extracting the number of charges) with the number of primaries. Output = histograms of amplified Fe spectra

- GetGain.C: takes the output files of Convolute.C to get the gain of the spectra

- GetIbf.C: take the output files of signal or feSignal to read the currents on the pad electrode and the drift electrode and compute the ibf

- PhotonConversionRatio.C: to draw the proportion of photons that converted in the detector as a function of the thickness of the detector, and compare it to the theoretical calculation using the cross section (not found! --> find the cross section of photon at 5.9keV in Ar to finish this)





Note:
_Utils.C and parameters.C have to be updated / at least checked when using another model than 1


In the folder unused, to compute gain and IBF by counting the electrons and ions :

To draw the gain:
1) ./gain
2) root -l -q Convolute.C
3) root -l -q Analyse.C

To compute the IBF:
1) ./ibf
2) root -l -q Analyse.C

- drawGeometry : draws a geometry of solid boxes (shape of the detector) in 3D, kind of useless

- feSignal: simulates the Fe source with a given rate of events, then the photons convert into primary electrons, each of which then starts an avalanche [takes too long]

- GetNumberOfCharges.C: takes in input the signal root files. It looks at the current = f(time); each time the current goes above a threshold, it'a a new event. For each event, it integrates the current above the thrshold to extract the number of charges (should correspond to the gain)

- Analyse.C : uses the previously generated root files to analyse the data
--> fit gain and ibf
--> draw gain and ibf curve
--> compare with real data
