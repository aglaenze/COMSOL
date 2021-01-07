#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(COMSOL)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Locate the ROOT package and define a number of useful targets and variables-
find_package(ROOT REQUIRED COMPONENTS Tree)

#---Build executables------------------------------------------------------------
add_executable(avalanche avalanche.C)
target_link_libraries(avalanche Garfield)

add_executable(plotField plotField.C)
target_link_libraries(plotField Garfield)

add_executable(spectrumFe55 spectrumFe55.C)
target_link_libraries(spectrumFe55 Garfield)

add_executable(signal signal.C)
target_link_libraries(signal Garfield ROOT::RIO ROOT::Tree)