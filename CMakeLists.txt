#---Setup the project-----------------------------------------------------------
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(COMSOL)

#---Find Garfield package-------------------------------------------------------
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

#---Locate the ROOT package and define a number of useful targets and variables-
#find_package(ROOT REQUIRED COMPONENTS Tree RIO)
if(NOT TARGET ROOT::Tree)
  find_package(ROOT Tree)
endif()

if(NOT TARGET ROOT::RIO)
  find_package(ROOT RIO)
endif()


#---Build executables------------------------------------------------------------

add_executable(avalanche avalanche.C)
target_link_libraries(avalanche Garfield::Garfield)

add_executable(plotField plotField.C)
target_link_libraries(plotField Garfield::Garfield)

add_executable(spectrumFe55 spectrumFe55.C)
target_link_libraries(spectrumFe55 Garfield::Garfield ROOT::RIO)

add_executable(signal signal.C)
target_link_libraries(signal Garfield::Garfield ROOT::RIO ROOT::Tree)


