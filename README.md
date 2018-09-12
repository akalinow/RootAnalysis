# This folder contains modifications to the analysis
by A. Kalinowski <Faculty of Physics, Univeristy of Warsaw>

## Installation instructions:

* Copy analysis from A. Kalinowski repository:
``` git clone https://github.com/akalinow/RootAnalysis -b Run2017 ```

* Follow the provided instructions to install the analysis. Make sure you can compile and run it.
* Copy folder RootAnalysis from this repository and swap files.
* Re-run cmake,  re-compile and re-run HTTAnalysis.
* You should obtain a file as output with TTree in it. The file should be in folder specified by htt_MuTau.ini

## Specifying parameters to extract

To specify which particle properties and for which particles to extract, see 
file HTTAnalysis/ml_Properties.ini

## IMPORTANT REMARK!

The main part goal of ML analysis is to produce a TTree output, which will be
later imported into python. However, the analysis framework DOES NOT produce
TTree output when running in multi-thread mode. That is why, ML analysis will
not be executed in multi-thread mode.
