## Installation instructions:
##The main tau->3mu LFV analysis package.

``` 
git clone https://github.com/akalinow/RootAnalysis -b TauLFV
cd RootAnalysis
mkdir build; cd build
cmake ../
make install -j 4
```
## Run instructions:
Before running update path to the data files in the ```TauLFV/Analysis/config/tau_lfv.ini``` file.
A test file is can be downloaded with following command:

```
cd RootAnalysis/build
wget http://akalinow.web.cern.ch/akalinow/TauLFV/DsToTau_To3Mu_MuFilter.root -P data
```

```
cd RootAnalysis/build
./bin/tauLFVAnalysis config/tau_lfv.ini
```

The resulting plots are stored in the fig_png directory.

```
display fig_png/h1DMass3Mu_DsToTau.png
```
