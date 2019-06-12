## Installation instructions:
##The main tau->3mu LFV analysis package.

``` 
git clone https://github.com/akalinow/RootAnalysis -b tauLFV
cd RootAnalysis
mkdir build; cd build
cmake ../
make install -j 4
```
## Run instructions:
Before running update path to the data files in the TauLFV/Analysis/config/tau_lfv.ini file

```
cd RootAnalysis/build
/bin/lfvAnalysis config/tau_lfv.ini
```

The resulting plots are stored in the fig_png directory.

