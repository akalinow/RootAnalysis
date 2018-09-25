## Installation instructions:
##The main H->tau tau analysis package.

``` 
git clone https://github.com/akalinow/RootAnalysis
cd -; cd RootAnalysis
git checkout relevant_tag
mkdir build; cd build
cmake ../
make install -j 4
```
## Run instructions:
Before running update path to the data files in the config/htt_MuTau.ini file
(or other, relevant for desired decay channel)

```
cd RootAnalysis/build
./bin/httAnalysis config/htt_MuTau.ini
```

The resulting plots are stored in the fig_png directory.

