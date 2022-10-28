## Installation instructions:

``` 
git clone https://github.com/akalinow/RootAnalysis
cd RootAnalysis
git checkout relevant_tag
mkdir build; cd build
cmake ../
make install -j 4
```


## Run instructions:

Before running update path to the data files in the config/omtf_emulator.ini file

```
cd RootAnalysis/build
./bin/omtfAnalysis config/omtf_emulator.ini
```

The resulting plots are stored in the fig_png directory.
