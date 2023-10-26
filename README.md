## Installation instructions:

``` 
git clone https://github.com/akalinow/RootAnalysis
cd RootAnalysis
git checkout relevant_tag
mkdir build; cd build
cmake ../
make install -j 4
```
If running on lxplus or a container a boost libraries might not be available by default.
In that case before running cmake please setup the LCG tool box environment:

```shell
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos8-gcc11-opt/setup.sh
```

## Run instructions:

Before running update path to the data files in the config/omtf_emulator.ini file

```
cd RootAnalysis/build
./bin/omtfAnalysis config/omtf_emulator.ini
```

The resulting plots are stored in the fig_png directory.
