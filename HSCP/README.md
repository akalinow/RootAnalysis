## Installation instructions:
##The main HSCP analysis package.

``` 
git clone https://github.com/akalinow/RootAnalysis
cd -; cd RootAnalysis
git checkout relevant_tag
mkdir build; cd build
cmake ../
make install -j 4
```
## Run instructions:
Before running update path to the data files in the config/hscp.ini file

```
cd RootAnalysis/build
./bin/hscpAnalysis config/hscp.ini
```

The resulting plots are stored in the fig_png directory.

