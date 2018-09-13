## Installation instructions:
##SVFit analysis requires downloading the SVFit package.

``` 
git clone https://github.com/akalinow/RootAnalysis
git clone https://github.com/akalinow/ClassicSVfit -b tau_sv_features_code_optimalization TauAnalysis/ClassicSVfit
git clone https://github.com/SVfit/SVfitTF TauAnalysis/SVfitTF

cd -; cd RootAnalysis
git checkout relevant_tag
mkdir build; cd build
cmake ../
make install -j 4
```
## Run instructions:
Before running update path to the data files in the svfit_MuTau.ini file.

```
cd RootAnalysis/build
./bin/svfitAnalysis config/svfit_MuTau.ini
```
