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
Before running update path to the data files in the config/svfit_MuTau.ini file.

```
cd RootAnalysis/build
./bin/svfitAnalysis config/svfit_MuTau.ini
```

Convert the plain ROOT TTRee into python numpy array

```
python python/root2pickle.py RootAnalysis_SVfitAnalysisMuTau.root
```

Train the NN (requires TensorFlow)

```
python python/neuralnetwork_svfit.py
```