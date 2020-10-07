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
See Pythia8Interface/README for instructions hot to generate events with Pythia8.

```
cd RootAnalysis/build
ln -s config/ml_Properties.ini ml_Properties.ini
./bin/svfitAnalysis config/svfit_MuTau_train.ini >& train.out &
./bin/svfitAnalysis config/svfit_MuTau_train_smearMET.ini >& train.out &
./bin/svfitAnalysis config/svfit_MuTau_test_ggH125.ini >& ggH125.out &
./bin/svfitAnalysis config/svfit_MuTau_test_DY.ini >& DY.out &	
```

Convert the plain ROOT TTRee into python numpy array

```
python python/root2pickle.py --input RootAnalysis_SVfitMLAnalysisMuTauTrain.root --output htt_features_train.pkl
python python/root2pickle.py --input RootAnalysis_SVfitMLAnalysisMuTau_ggH125.root --output htt_features_ggH125.pkl
python python/root2pickle.py --input RootAnalysis_SVfitMLAnalysisMuTau_DY.root --output htt_features_DY.pkl
```

Train a simple NN (requires TensorFlow)

```
python neuralnetwork_svfit.py
```

For more elaborate ML code see the [MachineLearning](https://github.com/akalinow/MachineLearning) repository