## Installation instructions:
The Pythia8Interface package requires downloading and compiling the Pythia8 package.
If Pythia is installed check if PYTHIA8_DIR exists.

``` 
git clone https://github.com/akalinow/RootAnalysis

cd -; cd RootAnalysis
git checkout relevant_tag
mkdir build; cd build
cmake ../
make install -j 4
```
## Run instructions:
Before running update path to the data files in the pythia8.ini file.
eventsToAnalyze should be set to 1.
eventsToGenerate controls the number of events per tau tau mass point.

```
cd RootAnalysis/build
./bin/generateEvents config/pythia8.ini
```
