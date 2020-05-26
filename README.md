# vlq-qWqW-NanoAOD
NanoAOD framework for selecting qWqW events and running HitFit for VLQ searches

# CMSSW setup

```
cmsrel CMSSW_10_2_16
cd CMSSW_10_2_16/src
cmsenv
git cms-init
git cms-merge-topic jmhogan:TopHitFit-TLorentzVectors
scramv1 b 
```

# NanoAOD analysis (interactive)

```
git clone https://github.com/jmhogan/vlq-qWqW-NanoAOD.git 
cd vlq-qWqW-NanoAOD
git clone -b 80X https://github.com/cms-jet/PuppiSoftdropMassCorrections.git
```
- Find a NanoAOD file name to use as a test file
- Study qWqW_NanoAnalysis.C!
```
root -l -b -q runAnalysis.C\(\"root://cmsxrootd.fnal.gov//store/mc/...path../inputfile.root\",\"testoutputfile.root\",0,5000\) # skip 0 events, run 5000 events
```

# Condor job submission 
- Edit runNanoAODJobs2018.py (or make a similar file) to have good input/output paths and NanoAOD samples. Currently the data is set to read from NanoAOD-tools output stored on the LPC EOS cluster. 
- Check runCondorAnalysis.sh
```
voms-proxy-init -voms cms -valid 168:00
python -u runNanoAODJobs2018.py >& submission.log &
```

# Missing analysis elements
Probably best added by running NanoAOD-tools crab jobs to create custom NanoAOD (add settings/configs to TopHitFit-TLorentzVectors branch above?). 
- Data lumi masking
- JEC/JER uncertainties
- b-tagging scale factor application and uncertainty
- Lepton scale factors and uncertainty

To be added here in the macro:
- Probably trigger scale factors
