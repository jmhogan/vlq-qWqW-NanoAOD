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
cd vlq-qWqW-NanoAOD.git
```
- Find a NanoAOD file name to use as a test file
- ROOT macro structure, add instructions for interactive test

# Condor job submission 
- add instructions for LPC condor jobs

# Missing analysis elements
Probably best added by running NanoAOD-tools crab jobs to create custom NanoAOD (add settings/configs to TopHitFit-TLorentzVectors branch above?). 
- Data lumi masking
- JEC/JER uncertainties
- b-tagging scale factor application and uncertainty
- Lepton scale factors and uncertainty

To be added here in the macro:
- Probably trigger scale factors
