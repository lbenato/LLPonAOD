# LLPonAOD

## Git instructions
In your working area, first set up the CMSSW release:
```bash
### CMSSW_8_0_26_patch1 is the release needed for 2016 data analyses
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src
cmsenv
git cms-init
```

then, clone the git repository:
```bash
### CLONE THE REPO
cd $CMSSW_BASE/src
mkdir Analyzer
cd Analyzer
git clone https://github.com/lbenato/LLPonAOD.git
cd LLPonAOD
```

and compile the code:
```bash
### COMPILE
cd $CMSSW_BASE/src
scram b -j 32
```

## Run the code
```bash
### COMPILE
cmsRun python/ConfFile_cfg.py
```