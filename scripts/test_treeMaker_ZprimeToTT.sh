#!/bin/bash

nohup cmsRun MyTools/EDAnalyzers/python/ConfFile_cfg.py \
outFileBaseName=ntupleTree_ZEE \
sourceFile=sourceFiles/ZprimeToTT_M2000_W20_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM/ZprimeToTT_M2000_W20_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM.txt \
maxEvents=100000 \
> logs/ConfFile_cfg_ZEE.log &
