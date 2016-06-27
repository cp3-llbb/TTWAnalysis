import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework

from cp3_llbb.TTWAnalysis.Configuration import addTTWAnalyzer, addTTWCandidatesAnalyzer, customizeProducers

globalTag_ = '76X_dataRun2_16Dec2015_v0'
processName_ = 'RECO'

framework = Framework.Framework(True, eras.Run2_25ns, globalTag=globalTag_, processName=processName_)

## ANALYZERS
addTTWAnalyzer          (framework) ## make candidates
addTTWCandidatesAnalyzer(framework) ## fill most branches

## PRODUCERS
customizeProducers(framework)

# framework.redoJEC()
# framework.smearJets()
# framework.doSystematics(['jec', 'jer'])
framework.doSystematics([])

process = framework.create()

process.source.fileNames = cms.untracked.vstring(
    '/store/data/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/00039A2E-D7A7-E511-98EE-3417EBE64696.root'
    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
