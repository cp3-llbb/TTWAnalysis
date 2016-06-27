import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework

from cp3_llbb.TTWAnalysis.Configuration import addTTWAnalyzer, addTTWCandidatesAnalyzer, customizeProducers

globalTag_ = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1'
processName_ = 'PAT'

framework = Framework.Framework(False, eras.Run2_25ns, globalTag=globalTag_, processName=processName_)

## ANALYZERS
addTTWAnalyzer          (framework) ## make candidates
addTTWCandidatesAnalyzer(framework) ## fill most branches

from cp3_llbb.TTWAnalysis.Configuration import lepton_WPs
framework.addAnalyzer('ttWTruth', cms.PSet(
        type = cms.string('ttwtruth_analyzer'),
        prefix = cms.string(''),
        enable = cms.bool(True),
        parameters = cms.PSet(
            TTWAnalyzer=cms.string("ttW"),
            JetWP=cms.vstring(*("ID{0}_Iso{1}".format(lID,lIso) for lID,lIso in lepton_WPs.iterkeys()))
            )
        ))

## PRODUCERS
customizeProducers(framework)

# framework.redoJEC()
# framework.smearJets()
# framework.doSystematics(['jec', 'jer'])
framework.doSystematics([])

process = framework.create()

process.source.fileNames = cms.untracked.vstring(
    '/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/04D51FB4-B2B8-E511-A399-047D7B881D6A.root'
    )

## process.source.fileNames = cms.untracked.vstring(
##     '/store/mc/RunIIFall15MiniAODv2/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/14C51DB0-D6B8-E511-8D9B-8CDCD4A9A484.root'
##     )

## process.source.fileNames = cms.untracked.vstring(
##     '/store/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0ED0E1CB-90BF-E511-B379-0025905C4432.root'
##     )

## Tricky gen event from /store/mc/RunIISpring15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/00000/0014DC94-DC5C-E511-82FB-7845C4FC39F5.root
## First one is g g -> t tbar with one W -> bbar c
## Second is b bar -> t tbar semi-leptonic
#process.source.eventsToProcess = cms.untracked.VEventRange(
#        '1:52386:13083444',
#        '1:34020:8496854'
#        )

## Other tricky gen events, with lots of ISR
## From file:/nfs/scratch/fynu/swertz/CMSSW_7_4_15/src/cp3_llbb/TTAnalysis/test/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_miniAODv2_oneFile.root
#process.source.eventsToProcess = cms.untracked.VEventRange(
#        '1:321521:80300260',
#        '1:357590:89308562',
#        '1:387992:96901374'
#        )

#process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
