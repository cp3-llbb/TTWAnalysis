import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework

globalTag_ = '76X_dataRun2_16Dec2015_v0'
processName_ = 'RECO'

framework = Framework.Framework(True, eras.Run2_25ns, globalTag=globalTag_, processName=processName_)

triggersPerChannel = {
      "ElEl" : ["HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*"]
    , "ElMu" : ["HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*", "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*"]
    , "MuMu" : ["HLT_Mu17_TrkIsoVVL_(Tk)?Mu8_TrkIsoVVL_DZ_v.*"]
    }
triggersPerChannel["MuEl"] = triggersPerChannel["ElMu"]

def makeCategoryParams():
    categs = dict()
    for l1 in ("El", "Mu"):
        for l2 in ("El", "Mu"):
            flav = "".join((l1,l2))
            base = cms.PSet(
                      NElectrons = cms.uint32(sum( 1 for l in (l1,l2) if l == "El" ))
                    , NMuons     = cms.uint32(sum( 1 for l in (l1,l2) if l == "Mu" ))
                    , Category   = cms.string("is{0}".format(flav))
                    , HLT        = cms.vstring(triggersPerChannel[flav])
                    , Cuts       = cms.VPSet(
                                        cms.PSet(Mll   = cms.string("p4.M > 20"))
                                      , cms.PSet(ZVeto = cms.string("( p4.M < 76 ) || ( p4.M > 116 )"))
                                      )
                    )
            categs["{0}OS".format(flav)]    = base.clone(Charge=cms.int32( 0), Category=cms.string("is{0} && isOS".format(flav)))
            categs["{0}Plus".format(flav)]  = base.clone(Charge=cms.int32( 1))
            categs["{0}Minus".format(flav)] = base.clone(Charge=cms.int32(-1))
    return cms.PSet(**categs)

framework.addAnalyzer('ttW', cms.PSet(
        type = cms.string('ttw_analyzer'),
        prefix = cms.string('ttW_'),
        enable = cms.bool(True),
        parameters = cms.PSet(
            electronsProducer = cms.string('electrons'),
            muonsProducer = cms.string('muons'),
            jetsProducer = cms.string('jets'),
            metProducer = cms.string('met'),

            electronPtCut = cms.untracked.double(20),
            electronEtaCut = cms.untracked.double(2.5),
            electronVetoIDName = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-veto'),
            electronLooseIDName = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-loose'),
            electronMediumIDName = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-medium'),
            electronTightIDName = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-tight'),

            muonPtCut = cms.untracked.double(20),
            muonEtaCut = cms.untracked.double(2.4),
            muonLooseIsoCut = cms.untracked.double(.25), # Loose cut recommended for dilepton analysis
            muonTightIsoCut = cms.untracked.double(.15),

            jetPtCut = cms.untracked.double(30),
            jetEtaCut = cms.untracked.double(2.5),
            #jetPUID = cms.untracked.double(-9999999),
            jetDRleptonCut = cms.untracked.double(0.3),
            jetID = cms.untracked.string('loose'), # not tightLeptonVeto since DeltaR(l,j) cut should be enough
            jetCSVv2Name = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
            jetCSVv2L = cms.untracked.double(0.460),
            jetCSVv2M = cms.untracked.double(0.8),
            jetCSVv2T = cms.untracked.double(0.935),

            hltDRCut = cms.untracked.double(0.3), # DeltaR cut for trigger matching
            hltDPtCut = cms.untracked.double(0.5), #Delta(Pt)/Pt cut for trigger matching
            ),
        categories_parameters = makeCategoryParams()
        )
    )

framework.removeProducer('fat_jets')

# framework.redoJEC()
# framework.smearJets()
# framework.doSystematics(['jec', 'jer'])
framework.doSystematics([])

process = framework.create()

process.source.fileNames = cms.untracked.vstring(
    '/store/data/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/00039A2E-D7A7-E511-98EE-3417EBE64696.root'
    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
