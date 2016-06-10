from collections import OrderedDict as odict
from itertools import product, tee, chain

import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework

globalTag_ = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1'
processName_ = 'PAT'

framework = Framework.Framework(False, eras.Run2_25ns, globalTag=globalTag_, processName=processName_)

triggersPerChannel = {
      "ElEl" : ["HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*"]
    , "ElMu" : ["HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*", "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*"]
    , "MuMu" : ["HLT_Mu17_TrkIsoVVL_(Tk)?Mu8_TrkIsoVVL_DZ_v.*"]
    }
triggersPerChannel["MuEl"] = triggersPerChannel["ElMu"]

def makeCategoryParams(llWPs=[], diLeptonTriggerMatch=False):
    categs = dict() ## take dilepton working points from input
    for l1 in ("El", "Mu"):
        for l2 in ("El", "Mu"):
            flav = "".join((l1,l2))
            base = cms.PSet(
                      NElectrons = cms.uint32(sum( 1 for l in (l1,l2) if l == "El" ))
                    , NMuons     = cms.uint32(sum( 1 for l in (l1,l2) if l == "Mu" ))
                    , Category   = cms.string("is{0}".format(flav))
                    , HLT        = cms.vstring(triggersPerChannel[flav]) if diLeptonTriggerMatch else cms.vstring()
                    , Cuts       = cms.VPSet(
                                        cms.PSet(Mll   = cms.string("p4.M > 20"))
                                      , cms.PSet(ZVeto = cms.string("( p4.M < 76 ) || ( p4.M > 116 )"))
                                      )
                    , WPs=cms.vstring(llWPs)
                    )
            categs["{0}OS".format(flav)]    = base.clone(Charge=cms.int32( 0), Category=cms.string("is{0} && isOS".format(flav)))
            categs["{0}Plus".format(flav)]  = base.clone(Charge=cms.int32( 1))
            categs["{0}Minus".format(flav)] = base.clone(Charge=cms.int32(-1))
    return cms.PSet(**categs)

## Lepton identification and isolation working points
el_ID_WPs  = odict((nm, "electronID('{0}')".format(sel))
                    for nm,sel in [
                      ("Veto"  , "cutBasedElectronID-Spring15-25ns-V1-standalone-veto")
                    , ("Loose" , "cutBasedElectronID-Spring15-25ns-V1-standalone-loose")
                    , ("Medium", "cutBasedElectronID-Spring15-25ns-V1-standalone-medium")
                    , ("Tight" , "cutBasedElectronID-Spring15-25ns-V1-standalone-tight")
                    ])
mu_ID_WPs  = odict(Loose ="isLooseMuon",
                   Medium="isMediumMuon",
                   Tight ="( userInt('tightMuonID') != 0 )")
mu_Iso_var = "( {iso}.sumChargedHadronPt + max(({iso}.sumNeutralHadronEt + {iso}.sumPhotonEt) - 0.5*({iso}.sumPUPt), 0.) ) / pt".format(iso="pfIsolationR04") # relativeIsoR04_deltaBeta
mu_Iso_WPs = odict(Loose ="( ( {0} ) < .25 )".format(mu_Iso_var),
                   Tight ="( ( {0} ) < .15 )".format(mu_Iso_var)
                  )
el_Iso_WPs = odict(Loose = "1 == 1",
                   Tight = "1 == 0"
                  )
## format: { (IDnm, ISOnm) : (el-cutStr, mu-cutStr) }
lepton_WPs = odict(((idKy[0], isoKy[0]), ("( {id} && {iso} )".format(id=el_ID_WPs[idKy], iso=el_Iso_WPs[isoKy]), "( {id} && {iso} )".format(id=mu_ID_WPs[idKy], iso=mu_Iso_WPs[isoKy]))) for idKy,isoKy in product(("Loose", "Medium", "Tight"), ("Loose", "Tight")))

## Jet ID: https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
JetIDDefs = dict((k,v.format(
                  NHF ="neutralHadronEnergyFraction"
                , NEMF="neutralEmEnergyFraction"
                , CHF ="chargedHadronEnergyFraction"
                , MUF ="muonEnergyFraction"
                , CEMF="chargedEmEnergyFraction"
                , CHM ="chargedMultiplicity"
                , NumNeutralParticles="neutralMultiplicity"
                , NumConst="(chargedMultiplicity+neutralMultiplicity)"
            )) for k,v in {
              "Loose" : ("( ( abs(eta) <= 3. ) && "
                             "( ( ({NHF}<0.99) && ({NEMF}<0.99) && ({NumConst}>1) ) "
                            "&& ( ( (abs(eta)<=2.4) && ({CHF}>0) && ({CHM}>0) && ({CEMF}<0.99) ) "
                                "|| (abs(eta)>2.4) ) ) ) "
                      "|| ( ( abs(eta) >  3. ) && "
                             "( ({NEMF}<0.90) && ({NumNeutralParticles}>10) ) )")
            , "Tight" : ("( ( abs(eta) <= 3. ) && "
                             "( ( ({NHF}<0.90) && ({NEMF}<0.90) && ({NumConst}>1) ) "
                            "&& ( ( (abs(eta)<=2.4) && ({CHF}>0) && ({CHM}>0) && ({CEMF}<0.99) ) "
                                "|| (abs(eta)>2.4) ) ) ) "
                      "|| ( ( abs(eta) >  3. ) && "
                             "( ({NEMF}<0.90) && ({NumNeutralParticles}>10) ) )")
            , "TightLeptonVeto" : ("( abs(eta) <= 3. ) && "
                             "( ( ({NHF}<0.90) && ({NEMF}<0.90) && ({NumConst}>1) && ({MUF}<0.8) ) "
                            "&& ( ( (abs(eta)<=2.4) && ({CHF}>0) && ({CHM}>0) && ({CEMF}<0.90) ) "
                                "|| (abs(eta)>2.4) ) ) ")

            }.iteritems())

bTagName = 'pfCombinedInclusiveSecondaryVertexV2BJetTags'
b_tag_WPs = odict((nm, "(abs(eta)<2.4) && (bDiscriminator('{0}')>{1:.5f})".format(bTagName, cutVal))
                  for nm,cutVal in [
                    ("Loose" , 0.460)
                  , ("Medium", 0.8  )
                  , ("Tight" , 0.935)
                  ])

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

            muonPtCut = cms.untracked.double(20),
            muonEtaCut = cms.untracked.double(2.4),

            jetPtCut = cms.untracked.double(30),
            jetEtaCut = cms.untracked.double(2.5),
            #jetPUID = cms.untracked.double(-9999999),
            jetDRleptonCut = cms.untracked.double(0.3),

            bTagName = cms.untracked.string(bTagName),

            hltDRCut = cms.untracked.double(0.3), # DeltaR cut for trigger matching
            hltDPtCut = cms.untracked.double(0.5), #Delta(Pt)/Pt cut for trigger matching

            ### new-style parameters
            ElectronWP = cms.PSet(),
            MuonWP     = cms.PSet(),
            LeptonWP   = cms.PSet(**dict(("ID{0}_Iso{1}".format(iID, iIso),
                                        cms.PSet(Electron=cms.string(elCut), Muon=cms.string(muCut)))
                                    for (iID, iIso), (elCut, muCut) in lepton_WPs.iteritems())),
            DiLeptonWP = cms.PSet(**dict(("ID{0}{1}_Iso{2}{3}".format(id1,id2,iso1,iso2),
                                        cms.PSet(Leading=cms.PSet(Electron=cms.string(el1), Muon=cms.string(mu1)),
                                                 SubLeading=cms.PSet(Electron=cms.string(el2), Muon=cms.string(mu2))))
                                    for ((id1,iso1), (el1,mu1)), ((id2,iso2), (el2,mu2)) in product(*tee(lepton_WPs.iteritems())))),
            JetID      = cms.string(JetIDDefs["Loose"]), # not tightLeptonVeto since DeltaR(l,j) cut should be enough
            BtagWP     = cms.PSet(**dict(
                ("ID{0}_Iso{1}_B{2}".format(lID,lIso, bWP[0]),
                    cms.PSet(LeptonWP=cms.string("ID{0}_Iso{1}".format(lID,lIso)), Selection=cms.string(wpSel)))
                for (lID,lIso),(bWP,wpSel) in product(lepton_WPs.iterkeys(), b_tag_WPs.iteritems()))),
            DiBtagWP   = cms.PSet(**dict(("ID{0}_Iso{1}_B{2}{3}".format(lID,lIso, nm1[0],nm2[0]),
                                        cms.PSet(LeptonWP  =cms.string("ID{0}_Iso{1}".format(lID,lIso)),
                                                 Leading   =cms.string(sel1),
                                                 SubLeading=cms.string(sel2)))
                                    for (lID, lIso), (nm1,sel1), (nm2,sel2) in product(lepton_WPs.iterkeys(), *tee(b_tag_WPs.iteritems())))),
            ),
        categories_parameters = makeCategoryParams(llWPs=["ID{0}{1}_Iso{2}{3}".format(id1,id2,iso1,iso2) for (id1,iso1), (id2,iso2) in product(*tee(lepton_WPs.iterkeys()))], diLeptonTriggerMatch=False)
        )
    )

framework.addAnalyzer('fillLists', cms.PSet(
    type = cms.string('ttw_delegatinganalyzer'),
    prefix = cms.string(''),
    enable = cms.bool(True),
    parameters = cms.PSet(
        Helpers    = cms.PSet(
            ### BASIC OBJECTS
            Electrons=cms.PSet(type=cms.string("ttw_electronsanalyzerhelper"), prefix=cms.string("electron_"),
                parameters=cms.PSet(TTWAnalyzer=cms.string("ttW"),
                    DictTools=cms.PSet(
                        Kin=cms.PSet(type=cms.string("ttw_electronKin"), parameters=cms.PSet()),
                        Gen=cms.PSet(type=cms.string("ttw_electronGenMatch"), parameters=cms.PSet()),
                        ID =cms.PSet(type=cms.string("ttw_electronHybridCuts"), parameters=cms.PSet(
                            Cuts=cms.PSet(**dict(("ID{0}".format(wpNm), cms.string(wpSel)) for wpNm, wpSel in el_ID_WPs.iteritems()))
                            )),
                        IDVars=cms.PSet(type=cms.string("ttw_electronIDVars"),parameters=cms.PSet()),
                        Iso=cms.PSet(type=cms.string("ttw_electronIso"), parameters=cms.PSet(
                            ea_R03 = cms.untracked.FileInPath("RecoEgamma/ElectronIdentification/data/PHYS14/effAreaElectrons_cone03_pfNeuHadronsAndPhotons.txt"),
                            ea_R04 = cms.untracked.FileInPath("cp3_llbb/Framework/data/effAreaElectrons_cone04_pfNeuHadronsAndPhotons.txt")
                            )),
                        MiniIso=cms.PSet(type=cms.string("ttw_electronMiniIso"), parameters=cms.PSet(
                            ea=cms.untracked.FileInPath("RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt"),
                            rho=cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                            packedCandidates=cms.InputTag("packedPFCandidates"),
                            )),
                        SF =cms.PSet(type=cms.string("ttw_electronSF" ), parameters=cms.PSet(scale_factors=cms.untracked.PSet(
                            id_veto   = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electron_CutBasedID_VetoWP_fromTemplates_withSyst_76X.json'),
                            id_loose  = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electron_CutBasedID_LooseWP_fromTemplates_withSyst_76X.json'),
                            id_medium = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electron_CutBasedID_MediumWP_fromTemplates_withSyst_76X.json'),
                            id_tight  = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electron_CutBasedID_TightWP_fromTemplates_withSyst_76X.json'),
                            hww_wp    = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electrons_HWW_CutBasedID_TightWP_fromTemplates_withSyst_v1_reco_id_iso.json'),
                            ))),
                        MVAttH=cms.PSet(type=cms.string("ttw_electronMVAttH"), parameters=cms.PSet(
                            ## for matching jet
                            JetLeptonDR=cms.double(.4),
                            Jets=cms.InputTag("slimmedJets"),
                            ## for mini-isolation
                            packedCandidates=cms.InputTag("packedPFCandidates"),
                            ea=cms.FileInPath("RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt"),
                            rho=cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                            ## BDT weights
                            WeightsFile=cms.FileInPath("cp3_llbb/TTWAnalysis/data/forMoriond16_el_sigTTZ_bkgTT_BDTG.weights.xml"),
                            )),
                        )
                    )),
            Muons    =cms.PSet(type=cms.string("ttw_muonsanalyzerhelper"), prefix=cms.string("muon_"),
                parameters=cms.PSet(TTWAnalyzer=cms.string("ttW"),
                    DictTools=cms.PSet(
                        Kin=cms.PSet(type=cms.string("ttw_muonKin"), parameters=cms.PSet()),
                        Gen=cms.PSet(type=cms.string("ttw_muonGenMatch"), parameters=cms.PSet()),
                        ID =cms.PSet(type=cms.string("ttw_muonHybridCuts"), parameters=cms.PSet(
                            Cuts=cms.PSet(**dict(("ID{0}".format(wpNm), cms.string(wpSel)) for wpNm, wpSel in mu_ID_WPs.iteritems()))
                            )),
                        IDVars=cms.PSet(type=cms.string("ttw_muonIDVars"),parameters=cms.PSet()),
                        Iso=cms.PSet(type=cms.string("ttw_muonIso"), parameters=cms.PSet(
                            ea_R03=cms.untracked.FileInPath("cp3_llbb/Framework/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons.txt"),
                            ea_R04=cms.untracked.FileInPath("cp3_llbb/Framework/data/effAreaMuons_cone04_pfNeuHadronsAndPhotons.txt"),
                            )),
                        MiniIso=cms.PSet(type=cms.string("ttw_muonMiniIso"), parameters=cms.PSet(
                            ea=cms.untracked.FileInPath("cp3_llbb/TTWAnalysis/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_Spring15_25ns.txt"),
                            rho=cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                            packedCandidates=cms.InputTag("packedPFCandidates"),
                            )),
                        SF =cms.PSet(type=cms.string("ttw_muonSF" ), parameters=cms.PSet(scale_factors=cms.untracked.PSet(
                            id_soft   = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_SoftID_genTracks_id.json'),
                            id_loose  = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_LooseID_genTracks_id.json'),
                            id_medium = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_MediumID_genTracks_id.json'),
                            id_tight  = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_TightIDandIPCut_genTracks_id.json'),
                            iso_loose_id_loose  = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_LooseRelIso_LooseID_iso.json'),
                            iso_loose_id_medium = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_LooseRelIso_MediumID_iso.json'),
                            iso_loose_id_tight  = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_LooseRelIso_TightID_iso.json'),
                            iso_tight_id_medium = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_TightRelIso_MediumID_iso.json'),
                            iso_tight_id_tight  = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_TightRelIso_TightID_iso.json'),
                            id_hww              = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_MediumID_Data_MC_25ns_PTvsETA_HWW_76.json'),
                            iso_tight_id_hww    = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_ISOTight_Data_MC_25ns_PTvsETA_HWW.json'),
                            ))),
                        MVAttH=cms.PSet(type=cms.string("ttw_muonMVAttH"), parameters=cms.PSet(
                            ## for matching jet
                            JetLeptonDR=cms.double(.4),
                            Jets=cms.InputTag("slimmedJets"),
                            ## for mini-isolation
                            packedCandidates=cms.InputTag("packedPFCandidates"),
                            ea=cms.FileInPath("cp3_llbb/TTWAnalysis/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_Spring15_25ns.txt"),
                            rho=cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                            ## BDT weights
                            WeightsFile=cms.FileInPath("cp3_llbb/TTWAnalysis/data/forMoriond16_mu_sigTTZ_bkgTT_BDTG.weights.xml"),
                            )),
                        )
                    )),
            Jets     =cms.PSet(type=cms.string("ttw_jetsanalyzerhelper"), prefix=cms.string("jet_"),
                parameters=cms.PSet(TTWAnalyzer=cms.string("ttW"),
                    DictTools=cms.PSet(
                        Kin=cms.PSet(type=cms.string("ttw_jetKin"), parameters=cms.PSet()),
                        Gen=cms.PSet(type=cms.string("ttw_jetGenMatch"), parameters=cms.PSet()),
                        ID =cms.PSet(type=cms.string("ttw_jetHybridCuts"), parameters=cms.PSet(Cuts=cms.PSet(
                            **dict(("ID{0}".format(wpNm), cms.string(wpSel)) for wpNm, wpSel in JetIDDefs.iteritems()))
                            )),
                        IDVars=cms.PSet(type=cms.string("ttw_jetIDVars"), parameters=cms.PSet(
                            bTaggers=cms.PSet(**dict((tagNm, cms.string(tagNm)) for tagNm in ("pfCombinedInclusiveSecondaryVertexV2BJetTags", "pfJetProbabilityBJetTags", "pfCombinedMVABJetTags")))
                            )),
                        BSF=cms.PSet(type=cms.string("ttw_jetSF"), parameters=cms.PSet(
                            scale_factors=cms.untracked.PSet(
                                **dict(("csvv2_{wp}".format(wp=wp),
                                    cms.untracked.PSet(
                                        algorithm    =cms.untracked.string("csvv2"),
                                        working_point=cms.untracked.string(wp),
                                        files=cms.untracked.VPSet(
                                            *(cms.untracked.PSet(
                                                flavor=cms.untracked.string(flav),
                                                file  =cms.untracked.FileInPath(
                                                        'cp3_llbb/Framework/data/ScaleFactors/BTagging_{wp}_{flav}_{src}_CSVv2.json'.format(
                                                            wp=wp, flav=flav, src=("incl" if flav=="lightjets" else "mujets") ))
                                                ) for flav in ("bjets", "cjets", "lightjets") )
                                            )
                                        )
                                    )
                                    for wp in ("loose", "medium", "tight") )
                                )
                            )),
                        )
                    )),

            ### SELECTED/COMBINED CANDIDATES
            Leptons  =cms.PSet(type=cms.string("ttw_leptonsanalyzerhelper"), prefix=cms.string("ttW_lepton_"),
                parameters=cms.PSet(TTWAnalyzer=cms.string("ttW"),
                    DictTools=cms.PSet(
                        Basic=cms.PSet(type=cms.string("ttw_leptonCandidate"), parameters=cms.PSet()),
                        HLT  =cms.PSet(type=cms.string("ttw_leptonHLTMatch"), parameters=cms.PSet(Selections=cms.PSet(
                            HLTMatch_SingleMu=cms.vstring(["HLT_IsoMu20_v.*", "HLT_IsoTkMu20_v.*"]),
                            HLTMatch_SingleEl=cms.vstring(["HLT_Ele23_WPLoose_Gsf_v*"]),
                            HLTMatch_SingleElMC=cms.vstring(["HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v.*"]),
                            ))),
                        HLT2 =cms.PSet(type=cms.string("ttw_leptonHLTMatchv2"), parameters=cms.PSet(
                            triggers=cms.untracked.FileInPath("cp3_llbb/TTWAnalysis/data/trigger.xml"),
                            Selections=cms.PSet(
                                HLTMatch2_SingleMu=cms.vstring(["HLT_IsoMu20_v.*", "HLT_IsoTkMu20_v.*"]),
                                HLTMatch2_SingleEl=cms.vstring(["HLT_Ele23_WPLoose_Gsf_v*"]),
                                HLTMatch2_SingleElMC=cms.vstring(["HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v.*"]),
                                ),
                            hltDRCut = cms.untracked.double(0.3), # DeltaR cut for trigger matching
                            hltDPtCut = cms.untracked.double(0.5), #Delta(Pt)/Pt cut for trigger matching
                            ))
                        )
                    )),
            DiLeptons=cms.PSet(type=cms.string("ttw_dileptonsanalyzerhelper"), prefix=cms.string("ttW_dilepton_"),
                parameters=cms.PSet(TTWAnalyzer=cms.string("ttW"),
                    DictTools=cms.PSet(
                        Basic=cms.PSet(type=cms.string("ttw_dileptonCandidate"), parameters=cms.PSet()),
                        HLT  =cms.PSet(type=cms.string("ttw_dileptonHLTMatch"), parameters=cms.PSet(Selections=cms.PSet(
                            **dict(("HLTMatch_{0}".format(ky), cms.vstring(*trigs)) for ky,trigs in triggersPerChannel.iteritems())
                            ))),
                        HLT2 =cms.PSet(type=cms.string("ttw_dileptonHLTMatchv2"), parameters=cms.PSet(
                            triggers=cms.untracked.FileInPath("cp3_llbb/TTWAnalysis/data/trigger.xml"),
                            Selections=cms.PSet(**dict(("HLTMatch2_{0}".format(ky), cms.vstring(*trigs)) for ky,trigs in triggersPerChannel.iteritems())),
                            hltDRCut = cms.untracked.double(0.3), # DeltaR cut for trigger matching
                            hltDPtCut = cms.untracked.double(0.5), #Delta(Pt)/Pt cut for trigger matching
                            ))
                        )
                    )),
            DiJets   =cms.PSet(type=cms.string("ttw_dijetsanalyzerhelper"), prefix=cms.string("ttW_dijet_"),
                parameters=cms.PSet(TTWAnalyzer=cms.string("ttW"),
                    DictTools=cms.PSet(
                        Basic=cms.PSet(type=cms.string("ttw_dijetCandidate"), parameters=cms.PSet()),
                        )
                    )),
            DiLeptonDiJets=cms.PSet(type=cms.string("ttw_dileptondijetsanalyzerhelper"), prefix=cms.string("ttW_dileptondijet_"),
                parameters=cms.PSet(TTWAnalyzer=cms.string("ttW"),
                    DictTools=cms.PSet(
                        Basic=cms.PSet(type=cms.string("ttw_dileptondijetCandidate"), parameters=cms.PSet()),
                        )
                    )),
            DiLeptonDiJetMets=cms.PSet(type=cms.string("ttw_dileptondijetmetsanalyzerhelper"), prefix=cms.string("ttW_dileptondijetmet_"),
                parameters=cms.PSet(TTWAnalyzer=cms.string("ttW"),
                    DictTools=cms.PSet(
                        Basic=cms.PSet(type=cms.string("ttw_dileptondijetmetCandidate"), parameters=cms.PSet()),
                        )
                    )),
            )
        )
    ))

framework.addAnalyzer('ttWTruth', cms.PSet(
        type = cms.string('ttwtruth_analyzer'),
        prefix = cms.string(''),
        enable = cms.bool(True),
        parameters = cms.PSet(
            TTWAnalyzer=cms.string("ttW"),
            JetWP=cms.vstring(*("ID{0}_Iso{1}".format(lID,lIso) for lID,lIso in lepton_WPs.iterkeys()))
            )
        ))


framework.getProducer("hlt").parameters.triggers = cms.untracked.FileInPath("cp3_llbb/TTWAnalysis/data/trigger.xml")

## replace default producers
framework.removeProducer("electrons")
framework.addProducer("electrons", cms.PSet(type=cms.string("ttw_electronproducer"), enable=cms.bool(True),
    prefix=cms.string("electron_"),
    parameters=cms.PSet(
        input=cms.InputTag("slimmedElectrons"))
    ))
framework.removeProducer("muons")
framework.addProducer("muons",     cms.PSet(type=cms.string("ttw_muonproducer")    , enable=cms.bool(True),
    prefix=cms.string("muon_"),
    parameters=cms.PSet(
        input=cms.InputTag("slimmedMuons"))
    ))
framework.removeProducer("jets")
framework.addProducer("jets",      cms.PSet(type=cms.string("ttw_jetproducer")     , enable=cms.bool(True),
    prefix=cms.string("jet_"),
    parameters=cms.PSet(
        input=cms.InputTag("slimmedJets"),
        cut=cms.untracked.string("pt > 10"))
    ))

framework.removeProducer('fat_jets')

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
