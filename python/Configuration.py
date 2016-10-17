"""
Configuration for ttW trees (shared between data and MC)
"""
__author__ = "Pieter David <pieter.david@uclouvain.be>"

from collections import OrderedDict as odict
from itertools import product, tee, chain

import FWCore.ParameterSet.Config as cms

dileptonTriggers = {
      "ElEl" : ["HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*"]
    , "ElMu" : ["HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*", "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*"]
    , "MuMu" : ["HLT_Mu17_TrkIsoVVL_(Tk)?Mu8_TrkIsoVVL_DZ_v.*"]
    }
dileptonTriggers["MuEl"] = dileptonTriggers["ElMu"]

def makeCategoryParams(llWPs=[], diLeptonTriggerMatch=False, addPassAll=False):
    categs = dict() ## take dilepton working points from input
    for l1 in ("El", "Mu"):
        for l2 in ("El", "Mu"):
            flav = "".join((l1,l2))
            base = cms.PSet(
                      NElectrons = cms.uint32(sum( 1 for l in (l1,l2) if l == "El" ))
                    , NMuons     = cms.uint32(sum( 1 for l in (l1,l2) if l == "Mu" ))
                    , Category   = cms.string("is{0}".format(flav))
                    , HLT        = cms.vstring(dileptonTriggers[flav]) if diLeptonTriggerMatch else cms.vstring()
                    , Cuts       = cms.VPSet(
                                        cms.PSet(Mll   = cms.string("p4.M > 20"))
                                      , cms.PSet(ZVeto = cms.string("( p4.M < 76 ) || ( p4.M > 116 )"))
                                      )
                    , WPs=cms.vstring(llWPs)
                    )
            categs["{0}OS".format(flav)]    = base.clone(Charge=cms.int32( 0), Category=cms.string("is{0} && isOS".format(flav)))
            categs["{0}Plus".format(flav)]  = base.clone(Charge=cms.int32( 1))
            categs["{0}Minus".format(flav)] = base.clone(Charge=cms.int32(-1))
    if addPassAll:
        categs["all"] = cms.PSet(
                  NElectrons = cms.uint32(0)
                , NMuons     = cms.uint32(0)
                , Category   = cms.string("")
                , HLT        = cms.vstring()
                , Cuts       = cms.VPSet()
                , WPs        = cms.vstring()
                , Charge     = cms.int32(0)
                )

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
lepton_WPs = odict(((idKy[0], isoKy[0]), ("( {id} && {iso} )".format(id=el_ID_WPs[idKy], iso=el_Iso_WPs[isoKy]), "( {id} && {iso} )".format(id=mu_ID_WPs[idKy], iso=mu_Iso_WPs[isoKy]))) for idKy,isoKy in product(("Loose", "Medium"), ("Loose",)))

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
                  #, ("Tight" , 0.935)
                  ])

def addTTWAnalyzer(framework, name="ttW", prefix="ttW_", applyFilter=True):
    framework.addAnalyzer(name, cms.PSet(
        type = cms.string('ttw_analyzer'),
        prefix = cms.string(prefix),
        enable = cms.bool(True),
        parameters = cms.PSet(
            electronsProducer = cms.string('electrons'),
            muonsProducer = cms.string('muons'),
            jetsProducer = cms.string('jets'),
            metProducer = cms.string('met'),

            electronPtCut = cms.untracked.double(5),
            electronEtaCut = cms.untracked.double(2.5),

            muonPtCut = cms.untracked.double(5),
            muonEtaCut = cms.untracked.double(2.4),

            jetPtCut = cms.untracked.double(25),
            jetEtaCut = cms.untracked.double(2.5),
            #jetPUID = cms.untracked.double(-9999999),
            jetDRleptonCut = cms.untracked.double(0.4),

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
        categories_parameters = makeCategoryParams(llWPs=["ID{0}{1}_Iso{2}{3}".format(id1,id2,iso1,iso2) for (id1,iso1), (id2,iso2) in product(*tee(lepton_WPs.iterkeys()))], diLeptonTriggerMatch=False, addPassAll=(not applyFilter))
        ))

def addTTWCandidatesAnalyzer(framework, name="fillLists", prefix=""):
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
                            MiniIso=cms.PSet(type=cms.string("ttw_electronHybridFunctions"), parameters=cms.PSet(
                                Functions=cms.PSet(**dict((nm, cms.string('userFloat("{0}")'.format(nm))) for nm in chain(
                                    ("miniIsoR", "miniAbsIsoCharged", "miniAbsIsoPho", "miniAbsIsoNHad", "miniAbsIsoPU"),
                                    ("_".join((chnabsrel, strat))
                                        for chnabsrel in ("miniAbsIsoNeutral", "miniAbsIso", "miniRelIso")
                                        for strat in ("weights", "raw", "rhoArea", "deltaBeta"))
                                    )))
                                )),
                            SF =cms.PSet(type=cms.string("ttw_electronSF" ), parameters=cms.PSet(scale_factors=cms.untracked.PSet(
                                id_veto   = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electron_CutBasedID_VetoWP_fromTemplates_withSyst_76X.json'),
                                id_loose  = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electron_CutBasedID_LooseWP_fromTemplates_withSyst_76X.json'),
                                id_medium = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electron_CutBasedID_MediumWP_fromTemplates_withSyst_76X.json'),
                                id_tight  = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electron_CutBasedID_TightWP_fromTemplates_withSyst_76X.json'),
                                hww_wp    = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electrons_HWW_CutBasedID_TightWP_76X_forHWW_Final.json'),
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
                                WeightsFile=cms.FileInPath("cp3_llbb/TTWAnalysis/data/el_BDTG.weights.xml"),
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
                            MiniIso=cms.PSet(type=cms.string("ttw_muonHybridFunctions"), parameters=cms.PSet(
                                Functions=cms.PSet(**dict((nm, cms.string('userFloat("{0}")'.format(nm))) for nm in chain(
                                    ("miniIsoR", "miniAbsIsoCharged", "miniAbsIsoPU"),
                                    ("_".join((chnabsrel, strat))
                                        for chnabsrel in ("miniAbsIsoNeutral", "miniAbsIso", "miniRelIso")
                                        for strat in ("weights", "raw", "rhoArea", "deltaBeta"))
                                    )))
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
                                WeightsFile=cms.FileInPath("cp3_llbb/TTWAnalysis/data/mu_BDTG.weights.xml"),
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
                                **dict(("HLTMatch_{0}".format(ky), cms.vstring(*trigs)) for ky,trigs in dileptonTriggers.iteritems())
                                ))),
                            HLT2 =cms.PSet(type=cms.string("ttw_dileptonHLTMatchv2"), parameters=cms.PSet(
                                triggers=cms.untracked.FileInPath("cp3_llbb/TTWAnalysis/data/trigger.xml"),
                                Selections=cms.PSet(**dict(("HLTMatch2_{0}".format(ky), cms.vstring(*trigs)) for ky,trigs in dileptonTriggers.iteritems())),
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

def customizeProducers(framework):
    ## add single lepton triggers to the default cp3_llbb HLT matching
    framework.getProducer("hlt").parameters.triggers = cms.untracked.FileInPath("cp3_llbb/TTWAnalysis/data/trigger.xml")

    ## replace default producers
    framework.removeProducer("electrons")
    framework.addProducer("electrons", cms.PSet(type=cms.string("ttw_electronproducer"), enable=cms.bool(True),
        prefix=cms.string("electron_"),
        parameters=cms.PSet(
            input=cms.InputTag("slimmedElectrons"),
            DictTools=cms.PSet(
                MiniIso=cms.PSet(type=cms.string("ttw_electronMiniIso"), parameters=cms.PSet(
                    ea=cms.untracked.FileInPath("RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt"),
                    rho=cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                    packedCandidates=cms.InputTag("packedPFCandidates"),
                    )),
                )
            )
        ))
    framework.removeProducer("muons")
    framework.addProducer("muons",     cms.PSet(type=cms.string("ttw_muonproducer")    , enable=cms.bool(True),
        prefix=cms.string("muon_"),
        parameters=cms.PSet(
            input=cms.InputTag("slimmedMuons"),
            DictTools=cms.PSet(
                MiniIso=cms.PSet(type=cms.string("ttw_muonMiniIso"), parameters=cms.PSet(
                    ea=cms.untracked.FileInPath("cp3_llbb/TTWAnalysis/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_Spring15_25ns.txt"),
                    rho=cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                    packedCandidates=cms.InputTag("packedPFCandidates"),
                    )),
                )
            )
        ))
    framework.removeProducer("jets")
    framework.addProducer("jets",      cms.PSet(type=cms.string("ttw_jetproducer")     , enable=cms.bool(True),
        prefix=cms.string("jet_"),
        parameters=cms.PSet(
            input=cms.InputTag("slimmedJets"),
            cut=cms.untracked.string("pt > 10"),
            DictTools=cms.PSet())
        ))

    framework.removeProducer('fat_jets')
