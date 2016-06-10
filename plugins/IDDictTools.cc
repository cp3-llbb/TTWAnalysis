#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "cp3_llbb/Framework/interface/Types.h"
#include "cp3_llbb/TTWAnalysis/interface/DictTool.h"

#include "Helpers.h"

namespace TTWAnalysis {
/**
 * Additional electron identification variables
 */
class DictElectronIDVars : public DictTool<pat::Electron>, protected DictPVHelper {
public:
  DictElectronIDVars(const edm::ParameterSet& config)
    : DictTool<pat::Electron>(config)
    , DictPVHelper(config)
  {}
  virtual ~DictElectronIDVars() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    DictPVHelper                ::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
  }

  virtual Dict evaluate(const pat::Electron& el,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    bool elValid{el.originalObjectRef().isNonnull()};
    double eta = elValid ? el.superCluster()->eta() : 0.;
    const reco::Vertex* pv = event ? getPV(event) : nullptr;

    Dict ret;
    ret.add("isEB", el.isEB());
    ret.add("isEE", el.isEE());
    // Same values used for cut-based electron ID. See:
    //     https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_15/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleDzCut.cc#L64
    //     https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_15/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleDxyCut.cc#L64
    ret.add("dxy", pv && elValid ? el.gsfTrack()->dxy(pv->position()) : 0. );
    ret.add("dz" , pv && elValid ? el.gsfTrack()->dz (pv->position()) : 0. );
    ret.add("dca", el.dB(pat::Electron::PV3D)/el.edB(pat::Electron::PV3D));

    // electron trigger emulation
    //     https://github.com/peruzzim/cmgtools-lite/blob/76X_for2016basis/TTHAnalysis/python/tools/functionsTTH.py#L10-L20
    bool passTrigEmu{false};
    if ( elValid ) {
      double eInvMinusPInv = el.ecalEnergy() > 0. ? (1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy()) : 9.e9;
      if ( std::abs(eta) < 1.479 ) {
        passTrigEmu = (
              ( el.hadronicOverEm() < 0.10 )
           && ( std::abs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.01 )
           && ( std::abs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.04 )
           && ( -0.05 < eInvMinusPInv ) && ( eInvMinusPInv < 0.01 )
           && ( el.full5x5_sigmaIetaIeta() < 0.011 )
           );
      } else {
        passTrigEmu = (
              ( el.hadronicOverEm() < 0.07 )
           && ( std::abs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.008 )
           && ( std::abs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.07 )
           && ( -0.05 < eInvMinusPInv ) && ( eInvMinusPInv < 0.005 )
           && ( el.full5x5_sigmaIetaIeta() < 0.030 )
           );
      }
    }
    ret.add("TrigEmu", passTrigEmu);

    ret.add("NMissInner", elValid ? el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) : 0 );
    ret.add("passConversionVeto", el.passConversionVeto());
    ret.add("ThreeChargeAgreement", el.isGsfCtfScPixChargeConsistent()); // ??

    return ret;
  }
};

/**
 * Additional muon identification variables
 */
class DictMuonIDVars : public DictTool<pat::Muon>, protected DictPVHelper {
public:
  DictMuonIDVars(const edm::ParameterSet& config)
    : DictTool<pat::Muon>(config)
    , DictPVHelper(config)
  {}
  virtual ~DictMuonIDVars() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    DictPVHelper            ::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
  }

  virtual Dict evaluate(const pat::Muon& mu,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    bool muValid{mu.originalObjectRef().isNonnull()};
    const reco::Vertex* pv = event ? getPV(event) : nullptr;

    Dict ret;
    // Same values used for cut-based muon ID. See:
    //     https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_15/DataFormats/MuonReco/src/MuonSelectors.cc#L756
    ret.add("dxy", pv && muValid && mu.muonBestTrack().isNonnull() ? mu.muonBestTrack()->dxy(pv->position()) : 0. );
    ret.add("dz" , pv && muValid && mu.muonBestTrack().isNonnull() ? mu.muonBestTrack()->dz (pv->position()) : 0. );
    ret.add("dca", mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D));

    ret.add("DPToPT", muValid && mu.innerTrack().isNonnull() ? mu.innerTrack()->ptError()/mu.innerTrack()->pt() : 0);

    return ret;
  }
};

/**
 * Additional jet identification variables
 */
class DictJetIDVars : public DictTool<pat::Jet> {
public:
  DictJetIDVars(const edm::ParameterSet& config)
    : DictTool<pat::Jet>(config)
  {
    const auto& bTagCfg = config.getParameter<edm::ParameterSet>("bTaggers");
    for ( const auto& bTagNm : bTagCfg.getParameterNames() ) {
      m_bTaggers.emplace(bTagNm, bTagCfg.getParameter<std::string>(bTagNm));
    }
  }
  virtual ~DictJetIDVars() {}

  virtual Dict evaluate(const pat::Jet& jet,
      const edm::Event* /**/, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    Dict ret;
    bool valid{jet.originalObjectRef().isNonnull()};
    ret.add("jecFactor"   , valid ? jet.jecFactor(0) : 0);
    ret.add("area"        , jet.jetArea());
    ret.add("partonFlavor", jet.partonFlavour());
    ret.add("hadronFlavor", jet.hadronFlavour());
    ret.add("puJetID"     , valid ? jet.userFloat("pileupJetId:fullDiscriminant") : 0);
    for ( const auto& varNm : m_bTaggers ) {
      ret.add(varNm.first, valid ? jet.bDiscriminator(varNm.second) : 0);
    }
    return ret;
  }
private:
  std::map<std::string,std::string> m_bTaggers;
};
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictElectronIDVars, "ttw_electronIDVars");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictMuonIDVars    , "ttw_muonIDVars");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Jet     >::factory, TTWAnalysis::DictJetIDVars     , "ttw_jetIDVars");
