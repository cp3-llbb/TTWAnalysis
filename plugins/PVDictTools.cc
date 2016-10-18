#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "cp3_llbb/Framework/interface/Types.h"
#include "cp3_llbb/TTWAnalysis/interface/DictTool.h"

#include "Helpers.h"

namespace TTWAnalysis {
/**
 * Impact parameter (significance) variables for electrons
 */
class DictElectronPVVars : public DictTool<pat::Electron>, protected DictPVHelper {
public:
  DictElectronPVVars(const edm::ParameterSet& config)
    : DictTool<pat::Electron>(config)
    , DictPVHelper(config)
  {}
  virtual ~DictElectronPVVars() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    DictPVHelper                ::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
  }

  virtual Dict evaluate(const pat::Electron& el,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    bool elValid{el.originalObjectRef().isNonnull()};
    const reco::Vertex* pv = event ? getPV(event) : nullptr;

    Dict ret;
    // Same values used for cut-based electron ID. See:
    //     https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_15/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleDzCut.cc#L64
    //     https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_15/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleDxyCut.cc#L64
    ret.add("dxy", pv && elValid ? el.gsfTrack()->dxy(pv->position()) : 0. );
    ret.add("dz" , pv && elValid ? el.gsfTrack()->dz (pv->position()) : 0. );
    ret.add("dca", el.dB(pat::Electron::PV3D)/el.edB(pat::Electron::PV3D));

    return ret;
  }
};
/**
 * Impact parameter (significance) variables for muons
 */
class DictMuonPVVars : public DictTool<pat::Muon>, protected DictPVHelper {
public:
  DictMuonPVVars(const edm::ParameterSet& config)
    : DictTool<pat::Muon>(config)
    , DictPVHelper(config)
  {}
  virtual ~DictMuonPVVars() {}

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

    return ret;
  }
};
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictElectronPVVars, "ttw_electronPVVars");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictMuonPVVars    , "ttw_muonPVVars");
