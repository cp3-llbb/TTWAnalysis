#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "cp3_llbb/Framework/interface/Types.h"
#include "cp3_llbb/TTWAnalysis/interface/DictTool.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

#include "Helpers.h"

namespace TTWAnalysis {
/**
 * Lepton iso (like LeptonsProducer)
 */
template<class Lepton>
class DictLeptonIso : public DictTool<Lepton> {
public:
  DictLeptonIso(const edm::ParameterSet& config)
    : DictTool<Lepton>(config)
  {
    m_ea.emplace("R04", EffectiveAreas(config.getUntrackedParameter<edm::FileInPath>("ea_R04").fullPath()));
    m_ea.emplace("R03", EffectiveAreas(config.getUntrackedParameter<edm::FileInPath>("ea_R03").fullPath()));
  }
  virtual ~DictLeptonIso() {}
  // no implementation for interface method
protected:
  void fillIsolations(float chargedHadronIso, float neutralHadronIso, float photonIso, float puChargedHadronIso,
                      float pt, float eta, float rho, const std::string& name, Dict& out) const
  {
    out.add("chargedHadronIso"  +name, chargedHadronIso  );
    out.add("neutralHadronIso"  +name, neutralHadronIso  );
    out.add("photonHadronIso"   +name, photonIso         );
    out.add("puChargedHadronIso"+name, puChargedHadronIso);
    out.add("relativeIso"       +name, (chargedHadronIso + neutralHadronIso + photonIso) / pt);
    out.add("relativeIso"+name+"_deltaBeta", (chargedHadronIso + std::max((neutralHadronIso + photonIso) - 0.5f * puChargedHadronIso, 0.0f)) / pt);

    float ea = m_ea.find(name)->second.getEffectiveArea(eta);
    out.add("relativeIso"+name+"_withEA"   , (chargedHadronIso + std::max((neutralHadronIso + photonIso) - rho * ea, 0.0f)) / pt);
    out.add("EA"+name, ea);
  }
private:
  // Effective areas
  std::map<std::string,EffectiveAreas> m_ea;
};

/**
 * Electron isolation
 */
class DictElectronIso : public DictLeptonIso<pat::Electron>, protected DictRhoHelper {
public:
  DictElectronIso(const edm::ParameterSet& config)
    : DictLeptonIso<pat::Electron>(config)
    , DictRhoHelper(config)
  {}
  virtual ~DictElectronIso() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    DictLeptonIso<pat::Electron>::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    DictRhoHelper               ::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
  }

  virtual Dict evaluate(edm::Ptr<pat::Electron> el,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    const double rho = event ? getRho(event) : 0.;
    const bool valid{el.isNonnull() && el->superCluster().isNonnull()};
    const double pt = valid ? el->pt() : -1.;
    const double eta = valid ? el->superCluster()->eta() : 0.;

    Dict ret{};
    reco::GsfElectron::PflowIsolationVariables pfIso = valid ? el->pfIsolationVariables() : reco::GsfElectron::PflowIsolationVariables{};
    fillIsolations(pfIso.sumChargedHadronPt, pfIso.sumNeutralHadronEt,
                   pfIso.sumPhotonEt, pfIso.sumPUPt,
                   pt, eta, rho,
                   "R03", ret);
    fillIsolations(valid ? el->chargedHadronIso() : -1., valid ? el->neutralHadronIso() : -1.,
                   valid ? el->photonIso() : -1., valid ? el->puChargedHadronIso() : -1.,
                   pt, eta, rho,
                   "R04", ret);

    ret.add("ecalPFClusterIso", valid ? el->ecalPFClusterIso() : -1.);
    ret.add("hcalPFClusterIso", valid ? el->hcalPFClusterIso() : -1.);
    ret.add("trackIso", valid ? el->trackIso() : -1.);

    return ret;
  }
};

class DictMuonIso : public DictLeptonIso<pat::Muon>, protected DictRhoHelper {
public:
  DictMuonIso(const edm::ParameterSet& config)
    : DictLeptonIso<pat::Muon>(config)
    , DictRhoHelper(config)
  {}
  virtual ~DictMuonIso() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    DictLeptonIso<pat::Muon>::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    DictRhoHelper           ::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
  }

  virtual Dict evaluate(edm::Ptr<pat::Muon> mu,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    const double rho = event ? getRho(event) : 0.;
    const bool valid{mu.isNonnull()};
    const double pt  = valid ? mu->pt() : -1.;
    const double eta = valid ? mu->eta() : -10.;

    Dict ret{};
    reco::MuonPFIsolation pfIso = valid ? mu->pfIsolationR03() : reco::MuonPFIsolation{};
    fillIsolations(pfIso.sumChargedHadronPt, pfIso.sumNeutralHadronEt,
                   pfIso.sumPhotonEt, pfIso.sumPUPt,
                   pt, eta, rho,
                   "R03", ret);
    if ( valid ) { pfIso = mu->pfIsolationR04(); }
    fillIsolations(pfIso.sumChargedHadronPt, pfIso.sumNeutralHadronEt,
                   pfIso.sumPhotonEt, pfIso.sumPUPt,
                   pt, eta, rho,
                   "R04", ret);

    return ret;
  }
};
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictElectronIso, "ttw_electronIso");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictMuonIso    , "ttw_muonIso");
