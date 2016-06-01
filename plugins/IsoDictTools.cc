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
 * Small helper to access rho
 */
class DictRhoHelper {
public:
  DictRhoHelper(const edm::ParameterSet& config) {}
  void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector)
  {
    m_rho_token = collector.consumes<double>(config.getUntrackedParameter<edm::InputTag>("rho",
                                             edm::InputTag("fixedGridRhoFastjetAll")));
  }
protected:
  double getRho(const edm::Event* event) const
  {
    edm::Handle<double> handle;
    event->getByToken(m_rho_token, handle);
    return *handle;
  }
private:
  edm::EDGetTokenT<double> m_rho_token;
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

  virtual Dict evaluate(const pat::Electron& el,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    double rho = event ? getRho(event) : 0.;
    bool elValid{el.originalObjectRef().isNonnull()};
    double pt = el.pt();
    double eta = elValid ? el.superCluster()->eta() : 0.;

    Dict ret{};
    reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
    fillIsolations(pfIso.sumChargedHadronPt, pfIso.sumNeutralHadronEt,
                   pfIso.sumPhotonEt, pfIso.sumPUPt,
                   pt, eta, rho,
                   "R03", ret);
    fillIsolations(el.chargedHadronIso(), el.neutralHadronIso(),
                   el.photonIso(), el.puChargedHadronIso(),
                   pt, eta, rho,
                   "R04", ret);

    ret.add("ecalPFClusterIso", el.ecalPFClusterIso());
    ret.add("hcalPFClusterIso", el.hcalPFClusterIso());
    ret.add("trackIso", el.trackIso());

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

  virtual Dict evaluate(const pat::Muon& mu,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    double rho = event ? getRho(event) : 0.;
    double pt = mu.pt();
    double eta = mu.eta();

    Dict ret{};
    reco::MuonPFIsolation pfIso = mu.pfIsolationR03();
    fillIsolations(pfIso.sumChargedHadronPt, pfIso.sumNeutralHadronEt,
                   pfIso.sumPhotonEt, pfIso.sumPUPt,
                   pt, eta, rho,
                   "R03", ret);
    pfIso = mu.pfIsolationR04();
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
