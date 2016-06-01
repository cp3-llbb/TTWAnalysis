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

namespace TTWAnalysis {
/**
 * Fill basic kinematic information: 4-momentum vector, rapidity and charge
 * Equivalent to what is filled by cp3_llbb/Framework's CandidateProducer
 */
template<class RecoCandidate>
class DictKinematic : public DictTool<RecoCandidate> {
public:
  DictKinematic(const edm::ParameterSet& config) : DictTool<RecoCandidate>(config) {}
  virtual ~DictKinematic() {}

  virtual Dict evaluate(const RecoCandidate& cand,
      const edm::Event* /**/, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    Dict ret{};
    ret.add("p4", LorentzVector(cand.pt(), cand.eta(), cand.phi(), cand.energy()));
    ret.add("y" , cand.rapidity());
    ret.add("charge", cand.charge());
    return ret;
  }
};
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictKinematic<pat::Electron>, "ttw_electronKin");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictKinematic<pat::Muon    >, "ttw_muonKin");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Jet     >::factory, TTWAnalysis::DictKinematic<pat::Jet     >, "ttw_jetKin");

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

namespace TTWAnalysis {
/**
 * Configurable selections
 */
template<class RecoCandidate>
class DictHybridCuts : public DictTool<RecoCandidate> {
public:
  DictHybridCuts(const edm::ParameterSet& config)
    : DictTool<RecoCandidate>(config)
  {
    const auto& cutsConfig = config.getParameter<edm::ParameterSet>("Cuts");
    for ( const auto& iName : cutsConfig.getParameterNames() ) {
      m_cuts.emplace(iName, HybridCut{cutsConfig.getParameter<std::string>(iName)});
    }
  }
  virtual ~DictHybridCuts() {}

  virtual Dict evaluate(const RecoCandidate& cand,
      const edm::Event* /**/, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    bool valid{cand.originalObjectRef().isNonnull()};
    Dict ret{};
    for ( const auto& nmAndCut : m_cuts ) {
      ret.add(nmAndCut.first, valid && nmAndCut.second(cand));
    }
    return ret;
  }
private:
  using HybridCut = StringCutObjectSelector<RecoCandidate>;
  std::map<std::string,HybridCut> m_cuts;
};
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictHybridCuts<pat::Electron>, "ttw_electronHybridCuts");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictHybridCuts<pat::Muon    >, "ttw_muonHybridCuts");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Jet     >::factory, TTWAnalysis::DictHybridCuts<pat::Jet     >, "ttw_jetHybridCuts");

#include "CommonTools/Utils/interface/StringObjectFunction.h"

namespace TTWAnalysis {
/**
 * Configurable functions
 */
template<class RecoCandidate>
class DictHybridFunctions : public DictTool<RecoCandidate> {
public:
  DictHybridFunctions(const edm::ParameterSet& config)
    : DictTool<RecoCandidate>(config)
  {
    const auto& funsConfig = config.getParameter<edm::ParameterSet>("Functions");
    for ( const auto& iName : funsConfig.getParameterNames() ) {
      m_funs.emplace(iName, HybridFunction{funsConfig.getParameter<std::string>(iName)});
    }
  }
  virtual ~DictHybridFunctions() {}

  virtual Dict evaluate(const RecoCandidate& cand,
      const edm::Event* /**/, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    Dict ret{};
    bool valid{cand.originalObjectRef().isNonnull()};
    for ( const auto& nmAndFun : m_funs ) {
      ret.add(nmAndFun.first, valid ? nmAndFun.second(cand) : 0.);
    }
    return ret;
  }
private:
  using HybridFunction = StringObjectFunction<RecoCandidate>;
  std::map<std::string,HybridFunction> m_funs;
};
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictHybridFunctions<pat::Electron>, "ttw_electronHybridFunctions");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictHybridFunctions<pat::Muon    >, "ttw_muonHybridFunctions");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Jet     >::factory, TTWAnalysis::DictHybridFunctions<pat::Jet     >, "ttw_jetHybridFunctions");

namespace TTWAnalysis {
/**
 * Gen match for leptons
 */
template<class Lepton>
class DictLeptonGenMatch : public DictTool<Lepton> {
public:
  DictLeptonGenMatch(const edm::ParameterSet& config) : DictTool<Lepton>(config) {}
  virtual ~DictLeptonGenMatch() {}
  virtual Dict evaluate(const Lepton& cand,
      const edm::Event* /**/, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    const auto gen = cand.genParticle();
    Dict ret{};
    ret.add("gen_p4"    , gen ? LorentzVector(gen->pt(), gen->eta(), gen->phi(), gen->energy()) : LorentzVector(0.,0.,0.,0.));
    ret.add("gen_y"     , gen ? gen->y() : 0.);
    ret.add("gen_charge", gen ? gen->charge() : 0);
    return ret;
  }
};
/**
 * Gen match for jets
 */
class DictJetGenMatch : public DictTool<pat::Jet> {
public:
  DictJetGenMatch(const edm::ParameterSet& config) : DictTool<pat::Jet>(config) {}
  virtual ~DictJetGenMatch() {}
  virtual Dict evaluate(const pat::Jet& cand,
      const edm::Event* /**/, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    const auto gen = cand.genJet();
    Dict ret{};
    ret.add("gen_p4"    , gen ? LorentzVector(gen->pt(), gen->eta(), gen->phi(), gen->energy()) : LorentzVector(0.,0.,0.,0.));
    ret.add("gen_y"     , gen ? gen->y() : 0.);
    ret.add("gen_charge", gen ? gen->charge() : 0);
    return ret;
  }
};
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictLeptonGenMatch<pat::Electron>, "ttw_electronGenMatch");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictLeptonGenMatch<pat::Muon    >, "ttw_muonGenMatch");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Jet     >::factory, TTWAnalysis::DictJetGenMatch                  , "ttw_jetGenMatch");

#include "cp3_llbb/TTWAnalysis/interface/NewTypes.h"

namespace TTWAnalysis {
// Basic indices of the different candidates
/**
 * TTWAnalysis::Lepton minimal info
 */
class LeptonCandidate : public DictTool<TTWAnalysis::Lepton> {
public:
  LeptonCandidate(const edm::ParameterSet& config)
    : DictTool<TTWAnalysis::Lepton>(config)
  {}
  virtual ~LeptonCandidate() {}

  virtual Dict evaluate(const TTWAnalysis::Lepton& lepton,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    Dict ret{};
    ret.add("isEl", lepton.isEl());
    ret.add("isMu", lepton.isMu());
    ret.add("idx" , lepton.idx);
    return ret;
  }
};
/**
 * TTWAnalysis::DiLepton minimal info
 */
class DiLeptonCandidate : public DictTool<TTWAnalysis::DiLepton> {
public:
  DiLeptonCandidate(const edm::ParameterSet& config)
    : DictTool<TTWAnalysis::DiLepton>(config)
  {}
  virtual ~DiLeptonCandidate() {}

  virtual Dict evaluate(const TTWAnalysis::DiLepton& ll,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    bool valid{ll.first && ll.second};
    Dict ret{};
    // essential
    ret.add("lidx1" , ll.lidxs.first );
    ret.add("lidx2" , ll.lidxs.second);
    // convenience
    ret.add("isElEl", valid && ll.isElEl());
    ret.add("isElMu", valid && ll.isElMu());
    ret.add("isMuEl", valid && ll.isMuEl());
    ret.add("isMuMu", valid && ll.isMuMu());
    ret.add("idx1"  , valid && ll.idxs().first );
    ret.add("idx2"  , valid && ll.idxs().second);
    ret.add("isOS"  , valid && ll.isOS());
    ret.add("isSF"  , valid && ll.isSF());
    return ret;
  }
};
/**
 * TTWAnalysis::DiJet minimal info
 */
class DiJetCandidate : public DictTool<TTWAnalysis::DiJet> {
public:
  DiJetCandidate(const edm::ParameterSet& config)
    : DictTool<TTWAnalysis::DiJet>(config)
  {}
  virtual ~DiJetCandidate() {}

  virtual Dict evaluate(const TTWAnalysis::DiJet& jj,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    Dict ret{};
    ret.add("idx1", jj.jidxs.first);
    ret.add("idx2", jj.jidxs.second);
    return ret;
  }
};
/**
 * TTWAnalysis::DiLeptonDiJet minimal info
 */
class DiLeptonDiJetCandidate : public DictTool<TTWAnalysis::DiLeptonDiJet> {
public:
  DiLeptonDiJetCandidate(const edm::ParameterSet& config)
    : DictTool<TTWAnalysis::DiLeptonDiJet>(config)
  {}
  virtual ~DiLeptonDiJetCandidate() {}

  virtual Dict evaluate(const TTWAnalysis::DiLeptonDiJet& lljj,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    Dict ret{};
    // essential
    ret.add("llidx", lljj.llIdx);
    ret.add("jjidx", lljj.jjIdx);
    // convenience
    bool valid{lljj.jj && lljj.ll};
    ret.add("jidx1", valid ? lljj.jj->jidxs.first  : -1);
    ret.add("jidx2", valid ? lljj.jj->jidxs.second : -1);
    ret.add("lidx1", valid ? lljj.ll->lidxs.first  : -1);
    ret.add("lidx2", valid ? lljj.ll->lidxs.second : -1);
    return ret;
  }
};
/**
 * TTWAnalysis::DiLeptonDiJetMet minimal info
 */
class DiLeptonDiJetMetCandidate : public DictTool<TTWAnalysis::DiLeptonDiJetMet> {
public:
  DiLeptonDiJetMetCandidate(const edm::ParameterSet& config)
    : DictTool<TTWAnalysis::DiLeptonDiJetMet>(config)
  {}
  virtual ~DiLeptonDiJetMetCandidate() {}

  virtual Dict evaluate(const TTWAnalysis::DiLeptonDiJetMet& lljjm,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    Dict ret{};
    // essential
    ret.add("lljjidx", lljjm.lljjIdx);
    // convenience
    bool valid{lljjm.lljj};
    ret.add("llidx", valid ? lljjm.lljj->llIdx : -1);
    ret.add("jjidx", valid ? lljjm.lljj->jjIdx : -1);
    valid = valid && ( lljjm.lljj->jj && lljjm.lljj->ll );
    ret.add("jidx1", valid ? lljjm.lljj->jj->jidxs.first  : -1);
    ret.add("jidx2", valid ? lljjm.lljj->jj->jidxs.second : -1);
    ret.add("lidx1", valid ? lljjm.lljj->ll->lidxs.first  : -1);
    ret.add("lidx2", valid ? lljjm.lljj->ll->lidxs.second : -1);
    return ret;
  }
};
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<TTWAnalysis::Lepton          >::factory, TTWAnalysis::LeptonCandidate          , "ttw_leptonCandidate");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<TTWAnalysis::DiLepton        >::factory, TTWAnalysis::DiLeptonCandidate        , "ttw_dileptonCandidate");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<TTWAnalysis::DiJet           >::factory, TTWAnalysis::DiJetCandidate           , "ttw_dijetCandidate");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<TTWAnalysis::DiLeptonDiJet   >::factory, TTWAnalysis::DiLeptonDiJetCandidate   , "ttw_dileptondijetCandidate");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<TTWAnalysis::DiLeptonDiJetMet>::factory, TTWAnalysis::DiLeptonDiJetMetCandidate, "ttw_dileptondijetmetCandidate");

#include "cp3_llbb/TTWAnalysis/interface/HLTMatch.h"
#include "cp3_llbb/Framework/interface/ProducersManager.h"
#include "cp3_llbb/Framework/interface/HLTProducer.h"

namespace TTWAnalysis {
template<class Candidate>
class DictHLTMatch : public DictTool<Candidate> {
public:
  DictHLTMatch(const edm::ParameterSet& config)
    : DictTool<Candidate>(config)
  {
    const auto& selsConfig = config.getParameter<edm::ParameterSet>("Selections");
    for ( const auto& iName : selsConfig.getParameterNames() ) {
      m_matchers.emplace(iName, HLTMatch<Candidate>{selsConfig.getParameter<std::vector<std::string>>(iName)});
    }
  }
  virtual ~DictHLTMatch() {}

  virtual Dict evaluate(const Candidate& cand,
      const edm::Event* /**/, const edm::EventSetup* /**/,
      const ProducersManager* producers, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    if ( producers ) {
      const HLTProducer& hlt = producers->get<HLTProducer>("hlt");
      for ( auto& nmMatcher : m_matchers ) { nmMatcher.second.setHLTProd(&hlt); }
    }

    Dict ret{};
    for ( const auto& nmMatcher : m_matchers ) {
      ret.add(nmMatcher.first, nmMatcher.second(cand));
    }
    return ret;
  }
private:
  mutable std::map<std::string,HLTMatch<Candidate>> m_matchers;
};
}
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<TTWAnalysis::Lepton  >::factory, TTWAnalysis::DictHLTMatch<TTWAnalysis::Lepton>  , "ttw_leptonHLTMatch");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<TTWAnalysis::DiLepton>::factory, TTWAnalysis::DictHLTMatch<TTWAnalysis::DiLepton>, "ttw_dileptonHLTMatch");
