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
 * Mini-isolation, see https://indico.cern.ch/event/388718/contributions/921752/attachments/777177/1065760/SUS_miniISO_4-21-15.pdf
 * Implementation from:
 * - https://github.com/CERN-PH-CMG/cmg-cmssw/blob/fc44bf40c03c475725d06d1d1f68c2e594c49eff/PhysicsTools/Heppy/python/analyzers/objects/LeptonAnalyzer.py
 * - https://github.com/CERN-PH-CMG/cmg-cmssw/blob/ebe350c7d742c6b29b55e60f5bf872bfc3c5afda/PhysicsTools/Heppy/interface/IsolationComputer.h
 * - https://github.com/CERN-PH-CMG/cmg-cmssw/blob/ebe350c7d742c6b29b55e60f5bf872bfc3c5afda/PhysicsTools/Heppy/src/IsolationComputer.cc
 */
class DictElectronMiniIsolation : public DictTool<pat::Electron>, protected DictRhoHelper {
public:
  DictElectronMiniIsolation(const edm::ParameterSet& config)
    : DictTool<pat::Electron>(config)
    , DictRhoHelper{config}
    , m_isoComp{0.4} // for "weights" case - only used with the "*Weighted" methods, so can use one IsolationComputer for all cases, see below
    , m_ea{EffectiveAreas{config.getUntrackedParameter<edm::FileInPath>("ea").fullPath()}} // R03
  {}

  virtual ~DictElectronMiniIsolation() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    DictRhoHelper::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    m_packedCandidates_token = collector.consumes<std::vector<pat::PackedCandidate>>(config.getParameter<edm::InputTag>("packedCandidates"));
    // TODO electrons and muons unless vetoing "None"
  }

  virtual Dict evaluate(const pat::Electron& cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  mutable edm::EventID m_lastEvent;
  mutable heppy::IsolationComputer m_isoComp;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> m_packedCandidates_token;
  EffectiveAreas m_ea;
};

class DictMuonMiniIsolation : public DictTool<pat::Muon>, protected DictRhoHelper {
public:
  DictMuonMiniIsolation(const edm::ParameterSet& config)
    : DictTool<pat::Muon>(config)
    , DictRhoHelper{config}
    , m_isoComp{0.4} // for "weights" case - only used with the "*Weighted" methods, so can use one IsolationComputer for all cases, see below
    , m_ea{EffectiveAreas{config.getUntrackedParameter<edm::FileInPath>("ea").fullPath()}} // R03
  {}

  virtual ~DictMuonMiniIsolation() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    DictRhoHelper::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    m_packedCandidates_token = collector.consumes<std::vector<pat::PackedCandidate>>(config.getParameter<edm::InputTag>("packedCandidates"));
    // TODO electrons and muons unless vetoing "None"
  }

  virtual Dict evaluate(const pat::Muon& cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  mutable edm::EventID m_lastEvent;
  mutable heppy::IsolationComputer m_isoComp;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> m_packedCandidates_token;
  EffectiveAreas m_ea;
};

}

// implementations

TTWAnalysis::Dict TTWAnalysis::DictElectronMiniIsolation::evaluate(const pat::Electron& cand,
    const edm::Event* event, const edm::EventSetup* /**/,
    const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  if ( event && ( m_lastEvent != event->id() ) ) {
    edm::Handle<std::vector<pat::PackedCandidate>> packedCandidatesHandle;
    event->getByToken(m_packedCandidates_token, packedCandidatesHandle);
    m_isoComp.setPackedCandidates(*(packedCandidatesHandle.product()));
    // TODO if vetoing "any": add all leptons from containers
    // TODO if vetoing "inclusive": add specifically pre-selected leptons (and for these add the miniIso)
    m_lastEvent = event->id();
  }
  using heppy::IsolationComputer;

  double rho = event ? getRho(event) : 0.;
  double eta = cand.originalObjectRef().isNonnull() ? cand.superCluster()->eta() : 0.;

  double outerR{10.0/std::min(std::max(cand.pt(), 50.),200.)};
  double innerRPh{(!cand.isEB()) ? 0.08 : 0.}; // inner rad for photons if eleE
  double innerRCh{(!cand.isEB()) ? 0.015 : 0.}; // inner rad for charged particles if eleE

  double absIsoCharged{m_isoComp.chargedAbsIso(cand, outerR, innerRCh, 0., IsolationComputer::selfVetoNone)};
  double absIsoPU{m_isoComp.puAbsIso(cand, outerR, innerRCh, 0.0, IsolationComputer::selfVetoNone)};
  double isoPhotRaw{m_isoComp.photonAbsIsoRaw(cand, outerR, innerRPh, 0., IsolationComputer::selfVetoNone)};
  double isoNHadRaw{m_isoComp.neutralHadAbsIsoRaw(cand, outerR, 0., 0., IsolationComputer::selfVetoNone)};
  double isoNeutralRaw{isoPhotRaw+isoNHadRaw};
  // puCorr "weights" case
  double isoNeutralWeights{
        m_isoComp.photonAbsIsoWeighted(cand, outerR, innerRPh, 0., IsolationComputer::selfVetoNone)
      + m_isoComp.neutralHadAbsIsoWeighted(cand, outerR, 0., 0., IsolationComputer::selfVetoNone)
      };
  // puCorr "rhoArea" case
  double isoNeutralRhoArea{std::max(0., isoNeutralRaw - rho * m_ea.getEffectiveArea(eta) * std::pow(outerR/0.3, 2))};
  // puCorr "deltaBeta" case
  double isoNeutralDeltaBeta{std::max(0., isoNeutralRaw - 0.5*absIsoPU)};
  //
  TTWAnalysis::Dict ret{};
  ret.add("miniIsoR", outerR);
  ret.add("miniAbsIsoCharged", absIsoCharged);
  ret.add("miniAbsIsoPho"    , isoPhotRaw);
  ret.add("miniAbsIsoNHad"   , isoNHadRaw);
  ret.add("miniAbsIsoPU"     , absIsoPU);
  //
  auto addNeutralAbsRel = [&ret,absIsoCharged,&cand] (double neuIsoAbs, std::string postfix)
  {
    ret.add("miniAbsIsoNeutral_"+postfix, neuIsoAbs);
    ret.add("miniAbsIso_"+postfix, absIsoCharged+neuIsoAbs);
    ret.add("miniRelIso_"+postfix, (absIsoCharged+neuIsoAbs)/cand.pt());
  };
  addNeutralAbsRel(isoNeutralWeights  , "weights"  );
  addNeutralAbsRel(isoNeutralRaw      , "raw"      );
  addNeutralAbsRel(isoNeutralRhoArea  , "rhoArea"  );
  addNeutralAbsRel(isoNeutralDeltaBeta, "deltaBeta");
  return ret;
}

TTWAnalysis::Dict TTWAnalysis::DictMuonMiniIsolation:: evaluate(const pat::Muon& cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  if ( event && ( m_lastEvent != event->id() ) ) {
    edm::Handle<std::vector<pat::PackedCandidate>> packedCandidatesHandle;
    event->getByToken(m_packedCandidates_token, packedCandidatesHandle);
    m_isoComp.setPackedCandidates(*(packedCandidatesHandle.product()));
    // TODO if vetoing "any": add all leptons from containers
    // TODO if vetoing "inclusive": add specifically pre-selected leptons (and for these add the miniIso)
    m_lastEvent = event->id();
  }
  using heppy::IsolationComputer;

  const double rho = event ? getRho(event) : 0.;

  const double outerR{10.0/std::min(std::max(cand.pt(), 50.),200.)};

  const double absIsoCharged{m_isoComp.chargedAbsIso(cand, outerR, 0.0001, 0.)};
  const double absIsoPU{m_isoComp.puAbsIso(cand, outerR, .01, .5)};
  const double isoNeutralRaw{m_isoComp.neutralAbsIsoRaw(cand, outerR, 0.01, 0.5)};
  // puCorr "weights" case
  const double isoNeutralWeights{m_isoComp.neutralAbsIsoWeighted(cand, outerR, 0.01, 0.5)};
  // puCorr "rhoArea" case
  const double isoNeutralRhoArea{std::max(0., isoNeutralRaw - rho * m_ea.getEffectiveArea(cand.eta()) * std::pow(outerR/0.3, 2))};
  // puCorr "deltaBeta" case
  const double isoNeutralDeltaBeta{std::max(0., isoNeutralRaw - 0.5*absIsoPU)};
  //
  TTWAnalysis::Dict ret{};
  ret.add("miniIsoR", outerR);
  ret.add("miniAbsIsoCharged", absIsoCharged);
  ret.add("miniAbsIsoPU"     , absIsoPU);
  //
  auto addNeutralAbsRel = [&ret,absIsoCharged,&cand] (double neuIsoAbs, std::string postfix)
  {
    ret.add("miniAbsIsoNeutral_"+postfix, neuIsoAbs);
    ret.add("miniAbsIso_"+postfix, absIsoCharged+neuIsoAbs);
    ret.add("miniRelIso_"+postfix, (absIsoCharged+neuIsoAbs)/cand.pt());
  };
  addNeutralAbsRel(isoNeutralWeights  , "weights"  );
  addNeutralAbsRel(isoNeutralRaw      , "raw"      );
  addNeutralAbsRel(isoNeutralRhoArea  , "rhoArea"  );
  addNeutralAbsRel(isoNeutralDeltaBeta, "deltaBeta");
  return ret;
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictElectronMiniIsolation, "ttw_electronMiniIso");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictMuonMiniIsolation    , "ttw_muonMiniIso");
