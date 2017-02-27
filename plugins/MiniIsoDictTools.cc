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
    m_isoComp.doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    // TODO electrons and muons unless vetoing "None"
  }

  virtual Dict evaluate(edm::Ptr<pat::Electron> cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  mutable heppy::IsolationComputer m_isoComp;
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
    m_isoComp.doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    // TODO electrons and muons unless vetoing "None"
  }

  virtual Dict evaluate(edm::Ptr<pat::Muon> cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  mutable heppy::IsolationComputer m_isoComp;
  EffectiveAreas m_ea;
};

}

// implementations

TTWAnalysis::Dict TTWAnalysis::DictElectronMiniIsolation::evaluate(edm::Ptr<pat::Electron> cand,
    const edm::Event* event, const edm::EventSetup* /**/,
    const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  using heppy::IsolationComputer;
  m_isoComp.updateEvent(event);

  const bool valid{cand.isNonnull() && cand->superCluster().isNonnull()};
  if ( ! valid ) {
    LogDebug("ttW-eleSync") << "Null electron: candidate is " << ( cand.isNonnull() ? "not " : "" ) << "null; supercluster ref is " << ( cand.isNonnull() ? ( cand->superCluster().isNonnull() ? "also not " : "" ) : "NA " ) << "null";
  }
  const double rho = event ? getRho(event) : 0.;
  const double eta = valid ? cand->superCluster()->eta() : 0.;

  const double outerR{valid ? 10.0/std::min(std::max(cand->pt(), 50.),200.) : 0.};

  double absIsoCharged{-1.}, isoPhotRaw{-1.}, isoNHadRaw{-1.}, absIsoPU{-1.}, isoNeutralRaw{-1.}, isoNeutralWeights{-1.}, isoNeutralRhoArea{-1.}, isoNeutralDeltaBeta{-1.};
  if ( valid ) {
    const pat::Electron& candR = *cand;
    const double innerRPh{(!candR.isEB()) ? 0.08 : 0.}; // inner rad for photons if eleE
    const double innerRCh{(!candR.isEB()) ? 0.015 : 0.}; // inner rad for charged particles if eleE

    absIsoCharged = m_isoComp.chargedAbsIso(candR, outerR, innerRCh, 0., IsolationComputer::selfVetoNone);
    absIsoPU = m_isoComp.puAbsIso(candR, outerR, innerRCh, 0.0, IsolationComputer::selfVetoNone);
    isoPhotRaw = m_isoComp.photonAbsIsoRaw(candR, outerR, innerRPh, 0., IsolationComputer::selfVetoNone);
    isoNHadRaw = m_isoComp.neutralHadAbsIsoRaw(candR, outerR, 0., 0., IsolationComputer::selfVetoNone);
    isoNeutralRaw = isoPhotRaw+isoNHadRaw;
    // puCorr "weights" case
    isoNeutralWeights =
          m_isoComp.photonAbsIsoWeighted(candR, outerR, innerRPh, 0., IsolationComputer::selfVetoNone)
        + m_isoComp.neutralHadAbsIsoWeighted(candR, outerR, 0., 0., IsolationComputer::selfVetoNone)
        ;
    // puCorr "rhoArea" case
    isoNeutralRhoArea = std::max(0., isoNeutralRaw - rho * m_ea.getEffectiveArea(eta) * std::pow(outerR/0.3, 2));
    // puCorr "deltaBeta" case
    isoNeutralDeltaBeta = std::max(0., isoNeutralRaw - 0.5*absIsoPU);
  }
  //
  TTWAnalysis::Dict ret{};
  ret.add("miniIso_R", outerR);
  ret.add("miniIso_AbsCharged", absIsoCharged);
  ret.add("miniIso_AbsPho"    , isoPhotRaw);
  ret.add("miniIso_AbsNHad"   , isoNHadRaw);
  ret.add("miniIso_AbsPU"     , absIsoPU);
  //
  const float candpt = valid ? cand->pt() : -1.;
  auto addNeutralAbsRel = [&ret,absIsoCharged,candpt] (double neuIsoAbs, std::string postfix)
  {
    ret.add("miniIso_AbsNeutral_"+postfix, neuIsoAbs);
    ret.add("miniIso_Abs_"+postfix, absIsoCharged+neuIsoAbs);
    ret.add("miniIso_Rel_"+postfix, (absIsoCharged+neuIsoAbs)/candpt);
  };
  addNeutralAbsRel(isoNeutralWeights  , "weights"  );
  addNeutralAbsRel(isoNeutralRaw      , "raw"      );
  addNeutralAbsRel(isoNeutralRhoArea  , "rhoArea"  );
  addNeutralAbsRel(isoNeutralDeltaBeta, "deltaBeta");
  return ret;
}

TTWAnalysis::Dict TTWAnalysis::DictMuonMiniIsolation:: evaluate(edm::Ptr<pat::Muon> cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  using heppy::IsolationComputer;
  m_isoComp.updateEvent(event);

  const bool valid{cand.isNonnull()};
  const double rho = event ? getRho(event) : 0.;

  const double outerR{valid ? 10.0/std::min(std::max(cand->pt(), 50.),200.) : 0.};

  double absIsoCharged{-1.}, absIsoPU{-1.}, isoNeutralRaw{-1.}, isoNeutralWeights{-1.}, isoNeutralRhoArea{-1.}, isoNeutralDeltaBeta{-1.};
  if ( valid ) {
    const pat::Muon& candR = *cand;
    absIsoCharged = m_isoComp.chargedAbsIso(candR, outerR, 0.0001, 0.);
    absIsoPU = m_isoComp.puAbsIso(candR, outerR, .01, .5);
    isoNeutralRaw = m_isoComp.neutralAbsIsoRaw(candR, outerR, 0.01, 0.5);
    // puCorr "weights" case
    isoNeutralWeights = m_isoComp.neutralAbsIsoWeighted(candR, outerR, 0.01, 0.5);
    // puCorr "rhoArea" case
    isoNeutralRhoArea = std::max(0., isoNeutralRaw - rho * m_ea.getEffectiveArea(candR.eta()) * std::pow(outerR/0.3, 2));
    // puCorr "deltaBeta" case
    isoNeutralDeltaBeta = std::max(0., isoNeutralRaw - 0.5*absIsoPU);
  }
  //
  TTWAnalysis::Dict ret{};
  ret.add("miniIso_R", outerR);
  ret.add("miniIso_AbsCharged", absIsoCharged);
  ret.add("miniIso_AbsPU"     , absIsoPU);
  //
  const float candpt = valid ? cand->pt() : -1.;
  auto addNeutralAbsRel = [&ret,absIsoCharged,candpt] (double neuIsoAbs, std::string postfix)
  {
    ret.add("miniIso_AbsNeutral_"+postfix, neuIsoAbs);
    ret.add("miniIso_Abs_"+postfix, absIsoCharged+neuIsoAbs);
    ret.add("miniIso_Rel_"+postfix, (absIsoCharged+neuIsoAbs)/candpt);
  };
  addNeutralAbsRel(isoNeutralWeights  , "weights"  );
  addNeutralAbsRel(isoNeutralRaw      , "raw"      );
  addNeutralAbsRel(isoNeutralRhoArea  , "rhoArea"  );
  addNeutralAbsRel(isoNeutralDeltaBeta, "deltaBeta");
  return ret;
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictElectronMiniIsolation, "ttw_electronMiniIso");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictMuonMiniIsolation    , "ttw_muonMiniIso");
