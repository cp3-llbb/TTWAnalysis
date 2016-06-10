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

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

#include "Helpers.h"
#include "Math/VectorUtil.h"
#include "TMVA/Reader.h"

namespace TTWAnalysis {
/**
 * (prompt) lepton multivariate ID
 *
 * See
 * - https://twiki.cern.ch/twiki/bin/viewauth/CMS/LeptonMVA
 * - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTH-Multileptons
 */
class DictTTHElectronMVA : public DictTool<pat::Electron>, protected DictRhoHelper, protected DictPVHelper {
public:
  DictTTHElectronMVA(const edm::ParameterSet& config)
    : DictTool<pat::Electron>{config}
    , DictRhoHelper{config}
    , DictPVHelper{config}
    , m_vars{
        {"LepGood_pt"                   , -1.}
      , {"LepGood_eta"                  , -1.}
      , {"LepGood_jetNDauChargedMVASel" , -1.}
      , {"LepGood_miniRelIsoCharged"    , -1.}
      , {"LepGood_miniRelIsoNeutral"    , -1.}
      , {"LepGood_jetPtRelv2"           , -1.}
      , {"min(LepGood_jetPtRatiov2,1.5)", -1.}
      , {"max(LepGood_jetBTagCSV,0)"    , -1.}
      , {"LepGood_sip3d"                , -1.}
      , {"log(abs(LepGood_dxy))"        , -1.}
      , {"log(abs(LepGood_dz))"         , -1.}
      , {"LepGood_mvaIdSpring15"        , -1.}
    }
    , m_ea{EffectiveAreas{config.getParameter<edm::FileInPath>("ea").fullPath()}}
    , m_maxDR{config.getParameter<double>("JetLeptonDR")}
  {
    for ( auto& entry : m_vars ) {
      m_tmvaReader.AddVariable(entry.first, &(entry.second));
    }
    m_tmvaReader.BookMVA("ttH-el-BDTG", config.getParameter<edm::FileInPath>("WeightsFile").fullPath());
  }

  virtual ~DictTTHElectronMVA() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    DictRhoHelper::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    DictPVHelper::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    m_packedCandidates_token = collector.consumes<std::vector<pat::PackedCandidate>>(config.getParameter<edm::InputTag>("packedCandidates"));
    m_jetsToken = collector.consumes<std::vector<pat::Jet>>(config.getParameter<edm::InputTag>("Jets"));
  }

  virtual Dict evaluate(const pat::Electron& cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  mutable edm::EventID m_lastEvent;
  mutable heppy::IsolationComputer m_isoComp;
  mutable std::vector<std::pair<std::string,float>> m_vars;
  float& getVar(const std::string& name) const { auto it = std::find_if(std::begin(m_vars), std::end(m_vars), [&name] ( const std::pair<std::string,float>& item ) { return name == item.first; }); return it->second; }
  mutable TMVA::Reader m_tmvaReader;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> m_packedCandidates_token;
  edm::EDGetTokenT<std::vector<pat::Jet>> m_jetsToken;
  EffectiveAreas m_ea;
  double m_maxDR;
};

class DictTTHMuonMVA : public DictTool<pat::Muon>, protected DictRhoHelper, protected DictPVHelper {
public:
  DictTTHMuonMVA(const edm::ParameterSet& config)
    : DictTool<pat::Muon>{config}
    , DictRhoHelper{config}
    , DictPVHelper{config}
    , m_vars{
        {"LepGood_pt"                   , -1.}
      , {"LepGood_eta"                  , -1.}
      , {"LepGood_jetNDauChargedMVASel" , -1.}
      , {"LepGood_miniRelIsoCharged"    , -1.}
      , {"LepGood_miniRelIsoNeutral"    , -1.}
      , {"LepGood_jetPtRelv2"           , -1.}
      , {"min(LepGood_jetPtRatiov2,1.5)", -1.}
      , {"max(LepGood_jetBTagCSV,0)"    , -1.}
      , {"LepGood_sip3d"                , -1.}
      , {"log(abs(LepGood_dxy))"        , -1.}
      , {"log(abs(LepGood_dz))"         , -1.}
      , {"LepGood_segmentCompatibility" , -1.}
    }
    , m_ea{EffectiveAreas{config.getParameter<edm::FileInPath>("ea").fullPath()}}
    , m_maxDR{config.getParameter<double>("JetLeptonDR")}
  {
    for ( auto& entry : m_vars ) {
      m_tmvaReader.AddVariable(entry.first, &(entry.second));
    }
    m_tmvaReader.BookMVA("ttH-mu-BDTG", config.getParameter<edm::FileInPath>("WeightsFile").fullPath());
  }

  virtual ~DictTTHMuonMVA() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    DictRhoHelper::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    DictPVHelper::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    m_packedCandidates_token = collector.consumes<std::vector<pat::PackedCandidate>>(config.getParameter<edm::InputTag>("packedCandidates"));
    m_jetsToken = collector.consumes<std::vector<pat::Jet>>(config.getParameter<edm::InputTag>("Jets"));
  }

  virtual Dict evaluate(const pat::Muon& cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  mutable edm::EventID m_lastEvent;
  mutable heppy::IsolationComputer m_isoComp;
  mutable std::vector<std::pair<std::string,float>> m_vars;
  float& getVar(const std::string& name) const { auto it = std::find_if(std::begin(m_vars), std::end(m_vars), [&name] ( const std::pair<std::string,float>& item ) { return name == item.first; }); return it->second; }
  mutable TMVA::Reader m_tmvaReader;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> m_packedCandidates_token;
  edm::EDGetTokenT<std::vector<pat::Jet>> m_jetsToken;
  EffectiveAreas m_ea;
  double m_maxDR;
};
}

// implementation
namespace {
  template<typename _Arg,typename _UnaryFunction>
  struct _LessOf {
    public:
      _LessOf(_UnaryFunction fun) : m_fun(fun) {}
      bool operator() ( const _Arg& __a, const _Arg& __b ) const
      {
        return m_fun(__a) < m_fun(__b);
      }
    private:
      _UnaryFunction m_fun;
  };
  // factory method
  template<typename _Arg,typename _UnaryFunction>
  _LessOf<_Arg,_UnaryFunction> LessOf( _UnaryFunction __fun )
  { return _LessOf<_Arg,_UnaryFunction>(std::forward<_UnaryFunction>(__fun)); }

  const pat::Jet* findMatchingJet( const reco::RecoCandidate& lepton, const std::vector<pat::Jet>& jets, double maxDR )
  {
    using namespace ROOT::Math;
    const auto it = std::min_element(std::begin(jets), std::end(jets),
          LessOf<const pat::Jet&>( [&lepton] ( const pat::Jet& j ) {
            return VectorUtil::DeltaR(j.p4(), lepton.p4());
          } ));
    if ( ( std::end(jets) != it ) && ( VectorUtil::DeltaR(lepton.p4(), it->p4()) < maxDR ) ) {
      return &*it;
    }
    return nullptr;
  }

  // Lepton-aware JEC
  reco::Candidate::LorentzVector jetLepAwareJEC(const pat::Jet* jet, const reco::Candidate::LorentzVector& lepton)
  {
    if ( ( ! jet ) || ( (jet->p4()*jet->jecFactor("Uncorrected")-lepton).Rho() < 1.e-4 ) ) {
      return lepton;
    } else {
      const float rawCorr{jet->jecFactor("Uncorrected")};
      const float l1Corr{jet->jecSetAvailable("L1") ? jet->jecFactor("L1") : float(1.)};
      return (jet->p4()*rawCorr-lepton*(1./l1Corr))/rawCorr+lepton;
    }
  }

  double ptRelv2(const pat::Jet* jet, const reco::Candidate::LorentzVector& lepton)
  {
    reco::Candidate::LorentzVector m{jetLepAwareJEC(jet, lepton)};
    return ( (m-lepton).Rho() < 1.e-4 ) ? 0. : lepton.Vect().Cross((m-lepton).Vect()).R();
  }

  int numberOfChargedDaughtersMVASel(const pat::Jet* jet, const reco::RecoCandidate& cand, const reco::Vertex* pv)
  {
    using namespace ROOT::Math;
    if ( ( ! jet ) || ( ! pv ) ) { return 0; }

    return std::count_if(std::begin(jet->daughterPtrVector()), std::end(jet->daughterPtrVector()),
        [&cand,pv] ( reco::CandidatePtr daug ) {
          if ( ( VectorUtil::DeltaR(daug->p4(), cand.p4()) < 0.4 ) && ( daug->charge() != 0 ) ) {
            const pat::PackedCandidate* pkCand = dynamic_cast<const pat::PackedCandidate*>(daug.get());
            if ( pkCand->fromPV() >= pat::PackedCandidate::PVTight ) {
              auto trk = pkCand->pseudoTrack();
              return ( ( trk.pt() > 1. )
                    && ( trk.hitPattern().numberOfValidHits() >= 8 )
                    && ( trk.hitPattern().numberOfValidPixelHits() >= 2 )
                    && ( trk.normalizedChi2() < 5. )
                    && ( std::fabs(trk.dxy(pv->position())) < 0.2 )
                    && ( std::fabs(trk.dz(pv->position())) < 17. )
                     ) ;
            }
          }
          return false;
        });
  }
}

TTWAnalysis::Dict TTWAnalysis::DictTTHElectronMVA::evaluate(const pat::Electron& cand,
    const edm::Event* event, const edm::EventSetup* /**/,
    const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  using heppy::IsolationComputer;

  const double rho = event ? getRho(event) : 0.;
  bool elValid{cand.originalObjectRef().isNonnull()};
  const double eta = elValid ? cand.superCluster()->eta() : 0.;
  const reco::Vertex* pv = event ? getPV(event) : nullptr;

  // TODO add all this mini isolation info to the pat candidates instead to void duplication

  const double outerR{10.0/std::min(std::max(cand.pt(), 50.),200.)};
  const double innerRPh{(!cand.isEB()) ? 0.08 : 0.}; // inner rad for photons if eleE
  const double innerRCh{(!cand.isEB()) ? 0.015 : 0.}; // inner rad for charged particles if eleE

  const double absIsoCharged{m_isoComp.chargedAbsIso(cand, outerR, innerRCh, 0., IsolationComputer::selfVetoNone)};
  const double isoPhotRaw{m_isoComp.photonAbsIsoRaw(cand, outerR, innerRPh, 0., IsolationComputer::selfVetoNone)};
  const double isoNHadRaw{m_isoComp.neutralHadAbsIsoRaw(cand, outerR, 0., 0., IsolationComputer::selfVetoNone)};
  const double isoNeutral{std::max(0., isoPhotRaw+isoNHadRaw - rho * m_ea.getEffectiveArea(eta) * std::pow(outerR/0.3, 2))};

  edm::Handle<std::vector<pat::Jet>> jets;
  const pat::Jet* jet{nullptr};
  if ( event ) {
    event->getByToken(m_jetsToken, jets);
    jet = findMatchingJet( cand, *jets, m_maxDR );
  }
  int jetNDauChargedMVASel = ( jet && pv ) ? numberOfChargedDaughtersMVASel(jet, cand, pv) : 0;

  getVar("LepGood_pt"                   ) = cand.pt();
  getVar("LepGood_eta"                  ) = eta;
  getVar("LepGood_jetNDauChargedMVASel" ) = jetNDauChargedMVASel;
  getVar("LepGood_miniRelIsoCharged"    ) = absIsoCharged/cand.pt();
  getVar("LepGood_miniRelIsoNeutral"    ) = isoNeutral/cand.pt();
  getVar("LepGood_jetPtRelv2"           ) = ptRelv2(jet, cand.p4());
  getVar("min(LepGood_jetPtRatiov2,1.5)") = std::min(cand.pt()/jetLepAwareJEC(jet, cand.p4()).Pt(), 1.5);
  getVar("max(LepGood_jetBTagCSV,0)"    ) = std::max(jet ? jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99., 0.);
  getVar("LepGood_sip3d"                ) = cand.dB(pat::Electron::PV3D)/cand.edB(pat::Electron::PV3D); // TODO check definition
  getVar("log(abs(LepGood_dxy))"        ) = pv && elValid ? std::log(std::abs(cand.gsfTrack()->dxy(pv->position()))) : 0.; // TODO check definition
  getVar("log(abs(LepGood_dz))"         ) = pv && elValid ? std::log(std::abs(cand.gsfTrack()->dz (pv->position()))) : 0.;
  getVar("LepGood_mvaIdSpring15"        ) = elValid ? cand.userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") : -1.;

  TTWAnalysis::Dict ret{};
  ret.add("LeptonMVA"     , elValid ? m_tmvaReader.EvaluateRegression("ttH-el-BDTG")[0] : -1.);
  // some more variables that are only calculated here
  ret.add("jetPtRelv2"    , getVar("LepGood_jetPtRelv2"));
  ret.add("jetPtRatio"    , getVar("min(LepGood_jetPtRatiov2,1.5)"));
  ret.add("jetNDaugMVASel", getVar("LepGood_jetNDauChargedMVASel"));
  ret.add("jetBTagCSV"    , getVar("max(LepGood_jetBTagCSV,0)"));
  ret.add("MVAIdNonTrig25nsSpring15", getVar("LepGood_mvaIdSpring15"));

  return ret;
}

TTWAnalysis::Dict TTWAnalysis::DictTTHMuonMVA::evaluate(const pat::Muon& cand,
    const edm::Event* event, const edm::EventSetup* /**/,
    const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  using heppy::IsolationComputer;

  const double rho = event ? getRho(event) : 0.;
  bool muValid{cand.originalObjectRef().isNonnull()};
  const reco::Vertex* pv = event ? getPV(event) : nullptr;

  const double outerR{10.0/std::min(std::max(cand.pt(), 50.),200.)};

  const double absIsoCharged{m_isoComp.chargedAbsIso(cand, outerR, 0.0001, 0.)};
  const double isoNeutralRaw{m_isoComp.neutralAbsIsoRaw(cand, outerR, 0.01, 0.5)};
  const double isoNeutral{std::max(0., isoNeutralRaw - rho * m_ea.getEffectiveArea(cand.eta()) * std::pow(outerR/0.3, 2))};

  edm::Handle<std::vector<pat::Jet>> jets;
  const pat::Jet* jet{nullptr};
  if ( event ) {
    event->getByToken(m_jetsToken, jets);
    jet = findMatchingJet( cand, *jets, m_maxDR );
  }
  int jetNDauChargedMVASel = ( jet && pv ) ? numberOfChargedDaughtersMVASel(jet, cand, pv) : 0;

  getVar("LepGood_pt"                   ) = cand.pt();
  getVar("LepGood_eta"                  ) = cand.eta();
  getVar("LepGood_jetNDauChargedMVASel" ) = jetNDauChargedMVASel;
  getVar("LepGood_miniRelIsoCharged"    ) = absIsoCharged/cand.pt();
  getVar("LepGood_miniRelIsoNeutral"    ) = isoNeutral/cand.pt();
  getVar("LepGood_jetPtRelv2"           ) = ptRelv2(jet, cand.p4());
  getVar("min(LepGood_jetPtRatiov2,1.5)") = std::min(cand.pt()/jetLepAwareJEC(jet, cand.p4()).Pt(), 1.5);
  getVar("max(LepGood_jetBTagCSV,0)"    ) = std::max(jet ? jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99., 0.);
  getVar("LepGood_sip3d"                ) = cand.dB(pat::Muon::PV3D)/cand.edB(pat::Muon::PV3D); // TODO check definition
  getVar("log(abs(LepGood_dxy))"        ) = pv && muValid && cand.muonBestTrack().isNonnull() ? std::log(std::abs(cand.muonBestTrack()->dxy(pv->position()))) : 0.; // TODO check definition
  getVar("log(abs(LepGood_dz))"         ) = pv && muValid && cand.muonBestTrack().isNonnull() ? std::log(std::abs(cand.muonBestTrack()->dz (pv->position()))) : 0.;
  getVar("LepGood_segmentCompatibility" ) = cand.segmentCompatibility();

  TTWAnalysis::Dict ret{};
  ret.add("LeptonMVA"     , muValid ? m_tmvaReader.EvaluateRegression("ttH-mu-BDTG")[0] : -1.);
  // some more variables that are only calculated here
  ret.add("jetPtRelv2"    , getVar("LepGood_jetPtRelv2"));
  ret.add("jetPtRatio"    , getVar("min(LepGood_jetPtRatiov2,1.5)"));
  ret.add("jetNDaugMVASel", getVar("LepGood_jetNDauChargedMVASel"));
  ret.add("jetBTagCSV"    , getVar("max(LepGood_jetBTagCSV,0)"));

  return ret;
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictTTHElectronMVA, "ttw_electronMVAttH");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictTTHMuonMVA, "ttw_muonMVAttH");
