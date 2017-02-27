#ifndef TTWANALYSIS_PLUGINS_HELPERS_H
#define TTWANALYSIS_PLUGINS_HELPERS_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

namespace TTWAnalysis {
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
 * Small helper to access the primary vertex
 */
class DictPVHelper {
public:
  DictPVHelper(const edm::ParameterSet& config) {}
  void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector)
  {
    m_vertices_token = collector.consumes<std::vector<reco::Vertex>>(config.getUntrackedParameter<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices")));
  }
protected:
  const reco::Vertex* getPV(const edm::Event* event) const
  {
    edm::Handle<std::vector<reco::Vertex>> handle;
    event->getByToken(m_vertices_token, handle);
    if ( ! handle->empty() ) {
      return &(handle->at(0));
    }
    return nullptr;
  }
private:
  edm::EDGetTokenT<std::vector<reco::Vertex>> m_vertices_token;
};

}


namespace heppy {
/* heppy::IsolationComputer
 * Copied and slightly modified version of
 * - https://github.com/CERN-PH-CMG/cmg-cmssw/blob/ebe350c7d742c6b29b55e60f5bf872bfc3c5afda/PhysicsTools/Heppy/interface/IsolationComputer.h
 * - https://github.com/CERN-PH-CMG/cmg-cmssw/blob/ebe350c7d742c6b29b55e60f5bf872bfc3c5afda/PhysicsTools/Heppy/src/IsolationComputer.cc
 */
class IsolationComputer {
public:
  /// Create the calculator; optionally specify a cone for computing deltaBeta weights
  IsolationComputer(float weightCone=-1) : weightCone_(weightCone) {}

  /// Self-veto policy
  enum SelfVetoPolicy {
      selfVetoNone=0, selfVetoAll=1, selfVetoFirst=2
  };
  /// Initialize with the list of packed candidates (note: clears also all vetos)
  void setPackedCandidates(const std::vector<pat::PackedCandidate> & all, int fromPV_thresh=1, float dz_thresh=9999., float dxy_thresh=9999., bool also_leptons=false) ;


  /// veto footprint from this candidate, for the isolation of all candidates and also for calculation of neutral weights (if used)
  void addVetos(const reco::Candidate &cand) ;

  /// clear all vetos
  void clearVetos() { vetos_.clear(); }

  /// Isolation from charged from the PV
  float chargedAbsIso(const reco::Candidate &cand, float dR, float innerR=0, float threshold=0, SelfVetoPolicy selfVeto=selfVetoAll) const
  { return isoSumRaw(charged_, cand, dR, innerR, threshold, selfVeto); }

  /// Isolation from charged from PU
  float puAbsIso(const reco::Candidate &cand, float dR, float innerR=0, float threshold=0, SelfVetoPolicy selfVeto=selfVetoAll) const
  { return isoSumRaw(pileup_, cand, dR, innerR, threshold, selfVeto); }

  /// Isolation from all neutrals (uncorrected)
  float neutralAbsIsoRaw(const reco::Candidate &cand, float dR, float innerR=0, float threshold=0, SelfVetoPolicy selfVeto=selfVetoAll) const
  { return isoSumRaw(neutral_, cand, dR, innerR, threshold, selfVeto); }

  /// Isolation from neutral hadrons (uncorrected)
  float neutralHadAbsIsoRaw(const reco::Candidate &cand, float dR, float innerR=0, float threshold=0, SelfVetoPolicy selfVeto=selfVetoAll) const
  { return isoSumRaw(neutral_, cand, dR, innerR, threshold, selfVeto, 130); }

  /// Isolation from photons (uncorrected)
  float photonAbsIsoRaw(const reco::Candidate &cand, float dR, float innerR=0, float threshold=0, SelfVetoPolicy selfVeto=selfVetoAll) const
  { return isoSumRaw(neutral_, cand, dR, innerR, threshold, selfVeto, 22); }

  /// Isolation from all neutrals (with weights)
  float neutralAbsIsoWeighted(const reco::Candidate &cand, float dR, float innerR=0, float threshold=0, SelfVetoPolicy selfVeto=selfVetoAll) const
  { return isoSumNeutralsWeighted(cand, dR, innerR, threshold, selfVeto); }

  /// Isolation from neutral hadrons (with weights)
  float neutralHadAbsIsoWeighted(const reco::Candidate &cand, float dR, float innerR=0, float threshold=0, SelfVetoPolicy selfVeto=selfVetoAll) const
  { return isoSumNeutralsWeighted(cand, dR, innerR, threshold, selfVeto, 130); }

  /// Isolation from photons (with weights)
  float photonAbsIsoWeighted(const reco::Candidate &cand, float dR, float innerR=0, float threshold=0, SelfVetoPolicy selfVeto=selfVetoAll) const
  { return isoSumNeutralsWeighted(cand, dR, innerR, threshold, selfVeto, 22); }

  void updateEvent( const edm::Event* event );

  void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector)
  {
    m_packedCandidates_token = collector.consumes<std::vector<pat::PackedCandidate>>(config.getParameter<edm::InputTag>("packedCandidates"));
  }

protected:
  edm::EventID m_lastEvent;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> m_packedCandidates_token;
  const std::vector<pat::PackedCandidate> * allcands_;
  float weightCone_;
  // collections of objects, sorted in eta
  std::vector<const pat::PackedCandidate *> charged_, neutral_, pileup_;
  mutable std::vector<float> weights_;
  std::vector<const reco::Candidate *> vetos_;

  float isoSumRaw(const std::vector<const pat::PackedCandidate *> & cands, const reco::Candidate &cand, float dR, float innerR, float threshold, SelfVetoPolicy selfVeto, int pdgId=-1) const ;
  float isoSumNeutralsWeighted(const reco::Candidate &cand, float dR, float innerR, float threshold, SelfVetoPolicy selfVeto, int pdgId=-1) const ;
};

}


#endif // TTWANALYSIS_PLUGINS_HELPERS_H
