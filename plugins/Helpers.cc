#include "DataFormats/Math/interface/deltaR.h"

#include "Helpers.h"

namespace {

struct ByEta {
  bool operator()(const pat::PackedCandidate *c1, const pat::PackedCandidate *c2) const {
    return c1->eta() < c2->eta();
  }
  bool operator()(float c1eta, const pat::PackedCandidate *c2) const {
    return c1eta < c2->eta();
  }
  bool operator()(const pat::PackedCandidate *c1, float c2eta) const {
    return c1->eta() < c2eta;
  }
};

/*
 * Small helper iterator class (to loop over source candidate pointers of a reco::Candidate
 */
class SourceCandidatePtrIterator
  : public boost::iterator_facade<SourceCandidatePtrIterator,
      reco::CandidatePtr, // value
      boost::random_access_traversal_tag,
      reco::CandidatePtr  // ref
    >
{
public:
  SourceCandidatePtrIterator( const reco::Candidate& base, std::size_t index )
    : m_base(base), m_idx(index)
  {}

  reco::CandidatePtr dereference() const { return m_base.sourceCandidatePtr(m_idx); }
  bool equal( const SourceCandidatePtrIterator& other ) const { return ( &m_base == &(other.m_base) ) && ( m_idx == other.m_idx ); }
  void increment() { ++m_idx; }
  void decrement() { --m_idx; }
  void advance(std::ptrdiff_t d) { m_idx += d; }
  std::ptrdiff_t distance_to(const SourceCandidatePtrIterator& other) const { return m_idx-other.m_idx; }
private:
  const reco::Candidate& m_base;
  std::size_t m_idx;
};
SourceCandidatePtrIterator begin_SourceCandidatePtrs( const reco::Candidate& base )
{ return SourceCandidatePtrIterator(base, 0); }
SourceCandidatePtrIterator end_SourceCandidatePtrs  ( const reco::Candidate& base )
{ return SourceCandidatePtrIterator(base, base.numberOfSourceCandidatePtrs()); }

}

void heppy::IsolationComputer::setPackedCandidates(const std::vector<pat::PackedCandidate> & all, int fromPV_thresh, float dz_thresh, float dxy_thresh, bool also_leptons) 
{
  allcands_ = &all;
  charged_.clear(); neutral_.clear(); pileup_.clear();

  for (const pat::PackedCandidate &p : all) {
    if (p.charge() == 0) {
      neutral_.push_back(&p);
    } else {
      if ( (abs(p.pdgId()) == 211 ) || ( also_leptons && ((abs(p.pdgId()) == 11 ) || (abs(p.pdgId()) == 13 )) ) )  {
        if (p.fromPV() > fromPV_thresh && fabs(p.dz()) < dz_thresh && fabs(p.dxy()) < dxy_thresh ) {
          charged_.push_back(&p);
        } else {
          pileup_.push_back(&p);
        }
      }
    }
  }
  if (weightCone_ > 0) weights_.resize(neutral_.size());
  std::fill(weights_.begin(), weights_.end(), -1.f);
  std::sort(charged_.begin(), charged_.end(), ByEta());
  std::sort(neutral_.begin(), neutral_.end(), ByEta());
  std::sort(pileup_.begin(),  pileup_.end(),  ByEta());
  clearVetos();
}

/// veto footprint from this candidate
void heppy::IsolationComputer::addVetos(const reco::Candidate &cand)
{
  for ( auto it = begin_SourceCandidatePtrs(cand); end_SourceCandidatePtrs(cand) != it; ++it ) {
    if ( it->isNonnull() && it->isAvailable() ) {
      vetos_.push_back(it->get());
    }
  }
}

float heppy::IsolationComputer::isoSumRaw(const std::vector<const pat::PackedCandidate *> & cands, const reco::Candidate &cand, float dR, float innerR, float threshold, SelfVetoPolicy selfVeto, int pdgId) const
{
  float dR2 = dR*dR, innerR2 = innerR*innerR;

  std::vector<const reco::Candidate *> vetos(vetos_);
  if ( selfVetoNone != selfVeto ) {
    for ( auto it = begin_SourceCandidatePtrs(cand); end_SourceCandidatePtrs(cand) != it; ++it ) {
      if ( it->isNonnull() && it->isAvailable() ) {
        vetos.push_back(it->get());
        if ( selfVetoFirst == selfVeto ) break;
      }
    }
  }

  using IT = std::vector<const pat::PackedCandidate *>::const_iterator;
  const IT candsbegin = std::lower_bound(cands.begin(), cands.end(), cand.eta() - dR, ByEta());
  const IT candsend = std::upper_bound(candsbegin, cands.end(), cand.eta() + dR, ByEta());

  double isosum = 0;
  for ( IT icharged{candsbegin}; candsend != icharged; ++icharged ) {
    // pdgId
    if (pdgId > 0 && abs((*icharged)->pdgId()) != pdgId) continue;
    // threshold
    if (threshold > 0 && (*icharged)->pt() < threshold) continue;
    // cone
    float mydr2 = reco::deltaR2(**icharged, cand);
    if (mydr2 >= dR2 || mydr2 < innerR2) continue;
    // veto
    if ( vetos.end() != std::find(vetos.begin(), vetos.end(), *icharged) ) {
      continue;
    }
    // add to sum
    isosum += (*icharged)->pt();
  }
  return isosum;
}

float heppy::IsolationComputer::isoSumNeutralsWeighted(const reco::Candidate &cand, float dR, float innerR, float threshold, SelfVetoPolicy selfVeto, int pdgId) const
{
  if (weightCone_ <= 0) throw cms::Exception("LogicError", "you must set a valid weight cone to use this method");
  float dR2 = dR*dR, innerR2 = innerR*innerR, weightCone2 = weightCone_*weightCone_;

  std::vector<const reco::Candidate *> vetos(vetos_);
  if ( selfVetoNone != selfVeto ) {
    for ( auto it = begin_SourceCandidatePtrs(cand); end_SourceCandidatePtrs(cand) != it; ++it ) {
      if ( it->isNonnull() && it->isAvailable() ) {
        vetos.push_back(it->get());
        if ( selfVetoFirst == selfVeto ) break;
      }
    }
  }

  using IT = std::vector<const pat::PackedCandidate *>::const_iterator;
  const IT charged_begin = std::lower_bound(charged_.begin(), charged_.end(), cand.eta() - dR - weightCone_, ByEta());
  const IT charged_end   = std::upper_bound(charged_begin, charged_.end(), cand.eta() + dR + weightCone_, ByEta());
  const IT pileup_begin  = std::lower_bound(pileup_.begin(), pileup_.end(), cand.eta() - dR - weightCone_, ByEta());
  const IT pileup_end    = std::upper_bound(pileup_begin, pileup_.end(), cand.eta() + dR + weightCone_, ByEta());
  const IT neutral_begin = std::lower_bound(neutral_.begin(), neutral_.end(), cand.eta() - dR, ByEta());
  const IT neutral_end   = std::upper_bound(neutral_begin, neutral_.end(), cand.eta() + dR, ByEta());

  double isosum = 0.0;
  for ( IT ineutral{neutral_begin}; neutral_end != ineutral; ++ineutral ) {
    // pdgId
    if (pdgId > 0 && abs((*ineutral)->pdgId()) != pdgId) continue;
    // threshold
    if (threshold > 0 && (*ineutral)->pt() < threshold) continue;
    // cone
    float mydr2 = reco::deltaR2(**ineutral, cand);
    if (mydr2 >= dR2 || mydr2 < innerR2) continue;
    // veto
    if ( vetos.end() != std::find(vetos.begin(), vetos.end(), *ineutral) ) {
      continue;
    }
    // weight
    float &w = weights_[ineutral-neutral_.begin()];
    if (w == -1.f) {
      double sumc = 0, sump = 0.0;
      for ( IT icharged{charged_begin}; charged_end != icharged; ++icharged ) {
        float hisdr2 = std::max<float>(reco::deltaR2(**icharged, **ineutral), 0.01f);
        if (hisdr2 > weightCone2) continue;
        if (std::find(vetos_.begin(), vetos_.end(), *icharged) != vetos_.end()) {
          continue;
        }
        sumc += std::log( (*icharged)->pt() / std::sqrt(hisdr2) );
      }
      for ( IT ipileup = pileup_begin; pileup_end != ipileup; ++ipileup ) {
        float hisdr2 = std::max<float>(reco::deltaR2(**ipileup, **ineutral), 0.01f);
        if (hisdr2 > weightCone2) continue;
        if (std::find(vetos_.begin(), vetos_.end(), *ipileup) != vetos_.end()) {
          continue;
        }
        sumc += std::log( (*ipileup)->pt() / std::sqrt(hisdr2) );
      }
      w = (sump == 0 ? 1 : sumc/(sump+sumc));
    }
    // add to sum
    isosum += w * (*ineutral)->pt();
  }
  return isosum;
}
