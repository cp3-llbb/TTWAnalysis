#ifndef TTWANALYSIS_NEWTYPES_H
#define TTWANALYSIS_NEWTYPES_H

#include <boost/variant.hpp>
// #include <boost/range/any_range.hpp>

#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <Math/LorentzVector.h>
#include "cp3_llbb/TTWAnalysis/interface/Types.h"

namespace TTWAnalysis {
  /*
  // TODO move to the boost::any_range version when boost version > 1.57
  template<typename ObjType>
  using Range = boost::any_range<ObjType,boost::random_access_traversal_tag,const ObjType&,std::ptrdiff_t>;
  */

  using ElectronCut = std::function<bool(const pat::Electron&)>;
  using MuonCut     = std::function<bool(const pat::Muon&)>;
  using JetCut      = std::function<bool(const pat::Jet&)>;

  using index_t = std::uint16_t;
  using sindex_t = std::int16_t;
  using indexlist_t = std::vector<index_t>;

  struct Lepton {
    bool m_elOrMu;
    index_t idx = -1;
    boost::variant<edm::Ptr<pat::Electron>,edm::Ptr<pat::Muon>> ptr;
    sindex_t hlt_idx = -1;
    Lepton() {}
    Lepton( index_t idx_, edm::Ptr<pat::Electron> ele ) : m_elOrMu(true) , idx(idx_), ptr(ele) {}
    Lepton( index_t idx_, edm::Ptr<pat::Muon>     mu  ) : m_elOrMu(false), idx(idx_), ptr(mu ) {}
    bool isEl() const { return m_elOrMu; }
    bool isMu() const { return !m_elOrMu; }
    // WARNING may throw
    edm::Ptr<pat::Electron> electron() const { return boost::get<edm::Ptr<pat::Electron>>(ptr); }
    edm::Ptr<pat::Muon    > muon    () const { return boost::get<edm::Ptr<pat::Muon    >>(ptr); }

    const reco::RecoCandidate* recoCand() const { if ( m_elOrMu ) return electron().get(); else return muon().get(); }
    const reco::RecoCandidate::LorentzVector& p4 () const { return recoCand()->p4(); }
    int pdgId() const { return recoCand()->pdgId(); }
    int charge() const { return recoCand()->charge(); }
  };

  class LeptonCut {
    public:
      LeptonCut() {}
      LeptonCut( ElectronCut elCut, MuonCut muCut ) : m_elCut(elCut), m_muCut(muCut) {}

      bool operator() ( const Lepton& l ) const { return l.isEl() ? m_elCut(*(l.electron().get())) : m_muCut(*(l.muon().get())); }
    private:
      ElectronCut m_elCut;
      MuonCut m_muCut;
  };

  // another helper : a pair of leptons, with their indices in the selected leptons container
  struct DiLepton {
    std::pair<index_t,index_t> lidxs;
    const Lepton* first  = nullptr;
    const Lepton* second = nullptr;
    DiLepton() {}
    DiLepton( index_t lidx1, const Lepton& first_, index_t lidx2, const Lepton& second_ )
      : lidxs(std::pair<index_t,index_t>(lidx1,lidx2)), first(&first_), second(&second_) {}
    // TODO maybe don't need all of these getters
    reco::RecoCandidate::LorentzVector p4 () const { return first->p4()+second->p4(); }
    std::pair<index_t,index_t> idxs () const { return std::pair<index_t,index_t>(first->idx,second->idx); }
    bool isElEl() const { return first->isEl() && second->isEl(); }
    bool isElMu() const { return first->isEl() && second->isMu(); }
    bool isMuEl() const { return first->isMu() && second->isEl(); }
    bool isMuMu() const { return first->isMu() && second->isMu(); }
    bool isOS  () const { return first->recoCand()->pdgId()*second->recoCand()->pdgId() < 0.; }// TODO check (above)
    bool isSF  () const { return first->m_elOrMu == second->m_elOrMu; }
    /*
    m_diLepton.DR = VectorUtil::DeltaR(l1.p4, l2.p4);
    m_diLepton.DEta = TTWAnalysis::DeltaEta(l1.p4, l2.p4);
    m_diLepton.DPhi = VectorUtil::DeltaPhi(l1.p4, l2.p4);
    */
  };

  class DiLeptonCut {
    public:
      DiLeptonCut() {}
      DiLeptonCut( LeptonCut cut1, LeptonCut cut2 ) : m_cut1(cut1), m_cut2(cut2) {}

      bool operator() ( const DiLepton& ll ) const { return m_cut1(*(ll.first)) && m_cut2(*(ll.second)); }
    private:
      LeptonCut m_cut1;
      LeptonCut m_cut2;
  };

  struct DiJet {
    std::pair<index_t,index_t> jidxs;
    edm::Ptr<pat::Jet> first;
    edm::Ptr<pat::Jet> second;
    DiJet() {}
    DiJet( index_t jidx1, edm::Ptr<pat::Jet> first_, index_t jidx2, edm::Ptr<pat::Jet> second_ )
      : jidxs(std::pair<index_t,index_t>(jidx1,jidx2)), first(first_), second(second_) {}
  };

  class DiJetCut {
    public:
      DiJetCut() {}
      DiJetCut( JetCut cut1, JetCut cut2 ) : m_cut1(cut1), m_cut2(cut2) {}

      bool operator() ( const DiJet& jj ) const { return m_cut1(*(jj.first.get())) && m_cut2(*(jj.second.get())); }
    private:
      JetCut m_cut1;
      JetCut m_cut2;
  };

  struct DiLeptonDiJet {
    index_t llIdx;
    index_t jjIdx;
    const DiLepton* ll = nullptr;
    const DiJet* jj    = nullptr;
    DiLeptonDiJet() {}
    DiLeptonDiJet( index_t llIdx_, const DiLepton& ll_, index_t jjIdx_, const DiJet& jj_ )
      : llIdx(llIdx_), jjIdx(jjIdx_), ll(&ll_), jj(&jj_) {}
  };

  class DiLeptonDiJetCut {
    public:
      DiLeptonDiJetCut() {}
      DiLeptonDiJetCut( DiLeptonCut llCut, DiJetCut jjCut ) : m_llCut(llCut), m_jjCut(jjCut) {}
      bool operator() ( const DiLeptonDiJet& lljj ) const { return m_llCut(*(lljj.ll)) && m_jjCut(*(lljj.jj)); }
    private:
      DiLeptonCut m_llCut;
      DiJetCut m_jjCut;
  };

  struct DiLeptonDiJetMet {
    index_t lljjIdx;
    myLorentzVector met;
    const DiLeptonDiJet* lljj = nullptr;
    DiLeptonDiJetMet() {}
    DiLeptonDiJetMet( index_t lljjIdx_, const DiLeptonDiJet& lljj_, const myLorentzVector& met_ )
      : lljjIdx(lljjIdx_), met(met_), lljj(&lljj_) {}
  };
}

#endif // TTWANALYSIS_NEWTYPES_H
