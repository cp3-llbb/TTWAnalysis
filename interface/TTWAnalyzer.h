#pragma once

#include <string>
#include <utility>
#include <vector>
#include <limits>
#include <boost/container/flat_map.hpp>

#include <cp3_llbb/Framework/interface/Analyzer.h>
class HLTProducer;

#include <cp3_llbb/TTWAnalysis/interface/NewTypes.h>
#include <cp3_llbb/TTWAnalysis/interface/Tools.h>

#include "cp3_llbb/TTWAnalysis/interface/CandidatesProducer.h"

template<class Map>
class IdxListForWPsHolder {
  public:
    using index_t = TTWAnalysis::index_t;
    using indexlist_t = TTWAnalysis::indexlist_t;

    IdxListForWPsHolder() : m_wpMap(nullptr) {} // TODO get rid of this

    IdxListForWPsHolder( const Map& wpMap, const ROOT::TreeGroup& tree, const std::string& prefix )
      : m_wpMap(&wpMap)
    {
      ROOT::TreeGroup tree_{tree};
      for ( const auto& iWP : *m_wpMap ) {
        std::string brName{prefix+iWP.first};
        indexlist_t& br = tree_[brName].write<indexlist_t>();
        this->m_indices.push_back(&br);
      }
    }

    class const_iterator
      : public boost::iterator_facade<const_iterator,
          std::pair<typename Map::const_reference,const indexlist_t*>,
          boost::random_access_traversal_tag,
          std::pair<typename Map::const_reference,const indexlist_t*>  // ref
        >
    {
    public:
      const_iterator( const IdxListForWPsHolder& parent, index_t idx ) : m_parent(parent), m_idx(idx) {}

      //std::pair<typename Map::const_reference,const indexlist_t*> dereference() const { return std::make_pair(this->m_parent.m_wpMap->nth(this->m_idx), this->m_parent.m_indices.at(this->m_idx)); }
      std::pair<typename Map::const_reference,const indexlist_t*> dereference() const { return std::pair<typename Map::const_reference,const indexlist_t*>(*(this->m_parent.m_wpMap->begin()+this->m_idx), this->m_parent.m_indices.at(this->m_idx)); }
      bool equal( const const_iterator& other ) const { return (&(this->m_parent)) == (&other.m_parent) && ( this->m_idx == other.m_idx); }
      void increment() { ++(this->m_idx); }
      void decrement() { --(this->m_idx); }
      void advance(std::ptrdiff_t d) { (this->m_idx) += d; }
      std::ptrdiff_t distance_to(const const_iterator& other) const { return this->m_idx-other.m_idx; }

      // const typename Map::key_type& name() const { return this->m_parent.m_wpMap->nth(this->m_idx).first ; }
      // const typename Map::mapped_type& cut () const { return this->m_parent.m_wpMap->nth(this->m_idx).second; }
      const typename Map::key_type& name() const { return (this->m_parent.m_wpMap->begin()+this->m_idx)->first ; }
      const typename Map::mapped_type& cut () const { return (this->m_parent.m_wpMap->begin()+this->m_idx)->second; }
      indexlist_t* idxList() { return this->m_parent.m_indices.at(this->m_idx); }
      index_t index() const { return this->m_idx; }

    private:
      const IdxListForWPsHolder& m_parent;
      index_t m_idx;
    };

    index_t size() const { return this->m_indices.size(); }
    const_iterator cbegin() const { return const_iterator(*this, 0); }
    const_iterator cend  () const { return const_iterator(*this, static_cast<index_t>(this->m_indices.size())); }
    const_iterator at( index_t idx ) const { return const_iterator(*this, idx); }
    //const_iterator find( const typename Map::key_type& ky ) const { return const_iterator(*this, this->m_wpMap->index_of(ky)); }
    const_iterator find( const typename Map::key_type& ky ) const { return const_iterator(*this, std::distance(this->m_wpMap->begin(), this->m_wpMap->find(ky))); }

    class iterator
      : public boost::iterator_facade<iterator,
          std::pair<typename Map::const_reference,indexlist_t*>,
          boost::random_access_traversal_tag,
          std::pair<typename Map::const_reference,indexlist_t*> // ref
        >
    {
    public:
      iterator( IdxListForWPsHolder& parent, index_t idx ) : m_parent(parent), m_idx(idx) {}

      //std::pair<typename Map::const_reference,indexlist_t*> dereference() { return std::make_pair(this->m_parent.m_wpMap->nth(this->m_idx), this->m_parent.m_indices.at(this->m_idx)); }
      std::pair<typename Map::const_reference,indexlist_t*> dereference() { return std::pair<typename Map::const_reference,indexlist_t*>(*(this->m_parent.m_wpMap->begin()+this->m_idx), this->m_parent.m_indices.at(this->m_idx)); }
      bool equal( const iterator& other ) const { return (&(this->m_parent)) == (&(other.m_parent)) && ( this->m_idx == other.m_idx); }
      void increment() { ++(this->m_idx); }
      void decrement() { --(this->m_idx); }
      void advance(std::ptrdiff_t d) { (this->m_idx) += d; }
      std::ptrdiff_t distance_to(const iterator& other) const { return this->m_idx-other.m_idx; }

      // const typename Map::key_type& name() const { return this->m_parent.m_wpMap->nth(this->m_idx).first ; }
      // const typename Map::mapped_type& cut () const { return this->m_parent.m_wpMap->nth(this->m_idx)->second; }
      const typename Map::key_type& name() const { return (this->m_parent.m_wpMap->begin()+this->m_idx).first ; }
      const typename Map::mapped_type& cut () const { return (this->m_parent.m_wpMap->begin()+this->m_idx)->second; }
      indexlist_t* idxList() { return this->m_parent.m_indices.at(this->m_idx); }
      index_t index() const { return this->m_idx; }

    private:
      IdxListForWPsHolder& m_parent;
      index_t m_idx;
    };
    iterator begin() { return iterator{*this, 0}; }
    iterator end() { return iterator{*this, static_cast<index_t>(this->m_indices.size())}; }
    iterator at( index_t idx ) { return iterator(*this, idx); }

  private:
    const Map* m_wpMap;
    std::vector<indexlist_t*> m_indices;
};

class TTWAnalyzer: public Framework::Analyzer {
    public:
        TTWAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config);

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const AnalyzersManager&, const CategoryManager&) override;
        virtual void doConsumes(const edm::ParameterSet&, edm::ConsumesCollector&& collector) override;

        virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet&) override;

        // TODO move to the boost::any_range version when boost version > 1.57
        // access to the containers
        // template<typename ObjType> TTWAnalysis::Range<ObjType> getList() const;
        template<typename patObjectType> const edm::PtrVector<patObjectType>& getPtrList() const;
        template<typename ttwObjectType> const std::vector   <ttwObjectType>& getObjList() const;

        // access to the indices of selected objects (by working point name)
        const TTWAnalysis::indexlist_t& selectedElectrons  ( const std::string& wpName ) const { return *(m_idxElWP.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedMuons      ( const std::string& wpName ) const { return *(m_idxMuWP.find(wpName)->second); }

        const TTWAnalysis::indexlist_t& selectedLeptons    ( const std::string& wpName ) const { return *(m_idxlWP.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedDiLeptons  ( const std::string& wpName ) const { return *(m_idxllWP.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedJets       ( const std::string& wpName ) const { return *(m_idxjDRWP.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedBJets      ( const std::string& wpName ) const { return *(m_idxbDRWP_PT.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedBJets_tag  ( const std::string& wpName ) const { return *(m_idxbDRWP_tag.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedDiJets     ( const std::string& wpName ) const { return *(m_idxjjDRWP.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedDiBJets    ( const std::string& wpName ) const { return *(m_idxbbDRWP_PT.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedDiBJets_tag( const std::string& wpName ) const { return *(m_idxbbDRWP_tag.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedDiLeptonDiJets     ( const std::string& wpName ) const { return *(m_idxlljjDRWP.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedDiLeptonDiBJets    ( const std::string& wpName ) const { return *(m_idxllbbDRWP_PT.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedDiLeptonDiBJets_tag( const std::string& wpName ) const { return *(m_idxllbbDRWP_tag.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedDiLeptonDiJetsMET     ( const std::string& wpName ) const { return *(m_idxlljjmDRWP.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedDiLeptonDiBJetsMET    ( const std::string& wpName ) const { return *(m_idxllbbmDRWP_PT.find(wpName)->second); }
        const TTWAnalysis::indexlist_t& selectedDiLeptonDiBJetsMET_tag( const std::string& wpName ) const { return *(m_idxllbbmDRWP_tag.find(wpName)->second); }

    private:

        // Producers name
        const std::string m_electrons_producer;
        const std::string m_muons_producer;
        const std::string m_jets_producer;
        const std::string m_met_producer;

        edm::EDGetTokenT<std::vector<reco::Vertex>> m_vertices_token;

        boost::container::flat_map<std::string,TTWAnalysis::ElectronCut> m_elWP;
        boost::container::flat_map<std::string,TTWAnalysis::MuonCut    > m_muWP;
        boost::container::flat_map<std::string,TTWAnalysis::LeptonCut  > m_lWP ;
        boost::container::flat_map<std::string,TTWAnalysis::DiLeptonCut> m_llWP;
        boost::container::flat_map<std::string,std::pair<std::string,TTWAnalysis::JetCut     >> m_bWP    ; // with the name of a j WP
        boost::container::flat_map<std::string,std::pair<std::string,TTWAnalysis::DiJetCut   >> m_bbWP   ; // with the name of a jj WP
        boost::container::flat_map<std::string,std::pair<std::string,TTWAnalysis::DiLeptonCut>> m_lljjWP ; // with the name of a lepton WP
        boost::container::flat_map<std::string,std::pair<std::string,TTWAnalysis::DiJetCut   >> m_llbbWP ; // with the name of a lljj WP
        boost::container::flat_map<std::string,std::pair<std::string,TTWAnalysis::DiJetCut   >> m_llbbmWP; // with the name of a lljj WP

        IdxListForWPsHolder<decltype(m_elWP  )> m_idxElWP;
        IdxListForWPsHolder<decltype(m_muWP  )> m_idxMuWP;
        IdxListForWPsHolder<decltype(m_lWP   )> m_idxlWP;
        IdxListForWPsHolder<decltype(m_llWP  )> m_idxllWP;
    public:
        decltype(m_idxllWP)::const_iterator begin_diLeptonSelections() const { return m_idxllWP.cbegin(); }
        decltype(m_idxllWP)::const_iterator end_diLeptonSelections  () const { return m_idxllWP.cend(); }
    private:
        IdxListForWPsHolder<decltype(m_lWP   )> m_idxjDRWP;
        IdxListForWPsHolder<decltype(m_bWP   )> m_idxbDRWP_PT;
        IdxListForWPsHolder<decltype(m_bWP   )> m_idxbDRWP_tag;
        IdxListForWPsHolder<decltype(m_lWP   )> m_idxjjDRWP;
        IdxListForWPsHolder<decltype(m_bbWP  )> m_idxbbDRWP_PT;
        IdxListForWPsHolder<decltype(m_bbWP  )> m_idxbbDRWP_tag;
        IdxListForWPsHolder<decltype(m_lljjWP)> m_idxlljjDRWP;
        IdxListForWPsHolder<decltype(m_llbbWP)> m_idxllbbDRWP_PT;
        IdxListForWPsHolder<decltype(m_llbbWP)> m_idxllbbDRWP_tag;
        IdxListForWPsHolder<decltype(m_lljjWP)> m_idxlljjmDRWP;
        IdxListForWPsHolder<decltype(m_llbbWP)> m_idxllbbmDRWP_PT;
        IdxListForWPsHolder<decltype(m_llbbWP)> m_idxllbbmDRWP_tag;

        const float m_electronPtCut, m_electronEtaCut;

        const float m_muonPtCut, m_muonEtaCut;

        const float m_jetPtCut, m_jetEtaCut, m_jetPUID, m_jetDRleptonCut;
        const std::string m_jetBTagName;
        TTWAnalysis::JetCut m_jetIDCut;

        const float m_hltDRCut, m_hltDPtCut;
        TTWAnalysis::sindex_t matchOfflineLepton( const TTWAnalysis::Lepton& lepton, const HLTProducer& hlt ) const;

        edm::PtrVector<pat::Electron>              m_electrons;
        edm::PtrVector<pat::Muon>                  m_muons;
        std::vector<TTWAnalysis::Lepton>           m_leptons;
        std::vector<TTWAnalysis::DiLepton>         m_dileptons;
        edm::PtrVector<pat::Jet>                   m_jets;
        std::vector<TTWAnalysis::DiJet>            m_dijets;
        std::vector<TTWAnalysis::DiLeptonDiJet>    m_dileptondijets;
        std::vector<TTWAnalysis::DiLeptonDiJetMet> m_dileptondijetmets;
};

// TODO move to the boost::any_range version when boost version > 1.57
//
//
template<> inline const edm::PtrVector<pat::Electron>& TTWAnalyzer::getPtrList<pat::Electron>() const { return m_electrons; }
template<> inline const edm::PtrVector<pat::Muon    >& TTWAnalyzer::getPtrList<pat::Muon    >() const { return m_muons    ; }
template<> inline const edm::PtrVector<pat::Jet     >& TTWAnalyzer::getPtrList<pat::Jet     >() const { return m_jets     ; }
template<> inline const std::vector<TTWAnalysis::Lepton          >& TTWAnalyzer::getObjList<TTWAnalysis::Lepton          >() const { return m_leptons          ; }
template<> inline const std::vector<TTWAnalysis::DiLepton        >& TTWAnalyzer::getObjList<TTWAnalysis::DiLepton        >() const { return m_dileptons        ; }
template<> inline const std::vector<TTWAnalysis::DiJet           >& TTWAnalyzer::getObjList<TTWAnalysis::DiJet           >() const { return m_dijets           ; }
template<> inline const std::vector<TTWAnalysis::DiLeptonDiJet   >& TTWAnalyzer::getObjList<TTWAnalysis::DiLeptonDiJet   >() const { return m_dileptondijets   ; }
template<> inline const std::vector<TTWAnalysis::DiLeptonDiJetMet>& TTWAnalyzer::getObjList<TTWAnalysis::DiLeptonDiJetMet>() const { return m_dileptondijetmets; }

/*
// #include "cp3_llbb/TTWAnalysis/interface/PtrVectorItr.h"

template<> TTWAnalysis::Range<edm::Ptr<pat::Electron>> TTWAnalyzer::getList<edm::Ptr<pat::Electron>>() const { return TTWAnalysis::Range<edm::Ptr<pat::Electron>>(TTWAnalysis::begin(m_electrons), TTWAnalysis::end(m_electrons)); }
template<> TTWAnalysis::Range<edm::Ptr<pat::Muon    >> TTWAnalyzer::getList<edm::Ptr<pat::Muon    >>() const { return TTWAnalysis::Range<edm::Ptr<pat::Muon    >>(TTWAnalysis::begin(m_muons    ), TTWAnalysis::end(m_muons)    ); }
template<> TTWAnalysis::Range<edm::Ptr<pat::Jet     >> TTWAnalyzer::getList<edm::Ptr<pat::Jet     >>() const { return TTWAnalysis::Range<edm::Ptr<pat::Jet     >>(TTWAnalysis::begin(m_jets     ), TTWAnalysis::end(m_jets)     ); }
template<> TTWAnalysis::Range<TTWAnalysis::Lepton          > TTWAnalyzer::getList<TTWAnalysis::Lepton          >() const { return TTWAnalysis::Range<TTWAnalysis::Lepton          >(m_leptons          ); }
template<> TTWAnalysis::Range<TTWAnalysis::DiLepton        > TTWAnalyzer::getList<TTWAnalysis::DiLepton        >() const { return TTWAnalysis::Range<TTWAnalysis::DiLepton        >(m_dileptons        ); }
template<> TTWAnalysis::Range<TTWAnalysis::DiJet           > TTWAnalyzer::getList<TTWAnalysis::DiJet           >() const { return TTWAnalysis::Range<TTWAnalysis::DiJet           >(m_dijets           ); }
template<> TTWAnalysis::Range<TTWAnalysis::DiLeptonDiJet   > TTWAnalyzer::getList<TTWAnalysis::DiLeptonDiJet   >() const { return TTWAnalysis::Range<TTWAnalysis::DiLeptonDiJet   >(m_dileptondijets   ); }
template<> TTWAnalysis::Range<TTWAnalysis::DiLeptonDiJetMet> TTWAnalyzer::getList<TTWAnalysis::DiLeptonDiJetMet>() const { return TTWAnalysis::Range<TTWAnalysis::DiLeptonDiJetMet>(m_dileptondijetmets); }
*/
