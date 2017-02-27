#include "cp3_llbb/Framework/interface/HLTProducer.h"
#include "cp3_llbb/TTWAnalysis/interface/TTWAnalyzer.h"
#include "cp3_llbb/TTWAnalysis/interface/CandidatesProducer.h"

#include "cp3_llbb/TTWAnalysis/interface/TTWDileptonCategory.h"

#include <cassert>
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>


using namespace TTWAnalysis;

GDileptonCategory::GDileptonCategory()
 : Category()
 , m_llIsInCateg("")
{}

GDileptonCategory::~GDileptonCategory() {}

void GDileptonCategory::configure(const edm::ParameterSet& conf)
{
  using NamedDiLeptonCut = std::pair<std::string,DiLeptonCut>;
  // pre-ana in category
  m_nReqEl = conf.getParameter<unsigned int>("NElectrons");
  m_nReqMu = conf.getParameter<unsigned int>("NMuons");
  m_reqCh  = conf.getParameter<int>("Charge");
  // working points
  m_llWPs = conf.getParameter<std::vector<std::string>>("WPs");
  // post-ana in category
  m_llIsInCateg = DiLeptonHybridCut(conf.getParameter<std::string>("Category"));
  // post-ana cuts
  m_postCuts.clear();
  m_postCuts.emplace_back(NamedDiLeptonCut("Category", m_llIsInCateg));
  const auto& triggers = conf.getParameter<std::vector<std::string>>("HLT");
  if ( ! triggers.empty() ) {
    m_hltMatcher.addRegex(triggers);
    m_postCuts.emplace_back(NamedDiLeptonCut("DiLeptonTriggerMatch", std::reference_wrapper<HLTMatch<TTWAnalysis::DiLepton>>(m_hltMatcher)));
  }
  for ( const auto& otherCutConfig : conf.getParameter<edm::VParameterSet>("Cuts") ) {
    assert(otherCutConfig.getParameterNames().size() == 1);
    for ( const std::string& otherCutName : otherCutConfig.getParameterNames() ) {
      m_postCuts.emplace_back(NamedDiLeptonCut(otherCutName, DiLeptonHybridCut(otherCutConfig.getParameter<std::string>(otherCutName))));
    }
  }
}

bool GDileptonCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const
{
  unsigned int nPass;
  if ( m_nReqEl > 0 ) {
    const auto& electrons = producers.get<CandidatesProducer<pat::Electron>>("electrons").selected();
    if ( m_reqCh != 0 ) {
      nPass = std::count_if(std::begin(electrons), std::end(electrons), [this] ( const edm::Ptr<pat::Electron> el ) { return m_reqCh == el->charge(); });
    } else {
      nPass = electrons.size();
    }
    if ( nPass < m_nReqEl ) {
      return false;
    }
  }
  if ( m_nReqMu > 0 ) {
    const auto& muons = producers.get<CandidatesProducer<pat::Muon>>("muons").selected();
    if ( m_reqCh != 0 ) {
      nPass = std::count_if(std::begin(muons), std::end(muons), [this] ( const edm::Ptr<pat::Muon> mu ) { return m_reqCh == mu->charge(); });
    } else {
      nPass = muons.size();
    }
    if ( nPass < m_nReqMu ) {
      return false;
    }
  }
  return true;
}

bool GDileptonCategory::diLeptonIsInCategory( const TTWAnalyzer* ttW, uint16_t dilepIdx ) const
{
  const auto& dilep = ttW->getObjList<TTWAnalysis::DiLepton>()[dilepIdx];
  const auto& leptons = ttW->getObjList<TTWAnalysis::Lepton>();
  if ( m_llIsInCateg(dilep) ) {
    if ( m_reqCh == 0 ) {
      return true;
    } else {
      if ( ( m_reqCh == leptons[dilep.lidxs.first].charge() )
        && ( m_reqCh == leptons[dilep.lidxs.second].charge() ) ) {
        return true;
      }
    }
  }
  return false;
}

bool GDileptonCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const
{
  const TTWAnalyzer& ttW = analyzers.get<TTWAnalyzer>("ttW");
  if ( m_llWPs.empty() ) {
    return true;
  }
  for ( const auto& sel : m_llWPs ) {
    const auto& catDiLeptons = ttW.selectedDiLeptons(sel);
    if ( ( ! catDiLeptons.empty() )
      && ( std::any_of(std::begin(catDiLeptons), std::end(catDiLeptons),
                  std::bind(&GDileptonCategory::diLeptonIsInCategory, this, &ttW, std::placeholders::_1)
                  ) ) ) {
      return true;
    }
  }
  return false;
}

void GDileptonCategory::register_cuts(CutManager& manager)
{
  for ( const auto& sel : m_llWPs ) {
    for ( const auto& cut : m_postCuts ) {
      manager.new_cut(cut.first+"_"+sel, cut.first+"_"+sel);
    }
  }
}

void GDileptonCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {

  const TTWAnalyzer& ttW = analyzers.get<TTWAnalyzer>("ttW");
  const auto& diLeptons = ttW.getObjList<TTWAnalysis::DiLepton>();
  m_hltMatcher.setHLTProd(&producers.get<HLTProducer>("hlt"));

  for ( const auto& sel : m_llWPs ) {
    if ( ! ttW.selectedDiLeptons(sel).empty() ) {
      const auto& sel_before_cat = ttW.selectedDiLeptons(sel);
      std::vector<uint16_t> sel_dileptons;
      sel_dileptons.reserve(sel_before_cat.size());
      std::copy_if(std::begin(sel_before_cat), std::end(sel_before_cat), std::back_inserter(sel_dileptons),
                   std::bind(&GDileptonCategory::diLeptonIsInCategory, this, &ttW, std::placeholders::_1));
      for ( const auto& cut : m_postCuts ) {
        std::vector<uint16_t> passSel_dileptons;
        passSel_dileptons.reserve(sel_dileptons.size());
        std::copy_if(std::begin(sel_dileptons), std::end(sel_dileptons), std::back_inserter(passSel_dileptons),
                    [&diLeptons,&cut] ( uint16_t idx ) -> bool { return cut.second(diLeptons[idx]); }
                );
        if ( ! passSel_dileptons.empty() ) {
          manager.pass_cut(cut.first+"_"+sel);
        }
      }
    }
  }
}
