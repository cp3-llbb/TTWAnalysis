#pragma once

#include <utility>
#include <vector>
#include <string>
#include <functional>

#include "cp3_llbb/Framework/interface/Category.h"
#include "cp3_llbb/TTWAnalysis/interface/Types.h"
#include "cp3_llbb/TTWAnalysis/interface/HLTMatch.h"

// fwd
class HLTProducer;
class TTWAnalyzer;

namespace TTWAnalysis {

class GDileptonCategory : public Category {
  public:
    GDileptonCategory();
    virtual ~GDileptonCategory();
    virtual void configure(const edm::ParameterSet& conf) override;
  private:
    using DiLeptonCut = std::function<bool(const TTWAnalysis::DiLepton&)>;
    using DiLeptonHybridCut = StringCutObjectSelector<TTWAnalysis::DiLepton>;
    // for in_category_pre_analyzer
    unsigned int m_nReqEl, m_nReqMu;
    int m_reqCh;
    // for in_category_post_analyzer
    DiLeptonHybridCut m_llIsInCateg;
    // list of dilepton working points
    std::vector<std::string> m_llWPs;
    // registered cuts
    std::vector<std::pair<std::string,std::function<bool(const TTWAnalysis::DiLepton&)>>> m_postCuts;
  public:
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
  private:
    mutable HLTMatch<TTWAnalysis::DiLepton> m_hltMatcher;

    bool diLeptonIsInCategory( const TTWAnalyzer* ttW, uint16_t dilepIdx ) const;
};

}
