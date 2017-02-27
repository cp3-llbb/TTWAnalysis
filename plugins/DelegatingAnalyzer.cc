#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "cp3_llbb/Framework/interface/Analyzer.h"
#include "cp3_llbb/TTWAnalysis/interface/AnalyzerHelper.h"

namespace TTWAnalysis {

class DelegatingAnalyzer: public Framework::Analyzer {
  public:
    DelegatingAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config);

    virtual ~DelegatingAnalyzer();

    virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const AnalyzersManager&, const CategoryManager&) override;
    virtual void doConsumes(const edm::ParameterSet&, edm::ConsumesCollector&& collector) override;

    // event handlers
    virtual void beginJob(MetadataManager&) override;
    virtual void endJob(MetadataManager&) override;
    virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
    virtual void endRun(const edm::Run&, const edm::EventSetup&) override;
    virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
    virtual void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
  private:
    std::vector<std::unique_ptr<TTWAnalysis::AnalyzerHelper>> m_helpers;
};

}

TTWAnalysis::DelegatingAnalyzer::DelegatingAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config)
  : Analyzer(name, tree_, config)
  , m_helpers{}
{
  LogDebug("ttW") << "Constructing helpers for delegating analyzer " << m_name;
  const edm::ParameterSet& helperConfigs = config.getParameter<edm::ParameterSet>("Helpers");
  for ( const std::string& helperName : helperConfigs.getParameterNames() ) {
    const auto& helperConfig = helperConfigs.getParameter<edm::ParameterSet>(helperName);
    LogDebug("ttW") << " - " << helperName;
    m_helpers.emplace_back(std::unique_ptr<TTWAnalysis::AnalyzerHelper>{
        TTWAnalyzerHelperFactory::get()->create(
          helperConfig.getParameter<std::string>("type")
        , helperName
        , tree_.group(helperConfig.getParameter<std::string>("prefix"))
        , helperConfig.getParameter<edm::ParameterSet>("parameters")
        )});
  }
  LogDebug("ttW") << "Constructed all helpers for delegating analyzer " << m_name;
}

TTWAnalysis::DelegatingAnalyzer::~DelegatingAnalyzer() {}

void TTWAnalysis::DelegatingAnalyzer::doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector)
{
  const edm::ParameterSet& helperConfigs = config.getParameter<edm::ParameterSet>("Helpers");
  for ( auto& hlp : m_helpers ) {
    hlp->doConsumes(helperConfigs.getParameter<edm::ParameterSet>(hlp->name()).getParameter<edm::ParameterSet>("parameters"), std::forward<edm::ConsumesCollector>(collector));
  }
}

void TTWAnalysis::DelegatingAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup, const ProducersManager& producers, const AnalyzersManager& analyzers, const CategoryManager& categories)
{
  for ( auto& hlp : m_helpers ) {
    LogDebug("ttW") << "-> calling analyze on helper " << hlp->name();
    hlp->analyze(event, setup, producers, analyzers, categories);
    LogDebug("ttW") << "<- end";
  }
}

// event handlers
void TTWAnalysis::DelegatingAnalyzer::beginJob(MetadataManager& mdMgr)
{
  for ( auto& hlp : m_helpers ) {
    hlp->beginJob(mdMgr);
  }
}
void TTWAnalysis::DelegatingAnalyzer::endJob(MetadataManager& mdMgr)
{
  for ( auto& hlp : m_helpers ) {
    hlp->endJob(mdMgr);
  }
}
void TTWAnalysis::DelegatingAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  for ( auto& hlp : m_helpers ) {
    hlp->beginRun(run, eventSetup);
  }
}
void TTWAnalysis::DelegatingAnalyzer::endRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  for ( auto& hlp : m_helpers ) {
    hlp->endRun(run, eventSetup);
  }
}
void TTWAnalysis::DelegatingAnalyzer::beginLuminosityBlock(const edm::LuminosityBlock& lumiBlock, const edm::EventSetup& eventSetup)
{
  for ( auto& hlp : m_helpers ) {
    hlp->beginLuminosityBlock(lumiBlock, eventSetup);
  }
}
void TTWAnalysis::DelegatingAnalyzer::endLuminosityBlock(const edm::LuminosityBlock& lumiBlock, const edm::EventSetup& eventSetup)
{
  for ( auto& hlp : m_helpers ) {
    hlp->endLuminosityBlock(lumiBlock, eventSetup);
  }
}

#include <FWCore/PluginManager/interface/PluginFactory.h>
DEFINE_EDM_PLUGIN(ExTreeMakerAnalyzerFactory, TTWAnalysis::DelegatingAnalyzer, "ttw_delegatinganalyzer");
