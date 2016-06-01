#ifndef TTWANALYSIS_ANALYZERHELPER_H
#define TTWANALYSIS_ANALYZERHELPER_H

namespace edm {
  class ParameterSet;
  class Run;
  class LuminosityBlock;
  class Event;
  class EventSetup;
  class ConsumesCollector;
};
class ProducersManager;
class AnalyzersManager;
class CategoryManager;
class MetadataManager;

#include "cp3_llbb/TreeWrapper/interface/TreeWrapper.h"

namespace TTWAnalysis {

  class AnalyzerHelper {
    public:
      AnalyzerHelper(const std::string& name, const ROOT::TreeGroup& tree, const edm::ParameterSet& config)
        : m_name(name)
        , m_tree(tree)
      {}

      virtual ~AnalyzerHelper();

      const std::string& name() const { return m_name; };
      const ROOT::TreeGroup& tree() const { return m_tree; }

      // interface method
      virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const AnalyzersManager&, const CategoryManager&) = 0;

      virtual void doConsumes(const edm::ParameterSet&, edm::ConsumesCollector&&) {}

      // event handlers
      virtual void beginJob(MetadataManager&) {};
      virtual void endJob(MetadataManager&) {};
      virtual void beginRun(const edm::Run&, const edm::EventSetup&) {};
      virtual void endRun(const edm::Run&, const edm::EventSetup&) {};
      virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) {};
      virtual void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) {};

    private:
      std::string m_name;
      ROOT::TreeGroup m_tree;
  };
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<TTWAnalysis::AnalyzerHelper* (const std::string&, const ROOT::TreeGroup&, const edm::ParameterSet&)> TTWAnalyzerHelperFactory;

#endif // TTWANALYSIS_ANALYZERHELPER_H
