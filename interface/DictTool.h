#ifndef TTWANALYSIS_DICTTOOL_H
#define TTWANALYSIS_DICTTOOL_H

#include <vector>
#include <string>
#include <boost/any.hpp>

#include "FWCore/PluginManager/interface/PluginFactory.h"

namespace edm {
  class ParameterSet;
  class ConsumesCollector;
  class Event;
  class EventSetup;
  template<class T> class Ptr;
};
class ProducersManager;
class AnalyzersManager;
class CategoryManager;
namespace pat {
  class Muon;
  class Electron;
  class Jet;
};

namespace TTWAnalysis {
  /*
   * Ordered list of key-value pairs, with heterogeneous values, implemented as a std::vector
   */
  class Dict {
    public:
      using entry = std::pair<std::string,boost::any>;
      using iterator = std::vector<entry>::const_iterator;

      // default constructors, assignments and destructors

      // read access
      std::size_t size() const { return m_data.size(); }
      iterator begin() const { return m_data.begin(); }
      iterator end() const { return m_data.end(); }
      const entry& operator[] ( std::size_t i ) const { return m_data[i]; }

      // write access
      template<class... Args>
      void add(Args&&... args) { m_data.emplace_back(std::forward<Args>(args)...); }

    private:
      std::vector<entry> m_data;
  };

  template<class T>
  class dict_argument_traits {
  public:
    using argument_type = const T&;
  };
  // specializations for pat:: objects
  template<> class dict_argument_traits<class pat::Muon>     { public: using argument_type = edm::Ptr<pat::Muon    >; };
  template<> class dict_argument_traits<class pat::Electron> { public: using argument_type = edm::Ptr<pat::Electron>; };
  template<> class dict_argument_traits<class pat::Jet>      { public: using argument_type = edm::Ptr<pat::Jet     >; };

  /*
   * Fill a list of values for an object
   *
   * The resulting dictionary should *always* have the same structure
   */
  template<class PatObject>
  class DictTool {
    public:
      using argument_type = typename dict_argument_traits<PatObject>::argument_type;

      explicit DictTool(const edm::ParameterSet& config) {}

      virtual ~DictTool() {}

      virtual void doConsumes(const edm::ParameterSet&, edm::ConsumesCollector&&) {}

      // interface method
      virtual Dict evaluate(argument_type, const edm::Event* = nullptr, const edm::EventSetup* = nullptr, const ProducersManager* = nullptr, const AnalyzersManager* = nullptr, const CategoryManager* = nullptr) const = 0;

      using factory = edmplugin::PluginFactory<DictTool<PatObject>* (const edm::ParameterSet&)>;
  };
};

#endif // TTWANALYSIS_DICTTOOL_H
