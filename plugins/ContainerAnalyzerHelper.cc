#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "cp3_llbb/Framework/interface/AnalyzersManager.h"
#include "cp3_llbb/TTWAnalysis/interface/AnalyzerHelper.h"
#include "cp3_llbb/TTWAnalysis/interface/DictTool.h"
#include "cp3_llbb/TTWAnalysis/interface/CandidatesProducer.h"
#include "cp3_llbb/TTWAnalysis/interface/TTWAnalyzer.h"

namespace TTWAnalysis {
  /*
   * Fill a list of variables (using a set of TTWAnalysis::DictTools) for all objects in a container
   * e.g. id variables, kinematics etc. for a pat::Electron collection
   * (make sure that the necessary DictTool<PatObject> factory is declared to the plugin service)
   */
  template<class PatObject,class Container>
  class ContainerAnalyzerHelper : public AnalyzerHelper {
    public:
      ContainerAnalyzerHelper(const std::string& name, const ROOT::TreeGroup& tree, const edm::ParameterSet& config)
        : AnalyzerHelper(name, tree, config)
        , m_tree(tree)
      {
        this->m_ttWName = config.getParameter<std::string>("TTWAnalyzer");

        LogDebug("ttW") << "DelegatingAnalyzer " << this->name() << " dict tools: ";
        const edm::ParameterSet& dictConfigs = config.getParameter<edm::ParameterSet>("DictTools");
        for ( const std::string& dictToolName : dictConfigs.getParameterNames() ) {
          const auto& dictToolConfig = dictConfigs.getParameter<edm::ParameterSet>(dictToolName);
          std::string dictType{dictToolConfig.getParameter<std::string>("type")};
          const auto& dictParams = dictToolConfig.getParameter<edm::ParameterSet>("parameters");
          this->m_dicts.emplace_back(
                  std::unique_ptr<DictTool<PatObject>>{DictTool<PatObject>::factory::get()->create(dictType, dictParams)}
                , dictParams);
          LogDebug("ttW") << m_dicts.size()-1 << " : " << dictToolName;
        }
        // initialize branch-ref-group for each dict based on that
        LogDebug("ttW") << "Container analyzer helper " << this->name() << ": initializing branches";
        std::transform(std::begin(this->m_dicts), std::end(this->m_dicts), std::back_inserter(this->m_branchesRefs),
            [this] ( const DictHolder& tool ) {
              PatObject tmp{};
              return BranchRefsHolder(m_tree, tool.first->evaluate(tmp));
            });
        LogDebug("ttW") << "Container analyzer helper " << this->name() << " initialized";
      }

      virtual ~ContainerAnalyzerHelper() {}

      virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
      {
        AnalyzerHelper::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
        for ( auto& idct : this->m_dicts ) {
          idct.first->doConsumes(idct.second, std::forward<edm::ConsumesCollector>(collector));
        }
      }

      virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup, const ProducersManager& prodMgr, const AnalyzersManager& anaMgr, const CategoryManager& catMgr) override
      {
        const TTWAnalyzer& ttW_ana = anaMgr.get<TTWAnalyzer>(this->m_ttWName);
        auto candAcc = _accessor<PatObject,Container>(ttW_ana);
        for ( std::size_t iO{0}; candAcc.size() != iO; ++iO ) {
          for ( std::size_t iT{0}; this->m_dicts.size() != iT; ++iT ) {
            LogDebug("ttW") << "Calling dict tool #" << iT;
            this->m_branchesRefs[iT].addDict(this->m_dicts[iT].first->evaluate(candAcc[iO], &event, &eventSetup, &prodMgr, &anaMgr, &catMgr));
          }
        }
      }

    private:
      ROOT::TreeGroup m_tree;
      std::string m_ttWName;
      using DictHolder = std::pair<std::unique_ptr<DictTool<PatObject>>,edm::ParameterSet>;
      std::vector<DictHolder> m_dicts;

      /*
       * Ref to a vector branch (interface taking boost::any)
       */
      class IBranchRef {
        public:
          virtual ~IBranchRef() {}

          virtual void addValue(const boost::any& val) = 0;
      };
      /*
       * Ref to a vector branch (generic implementation, cast boost::any)
       */
      template<typename T>
      class BranchRef : public IBranchRef
      {
        public:
          using value_type = T;

          explicit BranchRef(std::vector<T>& data) : m_data(data) {}

          virtual ~BranchRef() {}

          virtual void addValue(const boost::any& val) override
          { m_data.emplace_back(boost::any_cast<T>(val)); }

        private:
          std::vector<T>& m_data;
      };

      /*
       * Construct a list of references to vector branches (using a factory method)
       * and allow to fill the corresponding dictionar at once
       */
      class BranchRefsHolder {
        public:
          BranchRefsHolder( ROOT::TreeGroup& tree, const Dict& initDict )
          {
            for ( const auto& entry : initDict ) {
              const auto& etyp = entry.second.type();
              // add more supported types here
              if        ( typeid(bool)          == etyp ) {
                addBranch<bool>(tree, entry.first);
              } else if ( typeid(std::uint16_t) == etyp ) {
                addBranch<std::uint16_t>(tree, entry.first);
              } else if ( typeid(std::int16_t)  == etyp ) {
                addBranch<std::int16_t>(tree, entry.first);
              } else if ( typeid(std::int32_t)  == etyp ) {
                addBranch<std::int32_t>(tree, entry.first);
              } else if ( typeid(float)         == etyp ) {
                addBranch<float>(tree, entry.first);
              } else if ( typeid(double)        == etyp ) {
                addBranch<double>(tree, entry.first);
              } else if ( typeid(math::XYZPointD) == etyp ) {
                addBranch<math::XYZPointD>(tree, entry.first);
              } else if ( typeid(math::PtEtaPhiELorentzVectorF) == etyp ) {
                addBranch<math::PtEtaPhiELorentzVectorF>(tree, entry.first);
              } else {
                throw std::bad_cast{};
              }
            }
          }

          void addDict(const Dict& valDict)
          {
            assert(m_branches.size() == valDict.size());
            for( std::size_t i{0}; i < valDict.size(); ++i ) {
              m_branches[i]->addValue(valDict[i].second);
            }
          }
        private:
          std::vector<std::unique_ptr<IBranchRef>> m_branches;

          template<typename T>
          void addBranch(ROOT::TreeGroup& tree, const std::string& name)
          {
            m_branches.emplace_back(new BranchRef<T>(tree[name].write<std::vector<T>>()));
          }
      };

      std::vector<BranchRefsHolder> m_branchesRefs;

      // helper to access the different kinds of containers
      template<typename OBJ, typename CONT>
      class _accessor {
      public:
        explicit _accessor( const TTWAnalyzer& ttW );
        std::size_t size() const { return m_cont.size(); }
        const OBJ& operator[] ( std::size_t i ) const;
      private:
        const CONT& m_cont;
      };
      template<typename OBJ>
      class _accessor<OBJ,edm::PtrVector<OBJ>> {
      public:
        explicit _accessor( const TTWAnalyzer& ttW ) : m_cont(ttW.getPtrList<OBJ>()) {}
        std::size_t size() const { return m_cont.size(); }
        const OBJ& operator[] ( std::size_t i ) const { return *(m_cont[i].get()); }
      private:
        const edm::PtrVector<OBJ>& m_cont;
      };
      template<typename OBJ>
      class _accessor<OBJ,std::vector<OBJ>> {
      public:
        explicit _accessor( const TTWAnalyzer& ttW ) : m_cont(ttW.getObjList<OBJ>()) {}
        std::size_t size() const { return m_cont.size(); }
        const OBJ& operator[] ( std::size_t i ) const { return m_cont[i]; }
      private:
        const std::vector<OBJ>& m_cont;
      };
  };
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

namespace reco {
  class Vertex;
};
namespace pat {
  class Electron;
  class Muon;
  class Jet;
};

namespace TTWAnalysis {
  struct Lepton;
  struct DiLepton;
  struct DiJet;
  struct DiLeptonDiJet;
  struct DiLeptonDiJetMet;
};

namespace TTWAnalysis {
  using VerticesAnalyzerHelper  = ContainerAnalyzerHelper<reco::Vertex ,std::vector<reco::Vertex    >>;

  using ElectronsAnalyzerHelper = ContainerAnalyzerHelper<pat::Electron,edm::PtrVector<pat::Electron>>;
  using MuonsAnalyzerHelper     = ContainerAnalyzerHelper<pat::Muon    ,edm::PtrVector<pat::Muon    >>;
  using JetsAnalyzerHelper      = ContainerAnalyzerHelper<pat::Jet     ,edm::PtrVector<pat::Jet     >>;

  using LeptonsAnalyzerHelper           = ContainerAnalyzerHelper<TTWAnalysis::Lepton          ,std::vector<TTWAnalysis::Lepton          >>;
  using DiLeptonsAnalyzerHelper         = ContainerAnalyzerHelper<TTWAnalysis::DiLepton        ,std::vector<TTWAnalysis::DiLepton        >>;
  using DiJetsAnalyzerHelper            = ContainerAnalyzerHelper<TTWAnalysis::DiJet           ,std::vector<TTWAnalysis::DiJet           >>;
  using DiLeptonDiJetsAnalyzerHelper    = ContainerAnalyzerHelper<TTWAnalysis::DiLeptonDiJet   ,std::vector<TTWAnalysis::DiLeptonDiJet   >>;
  using DiLeptonDiJetMetsAnalyzerHelper = ContainerAnalyzerHelper<TTWAnalysis::DiLeptonDiJetMet,std::vector<TTWAnalysis::DiLeptonDiJetMet>>;
};

DEFINE_EDM_PLUGIN(TTWAnalyzerHelperFactory, TTWAnalysis::VerticesAnalyzerHelper , "ttw_verticesanalyzerhelper");

DEFINE_EDM_PLUGIN(TTWAnalyzerHelperFactory, TTWAnalysis::ElectronsAnalyzerHelper, "ttw_electronsanalyzerhelper");
DEFINE_EDM_PLUGIN(TTWAnalyzerHelperFactory, TTWAnalysis::MuonsAnalyzerHelper    , "ttw_muonsanalyzerhelper");
DEFINE_EDM_PLUGIN(TTWAnalyzerHelperFactory, TTWAnalysis::JetsAnalyzerHelper     , "ttw_jetsanalyzerhelper");

DEFINE_EDM_PLUGIN(TTWAnalyzerHelperFactory, TTWAnalysis::LeptonsAnalyzerHelper          , "ttw_leptonsanalyzerhelper");
DEFINE_EDM_PLUGIN(TTWAnalyzerHelperFactory, TTWAnalysis::DiLeptonsAnalyzerHelper        , "ttw_dileptonsanalyzerhelper");
DEFINE_EDM_PLUGIN(TTWAnalyzerHelperFactory, TTWAnalysis::DiJetsAnalyzerHelper           , "ttw_dijetsanalyzerhelper");
DEFINE_EDM_PLUGIN(TTWAnalyzerHelperFactory, TTWAnalysis::DiLeptonDiJetsAnalyzerHelper   , "ttw_dileptondijetsanalyzerhelper");
DEFINE_EDM_PLUGIN(TTWAnalyzerHelperFactory, TTWAnalysis::DiLeptonDiJetMetsAnalyzerHelper, "ttw_dileptondijetmetsanalyzerhelper");
