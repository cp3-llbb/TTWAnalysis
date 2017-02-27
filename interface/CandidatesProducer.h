#ifndef TTWANALYSIS_CANDIDATES_PRODUCER_H
#define TTWANALYSIS_CANDIDATES_PRODUCER_H

#include "DataFormats/Common/interface/PtrVector.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "cp3_llbb/Framework/interface/Producer.h"
#include "cp3_llbb/TTWAnalysis/interface/DictTool.h"

namespace TTWAnalysis {

  template<typename ObjectType>
  class CandidatesProducer : public Framework::Producer {
    public:
      CandidatesProducer(const std::string& name, const ROOT::TreeGroup& tree, const edm::ParameterSet& config)
        : Producer(name, tree, config)
        , m_cut(config.getUntrackedParameter<std::string>("cut", "1"))
      {
        LogDebug("ttW") << "Constructing candidates producer " << m_name;
        const edm::ParameterSet& dictConfigs = config.getParameter<edm::ParameterSet>("DictTools");
        for ( const std::string& dictToolName : dictConfigs.getParameterNames() ) {
          const auto& dictToolConfig = dictConfigs.getParameter<edm::ParameterSet>(dictToolName);
          std::string dictType{dictToolConfig.getParameter<std::string>("type")};
          const auto& dictParams = dictToolConfig.getParameter<edm::ParameterSet>("parameters");
          LogDebug("ttW") << "Adding dict tool " << dictToolName << " to CandidatesProducer " << m_name;
          this->m_dicts.emplace_back(
                  std::unique_ptr<DictTool<ObjectType>>{DictTool<ObjectType>::factory::get()->create(dictType, dictParams)}
                , dictParams);
        }
        LogDebug("ttW") << "Added all dict tools";
      }

      virtual ~CandidatesProducer() {}

      virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
      {
        m_input_token = collector.consumes<std::vector<ObjectType>>(config.getParameter<edm::InputTag>("input"));
        for ( auto& idct : this->m_dicts ) {
          idct.first->doConsumes(idct.second, std::forward<edm::ConsumesCollector>(collector));
        }
      }

      virtual void produce(edm::Event& event, const edm::EventSetup& eventSetup) override
      {
        edm::Handle<std::vector<ObjectType>> candidates;
        event.getByToken(m_input_token, candidates);

        m_selected.clear();
        for ( std::size_t i{0}; candidates->size() != i; ++i ) {
          if ( m_cut((*candidates)[i]) ) {
            m_selected.push_back(edm::Ptr<ObjectType>{candidates.product(), i});
            // calculate and add some values
            auto& cand = const_cast<ObjectType&>(candidates->at(i));
            LogDebug("ttW") << m_name << " evaluating dicts for selected candidate";
            for ( std::size_t iT{0}; this->m_dicts.size() != iT; ++iT ) {
              auto dct = this->m_dicts[iT].first->evaluate(cand, &event, &eventSetup);
              for ( const auto& elm : dct ) {
                cand.addUserFloat(elm.first, boost::any_cast<double>(elm.second));
                LogDebug("ttW") << m_name << " added user float " << elm.first << " : " << boost::any_cast<double>(elm.second);
              }
            }
          }
        }
      }

      const edm::PtrVector<ObjectType>& selected() const { return m_selected; }

    private:
      edm::EDGetTokenT<std::vector<ObjectType>> m_input_token;
      StringCutObjectSelector<ObjectType> m_cut;
      edm::PtrVector<ObjectType> m_selected;
      // floats to add
      using DictHolder = std::pair<std::unique_ptr<DictTool<ObjectType>>,edm::ParameterSet>;
      std::vector<DictHolder> m_dicts;
  };

};

#endif // TTWANALYSIS_CANDIDATES_PRODUCER_H
