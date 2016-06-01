#ifndef TTWANALYSIS_CANDIDATES_PRODUCER_H
#define TTWANALYSIS_CANDIDATES_PRODUCER_H

#include "DataFormats/Common/interface/PtrVector.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "cp3_llbb/Framework/interface/Producer.h"

namespace TTWAnalysis {

  template<typename ObjectType>
  class CandidatesProducer : public Framework::Producer {
    public:
      CandidatesProducer(const std::string& name, const ROOT::TreeGroup& tree, const edm::ParameterSet& config)
        : Producer(name, tree, config)
        , m_cut(config.getUntrackedParameter<std::string>("cut", "1"))
      { }

      virtual ~CandidatesProducer() {}

      virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
      {
        m_input_token = collector.consumes<std::vector<ObjectType>>(config.getParameter<edm::InputTag>("input"));
      }

      virtual void produce(edm::Event& event, const edm::EventSetup& eventSetup) override
      {
        edm::Handle<std::vector<ObjectType>> candidates;
        event.getByToken(m_input_token, candidates);

        m_selected.clear();
        for ( std::size_t i{0}; candidates->size() != i; ++i ) {
          if ( m_cut((*candidates)[i]) ) {
            m_selected.push_back(edm::Ptr<ObjectType>{candidates.product(), i});
          }
        }
      }

      const edm::PtrVector<ObjectType>& selected() const { return m_selected; }

    private:
      edm::EDGetTokenT<std::vector<ObjectType>> m_input_token;
      StringCutObjectSelector<ObjectType> m_cut;
      edm::PtrVector<ObjectType> m_selected;
  };

};

#endif // TTWANALYSIS_CANDIDATES_PRODUCER_H
