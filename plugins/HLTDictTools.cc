#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/Event.h"

#include "cp3_llbb/Framework/interface/Types.h"
#include "cp3_llbb/TTWAnalysis/interface/DictTool.h"
#include "cp3_llbb/TTWAnalysis/interface/NewTypes.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include <cp3_llbb/Framework/interface/HLTService.h>

#include <Math/VectorUtil.h>
using namespace ROOT::Math;

#include <boost/regex.hpp>

namespace TTWAnalysis {

template<class Candidate>
class DictHLTMatchv2 : public DictTool<Candidate> {
public:
  DictHLTMatchv2(const edm::ParameterSet& config)
    : DictTool<Candidate>(config)
    , m_hltDRCut( config.getUntrackedParameter<double>("hltDRCut", std::numeric_limits<float>::max()) )
    , m_hltDPtCut( config.getUntrackedParameter<double>("hltDPtCut", std::numeric_limits<float>::max()) )
  {
    if (config.exists("triggers")) {
      m_hlt_service.reset(new HLTService(config.getUntrackedParameter<edm::FileInPath>("triggers").fullPath()));
      std::cout << "HLT configuration: " << std::endl;
      m_hlt_service->print();
      std::cout << std::endl;
    }

    const auto& selsConfig = config.getParameter<edm::ParameterSet>("Selections");
    for ( const auto& iName : selsConfig.getParameterNames() ) {
      const auto& sels = selsConfig.getParameter<std::vector<std::string>>(iName);
      std::vector<boost::regex> patterns; patterns.reserve(sels.size());
      std::transform(std::begin(sels), std::end(sels), std::back_inserter(patterns), 
          [] ( const std::string& ireg ) { return boost::regex(ireg, boost::regex_constants::icase); });
      m_selections.emplace(iName, patterns);
    }

  }
  virtual ~DictHLTMatchv2() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override {
    m_hlt_token = collector.consumes<edm::TriggerResults>(config.getUntrackedParameter<edm::InputTag>("hlt", edm::InputTag("TriggerResults", "", "HLT")));
    // m_prescales_token = collector.consumes<pat::PackedTriggerPrescales>(config.getUntrackedParameter<edm::InputTag>("prescales", edm::InputTag("patTrigger")));
    m_trigger_objects_token = collector.consumes<pat::TriggerObjectStandAloneCollection>(config.getUntrackedParameter<edm::InputTag>("objects", edm::InputTag("selectedPatTrigger")));
  }

  using FilteredObject = std::pair<pat::TriggerObjectStandAlone,std::vector<std::string>>;

  virtual Dict evaluate(const Candidate& cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  std::map<std::string,std::vector<boost::regex>> m_selections;
  // Tokens
  edm::EDGetTokenT<edm::TriggerResults> m_hlt_token;
  // edm::EDGetTokenT<pat::PackedTriggerPrescales> m_prescales_token;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> m_trigger_objects_token;
  // Service
  std::shared_ptr<HLTService> m_hlt_service;

  const float m_hltDRCut, m_hltDPtCut;

  mutable std::vector<std::string> m_paths;
  // mutable std::vector<uint16_t> m_prescales;

  void readTriggerNamesAndPrescales( const edm::Event* event ) const
  {
    m_paths.clear(); // m_prescales.clear();

    edm::Handle<edm::TriggerResults> hlt;
    event->getByToken(m_hlt_token, hlt);
    // edm::Handle<pat::PackedTriggerPrescales> prescales;
    //event->getByToken(m_prescales_token, prescales);

    const edm::TriggerNames& triggerNames = event->triggerNames(*hlt);

    bool filter = m_hlt_service.get() != nullptr;
    const HLTService::PathVector* valid_paths = nullptr;
    if (filter) {
      valid_paths = &m_hlt_service->getPaths(event->id().run());
    }

    // get paths and prescales
    for ( size_t i = 0 ; hlt->size() != i; ++i ) {
      if (hlt->accept(i)) {
        std::string triggerName = triggerNames.triggerName(i);
        if (triggerName == "HLTriggerFinalPath")
          continue; // This one is pretty useless...
        if (triggerName[0] == 'A')
          continue; // Remove AlCa HLT paths

        if ( (!filter) || std::any_of(valid_paths->begin(), valid_paths->end(),
          [&triggerName] ( const HLTService::PathName& regex ) { return boost::regex_match(triggerName, regex); } ) )
        {
          m_paths.push_back(triggerName);
          // m_prescales.push_back(prescales.isValid() ? prescales->getPrescaleForIndex(i) : 1.);
        }
      }
    }
    std::sort(std::begin(m_paths), std::end(m_paths));
  }

  std::vector<FilteredObject> collectTriggerObjects( const edm::Event* event ) const
  {
    edm::Handle<edm::TriggerResults> hlt;
    event->getByToken(m_hlt_token, hlt);
    const edm::TriggerNames& triggerNames = event->triggerNames(*hlt);
    bool filter = m_hlt_service.get() != nullptr;
    edm::Handle<pat::TriggerObjectStandAloneCollection> objects;
    event->getByToken(m_trigger_objects_token, objects);

    std::vector<FilteredObject> result;
    if ( objects.isValid() ) {
      for ( pat::TriggerObjectStandAlone obj : *objects ) {
        obj.unpackPathNames(triggerNames);

        std::vector<std::string> object_paths = obj.pathNames(false);

        // Check if this object has triggered at least one of the path we are interesting in
        if ( (!filter) || std::any_of(std::begin(object_paths), std::end(object_paths),
          [this] ( const std::string& path ) { return std::end(m_paths) != std::find(std::begin(m_paths), std::end(m_paths), path); } ) )
        {
          std::sort(std::begin(object_paths), std::end(object_paths));

          std::vector<std::string> filtered_paths;
          std::set_intersection(m_paths.begin(), m_paths.end(), object_paths.begin(), object_paths.end(), std::back_inserter(filtered_paths));

          result.emplace_back(obj, filtered_paths);
        }
      }
    }
    return result;
  }
};

bool matchHLTCandidate( const Lepton& lepton, const DictHLTMatchv2<Lepton>::FilteredObject& trigObj, const boost::regex& sel, float maxDR, float maxDPT )
{
  return ( std::any_of(std::begin(trigObj.second), std::end(trigObj.second),
        [&sel] ( const std::string& trig ) {
          return boost::regex_match(trig, sel);
        })
       && ( VectorUtil::DeltaR(lepton.p4(), trigObj.first.p4()) <= maxDR )
       && ( std::abs(lepton.p4().Pt() - trigObj.first.p4().Pt()) / lepton.p4().Pt() <= maxDPT )
       );
}

template<class Candidate>
bool hasHLTCandidateMatch( const Candidate& cand, const std::vector<typename DictHLTMatchv2<Candidate>::FilteredObject>& trigObjects, const std::vector<boost::regex>& selections, float maxDR, float maxDPT );

template<>
bool hasHLTCandidateMatch<Lepton>( const Lepton& lepton, const std::vector<typename DictHLTMatchv2<Lepton>::FilteredObject>& trigObjects, const std::vector<boost::regex>& selections, float maxDR, float maxDPT )
{
  return std::any_of(std::begin(trigObjects), std::end(trigObjects),
      [&,maxDR,maxDPT] ( const DictHLTMatchv2<Lepton>::FilteredObject& obj ) {
        return std::any_of(std::begin(selections), std::end(selections),
          [&,maxDR,maxDPT] ( const boost::regex& sel ) {
            return matchHLTCandidate(lepton, obj, sel, maxDR, maxDPT);
          });
        });
}

template<>
bool hasHLTCandidateMatch<DiLepton>( const DiLepton& dilepton, const std::vector<typename DictHLTMatchv2<DiLepton>::FilteredObject>& trigObjects, const std::vector<boost::regex>& selections, float maxDR, float maxDPT )
{
  std::vector<std::string> mObj1, mObj2, mComm;
  for ( const auto& sel : selections ) {
    for ( const auto& obj : trigObjects ) {
      if ( dilepton.first  && matchHLTCandidate(*(dilepton.first), obj, sel, maxDR, maxDPT ) ) {
        std::copy_if(std::begin(obj.second), std::end(obj.second), std::back_inserter(mObj1),
            [&sel] ( const std::string& obj ) { return boost::regex_match(obj, sel); } );
      }
      if ( dilepton.second && matchHLTCandidate(*(dilepton.second), obj, sel, maxDR, maxDPT ) ) {
        std::copy_if(std::begin(obj.second), std::end(obj.second), std::back_inserter(mObj2),
            [&sel] ( const std::string& obj ) { return boost::regex_match(obj, sel); } );
      }
    }
  }
  std::set_intersection(std::begin(mObj1), std::end(mObj1), std::begin(mObj2), std::end(mObj2), std::back_inserter(mComm));
  return ! mComm.empty();
}

template<class Candidate>
Dict DictHLTMatchv2<Candidate>::evaluate(const Candidate& cand,
    const edm::Event* event, const edm::EventSetup* /**/,
    const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  std::vector<FilteredObject> triggerObjects;
  if ( event ) {
    // TODO room for optimisation
    readTriggerNamesAndPrescales(event);
    triggerObjects = collectTriggerObjects(event);
  }

  Dict ret{};
  for ( const auto& selection : m_selections ) {
    ret.add(selection.first, hasHLTCandidateMatch<Candidate>(cand, triggerObjects, selection.second, m_hltDRCut, m_hltDPtCut));
  }
  return ret;
}

}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<TTWAnalysis::Lepton  >::factory, TTWAnalysis::DictHLTMatchv2<TTWAnalysis::Lepton>  , "ttw_leptonHLTMatchv2");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<TTWAnalysis::DiLepton>::factory, TTWAnalysis::DictHLTMatchv2<TTWAnalysis::DiLepton>, "ttw_dileptonHLTMatchv2");
