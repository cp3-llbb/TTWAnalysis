#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "cp3_llbb/Framework/interface/Types.h"
#include "cp3_llbb/TTWAnalysis/interface/DictTool.h"

#include "cp3_llbb/Framework/interface/BinnedValuesJSONParser.h"

namespace TTWAnalysis {
/**
 * Scale factors as a function of PT and abs(ETA)
 */
template<class PatObject>
class DictScaleFactors : public DictTool<PatObject> {
public:
  DictScaleFactors(const edm::ParameterSet& config)
    : DictTool<PatObject>(config)
  {
    const auto& sfConfig = config.getUntrackedParameter<edm::ParameterSet>("scale_factors");
    for ( const auto& sfName : sfConfig.getParameterNames() ) {
      BinnedValuesJSONParser parser(sfConfig.getUntrackedParameter<edm::FileInPath>(sfName).fullPath());
      m_sf.emplace(sfName, std::move(parser.get_values()));
    }
  }
  virtual ~DictScaleFactors() {}

  virtual Dict evaluate(edm::Ptr<PatObject> obj,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    const bool valid{obj.isNonnull()};
    Dict ret{};
    const Parameters vars{ {BinningVariable::Eta, valid ? obj->eta() : -100.}, {BinningVariable::Pt, valid ? obj->pt() : -1.} };
    for ( const auto& sf : m_sf ) {
      const std::vector<float> sfWErrs = ( event && ( ! event->isRealData() ) ) ? sf.second.get(vars)
                                      : std::vector<float>({ 1., 0., 0. });
      ret.add("sf_"+sf.first       , sfWErrs[0]);
      ret.add("sf_"+sf.first+"_elo", sfWErrs[1]);
      ret.add("sf_"+sf.first+"_ehi", sfWErrs[2]);
    }
    return ret;
  }
private:
  std::map<std::string,BinnedValues> m_sf;
  // TODO can make generic in input variables as well by using functors
};
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictScaleFactors<pat::Electron>, "ttw_electronSF");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictScaleFactors<pat::Muon    >, "ttw_muonSF");

#include "cp3_llbb/Framework/interface/BTaggingScaleFactors.h"

namespace TTWAnalysis {
/**
 * B-tagging scale factors as a function of PT and abs(ETA)
 */
class DictJetScaleFactors : public DictTool<pat::Jet> {
public:
  DictJetScaleFactors(const edm::ParameterSet& config)
    : DictTool<pat::Jet>(config)
  {
    const auto& sfConfig = config.getUntrackedParameter<edm::ParameterSet>("scale_factors");
    for ( const auto& sfName : sfConfig.getParameterNames() ) {
      const auto& sfSet = sfConfig.getUntrackedParameter<edm::ParameterSet>(sfName);
      std::string algo{sfSet.getUntrackedParameter<std::string>("algorithm")};
      if ( Algorithm::UNKNOWN == BTaggingScaleFactors::string_to_algorithm(algo) ) {
#ifdef SF_DEBUG
        std::cout << "\tUnsupported b-tagging algorithm: " << algo << std::endl;
#endif
        continue;
      }
      std::string wp{sfSet.getUntrackedParameter<std::string>("working_point")};
      m_algos[algo].push_back(wp);
      std::vector<edm::ParameterSet> files = sfSet.getUntrackedParameter<std::vector<edm::ParameterSet>>("files");

      for (const auto& file_set: files) {
        std::string file = file_set.getUntrackedParameter<edm::FileInPath>("file").fullPath();
        std::string flavor = file_set.getUntrackedParameter<std::string>("flavor");

#ifdef SF_DEBUG
        std::cout << "\tAdding scale factor for algo: " << algo_str << "  wp: " << wp << "  flavor: " << flavor << " from file '" << file << "'" << std::endl;
#endif

        sf_key_type sf_key = std::make_tuple(BTaggingScaleFactors::string_to_algorithm(algo), BTaggingScaleFactors::string_to_flavor(flavor), wp);

        BinnedValuesJSONParser parser(file);
        m_sf.emplace(sf_key, std::move(parser.get_values()));
      }
#ifdef SF_DEBUG
      std::cout << std::endl;
#endif
    }
  }
  virtual ~DictJetScaleFactors() {}

  virtual Dict evaluate(edm::Ptr<pat::Jet> jet,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    const bool valid{jet.isNonnull()};
    Dict ret{};
    for ( const auto algoWithWPs : m_algos ) {
      const Algorithm algo = BTaggingScaleFactors::string_to_algorithm(algoWithWPs.first);
      const std::string algoShort{BTaggingScaleFactors::algorithm_to_string(algo)};
      for ( const auto& wp : algoWithWPs.second ) {
        const sf_key_type sf_key = std::make_tuple(algo, valid ? BTaggingScaleFactors::get_flavor(jet->hadronFlavour()) : Flavor::LIGHT, wp);
        const std::string sf_name{algoShort+"_"+wp};
        const Parameters vars{ {BinningVariable::Eta, valid ? jet->eta() : -100.}, {BinningVariable::Pt, valid ? jet->pt() : -1.}, {BinningVariable::BTagDiscri, valid ? jet->bDiscriminator(algoWithWPs.first) : -1.} };
        const std::vector<float> sfWErrs = ( event && ( ! event->isRealData() ) ) ? m_sf.find(sf_key)->second.get(vars)
                                        : std::vector<float>({ 1., 0., 0. });
        ret.add("sf_"+sf_name       , sfWErrs[0]);
        ret.add("sf_"+sf_name+"_elo", sfWErrs[1]);
        ret.add("sf_"+sf_name+"_ehi", sfWErrs[2]);
      }
    }
    return ret;
  }
private:
  using sf_key_type = BTaggingScaleFactors::sf_key_type;

  std::map<std::string, std::vector<std::string>> m_algos; // WPs per algorithm
  std::map<sf_key_type,BinnedValues> m_sf;
};
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Jet>::factory, TTWAnalysis::DictJetScaleFactors, "ttw_jetSF");
