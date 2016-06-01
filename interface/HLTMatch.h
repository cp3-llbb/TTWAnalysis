#ifndef TTWANALYSIS_HLTMATCH_H
#define TTWANALYSIS_HLTMATCH_H

#include <boost/regex.hpp>

#include "cp3_llbb/TTWAnalysis/interface/NewTypes.h"
class HLTProducer;

namespace TTWAnalysis {
  /*
   * Helper functor for HLT candidate matching (based on HLTProducer outputs)
   */
  template<class Candidate>
  class HLTMatch {
  public:
    HLTMatch() {}
    explicit HLTMatch(const std::vector<std::string>& regex) { addRegex(regex); }

    const HLTProducer* hltProd() const { return m_hltProd; }
    void setHLTProd( const HLTProducer* hltProd ) { m_hltProd = hltProd; }

    void addRegex( const std::vector<std::string>& regex )
    {
      std::transform(std::begin(regex), std::end(regex), std::back_inserter(m_regex),
          [] ( const std::string& ireg ) { return boost::regex(ireg, boost::regex_constants::icase); } );
    }

    bool operator() ( const Candidate& cand ) const;

  private:
    const HLTProducer* m_hltProd = nullptr;
    std::vector<boost::regex> m_regex;
  };

  template<>
  bool HLTMatch<Lepton  >::operator()( const Lepton& lepton ) const;
  template<>
  bool HLTMatch<DiLepton>::operator()( const DiLepton& diLepton ) const;
}
#endif // TTWANALYSIS_HLTMATCH_H
