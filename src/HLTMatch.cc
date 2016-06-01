#include "cp3_llbb/Framework/interface/HLTProducer.h"

#include "cp3_llbb/TTWAnalysis/interface/HLTMatch.h"

template<>
bool TTWAnalysis::HLTMatch<TTWAnalysis::Lepton>::operator() ( const TTWAnalysis::Lepton& lep ) const
{
  if ( ( ! m_hltProd )
    || ( lep.hlt_idx  < 0 ) || ( lep.hlt_idx  >= int(m_hltProd->object_paths.size()) ) )
  { return false; }

  return std::any_of( std::begin(m_regex), std::end(m_regex),
      [this,&lep] ( const boost::regex& paths ) -> bool {
        const auto& objp = m_hltProd->object_paths[lep.hlt_idx];
        return std::any_of(std::begin(objp), std::end(objp),
            [&paths] ( const std::string& obj ) {
              return boost::regex_match(obj, paths);
            } );
      } );
}

template<>
bool TTWAnalysis::HLTMatch<TTWAnalysis::DiLepton>::operator() ( const TTWAnalysis::DiLepton& diLep ) const
{
  if ( ( ! m_hltProd )
    || ( diLep.first->hlt_idx  < 0 ) || ( diLep.first->hlt_idx  >= int(m_hltProd->object_paths.size()) )
    || ( diLep.second->hlt_idx < 0 ) || ( diLep.second->hlt_idx >= int(m_hltProd->object_paths.size()) ) )
  { return false; }

  std::vector<std::string> mObj1, mObj2, mComm;
  for ( const auto& paths : m_regex ) {
    const auto& objp1 = m_hltProd->object_paths[diLep.first->hlt_idx];
    std::copy_if(std::begin(objp1), std::end(objp1), std::back_inserter(mObj1),
        [&paths] ( const std::string& obj ) { return boost::regex_match(obj, paths); } );
    const auto& objp2 = m_hltProd->object_paths[diLep.second->hlt_idx];
    std::copy_if(std::begin(objp2), std::end(objp2), std::back_inserter(mObj2),
        [&paths] ( const std::string& obj ) { return boost::regex_match(obj, paths); } );
  }
  std::set_intersection(std::begin(mObj1), std::end(mObj1), std::begin(mObj2), std::end(mObj2), std::back_inserter(mComm));
  return ! mComm.empty();
}
