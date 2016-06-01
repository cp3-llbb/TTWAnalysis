#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "cp3_llbb/TTWAnalysis/interface/CandidatesProducer.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

DEFINE_EDM_PLUGIN(ExTreeMakerProducerFactory, TTWAnalysis::CandidatesProducer<pat::Electron>, "ttw_electronproducer");
DEFINE_EDM_PLUGIN(ExTreeMakerProducerFactory, TTWAnalysis::CandidatesProducer<pat::Muon>, "ttw_muonproducer");
DEFINE_EDM_PLUGIN(ExTreeMakerProducerFactory, TTWAnalysis::CandidatesProducer<pat::Jet>, "ttw_jetproducer");
