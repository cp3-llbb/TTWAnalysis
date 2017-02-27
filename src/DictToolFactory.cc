#include "cp3_llbb/TTWAnalysis/interface/DictTool.h"

namespace reco {
  class Vertex;
}

EDM_REGISTER_PLUGINFACTORY(TTWAnalysis::DictTool<reco::Vertex>::factory, "TTWVertexDictToolFactory");

namespace pat {
  class Electron;
  class Muon;
  class Jet;
}

EDM_REGISTER_PLUGINFACTORY(TTWAnalysis::DictTool<pat::Electron>::factory, "TTWElectronDictToolFactory");
EDM_REGISTER_PLUGINFACTORY(TTWAnalysis::DictTool<pat::Muon>::factory, "TTWMuonDictToolFactory");
EDM_REGISTER_PLUGINFACTORY(TTWAnalysis::DictTool<pat::Jet>::factory, "TTWJetDictToolFactory");

namespace TTWAnalysis {
  struct Lepton;
  struct DiLepton;
  struct DiJet;
  struct DiLeptonDiJet;
  struct DiLeptonDiJetMet;
};

EDM_REGISTER_PLUGINFACTORY(TTWAnalysis::DictTool<TTWAnalysis::Lepton          >::factory, "TTWLeptonDictToolFactory");
EDM_REGISTER_PLUGINFACTORY(TTWAnalysis::DictTool<TTWAnalysis::DiLepton        >::factory, "TTWDiLeptonDictToolFactory");
EDM_REGISTER_PLUGINFACTORY(TTWAnalysis::DictTool<TTWAnalysis::DiJet           >::factory, "TTWDiJetDictToolFactory");
EDM_REGISTER_PLUGINFACTORY(TTWAnalysis::DictTool<TTWAnalysis::DiLeptonDiJet   >::factory, "TTWDiLeptonDiJetDictToolFactory");
EDM_REGISTER_PLUGINFACTORY(TTWAnalysis::DictTool<TTWAnalysis::DiLeptonDiJetMet>::factory, "TTWDiLeptonDiJetMetDictToolFactory");
