#include <cp3_llbb/TTWAnalysis/interface/Types.h>
#include <vector>
#include <utility>

#include "cp3_llbb/TTWAnalysis/interface/NewTypes.h"

namespace TTWAnalysis {
  struct dictionary {
    TTWAnalysis::BaseObject dummy;
    std::vector<TTWAnalysis::BaseObject> dummy2;
    std::vector<uint16_t> dummy13;
    std::vector<std::vector<uint16_t>> dummy14;
    std::vector<float> dummy17;
    TTWAnalysis::GenParticle dummy22;
    std::vector<TTWAnalysis::GenParticle> dummy23;
    // new types
    TTWAnalysis::Lepton dummy24;
    TTWAnalysis::DiLepton dummy25;
  };
}
