#include <vector>
#include "Analyzer/LLPonAOD/plugins/Objects.h"

namespace {
  struct dictionary {

    //Structures                                                                
    JetType dummy0;
    CaloJetType dummy1;
    GenPType dummy2;

    //Vector of structures                                                      
    std::vector<JetType> dummyVector0;
    std::vector<CaloJetType> dummyVector1;
    std::vector<GenPType> dummyVector2;
  };
}
