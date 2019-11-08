#include <vector>
#include "Analyzer/LLPonAOD/plugins/RecoObjects.h"
#include "Analyzer/LLPonAOD/plugins/Objects.h"

namespace {
  struct dictionary {

    //Structures                                                                
    RecoJetType dummy0;
    CaloJetType dummy1;
    GenPType dummy2;
    JetType dummy3;

    //Vector of structures                                                      
    std::vector<RecoJetType> dummyVector0;
    std::vector<CaloJetType> dummyVector1;
    std::vector<GenPType> dummyVector2;
    std::vector<JetType> dummyVector3;
  };
}
