#ifndef OBJECTSFORMAT_H
#define OBJECTSFORMAT_H

#include <string>
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "Objects.h"

class ObjectsFormat {

 public:
  ObjectsFormat() {};
  ~ObjectsFormat() {};

  static void FillJetType(JetType&, const reco::Jet*, bool);
  static void FillCaloJetType(CaloJetType&, const reco::CaloJet*, bool, bool);//new
  static void FillMEtType(MEtType&, const reco::PFMET*, bool);
  static void FillCandidateType(CandidateType&, reco::CompositeCandidate*, bool);
  static void FillGenPType(GenPType&, const reco::GenParticle*);

  static void ResetJetType(JetType&);
  static void ResetCaloJetType(CaloJetType&);//new
  static void ResetMEtType(MEtType&);
  static void ResetCandidateType(CandidateType&);
  static void ResetGenPType(GenPType&);

  static std::string ListJetType();
  static std::string ListCaloJetType();
  static std::string ListMEtType();
  static std::string ListCandidateType();
  static std::string ListGenPType();

 private:

};

#endif
