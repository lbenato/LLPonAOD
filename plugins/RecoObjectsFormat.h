#ifndef RECOOBJECTSFORMAT_H
#define RECOOBJECTSFORMAT_H

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

#include "RecoObjects.h"

class RecoObjectsFormat {

 public:
  RecoObjectsFormat() {};
  ~RecoObjectsFormat() {};

  static void FillRecoJetType(RecoJetType&, const reco::Jet*, bool);
  static void FillRecoElectronType(RecoLeptonType&, const reco::GsfElectron*, bool);
  //static void FillCaloJetType(CaloJetType&, const reco::CaloJet*, bool, bool);//new
  static void FillRecoMEtType(RecoMEtType&, const reco::PFMET*, bool);
  //static void FillCandidateType(CandidateType&, reco::CompositeCandidate*, bool);
  //static void FillGenPType(GenPType&, const reco::GenParticle*);

  static void ResetRecoJetType(RecoJetType&);
  static void ResetRecoElectronType(RecoLeptonType&);
  //static void ResetCaloJetType(CaloJetType&);//new
  static void ResetRecoMEtType(RecoMEtType&);
  //static void ResetCandidateType(CandidateType&);
  //static void ResetGenPType(GenPType&);

  static std::string ListRecoJetType();
  static std::string ListRecoElectronType();
  //static std::string ListCaloJetType();
  static std::string ListRecoMEtType();
  //static std::string ListCandidateType();
  //static std::string ListGenPType();

 private:

};

#endif
