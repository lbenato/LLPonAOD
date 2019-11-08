#ifndef UTILITIES_H
#define UTILITIES_H


// dR and dPhi
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefToBase.h"

// Jets and PF
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "RecoParticleFlow/PFProducer/interface/Utils.h"
// Transient Track and IP
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
// Gen Info
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
// Association
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

// Pat
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
// Trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
// KinFitter
//#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
//#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
//#include "PhysicsTools/KinFitter/interface/TFitParticlePtEtaPhi.h"
//#include "PhysicsTools/KinFitter/interface/TFitParticleEScaledMomDev.h"
//#include "PhysicsTools/KinFitter/interface/TKinFitter.h"

#include "Numbers.h"

#include "TLorentzVector.h"
#include "TVector3.h"


//bool SortByPt(pat::Jet j1, pat::Jet j2) {return(j1.pt() > j2.pt());}
//bool SortByCSV(pat::Jet j1, pat::Jet j2) {return(j1.bDiscriminator("combinedSecondaryVertexBJetTags") > j2.bDiscriminator("combinedSecondaryVertexBJetTags"));}
//bool SortByCSVV1(pat::Jet j1, pat::Jet j2) {return(j1.bDiscriminator("combinedSecondaryVertexV1BJetTags") > j2.bDiscriminator("combinedSecondaryVertexV1BJetTags"));}
//bool SortByCSVSL(pat::Jet j1, pat::Jet j2) {return(j1.bDiscriminator("combinedSecondaryVertexSoftPFLeptonV1BJetTags") > j2.bDiscriminator("combinedSecondaryVertexSoftPFLeptonV1BJetTags"));}
//bool SortByCSVV1Reshaped(pat::Jet j1, pat::Jet j2) {return(j1.bDiscriminator("combinedSecondaryVertexV1BJetTagsReshaped") > j2.bDiscriminator("combinedSecondaryVertexV1BJetTagsReshaped"));}

class Utilities {
  public:
    Utilities() {};
    ~Utilities() {};
    
    // Angular
    static float ReturnCosThetaStar(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
    static float ReturnCosTheta1(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
    static float ReturnCosTheta2(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
    static float ReturnPhi(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
    static float ReturnPhi1(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
    // Kinematics and reconstruction
    static float PerformKinematicFit(pat::Jet*, pat::Jet*, reco::Candidate::LorentzVector*, reco::Candidate::LorentzVector*, float);
    static float RecoverNeutrinoPz(const reco::Particle::LorentzVector*, const reco::Particle::LorentzVector*);
    // Miscellanea
    static pat::CompositeCandidate RecoilMassFormula(pat::CompositeCandidate&, pat::MET&);
    static int FindMotherId(const reco::GenParticle*);
    
  private:
    
};


#endif
