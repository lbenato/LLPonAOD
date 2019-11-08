#ifndef OBJECTSFORMAT_H
#define OBJECTSFORMAT_H

#include <string>
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Utilities.h"
#include "Objects.h"

class ObjectsFormat {

    public:
        ObjectsFormat() {};
        ~ObjectsFormat() {};

        static void FillElectronType(LeptonType&, const pat::Electron*, bool);
        static void FillMuonType(LeptonType&, const pat::Muon*, bool);
        static void FillPhotonType(PhotonType&, const pat::Photon*, bool);
        static void FillTauType(TauType&, const pat::Tau*, bool);
        static void FillJetType(JetType&, const pat::Jet*, bool);
        static void FillCaloJetType(CaloJetType&, const reco::CaloJet*, bool, bool);//new
        static void FillFatJetType(FatJetType&, const pat::Jet*, std::string, bool);
        static void FillCustomFatJetType(CustomFatJetType&, const pat::Jet*, std::string, bool);
        static void FillMEtType(MEtType&, const pat::MET*, bool);
//        static void FillMEtCorType(MEtCorType&, const pat::MET*, bool);
        static void FillMEtFullType(MEtFullType&, const pat::MET*, bool);
        static void FillCandidateType(CandidateType&, pat::CompositeCandidate*, bool);
        static void FillLorentzType(LorentzType&, const reco::Candidate::LorentzVector*);
        static void FillGenPType(GenPType&, const reco::GenParticle*);
        static void FillTriggerObjectType(TriggerObjectType&, const pat::TriggerObjectStandAlone*);
        static void FillPFCandidateType(PFCandidateType&, const pat::PackedCandidate*, int, int, int);
        static void FillPrimVertexType(VertexType&, const reco::Vertex*);
	static void FillSecVertexType(VertexType&, const reco::VertexCompositePtrCandidate*);
	static void FillBtagSecVertexType(VertexType&, const reco::CandSecondaryVertexTagInfo*, const reco::CandIPTagInfo*, unsigned int, float, int);
        static void FillSimplifiedJetType(SimplifiedJetType&, const pat::Jet*, bool);

        static void ResetLeptonType(LeptonType&);
        static void ResetPhotonType(PhotonType&);
        static void ResetTauType(TauType&);
        static void ResetJetType(JetType&);
        static void ResetCaloJetType(CaloJetType&);//new
        static void ResetFatJetType(FatJetType&);
        static void ResetCustomFatJetType(CustomFatJetType&);
        static void ResetMEtType(MEtType&);
//        static void ResetMEtCorType(MEtCorType&);
        static void ResetMEtFullType(MEtFullType&);
        static void ResetCandidateType(CandidateType&);
        static void ResetLorentzType(LorentzType&);
        static void ResetGenPType(GenPType&);
        static void ResetTriggerObjectType(TriggerObjectType&);
        static void ResetPFCandidateType(PFCandidateType&);
	static void ResetVertexType(VertexType&);
        static void ResetSimplifiedJetType(SimplifiedJetType&);

        static std::string ListLeptonType();
        static std::string ListPhotonType();
        static std::string ListTauType();
        static std::string ListJetType();
        static std::string ListCaloJetType();
        static std::string ListFatJetType();
        static std::string ListCustomFatJetType();
        static std::string ListMEtType();
//        static std::string ListMEtCorType();
        static std::string ListMEtFullType();
        static std::string ListCandidateType();
        static std::string ListLorentzType();
        static std::string ListGenPType();
        static std::string ListTriggerObjectType();
        static std::string ListPFCandidateType();
        static std::string ListVertexType();
        static std::string ListSimplifiedJetType();


    private:

};


#endif
