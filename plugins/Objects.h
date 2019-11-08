#ifndef OBJECTS_H
#define OBJECTS_H

struct LeptonType {
LeptonType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), inTrkPt(-1.), pfIso03(-1.), pfIso04(-1.), trkIso(-1.), miniIso(-1.), dxy(-99.), dz(-99.), ip3d(-99.), sip3d(-99.), nPixelHits(-1.), dPhi_met(-1.), charge(0), pdgId(0), isElectron(false), isMuon(false), isVeto(false), isLoose(false), isMedium(false), isTight(false), isHighPt(false), isTrackerHighPt(false), SSscale(-1.), SSsigma(-1.), SSscaleUnc(-1.), SSsigmaUncUp(-1.), SSsigmaUncDown(-1.), SScorr(-1.), energySScorr(-1.), energySScorrUncUp(-1.), energySScorrUncDown(-1.), ptSScorr(-1.), ptSScorrUncUp(-1.), ptSScorrUncDown(-1.), isMatched(false) {} // isHEEP(false), isMVANonTrigMedium(false), isMVANonTrigTight(false), isMVATrigMedium(false), isMVATrigTight(false)
    float pt;
    float eta;
    float phi;
    float mass;
    float energy;
    float inTrkPt;
    float pfIso03;
    float pfIso04;
    float trkIso;
    float miniIso;
    float dxy;
    float dz;
    float ip3d;
    float sip3d;
    float nPixelHits;
    float dPhi_met;
    int charge;
    int pdgId;
    bool isElectron;
    bool isMuon;
    bool isVeto;
    bool isLoose;
    bool isMedium;
    bool isTight;
//    bool isHEEP;
//    bool isMVANonTrigMedium;
//    bool isMVANonTrigTight;
//    bool isMVATrigMedium;
//    bool isMVATrigTight;
    bool isHighPt;
    bool isTrackerHighPt;

    float SSscale;
    float SSsigma;
    float SSscaleUnc;
    float SSsigmaUncUp;
    float SSsigmaUncDown;
    float SScorr;
    float energySScorr;
    float energySScorrUncUp;
    float energySScorrUncDown;
    float ptSScorr;
    float ptSScorrUncUp;
    float ptSScorrUncDown;

    bool isMatched;
};

struct PhotonType {
PhotonType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), pfIso(-1.), dz(-99.), charge(0), pdgId(0), isLoose(false), isMedium(false), isTight(false), isMVANonTrigMedium(false), isMatched(false) {}
    float pt;
    float eta;
    float phi;
    float mass;
    float energy;
    float pfIso;
    float dz;
    int charge;
    int pdgId;
    bool isLoose;
    bool isMedium;
    bool isTight;
    bool isMVANonTrigMedium;
    bool isMatched;
};

struct TauType {
TauType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), pfIso(-1.), dz(-99.), charge(0), pdgId(0), isLoose(false), isMedium(false), isTight(false), isMatched(false) {}
    float pt;
    float eta;
    float phi;
    float mass;
    float energy;
    float pfIso;
    float dz;
    int charge;
    int pdgId;
    bool isLoose;
    bool isMedium;
    bool isTight;
    bool isMatched;
};


struct JetType {
JetType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), energyRaw(-1.), ptRaw(-1.), ptUnc(-1.), CSV(-99.), CMVA(-99.),
//dPhi_met(-1.), dPhi_Jet1(-1.), puId(-1.),
//CSVR(-99.), CSVRUp(-99.), CSVRDown(-99.),
//CMVAR(-99.), CMVARUp(-99.), CMVARDown(-99.), QGLikelihood(-1.),
//energies
cHadE(-1.), nHadE (-1.), eleE(-1.), photonE(-1.), muE(-1.), nEmE(-1.), cEmE(-1.), cmuE(-1.),
//energy fractions
cHadEFrac(-1.), nHadEFrac(-1.), eleEFrac(-1.), photonEFrac(-1.), muEFrac(-1.), nEmEFrac(-1.), cEmEFrac(-1.), cmuEFrac(-1.),
//multiplicities
cHadMulti(-1.), nHadMulti(-1.), eleMulti(-1.), photonMulti(-1.), muMulti(-1.), cMulti(-1.), nMulti(-1.), npr(-1.),
//multiplicity fractions
cHadMultiFrac(-1.), nHadMultiFrac(-1.), eleMultiFrac(-1.), photonMultiFrac(-1.), muMultiFrac(-1.), cMultiFrac(-1.), nMultiFrac(-1.),

	   //emEFrac(-1.), emEinEB(-1.), emEinEE(-1.), emEinHF(-1.), EFracHad(-1.), hadEinHB(-1.), hadEinHE(-1.), hadEinHF(-1.), hadEinHO(-1.),
ptGenJ(-10.), etaGenJ(-4.), phiGenJ(-4.), massGenJ(-10.),ptGen(-10.), etaGen(-4.), phiGen(-4.), massGen(-10.), pdgIdGen(0.),
// ptLhe(-10.), etaLhe(-4.), phiLhe(-4.),
    partonFlavour(0), hadronFlavour(0), mother(0), isLoose(false), isMedium(false), isTight(false), isTightLepVeto(false), isCSVL(false), isCSVM(false), isCSVT(false),
//isMatched(false), dR_q1(1000), dR_q2(1000), dR_q3(1000), dR_q4(1000), m_q1(false), m_q2(false), m_q3(false), m_q4(false), dR_pi1(1000), dR_pi2(1000),
matchBquark(-1), matchLL(-1),
//original_jet_index(-1),
    isGenMatched(0), isVBFGenMatched(0),
    alphaMax(-100.), sigIP2DMedian(-100.), theta2DMedian(-100.), POCA_theta2DMedian(-100.), nPixelHitsMedian(-1.0), nHitsMedian(-1.0),
    hcalE(-100.), ecalE(-100.), FracCal(-100.), flightDist2d(-100.), flightDist2dError(-100.), flightDist3d(-100.), flightDist3dError(-100.), nSV(0), nSVCand(0), nVertexTracks(0), nSelectedTracks(0), dRSVJet(-100.), SV_x(-100.), SV_y(-100.), SV_z(-100.), SV_dx(-100.), SV_dy(-100.), SV_dz(-100.), nTracksSV(0), SV_mass(-100.),  isCaloTag(0),
//VBF_DisplacedJet40_VTightID_Hadronic_match(0), VBF_DisplacedJet40_VVTightID_Hadronic_match(0),
    ptJESUp (-1.), ptJESDown (-1.), ptJERUp (-1.), ptJERDown (-1.), tau1(-1.), tau2(-1.), tau3(-1.), nSubJets(-1), tau21(-1.), tau31(-1.), tau32(-1.), tau1_neutral(-1.), tau2_neutral(-1.), tau21_neutral(-1.), tau1_charged(-1.), tau2_charged(-1.), tau21_charged(-1.), TriggerMatched_VBFJet(0), TriggerMatched_DisplacedJet(0), TriggerMatched_TripleJet50(0),
    nConstituents (-1), nTrackConstituents (-1), nTracks0PixelHits(-1), nTracks1PixelHit(-1),nTracks2PixelHits(-1),nTracks3PixelHits(-1),nTracks4PixelHits(-1),nTracks5PixelHits(-1),nTracksAtLeast6PixelHits(-1),
    nTracksValidHitInBPix1(-1),nTracks0LostInnerHits(-1), nTracks1LostInnerHit(-1), nTracksAtLeast2LostInnerHits(-1), nMatchedGenBquarks(-1){}
    float pt;
    float eta;
    float phi;
    float mass;
    float energy;
    float energyRaw;
    float ptRaw;
    float ptUnc;
  //float dPhi_met;
  // float dPhi_Jet1;
  //float puId;
    float CSV;
  //float CSVR;
  //float CSVRUp;
  //float CSVRDown;
    float CMVA;
  //float CMVAR;
  //float CMVARUp;
  //float CMVARDown;
  //float QGLikelihood;
    float cHadE;
    float nHadE;
    float eleE;
    float photonE;
    float muE;
    float nEmE;
    float cEmE;
    float cmuE;
    float cHadEFrac;
    float nHadEFrac;
    float eleEFrac;
    float photonEFrac;
    float muEFrac;
    float nEmEFrac;
    float cEmEFrac;
    float cmuEFrac;
    float cHadMulti;
    float nHadMulti;
    float eleMulti;
    float photonMulti;
    float muMulti;
    float cMulti;
    float nMulti;
    float npr;
    float cHadMultiFrac;
    float nHadMultiFrac;
    float eleMultiFrac;
    float photonMultiFrac;
    float muMultiFrac;
    float cMultiFrac;
    float nMultiFrac;
  //float emEFrac;
  //float emEinEB;
  //float emEinEE;
  //float emEinHF;
  //float EFracHad;
  //float hadEinHB;
  //float hadEinHE;
  //float hadEinHF;
  //float hadEinHO;
    float ptGenJ;
    float etaGenJ;
    float phiGenJ;
    float massGenJ;
    float ptGen;
    float etaGen;
    float phiGen;
    float massGen;
    int pdgIdGen;
//    float ptLhe;
//    float etaLhe;
//    float phiLhe;
//    int chm;
//    int cm;
//    int nm;
    int partonFlavour;
    int hadronFlavour;
    int mother;
    bool isLoose;
    bool isMedium;
    bool isTight;
    bool isTightLepVeto;
    bool isCSVL;
    bool isCSVM;
    bool isCSVT;
  //bool isMatched;
  //float dR_q1;
  //float dR_q2;
  //float dR_q3;
  //float dR_q4;
  //bool m_q1;
  //bool m_q2;
  //bool m_q3;
  //bool m_q4;
  //float dR_pi1;
  //float dR_pi2;
    int matchBquark;
    int matchLL;
  //int original_jet_index;
    int isGenMatched;
    int isVBFGenMatched;
    float alphaMax;
    float sigIP2DMedian;
    float theta2DMedian;
    float POCA_theta2DMedian;
    float nPixelHitsMedian;
    float nHitsMedian;
    float hcalE;
    float ecalE;
    float FracCal;
    float flightDist2d;
    float flightDist2dError;
    float flightDist3d;
    float flightDist3dError;
    int nSV;
    int nSVCand;
    int nVertexTracks;
    int nSelectedTracks;
    float dRSVJet;
    float SV_x;
    float SV_y;
    float SV_z;
    float SV_dx;
    float SV_dy;
    float SV_dz;
    int nTracksSV;
    float SV_mass;
    int isCaloTag;
  //int VBF_DisplacedJet40_VTightID_Hadronic_match;
  //int VBF_DisplacedJet40_VVTightID_Hadronic_match;
    float ptJESUp;
    float ptJESDown;
    float ptJERUp;
    float ptJERDown;
    float tau1;
    float tau2;
    float tau3;
    int nSubJets;
    float tau21;
    float tau31;
    float tau32;
    float tau1_neutral;
    float tau2_neutral;
    float tau21_neutral;
    float tau1_charged;
    float tau2_charged;
    float tau21_charged;
    int TriggerMatched_VBFJet;
    int TriggerMatched_DisplacedJet;
    int TriggerMatched_TripleJet50;
    int nConstituents;
    int nTrackConstituents;
    int nTracks0PixelHits;
    int nTracks1PixelHit;
    int nTracks2PixelHits;
    int nTracks3PixelHits;
    int nTracks4PixelHits;
    int nTracks5PixelHits;
    int nTracksAtLeast6PixelHits;
    int nTracksValidHitInBPix1;
    int nTracks0LostInnerHits;
    int nTracks1LostInnerHit;
    int nTracksAtLeast2LostInnerHits;
    int nMatchedGenBquarks;
};



struct FatJetType {
FatJetType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), ptRaw(-1.), ptUnc(-1.), dPhi_met(-1.), dPhi_Jet1(-1.), puId(-1.), CSV(-99.), CSVR(-99.), CSVRUp(-99.), CSVRDown(-99.), pfBoostedDoubleSVAK8(-1.), CHSprunedMass(-1.), CHSsoftdropMass(-1.), softdropPuppiMass(-1.), CHSprunedMassCorr(-1.), CHSsoftdropMassCorr(-1.), softdropPuppiMassCorr(-1.), softdropPuppiMassCorrNotSmeared(-1.), nSoftDropSubJets(-1), pt1(-1.), eta1(-9.), phi1(-9.), mass1(-1.), CSV1(-99.), //CSVR1(-99.), CSVR1Up(-99.), CSVR1Down(-99.),
    CMVA1(-99.), //CMVAR1(-99.), CMVAR1Up(-99.), CMVAR1Down(-99.),
    flavour1(-1.), nSV1(-1), nVertexTracks1(-1), pt2(-1.), eta2(-9.), phi2(-9.), mass2(-1.), CSV2(-99.), //CSVR2(-99.), CSVR2Up(-99.), CSVR2Down(-99.),
    CMVA2(-99.), //CMVAR2(-99.), CMVAR2Up(-99.), CMVAR2Down(-99.),
    flavour2(-1.), nSV2(-1), nVertexTracks2(-1), dR(-1.), chsTau21(-1.), puppiTau21(-1.), BDSV(-1.),
//energies
cHadE(-1.), nHadE (-1.), eleE(-1.), photonE(-1.), muE(-1.), nEmE(-1.), cEmE(-1.), cmuE(-1.),
//energy fractions
cHadEFrac(-1.), nHadEFrac(-1.), eleEFrac(-1.), photonEFrac(-1.), muEFrac(-1.), nEmEFrac(-1.), cEmEFrac(-1.), cmuEFrac(-1.),
//multiplicities
cHadMulti(-1.), nHadMulti(-1.), eleMulti(-1.), photonMulti(-1.), muMulti(-1.), cMulti(-1.), nMulti(-1.), npr(-1.),
//multiplicity fractions
cHadMultiFrac(-1.), nHadMultiFrac(-1.), eleMultiFrac(-1.), photonMultiFrac(-1.), muMultiFrac(-1.), cMultiFrac(-1.), nMultiFrac(-1.),
partonFlavour(0), hadronFlavour(0), mother(0), isLoose(false), isMedium(false), isTight(false), isTightLepVeto(false), isMatched(false), JESUnc(-1.), ptJERUp(-1.), etaJERUp(-1.), phiJERUp(-9.), energyJERUp(-1.), ptJERDown(-1.), etaJERDown(-1.), phiJERDown(-9.), energyJERDown(-1.), smearFact(-1.), smearFactUp(-1.), smearFactDown(-1.), softdropPuppiMassCorrJMS(-1.), softdropPuppiMassCorrJMSUp(-1.), softdropPuppiMassCorrJMSDown(-1.), softdropPuppiMassCorrJMR(-1.), softdropPuppiMassCorrJMRUp(-1.), softdropPuppiMassCorrJMRDown(-1.), dR_q1(1000), dR_q2(1000), dR_q3(1000), dR_q4(1000), m_q1(false), m_q2(false), m_q3(false), m_q4(false), dR_pi1(1000), dR_pi2(1000), matchBquark(-1), matchLL(-1), isGenMatched(false),
//From PFCandidates
alphaMax(-100.), sigIP2DMedian(-100.), theta2DMedian(-100.), POCA_theta2DMedian(-100.), nPixelHitsMedian(-1.0), nHitsMedian(-1.0), nConstituents(-1), nTrackConstituents(-1),
nTracks0PixelHits(-1), nTracks1PixelHit(-1), nTracks2PixelHits(-1), nTracks3PixelHits(-1), nTracks4PixelHits(-1), nTracks5PixelHits(-1), nTracksAtLeast6PixelHits(-1),
nTracksValidHitInBPix1(-1), nTracks0LostInnerHits(-1), nTracks1LostInnerHit(-1), nTracksAtLeast2LostInnerHits(-1),
//
nMatchedGenBquarks(-1) {}
    float pt;
    float eta;
    float phi;
    float mass;
    float energy;
    float ptRaw;
    float ptUnc;
    float dPhi_met;
    float dPhi_Jet1;
    float puId;
    float CSV;
    float CSVR;
    float CSVRUp;
    float CSVRDown;
    float pfBoostedDoubleSVAK8;
    float CHSprunedMass;
    float CHSsoftdropMass;
    float softdropPuppiMass;
    float CHSprunedMassCorr;
    float CHSsoftdropMassCorr;
    float softdropPuppiMassCorr;
    float softdropPuppiMassCorrNotSmeared;
    int nSoftDropSubJets;
    float pt1;
    float eta1;
    float phi1;
    float mass1;
    float CSV1;
  //float CSVR1;
  //float CSVR1Up;
  //float CSVR1Down;
    float CMVA1;
  //float CMVAR1;
  //float CMVAR1Up;
  //float CMVAR1Down;
    float flavour1;
    int nSV1;
    int nVertexTracks1;
    float pt2;
    float eta2;
    float phi2;
    float mass2;
    float CSV2;
  //float CSVR2;
  //float CSVR2Up;
  //float CSVR2Down;
    float CMVA2;
  //float CMVAR2;
  //float CMVAR2Up;
  //float CMVAR2Down;
    float flavour2;
    int nSV2;
    int nVertexTracks2;
    float dR;
    float chsTau21;
    float puppiTau21;
    float ddtTau21;
    float BDSV;
    float cHadE;
    float nHadE;
    float eleE;
    float photonE;
    float muE;
    float nEmE;
    float cEmE;
    float cmuE;
    float cHadEFrac;
    float nHadEFrac;
    float eleEFrac;
    float photonEFrac;
    float muEFrac;
    float nEmEFrac;
    float cEmEFrac;
    float cmuEFrac;
    float cHadMulti;
    float nHadMulti;
    float eleMulti;
    float photonMulti;
    float muMulti;
    float cMulti;
    float nMulti;
    float npr;
    float cHadMultiFrac;
    float nHadMultiFrac;
    float eleMultiFrac;
    float photonMultiFrac;
    float muMultiFrac;
    float cMultiFrac;
    float nMultiFrac;
    int partonFlavour;
    int hadronFlavour;
    int mother;
    bool isLoose;
    bool isMedium;
    bool isTight;
    bool isTightLepVeto;
    bool isMatched;
    float JESUnc;
    float ptJERUp;
    float etaJERUp;
    float phiJERUp;
    float energyJERUp;
    float ptJERDown;
    float etaJERDown;
    float phiJERDown;
    float energyJERDown;
    float smearFact;
    float smearFactUp;
    float smearFactDown;
    float softdropPuppiMassCorrJMS;
    float softdropPuppiMassCorrJMSUp;
    float softdropPuppiMassCorrJMSDown;
    float softdropPuppiMassCorrJMR;
    float softdropPuppiMassCorrJMRUp;
    float softdropPuppiMassCorrJMRDown;
    float dR_q1;
    float dR_q2;
    float dR_q3;
    float dR_q4;
    bool m_q1;
    bool m_q2;
    bool m_q3;
    bool m_q4;
    float dR_pi1;
    float dR_pi2;
    int matchBquark;
    int matchLL;
    bool isGenMatched;
    float alphaMax;
    float sigIP2DMedian;
    float theta2DMedian;
    float POCA_theta2DMedian;
    float nPixelHitsMedian;
    float nHitsMedian;
    int nConstituents;
    int nTrackConstituents;
    int nTracks0PixelHits;
    int nTracks1PixelHit;
    int nTracks2PixelHits;
    int nTracks3PixelHits;
    int nTracks4PixelHits;
    int nTracks5PixelHits;
    int nTracksAtLeast6PixelHits;
    int nTracksValidHitInBPix1;
    int nTracks0LostInnerHits;
    int nTracks1LostInnerHit;
    int nTracksAtLeast2LostInnerHits;
    int nMatchedGenBquarks;
};


//struct MEtType {
//    MEtType(): pt(-1.), eta(-9.), phi(-9.), sign(-1.) {}
//    float pt;
//    float eta;
//    float phi;
//    float sign;
//};

struct MEtType {
  //MEtType(): pt(-1.), eta(-9.), phi(-9.), sign(-1.), ptRaw(-1.), phiRaw(-9.), ptType1(-1.), phiType1(-9.), ptGen(-1.), phiGen(-9.), ptScaleUp(-1.), ptScaleDown(-1.), ptResUp(-1.), ptResDown(-1.), ptCalo(-1.) {}
    MEtType(): pt(-1.), eta(-9.), phi(-9.), sign(-1.), ptShiftJetResUp(-1.), ptShiftJetResDown(-1.), ptShiftJetEnUp(-1.), ptShiftJetEnDown(-1.), ptShiftUnclusteredEnUp(-1.), ptShiftUnclusteredEnDown(-1.), ptShiftJetResUpSmear(-1.), ptShiftJetResDownSmear(-1.), ptRaw(-1.), phiRaw(-9.), ptGen(-1.), phiGen(-9.), ptCalo(-1.) {}
    float pt;
    float eta;
    float phi;
    float sign;
    float ptShiftJetResUp;
    float ptShiftJetResDown;
    float ptShiftJetEnUp;
    float ptShiftJetEnDown;
    float ptShiftUnclusteredEnUp;
    float ptShiftUnclusteredEnDown;
    float ptShiftJetResUpSmear;
    float ptShiftJetResDownSmear;
    float ptRaw;
    float phiRaw;
  //float ptType1;
  //float phiType1;
    float ptGen;
    float phiGen;
  //float ptScaleUp;
  //float ptScaleDown;
  //float ptResUp;
  //float ptResDown;
    float ptCalo;
};

struct MEtFullType {
    MEtFullType(): pt(-1.), eta(-9.), phi(-9.), sign(-1.), ptRaw(-1.), phiRaw(-9.), ptGen(-1.), phiGen(-9.), ptJERUp(-1.), ptJERDown(-1.), ptJERUpSmear(-1.), ptJERDownSmear(-1.), ptJESUp(-1.), ptJESDown(-1.), ptMUSUp(-1.), ptMUSDown(-1.), ptELSUp(-1.), ptELSDown(-1.), ptTAUUp(-1.), ptTAUDown(-1.), ptUNCUp(-1.), ptUNCDown(-1.), ptPHOUp(-1.), ptPHODown(-1.), phf(-1.), nhf(-1.), elf(-1.), chf(-1.), muf(-1.) {}
    float pt;
    float eta;
    float phi;
    float sign;
    float ptRaw;
    float phiRaw;
    float ptGen;
    float phiGen;
    float ptJERUp;
    float ptJERDown;
    float ptJERUpSmear;
    float ptJERDownSmear;
    float ptJESUp;
    float ptJESDown;
    float ptMUSUp;
    float ptMUSDown;
    float ptELSUp;
    float ptELSDown;
    float ptTAUUp;
    float ptTAUDown;
    float ptUNCUp;
    float ptUNCDown;
    float ptPHOUp;
    float ptPHODown;
    float phf;
    float nhf;
    float elf;
    float chf;
    float muf;
};


struct CandidateType {
    CandidateType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), tmass(-1.), dR(-1.), dEta(-1.), dPhi(-1.), twist(-1.) {} // angle(-1.), ptBalance(-9.), centrality(-1.), charge(0),
    float pt;
    float eta;
    float phi;
    float mass;
    float tmass;
//    float softdropMass;
//    float tmassScaleUp;
//    float tmassScaleDown;
//    float tmassResUp;
//    float tmassResDown;
    float dR;
    float dEta;
    float dPhi;
    float twist;
  //  float angle;
  //  float ptBalance;
  //  float centrality;
  //  int charge;
//    bool isMatched;
};





struct LorentzType {
  LorentzType(): pt(-1.), eta(-9.), phi(-9.), energy(-1.), mass(-1.) {}
  float pt;
  float eta;
  float phi;
  float energy;
  float mass;
};

struct EventType {
  EventType(): id1(-9), id2(-9), x1(-1.), x2(-1.), Q(-1.) {}
  int id1;
  int id2;
  float x1;
  float x2;
  float Q;
};

struct GenPType {
GenPType(): pt(-1.), eta(-9.), rapidity(-99999.), phi(-9.), mass(-1.), energy(-1.), charge(0), pdgId(0), status(0), radius(-1), motherid(0) {}
    float pt;
    float eta;
    float rapidity;
    float phi;
    float mass;
    float energy;
    int charge;
    int pdgId;
    int status;
    float radius;
    int motherid;
    float vx;
    float vy;
    float vz;
};

struct TriggerObjectType {
TriggerObjectType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), charge(0) {}
    float pt;
    float eta;
    float phi;
    float mass;
    float energy;
    int charge;
};


struct CustomFatJetType {
//CustomFatJetType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), ptRaw(-1.), ptUnc(-1.), dPhi_met(-1.), dPhi_Jet1(-1.), puId(-1.), CSV(-99.), CSVR(-99.), CSVRUp(-99.), CSVRDown(-99.), CHSprunedMass(-1.), CHSsoftdropMass(-1.), prunedMass(-1.), softdropMass(-1.), softdropPuppiMass(-1.), CHSprunedMassCorr(-1.), CHSsoftdropMassCorr(-1.), prunedMassCorr(-1.), softdropMassCorr(-1.), softdropPuppiMassCorr(-1.), softdropPuppiMassCorrNotSmeared(-1.), pt1(-1.), eta1(-9.), phi1(-9.), mass1(-1.), CSV1(-99.), CSVR1(-99.), CSVR1Up(-99.), CSVR1Down(-99.), CMVA1(-99.), CMVAR1(-99.), CMVAR1Up(-99.), CMVAR1Down(-99.), flavour1(-1.), pt2(-1.), eta2(-9.), phi2(-9.), mass2(-1.), CSV2(-99.), CSVR2(-99.), CSVR2Up(-99.), CSVR2Down(-99.), CMVA2(-99.), CMVAR2(-99.), CMVAR2Up(-99.), CMVAR2Down(-99.), flavour2(-1.), dR(-1.), Tau21(-1.), puppiTau21(-1.), BDSV(-1.), chf(-1.), nhf(-1.), phf(-1.), elf(-1.), muf(-1.), chm(-1), npr(-1), flavour(0), mother(0), isLoose(false), isMedium(false), isTight(false), isTightLepVeto(false), isMatched(false), JESUnc(-1.), ptJERUp(-1.), etaJERUp(-1.), phiJERUp(-9.), energyJERUp(-1.), ptJERDown(-1.), etaJERDown(-1.), phiJERDown(-9.), energyJERDown(-1.), smearFact(-1.), smearFactUp(-1.), smearFactDown(-1.), softdropPuppiMassCorrJMS(-1.), softdropPuppiMassCorrJMSUp(-1.), softdropPuppiMassCorrJMSDown(-1.), softdropPuppiMassCorrJMR(-1.), softdropPuppiMassCorrJMRUp(-1.), softdropPuppiMassCorrJMRDown(-1.), dR_q1(1000), dR_q2(1000), dR_q3(1000), dR_q4(1000), m_q1(false), m_q2(false), m_q3(false), m_q4(false), dR_pi1(1000), dR_pi2(1000), matchBquark(-1), matchLL(-1) {}
//simplified version
CustomFatJetType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), ptRaw(-1.), ptUnc(-1.), dPhi_met(-1.), dPhi_Jet1(-1.), puId(-1.), CSV(-99.), CHSprunedMass(-1.), CHSsoftdropMass(-1.), prunedMass(-1.), softdropMass(-1.), softdropPuppiMass(-1.), CHSprunedMassCorr(-1.), CHSsoftdropMassCorr(-1.), prunedMassCorr(-1.), softdropMassCorr(-1.), softdropPuppiMassCorr(-1.), softdropPuppiMassCorrNotSmeared(-1.), pt1(-1.), eta1(-9.), phi1(-9.), mass1(-1.), CSV1(-99.), CMVA1(-99.), flavour1(-1.), pt2(-1.), eta2(-9.), phi2(-9.), mass2(-1.), CSV2(-99.), CMVA2(-99.), flavour2(-1.), pt1SDP(-1.), eta1SDP(-9.), phi1SDP(-9.), mass1SDP(-1.), CSV1SDP(-99.), CMVA1SDP(-99.), flavour1SDP(-1.), pt2SDP(-1.), eta2SDP(-9.), phi2SDP(-9.), mass2SDP(-1.), CSV2SDP(-99.), CMVA2SDP(-99.), flavour2SDP(-1.), dR(-1.), Tau21(-1.), puppiTau21(-1.), BDSV(-1.), chf(-1.), nhf(-1.), phf(-1.), elf(-1.), muf(-1.), chm(-1), npr(-1), partonFlavour(0), hadronFlavour(0), mother(0), isLoose(false), isMedium(false), isTight(false), isTightLepVeto(false), dR_q1(1000), dR_q2(1000), dR_q3(1000), dR_q4(1000), m_q1(false), m_q2(false), m_q3(false), m_q4(false), dR_pi1(1000), dR_pi2(1000), matchBquark(-1), matchLL(-1), dR_q1_sj1(1000), dR_q2_sj1(1000), dR_q3_sj1(1000), dR_q4_sj1(1000), dR_pi1_sj1(1000), dR_pi2_sj1(1000), dR_q1_sj2(1000), dR_q2_sj2(1000), dR_q3_sj2(1000), dR_q4_sj2(1000), dR_pi1_sj2(1000), dR_pi2_sj2(1000) {}
    float pt;
    float eta;
    float phi;
    float mass;
    float energy;
    float ptRaw;
    float ptUnc;
    float dPhi_met;
    float dPhi_Jet1;
    float puId;
    float CSV;
  //float CSVR;
  //float CSVRUp;
  //float CSVRDown;
    float CHSprunedMass;
    float CHSsoftdropMass;
    float prunedMass;
    float softdropMass;
    float softdropPuppiMass;
    float CHSprunedMassCorr;
    float CHSsoftdropMassCorr;
    float prunedMassCorr;
    float softdropMassCorr;
    float softdropPuppiMassCorr;
    float softdropPuppiMassCorrNotSmeared;
    float pt1;
    float eta1;
    float phi1;
    float mass1;
    float CSV1;
  //float CSVR1;
  //float CSVR1Up;
  //float CSVR1Down;
    float CMVA1;
  //float CMVAR1;
  //float CMVAR1Up;
  //float CMVAR1Down;
    float flavour1;
    float pt2;
    float eta2;
    float phi2;
    float mass2;
    float CSV2;
  //float CSVR2;
  //float CSVR2Up;
  //float CSVR2Down;
    float CMVA2;
  //float CMVAR2;
  //float CMVAR2Up;
  //float CMVAR2Down;
    float flavour2;
    float pt1SDP;
    float eta1SDP;
    float phi1SDP;
    float mass1SDP;
    float CSV1SDP;
    float CMVA1SDP;
    float flavour1SDP;
    float pt2SDP;
    float eta2SDP;
    float phi2SDP;
    float mass2SDP;
    float CSV2SDP;
    float CMVA2SDP;
    float flavour2SDP;
    float dR;
    float Tau21;
    float puppiTau21;
    float ddtTau21;
    float BDSV;
    float chf;
    float nhf;
    float phf;
    float elf;
    float muf;
    int chm;
    int npr;
    int partonFlavour;
    int hadronFlavour;
    int mother;
    bool isLoose;
    bool isMedium;
    bool isTight;
    bool isTightLepVeto;
  //bool isMatched;
  //float JESUnc;
  //float ptJERUp;
  //float etaJERUp;
  //float phiJERUp;
  //float energyJERUp;
  //float ptJERDown;
  //float etaJERDown;
  //float phiJERDown;
  //float energyJERDown;
  //float smearFact;
  //float smearFactUp;
  //float smearFactDown;
  //float softdropPuppiMassCorrJMS;
  //float softdropPuppiMassCorrJMSUp;
  //float softdropPuppiMassCorrJMSDown;
  //float softdropPuppiMassCorrJMR;
  //float softdropPuppiMassCorrJMRUp;
  //float softdropPuppiMassCorrJMRDown;
    float dR_q1;
    float dR_q2;
    float dR_q3;
    float dR_q4;
    bool m_q1;
    bool m_q2;
    bool m_q3;
    bool m_q4;
    float dR_pi1;
    float dR_pi2;
    int matchBquark;
    int matchLL;
    float dR_q1_sj1;
    float dR_q2_sj1;
    float dR_q3_sj1;
    float dR_q4_sj1;
    float dR_pi1_sj1;
    float dR_pi2_sj1;
    float dR_q1_sj2;
    float dR_q2_sj2;
    float dR_q3_sj2;
    float dR_q4_sj2;
    float dR_pi1_sj2;
    float dR_pi2_sj2;
};



struct CaloJetType {
CaloJetType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), emEnergyFraction(-1.), emEnergyInEB(-1.), emEnergyInEE(-1.), emEnergyInHF(-1.), energyFractionHadronic(-1.), hadEnergyInHB(-1.), hadEnergyInHE(-1.), hadEnergyInHF(-1.), hadEnergyInHO(-1.), longLived(false), maxEInEmTowers(-1.), maxEInHadTowers(-1.), n60(-1), n90(-1), nConstituents(-1), nPasses(-1), isGenMatched(false)  {}

    float pt;
    float eta;
    float phi;
    float mass;
    float energy;
    float emEnergyFraction;
    float emEnergyInEB;
    float emEnergyInEE;
    float emEnergyInHF;
    float energyFractionHadronic;
    float hadEnergyInHB;
    float hadEnergyInHE;
    float hadEnergyInHF;
    float hadEnergyInHO;
    bool longLived;
    float maxEInEmTowers;
    float maxEInHadTowers;
    int n60;
    int n90;
    int nConstituents;
    int nPasses;
    bool isGenMatched;
};


struct PFCandidateType {
PFCandidateType():
    charge(-99), pt(-1.), eta(-9.), phi(-9.),
    energy(-1.), mass(-1.), px(-99.), py(-99.), pz(-99.),
    pdgId(0), isTrack(false), jetIndex (-1), fatJetIndex(-1), pvIndex(-1),
    dxy(-99.), dz(-99.), POCA_x(-99.), POCA_y(-99.), POCA_z(-99.), POCA_phi(-9.),
    ptError(-1.), etaError(-1.), phiError(-1.), dxyError(-99.), dzError(-99.), theta(-9.), thetaError(-1.),
    chi2(-1.), ndof(-1), normalizedChi2(-1.), nHits(-1), nPixelHits(-1), lostInnerHits(-9) {}
    //hcalFraction(-1.), longLived(-1)
    //innerDetId(-1), innerPosition_x(-9999.), innerPosition_y(-9999.), innerPosition_z(-9999.),
    //innerMomentum_x(-9999.), innerMomentum_y(-9999.), innerMomentum_z(-9999.),

    // General
    int   charge;
    float pt;
    float eta;
    float phi;

    float energy;
    float mass;
    float px;
    float py;
    float pz;

    int pdgId;
    bool isTrack;
    int jetIndex;
    int fatJetIndex;
    int pvIndex;

    // Tracking:
    float dxy;              // Transverse impact parameter
    float dz;               // Longitudinal impact parameter
    float POCA_x;           // Point Of Closest Approach coordinates
    float POCA_y;
    float POCA_z;
    float POCA_phi;         // Momentum direction (from tracking-only information)

    // Additional tracking information (for pT > 0.95 GeV)
    float ptError;
    float etaError;
    float phiError;
    float dxyError;
    float dzError;
    float theta;
    float thetaError;

    float chi2;
    int ndof;
    float normalizedChi2;

    int   nHits;
    int   nPixelHits;
    int   lostInnerHits;

    // Not used:
    //float hcalFraction;
    //int longLived;

    // Not stored in miniAOD:
    //int innerDetId;         // Detector with the innermost hit
    //float innerPosition_x;  // Position of the innermost hit
    //float innerPosition_y;
    //float innerPosition_z;
    //float innerMomentum_x;  // Momentum direction at the innermost hit
    //float innerMomentum_y;
    //float innerMomentum_z;
};


struct VertexType {
VertexType(): chi2(-1.), ndof(-1.), x(-100.), y(-100.), z(-100.), px(-100.), py(-100.), pz(-100.), mass(-1.), isValid(false), isFake(false), pt(-1.), flightDist2d(-100.), flightDist2dErr(-100.), flightDist2dSig(-100.), flightDist3d(-100.), flightDist3dErr(-100.), flightDist3dSig(-100.), nVertexTracks(-1), dRSVJet(-100.), jetIndex(-1), vertexMass(-1.), trackEta(-100.), trackMomentum(-100.), trackPhi(-100.), trackPtRel(-100.), trackPPar(-100.), trackEtaRel(-100.), trackDeltaR(-100.), trackPtRatio(-100.), trackPParRatio(-100.), trackSip2dVal(-100.), trackSip2dSig(-100.), trackSip3dVal(-100.), trackSip3dSig(-100.), trackDecayLenVal(-100.), trackDecayLenSig(-100.), trackJetDistVal(-100.), trackJetDistSig(-100.), trackChi2(-100.), trackNTotalHits(-1), trackNPixelHits(-1), trackJetPt(-100.)  {}
  //, trackSip3dSig_3(-100.), trackSip3dSig_2(-100.), trackSip3dSig_1(-100.), trackSip3dSig_0(-100.), trackSip2dSigAboveCharm_0(-100.), trackSip2dSigAboveBottom_0(-100.) , trackSip2dSigAboveBottom_1(-100.)
  float chi2;
  float ndof;
  float x;
  float y;
  float z;
  float px;
  float py;
  float pz;
  float mass;
  bool isValid;
  bool isFake;
  float pt;
  float flightDist2d;
  float flightDist2dErr;
  float flightDist2dSig;
  float flightDist3d;
  float flightDist3dErr;
  float flightDist3dSig;
  int nVertexTracks;
  float dRSVJet;
  int jetIndex;
  float vertexMass;
  float trackEta;
  float trackMomentum;
  float trackPhi;
  float trackPtRel;
  float trackPPar;
  float trackEtaRel;
  float trackDeltaR;
  float trackPtRatio;
  float trackPParRatio;
  float trackSip2dVal;
  float trackSip2dSig;
  float trackSip3dVal;
  float trackSip3dSig;
  float trackDecayLenVal;
  float trackDecayLenSig;
  float trackJetDistVal;
  float trackJetDistSig;
  float trackChi2;
  int trackNTotalHits;
  int trackNPixelHits;
  float trackJetPt;
  /* float trackSip3dSig_3; */
  /* float trackSip3dSig_2; */
  /* float trackSip3dSig_1; */
  /* float trackSip3dSig_0; */
  /* float trackSip2dSigAboveCharm_0; */
  /* float trackSip2dSigAboveBottom_0; */
  /* float trackSip2dSigAboveBottom_1; */

};


struct SimplifiedJetType {
SimplifiedJetType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.),cHadE(-1.), nHadE (-1.), cHadEFrac(-1.), nHadEFrac(-1.), nEmE(-1.), nEmEFrac(-1.), cEmE(-1.), cEmEFrac(-1.), cmuE(-1.), cmuEFrac(-1.), muE(-1.), muEFrac(-1.), eleE(-1.), eleEFrac(-1.), eleMulti(-1.), photonE(-1.), photonEFrac(-1.), photonMulti(-1.), cHadMulti(-1.), nHadMulti(-1.), npr(-1.), cMulti(-1.), nMulti(-1.), isLoose(false), isMedium(false), isTight(false), isGenMatched(0), nSelectedTracks(0), nConstituents (-1), TriggerMatched_VBFJet(0), TriggerMatched_DisplacedJet(0), TriggerMatched_TripleJet50(0){}
    float pt;
    float eta;
    float phi;
    float mass;
    float energy;
    float cHadE;
    float nHadE;
    float cHadEFrac;
    float nHadEFrac;
    float nEmE;
    float nEmEFrac;
    float cEmE;
    float cEmEFrac;
    float cmuE;
    float cmuEFrac;
    float muE;
    float muEFrac;
    float eleE;
    float eleEFrac;
    float eleMulti;
    float photonE;
    float photonEFrac;
    float photonMulti;
    float cHadMulti;
    float nHadMulti;
    float npr;
    float cMulti;
    float nMulti;
    bool isLoose;
    bool isMedium;
    bool isTight;
    int isGenMatched;
    int nSelectedTracks;
    int nConstituents;
    int TriggerMatched_VBFJet;
    int TriggerMatched_DisplacedJet;
    int TriggerMatched_TripleJet50;

};


#endif
