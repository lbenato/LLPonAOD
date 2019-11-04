#ifndef OBJECTS_H
#define OBJECTS_H

struct LeptonType {
LeptonType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), inTrkPt(-1.), pfIso03(-1.), pfIso04(-1.), trkIso(-1.), miniIso(-1.), dxy(-99.), dz(-99.), ip3d(-99.), sip3d(-99.), nPixelHits(-1.), dPhi_met(-1.), charge(0), pdgId(0), isElectron(false), isMuon(false), isVeto(false), isLoose(false), isMedium(false), isTight(false), isHighPt(false), isTrackerHighPt(false), SSscale(-1.), SSsigma(-1.), SSscaleUnc(-1.), SSsigmaUncUp(-1.), SSsigmaUncDown(-1.), SScorr(-1.), energySScorr(-1.), energySScorrUncUp(-1.), energySScorrUncDown(-1.), ptSScorr(-1.), ptSScorrUncUp(-1.), ptSScorrUncDown(-1.), isMatched(false) {}
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
JetType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.) {}

  float pt;
  float eta;
  float phi;
  float mass;
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

struct MEtType {
    MEtType(): pt(-1.), eta(-9.), phi(-9.), sign(-1.), ptCalo(-1.) {}
    float pt;
    float eta;
    float phi;
    float sign;
    float ptCalo;
};

struct CandidateType {
    CandidateType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), tmass(-1.), dR(-1.), dEta(-1.), dPhi(-1.), twist(-1.) {} // angle(-1.), ptBalance(-9.), centrality(-1.), charge(0),
    float pt;
    float eta;
    float phi;
    float mass;
    float tmass;
    float dR;
    float dEta;
    float dPhi;
    float twist;
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


#endif
