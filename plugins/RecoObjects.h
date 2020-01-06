#ifndef RECOOBJECTS_H
#define RECOOBJECTS_H

struct RecoLeptonType {
RecoLeptonType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.) {}

  float pt;
  float eta;
  float phi;
  float mass;
};


struct RecoJetType {
RecoJetType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.) {}

  float pt;
  float eta;
  float phi;
  float mass;
};
/*
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
*/
struct RecoMEtType {
    RecoMEtType(): pt(-1.), eta(-9.), phi(-9.) {}
    float pt;
    float eta;
    float phi;
};

/*
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
*/

#endif
