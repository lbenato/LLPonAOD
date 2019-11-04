#ifndef GENANALYZER_H
#define GENANALYZER_H

#include <iostream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TFile.h"
#include "TH3.h"
#include "TF1.h"
#include "TKey.h"

class GenAnalyzer {
    public:
        GenAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~GenAnalyzer();
        virtual std::map<int, float> FillWeightsMap(const edm::Event&);
        virtual std::map<std::string, float> FillLheMap(const edm::Event&);
        virtual std::vector<reco::GenParticle> FillGenVector(const edm::Event&);
        virtual std::vector<reco::GenParticle> FillGenVectorByIdAndStatus(const edm::Event&, int, int);
        virtual std::vector<reco::GenParticle> FillGenVectorByIdListAndStatusAndMotherList(const edm::Event&, std::vector<int>, int, std::vector<int>);
        virtual std::vector<reco::GenParticle> FillVBFGenVector(const edm::Event&);//, std::vector<int>, int, std::vector<int>, std::vector<int>);
        virtual std::vector<reco::GenParticle> FillGenVectorByIdStatusAndMother(const edm::Event&, int, int, int);

        virtual std::vector<reco::GenParticle> FillGenVectorByIdAndStatusAndKin(const edm::Event&, int, int, float, float);
        virtual std::vector<reco::GenParticle> FillGenVectorByIdStatusAndMotherAndKin(const edm::Event&, int, int, int, float, float);

        virtual reco::Candidate* FindGenParticle(std::vector<reco::GenParticle>&, int);
        virtual reco::Candidate* FindLastDaughter(reco::Candidate*);
        virtual reco::GenParticle* FindGenParticleGenByIds(std::vector<reco::GenParticle>&, std::vector<int>, int=-1);
        virtual reco::GenParticle* FindLastDaughterGen(reco::GenParticle*);
        virtual const reco::GenParticle* FindLastDaughterGen(const reco::GenParticle*);
        virtual const reco::Candidate* FindMother(reco::GenParticle*);
        virtual reco::Candidate* FindGenParticleByIdAndStatus(std::vector<reco::GenParticle>&, int, int);
        virtual float GetStitchWeight(std::map<std::string, float>);
        virtual float GetZewkWeight(float);
        virtual float GetWewkWeight(float);
        virtual float GetTopPtWeight(float );
        virtual float GetPUWeight(const edm::Event&);
    //    virtual float GetPDFWeight(const edm::Event&);
        virtual std::pair<float, float> GetQ2Weight(const edm::Event&);
        virtual std::vector<reco::GenParticle> PartonsFromDecays(const std::vector<int> & pdgIds);
        virtual std::vector<reco::GenParticle> PartonsFromDecays(const std::vector<int> & pdgIds, std::vector<reco::GenParticle> & genDecay );
        virtual std::vector<reco::GenParticle> FirstNGenParticles(const std::vector<int> & pdgIds, std::size_t n); 

      
    private:
        edm::EDGetTokenT<GenEventInfoProduct> GenToken;
        edm::EDGetTokenT<LHEEventProduct> LheToken;
        edm::EDGetTokenT<std::vector<reco::GenParticle> > GenParticlesToken;
        std::vector<int> ParticleList;
        std::vector<int> ParticleStatus;
        std::vector<std::string> SampleDYJetsToLL;
        std::vector<std::string> SampleZJetsToNuNu;
        std::vector<std::string> SampleWJetsToLNu;
        std::string SampleDir;
        std::string Sample;
        
        std::map<std::string, TFile*> Files;
        std::map<std::string, TH1F*> hPartons;
        std::map<std::string, TH1F*> hBPartons;
        std::map<std::string, TH1F*> hHT;
        std::map<std::string, TH1F*> hPtV;
        
        std::string EWKFileName;
        bool ApplyEWK;
        bool ApplyTopPtReweigth;
        bool PythiaLOSample;
        bool isRealData;

        
        TFile* EWKFile;
        TF1* fZEWK;
        TF1* fWEWK;
        
        edm::LumiReWeighting* LumiWeights;

        // new member to hold GenCollection and about copy overheads
        edm::Handle<std::vector<reco::GenParticle> > GenCollection;

};


//namespace LHAPDF {
//  void initPDFSet(int nset, int setid, int member=0);
//  void initPDFSet(int nset, const std::string& filename, int member=0);
//  int numberPDF(int nset);
//  void usePDFMember(int nset, int member);
//  double xfx(int nset, double x, double Q, int fl);
//  double getXmin(int nset, int member);
//  double getXmax(int nset, int member);
//  double getQ2min(int nset, int member);
//  double getQ2max(int nset, int member);
//  void extrapolate(bool extrapolate=true);
//}


#endif
