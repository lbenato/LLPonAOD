#ifndef ELECTRONANALYZER_H
#define ELECTRONANALYZER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "TFile.h"
#include "TH2.h"
#include "TRandom3.h"

#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"

class EnergyScaleCorrection_class;

class ElectronAnalyzer {
    public:
        ElectronAnalyzer(const edm::ParameterSet&, edm::ConsumesCollector&&);
        ~ElectronAnalyzer();
        virtual std::vector<pat::Electron> FillElectronVector(const edm::Event&);
        virtual void AddVariables(std::vector<pat::Electron>&, pat::MET&);
        virtual float GetDoubleElectronTriggerSF(pat::Electron&, pat::Electron&);
        //virtual float GetLooseElectronSF(pat::Electron&);
        virtual float GetElectronIdSF(pat::Electron&, int);
        virtual float GetElectronIdSFError(pat::Electron&, int);
        virtual float GetElectronRecoEffSF(pat::Electron&);
        virtual float GetElectronRecoEffSFError(pat::Electron&);
        virtual float GetElectronTriggerSFEle105(pat::Electron&);
        virtual float GetElectronTriggerSFErrorEle105(pat::Electron&);
        virtual float GetElectronTriggerSFEle27Tight(pat::Electron&);
	virtual float GetElectronTriggerSFErrorEle27Tight(pat::Electron&);
        
      
    private:
      
        edm::EDGetTokenT<std::vector<pat::Electron> > ElectronToken;
        edm::EDGetTokenT<reco::VertexCollection> VertexToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> EleVetoIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> EleLooseIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> EleMediumIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> EleTightIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> EleHEEPIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> EleMVANonTrigMediumIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> EleMVANonTrigTightIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> EleMVATrigMediumIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> EleMVATrigTightIdMapToken;
        edm::EDGetTokenT<EcalRecHitCollection> EleEcalRecHitCollectionToken;
        std::string EleSingleTriggerIsoFileName;
        std::string EleSingleTriggerFileName;
        std::string EleVetoIdFileName;
        std::string EleLooseIdFileName;
        std::string EleMediumIdFileName;
        std::string EleTightIdFileName;
        std::string EleMVATrigMediumIdFileName;
        std::string EleMVATrigTightIdFileName;
        std::string EleRecoEffFileName;
        std::string EleScaleSmearCorrectionName;

        int Electron1Id, Electron2Id;// Electron1Iso, Electron2Iso;
        float Electron1Pt, Electron2Pt;
        float EleTriggerPtMax;
        
        bool isEleVetoIdFile, isEleLooseIdFile, isEleMediumIdFile, isEleTightIdFile, isEleMVATrigMediumIdFile, isEleMVATrigTightIdFile, isEleTriggerFile, isEleTriggerIsoFile, isEleSingleTriggerFile, isEleSingleTriggerIsoFile, isEleRecoEffFile;
        
        TFile* EleTriggerIsoFile;
        TFile* EleTriggerFile;
        TFile* EleSingleTriggerIsoFile;
        TFile* EleSingleTriggerFile;
        TFile* EleVetoIdFile;
        TFile* EleLooseIdFile;
        TFile* EleMediumIdFile;
        TFile* EleTightIdFile;
        TFile* EleMVATrigMediumIdFile;
        TFile* EleMVATrigTightIdFile;
        TFile* EleRecoEffFile;
        
        TH2F* EleTriggerDATAHighLeg;
        TH2F* EleTriggerDATALowLeg;
        TH2F* EleTriggerMCHighLeg;
        TH2F* EleTriggerMCLowLeg;
        TH2F* ElectronTriggerEle105;
        TH2F* ElectronTriggerEle27Tight;
        TH2F* ElectronTriggerErrorEle27Tight;
        TH2F* ElectronIdVeto;
        TH2F* ElectronIdLoose;
        TH2F* ElectronIdMedium;
        TH2F* ElectronIdTight;
        TH2F* ElectronIdMVATrigMedium;
        TH2F* ElectronIdMVATrigTight;
        TH2F* ElectronRecoEff;

        EnergyScaleCorrection_class * eScaleSmearer;

};


#endif
