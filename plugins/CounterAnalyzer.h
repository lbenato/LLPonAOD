// -*- C++ -*-
//
// Package:    Analysis/CounterAnalyzer
// Class:      CounterAnalyzer
// 
/**\class CounterAnalyzer CounterAnalyzer.cc Analysis/CounterAnalyzer/plugins/CounterAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alberto Zucchetta
//         Created:  Tue, 31 May 2016 12:47:56 GMT
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"

#include "TH1.h"
#include "TH3.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class CounterAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit CounterAnalyzer(const edm::ParameterSet&);
      ~CounterAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
      //edm::EDGetToken LheToken;
      edm::EDGetToken GenToken;
      edm::Service<TFileService> fs;
      TH1F* Weight;
      //TH1F* NPartons;
      //TH1F* NBPartons;
      //TH1F* LheHT;
      //TH1F* LhePtZ;
      //TH3F* Bin;
      bool PythiaLOSample;
};

