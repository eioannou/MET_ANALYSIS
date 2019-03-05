// -*- C++ -*-
//
// Package:    MET_ANALYSIS/MissingEtAnalyzer
// Class:      METUncertainties
// 
/**\class METUncertainties METUncertainties.cc MET_ANALYSIS/MissingEtAnalyzer/plugins/METUncertainties.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Aimilios Ioannou
//         Created:  Mon, 11 Feb 2019 16:09:33 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <TROOT.h>
#include "TFile.h"
#include "TH1.h"
#include "TLorentzVector.h"


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "Math/GenVector/LorentzVector.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class METUncertainties : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit METUncertainties(const edm::ParameterSet&);
      ~METUncertainties();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  edm::InputTag metSrcTag_;
  edm::InputTag metmodifiedSrcTag_;
  edm::InputTag metPuppiSrcTag_;
  edm::InputTag metPuppiRecorrectSrcTag_;
  edm::EDGetTokenT<edm::View<pat::MET> > metSrcToken_;
  edm::EDGetTokenT<edm::View<pat::MET> > metmodifiedSrcToken_;
  edm::EDGetTokenT<edm::View<pat::MET> > metPuppiSrcToken_;
  edm::EDGetTokenT<edm::View<pat::MET> > metPuppiRecorrectSrcToken_;

  int nEvent;

  // MET belongs to MiniAOD
  TH1F *hMET_pt, *hMET_phi;
  TH1F *hMET_pt_jecup,  *hMET_phi_jecup;
  TH1F *hMET_pt_jecdn,  *hMET_phi_jecdn;
  TH1F *hMET_pt_muonup, *hMET_phi_muonup;
  TH1F *hMET_pt_muondn, *hMET_phi_muondn;
  TH1F *hMET_pt_eleup,  *hMET_phi_eleup;
  TH1F *hMET_pt_eledn,  *hMET_phi_eledn;
  TH1F *hMET_pt_tauup,  *hMET_phi_tauup;
  TH1F *hMET_pt_taudn,  *hMET_phi_taudn;
  TH1F *hMET_pt_uncup,  *hMET_phi_uncup;
  TH1F *hMET_pt_uncdn,  *hMET_phi_uncdn;
  // MET belongs to Modified MiniAOD from Met tool
  TH1F *hMETmod_pt, *hMETmod_phi;
  TH1F *hMETmod_pt_jecup,  *hMETmod_phi_jecup;
  TH1F *hMETmod_pt_jecdn,  *hMETmod_phi_jecdn;
  TH1F *hMETmod_pt_muonup, *hMETmod_phi_muonup;
  TH1F *hMETmod_pt_muondn, *hMETmod_phi_muondn;
  TH1F *hMETmod_pt_eleup,  *hMETmod_phi_eleup;
  TH1F *hMETmod_pt_eledn,  *hMETmod_phi_eledn;
  TH1F *hMETmod_pt_tauup,  *hMETmod_phi_tauup;
  TH1F *hMETmod_pt_taudn,  *hMETmod_phi_taudn;
  TH1F *hMETmod_pt_uncup,  *hMETmod_phi_uncup;
  TH1F *hMETmod_pt_uncdn,  *hMETmod_phi_uncdn;
  // MET Puppi
  TH1F *hMETPuppi_pt, *hMETPuppi_phi;
  TH1F *hMETPuppi_pt_jecup,  *hMETPuppi_phi_jecup;
  TH1F *hMETPuppi_pt_jecdn,  *hMETPuppi_phi_jecdn;
  TH1F *hMETPuppi_pt_muonup, *hMETPuppi_phi_muonup;
  TH1F *hMETPuppi_pt_muondn, *hMETPuppi_phi_muondn;
  TH1F *hMETPuppi_pt_eleup,  *hMETPuppi_phi_eleup;
  TH1F *hMETPuppi_pt_eledn,  *hMETPuppi_phi_eledn;
  TH1F *hMETPuppi_pt_tauup,  *hMETPuppi_phi_tauup;
  TH1F *hMETPuppi_pt_taudn,  *hMETPuppi_phi_taudn;
  TH1F *hMETPuppi_pt_uncup,  *hMETPuppi_phi_uncup;
  TH1F *hMETPuppi_pt_uncdn,  *hMETPuppi_phi_uncdn;
  TH1F *hMETPuppi_px_jecup;
  TH1F *hMETPuppi_px_jecdn;
  TH1F *hMETPuppi_px_muonup;
  TH1F *hMETPuppi_px_muondn;
  TH1F *hMETPuppi_px_eleup;
  TH1F *hMETPuppi_px_eledn;
  TH1F *hMETPuppi_px_tauup;
  TH1F *hMETPuppi_px_taudn;
  TH1F *hMETPuppi_px_uncup;
  TH1F *hMETPuppi_px_uncdn;
  // MET Puppi Recorrect
  TH1F *hMETPuppiRecorrect_pt, *hMETPuppiRecorrect_phi;
  TH1F *hMETPuppiRecorrect_pt_jecup,  *hMETPuppiRecorrect_phi_jecup;
  TH1F *hMETPuppiRecorrect_pt_jecdn,  *hMETPuppiRecorrect_phi_jecdn;
  TH1F *hMETPuppiRecorrect_pt_muonup, *hMETPuppiRecorrect_phi_muonup;
  TH1F *hMETPuppiRecorrect_pt_muondn, *hMETPuppiRecorrect_phi_muondn;
  TH1F *hMETPuppiRecorrect_pt_eleup,  *hMETPuppiRecorrect_phi_eleup;
  TH1F *hMETPuppiRecorrect_pt_eledn,  *hMETPuppiRecorrect_phi_eledn;
  TH1F *hMETPuppiRecorrect_pt_tauup,  *hMETPuppiRecorrect_phi_tauup;
  TH1F *hMETPuppiRecorrect_pt_taudn,  *hMETPuppiRecorrect_phi_taudn;
  TH1F *hMETPuppiRecorrect_pt_uncup,  *hMETPuppiRecorrect_phi_uncup;
  TH1F *hMETPuppiRecorrect_pt_uncdn,  *hMETPuppiRecorrect_phi_uncdn;
  TH1F *hMETPuppiRecorrect_px_jecup;
  TH1F *hMETPuppiRecorrect_px_jecdn;
  TH1F *hMETPuppiRecorrect_px_muonup;
  TH1F *hMETPuppiRecorrect_px_muondn;
  TH1F *hMETPuppiRecorrect_px_eleup;
  TH1F *hMETPuppiRecorrect_px_eledn;
  TH1F *hMETPuppiRecorrect_px_tauup;
  TH1F *hMETPuppiRecorrect_px_taudn;
  TH1F *hMETPuppiRecorrect_px_uncup;
  TH1F *hMETPuppiRecorrect_px_uncdn;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
METUncertainties::METUncertainties(const edm::ParameterSet& iConfig)

{
  metSrcTag_ = iConfig.getUntrackedParameter<edm::InputTag>("metSrc");
  metSrcToken_ = consumes<edm::View<pat::MET> >(metSrcTag_);

  metmodifiedSrcTag_ = iConfig.getUntrackedParameter<edm::InputTag>("metmodifiedSrc");
  metmodifiedSrcToken_ = consumes<edm::View<pat::MET> >(metmodifiedSrcTag_);

  metPuppiSrcTag_ = iConfig.getUntrackedParameter<edm::InputTag>("metPuppiSrc");
  metPuppiSrcToken_ = consumes<edm::View<pat::MET> >(metPuppiSrcTag_);

  metPuppiRecorrectSrcTag_ = iConfig.getUntrackedParameter<edm::InputTag>("metPuppiRecorrectSrc");
  metPuppiRecorrectSrcToken_ = consumes<edm::View<pat::MET> >(metPuppiRecorrectSrcTag_);

  nEvent = 0;

  edm::Service<TFileService> fs;
  hMET_pt         = fs -> make<TH1F>("hMET_pt",         "hMET_pt",         100, 0.0, 200.0);
  hMET_phi        = fs -> make<TH1F>("hMET_phi",        "hMET_phi",        100, 0.0, 200.0);
  hMET_pt_jecup   = fs -> make<TH1F>("hMET_pt_jecup",   "hMET_pt_jecup",   100, 0.0, 200.0);
  hMET_phi_jecup  = fs -> make<TH1F>("hMET_phi_jecup",  "hMET_phi_jecup",  100, 0.0, 200.0);
  hMET_pt_jecdn   = fs -> make<TH1F>("hMET_pt_jecdn",   "hMET_pt_jecdn",   100, 0.0, 200.0);
  hMET_phi_jecdn  = fs -> make<TH1F>("hMET_phi_jecdn",  "hMET_phi_jecdn",  100, 0.0, 200.0);
  hMET_pt_muonup  = fs -> make<TH1F>("hMET_pt_muonup",  "hMET_pt_muonup",  100, 0.0, 200.0);
  hMET_phi_muonup = fs -> make<TH1F>("hMET_phi_muonup", "hMET_phi_muonup", 100, 0.0, 200.0);
  hMET_pt_muondn  = fs -> make<TH1F>("hMET_pt_muondn",  "hMET_pt_muondn",  100, 0.0, 200.0);
  hMET_phi_muondn = fs -> make<TH1F>("hMET_phi_muondn", "hMET_phi_muondn", 100, 0.0, 200.0);
  hMET_pt_eleup   = fs -> make<TH1F>("hMET_pt_eleup",   "hMET_pt_eleup",   100, 0.0, 200.0);
  hMET_phi_eleup  = fs -> make<TH1F>("hMET_phi_eleup",  "hMET_phi_eleup",  100, 0.0, 200.0);
  hMET_pt_eledn   = fs -> make<TH1F>("hMET_pt_eledn",   "hMET_pt_eledn",   100, 0.0, 200.0);
  hMET_phi_eledn  = fs -> make<TH1F>("hMET_phi_eledn",  "hMET_phi_eledn",  100, 0.0, 200.0);
  hMET_pt_tauup   = fs -> make<TH1F>("hMET_pt_tauup",   "hMET_pt_tauup",   100, 0.0, 200.0);
  hMET_phi_tauup  = fs -> make<TH1F>("hMET_phi_tauup",  "hMET_phi_tauup",  100, 0.0, 200.0);
  hMET_pt_taudn   = fs -> make<TH1F>("hMET_pt_taudn",   "hMET_pt_taudn",   100, 0.0, 200.0);
  hMET_phi_taudn  = fs -> make<TH1F>("hMET_phi_taudn",  "hMET_phi_taudn",  100, 0.0, 200.0);
  hMET_pt_uncup   = fs -> make<TH1F>("hMET_pt_uncup",   "hMET_pt_uncup",   100, 0.0, 200.0);
  hMET_phi_uncup  = fs -> make<TH1F>("hMET_phi_uncup",  "hMET_phi_uncup",  100, 0.0, 200.0);
  hMET_pt_uncdn   = fs -> make<TH1F>("hMET_pt_uncdn",   "hMET_phi_uncdn",  100, 0.0, 200.0);
  hMET_phi_uncdn  = fs -> make<TH1F>("hMET_phi_uncdn",  "hMET_phi_uncdn",  100, 0.0, 200.0);
  //----------------------------------------------------------------------------------------
  hMETmod_pt         = fs -> make<TH1F>("hMETmod_pt",         "hMETmod_pt",         100, 0.0, 200.0);
  hMETmod_phi        = fs -> make<TH1F>("hMETmod_phi",        "hMETmod_phi",        100, 0.0, 200.0);
  hMETmod_pt_jecup   = fs -> make<TH1F>("hMETmod_pt_jecup",   "hMETmod_pt_jecup",   100, 0.0, 200.0);
  hMETmod_phi_jecup  = fs -> make<TH1F>("hMETmod_phi_jecup",  "hMETmod_phi_jecup",  100, 0.0, 200.0);
  hMETmod_pt_jecdn   = fs -> make<TH1F>("hMETmod_pt_jecdn",   "hMETmod_pt_jecdn",   100, 0.0, 200.0);
  hMETmod_phi_jecdn  = fs -> make<TH1F>("hMETmod_phi_jecdn",  "hMETmod_phi_jecdn",  100, 0.0, 200.0);
  hMETmod_pt_muonup  = fs -> make<TH1F>("hMETmod_pt_muonup",  "hMETmod_pt_muonup",  100, 0.0, 200.0);
  hMETmod_phi_muonup = fs -> make<TH1F>("hMETmod_phi_muonup", "hMETmod_phi_muonup", 100, 0.0, 200.0);
  hMETmod_pt_muondn  = fs -> make<TH1F>("hMETmod_pt_muondn",  "hMETmod_pt_muondn",  100, 0.0, 200.0);
  hMETmod_phi_muondn = fs -> make<TH1F>("hMETmod_phi_muondn", "hMETmod_phi_muondn", 100, 0.0, 200.0);
  hMETmod_pt_eleup   = fs -> make<TH1F>("hMETmod_pt_eleup",   "hMETmod_pt_eleup",   100, 0.0, 200.0);
  hMETmod_phi_eleup  = fs -> make<TH1F>("hMETmod_phi_eleup",  "hMETmod_phi_eleup",  100, 0.0, 200.0);
  hMETmod_pt_eledn   = fs -> make<TH1F>("hMETmod_pt_eledn",   "hMETmod_pt_eledn",   100, 0.0, 200.0);
  hMETmod_phi_eledn  = fs -> make<TH1F>("hMETmod_phi_eledn",  "hMETmod_phi_eledn",  100, 0.0, 200.0);
  hMETmod_pt_tauup   = fs -> make<TH1F>("hMETmod_pt_tauup",   "hMETmod_pt_tauup",   100, 0.0, 200.0);
  hMETmod_phi_tauup  = fs -> make<TH1F>("hMETmod_phi_tauup",  "hMETmod_phi_tauup",  100, 0.0, 200.0);
  hMETmod_pt_taudn   = fs -> make<TH1F>("hMETmod_pt_taudn",   "hMETmod_pt_taudn",   100, 0.0, 200.0);
  hMETmod_phi_taudn  = fs -> make<TH1F>("hMETmod_phi_taudn",  "hMETmod_phi_taudn",  100, 0.0, 200.0);
  hMETmod_pt_uncup   = fs -> make<TH1F>("hMETmod_pt_uncup",   "hMETmod_pt_uncup",   100, 0.0, 200.0);
  hMETmod_phi_uncup  = fs -> make<TH1F>("hMETmod_phi_uncup",  "hMETmod_phi_uncup",  100, 0.0, 200.0);
  hMETmod_pt_uncdn   = fs -> make<TH1F>("hMETmod_pt_uncdn",   "hMETmod_pt_uncdn",   100, 0.0, 200.0);
  hMETmod_phi_uncdn  = fs -> make<TH1F>("hMETmod_phi_uncdn",  "hMETmod_phi_uncdn",  100, 0.0, 200.0);
  //-------------------------------------------------------------------------------------------------
  hMETPuppi_pt         = fs -> make<TH1F>("hMETPuppi_pt",         "hMETPuppi_pt",         100, 0.0, 200.0);
  hMETPuppi_phi        = fs -> make<TH1F>("hMETPuppi_phi",        "hMETPuppi_phi",        100, 0.0, 200.0);
  hMETPuppi_pt_jecup   = fs -> make<TH1F>("hMETPuppi_pt_jecup",   "hMETPuppi_pt_jecup",   100, 0.0, 200.0);
  hMETPuppi_phi_jecup  = fs -> make<TH1F>("hMETPuppi_phi_jecup",  "hMETPuppi_phi_jecup",  100, 0.0, 200.0);
  hMETPuppi_pt_jecdn   = fs -> make<TH1F>("hMETPuppi_pt_jecdn",   "hMETPuppi_pt_jecdn",   100, 0.0, 200.0);
  hMETPuppi_phi_jecdn  = fs -> make<TH1F>("hMETPuppi_phi_jecdn",  "hMETPuppi_phi_jecdn",  100, 0.0, 200.0);
  hMETPuppi_pt_muonup  = fs -> make<TH1F>("hMETPuppi_pt_muonup",  "hMETPuppi_pt_muonup",  100, 0.0, 200.0);
  hMETPuppi_phi_muonup = fs -> make<TH1F>("hMETPuppi_phi_muonup", "hMETPuppi_phi_muonup", 100, 0.0, 200.0);
  hMETPuppi_pt_muondn  = fs -> make<TH1F>("hMETPuppi_pt_muondn",  "hMETPuppi_pt_muondn",  100, 0.0, 200.0);
  hMETPuppi_phi_muondn = fs -> make<TH1F>("hMETPuppi_phi_muondn", "hMETPuppi_phi_muondn", 100, 0.0, 200.0);
  hMETPuppi_pt_eleup   = fs -> make<TH1F>("hMETPuppi_pt_eleup",   "hMETPuppi_pt_eleup",   100, 0.0, 200.0);
  hMETPuppi_phi_eleup  = fs -> make<TH1F>("hMETPuppi_phi_eleup",  "hMETPuppi_phi_eleup",  100, 0.0, 200.0);
  hMETPuppi_pt_eledn   = fs -> make<TH1F>("hMETPuppi_pt_eledn",   "hMETPuppi_pt_eledn",   100, 0.0, 200.0);
  hMETPuppi_phi_eledn  = fs -> make<TH1F>("hMETPuppi_phi_eledn",  "hMETPuppi_phi_eledn",  100, 0.0, 200.0);
  hMETPuppi_pt_tauup   = fs -> make<TH1F>("hMETPuppi_pt_tauup",   "hMETPuppi_pt_tauup",   100, 0.0, 200.0);
  hMETPuppi_phi_tauup  = fs -> make<TH1F>("hMETPuppi_phi_tauup",  "hMETPuppi_phi_tauup",  100, 0.0, 200.0);
  hMETPuppi_pt_taudn   = fs -> make<TH1F>("hMETPuppi_pt_taudn",   "hMETPuppi_pt_taudn",   100, 0.0, 200.0);
  hMETPuppi_phi_taudn  = fs -> make<TH1F>("hMETPuppi_phi_taudn",  "hMETPuppi_phi_taudn",  100, 0.0, 200.0);
  hMETPuppi_pt_uncup   = fs -> make<TH1F>("hMETPuppi_pt_uncup",   "hMETPuppi_pt_uncup",   100, 0.0, 200.0);
  hMETPuppi_phi_uncup  = fs -> make<TH1F>("hMETPuppi_phi_uncup",  "hMETPuppi_phi_uncup",  100, 0.0, 200.0);
  hMETPuppi_pt_uncdn   = fs -> make<TH1F>("hMETPuppi_pt_uncdn",   "hMETPuppi_pt_uncdn",   100, 0.0, 200.0);
  hMETPuppi_phi_uncdn  = fs -> make<TH1F>("hMETPuppi_phi_uncdn",  "hMETPuppi_phi_uncdn",  100, 0.0, 200.0);
  hMETPuppi_px_jecup   = fs -> make<TH1F>("hMETPuppi_px_jecup",   "hMETPuppi_px_jecup",   100, 0.0, 200.0);
  hMETPuppi_px_jecdn   = fs -> make<TH1F>("hMETPuppi_px_jecdn",   "hMETPuppi_px_jecdn",   100, 0.0, 200.0);
  hMETPuppi_px_muonup  = fs -> make<TH1F>("hMETPuppi_px_muonup",  "hMETPuppi_px_muonup",  100, 0.0, 200.0);
  hMETPuppi_px_muondn  = fs -> make<TH1F>("hMETPuppi_px_muondn",  "hMETPuppi_px_muondn",  100, 0.0, 200.0);
  hMETPuppi_px_eleup   = fs -> make<TH1F>("hMETPuppi_px_eleup",   "hMETPuppi_px_eleup",   100, 0.0, 200.0);
  hMETPuppi_px_eledn   = fs -> make<TH1F>("hMETPuppi_px_eledn",   "hMETPuppi_px_eledn",   100, 0.0, 200.0);
  hMETPuppi_px_tauup   = fs -> make<TH1F>("hMETPuppi_px_tauup",   "hMETPuppi_px_tauup",   100, 0.0, 200.0);
  hMETPuppi_px_taudn   = fs -> make<TH1F>("hMETPuppi_px_taudn",   "hMETPuppi_px_taudn",   100, 0.0, 200.0);
  hMETPuppi_px_uncup   = fs -> make<TH1F>("hMETPuppi_px_uncup",   "hMETPuppi_px_uncup",   100, 0.0, 200.0);
  hMETPuppi_px_uncdn   = fs -> make<TH1F>("hMETPuppi_px_uncdn",   "hMETPuppi_px_uncdn",   100, 0.0, 200.0);
  //--------------------------------------------------------------------------------------------------------
  hMETPuppiRecorrect_pt         = fs -> make<TH1F>("hMETPuppiRecorrect_pt",         "hMETPuppiRecorrect_pt",         100, 0.0, 200.0);
  hMETPuppiRecorrect_phi        = fs -> make<TH1F>("hMETPuppiRecorrect_phi",        "hMETPuppiRecorrect_phi",        100, 0.0, 200.0);
  hMETPuppiRecorrect_pt_jecup   = fs -> make<TH1F>("hMETPuppiRecorrect_pt_jecup",   "hMETPuppiRecorrect_pt_jecup",   100, 0.0, 200.0);
  hMETPuppiRecorrect_phi_jecup  = fs -> make<TH1F>("hMETPuppiRecorrect_phi_jecup",  "hMETPuppiRecorrect_phi_jecup",  100, 0.0, 200.0);
  hMETPuppiRecorrect_pt_jecdn   = fs -> make<TH1F>("hMETPuppiRecorrect_pt_jecdn",   "hMETPuppiRecorrect_pt_jecdn",   100, 0.0, 200.0);
  hMETPuppiRecorrect_phi_jecdn  = fs -> make<TH1F>("hMETPuppiRecorrect_phi_jecdn",  "hMETPuppiRecorrect_phi_jecdn",  100, 0.0, 200.0);
  hMETPuppiRecorrect_pt_muonup  = fs -> make<TH1F>("hMETPuppiRecorrect_pt_muonup",  "hMETPuppiRecorrect_pt_muonup",  100, 0.0, 200.0);
  hMETPuppiRecorrect_phi_muonup = fs -> make<TH1F>("hMETPuppiRecorrect_phi_muonup", "hMETPuppiRecorrect_phi_muonup", 100, 0.0, 200.0);
  hMETPuppiRecorrect_pt_muondn  = fs -> make<TH1F>("hMETPuppiRecorrect_pt_muondn",  "hMETPuppiRecorrect_pt_muondn",  100, 0.0, 200.0);
  hMETPuppiRecorrect_phi_muondn = fs -> make<TH1F>("hMETPuppiRecorrect_phi_muondn", "hMETPuppiRecorrect_phi_muondn", 100, 0.0, 200.0);
  hMETPuppiRecorrect_pt_eleup   = fs -> make<TH1F>("hMETPuppiRecorrect_pt_eleup",   "hMETPuppiRecorrect_pt_eleup",   100, 0.0, 200.0);
  hMETPuppiRecorrect_phi_eleup  = fs -> make<TH1F>("hMETPuppiRecorrect_phi_eleup",  "hMETPuppiRecorrect_phi_eleup",  100, 0.0, 200.0);
  hMETPuppiRecorrect_pt_eledn   = fs -> make<TH1F>("hMETPuppiRecorrect_pt_eledn",   "hMETPuppiRecorrect_pt_eledn",   100, 0.0, 200.0);
  hMETPuppiRecorrect_phi_eledn  = fs -> make<TH1F>("hMETPuppiRecorrect_phi_eledn",  "hMETPuppiRecorrect_phi_eledn",  100, 0.0, 200.0);
  hMETPuppiRecorrect_pt_tauup   = fs -> make<TH1F>("hMETPuppiRecorrect_pt_tauup",   "hMETPuppiRecorrect_pt_tauup",   100, 0.0, 200.0);
  hMETPuppiRecorrect_phi_tauup  = fs -> make<TH1F>("hMETPuppiRecorrect_phi_tauup",  "hMETPuppiRecorrect_phi_tauup",  100, 0.0, 200.0);
  hMETPuppiRecorrect_pt_taudn   = fs -> make<TH1F>("hMETPuppiRecorrect_pt_taudn",   "hMETPuppiRecorrect_pt_taudn",   100, 0.0, 200.0);
  hMETPuppiRecorrect_phi_taudn  = fs -> make<TH1F>("hMETPuppiRecorrect_phi_taudn",  "hMETPuppiRecorrect_phi_taudn",  100, 0.0, 200.0);
  hMETPuppiRecorrect_pt_uncup   = fs -> make<TH1F>("hMETPuppiRecorrect_pt_uncup",   "hMETPuppiRecorrect_pt_uncup",   100, 0.0, 200.0);
  hMETPuppiRecorrect_phi_uncup  = fs -> make<TH1F>("hMETPuppiRecorrect_phi_uncup",  "hMETPuppiRecorrect_phi_uncup",  100, 0.0, 200.0);
  hMETPuppiRecorrect_pt_uncdn   = fs -> make<TH1F>("hMETPuppiRecorrect_pt_uncdn",   "hMETPuppiRecorrect_pt_uncdn",   100, 0.0, 200.0);
  hMETPuppiRecorrect_phi_uncdn  = fs -> make<TH1F>("hMETPuppiRecorrect_phi_uncdn",  "hMETPuppiRecorrect_phi_uncdn",  100, 0.0, 200.0);
  hMETPuppiRecorrect_px_jecup   = fs -> make<TH1F>("hMETPuppiRecorrect_px_jecup",   "hMETPuppiRecorrect_px_jecup",   100, 0.0, 200.0);
  hMETPuppiRecorrect_px_jecdn   = fs -> make<TH1F>("hMETPuppiRecorrect_px_jecdn",   "hMETPuppiRecorrect_px_jecdn",   100, 0.0, 200.0);
  hMETPuppiRecorrect_px_muonup  = fs -> make<TH1F>("hMETPuppiRecorrect_px_muonup",  "hMETPuppiRecorrect_px_muonup",  100, 0.0, 200.0);
  hMETPuppiRecorrect_px_muondn  = fs -> make<TH1F>("hMETPuppiRecorrect_px_muondn",  "hMETPuppiRecorrect_px_muondn",  100, 0.0, 200.0);
  hMETPuppiRecorrect_px_eleup   = fs -> make<TH1F>("hMETPuppiRecorrect_px_eleup",   "hMETPuppiRecorrect_px_eleup",   100, 0.0, 200.0);
  hMETPuppiRecorrect_px_eledn   = fs -> make<TH1F>("hMETPuppiRecorrect_px_eledn",   "hMETPuppiRecorrect_px_eledn",   100, 0.0, 200.0);
  hMETPuppiRecorrect_px_tauup   = fs -> make<TH1F>("hMETPuppiRecorrect_px_tauup",   "hMETPuppiRecorrect_px_tauup",   100, 0.0, 200.0);
  hMETPuppiRecorrect_px_taudn   = fs -> make<TH1F>("hMETPuppiRecorrect_px_taudn",   "hMETPuppiRecorrect_px_taudn",   100, 0.0, 200.0);
  hMETPuppiRecorrect_px_uncup   = fs -> make<TH1F>("hMETPuppiRecorrect_px_uncup",   "hMETPuppiRecorrect_px_uncup",   100, 0.0, 200.0);
  hMETPuppiRecorrect_px_uncdn   = fs -> make<TH1F>("hMETPuppiRecorrect_px_uncdn",   "hMETPuppiRecorrect_px_uncdn",   100, 0.0, 200.0);


   //now do what ever initialization is needed
   usesResource("TFileService");

}


METUncertainties::~METUncertainties()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
METUncertainties::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   edm::Handle<edm::View<pat::MET> > meth;
   iEvent.getByToken(metSrcToken_, meth);
   const pat::MET &met = meth -> front();

   edm::Handle<edm::View<pat::MET> > metmodh;
   iEvent.getByToken(metmodifiedSrcToken_, metmodh);
   const pat::MET &metmod = metmodh -> front();

   edm::Handle<edm::View<pat::MET> > metpuppih;
   iEvent.getByToken(metPuppiSrcToken_, metpuppih);
   const pat::MET &metpuppi = metpuppih -> front();

   edm::Handle<edm::View<pat::MET> > metpuppirecorrecth;
   iEvent.getByToken(metPuppiRecorrectSrcToken_, metpuppirecorrecth);
   const pat::MET &metpuppirecorrect = metpuppirecorrecth -> front();

   // Fill histograms
   //----------------
   hMET_pt  -> Fill(met.pt());
   hMET_phi -> Fill(met.phi());
   hMET_pt_jecup   -> Fill(met.shiftedPt(pat::MET::JetEnUp));
   hMET_phi_jecup  -> Fill(met.shiftedPhi(pat::MET::JetEnUp));
   hMET_pt_jecdn   -> Fill(met.shiftedPt(pat::MET::JetEnDown));
   hMET_phi_jecdn  -> Fill(met.shiftedPhi(pat::MET::JetEnDown));
   hMET_pt_muonup  -> Fill(met.shiftedPt(pat::MET::MuonEnUp));
   hMET_phi_muonup -> Fill(met.shiftedPhi(pat::MET::MuonEnUp));
   hMET_pt_muondn  -> Fill(met.shiftedPt(pat::MET::MuonEnDown));
   hMET_phi_muondn -> Fill(met.shiftedPhi(pat::MET::MuonEnDown));
   hMET_pt_eleup   -> Fill(met.shiftedPt(pat::MET::ElectronEnUp)); 
   hMET_phi_eleup  -> Fill(met.shiftedPhi(pat::MET::ElectronEnUp));
   hMET_pt_eledn   -> Fill(met.shiftedPt(pat::MET::ElectronEnDown));
   hMET_phi_eledn  -> Fill(met.shiftedPhi(pat::MET::ElectronEnDown));
   hMET_pt_tauup   -> Fill(met.shiftedPt(pat::MET::TauEnUp));
   hMET_phi_tauup  -> Fill(met.shiftedPhi(pat::MET::TauEnUp));
   hMET_pt_taudn   -> Fill(met.shiftedPt(pat::MET::TauEnDown));
   hMET_phi_taudn  -> Fill(met.shiftedPhi(pat::MET::TauEnDown));
   hMET_pt_uncup   -> Fill(met.shiftedPt(pat::MET::UnclusteredEnUp));
   hMET_pt_uncdn   -> Fill(met.shiftedPt(pat::MET::UnclusteredEnDown));
   hMET_phi_uncup  -> Fill(met.shiftedPhi(pat::MET::UnclusteredEnUp));
   hMET_pt_uncdn   -> Fill(met.shiftedPt(pat::MET::UnclusteredEnDown));
   hMET_phi_uncdn  -> Fill(met.shiftedPhi(pat::MET::UnclusteredEnDown));
   //--------------------------------------------------------------------------------------
   hMETmod_pt  -> Fill(metmod.pt());
   hMETmod_phi -> Fill(metmod.phi());
   hMETmod_pt_jecup   -> Fill(metmod.shiftedPt(pat::MET::JetEnUp));
   hMETmod_phi_jecup  -> Fill(metmod.shiftedPhi(pat::MET::JetEnUp));
   hMETmod_pt_jecdn   -> Fill(metmod.shiftedPt(pat::MET::JetEnDown));
   hMETmod_phi_jecdn  -> Fill(metmod.shiftedPhi(pat::MET::JetEnDown));
   hMETmod_pt_muonup  -> Fill(metmod.shiftedPt(pat::MET::MuonEnUp));
   hMETmod_phi_muonup -> Fill(metmod.shiftedPhi(pat::MET::MuonEnUp));
   hMETmod_pt_muondn  -> Fill(metmod.shiftedPt(pat::MET::MuonEnDown));
   hMETmod_phi_muondn -> Fill(metmod.shiftedPhi(pat::MET::MuonEnDown));
   hMETmod_pt_eleup   -> Fill(metmod.shiftedPt(pat::MET::ElectronEnUp)); 
   hMETmod_phi_eleup  -> Fill(metmod.shiftedPhi(pat::MET::ElectronEnUp));
   hMETmod_pt_eledn   -> Fill(metmod.shiftedPt(pat::MET::ElectronEnDown));
   hMETmod_phi_eledn  -> Fill(metmod.shiftedPhi(pat::MET::ElectronEnDown));
   hMETmod_pt_tauup   -> Fill(metmod.shiftedPt(pat::MET::TauEnUp));
   hMETmod_phi_tauup  -> Fill(metmod.shiftedPhi(pat::MET::TauEnUp));
   hMETmod_pt_taudn   -> Fill(metmod.shiftedPt(pat::MET::TauEnDown));
   hMETmod_phi_taudn  -> Fill(metmod.shiftedPhi(pat::MET::TauEnDown));
   hMETmod_pt_uncup   -> Fill(metmod.shiftedPt(pat::MET::UnclusteredEnUp));
   hMETmod_phi_uncup  -> Fill(metmod.shiftedPhi(pat::MET::UnclusteredEnUp));
   hMETmod_pt_uncdn   -> Fill(metmod.shiftedPt(pat::MET::UnclusteredEnDown));
   hMETmod_phi_uncdn  -> Fill(metmod.shiftedPhi(pat::MET::UnclusteredEnDown));
   //---------------------------------------------------------------------------------------------
   hMETPuppi_pt  -> Fill(metpuppi.pt());
   hMETPuppi_phi -> Fill(metpuppi.phi());
   hMETPuppi_pt_jecup   -> Fill(metpuppi.shiftedPt(pat::MET::JetEnUp));
   hMETPuppi_phi_jecup  -> Fill(metpuppi.shiftedPhi(pat::MET::JetEnUp));
   hMETPuppi_pt_jecdn   -> Fill(metpuppi.shiftedPt(pat::MET::JetEnDown));
   hMETPuppi_phi_jecdn  -> Fill(metpuppi.shiftedPhi(pat::MET::JetEnDown));
   hMETPuppi_pt_muonup  -> Fill(metpuppi.shiftedPt(pat::MET::MuonEnUp));
   hMETPuppi_phi_muonup -> Fill(metpuppi.shiftedPhi(pat::MET::MuonEnUp));
   hMETPuppi_pt_muondn  -> Fill(metpuppi.shiftedPt(pat::MET::MuonEnDown));
   hMETPuppi_phi_muondn -> Fill(metpuppi.shiftedPhi(pat::MET::MuonEnDown));
   hMETPuppi_pt_eleup   -> Fill(metpuppi.shiftedPt(pat::MET::ElectronEnUp)); 
   hMETPuppi_phi_eleup  -> Fill(metpuppi.shiftedPhi(pat::MET::ElectronEnUp));
   hMETPuppi_pt_eledn   -> Fill(metpuppi.shiftedPt(pat::MET::ElectronEnDown));
   hMETPuppi_phi_eledn  -> Fill(metpuppi.shiftedPhi(pat::MET::ElectronEnDown));
   hMETPuppi_pt_tauup   -> Fill(metpuppi.shiftedPt(pat::MET::TauEnUp));
   hMETPuppi_phi_tauup  -> Fill(metpuppi.shiftedPhi(pat::MET::TauEnUp));
   hMETPuppi_pt_taudn   -> Fill(metpuppi.shiftedPt(pat::MET::TauEnDown));
   hMETPuppi_phi_taudn  -> Fill(metpuppi.shiftedPhi(pat::MET::TauEnDown));
   hMETPuppi_pt_uncup   -> Fill(metpuppi.shiftedPt(pat::MET::UnclusteredEnUp));
   hMETPuppi_pt_uncdn   -> Fill(metpuppi.shiftedPt(pat::MET::UnclusteredEnDown));
   hMETPuppi_phi_uncup  -> Fill(metpuppi.shiftedPhi(pat::MET::UnclusteredEnUp));
   hMETPuppi_pt_uncdn   -> Fill(metpuppi.shiftedPt(pat::MET::UnclusteredEnDown));
   hMETPuppi_phi_uncdn  -> Fill(metpuppi.shiftedPhi(pat::MET::UnclusteredEnDown));
   hMETPuppi_px_jecup   -> Fill(metpuppi.shiftedPx(pat::MET::JetEnUp));
   hMETPuppi_px_jecdn   -> Fill(metpuppi.shiftedPx(pat::MET::JetEnDown));
   hMETPuppi_px_muonup  -> Fill(metpuppi.shiftedPx(pat::MET::MuonEnUp));
   hMETPuppi_px_muondn  -> Fill(metpuppi.shiftedPx(pat::MET::MuonEnDown));
   hMETPuppi_px_eleup   -> Fill(metpuppi.shiftedPx(pat::MET::ElectronEnUp));
   hMETPuppi_px_eledn   -> Fill(metpuppi.shiftedPx(pat::MET::ElectronEnDown));
   hMETPuppi_px_tauup   -> Fill(metpuppi.shiftedPx(pat::MET::TauEnUp));
   hMETPuppi_px_taudn   -> Fill(metpuppi.shiftedPx(pat::MET::TauEnDown));
   hMETPuppi_px_uncup   -> Fill(metpuppi.shiftedPx(pat::MET::UnclusteredEnUp));
   hMETPuppi_px_uncdn   -> Fill(metpuppi.shiftedPx(pat::MET::UnclusteredEnDown));
   //--------------------------------------------------------------------------------------------------
   hMETPuppiRecorrect_pt  -> Fill(metpuppirecorrect.pt());
   hMETPuppiRecorrect_phi -> Fill(metpuppirecorrect.phi());
   hMETPuppiRecorrect_pt_jecup   -> Fill(metpuppirecorrect.shiftedPt(pat::MET::JetEnUp));
   hMETPuppiRecorrect_phi_jecup  -> Fill(metpuppirecorrect.shiftedPhi(pat::MET::JetEnUp));
   hMETPuppiRecorrect_pt_jecdn   -> Fill(metpuppirecorrect.shiftedPt(pat::MET::JetEnDown));
   hMETPuppiRecorrect_phi_jecdn  -> Fill(metpuppirecorrect.shiftedPhi(pat::MET::JetEnDown));
   hMETPuppiRecorrect_pt_muonup  -> Fill(metpuppirecorrect.shiftedPt(pat::MET::MuonEnUp));
   hMETPuppiRecorrect_phi_muonup -> Fill(metpuppirecorrect.shiftedPhi(pat::MET::MuonEnUp));
   hMETPuppiRecorrect_pt_muondn  -> Fill(metpuppirecorrect.shiftedPt(pat::MET::MuonEnDown));
   hMETPuppiRecorrect_phi_muondn -> Fill(metpuppirecorrect.shiftedPhi(pat::MET::MuonEnDown));
   hMETPuppiRecorrect_pt_eleup   -> Fill(metpuppirecorrect.shiftedPt(pat::MET::ElectronEnUp)); 
   hMETPuppiRecorrect_phi_eleup  -> Fill(metpuppirecorrect.shiftedPhi(pat::MET::ElectronEnUp));
   hMETPuppiRecorrect_pt_eledn   -> Fill(metpuppirecorrect.shiftedPt(pat::MET::ElectronEnDown));
   hMETPuppiRecorrect_phi_eledn  -> Fill(metpuppirecorrect.shiftedPhi(pat::MET::ElectronEnDown));
   hMETPuppiRecorrect_pt_tauup   -> Fill(metpuppirecorrect.shiftedPt(pat::MET::TauEnUp));
   hMETPuppiRecorrect_phi_tauup  -> Fill(metpuppirecorrect.shiftedPhi(pat::MET::TauEnUp));
   hMETPuppiRecorrect_pt_taudn   -> Fill(metpuppirecorrect.shiftedPt(pat::MET::TauEnDown));
   hMETPuppiRecorrect_phi_taudn  -> Fill(metpuppirecorrect.shiftedPhi(pat::MET::TauEnDown));
   hMETPuppiRecorrect_pt_uncup   -> Fill(metpuppirecorrect.shiftedPt(pat::MET::UnclusteredEnUp));
   hMETPuppiRecorrect_pt_uncdn   -> Fill(metpuppirecorrect.shiftedPt(pat::MET::UnclusteredEnDown));
   hMETPuppiRecorrect_phi_uncup  -> Fill(metpuppirecorrect.shiftedPhi(pat::MET::UnclusteredEnUp));
   hMETPuppiRecorrect_pt_uncdn   -> Fill(metpuppirecorrect.shiftedPt(pat::MET::UnclusteredEnDown));
   hMETPuppiRecorrect_phi_uncdn  -> Fill(metpuppirecorrect.shiftedPhi(pat::MET::UnclusteredEnDown));
   hMETPuppiRecorrect_px_jecup   -> Fill(metpuppirecorrect.shiftedPx(pat::MET::JetEnUp));
   hMETPuppiRecorrect_px_jecdn   -> Fill(metpuppirecorrect.shiftedPx(pat::MET::JetEnDown));
   hMETPuppiRecorrect_px_muonup  -> Fill(metpuppirecorrect.shiftedPx(pat::MET::MuonEnUp));
   hMETPuppiRecorrect_px_muondn  -> Fill(metpuppirecorrect.shiftedPx(pat::MET::MuonEnDown));
   hMETPuppiRecorrect_px_eleup   -> Fill(metpuppirecorrect.shiftedPx(pat::MET::ElectronEnUp));
   hMETPuppiRecorrect_px_eledn   -> Fill(metpuppirecorrect.shiftedPx(pat::MET::ElectronEnDown));
   hMETPuppiRecorrect_px_tauup   -> Fill(metpuppirecorrect.shiftedPx(pat::MET::TauEnUp));
   hMETPuppiRecorrect_px_taudn   -> Fill(metpuppirecorrect.shiftedPx(pat::MET::TauEnDown));
   hMETPuppiRecorrect_px_uncup   -> Fill(metpuppirecorrect.shiftedPx(pat::MET::UnclusteredEnUp));
   hMETPuppiRecorrect_px_uncdn   -> Fill(metpuppirecorrect.shiftedPx(pat::MET::UnclusteredEnDown));


   nEvent ++;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
METUncertainties::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METUncertainties::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
METUncertainties::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(METUncertainties);
