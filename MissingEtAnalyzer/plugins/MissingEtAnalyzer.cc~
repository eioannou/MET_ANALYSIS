// -*- C++ -*-
//
// Package:    MET_Analyzer/MissingEtAnalyzer
// Class:      MissingEtAnalyzer
// 
/**\class MissingEtAnalyzer MissingEtAnalyzer.cc MET_Analyzer/MissingEtAnalyzer/plugins/MissingEtAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Aimilios Ioannou
//         Created:  Wed, 18 Jul 2018 05:39:08 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"

// MET stuff
#include "DataFormats/PatCandidates/interface/MET.h"

// Muon stuff
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

// Jet stuff
#include "DataFormats/PatCandidates/interface/Jet.h"

// Vertex stuff
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Trigger stuff
#include "DataFormats/Common/interface/TriggerResults.h"

// TFile Service stuff
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "TTree.h"
#include "vector" 
 
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MissingEtAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MissingEtAnalyzer(const edm::ParameterSet&);
      ~MissingEtAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

  edm::InputTag metSrcTag_;
  edm::InputTag metmodifiedSrcTag_;
  edm::InputTag muonSrcTag_;
  edm::InputTag jetSrcTag_;
  edm::InputTag verticesSrcTag_;
  edm::EDGetTokenT <edm::View<pat::MET> > metSrcToken_;
  edm::EDGetTokenT <edm::View<pat::MET> > metmodifiedSrcToken_;
  edm::EDGetTokenT <edm::View<pat::Muon> > muonSrcToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetSrcToken_;
  edm::EDGetTokenT <edm::TriggerResults> noiseFilterTag_;
  edm::EDGetTokenT <reco::VertexCollection> verticesSrcToken_;

  std::string HBHENoiseIsoFilter_Selector;
  std::string HBHENoiseFilter_Selector;
  std::string EEBadScNoiseFilter_Selector;
  std::string GoodVtxNoiseFilter_Selector;
  std::string GlobalTightHalo2016NoiseFilter_Selector;
  std::string EcalDeadCellTriggerPrimitiveNoiseFilter_Selector;
  std::string BadPFMuonFilter_Selector;
  std::string BadChargedCandidateFilter_Selector;
  std::string ecalBadCalibFilter_Selector;





  int nEvent;

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
MissingEtAnalyzer::MissingEtAnalyzer(const edm::ParameterSet& iConfig):
  noiseFilterTag_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilterTag")))

{

  metSrcTag_    = iConfig.getUntrackedParameter<edm::InputTag>("metSrc");
  metSrcToken_  = consumes<edm::View<pat::MET> >(metSrcTag_);

  metmodifiedSrcTag_ = iConfig.getUntrackedParameter<edm::InputTag>("metmodifiedSrc");
  metmodifiedSrcToken_ = consumes<edm::View<pat::MET> > (metmodifiedSrcTag_);
  
  muonSrcTag_   = iConfig.getUntrackedParameter<edm::InputTag>("muonSrc");
  muonSrcToken_ = consumes<edm::View<pat::Muon> >(muonSrcTag_);
  
  jetSrcTag_   = iConfig.getUntrackedParameter<edm::InputTag>("jetSrc");
  jetSrcToken_ = consumes<edm::View<pat::Jet> >(jetSrcTag_); 

  verticesSrcTag_   = iConfig.getUntrackedParameter<edm::InputTag>("verticesSrc");
  verticesSrcToken_ = consumes<reco::VertexCollection>(verticesSrcTag_);
  
  noiseFilterTag_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilterTag"));

  HBHENoiseFilter_Selector     = iConfig.getParameter<std::string>("HBHENoiseFilter_Selector");
  EEBadScNoiseFilter_Selector = iConfig.getParameter<std::string>("EEBadScNoiseFilter_Selector");
  GoodVtxNoiseFilter_Selector  = iConfig.getParameter<std::string>("GoodVtxNoiseFilter_Selector");
  GlobalTightHalo2016NoiseFilter_Selector = iConfig.getParameter<std::string>("GlobalTightHalo2016NoiseFilter_Selector");
  HBHENoiseIsoFilter_Selector = iConfig.getParameter<std::string>("HBHENoiseIsoFilter_Selector");
  EcalDeadCellTriggerPrimitiveNoiseFilter_Selector = iConfig.getParameter<std::string>("EcalDeadCellTriggerPrimitiveNoiseFilter_Selector");
  BadPFMuonFilter_Selector = iConfig.getParameter<std::string>("BadPFMuonFilter_Selector");
  BadChargedCandidateFilter_Selector = iConfig.getParameter<std::string>("BadChargedCandidateFilter_Selector");
  ecalBadCalibFilter_Selector = iConfig.getParameter<std::string>("ecalBadCalibFilter_Selector");


  //-- Initialize ---                                                                            
  nEvent = 0;


  //now do what ever initialization is needed
  usesResource("TFileService");
  
}


MissingEtAnalyzer::~MissingEtAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MissingEtAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle<edm::View<pat::MET> > meth;
   iEvent.getByToken(metSrcToken_, meth);
   const pat::MET &met = meth -> front();

   edm::Handle<edm::View<pat::MET> > metmodifiedh;
   iEvent.getByToken(metmodifiedSrcToken_, metmodifiedh);
   const pat::MET &metmodified =metmodifiedh -> front();

   edm::Handle<edm::View<pat::Jet> > jeth;
   iEvent.getByToken(jetSrcToken_, jeth); 

   auto cov = met.getSignificanceMatrix();
   double covXX = cov[0][0];
   double covXY = cov[0][1];
   double covYY = cov[1][1];
   double met_phi = met.phi();
   double met_pt  = met.corPt(pat::MET::Type1);

   double met_significance = met.significance();
   double met_sumEt = met.sumEt();

   // Modified MET:
   auto covmod = metmodified.getSignificanceMatrix();
   double covXXmod = covmod[0][0];
   double covXYmod = covmod[0][1];
   double covYYmod = covmod[1][1];
   double metmod_phi = metmodified.phi();
   double metmod_pt  = metmodified.pt();
   double metmod_significance = metmodified.significance();
   double metmod_sumEt = metmodified.sumEt();

   //---Print event information --//       
   std::cout << "=============================================================================" << std::endl;
   std::cout << "RUN: " << iEvent.id().run() << " LUMI: " << iEvent.id().luminosityBlock() << " Event: " << iEvent.id().event() << std::endl;
   std::cout << "****** Information of Event ******************" << std::endl;
   std::cout << "MET_covXX    = " << covXX << " MET_covXY    = " << covXY << " MET_covYY    = " << covYY << " MET_phi    = " << met_phi << " MET_pt    = " << met_pt << " MET_significance    = " << met_significance << " MET_sumEt    = " << met_sumEt << std::endl;
   std::cout << "METMod_covXX = " << covXXmod << " METMod_covXY = " << covXYmod << " METMod_covYY = " << covYYmod << " METMod_phi = " << metmod_phi << " METMod_pt = " << metmod_pt << " METMod_significance = " << metmod_significance << " METMod_sumEt = " << metmod_sumEt << std::endl;



   //for( const pat::Jet &jet : *jeth ){
   //std::cout << " Jet Pt =" << jet.pt() << " Jet Eta =" << jet.eta() << std::endl;
   //}
   //std::cout << "MISSING ET = " << met.pt() << std::endl;
   // std::cout << "Type 1 MET = " << met.corPt(pat::MET::Type1) << std::endl;

   //std::cout << "=========================================================================================" << std::endl;

   //std::cout << "Shifted Up = " << met.shiftedPt(pat::MET::JetEnUp) << " Shifted Down = " << met.shiftedPt(pat::MET::JetEnDown) << std::endl; 

   //std::cout << " Calo MET = " << met.caloMETPt() << std::endl;
   //std::cout << " Calo MET Phi = " << met.caloMETPhi() << " Calo MET Et = " << met.caloMETSumEt() << std::endl;


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
MissingEtAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MissingEtAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MissingEtAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MissingEtAnalyzer);
