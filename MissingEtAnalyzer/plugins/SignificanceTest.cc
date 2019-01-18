// -*- C++ -*-
//
// Package:    MET_ANALYSIS/MissingEtAnalyzer
// Class:      SignificanceTest
//
/**\class SignificanceTest SignificanceTest.cc MET_ANALYSIS/MissingEtAnalyzer/plugins/SignificanceTest.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Aimilios Ioannou
//         Created:  Wed, 16 Jan 2019 08:50:08 GMT
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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/MET.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class SignificanceTest : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit SignificanceTest(const edm::ParameterSet&);
      ~SignificanceTest();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  edm::InputTag metSrcTag_;
  edm::InputTag metmodifiedSrcTag_;
  edm::EDGetTokenT <edm::View<pat::MET> > metSrcToken_;
  edm::EDGetTokenT <edm::View<pat::MET> > metmodifiedSrcToken_;

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
SignificanceTest::SignificanceTest(const edm::ParameterSet& iConfig)

{

  metSrcTag_   = iConfig.getUntrackedParameter<edm::InputTag>("metSrc");
  metSrcToken_ = consumes<edm::View<pat::MET> >(metSrcTag_);

  metmodifiedSrcTag_   = iConfig.getUntrackedParameter<edm::InputTag>("metmodifiedSrc");
  metmodifiedSrcToken_ = consumes<edm::View<pat::MET> >(metmodifiedSrcTag_);
   //now do what ever initialization is needed


}


SignificanceTest::~SignificanceTest()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SignificanceTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   edm::Handle<edm::View<pat::MET> > meth;
   iEvent.getByToken(metSrcToken_, meth);
   const pat::MET &met = meth -> front();
   
   edm::Handle<edm::View<pat::MET> > metmodifiedh;
   iEvent.getByToken(metmodifiedSrcToken_, metmodifiedh);
   const pat::MET &metmodified = metmodifiedh -> front();

   // MiniAOD MET:
   auto cov = met.getSignificanceMatrix();
   double covXX = cov[0][0];
   double covXY = cov[0][1];
   double covYY = cov[1][1];
   double met_phi = met.phi();
   double met_pt = met.corPt(pat::MET::Type1);
   double met_rawpt = met.corPt(pat::MET::Raw);
   double met_significance = met.significance();
   double met_sumEt = met.sumEt();

   // Modified MET:
   auto covmod = metmodified.getSignificanceMatrix();
   double covXXmod = covmod[0][0];
   double covXYmod = covmod[0][1];
   double covYYmod = covmod[1][1];
   double metmod_phi = metmodified.phi();
   double metmod_pt = metmodified.corPt(pat::MET::Type1);
   double metmod_rawpt = metmodified.corPt(pat::MET::Raw);
   double metmod_significance = metmodified.significance();
   double metmod_sumEt = metmodified.sumEt();


   //--- Print event information ----//
   std::cout << "----------------------------------------------------------" << std::endl;
   
   std::cout << "Run: " << iEvent.id().run() << " Lumi: " << iEvent.id().luminosityBlock() << " Event: " << iEvent.id().event() << std::endl;
   
   std::cout << "---------------- MET Information ------------------------" << std::endl;
   
   std::cout << "MINIAOD:      " << " covXX = " << covXX << " covXY = " << covXY << " covYY = " << covYY << " Significance = " << met_significance  <<" Phi = " << met_phi << " Type1 = " << met_pt << " Raw = " << met_rawpt << " SumEt = " << met_sumEt << std::endl;
   
   std::cout << "ReappliedJEC: " << " covXX = " << covXXmod << " covXY = " << covXYmod << " covYY = " << covYYmod << " Significance = " << metmod_significance  <<" Phi = " << metmod_phi << " Type1 = " << metmod_pt << " Raw = " << metmod_rawpt << " SumEt = " << metmod_sumEt << std::endl;
   

 
   

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
SignificanceTest::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SignificanceTest::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SignificanceTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SignificanceTest);
