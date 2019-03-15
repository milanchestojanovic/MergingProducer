// -*- C++ -*-
//
// Package:    MergingProducer/generalAndHiPixelTracks
// Class:      generalAndHiPixelTracks
// 
/**\class generalAndHiPixelTracks generalAndHiPixelTracks.cc MergingProducer/generalAndHiPixelTracks/plugins/generalAndHiPixelTracks.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Milan Stojanovic
//         Created:  Tue, 26 Feb 2019 13:22:38 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Utilities/interface/EDGetToken.h"

#include <vector>
#include <iostream>
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "TVector3.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include <DataFormats/VertexReco/interface/Vertex.h>
//#include <DataFormats/VertexReco/interface/VertexFwd.h>
//#include <DataFormats/TrackReco/interface/Track.h>
//#include <DataFormats/TrackReco/interface/TrackFwd.h>

//
// class declaration
//

class generalAndHiPixelTracks : public edm::stream::EDProducer<> {
   public:
      explicit generalAndHiPixelTracks(const edm::ParameterSet&);
      ~generalAndHiPixelTracks();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;



      edm::EDGetTokenT<reco::TrackCollection>  genTrackSrc_;
      edm::EDGetTokenT<reco::TrackCollection>  pixTrackSrc_;

      edm::EDGetTokenT<int> centralitySrc_;

//      typedef math::XYZPointD Point;
//      typedef std::vector<Point> PointCollection;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

     std::vector<double> genEta;
     std::vector<double> genPhi;

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
generalAndHiPixelTracks::generalAndHiPixelTracks(const edm::ParameterSet& iConfig):
  genTrackSrc_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("genTrackSrc"))),
  pixTrackSrc_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("pixTrackSrc"))),
  centralitySrc_(consumes<int>(iConfig.getParameter<edm::InputTag>("centralitySrc")))
{
   produces<std::vector<reco::Track> >(); 
}


generalAndHiPixelTracks::~generalAndHiPixelTracks()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
generalAndHiPixelTracks::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   using namespace reco;
   using namespace std;

  unique_ptr<vector<reco::Track> > mrgTcol(new vector<reco::Track>());

  Handle<reco::TrackCollection> genTcol;
  iEvent.getByToken(genTrackSrc_, genTcol);

  Handle<reco::TrackCollection> pixTcol;
  iEvent.getByToken(pixTrackSrc_, pixTcol);

  int occ=10, cbin=0;
  edm::Handle<int> centralityBin;
  iEvent.getByToken(centralitySrc_, centralityBin);
  cbin = *centralityBin;
  occ = cbin;

  double ptCut = 0.6;
  if (occ < 20) ptCut = 1.0;


  for(TrackCollection::const_iterator tr = genTcol->begin(); tr != genTcol->end(); tr++) {
    if ( tr->pt() < ptCut ) continue;
    if (tr->pt() < (ptCut + 0.2)  ) {
	genEta.push_back(tr->eta());
        genPhi.push_back(tr->phi());
    }
    mrgTcol->push_back(*tr);
  } 

  for(TrackCollection::const_iterator tr = pixTcol->begin(); tr != pixTcol->end(); tr++) {//COUNTING MULTIPLICITY
     int nrec=0;
     if (tr->pt() >= ptCut ) continue;
     if ( tr->pt() > (ptCut - 0.2) ) {
           int Ngen = genEta.size();
           for(int i=0; i< Ngen; i++){

             if (fabs(tr->eta()- genEta[i]) < 0.02 && fabs(tr->phi()- genPhi[i]) < 0.02 ) {
               nrec++;
               
             }

           }
     }
     if (nrec>0) continue;
     mrgTcol->push_back(*tr);

  }
  genEta.resize(0);
  genPhi.resize(0);

  iEvent.put(std::move(mrgTcol));// (xxxxgeneralAndHiPixelTracks));

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
generalAndHiPixelTracks::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
generalAndHiPixelTracks::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
generalAndHiPixelTracks::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
generalAndHiPixelTracks::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
generalAndHiPixelTracks::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
generalAndHiPixelTracks::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
generalAndHiPixelTracks::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(generalAndHiPixelTracks);
