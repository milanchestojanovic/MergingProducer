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
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
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

      bool passesGeneralTrackCuts(const reco::Track & track, const reco::VertexCollection&, double chi2n);
      bool passesPixelTrackCuts(const reco::Track & track, const reco::VertexCollection&);


      edm::EDGetTokenT<reco::TrackCollection>  genTrackSrc_;
      //edm::EDGetTokenT<reco::TrackCollection>  pixTrackSrc_;
      edm::EDGetTokenT<edm::View<pat::PackedCandidate>  > pixTrackSrc_; 
      edm::EDGetTokenT<edm::View<pat::PackedCandidate>  > srcPFcand_;
      edm::InputTag chi2MapLabel_;
      edm::EDGetTokenT<edm::ValueMap<float>> chi2Map_;

      edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
      edm::EDGetTokenT<int> centralitySrc_;
      edm::EDGetTokenT<std::vector<float>> mvaSrc_;

      double cutWidth_;
      double trkRes_;
      std::string qualityString_;
      double dxyErrMax_;
      double dzErrMax_;
      double ptErrMax_;
      int    nhitsMin_;
      double chi2nMax_;
      double chi2nMaxPixel_;
      double dzErrMaxPixel_;
      double dxyErrMaxPixel_;
      float trkMva;

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
 //pixTrackSrc_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("pixTrackSrc"))),
 pixTrackSrc_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("pixTrackSrc"))),
 srcPFcand_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("srcPFcand"))),
 chi2MapLabel_(iConfig.getParameter<edm::InputTag>("chi2Map")),
 chi2Map_(consumes<edm::ValueMap<float>>(chi2MapLabel_)),
 vertexSrc_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
 centralitySrc_(consumes<int>(iConfig.getParameter<edm::InputTag>("centralitySrc"))),
 mvaSrc_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("mvaSrc"))),
 cutWidth_(iConfig.getParameter<double>("cutWidth")),
 trkRes_(iConfig.getParameter<double>("trkRes")),
 qualityString_(iConfig.getParameter<std::string>("qualityString")),
 dxyErrMax_(iConfig.getParameter<double>("dxyErrMax")),
 dzErrMax_(iConfig.getParameter<double>("dzErrMax")),
 ptErrMax_(iConfig.getParameter<double>("ptErrMax")),
 nhitsMin_(iConfig.getParameter<int>("nhitsMin")),
 chi2nMax_(iConfig.getParameter<double>("chi2nMax")),
 chi2nMaxPixel_(iConfig.getParameter<double>("chi2nMaxPixel")),
 dzErrMaxPixel_(iConfig.getParameter<double>("dzErrMaxPixel")),
 dxyErrMaxPixel_(iConfig.getParameter<double>("dxyErrMaxPixel"))
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

  Handle<reco::VertexCollection> vertex;
  iEvent.getByToken(vertexSrc_, vertex);

  VertexCollection vsorted = *vertex;
  const VertexCollection *recoV = vertex.product();


  Handle<reco::TrackCollection> genTcol;
  iEvent.getByToken(genTrackSrc_, genTcol);

  //Handle<reco::TrackCollection> pixTcol;
  edm::Handle<edm::View<pat::PackedCandidate>> pixTcol;
  iEvent.getByToken(pixTrackSrc_, pixTcol);

  edm::Handle<edm::View<pat::PackedCandidate>> pfCandidates;
  iEvent.getByToken(srcPFcand_,pfCandidates);

  edm::Handle<edm::ValueMap<float>> chi2Map;
  iEvent.getByToken(chi2Map_,chi2Map);

  edm::Handle<std::vector<float>> mvaoutput;
  iEvent.getByToken(mvaSrc_, mvaoutput);


  int occ=10, cbin=0;
  edm::Handle<int> centralityBin;
  iEvent.getByToken(centralitySrc_, centralityBin);
  cbin = *centralityBin;
  occ = cbin;

  double ptCut = 0.6;
  if (occ < 20) ptCut = 1.0;

  //int i_tr=-1;
  for(unsigned int i = 0, n = pfCandidates->size(); i < n; ++i){
     const pat::PackedCandidate &pf = (*pfCandidates)[i];
      if(!(pf.hasTrackDetails()))continue;
       const reco::Track  &trk = pf.pseudoTrack();
       const reco::Track * tr = (&trk);
       //const reco::Track tr = *trk;
    //trkMva = (*mvaoutput)[i_tr]; 
    //i_tr++;
    if ( tr->pt() < ptCut ) continue;
    double chi2n=(*chi2Map)[pfCandidates->ptrAt(i)];
    if ( ! passesGeneralTrackCuts(*tr, vsorted, chi2n) ) continue;
    if (tr->pt() < (ptCut + cutWidth_)  ) {
	genEta.push_back(tr->eta());
        genPhi.push_back(tr->phi());
    }
    mrgTcol->push_back(*tr);
  } 

  for(unsigned int i = 0, n = pixTcol->size(); i < n; ++i){
     const pat::PackedCandidate &pf = (*pixTcol)[i];
      if(!(pf.hasTrackDetails()))continue;
       const reco::Track  &trk = pf.pseudoTrack();
       const reco::Track * tr = (&trk);

     int nrec=0;
     if (tr->pt() >= ptCut ) continue;
     if ( ! passesPixelTrackCuts(*tr, vsorted) ) continue;
     if ( tr->pt() > (ptCut - cutWidth_) ) {
           int Ngen = genEta.size();
           for(int i=0; i< Ngen; i++){

             if (fabs(tr->eta()- genEta[i]) < trkRes_ && fabs(tr->phi()- genPhi[i]) < trkRes_ ) {
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

bool
generalAndHiPixelTracks::passesGeneralTrackCuts(const reco::Track & track, const reco::VertexCollection & vertex, double chi2n)
{

   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
   int primaryvtx = 0;
   vtxPoint=vertex[primaryvtx].position();
   vzErr=vertex[primaryvtx].zError();
   vxErr=vertex[primaryvtx].xError();
   vyErr=vertex[primaryvtx].yError();

   double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
   dxy = track.dxy(vtxPoint);
   dz = track.dz(vtxPoint);
   dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
   dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);

   //double chi2n = track.normalizedChi2();
   double nlayers = track.hitPattern().trackerLayersWithMeasurement();
   chi2n = chi2n/nlayers;
   int nhits = track.numberOfValidHits();
   int algo  = track.algo();


   if(track.quality(reco::TrackBase::qualityByName(qualityString_)) != 1)
       return false;

   if(fabs(dxy/dxysigma) > dxyErrMax_) return false;
   if(fabs(dz/dzsigma) > dzErrMax_) return false;

   if(fabs(track.ptError()) / track.pt() > ptErrMax_) return false;

   if(nhits < nhitsMin_ ) return false;
   if(chi2n > chi2nMax_ ) return false;

   if ( algo == 6 && trkMva < 0.98 ) return false;

   return true;
}

bool
generalAndHiPixelTracks::passesPixelTrackCuts(const reco::Track & track, const reco::VertexCollection & vertex)
{

   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
   int primaryvtx = 0;
   vtxPoint=vertex[primaryvtx].position();
   vzErr=vertex[primaryvtx].zError();
   vxErr=vertex[primaryvtx].xError();
   vyErr=vertex[primaryvtx].yError();

   double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
   dxy = track.dxy(vtxPoint);
   dz = track.dz(vtxPoint);
   dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
   dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);

   double chi2n = track.normalizedChi2();
   double nlayers = track.hitPattern().trackerLayersWithMeasurement();
   chi2n = chi2n/nlayers;

   if (chi2n > chi2nMaxPixel_ ) return false;
   if (fabs(dz/dzsigma) > dzErrMaxPixel_ ) return false;
   if (fabs(dxy/dxysigma) > dxyErrMaxPixel_ ) return false;

   return true;
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
