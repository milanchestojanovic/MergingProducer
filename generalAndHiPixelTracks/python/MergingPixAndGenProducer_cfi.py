import FWCore.ParameterSet.Config as cms

generalAndHiPixelTracks = cms.EDProducer('generalAndHiPixelTracks',
  genTrackSrc = cms.InputTag("generalTracks"),
  pixTrackSrc = cms.InputTag("hiConformalPixelTracks"),
  centralitySrc = cms.InputTag("centralityBin","HFtowers")

)

