import FWCore.ParameterSet.Config as cms

generalAndHiPixelTracks = cms.EDProducer('generalAndHiPixelTracks',
  srcPFcand = cms.InputTag("packedPFCandidates"),
  chi2Map = cms.InputTag("packedPFCandidateTrackChi2"),
  genTrackSrc = cms.InputTag("generalTracks"),
  pixTrackSrc = cms.InputTag("hiPixelTracks"),
  vertexSrc = cms.InputTag("offlineSlimmedPrimaryVerticesRecovery"),
  centralitySrc = cms.InputTag("centralityBin","HFtowers"),
  mvaSrc = cms.InputTag('generalTracks','MVAValues'),
  cutWidth = cms.double(0.2),
  trkRes = cms.double(0.02),
  qualityString = cms.string("highPurity"),
  dxyErrMax = cms.double(3.0),
  dzErrMax = cms.double(3.0),
  ptErrMax = cms.double(0.1),
  nhitsMin = cms.int32(11),
  chi2nMax = cms.double(0.18),
  chi2nMaxPixel = cms.double(9999.0),
  dzErrMaxPixel = cms.double(20.0),
  dxyErrMaxPixel = cms.double(7.0)
)

