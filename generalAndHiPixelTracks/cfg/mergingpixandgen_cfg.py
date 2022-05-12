import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackPoducer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("HeavyIonsAnalysis.Configuration.hfCoincFilter_cff")
process.load('MergingProducer.generalAndHiPixelTracks.MergingPixAndGenProducer_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

       #/store/relval/CMSSW_3_3_1/RelValTTbar/GEN-SIM-RECO/MC_31X_V9-v1/0002/4017BFFA-E9BF-DE11-AC36-001617C3B76E.root
      #'file:/eos/cms/store/hidata/HIRun2018A/HIMinimumBias13/AOD/PromptReco-v1/000/326/623/00000/5EC2FCC1-2B5D-9445-B39C-54B0632A29EB.root'
      #'/store/hidata/HIRun2018A/HIMinimumBias18/AOD/04Apr2019-v1/110000/0271013E-1B89-E24F-8577-52EDD2E6106E.root'
      '/store/hidata/HIRun2018A/HIMinimumBias10/MINIAOD/PbPb18_MiniAODv1-v1/100000/02620e1d-2291-4d67-bd92-81e06f2a7159.root'
    )
)

#from ProdTutorial.TrackAndPointsProducer.trackandpointsproducer_cfi import *
#process.MergingPixAndGen = cms.EDProducer('MergingPixAndGen',
#  genTrackSrc = cms.InputTag("generalTracks"),
#  pixTrackSrc = cms.InputTag("hiConformalPixelTracks"),
#  centralitySrc = cms.InputTag("centralityBin","HFtowers")
#
#)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
    ,outputCommands = cms.untracked.vstring( 'keep *',#'drop *',
      "drop *_generalTracks_*_*",
      "drop *_hiConformalPixelTracks_*_*") 
)

process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

#centrality
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '112X_dataRun2_PromptLike_HI_v3', '')
#process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

"""
print('\n\033[31m~*~ USING CENTRALITY TABLE FOR PbPb 2018 ~*~\033[0m\n')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1031x01_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
"""
 
#process.p = cms.Path(process.MuonTrackPoints*process.TrackTrackPoints)
process.p = cms.Path(
	#process.centralityBin *
	process.generalAndHiPixelTracks)

process.e = cms.EndPath(process.out)

