import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.parseArguments()

process = cms.Process("Demo")

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    allowUnscheduled = cms.untracked.bool(True),
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/pnfs/desy.de/cms/tier2/store/user/lbenato/VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000_Summer16_AODSIM_Tranche2/VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_Tranche2_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_Tranche2_AODSIM/181214_110750/0000/aodsim_1.root'
    )
)


process.TFileService = cms.Service( "TFileService",
    fileName = cms.string('output.root' if len(options.outputFile)==0 else options.outputFile),
    closeFileFast = cms.untracked.bool(True),
)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

#-----------------------#
#    VERTEX FILTER      #
#-----------------------#

import RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi
process.primaryVertexFilter = cms.EDFilter('GoodVertexFilter',
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)


#-----------------------#
#       ANALYZER        #
#-----------------------#

process.demo = cms.EDAnalyzer('Demo'
)

process.seq = cms.Sequence(
    process.demo
)

process.p = cms.Path(process.demo)
