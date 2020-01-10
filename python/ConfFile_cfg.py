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
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

        #VBF:
        #'file:/pnfs/desy.de/cms/tier2/store/user/lbenato/VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000_Summer16_AODSIM_Tranche2/VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_Tranche2_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_Tranche2_AODSIM/181214_110750/0000/aodsim_1.root',
        #'file:/pnfs/desy.de/cms/tier2/store/user/lbenato/VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000_Summer16_AODSIM_Tranche2/VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_Tranche2_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_Tranche2_AODSIM/181214_110750/0000/aodsim_2.root',
        #'file:/pnfs/desy.de/cms/tier2/store/user/lbenato/VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000_Summer16_AODSIM_Tranche2/VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_Tranche2_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_Tranche2_AODSIM/181214_110750/0000/aodsim_3.root',
        #'file:/pnfs/desy.de/cms/tier2/store/user/lbenato/VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000_Summer16_AODSIM_Tranche2/VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_Tranche2_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_Tranche2_AODSIM/181214_110750/0000/aodsim_4.root',

        #ggH:
        #'file:/pnfs/desy.de/cms/tier2/store/user/lbenato/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1000_Summer16_AODSIM/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_AODSIM/181128_153659/0000/aodsim_1.root',
        #'file:/pnfs/desy.de/cms/tier2/store/user/lbenato/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1000_Summer16_AODSIM/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_AODSIM/181128_153659/0000/aodsim_2.root',
        #'file:/pnfs/desy.de/cms/tier2/store/user/lbenato/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1000_Summer16_AODSIM/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_AODSIM/181128_153659/0000/aodsim_3.root',
        #'file:/pnfs/desy.de/cms/tier2/store/user/lbenato/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1000_Summer16_AODSIM/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_AODSIM/181128_153659/0000/aodsim_4.root',

        #ZH
        'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISummer16DR80Premix/ZH_HToSSTobbbb_ZToLL_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/D8333AE5-5ECE-E611-861C-001E674FBA1D.root',

        #'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISummer16DR80Premix/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/FAA214B2-3CB8-E611-AD5D-00304867FDCF.root',
        ##'/store/mc/RunIISummer16DR80Premix/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/40D4FD22-CDBA-E611-AE2A-02163E016539.root',
        #'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/38A97D70-9DBE-E611-9811-FA163E90EFC1.root',
        #'/store/mc/RunIISummer16DR80Premix/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/FE441877-55B3-E611-AB31-A0369F7FCDF4.root',
        #'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISummer16DR80Premix/WW_TuneCUETP8M1_13TeV-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/4EC89234-9ED7-E611-91FE-0CC47A4C8E0E.root',
        #'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISummer16DR80Premix/ZJetsToNuNu_HT-200To400_13TeV-madgraph/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/F642326F-88C8-E611-81C9-0CC47AA992B4.root',
        #'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/1C5C7DB3-66BC-E611-ADA2-0242AC130004.root',
        #DYJets aMCatnlo miniaod:
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToNuNu_PtZ-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/FACE8698-A7CD-E611-A684-02163E019C06.root'
        #'file:/pnfs/desy.de/cms/tier2/store/data/Run2016G/MET/AOD/07Aug17-v1/110000/3C4239F2-E9A0-E711-82F7-02163E014117.root',
    )
)

noLHEinfo = True if ('WW_TuneCUETP8M1_13TeV-pythia8' or 'WZ_TuneCUETP8M1_13TeV-pythia8' or 'ZZ_TuneCUETP8M1_13TeV-pythia8') in process.source.fileNames[0] else False

process.TFileService = cms.Service( "TFileService",
    fileName = cms.string('output.root' if len(options.outputFile)==0 else options.outputFile),
    closeFileFast = cms.untracked.bool(True),
)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '') 

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

process.demo = cms.EDAnalyzer('DemoMini' if 'MINIAODSIM' in process.source.fileNames[0] else 'Demo',
    pythiaLOSample = cms.bool(True if noLHEinfo else False),
)

process.seq = cms.Sequence(
    process.demo
)

process.p = cms.Path(process.demo)
