import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')


## Var parsing, for CRAB
#options.register(
#    "runLocal", True,
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.bool,
#    "Decide if CRAB should overwrite variables"
#)
options.register(
    "runLocal", True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Decide if CRAB should overwrite variables"
)

options.register(
    "PisData", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "isData parser flag"
)
options.register(
    "PisReHLT", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "isReHLT parser flag"
)
options.register(
    "PisReReco", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "isReReco parser flag"
)
options.register(
    "PisReMiniAod", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "isReMiniAod parser flag"
)
options.register(
    "PisPromptReco", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "isPromptReco parser flag"
)
options.register(
    "PnoLHEinfo", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "noLHEinfo parser flag"
)
options.register(
    "PisbbH", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "isbbH parser flag"
)
options.register(
    "PisSignal", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "isSignal parser flag"
)
options.register(
    "PGT", "",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "GT parser flag"
)
options.register(
    "PJECstring", "",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "JECstring parser flag"
)
options.register(
    "PjsonName", "",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    " jsonName parser flag"
)
options.register(
    "PtriggerTag", "",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "triggerTag parser flag"
)
options.register(
    "PfilterString", "",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "filterString parser flag"
)
options.register(
    "Pcalo", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "calo parser flag"
)

options.parseArguments()

## Important: decide if keeping local options register or CRAB options register
#RunLocalCMS = cms.bool( options.runLocal )
RunLocal = cms.bool( options.runLocal )

if RunLocal:
   print "Taking configurations from cfg file; not designed for CRAB submission!"

process = cms.Process("ntuple")

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    allowUnscheduled = cms.untracked.bool(True),
)

#Enable multithreading!
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Messagge logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

## Input files
if len(options.inputFiles) == 0:

    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
            'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISummer16DR80Premix/ZH_HToSSTobbbb_ZToLL_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/D8333AE5-5ECE-E611-861C-001E674FBA1D.root',
            #'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISummer16DR80Premix/ZH_HToSSTobbbb_ZToLL_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/02A3FAE1-67CE-E611-A13F-001E67457107.root',
            #'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISummer16DR80Premix/ZH_HToSSTobbbb_ZToLL_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/8606DF62-61CE-E611-8FA4-ECF4BBE1DD48.root',
            #'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISummer16DR80Premix/ZH_HToSSTobbbb_ZToLL_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0ADCA7A5-65CE-E611-85E8-0CC47A546E5E.root',
            #'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISummer16DR80Premix/ZH_HToSSTobbbb_ZToLL_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/D6CC71A9-62CE-E611-A892-001E67456EBE.root',
            #'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISummer16DR80Premix/ZH_HToSSTobbbb_ZToLL_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/6ADE619E-64CE-E611-A504-D067E5F914D3.root',
            #'file:/pnfs/desy.de/cms/tier2/store/user/lbenato/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1000_Summer16_AODSIM/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_AODSIM/181128_153659/0000/aodsim_1.root'
            #'file:/pnfs/desy.de/cms/tier2/store/user/lbenato/VBFH_HToSSTobbbb_MH-125_MS-15_ctauS-5000_Summer16_AODSIM_Tranche2/VBFH_HToSSTobbbb_MH-125_MS-15_ctauS-5000_TuneCUETP8M1_13TeV-powheg-pythia8_Tranche2_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_Tranche2_AODSIM/181214_110243/0000/aodsim_1.root'
            #'/store/mc/RunIISummer16DR80Premix/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/FE57DDB4-DDBA-E611-A344-0025905A6064.root',
            #'file:/pnfs/desy.de/cms/tier2/store/data/Run2016G/MET/AOD/07Aug17-v1/110000/3C4239F2-E9A0-E711-82F7-02163E014117.root' 
        )
    )

if RunLocal:
    isData = ('/store/data/' in process.source.fileNames[0])
else:
    isData = options.PisData

process.TFileService = cms.Service( "TFileService",
    fileName = cms.string('output.root' if len(options.outputFile)==0 else options.outputFile),
    closeFileFast = cms.untracked.bool(True),
)


#-----------------------#
#     DATA FLAGS        #
#-----------------------#

if RunLocal:
    isData            = ('/store/data/' in process.source.fileNames[0])
    isReHLT           = ('_reHLT_' in process.source.fileNames[0])
    isReReco          = ('23Sep2016' in process.source.fileNames[0])
    isReMiniAod       = ('03Feb2017' in process.source.fileNames[0])
    isPromptReco      = ('PromptReco' in process.source.fileNames[0])
    noLHEinfo         = True if ('WW_TuneCUETP8M1_13TeV-pythia8' or 'WZ_TuneCUETP8M1_13TeV-pythia8' or 'ZZ_TuneCUETP8M1_13TeV-pythia8') in process.source.fileNames[0] else False #check for PythiaLO samples
    isbbH             = True if ('bbHToBB_M-125_4FS_yb2_13TeV_amcatnlo' in process.source.fileNames[0]) else False #bbH has a different label in LHEEventProduct
    isSignal          = True if ('HToSSTobbbb_MH-125' in process.source.fileNames[0]) else False
    isCalo            = True #HERE for calo analyses!!!

else:
    isData            = options.PisData
    isReHLT           = options.PisReHLT
    isReReco          = options.PisReReco
    isReMiniAod       = options.PisReMiniAod
    isPromptReco      = options.PisPromptReco
    noLHEinfo         = options.PnoLHEinfo
    isbbH             = options.PisbbH
    isSignal          = options.PisSignal
    isCalo            = options.Pcalo

theRunBCD = ['Run2016B','Run2016C','Run2016D']
theRunEF  = ['Run2016E','Run2016F']
theRunG   = ['Run2016G']
theRunH   = ['Run2016H']

print 'isData',isData
print 'isReHLT',isReHLT
print 'isReReco',isReReco
print 'isReMiniAod',isReMiniAod
print 'isPromptReco',isPromptReco
print 'isSignal', isSignal

if isCalo:
    print "\n"
    print "***************************************"
    print "***************************************"
    print "***************************************"
    print "\n"
    print "Performing analysis for CALO LIFETIMES!"
    print "\n"
    print "***************************************"
    print "***************************************"
    print "***************************************"
    print "\n"

#-----------------------#
#     GLOBAL TAG        #
#-----------------------#

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
GT = ''

if RunLocal:
    if isData:
        if isReMiniAod and any(s in process.source.fileNames[0] for s in theRunH): GT = '80X_dataRun2_Prompt_v16'
        else: GT = '80X_dataRun2_2016SeptRepro_v7'#'auto:run2_data'#'80X_dataRun2_2016LegacyRepro_v4'#
    elif not(isData):
        GT = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'#Moriond17 GT
else:
    GT = options.PGT

process.GlobalTag = GlobalTag(process.GlobalTag, GT)
print 'GlobalTag loaded: ', GT

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
#         JEC           #
#-----------------------#

# Jet corrector https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrOnTheFly
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')


JECstring = ''
if RunLocal:
    if isData:# and (isReReco or isReMiniAod):
      if any(s in process.source.fileNames[0] for s in theRunBCD):
        JECstring = "Summer16_23Sep2016BCDV3_DATA" #if isReMiniAod else "Summer16_23Sep2016BCDV3_DATA"
      if any(s in process.source.fileNames[0] for s in theRunEF):
        JECstring = "Summer16_23Sep2016EFV3_DATA" #if isReMiniAod else "Summer16_23Sep2016EFV3_DATA"
      if any(s in process.source.fileNames[0] for s in theRunG):
        JECstring = "Summer16_23Sep2016GV3_DATA" #if isReMiniAod else "Summer16_23Sep2016GV3_DATA"
      if any(s in process.source.fileNames[0] for s in theRunH):
        JECstring = "Summer16_23Sep2016HV3_DATA" #if isReMiniAod else "Summer16_23Sep2016HV3_DATA"
    elif isData and isPromptReco:
        JECstring = "Spring16_25nsV6_DATA"
    elif not isData:
        JECstring = "Summer16_23Sep2016V3_MC"

else:
    JECstring = options.PJECstring
print "JEC ->",JECstring

#-----------------------#
#       COUNTER         #
#-----------------------#
process.counter = cms.EDAnalyzer('CounterAnalyzer',
    lheProduct = cms.InputTag('externalLHEProducer' if not isbbH else 'source'),
    pythiaLOSample = cms.bool(True if noLHEinfo else False),
)


#-----------------------#
#     PAT OBJECTS       #
#-----------------------#

#Transient track builder needed for vertices
process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

#Load Pat sequences
#process.load("PhysicsTools.PatAlgos.patSequences_cff")
#from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *

#Puppi
process.load('CommonTools.PileupAlgos.Puppi_cff')
process.load('PhysicsTools.PatAlgos.slimming.puppiForMET_cff')
process.pfNoLepPUPPI = cms.EDFilter("PdgIdCandViewSelector",
                                    src = cms.InputTag("particleFlow"),
                                    pdgId = cms.vint32( 1,2,22,111,130,310,2112,211,-211,321,-321,999211,2212,-2212 )
                                    )

process.puppiNoLep = process.puppi.clone()
process.puppiNoLep.candName = cms.InputTag('pfNoLepPUPPI')

#packedPFCandidates
process.load('PhysicsTools.PatAlgos.slimming.packedPFCandidates_cff')

#lostTracks
process.load('PhysicsTools.PatAlgos.slimming.lostTracks_cfi')

#Vertex association and slimmed vertices
process.load('PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi')
process.primaryVertexAssociation.jets = cms.InputTag("ak4PFJets")
process.load('PhysicsTools.PatAlgos.slimming.offlineSlimmedPrimaryVertices_cfi')

#genParticles
process.load('PhysicsTools.PatAlgos.slimming.genParticles_cff')
#trigger
process.load('PhysicsTools.PatAlgos.slimming.selectedPatTrigger_cfi')

#slimmedMuons
process.load('PhysicsTools.PatAlgos.producersLayer1.muonProducer_cff')
process.patMuons.addGenMatch = False if isData else True
process.patMuons.embedGenMatch = False if isData else True
process.load('PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi')
process.selectedPatMuons.cut = cms.string("pt > 5 || isPFMuon || (pt > 3 && (isGlobalMuon || isStandAloneMuon || numberOfMatches > 0 || muonID(\'RPCMuLoose\')))")
process.load('PhysicsTools.PatAlgos.slimming.slimmedMuons_cfi')

#reducedEgamma
process.load('RecoEgamma.EgammaPhotonProducers.reducedEgamma_cfi')
process.load('RecoEgamma.EgammaTools.egammaObjectModificationsInMiniAOD_cff')

#patElectrons
process.load('PhysicsTools.PatAlgos.producersLayer1.electronProducer_cff')
process.patElectrons.addGenMatch = False if isData else True
process.patElectrons.embedGenMatch = False if isData else True
#patPhotons
process.load('PhysicsTools.PatAlgos.producersLayer1.photonProducer_cff')
process.patPhotons.addGenMatch = False if isData else True
process.patPhotons.embedGenMatch = False if isData else True

#slimmedTaus
process.load('PhysicsTools.PatAlgos.producersLayer1.tauProducer_cff')
process.patTaus.addGenMatch = False if isData else True
process.patTaus.addGenJetMatch = False if isData else True
process.patTaus.embedGenMatch = False if isData else True
process.patTaus.embedGenJetMatch = False if isData else True
process.load('PhysicsTools.PatAlgos.selectionLayer1.tauSelector_cfi')
process.selectedPatTaus.cut = cms.string("pt > 18. && tauID(\'decayModeFindingNewDMs\')> 0.5")
process.load('PhysicsTools.PatAlgos.slimming.slimmedTaus_cfi')


'''
#slimmedJets
process.load('PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff')
##we might need to add stuff, cross-check
process.load('PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi')
process.selectedPatJets.cut = cms.string("pt > 5.")#default is 10!
#here: AK8 missing, work on that
#process.load('PhysicsTools.PatAlgos.slimming.slimmedJets_cfi')
'''

#-----------------------#
#  E-MU-GAMMA MODULES   #
#-----------------------#

#electron/photon ID maps
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.AOD)
ele_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                  'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                  'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff']

for ele_idmod in ele_id_modules:
    setupAllVIDIdsInModule(process,ele_idmod,setupVIDElectronSelection)

process.heepIDVarValueMaps.dataFormat = cms.int32(1)#0 = auto detection, 1 = AOD, 2 = miniAOD
process.load('PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi')
process.load('PhysicsTools.PatAlgos.slimming.slimmedElectrons_cfi')

#photons upstream modules
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDPhotonIdProducer(process, DataFormat.AOD)

ph_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff',
                 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']

for ph_idmod in ph_id_modules:
    setupAllVIDIdsInModule(process,ph_idmod,setupVIDPhotonSelection)

process.load('PhysicsTools.PatAlgos.selectionLayer1.photonSelector_cfi')
process.load('PhysicsTools.PatAlgos.slimming.slimmedPhotons_cfi')


#electron/photon regression modules
'''
Do we really need that? It looks like slimmedElectrons are needed as input...

from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
  EGMSmearerElectrons = cms.PSet(
    initialSeed = cms.untrackesd.uint32(8675389),
    engineName = cms.untracked.string('TRandom3')
  ),
  EGMSmearerPhotons = cms.PSet(
    initialSeed = cms.untracked.uint32(8675389),
    engineName = cms.untracked.string('TRandom3')
  )
)

process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')
process.EGMRegression       = cms.Sequence(process.regressionApplication)

process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')
from EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi import *
process.EGMSmearerElectrons = calibratedPatElectrons.clone(
  isMC = cms.bool(False if isData else True)
)

process.load('EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi')
from EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi import *
process.EGMSmearerPhotons = calibratedPatPhotons.clone(
  isMC = cms.bool(False if isData else True)
)
'''


#muons upstream modules
process.cleanedMuons = cms.EDProducer('PATMuonCleanerBySegments',
                                      src = cms.InputTag('slimmedMuons'),#('calibratedMuons'),#
                                      preselection = cms.string('track.isNonnull'),
                                      passthrough = cms.string('isGlobalMuon && numberOfMatches >= 2'),
                                      fractionOfSharedSegments = cms.double(0.499)
                                      )




#-----------------------#
#        FILTERS        #
#-----------------------#

# JSON filter
if isData:
    import FWCore.PythonUtilities.LumiList as LumiList
    jsonName = "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON"#"Cert_294927-305364_13TeV_PromptReco_Collisions17_JSON"#"Cert_294927-301567_13TeV_PromptReco_Collisions17_JSON" #golden json
    process.source.lumisToProcess = LumiList.LumiList(filename = 'data/JSON/'+jsonName+'.txt').getVLuminosityBlockRange()
    print "JSON file loaded: ", jsonName

if RunLocal:
    # Trigger filter
    triggerTag = 'HLT2' if isReHLT else 'HLT'

    # MET filters string
    filterString = "RECO"
    if isData:
        filterString = "RECO"
    else:
        filterString = "SIM"#"PAT"
else:
    triggerTag = options.PtriggerTag
    filterString = options.PfilterString


## MET filters, not available on AOD? TODO
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag('slimmedMuons')
process.BadPFMuonFilter.PFCandidates = cms.InputTag('packedPFCandidates')

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag('slimmedMuons')
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag('packedPFCandidates')

#-----------------------#
#  b TAGGING tagInfos   #
#-----------------------#

bTagInfos = [
    'pfImpactParameterTagInfos'
   ,'pfSecondaryVertexTagInfos'
   ,'pfInclusiveSecondaryVertexFinderTagInfos'
]

bTagDiscriminators = [
   'pfCombinedInclusiveSecondaryVertexV2BJetTags',
   'pfCombinedSecondaryVertexV2BJetTags',
   'pfBoostedDoubleSecondaryVertexAK8BJetTags'
   ]


# # taken from here: https://github.com/cms-sw/cmssw/blob/02d4198c0b6615287fd88e9a8ff650aea994412e/RecoBTag/ImpactParameter/python/impactParameterTagInfos_cfi.py
process.load("RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi")
process.pfImpactParameterTagInfos.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pfImpactParameterTagInfos.maximumChiSquared = cms.double(99999.9)
process.pfImpactParameterTagInfos.maximumLongitudinalImpactParameter = cms.double(99999.9)
process.pfImpactParameterTagInfos.maximumTransverseImpactParameter = cms.double(99999.9)
process.pfImpactParameterTagInfos.minimumNumberOfHits = cms.int32(0)
process.pfImpactParameterTagInfos.minimumNumberOfPixelHits = cms.int32(1)#at least 1, otherwise tracking issues!13Dec
process.pfImpactParameterTagInfos.minimumTransverseMomentum = cms.double(1.0)
process.pfImpactParameterTagInfos.computeGhostTrack = cms.bool(True)

process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

process.load("RecoBTag.SecondaryVertex.pfSecondaryVertexTagInfos_cfi")
process.load('RecoBTag/SecondaryVertex/pfInclusiveSecondaryVertexFinderTagInfos_cfi')
process.load("RecoBTag.SecondaryVertex.pfSecondaryVertexNegativeTagInfos_cfi")
process.load('RecoBTag/SecondaryVertex/pfInclusiveSecondaryVertexFinderNegativeTagInfos_cfi')
process.load('RecoBTag/CTagging/pfInclusiveSecondaryVertexFinderCvsLTagInfos_cfi')
process.load('RecoBTag/CTagging/pfInclusiveSecondaryVertexFinderNegativeCvsLTagInfos_cfi')

for n in ['pfSecondaryVertexTagInfos', 'pfInclusiveSecondaryVertexFinderTagInfos', 'pfSecondaryVertexNegativeTagInfos', 'pfInclusiveSecondaryVertexFinderNegativeTagInfos', 'pfInclusiveSecondaryVertexFinderCvsLTagInfos', 'pfInclusiveSecondaryVertexFinderNegativeCvsLTagInfos']:
    setattr( getattr(process,n), 'useExternalSV', cms.bool(False) )
    setattr( getattr(process,n), 'extSVCollection', cms.InputTag("") )
    setattr( getattr(process,n), 'useSVClustering', cms.bool(True) )
    setattr( getattr(process,n), 'jetAlgorithm', cms.string("AntiKt") )
    setattr( getattr(process,n), 'rParam', cms.double(0.4) )
    setattr( getattr(process,n), 'trackSelection', cms.PSet(
            max_pT = cms.double(99999.9),
            max_pT_dRcut = cms.double(99999.9),
            max_pT_trackPTcut = cms.double(99999.9),
            min_pT = cms.double(-99999.9),
            min_pT_dRcut = cms.double(-99999.9),
            pixelHitsMin = cms.uint32(0),
            totalHitsMin = cms.uint32(0)
            ))
    setattr( getattr(process,n), 'vertexCuts', cms.PSet(
            distVal2dMax = cms.double(99999.9),
            distSig2dMin = cms.double(-99999.9),
            distSig2dMax = cms.double(99999.9),
            distVal2dMin = cms.double(-99999.9),
            minimumTrackWeight = cms.double(-99999.9),
            massMax = cms.double(99999.9)
))

process.load('RecoBTag/SoftLepton/softPFMuonTagInfos_cfi')
process.load('RecoBTag/SoftLepton/softPFElectronTagInfos_cfi')

#-----------------------#
#       AK4 JETS        #
#-----------------------#
pt_AK8 = 170 #dummy
#Jet labels
if isCalo:
   chosen_JEC = "AK4PFchs"
   chosen_jet_source = 'ak4PFJetsCHS'
   chosen_label = ''#'Reclustered'
   chosen_pfcand = 'pfCHS'
   chosen_jets = "patJets"+ chosen_label
   pt_AK4 = 5
else:
   #chosen_jets = "slimmedJets"#no longer existing in AOD
   pt_AK4 = 15

## packedPFCandidates with CHS are used by both AK4 and AK8
#process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
#equivalent to:
#process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV() || charge()==0"))

#print "\n"
#print "CLUSTERING ONLY NEUTRAL STUFF!!!"
#print "\n"
#process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("charge==0"))
process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV && (abs(dz)>0.05 || charge==0)"))

## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))

## Define GenJets
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsNoNu = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')

if isCalo:

   print "% % % % % % % % % % % % % % % % % % % % % % % % % %"
   print "Performing AK4PFchs jet reclustering, pT = "+str(pt_AK4)+" GeV"
   print "% % % % % % % % % % % % % % % % % % % % % % % % % %"

   ##Recluster reco jets
   from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
   ##non CHS
   ##process.ak4PFJets = ak4PFJets.clone(src = 'packedPFCandidates', doAreaFastjet = True, jetPtMin = pt_AK4)
   process.ak4PFJetsCHSCustom = ak4PFJets.clone(src = "pfCHS", doAreaFastjet = True, jetPtMin = pt_AK4)
   chosen_jet_source = 'ak4PFJetsCHSCustom'

   from PhysicsTools.PatAlgos.tools.jetTools import *
   addJetCollection(
      process,
      labelName = chosen_label,#'Reclustered',
      jetSource = cms.InputTag(chosen_jet_source),#reco jets
      pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
      pfCandidates = cms.InputTag(chosen_pfcand),#pfchs substracted
      svSource = cms.InputTag('slimmedSecondaryVertices'),
      btagDiscriminators = list(bTagDiscriminators),#btagging
      btagInfos = bTagInfos,
      jetCorrections = (chosen_JEC, ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], 'None'),#correct JEC
      genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
      genParticles = cms.InputTag('prunedGenParticles'),
      algo = 'AK',
      rParam = 0.4
      )

process.patJets.addGenJetMatch = False if isData else True
process.patJets.addGenPartonMatch = False if isData else True
process.patJets.embedGenJetMatch = False if isData else True
process.patJets.embedGenPartonMatch = False if isData else True
process.patJets.addPartonJetMatch = False if isData else True
process.patJets.addJetFlavourInfo = False if isData else True
process.patJets.getJetMCFlavour = False if isData else True


#slimmedMETs
process.load('PhysicsTools.PatAlgos.producersLayer1.metProducer_cff')
process.patMETs.computeMETSignsificance = cms.bool(True)
process.patMETs.srcJets = cms.InputTag("patJets")
process.patMETs.addGenMET = cms.bool(False if isData else True)
process.load('PhysicsTools.PatAlgos.selectionLayer1.metSelector_cfi')

#Calculating met corrections, uncertainties and slimmed mets
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import *
#runMetCorAndUncForMiniAODProduction(process,
#                                    jetSelection="pt>543 && abs(eta)<9.9",#this is important!
#                                    )

runMETCorrectionsAndUncertainties = RunMETCorrectionsAndUncertainties()
##
metType="PF"
#jetCollUnskimmed="patJets"
jetCollUnskimmed="patJets"
photonColl="selectedPatPhotons"
electronColl="selectedPatElectrons"
muonColl="selectedPatMuons"
tauColl="selectedPatTaus"
pfCandColl = "particleFlow"
jetCleaning="LepClean"
jetSelection="pt>543 && abs(eta)<9.9"
jecUnFile=""
jetFlavor="AK4PFchs"
recoMetFromPFCs=False
postfix=""

runMETCorrectionsAndUncertainties(process,
                                  metType=metType,
                                  correctionLevel=["T0","T1","T2","Smear","Txy"],
                                  computeUncertainties=False,
                                  produceIntermediateCorrections=True,
                                  addToPatDefaultSequence=False,
                                  jetCollectionUnskimmed=jetCollUnskimmed,
                                  photonCollection=photonColl,
                                  electronCollection=electronColl,
                                  muonCollection=muonColl,
                                  tauCollection=tauColl,
                                  pfCandCollection =pfCandColl,
                                  autoJetCleaning=jetCleaning,
                                  jecUncertaintyFile=jecUnFile,
                                  jetSelection=jetSelection,
                                  jetFlavor=jetFlavor,
                                  recoMetFromPFCs=recoMetFromPFCs,
                                  postfix=postfix,
                                  runOnData = True if isData else False
                                  )

#MET T1 uncertainties
runMETCorrectionsAndUncertainties(process,
                                  metType=metType,
                                  correctionLevel=["T1"],
                                  computeUncertainties=True,
                                  produceIntermediateCorrections=False,
                                  addToPatDefaultSequence=False,
                                  jetCollectionUnskimmed=jetCollUnskimmed,
                                  photonCollection=photonColl,
                                  electronCollection=electronColl,
                                  muonCollection=muonColl,
                                  tauCollection=tauColl,
                                  pfCandCollection =pfCandColl,
                                  autoJetCleaning=jetCleaning,
                                  jecUncertaintyFile=jecUnFile,
                                  jetSelection=jetSelection,
                                  jetFlavor=jetFlavor,
                                  recoMetFromPFCs=recoMetFromPFCs,
                                  postfix=postfix,
                                  runOnData = True if isData else False
                                  )

#MET T1 Smeared JER uncertainties
runMETCorrectionsAndUncertainties(process,
                                  metType=metType,
                                  correctionLevel=["T1","Smear"],
                                  computeUncertainties=True,
                                  produceIntermediateCorrections=False,
                                  addToPatDefaultSequence=False,
                                  jetCollectionUnskimmed=jetCollUnskimmed,
                                  photonCollection=photonColl,
                                  electronCollection=electronColl,
                                  muonCollection=muonColl,
                                  tauCollection=tauColl,
                                  pfCandCollection =pfCandColl,
                                  autoJetCleaning=jetCleaning,
                                  jecUncertaintyFile=jecUnFile,
                                  jetSelection=jetSelection,
                                  jetFlavor=jetFlavor,
                                  recoMetFromPFCs=recoMetFromPFCs,
                                  postfix=postfix,
                                  runOnData = True if isData else False
                                  )

#import PhysicsTools.PatAlgos.tools.helpers as configtools
#postfix=''
#configtools.removeIfInSequence(process, "selectedPatJetsForMetT1T2Corr", "patPFMetT1T2CorrSequence", postfix )
process.patSmearedJets.skipGenMatching = cms.bool(True if isData else False)

from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs
process.patCaloMet = patMETs.clone(
    metSource = cms.InputTag('caloMetM'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue'),
    addGenMET = cms.bool(False if isData else True)
)
process.load('PhysicsTools.PatAlgos.slimming.slimmedMETs_cfi')


#-----------------------#
#       B-Tag           #
#-----------------------#

from PhysicsTools.PatAlgos.tools.jetTools import *

#Seth
jetSource = chosen_jets
jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
pfCandidates = 'packedPFCandidates'
pvSource = 'offlinePrimaryVertices'#'offlineSlimmedPrimaryVertices'
svSource = 'inclusiveCandidateSecondaryVertices'#'slimmedSecondaryVertices'
muSource = 'slimmedMuons'
elSource = 'selectedPatElectrons'#'slimmedElectrons'

bTagInfos = [
    'pfImpactParameterTagInfos'
   ,'pfSecondaryVertexTagInfos'
   ,'pfInclusiveSecondaryVertexFinderTagInfos'
   ,'pfSecondaryVertexNegativeTagInfos'
   ,'pfInclusiveSecondaryVertexFinderNegativeTagInfos'
   ,'softPFMuonsTagInfos'
   ,'softPFElectronsTagInfos'
   ,'pfInclusiveSecondaryVertexFinderCvsLTagInfos'
   ,'pfInclusiveSecondaryVertexFinderNegativeCvsLTagInfos'
#   ,'pfDeepFlavourTagInfos' # not available for 2016
]


bTagDiscriminators = set([
    'pfJetBProbabilityBJetTags'
   ,'pfJetProbabilityBJetTags'
   ,'pfPositiveOnlyJetBProbabilityBJetTags'
   ,'pfPositiveOnlyJetProbabilityBJetTags'
   ,'pfNegativeOnlyJetBProbabilityBJetTags'
   ,'pfNegativeOnlyJetProbabilityBJetTags'
   ,'pfTrackCountingHighPurBJetTags'
   ,'pfTrackCountingHighEffBJetTags'
   ,'pfNegativeTrackCountingHighPurBJetTags'
   ,'pfNegativeTrackCountingHighEffBJetTags'
   ,'pfSimpleSecondaryVertexHighEffBJetTags'
   ,'pfSimpleSecondaryVertexHighPurBJetTags'
   ,'pfNegativeSimpleSecondaryVertexHighEffBJetTags'
   ,'pfNegativeSimpleSecondaryVertexHighPurBJetTags'
   ,'pfCombinedSecondaryVertexV2BJetTags'
   ,'pfPositiveCombinedSecondaryVertexV2BJetTags'
   ,'pfNegativeCombinedSecondaryVertexV2BJetTags'
   ,'pfCombinedInclusiveSecondaryVertexV2BJetTags'
   ,'pfPositiveCombinedInclusiveSecondaryVertexV2BJetTags'
   ,'pfNegativeCombinedInclusiveSecondaryVertexV2BJetTags'
   ,'softPFMuonBJetTags'
   ,'positiveSoftPFMuonBJetTags'
   ,'negativeSoftPFMuonBJetTags'
   ,'softPFElectronBJetTags'
   ,'positiveSoftPFElectronBJetTags'
   ,'negativeSoftPFElectronBJetTags'
   ,'pfCombinedMVAV2BJetTags'
   ,'pfNegativeCombinedMVAV2BJetTags'
   ,'pfPositiveCombinedMVAV2BJetTags'
   ,'pfCombinedCvsBJetTags'
   ,'pfNegativeCombinedCvsBJetTags'
   ,'pfPositiveCombinedCvsBJetTags'
   ,'pfCombinedCvsLJetTags'
   ,'pfNegativeCombinedCvsLJetTags'
   ,'pfPositiveCombinedCvsLJetTags'
])


useExplicitJTA = False #? try true also!
postfix = "" #"Update"


#-----------------------#
#       Vertices        #
#-----------------------#

# taken from: https://github.com/cms-sw/cmssw/blob/02d4198c0b6615287fd88e9a8ff650aea994412e/PhysicsTools/PatAlgos/test/btag-from-packedPat.py
#Info: these two producers seem not to be used
postfix = 'Final'#TODO


updateJetCollection(
    process,
    jetSource = cms.InputTag(jetSource),
    jetCorrections = jetCorrectionsAK4,
    pfCandidates = cms.InputTag(pfCandidates),
    pvSource = cms.InputTag(pvSource),
    svSource = cms.InputTag(svSource),
    muSource = cms.InputTag(muSource),
    elSource = cms.InputTag(elSource),
    btagInfos = bTagInfos,
    btagDiscriminators = list(bTagDiscriminators),
    explicitJTA = useExplicitJTA,
    postfix = postfix#,
    )

for m in ['updatedPatJets'+postfix, 'updatedPatJetsTransientCorrected'+postfix]:
    setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )

jets_after_btag_tools = 'updatedPatJetsTransientCorrected'+postfix

print("??????????????????????")
#-----------------------#
#       ANALYZER        #
#-----------------------#

process.ntuple = cms.EDAnalyzer('ZHNtuplizer',
    genSet = cms.PSet(
        genProduct = cms.InputTag('generator'),
        lheProduct = cms.InputTag('externalLHEProducer'),
        genParticles = cms.InputTag('genParticles'),
        pdgId = cms.vint32(5,9000006,23,24,25),#(1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 21, 23, 24, 25, 36, 39, 1000022, 9100000, 9000001, 9000002, 9100012, 9100022, 9900032, 1023),
        status = cms.vint32(22,23),
        samplesDYJetsToLL = cms.vstring(),
        samplesZJetsToNuNu = cms.vstring(),
        samplesWJetsToLNu = cms.vstring(),
        samplesDir = cms.string('data/Stitch/'),
        sample = cms.string("" ), #( sample )
        ewkFile = cms.string('data/scalefactors_v4.root'),
        applyEWK = cms.bool(False),#(True if sample.startswith('DYJets') or sample.startswith('WJets') else False),
        applyTopPtReweigth = cms.bool(False),#(True if sample.startswith('TT_') else False),
        pythiaLOSample = cms.bool(True if noLHEinfo else False),#(True if isDibosonInclusive else False),
    ),

    pileupSet = cms.PSet(
        pileup = cms.InputTag('addPileupInfo'),
        vertices = cms.InputTag('offlinePrimaryVertices'),
        dataFileName     = cms.string('data/PU_69200_ReReco.root'),#updated
        dataFileNameUp   = cms.string('data/PU_72380_ReReco.root'),#updated
        dataFileNameDown = cms.string('data/PU_66020_ReReco.root'),#updated
        mcFileName = cms.string('data/PU_MC_Moriond17.root'),#updated
        dataName = cms.string('pileup'),
        mcName = cms.string('2016_25ns_Moriond17MC_PoissonOOTPU'),#updated
    ),
    triggerSet = cms.PSet(
        trigger = cms.InputTag('TriggerResults', '', triggerTag),
        paths = cms.vstring(
*[
##lepton
'HLT_IsoMu24_v', 'HLT_Ele27_WPTight_Gsf_v',
##others:
'HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140_v','HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63_v','HLT_PFHT300_PFMET110_v','HLT_RsqMR270_Rsq0p09_MR200_v','HLT_RsqMR270_Rsq0p09_MR200_4jet_v','HLT_HT650_v',
### b-like
#'HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240_v', 'HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500_v', 'HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v', 'HLT_QuadJet45_TripleBTagCSV_p087_v', 'HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6_v', 'HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172_v',
### displaced tracks
'HLT_VBF_DisplacedJet40_DisplacedTrack_v', 'HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5_v', 'HLT_HT350_DisplacedDijet40_DisplacedTrack_v', 'HLT_HT350_DisplacedDijet80_DisplacedTrack_v', 'HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack_v', 'HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack_v', 'HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack_v', 'HLT_HT650_DisplacedDijet80_Inclusive_v', 'HLT_HT750_DisplacedDijet80_Inclusive_v',
### calo lifetimes
'HLT_VBF_DisplacedJet40_VTightID_Hadronic_v', 'HLT_VBF_DisplacedJet40_VVTightID_Hadronic_v',
###
#'HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067_v', 'HLT_MET200_v', 'HLT_MET250_v', 'HLT_MET75_IsoTrk50_v', 'HLT_MET90_IsoTrk50_v', 'HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_v', 'HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v', 'HLT_PFMET110_PFMHT110_IDTight_v', 'HLT_PFMET120_PFMHT120_IDTight_v', 'HLT_PFMET170_HBHECleaned_v', 'HLT_PFMET300_v', 'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v', 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v',
###production for MET
'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v',
### All studied triggers:
#
###Control paths for VBF, prescaled
#'HLT_L1_TripleJet_VBF_v', 'HLT_QuadPFJet_VBF_v','HLT_DiPFJetAve40_v','HLT_DiPFJetAve60_v','HLT_DiPFJetAve80_v','HLT_PFJet40_v','HLT_PFJet60_v','HLT_PFJet80_v',
###TEST
##'HLT_AK8PFJet450_v',###########TEST
#'HLT_VBF_DisplacedJet40_VTightID_Hadronic_v', 'HLT_VBF_DisplacedJet40_VVTightID_Hadronic_v'#,'HLT_AK4PFJet30_v7'
]
        ),
        metfilters = cms.InputTag('TriggerResults', '', filterString),
        metpaths = cms.vstring('Flag_HBHENoiseFilter', 'Flag_HBHENoiseIsoFilter', 'Flag_EcalDeadCellTriggerPrimitiveFilter', 'Flag_goodVertices', 'Flag_eeBadScFilter', 'Flag_globalTightHalo2016Filter','Flag_badMuons','Flag_duplicateMuons','Flag_noBadMuons') if isReMiniAod else cms.vstring('Flag_HBHENoiseFilter', 'Flag_HBHENoiseIsoFilter', 'Flag_EcalDeadCellTriggerPrimitiveFilter', 'Flag_goodVertices', 'Flag_eeBadScFilter', 'Flag_globalTightHalo2016Filter'),
        #prescales = cms.InputTag('patTrigger','','PAT'),
        #l1Minprescales = cms.InputTag('patTrigger','l1min','PAT'),
        #l1Maxprescales = cms.InputTag('patTrigger','l1max','PAT'),
        #objects = cms.InputTag('selectedPatTrigger','','PAT'),
        badPFMuonFilter = cms.InputTag("BadPFMuonFilter"),
        badChCandFilter = cms.InputTag("BadChargedCandidateFilter"),
        l1Gt = cms.InputTag("gtStage2Digis"),
        #l1filters = cms.vstring('hltL1sTripleJet846848VBFIorTripleJet887256VBFIorTripleJet927664VBFIorHTT300','hltL1sDoubleJetC112','hltL1sQuadJetC50IorQuadJetC60IorHTT280IorHTT300IorHTT320IorTripleJet846848VBFIorTripleJet887256VBFIorTripleJet927664VBF','hltL1sTripleJetVBFIorHTTIorDoubleJetCIorSingleJet','hltL1sSingleMu22','hltL1sV0SingleMu22IorSingleMu25','hltL1sZeroBias','hltL1sSingleJet60','hltL1sSingleJet35','hltTripleJet50','hltDoubleJet65','hltSingleJet80','hltVBFFilterDisplacedJets'),
    ),
    chsJetSet = cms.PSet(
        jets = cms.InputTag(jets_after_btag_tools),#('ak4PFJetsCHS'),#('updatedPatJetsTransientCorrected'+postfix),
        jetid = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight
        jet1pt = cms.double(5),
        jet2pt = cms.double(5),
        jeteta = cms.double(2.4),
        addQGdiscriminator = cms.bool(False),
        recalibrateJets = cms.bool(True),#(True),
        recalibrateMass = cms.bool(False),
        recalibratePuppiMass = cms.bool(False),
        softdropPuppiMassString = cms.string("ak8PFJetsPuppiValueMap:ak8PFJetsPuppiSoftDropMass" if pt_AK8<170 else "ak8PFJetsPuppiSoftDropMass"),
        smearJets = cms.bool(False),#just for comparison!
        vertices = cms.InputTag('offlinePrimaryVertices'),
        rho = cms.InputTag('fixedGridRhoFastjetAll'),
        jecUncertaintyDATA = cms.string('data/%s/%s_Uncertainty_AK4PFchs.txt' % (JECstring, JECstring)),#updating
        jecUncertaintyMC = cms.string('data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt'),#updating
        jecCorrectorDATA = cms.vstring(#updating
            'data/%s/%s_L1FastJet_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L2Relative_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK4PFchs.txt' % (JECstring, JECstring),
        ),
        jecCorrectorMC = cms.vstring(#updating!!!
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt',
        ),
        massCorrectorDATA = cms.vstring(#updating!!!
            'data/%s/%s_L2Relative_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK4PFchs.txt' % (JECstring, JECstring),
        ),
        massCorrectorMC = cms.vstring(#updating!!!
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt',
        ),
        massCorrectorPuppi = cms.string('data/puppiCorrSummer16.root'),#updating
        reshapeBTag = cms.bool(True),
        btag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
        btagDB = cms.string('data/CSVv2_Moriond17_B_H.csv'),
        jet1btag = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight
        jet2btag = cms.int32(0),
        met = cms.InputTag('slimmedMETs'),#('patMETs'),# if isReMiniAod else cms.InputTag('slimmedMETs', '', ''),# 'LLP'
        recomet = cms.InputTag('pfMet'),
        metRecoil = cms.bool(False),
        metRecoilMC = cms.string('data/recoilfit_gjetsMC_Zu1_pf_v5.root'),
        metRecoilData = cms.string('data/recoilfit_gjetsData_Zu1_pf_v5.root'),
        metTriggerFileName = cms.string('data/MET_trigger_eff_data_SingleMuRunBH.root'),
        jerNameRes = cms.string('data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt'),#v10 is the latest
        jerNameSf = cms.string('data/JER/Spring16_25nsV10_MC_SF_AK4PFchs.txt'),#v10 is the latest
    ),

    caloJetSet = cms.PSet(
        jets = cms.InputTag('ak4CaloJets'),
        jet1pt = cms.double(10.),
        jet2pt = cms.double(10.),
        jeteta = cms.double(2.4),
        recalibrateJets = cms.bool(True),
        recalibrateMass = cms.bool(False),
        smearJets = cms.bool(False),
        vertices = cms.InputTag('offlinePrimaryVertices'),
        rho = cms.InputTag('fixedGridRhoFastjetAll'),
        jecUncertaintyDATA = cms.string('data/%s/%s_Uncertainty_AK4Calo.txt' % (JECstring, JECstring)),
        jecUncertaintyMC = cms.string('data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_Uncertainty_AK4Calo.txt'),
        jecCorrectorDATA = cms.vstring(
            'data/%s/%s_L1FastJet_AK4Calo.txt' % (JECstring, JECstring),
            'data/%s/%s_L2Relative_AK4Calo.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK4Calo.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK4Calo.txt' % (JECstring, JECstring),
        ),
        jecCorrectorMC = cms.vstring(
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4Calo.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4Calo.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4Calo.txt',
        ),
        massCorrectorDATA = cms.vstring(
            'data/%s/%s_L2Relative_AK4Calo.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK4Calo.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK4Calo.txt' % (JECstring, JECstring),
        ),
        massCorrectorMC = cms.vstring(                                                         #
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4Calo.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4Calo.txt',
        ),
        jerNameRes = cms.string('data/JER/Spring16_25nsV10_MC_PtResolution_AK8PF.txt'),#NOT PROVIDED FOR CALO JETS
        jerNameSf = cms.string('data/JER/Spring16_25nsV10_MC_SF_AK8PF.txt'),#NOT PROVIDED FOR CALO JETS
    ),

    vbfJetSet = cms.PSet(
        jets = cms.InputTag(jets_after_btag_tools),#('ak4PFJetsCHS'),#('updatedPatJetsTransientCorrected'+postfix),
        jetid = cms.int32(3), # 0: no selection, 1: loose, 2: medium, 3: tight
        ##jet1pt = cms.double(30.),#https://indico.desy.de/indico/event/20983/contribution/0/material/slides/0.pdf
        ##jet2pt = cms.double(30.),#https://indico.desy.de/indico/event/20983/contribution/0/material/slides/0.pdf
        #new cut, motivated by calo-lifetimes trigger path
        #still to be optimized!
        jet1pt = cms.double(20.),
        jet2pt = cms.double(20.),
        jeteta = cms.double(5.2),
        addQGdiscriminator = cms.bool(False),
        recalibrateJets = cms.bool(True),
        recalibrateMass = cms.bool(False),
        recalibratePuppiMass = cms.bool(False),
        softdropPuppiMassString = cms.string("ak8PFJetsPuppiValueMap:ak8PFJetsPuppiSoftDropMass" if pt_AK8<170 else "ak8PFJetsPuppiSoftDropMass"),
        smearJets = cms.bool(False),#just for comparison!
        vertices = cms.InputTag('offlinePrimaryVertices'),
        rho = cms.InputTag('fixedGridRhoFastjetAll'),
        jecUncertaintyDATA = cms.string('data/%s/%s_Uncertainty_AK4PFchs.txt' % (JECstring, JECstring)),#updating
        jecUncertaintyMC = cms.string('data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt'),#updating
        jecCorrectorDATA = cms.vstring(#updating
            'data/%s/%s_L1FastJet_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L2Relative_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK4PFchs.txt' % (JECstring, JECstring),
        ),
        jecCorrectorMC = cms.vstring(#updating!!!
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt',
        ),
        massCorrectorDATA = cms.vstring(#updating!!!
            'data/%s/%s_L2Relative_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK4PFchs.txt' % (JECstring, JECstring),
        ),
        massCorrectorMC = cms.vstring(#updating!!!
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt',
        ),
        massCorrectorPuppi = cms.string('data/puppiCorrSummer16.root'),#updating
        reshapeBTag = cms.bool(True),
        btag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
        btagDB = cms.string('data/CSVv2_Moriond17_B_H.csv'),
        jet1btag = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight
        jet2btag = cms.int32(0),
        met = cms.InputTag('patMETs'),#('patMETs'),#('pfMet'),# if isReMiniAod else cms.InputTag('slimmedMETs', '', ''),# 'LLP'
        recomet = cms.InputTag('pfMet'),
        metRecoil = cms.bool(False),
        metRecoilMC = cms.string('data/recoilfit_gjetsMC_Zu1_pf_v5.root'),
        metRecoilData = cms.string('data/recoilfit_gjetsData_Zu1_pf_v5.root'),
        metTriggerFileName = cms.string('data/MET_trigger_eff_data_SingleMuRunBH.root'),
        jerNameRes = cms.string('data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt'),#v10 is the latest
        jerNameSf = cms.string('data/JER/Spring16_25nsV10_MC_SF_AK4PFchs.txt'),#v10 is the latest
    ),
    electronSet = cms.PSet(
        electrons = cms.InputTag('gedGsfElectrons'),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
        eleVetoIdMap = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto'),
        eleLooseIdMap = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose'),
        eleMediumIdMap = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium'),
        eleTightIdMap = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight'),
        eleHEEPIdMap = cms.InputTag('egmGsfElectronIDs:heepElectronID-HEEPV70'),
        eleMVANonTrigMediumIdMap = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90'),
        eleMVANonTrigTightIdMap = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80'),
        eleMVATrigMediumIdMap = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90'), ### NOTE -> SAME AS NON-TRIG IN 2017
        eleMVATrigTightIdMap = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80'), ### NOTE -> SAME AS NON-TRIG IN 2017
        eleEcalRecHitCollection = cms.InputTag("reducedEgamma:reducedEBRecHits"),
        eleSingleTriggerIsoFileName = cms.string('data/SingleEleTriggerEff.root'),
        eleSingleTriggerFileName = cms.string('data/eleTriggerEff_MORIOND17.root'),
        eleVetoIdFileName = cms.string('data/eleVetoIDSF_MORIOND17.root'),
        eleLooseIdFileName = cms.string('data/eleLooseIDSF_MORIOND17.root'),
        eleMediumIdFileName = cms.string('data/eleMediumIDSF_MORIOND17.root'),
        eleTightIdFileName = cms.string('data/eleTightIDSF_MORIOND17.root'),
        eleMVATrigMediumIdFileName = cms.string('data/eleMVA90IDSF_MORIOND17.root'),
        eleMVATrigTightIdFileName = cms.string('data/eleMVA80IDSF_MORIOND17.root'),
        eleRecoEffFileName = cms.string('data/eleRecoSF_MORIOND17.root'),
        eleScaleSmearCorrectionName = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_23Jan_ele'),
        electron1id = cms.int32(3), # 0: veto, 1: loose, 2: medium, 3: tight, 4: HEEP, 5: MVA medium nonTrig, 6: MVA tight nonTrig, 7: MVA medium Trig, 8: MVA tight Trig
        electron2id = cms.int32(3),
        electron1pt = cms.double(30),
        electron2pt = cms.double(30),
    ),
    muonSet = cms.PSet(
        muons = cms.InputTag('cleanedMuons'),#('slimmedMuons'),#
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
        muonTrkFileName = cms.string('data/MuonTrkEfficienciesAndSF_MORIOND17.root'),
        muonIdFileName = cms.string('data/MuonIdEfficienciesAndSF_MORIOND17.root'),
        muonIsoFileName = cms.string('data/MuonIsoEfficienciesAndSF_MORIOND17.root'),
        muonTrkHighptFileName = cms.string('data/tkhighpt_2016full_absetapt.root'),
        muonTriggerFileName = cms.string('data/MuonTrigEfficienciesAndSF_MORIOND17.root'),
        doubleMuonTriggerFileName = cms.string('data/MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.root'),#FIXME -> obsolete
        muon1id = cms.int32(3), # 0: tracker high pt muon id, 1: loose, 2: medium, 3: tight, 4: high pt
        muon2id = cms.int32(3),
        muon1iso = cms.int32(2), # 0: trk iso (<0.1), 1: loose (<0.25), 2: tight (<0.15) (pfIso in cone 0.4)
        muon2iso = cms.int32(2),
        muon1pt = cms.double(27.),
        muon2pt = cms.double(27.),
        useTuneP = cms.bool(False),
        doRochester = cms.bool(False),
    ),
    tauSet = cms.PSet(
        taus = cms.InputTag('slimmedTaus'),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
        taupt = cms.double(18.),
        taueta = cms.double(2.3),
        tauIdByDecayMode = cms.int32(1),# 0: not set, 1: old, 2: new
        tauIdByDeltaBetaIso = cms.int32(1),# 0: not set, 1: loose, 2: medium, 3: tight
        tauIdByMVAIso = cms.int32(0),# 0: not set, 1: V loose, 2: loose, 3: medium, 4: tight, 5: V tight
        tauIdByMuonRejection = cms.int32(0),# 0: not set, 1: loose, 2: tight
        tauIdByElectronRejection = cms.int32(0),# 0: not set, 1: V loose, 2: loose, 3: medium, 4: tight
    ),
    photonSet = cms.PSet(
        photons = cms.InputTag('gedPhotons'),
        vertices = cms.InputTag('offlinePrimaryVertices'),
        phoLooseIdMap = cms.InputTag('egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose'),
        phoMediumIdMap = cms.InputTag('egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium'),
        phoTightIdMap = cms.InputTag('egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight'),
        phoMVANonTrigMediumIdMap = cms.InputTag('egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90'),
        phoEcalRecHitCollection = cms.InputTag("reducedEgamma:reducedEBRecHits"),
        phoLooseIdFileName = cms.string('data/phoLooseIDSF_MORIOND17.root'),
        phoMediumIdFileName = cms.string('data/phoMediumIDSF_MORIOND17.root'),
        phoTightIdFileName = cms.string('data/phoTightIDSF_MORIOND17.root'),
        phoMVANonTrigMediumIdFileName = cms.string('data/phoMVA90IDSF_MORIOND17.root'),
        photonid = cms.int32(1), # 1: loose, 2: medium, 3: tight, 4:MVA NonTrig medium
        photonpt = cms.double(15.),
    ),
    vertexSet = cms.PSet(
        primaryVertices = cms.InputTag('offlinePrimaryVertices'),
        secondaryVertices =  cms.InputTag('inclusiveCandidateSecondaryVertices'),
    ),
    pfCandidateSet = cms.PSet(
        pfCandidates = cms.InputTag('packedPFCandidates'),
        lostTracks = cms.InputTag('lostTracks'),
    ),
    dtSet = cms.PSet(
        dtsegments = cms.InputTag('dt4DSegments')
        ),
    cscSet = cms.PSet(
        cscsegments = cms.InputTag('cscSegments')
        ),
    standaloneMuonsSet = cms.PSet(
        standaloneMuons = cms.InputTag('standAloneMuons')
        ),
    displacedStandaloneMuonsSet = cms.PSet(
        standaloneMuons = cms.InputTag('displacedStandAloneMuons')
        ),
    minGenBpt = cms.double(0.),#(15.),#gen b quarks in acceptance
    maxGenBeta = cms.double(999.),#(2.4),#gen b quarks in acceptance
    minGenBradius2D = cms.double(129.),#new!! in cm
    maxGenBradius2D = cms.double(402.),#new!! in cm
    minGenBetaAcc = cms.double(0.),#(2.4),#
    maxGenBetaAcc = cms.double(1.1),#(2.4),#
    #invmassVBF = cms.double(400.?),#https://indico.desy.de/indico/event/20983/contribution/0/material/slides/0.pdf
    #new cut, motivated by calo-lifetimes trigger path
    invmassVBF = cms.double(250.),
    #detaVBF = cms.double(3.0?),#https://indico.desy.de/indico/event/20983/contribution/0/material/slides/0.pdf
    #new cut, motivated by calo-lifetimes trigger path
    detaVBF = cms.double(2.5),
    writeGenVBFquarks = cms.bool(True),
    writeGenHiggs = cms.bool(True),
    writeGenBquarks = cms.bool(True), #Acceptance cuts a few lines above!
    writeGenLLPs = cms.bool(True),
    writeOnlyTriggerEvents = cms.bool(True),#slims down ntuples a lot
    writeOnlyisVBFEvents = cms.bool(False),#slims down ntuples a lot
    ## PFCandidates:
    writeAK4JetPFCandidates = cms.bool(False), #Matched to AK4 only!
    writeAK8JetPFCandidates = cms.bool(False), #Matched to AK8 only!
    writeAllJetPFCandidates = cms.bool(False), #Matched to either AK4 or AK8
    writeAllPFCandidates = cms.bool(False), #All PFCandidates. Large collection: Please write only if needed!
    performPreFiringStudies = cms.bool(True if ('unprefirable' in process.source.fileNames[0]) else False),
    verbose = cms.bool(False),
)

process.demo = cms.EDAnalyzer('Demo'
)

process.seq = cms.Sequence(
    #process.packedPFCandidates *
    #process.packedCandsForTkIso *
    #process.lostTracksForTkIso *
    #process.lostTracks *
    #process.reducedEgamma *
    #process.heepIDVarValueMaps *
    #process.egmPhotonIDs *
    #process.egmGsfElectronIDs *
    #process.patPhotons *
    #process.patElectrons *
    #process.selectedPatPhotons *
    #process.selectedPatElectrons *
    ##process.slimmedPhotons *
    ##process.slimmedElectrons *
    ####process.fullPatMetSequenceTEST *#leading to segfault
    process.counter *
    process.ntuple
    #process.demo
)

process.p = cms.Path(process.seq)

outFile = open("tmpConfig_ZHNtuplizer.py","w")
outFile.write(process.dumpPython())
outFile.close()
