#Here: standard crab config file

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys
config = config()

config.General.workArea = 'crab_projects_LLP'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/AODNtuplizer.py'
config.JobType.inputFiles = ['data']
#config.JobType.pyCfgParams = [string_runLocal, string_isData, string_isREHLT, string_isReReco, string_isReMiniAod, string_isPromptReco,string_noLHEinfo, string_isbbH, string_GT, string_JECstring, string_jsonName, string_triggerTag, string_filterString]

config.General.requestName = 'test'

config.Data.inputDataset =  '/VBFH_HToSSTobbbb_MH-125_MS-15_ctauS-5000_TuneCUETP8M1_13TeV-powheg-pythia8_Tranche2_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_Tranche2_AODSIM-b1a4edca9adfa7a2e4059536bf605cd7/USER'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'

config.Data.unitsPerJob = 10000
config.Data.outLFNDirBase = '/store/user/lbenato/choose_a_folder_name'
config.Data.publication = False

config.Site.storageSite = 'T2_DE_DESY'
config.Site.blacklist   = ['T2_FR_IPHC']
            

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)


    # Selection of samples via python lists
    import os
    from Analyzer.LLPonAOD.samples import sample, samples

    list_of_samples = ["SM_Higgs","VV","WJetsToQQ","WJetsToLNu","WJetsToLNu_Pt","DYJetsToQQ","DYJetsToNuNu","DYJetsToLL","ST","TTbar","QCD","signal_VBF","signal_ggH","all","data_obs","ZJetsToNuNu", "DYJets", "WJets", "signal_ZH"]#,"data_obs"
    print "Possible subgroups of samples:"
    for a in list_of_samples:
        print a
    print "---------------"

    ########parser#######
    import optparse
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option("-a", "--crabaction", action="store", type="string", dest="crabaction", default="test")
    parser.add_option("-l", "--lists", action="store", type="string", dest="lists", default="")
    parser.add_option("-g", "--groupofsamples", action="store", type="string", dest="groupofsamples", default="")
    parser.add_option("-c", "--calo", action="store_true", dest="calo", default=False)
    (options, args) = parser.parse_args()


    ####################
    #crabConfig = ''
    if options.calo:
       isCalo=True
    else:
       isCalo=False

    folder = ''
    pset = ''
    workarea = ''
    if  options.lists == "v0_calo_AOD":
        from Analyzer.LLPonAOD.crab_requests_lists_AOD_calo import * #This list is fine for us!
        #crabConfig = 'crabConfig.py'
        pset = "AODNtuplizer.py"
        folder = "v0_calo_AOD" #CHANGE here your crab folder name
        outLFNDirBase = "/store/user/lbenato/"+folder #CHANGE here according to your username!
        workarea = "/nfs/dust/cms/user/lbenato/" + folder #CHANGE here according to your username!
        isCalo=True
    else:
        print "No list indicated, aborting!"
        exit()

    selected_requests = {}
    if options.groupofsamples not in list_of_samples:
        print "Invalid subgroup of samples, aborting!"
        exit()

    for b, k in enumerate(requests.keys()):
        if options.groupofsamples=="signal_VBF":
            if "VBFH_HToSSTobb" in k:
                print k
                selected_requests[k] = requests[k]
        elif options.groupofsamples=="signal_ggH":
            if "GluGluH_HToSSTobb" in k:
                print k
                selected_requests[k] = requests[k]
        elif options.groupofsamples=="signal_ZH":
            if "ZH_HToSSTobb" in k:
                print k
                selected_requests[k] = requests[k]
        elif options.groupofsamples=="all":
            print "All samples considered"
            selected_requests[k] = requests[k]
        else:
            if k in samples[options.groupofsamples]["files"]:
                print k
                selected_requests[k] = requests[k]

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
    
    for a, j in enumerate(selected_requests):
        print "#"*65
	print "Dataset: ", j
        # Here: determines every needed parameter as per config file
        isData = True if ('SingleMuon' in j or 'SingleElectron' in j or 'JetHT' in j or 'BTagCSV' in j or 'DisplacedJet' in j or 'MET' in j) else False
        print "isData?", isData
        isReHLT = False
        isReReco          = True if ('23Sep2016' in j) else False
        isReMiniAod       = True if ('03Feb2017' in j) else False
        isPromptReco      = True if ('PromptReco' in j) else False
        theRunBCD = ['Run2016B','Run2016C','Run2016D']    
        theRunEF  = ['Run2016E','Run2016F']
        theRunG   = ['Run2016G']
        theRunH   = ['Run2016H']
        noLHEinfo = True if ('WW_TuneCUETP8M1_13TeV-pythia8' in j or 'WZ_TuneCUETP8M1_13TeV-pythia8' in j or 'ZZ_TuneCUETP8M1_13TeV-pythia8' in j) else False #check for PythiaLO samples
        isbbH = True if ('bbHToBB_M-125_4FS_yb2_13TeV_amcatnlo' in j) else False #bbH has a different label in LHEEventProduct
        isSignal = True if ('HToSSTobbbb_MH-125' in j) else False
        GT = ''
        if isData:
            if isReMiniAod and any(s in j for s in theRunH): GT = '80X_dataRun2_Prompt_v16' 
            else: GT = '80X_dataRun2_2016SeptRepro_v7'
        elif not(isData):                                       
            GT = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'#Moriond17 GT
        print "GT ->", GT

        JECstring = ''
        if isData:# and (isReReco or isReMiniAod):
          if any(s in j for s in theRunBCD):
            JECstring = "Summer16_23Sep2016BCDV3_DATA" #if isReMiniAod else "Summer16_23Sep2016BCDV3_DATA"
          if any(s in j for s in theRunEF):
            JECstring = "Summer16_23Sep2016EFV3_DATA" #if isReMiniAod else "Summer16_23Sep2016EFV3_DATA"
          if any(s in j for s in theRunG):
            JECstring = "Summer16_23Sep2016GV3_DATA" #if isReMiniAod else "Summer16_23Sep2016GV3_DATA"
          if any(s in j for s in theRunH):
            JECstring = "Summer16_23Sep2016HV3_DATA" #if isReMiniAod else "Summer16_23Sep2016HV3_DATA"
        elif isData and isPromptReco:
           JECstring = "Spring16_25nsV6_DATA"
        elif not isData:
           JECstring = "Summer16_23Sep2016V3_MC"

        print "JEC ->",JECstring

        # JSON filter
        jsonName = "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON"

        # Trigger filter
        triggerTag = 'HLT2' if isReHLT else 'HLT'

        filterString = "RECO"
        #if isData:
        #    filterString = "RECO"
        #else:
        #    filterString = "PAT"

        #Prepare inputstrings for pyCfg
        string_runLocal = 'runLocal=False'
        string_isData = 'PisData='+str(isData)
        string_isREHLT = 'PisReHLT='+str(isReHLT)
        string_isReReco = 'PisReReco='+str(isReReco)
        string_isReMiniAod = 'PisReMiniAod='+str(isReMiniAod)
        string_isPromptReco = 'PisPromptReco='+str(isPromptReco)
        string_noLHEinfo = 'PnoLHEinfo='+str(noLHEinfo)
        string_isbbH = 'PisbbH='+str(isbbH)
        string_isSignal = 'PisSignal='+str(isSignal)
        string_GT = 'PGT='+str(GT)
        string_JECstring = 'PJECstring='+str(JECstring)
        string_jsonName = 'PjsonName='+str(jsonName)
        string_triggerTag = 'PtriggerTag='+str(triggerTag)
        string_filterString = 'PfilterString='+str(filterString)
        string_calo = 'Pcalo=True' if isCalo else 'Pcalo=False'


        # submission of the python config
        if options.crabaction=="submit":
            if "VBFH_HToSS" in j:
                #automatic implementation of the choice bewteen inputDBS global/phys03
                config.Data.inputDBS = "phys03"
            elif "GluGluH_HToSS" in j:
                #automatic implementation of the choice bewteen inputDBS global/phys03
                config.Data.inputDBS = "phys03"
            else:
                config.Data.inputDBS = "global"

            os.system('echo submitting this config...\n')
            #modify parameters here
            config.General.requestName = j
            config.Data.inputDataset = selected_requests[j]
            config.JobType.psetName = "python/" + pset
            config.Data.outLFNDirBase = outLFNDirBase
            config.General.workArea= workarea
            if isData:
                config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
                #config.Data.splitting = 'Automatic'
                config.Data.unitsPerJob = 100000
            #config.JobType.pyCfgParams = ['runLocal=False']
            config.JobType.pyCfgParams = [string_runLocal, string_isData, string_isREHLT, string_isReReco, string_isReMiniAod, string_isPromptReco,string_noLHEinfo, string_isbbH, string_isSignal, string_GT, string_JECstring, string_jsonName, string_triggerTag, string_filterString, string_calo]
            print config
            submit(config)

        elif options.crabaction=="status":
            os.system('echo status -d ' + workarea + '/crab_'+j+'\n')
            os.system('crab status -d ' + workarea + '/crab_'+j+'\n')
            os.system('echo ----------------------------------------------------\n') 
        elif options.crabaction=="resubmit":
            os.system('echo resubmit -d ' + workarea + '/crab_'+j+'\n')
            os.system('crab resubmit -d ' + workarea + '/crab_'+j+'\n')
        elif options.crabaction=="getoutput":
            os.system('echo getoutput -d ' + workarea + '/crab_'+j+'\n')
            os.system('crab getoutput -d ' + workarea + '/crab_'+j+'\n')
        elif options.crabaction=="kill":
            os.system('echo kill -d ' + workarea + '/crab_'+j+'\n')
            os.system('crab kill -d ' + workarea + '/crab_'+j+'\n')
        elif options.crabaction=="report":
            os.system('echo report -d ' + workarea + '/crab_'+j+'\n')
            os.system('crab report -d ' + workarea + '/crab_'+j+'\n')
        elif options.crabaction=="test":
            if "VBFH_HToSS" in j:
                #automatic implementation of the choice bewteen inputDBS global/phys03
                config.Data.inputDBS = "phys03"
            else:
                config.Data.inputDBS = "global"

            os.system('echo submitting this config...\n')
            #modify parameters here
            config.General.requestName = j
            config.Data.inputDataset = selected_requests[j]
            config.JobType.psetName = "python/" + pset
            config.Data.outLFNDirBase = outLFNDirBase
            config.General.workArea= workarea
            if isData:
                config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
                #config.Data.splitting = 'Automatic'
                config.Data.unitsPerJob = 100000
            config.JobType.pyCfgParams = [string_runLocal, string_isData, string_isREHLT, string_isReReco, string_isReMiniAod, string_isPromptReco,string_noLHEinfo, string_isbbH, string_isSignal, string_GT, string_JECstring, string_jsonName, string_triggerTag, string_filterString, string_calo]
            print config
        else:
            print "Invalid crab action. Please type: -a submit/status/resubmit/getoutput/kill"
            exit()
    os.system('echo -%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-\n') 





