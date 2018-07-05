#//////////////////////////////////////////////////
#                                                 #
# Maximum likelyhood fit at 13 TeV 		  #
#                                                 #
#//////////////////////////////////////////////////

import os
import sys
import os.path

def execme(command):
    #print "\033[01;32m"+ "Excecuting: "+ "\033[00m",  command
    os.system(command)

#------------------------------------------------
#function to 'copy text to workspace' root files
#------------------------------------------------
def copyWSRootFile(LIMIT_DIR, CHANNEL_NAME, HIST_NAME, MASS):
    fullPath = LIMIT_DIR+"/"+CHANNEL_NAME+"/"+HIST_NAME
    fileName = 't2w_datacard_csbar_'+CHANNEL_NAME+'_'+HIST_NAME+'_13TeV_mH'+str(MASS)+'.root'
    execme('mkdir -p '+fullPath)
    filePath = '../'+fullPath+'/'+fileName
    if(os.path.exists(filePath)):
        #print 'copyWSRootFile: ' + filePath
        execme('cp '+filePath+' '+fullPath)
    else:
        print "ERROR: "+filePath+' does not exit !!!'

#---------------------------------------------
#function to do maximum likelyhood fit
#---------------------------------------------
def doMaxLikeFit(LIMIT_DIR, CHANNEL_NAME, HIST_NAME, MASS):
    copyWSRootFile(LIMIT_DIR, CHANNEL_NAME, HIST_NAME, MASS)
    fullPath = LIMIT_DIR+"/"+CHANNEL_NAME+"/"+HIST_NAME
    fileName = 't2w_datacard_csbar_'+CHANNEL_NAME+'_'+HIST_NAME+'_13TeV_mH'+str(MASS)+'.root'
    combineLogFile = fullPath+"/log_fitDiagnostics_"+CHANNEL_NAME+'_'+HIST_NAME+'_13TeV_mH'+str(MASS)+'.txt'
    fitDiagnosticsFile = fullPath+"/fitDiagnostics_"+CHANNEL_NAME+'_'+HIST_NAME+'_13TeV_mH'+str(MASS)+'.root'
    filePath = fullPath+'/'+fileName
    if(os.path.exists(filePath)):
        #print 'doMaxLikeFit on: ' + filePath
        #execme('combine '+fullPath+'/'+fileName+' -M MaxLikelihoodFit -t -1 --expectSignal 0 >'+combineLogFile)
        execme('combine '+fullPath+'/'+fileName+' -M FitDiagnostics -t -1 --expectSignal 0 >'+combineLogFile)
        execme('mv fitDiagnostics.root '+fitDiagnosticsFile)
    else:
        print "ERROR: "+filePath+' does not exit !!!'

#---------------------------------------------
#function to get the nuisance parameters
#---------------------------------------------
def diffNuisances(LIMIT_DIR, CHANNEL_NAME, HIST_NAME, MASS):
    print "\n/------------------------------------------------------------/"
    print "CHANNEL: %s, HIST: %s, MASS: %s" %(CHANNEL_NAME,HIST_NAME,MASS)
    print "/------------------------------------------------------------/"
    doMaxLikeFit(LIMIT_DIR, CHANNEL_NAME, HIST_NAME, MASS)
    fullPath = LIMIT_DIR+"/"+CHANNEL_NAME+"/"+HIST_NAME
    fitDiagnosticsFile = fullPath+"/fitDiagnostics_"+CHANNEL_NAME+'_'+HIST_NAME+'_13TeV_mH'+str(MASS)+'.root'
    diffNuFilePath = '${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py'
    diffNuLogFile = fullPath+"/log_diffNuisances_"+CHANNEL_NAME+'_'+HIST_NAME+'_13TeV_mH'+str(MASS)+'.txt'
    diffNuOutFile = fullPath+"/out_diffNuisances_"+CHANNEL_NAME+'_'+HIST_NAME+'_13TeV_mH'+str(MASS)+'.root'
    if(os.path.exists(fitDiagnosticsFile)):
        #print 'diffNuisances: '+fitDiagnosticsFile
        fitParam = "nuisancs"
        #fitParam = "post_fit_errs"
        #fitParam = "asdf"
        execme('python '+diffNuFilePath+' -a '+fitDiagnosticsFile+' --poi=BR --all -g'+ diffNuOutFile+' >'+diffNuLogFile)
        execme('root -l -q \"plotNuisanceParam_13TeV.C(\\\"'+diffNuOutFile+'\\\", \\\"'+fitParam+'\\\",\\\"'+CHANNEL_NAME+'\\\",\\\"'+HIST_NAME+'\\\",\\\"'+str(MASS)+'\\\")\"')
    else:
        print "ERROR: "+diffNuFilePath+' does not exit !!!'
    
#---------------------------------------------
#USERS INPUTS
#---------------------------------------------
mass_array = [80]
#mass_array = [80, 90, 100, 120, 140, 150, 155, 160]

###########################################
# do ML fit, get nuisance parameters
###########################################
for mass in mass_array:
    diffNuisances("limit", "mu", "mjj_kfit", mass)
    '''
    diffNuisances("limit", "ele", "mjj_kfit", mass)
    diffNuisances("limit", "mu_ele", "mjj_kfit", mass)

    diffNuisances("limit", "mu", "mjj_kfit_CTagL_SF", mass)
    diffNuisances("limit", "ele", "mjj_kfit_CTagL_SF", mass)
    diffNuisances("limit", "mu_ele", "mjj_kfit_CTagL_SF", mass)
    diffNuisances("limit", "mu", "mjj_kfit_CTagM_SF", mass)
    diffNuisances("limit", "ele", "mjj_kfit_CTagM_SF", mass)
    diffNuisances("limit", "mu_ele", "mjj_kfit_CTagM_SF", mass)

    diffNuisances("limit", "mu", "mjj_kfit_CTagT_SF", mass)
    diffNuisances("limit", "ele", "mjj_kfit_CTagT_SF", mass)
    diffNuisances("limit", "mu_ele", "mjj_kfit_CTagT_SF", mass)
    '''
