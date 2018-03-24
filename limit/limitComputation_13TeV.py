#//////////////////////////////////////////////////
#                                                 #
# 	Limit computation at 13 TeV 		  #
#                                                 #
#//////////////////////////////////////////////////

import os
import sys

def execme(command):
    #print "\033[01;32m"+ "Excecuting: "+ "\033[00m",  command
    os.system(command)
execme("mkdir -p ~/ARCHIVE/limit")
execme("cp *.C ~/ARCHIVE/limit")

#---------------------------------------------
#function to prepare data cards
#---------------------------------------------
def makeDataCards(CHANNEL_NAME, IN_FILE_DIR, HIST_NAME):
    execme("cp ~/ARCHIVE/limit/*.C .")
    execme("cp makeHPlusDataCard_13TeV.C ../")
    execme("cp template_datacard_csbar.txt ../")
    execme('sed -i '+'s#inFileDir_#'+IN_FILE_DIR+'# makeHPlusDataCard_13TeV.C')
    execme('sed -i '+'s#inShapeHisto_#'+HIST_NAME+'# makeHPlusDataCard_13TeV.C')
    if(CHANNEL_NAME=="mu"):
        execme('sed -i '+'s#isMuonChannel_#true# makeHPlusDataCard_13TeV.C')
        execme('sed -i '+'s#isEleChannel_#false# makeHPlusDataCard_13TeV.C')
        execme('cp '+muon_file_dir+'/all_muData.root '+muon_file_dir+'/all_Data.root')
        execme('sed -i '+'s#LUMI#'+str(muon_lumi)+'# makeHPlusDataCard_13TeV.C')
    if(CHANNEL_NAME=="ele"):
        execme('sed -i '+'s#isMuonChannel_#false# makeHPlusDataCard_13TeV.C')
        execme('sed -i '+'s#isEleChannel_#true# makeHPlusDataCard_13TeV.C')
        execme('cp '+ele_file_dir+'/all_EleData.root '+ele_file_dir+'/all_Data.root')
        execme('sed -i '+'s#LUMI#'+str(ele_lumi)+'# makeHPlusDataCard_13TeV.C')
    execme('sed -i '+'s#CHANNEL#'+CHANNEL_NAME+'# makeHPlusDataCard_13TeV.C')
    execme('sed -i '+'s#CHANNEL#'+CHANNEL_NAME+'# template_datacard_csbar.txt')
    execme('root -b -q makeHPlusDataCard_13TeV.C')
    execme("mv ../makeHPlusDataCard_13TeV.C . ")
    execme("mv ../template_datacard_csbar.txt . ")

#---------------------------------------------
#function to plot limits
#---------------------------------------------
def plotLimts(CHANNEL_NAME, HIST_NAME):
    execme("cp ~/ARCHIVE/limit/*.C .")
    execme("cp plotLimits_13TeV.C ../")
    execme('sed -i '+'s#HISTDIR#'+CHANNEL_NAME+'# plotLimits_13TeV.C')
    execme('sed -i '+'s#CHANNEL_HIST#'+CHANNEL_NAME+'_'+HIST_NAME+'# plotLimits_13TeV.C')
    execme('sed -i '+'s#HISTNAME#'+HIST_NAME+'# plotLimits_13TeV.C')
    execme('sed -i '+'s#CHANNELNAME#'+CHANNEL_NAME+'+jets# plotLimits_13TeV.C')
    execme('sed -i '+'s#CHANNELNAME#'+CHANNEL_NAME+'+jets# plotLimits_13TeV.C')
    if(CHANNEL_NAME=="mu"):
        execme('sed -i '+'s#LUMI#'+str(muon_lumi)+'# plotLimits_13TeV.C')
    if(CHANNEL_NAME=="ele"):
        execme('sed -i '+'s#LUMI#'+str(ele_lumi)+'# plotLimits_13TeV.C')
    execme('root -b -q plotLimits_13TeV.C')
    execme('mv limit_* '+CHANNEL_NAME+"/"+HIST_NAME)
    execme('cp *.C '+CHANNEL_NAME+"/"+HIST_NAME)
    execme("mv ../plotLimits_13TeV.C .")

#---------------------------------------------
#function to calculate limits
#---------------------------------------------
def calcLimits(CHANNEL_NAME, HIST_NAME, MASS_ARRAY):
    for MASS_POINT in MASS_ARRAY:
        DATACARD = 'datacard_csbar_'+CHANNEL_NAME+'_'+HIST_NAME+'_13TeV_mH'+str(MASS_POINT)
        execme('text2workspace.py '+DATACARD+'.txt -P HiggsAnalysis.CombinedLimit.ChargedHiggs:brChargedHiggs -o t2w_'+DATACARD+'.root')
        execme('combine --rAbsAcc 0.000001 t2w_'+DATACARD+'.root -M Asymptotic --mass '+str(MASS_POINT)+' --name ChargedHiggs_'+CHANNEL_NAME+'_'+HIST_NAME)
    #MOVE OUTPUT FILES TO THE LIMIT DIR
    LIMIT_DIR = CHANNEL_NAME+"/"+HIST_NAME
    execme('mkdir -p '+LIMIT_DIR)
    execme('cp *.C '+LIMIT_DIR)
    execme('mv t2w* '+LIMIT_DIR)
    execme('mv datacard_* '+LIMIT_DIR)
    execme('mv HplusShapes* '+LIMIT_DIR)
    execme('mv higgsCombine* '+LIMIT_DIR)

#---------------------------------------------------
#function to get datacards to be combined
#---------------------------------------------------
def getCardsToBeCombined(CHANNEL_ARRAY, IN_FILE_DIR_ARRAY, HIST_ARRAY, MASS_ARRAY):
    #make separate cards first
    for CH in range(len(CHANNEL_ARRAY)):
        for HIST in range(len(HIST_ARRAY)):
 	        makeDataCards(CHANNEL_ARRAY[CH], IN_FILE_DIR_ARRAY[CH], HIST_ARRAY[HIST])
    #store separate cards in an array
    COMB_CARD_CHANNEL_HIST_MASS = []
    for CH in range(len(CHANNEL_ARRAY)):
        COMB_CARD_HIST_MASS = []
        for HIST in range(len(HIST_ARRAY)):
            COMB_CARD_MASS = []
            for MASS in range(len(MASS_ARRAY)):
                COMB_CARD_MASS.append('datacard_csbar_'+CHANNEL_ARRAY[CH]+'_'+HIST_ARRAY[HIST]+'_13TeV_mH'+str(MASS_ARRAY[MASS])+'.txt')
            COMB_CARD_HIST_MASS.append(COMB_CARD_MASS)
        COMB_CARD_CHANNEL_HIST_MASS.append(COMB_CARD_HIST_MASS)    
    return COMB_CARD_CHANNEL_HIST_MASS

#---------------------------------------------------
#function to arrange datacards for combined limits
#---------------------------------------------------
def sortCardsForCombine(COMB_CARD_CHANNEL_HIST_MASS_ARRAY, CHANNEL_ARRAY, HIST_ARRAY, MASS):
    SORT_CARD = ' '
    COMB_CARD_MASS = []
    for CH in range(len(CHANNEL_ARRAY)):
        for HIST in range(len(HIST_ARRAY)):
            COMB_CARD_MASS.append(COMB_CARD_CHANNEL_HIST_MASS_ARRAY[CH][HIST][MASS])
    for STR in COMB_CARD_MASS:
        SORT_CARD = SORT_CARD+STR+' '
    return SORT_CARD

#---------------------------------------------------
#function to calculate combined limits
#---------------------------------------------------
def calcPlotCombinedLimits(CHANNEL_ARRAY, IN_FILE_DIR_ARRAY, HIST_ARRAY, MASS_ARRAY):
    COMB_CHANNEL_NAME = '_'.join(CHANNEL_ARRAY)
    COMB_HIST_NAME = '_'.join(HIST_ARRAY)
    getCardsToBeCombined_ = getCardsToBeCombined(CHANNEL_ARRAY, IN_FILE_DIR_ARRAY, HIST_ARRAY, MASS_ARRAY)
    for MASS in range(len(MASS_ARRAY)):
        sortCardsForCombine_ = sortCardsForCombine(getCardsToBeCombined_, CHANNEL_ARRAY, HIST_ARRAY, MASS)
        print sortCardsForCombine_
        COMB_DATACARD_NAME = 'datacard_csbar_'+COMB_CHANNEL_NAME+'_'+COMB_HIST_NAME+'_13TeV_mH'+str(MASS_ARRAY[MASS])
        if len(CHANNEL_ARRAY)>1 or len(HIST_ARRAY)>1:
            execme('combineCards.py '+sortCardsForCombine_+' > '+COMB_DATACARD_NAME+'.txt')
    calcLimits(COMB_CHANNEL_NAME, COMB_HIST_NAME, MASS_ARRAY)
    plotLimts(COMB_CHANNEL_NAME, COMB_HIST_NAME)

#---------------------------------------------
#USERS INPUTS
#---------------------------------------------
isSepLimit = False
isCombLimit = True
muon_lumi = 35.45
ele_lumi = 35.49
muon_file_dir="stack_20180116_Mu_sys"
ele_file_dir="stack_20180116_Ele_sys"
#channel_array = ["mu"]
#file_dir_array = [muon_file_dir]

channel_array = ["ele"]
file_dir_array = [ele_file_dir]

#channel_array = ["mu", "ele"]
#file_dir_array = [muon_file_dir, ele_file_dir]

hist_array = []
#hist_array.append("mjj_kfit")
hist_array.append("pt_bjetH")
hist_array.append("mjj_kfit_CTagL_SF_Cat")
hist_array.append("mjj_kfit_CTagM_SF_Cat")
hist_array.append("mjj_kfit_CTagT_SF_Cat")
mass_array = [80, 90, 100, 120, 140, 150, 155, 160]

###########################################
# calculate and plot separate limits
###########################################
if(isSepLimit):
    for shape in hist_array:
        #make datacards
        #makeDataCards("mu", muon_file_dir, shape)
        makeDataCards("ele", ele_file_dir, shape)
        
        #calc limits
        #calcLimits("mu",  shape, mass_array)
        calcLimits("ele", shape, mass_array)
        
        #plot limits
        plotLimts("ele", shape)

###########################################
# calculate, plot combined limits
###########################################
if(isCombLimit):
    calcPlotCombinedLimits(channel_array, file_dir_array, hist_array, mass_array)

execme("cp ~/ARCHIVE/limit/*.C .")
