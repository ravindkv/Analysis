
#//////////////////////////////////////////////////
#                                                 #
# 	Limit computation at 13 TeV 		  #
#                                                 #
#//////////////////////////////////////////////////

import os
import sys
import datetime

#USERS INPUTS
isMuChannel = True
isEleChannel= False
isLepChannel= False

muon_file_dir="stack_20171115_Mu_sys"
current_dir = os.getcwd()
shapeVar = ["mjj_kfit", "mjj_kfit_CTagL", "mjj_kfit_CTagM", "mjj_kfit_CTagT"]
massRange = [80, 90, 100, 120, 140, 150, 155, 160]
#-------------------------------
my_date = str(datetime.date.today()).replace("-","")
def execme(command):
    #print ""
    #print "\033[01;32m"+ "Excecuting: "+ "\033[00m",  command
    os.system(command)

isMakeDataCard= False
isCalculateLimit= False
isPlotLimit= True
if isMuChannel:
    muon_limit_dir = str("muon_limit_"+my_date)
    execme("mkdir -p "+muon_limit_dir)
    for shape in range(len(shapeVar)):
	shape_dir =shapeVar[shape]
        execme("mkdir -p "+muon_limit_dir+"/"+shape_dir)
        
	#---------------------------------------------
	#Prepare data cards
        #---------------------------------------------
	if(isMakeDataCard):
	    execme("cp makeHPlusDataCard_13TeV.C ../")
            execme("cp template_datacard_csbar_mu_8TeV.txt ../")
            execme('sed -i '+'s#inFileDir_#'+current_dir+"/"+muon_file_dir+'# makeHPlusDataCard_13TeV.C')
            execme('sed -i '+'s#inShapeHisto_#'+shape_dir+'# makeHPlusDataCard_13TeV.C')
            execme('sed -i '+'s#isMuonChannel_#true# makeHPlusDataCard_13TeV.C')
            execme('sed -i '+'s#isEleChannel_#false# makeHPlusDataCard_13TeV.C')
            execme('sed -i '+'s#CHANNEL#mu# makeHPlusDataCard_13TeV.C')
            execme('sed -i '+'s#PATH#'+current_dir+"/"+muon_limit_dir+"/"+shape_dir+'# template_datacard_csbar_mu_8TeV.txt')
            execme('root -l -q makeHPlusDataCard_13TeV.C')
            execme('mv datacard_csbar_* '+muon_limit_dir+"/"+shape_dir)
            execme('mv HplusShapes* '+muon_limit_dir+"/"+shape_dir)
            execme("mv ../makeHPlusDataCard_13TeV.C . ")
            execme("mv ../template_datacard_csbar_mu_8TeV.txt . ")
	
	#---------------------------------------------
	#Calculate limits
        #---------------------------------------------
	if(isCalculateLimit):
	    for mass in massRange:
	        execme('text2workspace.py '+muon_limit_dir+'/'+shape_dir+'/datacard_csbar_mu_'+shape_dir+'_13TeV_mH'+
				str(mass)+'.txt -P HiggsAnalysis.CombinedLimit.ChargedHiggs:brChargedHiggs -o comb_mu_'+shape_dir+'_mH'+str(mass)+'.root')
	        execme('combine comb_mu_'+shape_dir+'_mH'+str(mass)+'.root -M Asymptotic --mass '+str(mass)+' --name ChargedHiggs_mu_'+shape_dir)
            execme('mv comb* '+muon_limit_dir+"/"+shape_dir)
            execme('mv higgsCombine* '+muon_limit_dir+"/"+shape_dir)
	
	#---------------------------------------------
	#Plot limits
        #---------------------------------------------
	if(isPlotLimit):
	    execme("cp makeLimitPlot_13TeV.C ../")
            execme('sed -i '+'s#CHANNEL_HIST#mu_'+shape_dir+'# makeLimitPlot_13TeV.C')
            execme('sed -i '+'s#HISTDIR#'+muon_limit_dir+'# makeLimitPlot_13TeV.C')
            execme('sed -i '+'s#HISTNAME#'+shape_dir+'# makeLimitPlot_13TeV.C')
            execme('root -b -q makeLimitPlot_13TeV.C')
            execme('mv limit_* '+muon_limit_dir+"/"+shape_dir)
	    execme("mv ../makeLimitPlot_13TeV.C .")

if(isEleChannel):
    for shape in range(len(shapeVar)):
        execme("hadd -k "+str(data[samp])+"_Merged.root "+str(data[samp])+"*")
    execme("hadd -k all_EleData.root EleRun*_Merged.root")

