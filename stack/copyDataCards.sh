
#------------------------------------------
# Global paths
#------------------------------------------
limitPath="/home/rverma/combine/CMSSW_8_1_0/src/HiggsAnalysis/LimitComputationHtoCS/limit"
limitBinPath="/home/rverma/t3store2/condor_out/limit_condor/limit_out"

#------------------------------------------
# Copy data cards for muon + jets channel
#------------------------------------------
muCat1Path="mu/mjj_kfit"
muCat2Path="MuChannel/bJetPtCat"
muCat3Path="mu/mjj_kfit_CTagL_SF"
muCat4Path="MuChannel/CTagCat"
muCat5Path="MuChannel/BothCat"
echo "Copying muon data card ..."
mkdir -p MuChannel/Cat1_Inc
cp $limitPath/$muCat1Path/HplusShapes* MuChannel/Cat1_Inc 
cp $limitPath/$muCat1Path/*.txt MuChannel/Cat1_Inc 
mkdir -p MuChannel/Cat2_bJetPtBin
cp $limitBinPath/$muCat2Path/HplusShapes* MuChannel/Cat2_bJetPtBin
cp $limitBinPath/$muCat2Path/*.txt MuChannel/Cat2_bJetPtBin
mkdir -p MuChannel/Cat3_cTagInc
cp $limitPath/$muCat3Path/HplusShapes* MuChannel/Cat3_cTagInc 
cp $limitPath/$muCat3Path/*.txt MuChannel/Cat3_cTagInc 
mkdir -p MuChannel/Cat4_cTagEx
cp $limitBinPath/$muCat4Path/HplusShapes* MuChannel/Cat4_cTagEx
cp $limitBinPath/$muCat4Path/*.txt MuChannel/Cat4_cTagEx
mkdir -p MuChannel/Cat5_cTagEx_bJetPtBin
cp $limitBinPath/$muCat5Path/HplusShapes* MuChannel/Cat5_cTagEx_bJetPtBin
cp $limitBinPath/$muCat5Path/*.txt MuChannel/Cat5_cTagEx_bJetPtBin
echo ""

#------------------------------------------
# Copy data cards for electron + jets channel
#------------------------------------------
eleCat1Path="ele/mjj_kfit"
eleCat2Path="EleChannel/bJetPtCat"
eleCat3Path="ele/mjj_kfit_CTagL_SF"
eleCat4Path="EleChannel/CTagCat"
eleCat5Path="EleChannel/BothCat"
echo "Copying electron data card ..."
mkdir -p EleChannel/Cat1_Inc
cp $limitPath/$eleCat1Path/HplusShapes* EleChannel/Cat1_Inc 
cp $limitPath/$eleCat1Path/*.txt EleChannel/Cat1_Inc 
mkdir -p EleChannel/Cat2_bJetPtBin
cp $limitBinPath/$eleCat2Path/HplusShapes* EleChannel/Cat2_bJetPtBin
cp $limitBinPath/$eleCat2Path/*.txt EleChannel/Cat2_bJetPtBin
mkdir -p EleChannel/Cat3_cTagInc
cp $limitPath/$eleCat3Path/HplusShapes* EleChannel/Cat3_cTagInc 
cp $limitPath/$eleCat3Path/*.txt EleChannel/Cat3_cTagInc 
mkdir -p EleChannel/Cat4_cTagEx
cp $limitBinPath/$eleCat4Path/HplusShapes* EleChannel/Cat4_cTagEx
cp $limitBinPath/$eleCat4Path/*.txt EleChannel/Cat4_cTagEx
mkdir -p EleChannel/Cat5_cTagEx_bJetPtBin
cp $limitBinPath/$eleCat5Path/HplusShapes* EleChannel/Cat5_cTagEx_bJetPtBin
cp $limitBinPath/$eleCat5Path/*.txt EleChannel/Cat5_cTagEx_bJetPtBin
echo ""

#------------------------------------------
# Copy data cards for lepton + jets channel
#------------------------------------------
lepCat1Path="mu_ele/mjj_kfit"
lepCat2Path="LepChannel/bJetPtCat"
lepCat3Path="mu_ele/mjj_kfit_CTagL_SF"
lepCat4Path="LepChannel/CTagCat"
lepCat5Path="LepChannel/BothCat"
echo "Copying lepton data card ..."
mkdir -p LepChannel/Cat1_Inc
cp $limitPath/$lepCat1Path/HplusShapes* LepChannel/Cat1_Inc 
cp $limitPath/$lepCat1Path/*.txt LepChannel/Cat1_Inc 
mkdir -p LepChannel/Cat2_bJetPtBin
cp $limitBinPath/$lepCat2Path/HplusShapes* LepChannel/Cat2_bJetPtBin
cp $limitBinPath/$lepCat2Path/*.txt LepChannel/Cat2_bJetPtBin
mkdir -p LepChannel/Cat3_cTagInc
cp $limitPath/$lepCat3Path/HplusShapes* LepChannel/Cat3_cTagInc 
cp $limitPath/$lepCat3Path/*.txt LepChannel/Cat3_cTagInc 
mkdir -p LepChannel/Cat4_cTagEx
cp $limitBinPath/$lepCat4Path/HplusShapes* LepChannel/Cat4_cTagEx
cp $limitBinPath/$lepCat4Path/*.txt LepChannel/Cat4_cTagEx
mkdir -p LepChannel/Cat5_cTagEx_bJetPtBin
cp $limitBinPath/$lepCat5Path/HplusShapes* LepChannel/Cat5_cTagEx_bJetPtBin
cp $limitBinPath/$lepCat5Path/*.txt LepChannel/Cat5_cTagEx_bJetPtBin


#-------------------------------------------------------
# MuChannel: Clean
#-------------------------------------------------------
cd MuChannel/Cat1_Inc
mv datacard_csbar_mu_mjj_kfit_13TeV_mH90.txt datacard_csbar_muCat1_13TeV_mH90.txt
mv datacard_csbar_mu_mjj_kfit_13TeV_mH100.txt datacard_csbar_muCat1_13TeV_mH100.txt
mv datacard_csbar_mu_mjj_kfit_13TeV_mH120.txt datacard_csbar_muCat1_13TeV_mH120.txt
mv datacard_csbar_mu_mjj_kfit_13TeV_mH140.txt datacard_csbar_muCat1_13TeV_mH140.txt
mv datacard_csbar_mu_mjj_kfit_13TeV_mH150.txt datacard_csbar_muCat1_13TeV_mH150.txt
mv datacard_csbar_mu_mjj_kfit_13TeV_mH155.txt datacard_csbar_muCat1_13TeV_mH155.txt
mv datacard_csbar_mu_mjj_kfit_13TeV_mH160.txt datacard_csbar_muCat1_13TeV_mH160.txt
rm *mH80* && cd ../../

cd MuChannel/Cat2_bJetPtBin
mv combine_13TeV_mH90.txt  datacard_csbar_muCat2_13TeV_mH90.txt
mv combine_13TeV_mH100.txt datacard_csbar_muCat2_13TeV_mH100.txt
mv combine_13TeV_mH120.txt datacard_csbar_muCat2_13TeV_mH120.txt
mv combine_13TeV_mH140.txt datacard_csbar_muCat2_13TeV_mH140.txt
mv combine_13TeV_mH150.txt datacard_csbar_muCat2_13TeV_mH150.txt
mv combine_13TeV_mH155.txt datacard_csbar_muCat2_13TeV_mH155.txt
mv combine_13TeV_mH160.txt datacard_csbar_muCat2_13TeV_mH160.txt
rm *mH80* && rm datacard*mjj_kfit*.txt && cd ../../

cd MuChannel/Cat3_cTagInc
mv datacard_csbar_mu_mjj_kfit_CTagL_SF_13TeV_mH90.txt datacard_csbar_muCat3_13TeV_mH90.txt
mv datacard_csbar_mu_mjj_kfit_CTagL_SF_13TeV_mH100.txt datacard_csbar_muCat3_13TeV_mH100.txt
mv datacard_csbar_mu_mjj_kfit_CTagL_SF_13TeV_mH120.txt datacard_csbar_muCat3_13TeV_mH120.txt
mv datacard_csbar_mu_mjj_kfit_CTagL_SF_13TeV_mH140.txt datacard_csbar_muCat3_13TeV_mH140.txt
mv datacard_csbar_mu_mjj_kfit_CTagL_SF_13TeV_mH150.txt datacard_csbar_muCat3_13TeV_mH150.txt
mv datacard_csbar_mu_mjj_kfit_CTagL_SF_13TeV_mH155.txt datacard_csbar_muCat3_13TeV_mH155.txt
mv datacard_csbar_mu_mjj_kfit_CTagL_SF_13TeV_mH160.txt datacard_csbar_muCat3_13TeV_mH160.txt
rm *mH80* && cd ../../

cd MuChannel/Cat4_cTagEx
mv combine_13TeV_mH90.txt  datacard_csbar_muCat4_13TeV_mH90.txt
mv combine_13TeV_mH100.txt datacard_csbar_muCat4_13TeV_mH100.txt
mv combine_13TeV_mH120.txt datacard_csbar_muCat4_13TeV_mH120.txt
mv combine_13TeV_mH140.txt datacard_csbar_muCat4_13TeV_mH140.txt
mv combine_13TeV_mH150.txt datacard_csbar_muCat4_13TeV_mH150.txt
mv combine_13TeV_mH155.txt datacard_csbar_muCat4_13TeV_mH155.txt
mv combine_13TeV_mH160.txt datacard_csbar_muCat4_13TeV_mH160.txt
rm *mH80* && rm datacard*mjj_kfit*.txt && cd ../../

cd MuChannel/Cat5_cTagEx_bJetPtBin
mv combine_13TeV_mH90.txt  datacard_csbar_muCat5_13TeV_mH90.txt
mv combine_13TeV_mH100.txt datacard_csbar_muCat5_13TeV_mH100.txt
mv combine_13TeV_mH120.txt datacard_csbar_muCat5_13TeV_mH120.txt
mv combine_13TeV_mH140.txt datacard_csbar_muCat5_13TeV_mH140.txt
mv combine_13TeV_mH150.txt datacard_csbar_muCat5_13TeV_mH150.txt
mv combine_13TeV_mH155.txt datacard_csbar_muCat5_13TeV_mH155.txt
mv combine_13TeV_mH160.txt datacard_csbar_muCat5_13TeV_mH160.txt
rm *mH80* && rm datacard*mjj_kfit*.txt && cd ../../


#-------------------------------------------------------
# EleChannel: Clean
#-------------------------------------------------------
cd EleChannel/Cat1_Inc
mv datacard_csbar_ele_mjj_kfit_13TeV_mH90.txt datacard_csbar_eleCat1_13TeV_mH90.txt
mv datacard_csbar_ele_mjj_kfit_13TeV_mH100.txt datacard_csbar_eleCat1_13TeV_mH100.txt
mv datacard_csbar_ele_mjj_kfit_13TeV_mH120.txt datacard_csbar_eleCat1_13TeV_mH120.txt
mv datacard_csbar_ele_mjj_kfit_13TeV_mH140.txt datacard_csbar_eleCat1_13TeV_mH140.txt
mv datacard_csbar_ele_mjj_kfit_13TeV_mH150.txt datacard_csbar_eleCat1_13TeV_mH150.txt
mv datacard_csbar_ele_mjj_kfit_13TeV_mH155.txt datacard_csbar_eleCat1_13TeV_mH155.txt
mv datacard_csbar_ele_mjj_kfit_13TeV_mH160.txt datacard_csbar_eleCat1_13TeV_mH160.txt
rm *mH80* && cd ../../

cd EleChannel/Cat2_bJetPtBin
mv combine_13TeV_mH90.txt  datacard_csbar_eleCat2_13TeV_mH90.txt
mv combine_13TeV_mH100.txt datacard_csbar_eleCat2_13TeV_mH100.txt
mv combine_13TeV_mH120.txt datacard_csbar_eleCat2_13TeV_mH120.txt
mv combine_13TeV_mH140.txt datacard_csbar_eleCat2_13TeV_mH140.txt
mv combine_13TeV_mH150.txt datacard_csbar_eleCat2_13TeV_mH150.txt
mv combine_13TeV_mH155.txt datacard_csbar_eleCat2_13TeV_mH155.txt
mv combine_13TeV_mH160.txt datacard_csbar_eleCat2_13TeV_mH160.txt
rm *mH80* && rm datacard*mjj_kfit*.txt && cd ../../

cd EleChannel/Cat3_cTagInc
mv datacard_csbar_ele_mjj_kfit_CTagL_SF_13TeV_mH90.txt datacard_csbar_eleCat3_13TeV_mH90.txt
mv datacard_csbar_ele_mjj_kfit_CTagL_SF_13TeV_mH100.txt datacard_csbar_eleCat3_13TeV_mH100.txt
mv datacard_csbar_ele_mjj_kfit_CTagL_SF_13TeV_mH120.txt datacard_csbar_eleCat3_13TeV_mH120.txt
mv datacard_csbar_ele_mjj_kfit_CTagL_SF_13TeV_mH140.txt datacard_csbar_eleCat3_13TeV_mH140.txt
mv datacard_csbar_ele_mjj_kfit_CTagL_SF_13TeV_mH150.txt datacard_csbar_eleCat3_13TeV_mH150.txt
mv datacard_csbar_ele_mjj_kfit_CTagL_SF_13TeV_mH155.txt datacard_csbar_eleCat3_13TeV_mH155.txt
mv datacard_csbar_ele_mjj_kfit_CTagL_SF_13TeV_mH160.txt datacard_csbar_eleCat3_13TeV_mH160.txt
rm *mH80* && cd ../../

cd EleChannel/Cat4_cTagEx
mv combine_13TeV_mH90.txt  datacard_csbar_eleCat4_13TeV_mH90.txt
mv combine_13TeV_mH100.txt datacard_csbar_eleCat4_13TeV_mH100.txt
mv combine_13TeV_mH120.txt datacard_csbar_eleCat4_13TeV_mH120.txt
mv combine_13TeV_mH140.txt datacard_csbar_eleCat4_13TeV_mH140.txt
mv combine_13TeV_mH150.txt datacard_csbar_eleCat4_13TeV_mH150.txt
mv combine_13TeV_mH155.txt datacard_csbar_eleCat4_13TeV_mH155.txt
mv combine_13TeV_mH160.txt datacard_csbar_eleCat4_13TeV_mH160.txt
rm *mH80* && rm datacard*mjj_kfit*.txt && cd ../../

cd EleChannel/Cat5_cTagEx_bJetPtBin
mv combine_13TeV_mH90.txt  datacard_csbar_eleCat5_13TeV_mH90.txt
mv combine_13TeV_mH100.txt datacard_csbar_eleCat5_13TeV_mH100.txt
mv combine_13TeV_mH120.txt datacard_csbar_eleCat5_13TeV_mH120.txt
mv combine_13TeV_mH140.txt datacard_csbar_eleCat5_13TeV_mH140.txt
mv combine_13TeV_mH150.txt datacard_csbar_eleCat5_13TeV_mH150.txt
mv combine_13TeV_mH155.txt datacard_csbar_eleCat5_13TeV_mH155.txt
mv combine_13TeV_mH160.txt datacard_csbar_eleCat5_13TeV_mH160.txt
rm *mH80* && rm datacard*mjj_kfit*.txt && cd ../../


#-------------------------------------------------------
# LepChannel: Clean 
#-------------------------------------------------------
cd LepChannel/Cat1_Inc
mv datacard_csbar_mu_ele_mjj_kfit_13TeV_mH90.txt datacard_csbar_lepCat1_13TeV_mH90.txt
mv datacard_csbar_mu_ele_mjj_kfit_13TeV_mH100.txt datacard_csbar_lepCat1_13TeV_mH100.txt
mv datacard_csbar_mu_ele_mjj_kfit_13TeV_mH120.txt datacard_csbar_lepCat1_13TeV_mH120.txt
mv datacard_csbar_mu_ele_mjj_kfit_13TeV_mH140.txt datacard_csbar_lepCat1_13TeV_mH140.txt
mv datacard_csbar_mu_ele_mjj_kfit_13TeV_mH150.txt datacard_csbar_lepCat1_13TeV_mH150.txt
mv datacard_csbar_mu_ele_mjj_kfit_13TeV_mH155.txt datacard_csbar_lepCat1_13TeV_mH155.txt
mv datacard_csbar_mu_ele_mjj_kfit_13TeV_mH160.txt datacard_csbar_lepCat1_13TeV_mH160.txt
rm *mH80* && rm datacard_*mjj_kfit*.txt && cd ../../

cd LepChannel/Cat2_bJetPtBin
mv combine_13TeV_mH90.txt  datacard_csbar_lepCat2_13TeV_mH90.txt
mv combine_13TeV_mH100.txt datacard_csbar_lepCat2_13TeV_mH100.txt
mv combine_13TeV_mH120.txt datacard_csbar_lepCat2_13TeV_mH120.txt
mv combine_13TeV_mH140.txt datacard_csbar_lepCat2_13TeV_mH140.txt
mv combine_13TeV_mH150.txt datacard_csbar_lepCat2_13TeV_mH150.txt
mv combine_13TeV_mH155.txt datacard_csbar_lepCat2_13TeV_mH155.txt
mv combine_13TeV_mH160.txt datacard_csbar_lepCat2_13TeV_mH160.txt
rm *mH80* && rm datacard*mjj_kfit*.txt && cd ../../

cd LepChannel/Cat3_cTagInc
mv datacard_csbar_mu_ele_mjj_kfit_CTagL_SF_13TeV_mH90.txt datacard_csbar_lepCat3_13TeV_mH90.txt
mv datacard_csbar_mu_ele_mjj_kfit_CTagL_SF_13TeV_mH100.txt datacard_csbar_lepCat3_13TeV_mH100.txt
mv datacard_csbar_mu_ele_mjj_kfit_CTagL_SF_13TeV_mH120.txt datacard_csbar_lepCat3_13TeV_mH120.txt
mv datacard_csbar_mu_ele_mjj_kfit_CTagL_SF_13TeV_mH140.txt datacard_csbar_lepCat3_13TeV_mH140.txt
mv datacard_csbar_mu_ele_mjj_kfit_CTagL_SF_13TeV_mH150.txt datacard_csbar_lepCat3_13TeV_mH150.txt
mv datacard_csbar_mu_ele_mjj_kfit_CTagL_SF_13TeV_mH155.txt datacard_csbar_lepCat3_13TeV_mH155.txt
mv datacard_csbar_mu_ele_mjj_kfit_CTagL_SF_13TeV_mH160.txt datacard_csbar_lepCat3_13TeV_mH160.txt
rm *mH80* && rm datacard_*mjj_kfit*.txt && cd ../../

cd LepChannel/Cat4_cTagEx
mv combine_13TeV_mH90.txt  datacard_csbar_lepCat4_13TeV_mH90.txt
mv combine_13TeV_mH100.txt datacard_csbar_lepCat4_13TeV_mH100.txt
mv combine_13TeV_mH120.txt datacard_csbar_lepCat4_13TeV_mH120.txt
mv combine_13TeV_mH140.txt datacard_csbar_lepCat4_13TeV_mH140.txt
mv combine_13TeV_mH150.txt datacard_csbar_lepCat4_13TeV_mH150.txt
mv combine_13TeV_mH155.txt datacard_csbar_lepCat4_13TeV_mH155.txt
mv combine_13TeV_mH160.txt datacard_csbar_lepCat4_13TeV_mH160.txt
rm *mH80* && rm datacard*mjj_kfit*.txt && cd ../../

cd LepChannel/Cat5_cTagEx_bJetPtBin
mv combine_13TeV_mH90.txt  datacard_csbar_lepCat5_13TeV_mH90.txt
mv combine_13TeV_mH100.txt datacard_csbar_lepCat5_13TeV_mH100.txt
mv combine_13TeV_mH120.txt datacard_csbar_lepCat5_13TeV_mH120.txt
mv combine_13TeV_mH140.txt datacard_csbar_lepCat5_13TeV_mH140.txt
mv combine_13TeV_mH150.txt datacard_csbar_lepCat5_13TeV_mH150.txt
mv combine_13TeV_mH155.txt datacard_csbar_lepCat5_13TeV_mH155.txt
mv combine_13TeV_mH160.txt datacard_csbar_lepCat5_13TeV_mH160.txt
rm *mH80* && rm datacard*mjj_kfit*.txt && cd ../../

