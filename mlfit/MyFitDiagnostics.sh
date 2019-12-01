t2wDataCard=$1
mass=$2
combine $t2wDataCard --mass $mass -M FitDiagnostics --plots --saveShapes --saveWithUncertainties  --saveNormalizations --initFromBonly --redefineSignalPOIs BR --setParameterRanges BR=0,0.1 | tee FitDiagnostics.log

#combine $t2wDataCard --mass $mass -M FitDiagnostics --plots --saveShapes --saveWithUncertainties --saveNormalizations --redefineSignalPOIs BR  --setParameterRanges BR=0,0.1  --robustFit=1 --setRobustFitAlgo=Minuit2 --setRobustFitStrategy=0 --setRobustFitTolerance=0.2 --saveOverallShapes 
#combine $t2wDataCard --mass $mass -M FitDiagnostics --saveNormalizations  --robustFit=1 --setRobustFitAlgo=Minuit2 --setRobustFitStrategy=0 --setRobustFitTolerance=0.2 --forceRecreateNLL 
#python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a fitDiagnostics.root -g fitDiag.pdf | tee diffNuisances.log
python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py --all fitDiagnostics.root  --poi=BR -g fitDiag.pdf | tee diffNuisances.log


