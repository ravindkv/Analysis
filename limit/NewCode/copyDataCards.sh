
#------------------------------------------
# Global paths
#------------------------------------------
limitPath="/home/rverma/t3store/AN-18-061/ExclusionLimit/CMSSW_8_0_25/src/HiggsAnalysis/HplusTocs13TeVLimit/Limit"
cp -r $limitPath/limit .
cd limit 
rm -r png_file/
rm -r root_file/
rm *.tex
rm */*/*.png
rm */*/*.root

rm */*/*/*.C
rm */*/*/*.h
rm */*/*/*.py
rm */*/*/MyTe*.txt
rm */*/*/data*.txt
rm */*/*/higg*.root

