# Analysis
   
### Download and compile the package  ###  
* cmsrel CMSSW_7_2_3
* cd CMSSW_7_2_3/src
* cmsenv
* git clone https://github.com/ravindkv/Analysis.git 
* cd Analysis/src
* make clean 
* make
* cd .. 

### Compile and run, in one go ### 
* root -l 'runme.C("hplusAnalyzer")'

### Submit condor batch jobs  ###

* voms-proxy-init -voms cms
* cd condor
* ./hplusRunCond.sh
