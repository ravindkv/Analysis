#pushd /home/gouranga/Analysis/ChargedHiggs/8TeV/JERStusy/CMSSW_5_3_5/src
#setcms
#eval `scramv1 runtime -csh`
#popd
export MY_PATH=${PWD}/src/

#if[ $?LD_LIBRARY_PATH ] 
#then
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MY_PATH}
#else
#  export LD_LIBRARY_PATH=${MY_PATH}
#fi

