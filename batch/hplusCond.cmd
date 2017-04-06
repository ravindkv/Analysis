Universe              = vanilla
#Initialdir            = ../../../
Executable            = hplusAnalyzer.sh
Arguments             = FNAME OUTPUTFILE OUTPUTDIR $(Process) $(Cluster)
Output                = condor_out_$(Process)_$(Cluster).stdout
Error                 = condor_out_$(Process)_$(Cluster).stderr
Log                   = condor_out_$(Process)_$(Cluster).log
Notification          = Error
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
#transfer_input_files  = runMe.C, runMe.sh, hplusAnalyzer.C
Queue

