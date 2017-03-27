#include "TString.h"
#include "TROOT.h"

void runme(char *arg){
  gROOT->ProcessLine(TString::Format(".L %s.C+",arg));
  gROOT->ProcessLine(TString::Format("%s t",arg));
  gROOT->ProcessLine("t.processEvents()");
}
