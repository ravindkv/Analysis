#include <iostream>
#include <fstream>
#include <iomanip>

//hadd all_MC.root all_DY.root all_ST.root all_TTJetsP.root all_VV.root all_WJets.root

TH1F* getHisto(TFile *histFile, TString histPath){
  TH1F* hist; 
  TString inFile("$PWD/");
  TFile *fTT= new TFile(inFile+"all_TTJetsP.root"); 
  if(!(histFile->Get(histPath))){
    hist = (TH1F*)(fTT->Get(histPath));
    hist->Add(hist, -1);
  }else hist = (TH1F*)(histFile->Get(histPath));
  return hist;
}

void getABCDNumbers(ofstream &outFile, TFile *f, string proc, TString A, TString B, TString C, TString D){
  TH1F * hA = getHisto(f, A);
  TH1F * hB = getHisto(f, B);
  TH1F * hC = getHisto(f, C);
  TH1F * hD = getHisto(f, D);
  outFile<< proc<<" & "<< hA->Integral()<<" & "<< hB->Integral()<<" & "<< hC->Integral()<<" & "<< hD->Integral()<<" \\\\ "<<endl;
}

void makeABCDTable_13TeV(){  
  TString inFile("$PWD/");
  TFile *ttbar    		= new TFile(inFile+"all_TTJetsP.root"); 
  TFile *wjet  			= new TFile(inFile+"all_WJets.root"); 
  TFile *zjet  			= new TFile(inFile+"all_DY.root");
  TFile *stop  			= new TFile(inFile+"all_ST.root");
  TFile *diboson 		= new TFile(inFile+"all_VV.root");
  TFile *allMC 			= new TFile(inFile+"all_MC.root");
  
  //TFile *data = new TFile(inFile+"all_muData.root");
  TFile *data = new TFile(inFile+"all_EleData.root");
  ofstream outFile; 
  outFile.open("qcdABCDTable.tex"); 
  //outFile<<"\\documentclass[landscape,letterpaper]{article}"<<endl;  
  outFile<<"\\documentclass[landscape]{article}"<<endl;  
  outFile<<"\\pagestyle{empty}"<<endl;  
  outFile<<"\\usepackage{epsfig}"<<endl;  
  outFile<<"\\usepackage{amsmath}"<<endl;  
  outFile<<"\\usepackage{array}"<<endl;  
  outFile<<"\\usepackage{multirow}"<<endl;  
  outFile<<"\\usepackage[cm]{fullpage}"<<endl;  
  outFile<<"\\textheight = 8.in"<<endl;  
  outFile<<"\\textwidth 7.0in"<<endl;  
  outFile<<"\\begin{document}"<<endl;  
  outFile<<""<<endl;
  outFile<<"\\begin{table}"<<endl; 
  outFile<<"\\begin{center}"<<endl; 
  outFile<<"\\begin{tabular}{ |c|c|c|c|c| }"<<endl; 
  outFile<<"\\multicolumn{5}{c}{ } \\\\"<<endl; 
  outFile<<"\\hline "<<endl;
  outFile<< "\\bf{Process}"<<" & "<< "Region-A & "<< "Region-B & "<< "Region-C & "<< "Region-D "<<" \\\\ "<<endl;
  outFile<<"" <<" & "<< "(Iso, high MET) & "<< "(Non iso, high MET) & "<< "(Non iso, low MET)& "<< "(Iso, low MET)"<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  //Add another table with JESUP
   getABCDNumbers(outFile, ttbar, "$t\\bar{t}$ + jets", "base/Iso/KinFit/mjj_kfit", "base/NonIso/KinFit/mjj_kfit", "baseLowMET/NonIso/KinFit/mjj_kfit", "baseLowMET/Iso/KinFit/mjj_kfit");
  outFile<<"\\hline "<<endl;
  
  getABCDNumbers(outFile, stop, "Single ~t",            "base/Iso/KinFit/mjj_kfit", "base/NonIso/KinFit/mjj_kfit", "baseLowMET/NonIso/KinFit/mjj_kfit", "baseLowMET/Iso/KinFit/mjj_kfit");
  outFile<<"\\hline "<<endl;

   getABCDNumbers(outFile, wjet, " W + jets",           "base/Iso/KinFit/mjj_kfit", "base/NonIso/KinFit/mjj_kfit", "baseLowMET/NonIso/KinFit/mjj_kfit", "baseLowMET/Iso/KinFit/mjj_kfit");
  outFile<<"\\hline "<<endl;

   getABCDNumbers(outFile, zjet, "$Z/\\gamma$ + jets",  "base/Iso/KinFit/mjj_kfit", "base/NonIso/KinFit/mjj_kfit", "baseLowMET/NonIso/KinFit/mjj_kfit", "baseLowMET/Iso/KinFit/mjj_kfit");
  outFile<<"\\hline "<<endl;
 
  getABCDNumbers(outFile, diboson, "VV",                "base/Iso/KinFit/mjj_kfit", "base/NonIso/KinFit/mjj_kfit", "baseLowMET/NonIso/KinFit/mjj_kfit", "baseLowMET/Iso/KinFit/mjj_kfit");
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  getABCDNumbers(outFile, allMC, "Bkg",                 "base/Iso/KinFit/mjj_kfit", "base/NonIso/KinFit/mjj_kfit", "baseLowMET/NonIso/KinFit/mjj_kfit", "baseLowMET/Iso/KinFit/mjj_kfit");
  outFile<<"\\hline "<<endl;

   getABCDNumbers(outFile, data, "Data",                "base/Iso/KinFit/mjj_kfit", "base/NonIso/KinFit/mjj_kfit", "baseLowMET/NonIso/KinFit/mjj_kfit", "baseLowMET/Iso/KinFit/mjj_kfit");
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\end{tabular}"<<endl; 
  outFile<<"\\end{center}"<<endl; 
  outFile<<"\\end{table}"<<endl; 
  outFile<<"\\end{document}"<<endl;  
  outFile.close(); 
} 
  
