#include <iostream>
#include <fstream>
#include <iomanip>

void makeCutFlowTable(bool NOTANPAS = false)
{
  TString inFile("$PWD/");

  cout << "inFile  " << inFile << endl;

  TFile *ttbar = new TFile(inFile+"all_TTJetsP.root");
  TH1F* h_ttbar = ((TH1F*)ttbar->Get("base/cutflow") )->Clone("h_ttbar");
  //TH1F* h_ttbar_toppt =((TH1F*)ttbar->Get("base/AvTopPtWeight") )->Clone("h_ttbar_toppt");

  double lumi = 19.703164;

  double wh_scale = 0.32; // 2x(1-x) assuming x = 0.2  

  //h_ttbar->Scale(lumi); 
  //h_ttbar->Scale(1.0/h_ttbar_toppt->GetBinContent(2));


  ofstream outFile;
  outFile.open("electron_cutflow.tex");
  outFile<< fixed << showpoint <<setprecision(1);
  if(NOTANPAS){
    outFile<<"\\documentclass[landscape,letterpaper]{article}"<<endl; 
    outFile<<"\\pagestyle{empty}"<<endl; 
    outFile<<"\\usepackage{epsfig}"<<endl; 
    outFile<<"\\usepackage{amsmath}"<<endl; 
    outFile<<"\\usepackage{array}"<<endl; 
    outFile<<"\\usepackage{multirow}"<<endl; 
    outFile<<"\\usepackage[cm]{fullpage}"<<endl; 
    outFile<<"\\textheight = 8.in"<<endl; 
    outFile<<"\\textwidth 7.0in"<<endl; 
    outFile<<"\\setlength{\\topmargin}{-.5in}"<<endl;
    outFile<<"\\setlength{\\textheight}{9in}"<<endl;
    outFile<<"\\setlength{\\oddsidemargin}{-.875in}"<<endl;
    
    outFile<<"\\begin{document}"<<endl; 
    outFile<<"Electron $p_T$ $\\ge$ 25 GeV,  $|\\eta|$ $<$ 2.5 ,$~~$ Jet  $p_T$ $\\ge$ 30 GeV, $|\\eta|$ $<$ 2.5 and 4 jets"<<" \\\\ "<<endl;
    outFile<<"$\\not\\!\\!E_T \\ge 20GeV $"<<" \\\\ "<< endl;
    outFile<<" With pt $\\ge$ 25 GeV in Kinematic fit and also applied MET px, py correction" << " \\\\ "<< endl;
    outFile<<"2 btag jet(2 Medium )  working point "<<" \\\\ "<<endl;
    outFile<<" Luminosity for full 2012 is 19.7 inv fb"<<" \\\\ "<<endl;
  }
  outFile<<"\\begin{center}"<<endl; 
  outFile<<"%\\begin{LARGE}"<<endl; 
  outFile<<"\\begin{tabular}{ | c| c| c| c| c|}"<<endl; 
  outFile<<"\\multicolumn{5}{c}{ } \\\\"<<endl; 
  outFile<<"\\hline "<<endl;
  outFile<<"\\multicolumn{1}{| c|}{ } & \\multicolumn{1}{ c|}{ $N_{electron}=1$ } & \\multicolumn{1}{ c|}{ $N_{jets}\\ge 4$ } & \\multicolumn{1}{ c|}{ $\\not\\!\\!E_T \\ge 20GeV$ }  &\\multicolumn{1}{ c |}{ $\\ge$ 2btag } \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;

  outFile<<"SM $t\\bar{t}$"<<" & "<<h_ttbar->GetBinContent(2)<<" & "<<h_ttbar->GetBinContent(3)<<" & "<<h_ttbar->GetBinContent(4)<<" & "<<h_ttbar->GetBinContent(5)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  
  outFile<<"\\end{tabular}"<<endl;
  outFile<<"%\\end{LARGE}"<<endl; 
  outFile<<"\\end{center}"<<endl; 
  if(NOTANPAS){
    outFile<<"\\end{document}"<<endl;
  }
  outFile.close();

}


void makeSummaryTableWithSys() 
{ 
  TString inFile("/afs/cern.ch/work/g/gkole/chargedHiggs/8TeV/JERStusy/uncertainty/ExampleAnalysis/8TeV/2M/kinfit11/v3/test5/");
  //TString inFile("/afs/cern.ch/work/g/gkole/chargedHiggs/8TeV/JERStusy/uncertainty/ExampleAnalysis/8TeV/2M/kinfit11/v3/test5/");
 
  TFile *wh = new TFile(inFile+"all_Hplus.root"); 
  TH1F* h_wh_base = ((TH1F*)wh->Get("base/cutflow") )->Clone("h_wh_base"); 
  TH1F* h_wh_JESPlus = ((TH1F*)wh->Get("JESPlus/cutflow") )->Clone("h_wh_JESPlus");
  TH1F* h_wh_JESMinus = ((TH1F*)wh->Get("JESMinus/cutflow") )->Clone("h_wh_JESMinus");
  TH1F* h_wh_JERPlus = ((TH1F*)wh->Get("JERPlus/cutflow") )->Clone("h_wh_JERPlus");
  TH1F* h_wh_JERMinus = ((TH1F*)wh->Get("JERMinus/cutflow") )->Clone("h_wh_JERMinus");
  TH1F* h_wh_METUCPlus = ((TH1F*)wh->Get("METUCPlus/cutflow") )->Clone("h_wh_METUCPlus");
  TH1F* h_wh_METUCMinus = ((TH1F*)wh->Get("METUCMinus/cutflow") )->Clone("h_wh_METUCMinus");
  TH1F* h_wh_bTagPlus = ((TH1F*)wh->Get("bTagPlus/cutflow") )->Clone("h_wh_bTagPlus");
  TH1F* h_wh_bTagMinus = ((TH1F*)wh->Get("bTagMinus/cutflow") )->Clone("h_wh_bTagMinus");
 
  TFile *ttbar = new TFile(inFile+"all_TTJetsP.root"); 
  TH1F* h_ttbar_base = ((TH1F*)ttbar->Get("base/cutflow") )->Clone("h_ttbar_base");  
  TH1F* h_ttbar_JESPlus = ((TH1F*)ttbar->Get("JESPlus/cutflow") )->Clone("h_ttbar_JESPlus"); 
  TH1F* h_ttbar_JESMinus = ((TH1F*)ttbar->Get("JESMinus/cutflow") )->Clone("h_ttbar_JESMinus"); 
  TH1F* h_ttbar_JERPlus = ((TH1F*)ttbar->Get("JERPlus/cutflow") )->Clone("h_ttbar_JERPlus"); 
  TH1F* h_ttbar_JERMinus = ((TH1F*)ttbar->Get("JERMinus/cutflow") )->Clone("h_ttbar_JERMinus"); 
  TH1F* h_ttbar_METUCPlus = ((TH1F*)ttbar->Get("METUCPlus/cutflow") )->Clone("h_ttbar_METUCPlus"); 
  TH1F* h_ttbar_METUCMinus = ((TH1F*)ttbar->Get("METUCMinus/cutflow") )->Clone("h_ttbar_METUCMinus"); 
  TH1F* h_ttbar_bTagPlus = ((TH1F*)ttbar->Get("bTagPlus/cutflow") )->Clone("h_ttbar_bTagPlus"); 
  TH1F* h_ttbar_bTagMinus = ((TH1F*)ttbar->Get("bTagMinus/cutflow") )->Clone("h_ttbar_bTagMinus"); 

   
  TFile *wjet = new TFile(inFile+"all_WJets.root"); 
  TH1F* h_wjet_base = ((TH1F*)wjet->Get("base/cutflow") )->Clone("h_wjet_base");  
  TH1F* h_wjet_JESPlus = ((TH1F*)wjet->Get("JESPlus/cutflow") )->Clone("h_wjet_JESPlus"); 
  TH1F* h_wjet_JESMinus = ((TH1F*)wjet->Get("JESMinus/cutflow") )->Clone("h_wjet_JESMinus"); 
  TH1F* h_wjet_JERPlus = ((TH1F*)wjet->Get("JERPlus/cutflow") )->Clone("h_wjet_JERPlus"); 
  TH1F* h_wjet_JERMinus = ((TH1F*)wjet->Get("JERMinus/cutflow") )->Clone("h_wjet_JERMinus"); 
  TH1F* h_wjet_METUCPlus = ((TH1F*)wjet->Get("METUCPlus/cutflow") )->Clone("h_wjet_METUCPlus"); 
  TH1F* h_wjet_METUCMinus = ((TH1F*)wjet->Get("METUCMinus/cutflow") )->Clone("h_wjet_METUCMinus"); 
  TH1F* h_wjet_bTagPlus = ((TH1F*)wjet->Get("bTagPlus/cutflow") )->Clone("h_wjet_bTagPlus"); 
  TH1F* h_wjet_bTagMinus = ((TH1F*)wjet->Get("bTagMinus/cutflow") )->Clone("h_wjet_bTagMinus"); 

  TFile *zjet = new TFile(inFile+"all_DY.root");
  TH1F* h_zjet_base = ((TH1F*)zjet->Get("base/cutflow") )->Clone("h_zjet_base");
  TH1F* h_zjet_JESPlus = ((TH1F*)zjet->Get("JESPlus/cutflow") )->Clone("h_zjet_JESPlus");
  TH1F* h_zjet_JESMinus = ((TH1F*)zjet->Get("JESMinus/cutflow") )->Clone("h_zjet_JESMinus");
  TH1F* h_zjet_JERPlus = ((TH1F*)zjet->Get("JERPlus/cutflow") )->Clone("h_zjet_JERPlus");
  TH1F* h_zjet_JERMinus = ((TH1F*)zjet->Get("JERMinus/cutflow") )->Clone("h_zjet_JERMinus");
  TH1F* h_zjet_METUCPlus = ((TH1F*)zjet->Get("METUCPlus/cutflow") )->Clone("h_zjet_METUCPlus");
  TH1F* h_zjet_METUCMinus = ((TH1F*)zjet->Get("METUCMinus/cutflow") )->Clone("h_zjet_METUCMinus");
  TH1F* h_zjet_bTagPlus = ((TH1F*)zjet->Get("bTagPlus/cutflow") )->Clone("h_zjet_bTagPlus");
  TH1F* h_zjet_bTagMinus = ((TH1F*)zjet->Get("bTagMinus/cutflow") )->Clone("h_zjet_bTagMinus");

  TFile *qcd = new TFile(inFile+"all_QCD.root");
  TH1F* h_qcd_base = ((TH1F*)qcd->Get("base/cutflow") )->Clone("h_qcd_base");
  TH1F* h_qcd_JESPlus = ((TH1F*)qcd->Get("JESPlus/cutflow") )->Clone("h_qcd_JESPlus");
  TH1F* h_qcd_JESMinus = ((TH1F*)qcd->Get("JESMinus/cutflow") )->Clone("h_qcd_JESMinus");
  TH1F* h_qcd_JERPlus = ((TH1F*)qcd->Get("JERPlus/cutflow") )->Clone("h_qcd_JERPlus");
  TH1F* h_qcd_JERMinus = ((TH1F*)qcd->Get("JERMinus/cutflow") )->Clone("h_qcd_JERMinus");
  TH1F* h_qcd_METUCPlus = ((TH1F*)qcd->Get("METUCPlus/cutflow") )->Clone("h_qcd_METUCPlus");
  TH1F* h_qcd_METUCMinus = ((TH1F*)qcd->Get("METUCMinus/cutflow") )->Clone("h_qcd_METUCMinus");
  TH1F* h_qcd_bTagPlus = ((TH1F*)qcd->Get("bTagPlus/cutflow") )->Clone("h_qcd_bTagPlus");
  TH1F* h_qcd_bTagMinus = ((TH1F*)qcd->Get("bTagMinus/cutflow") )->Clone("h_qcd_bTagMinus");

  TFile *stop = new TFile(inFile+"all_ST.root");
  TH1F* h_stop_base = ((TH1F*)stop->Get("base/cutflow") )->Clone("h_stop_base");
  TH1F* h_stop_JESPlus = ((TH1F*)stop->Get("JESPlus/cutflow") )->Clone("h_stop_JESPlus");
  TH1F* h_stop_JESMinus = ((TH1F*)stop->Get("JESMinus/cutflow") )->Clone("h_stop_JESMinus");
  TH1F* h_stop_JERPlus = ((TH1F*)stop->Get("JERPlus/cutflow") )->Clone("h_stop_JERPlus");
  TH1F* h_stop_JERMinus = ((TH1F*)stop->Get("JERMinus/cutflow") )->Clone("h_stop_JERMinus");
  TH1F* h_stop_METUCPlus = ((TH1F*)stop->Get("METUCPlus/cutflow") )->Clone("h_stop_METUCPlus");
  TH1F* h_stop_METUCMinus = ((TH1F*)stop->Get("METUCMinus/cutflow") )->Clone("h_stop_METUCMinus");
  TH1F* h_stop_bTagPlus = ((TH1F*)stop->Get("bTagPlus/cutflow") )->Clone("h_stop_bTagPlus");
  TH1F* h_stop_bTagMinus = ((TH1F*)stop->Get("bTagMinus/cutflow") )->Clone("h_stop_bTagMinus");

  TFile *diboson = new TFile(inFile+"all_VV.root");
  TH1F* h_diboson_base = ((TH1F*)diboson->Get("base/cutflow") )->Clone("h_diboson_base");
  TH1F* h_diboson_JESPlus = ((TH1F*)diboson->Get("JESPlus/cutflow") )->Clone("h_diboson_JESPlus");
  TH1F* h_diboson_JESMinus = ((TH1F*)diboson->Get("JESMinus/cutflow") )->Clone("h_diboson_JESMinus");
  TH1F* h_diboson_JERPlus = ((TH1F*)diboson->Get("JERPlus/cutflow") )->Clone("h_diboson_JERPlus");
  TH1F* h_diboson_JERMinus = ((TH1F*)diboson->Get("JERMinus/cutflow") )->Clone("h_diboson_JERMinus");
  TH1F* h_diboson_METUCPlus = ((TH1F*)diboson->Get("METUCPlus/cutflow") )->Clone("h_diboson_METUCPlus");
  TH1F* h_diboson_METUCMinus = ((TH1F*)diboson->Get("METUCMinus/cutflow") )->Clone("h_diboson_METUCMinus");
  TH1F* h_diboson_bTagPlus = ((TH1F*)diboson->Get("bTagPlus/cutflow") )->Clone("h_diboson_bTagPlus");
  TH1F* h_diboson_bTagMinus = ((TH1F*)diboson->Get("bTagMinus/cutflow") )->Clone("h_diboson_bTagMinus");


  h_wh_base->Scale(1.0/0.98637);
  h_wh_JESPlus->Scale(1.0/0.98637);
  h_wh_JESMinus->Scale(1.0/0.98637);
  h_wh_JERPlus->Scale(1.0/0.98637);
  h_wh_JERMinus->Scale(1.0/0.98637);
  h_wh_METUCPlus->Scale(1.0/0.98637);
  h_wh_METUCMinus->Scale(1.0/0.98637);
  h_wh_bTagPlus->Scale(1.0/0.98637);
  h_wh_bTagMinus->Scale(1.0/0.98637);

  h_ttbar_base->Scale(1.0/0.989606);
  h_ttbar_JESPlus->Scale(1.0/0.989606);
  h_ttbar_JESMinus->Scale(1.0/0.989606);
  h_ttbar_JERPlus->Scale(1.0/0.989606);
  h_ttbar_JERMinus->Scale(1.0/0.989606);
  h_ttbar_METUCPlus->Scale(1.0/0.989606);
  h_ttbar_METUCMinus->Scale(1.0/0.989606);
  h_ttbar_bTagPlus->Scale(1.0/0.989606);
  h_ttbar_bTagMinus->Scale(1.0/0.989606);


  TH1F* h_TotalBkg = h_wjet_base->Clone("h_TotalBkg");
  h_TotalBkg->Reset();
  h_TotalBkg->Add(h_ttbar_base);
  h_TotalBkg->Add(h_wjet_base);
  h_TotalBkg->Add(h_zjet_base);
  h_TotalBkg->Add(h_qcd_base);
  h_TotalBkg->Add(h_stop_base);
  h_TotalBkg->Add(h_diboson_base);


  TFile *data = new TFile(inFile+"data_selection_.root");
  TH1F* h_data = ((TH1F*)data->Get("base/cutflow") )->Clone("h_data");


  ofstream outFile; 
  outFile.open("summaryWithSyst.tex"); 
   
  outFile<<"\\documentclass[landscape,letterpaper]{article}"<<endl;  
  outFile<<"\\pagestyle{empty}"<<endl;  
  outFile<<"\\usepackage{epsfig}"<<endl;  
  outFile<<"\\usepackage{amsmath}"<<endl;  
  outFile<<"\\usepackage{array}"<<endl;  
  outFile<<"\\usepackage{multirow}"<<endl;  
  outFile<<"\\usepackage[cm]{fullpage}"<<endl;  
  outFile<<"\\textheight = 8.in"<<endl;  
  outFile<<"\\textwidth 7.0in"<<endl;  
  outFile<<"\\begin{document}"<<endl;  
  outFile<<"\\begin{center}"<<endl;  
  outFile<<"\\begin{LARGE}"<<endl;  
  outFile<<"\\begin{tabular}{ | c| c| }"<<endl;  
  outFile<<"\\multicolumn{2}{c}{ } \\\\"<<endl;  
  outFile<<"\\hline "<<endl; 
  outFile<<"\\multicolumn{1}{|c|}{Source} & \\multicolumn{1}{|c|}{N$_{\rm events}$ $\\pm$ MC stat $\\pm$ JES/MET scale $\\pm$ bTag } \\\\"<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<"\\hline "<<endl;
  outFile<<"HW, $M_{H}=120~GeV/c^{2}$"<<" & "<<h_wh_base->GetBinContent(7)<<" $\\pm$ "<<h_wh_base->GetBinError(7)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_wh_JESPlus->GetBinContent(7) - h_wh_base->GetBinContent(7)), fabs(h_wh_base->GetBinContent(7) - h_wh_JESMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_wh_JERPlus->GetBinContent(7) - h_wh_base->GetBinContent(7)), fabs(h_wh_base->GetBinContent(7) - h_wh_JERMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_wh_METUCPlus->GetBinContent(7) - h_wh_base->GetBinContent(7)), fabs(h_wh_base->GetBinContent(7) - h_wh_METUCMinus->GetBinContent(7))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_wh_bTagPlus->GetBinContent(7) - h_wh_base->GetBinContent(7)), fabs(h_wh_base->GetBinContent(7) - h_wh_bTagMinus->GetBinContent(7))) <<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;  
  outFile<<"SM $t\\bar{t}$"<<" & "<<h_ttbar_base->GetBinContent(7)<<" $\\pm$ "<<h_ttbar_base->GetBinError(7)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_ttbar_JESPlus->GetBinContent(7) - h_ttbar_base->GetBinContent(7)), fabs(h_ttbar_base->GetBinContent(7) - h_ttbar_JESMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_ttbar_JERPlus->GetBinContent(7) - h_ttbar_base->GetBinContent(7)), fabs(h_ttbar_base->GetBinContent(7) - h_ttbar_JERMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_ttbar_METUCPlus->GetBinContent(7) - h_ttbar_base->GetBinContent(7)), fabs(h_ttbar_base->GetBinContent(7) - h_ttbar_METUCMinus->GetBinContent(7))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_ttbar_bTagPlus->GetBinContent(7) - h_ttbar_base->GetBinContent(7)), fabs(h_ttbar_base->GetBinContent(7) - h_ttbar_bTagMinus->GetBinContent(7))) <<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;

  outFile<<"W+Jets"<<" & "<<h_wjet_base->GetBinContent(7)<<" $\\pm$ "<<h_wjet_base->GetBinError(7)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_wjet_JESPlus->GetBinContent(7) - h_wjet_base->GetBinContent(7)), fabs(h_wjet_base->GetBinContent(7) - h_wjet_JESMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_wjet_JERPlus->GetBinContent(7) - h_wjet_base->GetBinContent(7)), fabs(h_wjet_base->GetBinContent(7) - h_wjet_JERMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_wjet_METUCPlus->GetBinContent(7) - h_wjet_base->GetBinContent(7)), fabs(h_wjet_base->GetBinContent(7) - h_wjet_METUCMinus->GetBinContent(7))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_wjet_bTagPlus->GetBinContent(7) - h_wjet_base->GetBinContent(7)), fabs(h_wjet_base->GetBinContent(7) - h_wjet_bTagMinus->GetBinContent(7))) <<" \\\\ "<<endl;  
  outFile<<"\\hline "<<endl;  
  
  outFile<<"Z+Jets"<<" & "<<h_zjet_base->GetBinContent(7)<<" $\\pm$ "<<h_zjet_base->GetBinError(7)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_zjet_JESPlus->GetBinContent(7) - h_zjet_base->GetBinContent(7)), fabs(h_zjet_base->GetBinContent(7) - h_zjet_JESMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_zjet_JERPlus->GetBinContent(7) - h_zjet_base->GetBinContent(7)), fabs(h_zjet_base->GetBinContent(7) - h_zjet_JERMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_zjet_METUCPlus->GetBinContent(7) - h_zjet_base->GetBinContent(7)), fabs(h_zjet_base->GetBinContent(7) - h_zjet_METUCMinus->GetBinContent(7))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_zjet_bTagPlus->GetBinContent(7) - h_zjet_base->GetBinContent(7)), fabs(h_zjet_base->GetBinContent(7) - h_zjet_bTagMinus->GetBinContent(7))) <<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;

  outFile<<"QCD"<<" & "<<h_qcd_base->GetBinContent(7)<<" $\\pm$ "<<h_qcd_base->GetBinError(7)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_qcd_JESPlus->GetBinContent(7) - h_qcd_base->GetBinContent(7)), fabs(h_qcd_base->GetBinContent(7) - h_qcd_JESMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_qcd_JERPlus->GetBinContent(7) - h_qcd_base->GetBinContent(7)), fabs(h_qcd_base->GetBinContent(7) - h_qcd_JERMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_qcd_METUCPlus->GetBinContent(7) - h_qcd_base->GetBinContent(7)), fabs(h_qcd_base->GetBinContent(7) - h_qcd_METUCMinus->GetBinContent(7))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_qcd_bTagPlus->GetBinContent(7) - h_qcd_base->GetBinContent(7)), fabs(h_qcd_base->GetBinContent(7) - h_qcd_bTagMinus->GetBinContent(7))) <<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;

  outFile<<"SingleTop"<<" & "<<h_stop_base->GetBinContent(7)<<" $\\pm$ "<<h_stop_base->GetBinError(7)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_stop_JESPlus->GetBinContent(7) - h_stop_base->GetBinContent(7)), fabs(h_stop_base->GetBinContent(7) - h_stop_JESMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_stop_JERPlus->GetBinContent(7) - h_stop_base->GetBinContent(7)), fabs(h_stop_base->GetBinContent(7) - h_stop_JERMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_stop_METUCPlus->GetBinContent(7) - h_stop_base->GetBinContent(7)), fabs(h_stop_base->GetBinContent(7) - h_stop_METUCMinus->GetBinContent(7))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_stop_bTagPlus->GetBinContent(7) - h_stop_base->GetBinContent(7)), fabs(h_stop_base->GetBinContent(7) - h_stop_bTagMinus->GetBinContent(7))) <<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;

  outFile<<"Dibosons"<<" & "<<h_diboson_base->GetBinContent(7)<<" $\\pm$ "<<h_diboson_base->GetBinError(7)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_diboson_JESPlus->GetBinContent(7) - h_diboson_base->GetBinContent(7)), fabs(h_diboson_base->GetBinContent(7) - h_diboson_JESMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_diboson_JERPlus->GetBinContent(7) - h_diboson_base->GetBinContent(7)), fabs(h_diboson_base->GetBinContent(7) - h_diboson_JERMinus->GetBinContent(7))), 2) + pow(TMath::Max(fabs(h_diboson_METUCPlus->GetBinContent(7) - h_diboson_base->GetBinContent(7)), fabs(h_diboson_base->GetBinContent(7) - h_diboson_METUCMinus->GetBinContent(7))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_diboson_bTagPlus->GetBinContent(7) - h_diboson_base->GetBinContent(7)), fabs(h_diboson_base->GetBinContent(7) - h_diboson_bTagMinus->GetBinContent(7))) <<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"Total Bkg"<<" & "<<h_TotalBkg->GetBinContent(7)<<" $\\\pm$ "<<h_TotalBkg->GetBinError(7)<<" $\\\pm$ "<<" -- "<<" $\\\pm$ "<<" -- "<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;   
  outFile<<"\\hline "<<endl;  
  outFile<<" Data "<<" & "<<h_data->GetBinContent(7)<<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;    
  outFile<<"\\hline "<<endl;   
  outFile<<"\\end{tabular}"<<endl; 
  outFile<<"\\end{LARGE}"<<endl;  
  outFile<<"\\end{center}"<<endl;  
  outFile<<"\\end{document}"<<endl; 
 
  outFile.close(); 
} 

void makeSummaryTable(TString hist1="mjj_kfit_Id_svCat1", TString hist2="mjj_kfit_Id_svCat2", bool NOTANPAS = false)  
{  
  TString inFile("/afs/cern.ch/work/g/gkole/chargedHiggs/8TeV/JERStusy/uncertainty/ExampleAnalysis/8TeV/2M/kinfit11/v3/test5/");

  //declare fixed uncertainties on efficiency and scale factors
  // these numbers are the statistical uncertainties
  float effUncWH = 0.10;
  float effUncTTbar = 0.10;
  float effUncWJets = 0.05;
  float effUncZJets = 0.04;
  float effUncQcd = 0.50;
  float effUncSingletop = 0.06;
  float effUncDiboson = 0.10;
  
  TFile *wh = new TFile(inFile+"all_Hplus.root");  
  TH1F* h1_wh_toppt =((TH1F*)wh->Get("base/AvTopPtWeight") )->Clone("h1_wh_toppt");
  TH1F* h1_wh_base = ((TH1F*)wh->Get("base/"+hist1))->Clone("h1_wh_base");  
  TH1F* h1_wh_JESPlus = ((TH1F*)wh->Get("JESPlus/"+hist1) )->Clone("h1_wh_JESPlus"); 
  TH1F* h1_wh_JESMinus = ((TH1F*)wh->Get("JESMinus/"+hist1) )->Clone("h1_wh_JESMinus"); 
  TH1F* h1_wh_JERPlus = ((TH1F*)wh->Get("JERPlus/"+hist1) )->Clone("h1_wh_JERPlus"); 
  TH1F* h1_wh_JERMinus = ((TH1F*)wh->Get("JERMinus/"+hist1) )->Clone("h1_wh_JERMinus"); 
  TH1F* h1_wh_METUCPlus = ((TH1F*)wh->Get("METUCPlus/"+hist1) )->Clone("h1_wh_METUCPlus"); 
  TH1F* h1_wh_METUCMinus = ((TH1F*)wh->Get("METUCMinus/"+hist1) )->Clone("h1_wh_METUCMinus"); 
  TH1F* h1_wh_bTagPlus = ((TH1F*)wh->Get("bTagPlus/"+hist1) )->Clone("h1_wh_bTagPlus"); 
  TH1F* h1_wh_bTagMinus = ((TH1F*)wh->Get("bTagMinus/"+hist1) )->Clone("h1_wh_bTagMinus"); 
  
  TH1F* h2_wh_toppt =((TH1F*)wh->Get("base/AvTopPtWeight") )->Clone("h2_wh_toppt");
  TH1F* h2_wh_base = ((TH1F*)wh->Get("base/"+hist2))->Clone("h2_wh_base");  
  TH1F* h2_wh_JESPlus = ((TH1F*)wh->Get("JESPlus/"+hist2) )->Clone("h2_wh_JESPlus"); 
  TH1F* h2_wh_JESMinus = ((TH1F*)wh->Get("JESMinus/"+hist2) )->Clone("h2_wh_JESMinus"); 
  TH1F* h2_wh_JERPlus = ((TH1F*)wh->Get("JERPlus/"+hist2) )->Clone("h2_wh_JERPlus"); 
  TH1F* h2_wh_JERMinus = ((TH1F*)wh->Get("JERMinus/"+hist2) )->Clone("h2_wh_JERMinus"); 
  TH1F* h2_wh_METUCPlus = ((TH1F*)wh->Get("METUCPlus/"+hist2) )->Clone("h2_wh_METUCPlus"); 
  TH1F* h2_wh_METUCMinus = ((TH1F*)wh->Get("METUCMinus/"+hist2) )->Clone("h2_wh_METUCMinus"); 
  TH1F* h2_wh_bTagPlus = ((TH1F*)wh->Get("bTagPlus/"+hist2) )->Clone("h2_wh_bTagPlus"); 
  TH1F* h2_wh_bTagMinus = ((TH1F*)wh->Get("bTagMinus/"+hist2) )->Clone("h2_wh_bTagMinus"); 
  
  TFile *ttbar = new TFile(inFile+"all_TTJetsP.root");  
  TH1F* h1_ttbar_toppt =((TH1F*)ttbar->Get("base/AvTopPtWeight") )->Clone("h1_ttbar_toppt");
  TH1F* h1_ttbar_base = ((TH1F*)ttbar->Get("base/"+hist1) )->Clone("h1_ttbar_base");   
  TH1F* h1_ttbar_JESPlus = ((TH1F*)ttbar->Get("JESPlus/"+hist1) )->Clone("h1_ttbar_JESPlus");  
  TH1F* h1_ttbar_JESMinus = ((TH1F*)ttbar->Get("JESMinus/"+hist1) )->Clone("h1_ttbar_JESMinus");  
  TH1F* h1_ttbar_JERPlus = ((TH1F*)ttbar->Get("JERPlus/"+hist1) )->Clone("h1_ttbar_JERPlus");  
  TH1F* h1_ttbar_JERMinus = ((TH1F*)ttbar->Get("JERMinus/"+hist1) )->Clone("h1_ttbar_JERMinus");  
  TH1F* h1_ttbar_METUCPlus = ((TH1F*)ttbar->Get("METUCPlus/"+hist1) )->Clone("h1_ttbar_METUCPlus");  
  TH1F* h1_ttbar_METUCMinus = ((TH1F*)ttbar->Get("METUCMinus/"+hist1) )->Clone("h1_ttbar_METUCMinus");  
  TH1F* h1_ttbar_bTagPlus = ((TH1F*)ttbar->Get("bTagPlus/"+hist1) )->Clone("h1_ttbar_bTagPlus");  
  TH1F* h1_ttbar_bTagMinus = ((TH1F*)ttbar->Get("bTagMinus/"+hist1) )->Clone("h1_ttbar_bTagMinus");  
 
  TH1F* h2_ttbar_toppt =((TH1F*)ttbar->Get("base/AvTopPtWeight") )->Clone("h2_ttbar_toppt");
  TH1F* h2_ttbar_base = ((TH1F*)ttbar->Get("base/"+hist2) )->Clone("h2_ttbar_base");   
  TH1F* h2_ttbar_JESPlus = ((TH1F*)ttbar->Get("JESPlus/"+hist2) )->Clone("h2_ttbar_JESPlus");  
  TH1F* h2_ttbar_JESMinus = ((TH1F*)ttbar->Get("JESMinus/"+hist2) )->Clone("h2_ttbar_JESMinus");  
  TH1F* h2_ttbar_JERPlus = ((TH1F*)ttbar->Get("JERPlus/"+hist2) )->Clone("h2_ttbar_JERPlus");  
  TH1F* h2_ttbar_JERMinus = ((TH1F*)ttbar->Get("JERMinus/"+hist2) )->Clone("h2_ttbar_JERMinus");  
  TH1F* h2_ttbar_METUCPlus = ((TH1F*)ttbar->Get("METUCPlus/"+hist2) )->Clone("h2_ttbar_METUCPlus");  
  TH1F* h2_ttbar_METUCMinus = ((TH1F*)ttbar->Get("METUCMinus/"+hist2) )->Clone("h2_ttbar_METUCMinus");  
  TH1F* h2_ttbar_bTagPlus = ((TH1F*)ttbar->Get("bTagPlus/"+hist2) )->Clone("h2_ttbar_bTagPlus");  
  TH1F* h2_ttbar_bTagMinus = ((TH1F*)ttbar->Get("bTagMinus/"+hist2) )->Clone("h2_ttbar_bTagMinus");  
 
  TFile *wjet = new TFile(inFile+"all_WJets.root");  
  TH1F* h1_wjet_base = ((TH1F*)wjet->Get("base/"+hist1) )->Clone("h1_wjet_base");   
  TH1F* h1_wjet_JESPlus = ((TH1F*)wjet->Get("JESPlus/"+hist1) )->Clone("h1_wjet_JESPlus");  
  TH1F* h1_wjet_JESMinus = ((TH1F*)wjet->Get("JESMinus/"+hist1) )->Clone("h1_wjet_JESMinus");  
  TH1F* h1_wjet_JERPlus = ((TH1F*)wjet->Get("JERPlus/"+hist1) )->Clone("h1_wjet_JERPlus");  
  TH1F* h1_wjet_JERMinus = ((TH1F*)wjet->Get("JERMinus/"+hist1) )->Clone("h1_wjet_JERMinus");  
  TH1F* h1_wjet_METUCPlus = ((TH1F*)wjet->Get("METUCPlus/"+hist1) )->Clone("h1_wjet_METUCPlus");  
  TH1F* h1_wjet_METUCMinus = ((TH1F*)wjet->Get("METUCMinus/"+hist1) )->Clone("h1_wjet_METUCMinus");  
  TH1F* h1_wjet_bTagPlus = ((TH1F*)wjet->Get("bTagPlus/"+hist1) )->Clone("h1_wjet_bTagPlus");  
  TH1F* h1_wjet_bTagMinus = ((TH1F*)wjet->Get("bTagMinus/"+hist1) )->Clone("h1_wjet_bTagMinus");  
 
  TH1F* h2_wjet_base = ((TH1F*)wjet->Get("base/"+hist2) )->Clone("h2_wjet_base");   
  TH1F* h2_wjet_JESPlus = ((TH1F*)wjet->Get("JESPlus/"+hist2) )->Clone("h2_wjet_JESPlus");  
  TH1F* h2_wjet_JESMinus = ((TH1F*)wjet->Get("JESMinus/"+hist2) )->Clone("h2_wjet_JESMinus");  
  TH1F* h2_wjet_JERPlus = ((TH1F*)wjet->Get("JERPlus/"+hist2) )->Clone("h2_wjet_JERPlus");  
  TH1F* h2_wjet_JERMinus = ((TH1F*)wjet->Get("JERMinus/"+hist2) )->Clone("h2_wjet_JERMinus");  
  TH1F* h2_wjet_METUCPlus = ((TH1F*)wjet->Get("METUCPlus/"+hist2) )->Clone("h2_wjet_METUCPlus");  
  TH1F* h2_wjet_METUCMinus = ((TH1F*)wjet->Get("METUCMinus/"+hist2) )->Clone("h2_wjet_METUCMinus");  
  TH1F* h2_wjet_bTagPlus = ((TH1F*)wjet->Get("bTagPlus/"+hist2) )->Clone("h2_wjet_bTagPlus");  
  TH1F* h2_wjet_bTagMinus = ((TH1F*)wjet->Get("bTagMinus/"+hist2) )->Clone("h2_wjet_bTagMinus");  
 
  TFile *zjet = new TFile(inFile+"all_DY.root");
  TH1F* h1_zjet_base = ((TH1F*)zjet->Get("base/"+hist1) )->Clone("h1_zjet_base");
  TH1F* h1_zjet_JESPlus = ((TH1F*)zjet->Get("JESPlus/"+hist1) )->Clone("h1_zjet_JESPlus");
  TH1F* h1_zjet_JESMinus = ((TH1F*)zjet->Get("JESMinus/"+hist1) )->Clone("h1_zjet_JESMinus");
  TH1F* h1_zjet_JERPlus = ((TH1F*)zjet->Get("JERPlus/"+hist1) )->Clone("h1_zjet_JERPlus");
  TH1F* h1_zjet_JERMinus = ((TH1F*)zjet->Get("JERMinus/"+hist1) )->Clone("h1_zjet_JERMinus");
  TH1F* h1_zjet_METUCPlus = ((TH1F*)zjet->Get("METUCPlus/"+hist1) )->Clone("h1_zjet_METUCPlus");
  TH1F* h1_zjet_METUCMinus = ((TH1F*)zjet->Get("METUCMinus/"+hist1) )->Clone("h1_zjet_METUCMinus");
  TH1F* h1_zjet_bTagPlus = ((TH1F*)zjet->Get("bTagPlus/"+hist1) )->Clone("h1_zjet_bTagPlus");
  TH1F* h1_zjet_bTagMinus = ((TH1F*)zjet->Get("bTagMinus/"+hist1) )->Clone("h1_zjet_bTagMinus");

  TH1F* h2_zjet_base = ((TH1F*)zjet->Get("base/"+hist2) )->Clone("h2_zjet_base");
  TH1F* h2_zjet_JESPlus = ((TH1F*)zjet->Get("JESPlus/"+hist2) )->Clone("h2_zjet_JESPlus");
  TH1F* h2_zjet_JESMinus = ((TH1F*)zjet->Get("JESMinus/"+hist2) )->Clone("h2_zjet_JESMinus");
  TH1F* h2_zjet_JERPlus = ((TH1F*)zjet->Get("JERPlus/"+hist2) )->Clone("h2_zjet_JERPlus");
  TH1F* h2_zjet_JERMinus = ((TH1F*)zjet->Get("JERMinus/"+hist2) )->Clone("h2_zjet_JERMinus");
  TH1F* h2_zjet_METUCPlus = ((TH1F*)zjet->Get("METUCPlus/"+hist2) )->Clone("h2_zjet_METUCPlus");
  TH1F* h2_zjet_METUCMinus = ((TH1F*)zjet->Get("METUCMinus/"+hist2) )->Clone("h2_zjet_METUCMinus");
  TH1F* h2_zjet_bTagPlus = ((TH1F*)zjet->Get("bTagPlus/"+hist2) )->Clone("h2_zjet_bTagPlus");
  TH1F* h2_zjet_bTagMinus = ((TH1F*)zjet->Get("bTagMinus/"+hist2) )->Clone("h2_zjet_bTagMinus");

  TFile *qcd = new TFile(inFile+"all_QCD.root");
  TH1F* h1_qcd_base = ((TH1F*)qcd->Get("base/"+hist1) )->Clone("h1_qcd_base");
  TH1F* h1_qcd_JESPlus = ((TH1F*)qcd->Get("JESPlus/"+hist1) )->Clone("h1_qcd_JESPlus");
  TH1F* h1_qcd_JESMinus = ((TH1F*)qcd->Get("JESMinus/"+hist1) )->Clone("h1_qcd_JESMinus");
  TH1F* h1_qcd_JERPlus = ((TH1F*)qcd->Get("JERPlus/"+hist1) )->Clone("h1_qcd_JERPlus");
  TH1F* h1_qcd_JERMinus = ((TH1F*)qcd->Get("JERMinus/"+hist1) )->Clone("h1_qcd_JERMinus");
  TH1F* h1_qcd_METUCPlus = ((TH1F*)qcd->Get("METUCPlus/"+hist1) )->Clone("h1_qcd_METUCPlus");
  TH1F* h1_qcd_METUCMinus = ((TH1F*)qcd->Get("METUCMinus/"+hist1) )->Clone("h1_qcd_METUCMinus");
  TH1F* h1_qcd_bTagPlus = ((TH1F*)qcd->Get("bTagPlus/"+hist1) )->Clone("h1_qcd_bTagPlus");
  TH1F* h1_qcd_bTagMinus = ((TH1F*)qcd->Get("bTagMinus/"+hist1) )->Clone("h1_qcd_bTagMinus");

  TH1F* h2_qcd_base = ((TH1F*)qcd->Get("base/"+hist2) )->Clone("h2_qcd_base");
  TH1F* h2_qcd_JESPlus = ((TH1F*)qcd->Get("JESPlus/"+hist2) )->Clone("h2_qcd_JESPlus");
  TH1F* h2_qcd_JESMinus = ((TH1F*)qcd->Get("JESMinus/"+hist2) )->Clone("h2_qcd_JESMinus");
  TH1F* h2_qcd_JERPlus = ((TH1F*)qcd->Get("JERPlus/"+hist2) )->Clone("h2_qcd_JERPlus");
  TH1F* h2_qcd_JERMinus = ((TH1F*)qcd->Get("JERMinus/"+hist2) )->Clone("h2_qcd_JERMinus");
  TH1F* h2_qcd_METUCPlus = ((TH1F*)qcd->Get("METUCPlus/"+hist2) )->Clone("h2_qcd_METUCPlus");
  TH1F* h2_qcd_METUCMinus = ((TH1F*)qcd->Get("METUCMinus/"+hist2) )->Clone("h2_qcd_METUCMinus");
  TH1F* h2_qcd_bTagPlus = ((TH1F*)qcd->Get("bTagPlus/"+hist2) )->Clone("h2_qcd_bTagPlus");
  TH1F* h2_qcd_bTagMinus = ((TH1F*)qcd->Get("bTagMinus/"+hist2) )->Clone("h2_qcd_bTagMinus");

  TFile *stop = new TFile(inFile+"all_ST.root");
  TH1F* h1_stop_base = ((TH1F*)stop->Get("base/"+hist1) )->Clone("h1_stop_base");
  TH1F* h1_stop_JESPlus = ((TH1F*)stop->Get("JESPlus/"+hist1) )->Clone("h1_stop_JESPlus");
  TH1F* h1_stop_JESMinus = ((TH1F*)stop->Get("JESMinus/"+hist1) )->Clone("h1_stop_JESMinus");
  TH1F* h1_stop_JERPlus = ((TH1F*)stop->Get("JERPlus/"+hist1) )->Clone("h1_stop_JERPlus");
  TH1F* h1_stop_JERMinus = ((TH1F*)stop->Get("JERMinus/"+hist1) )->Clone("h1_stop_JERMinus");
  TH1F* h1_stop_METUCPlus = ((TH1F*)stop->Get("METUCPlus/"+hist1) )->Clone("h1_stop_METUCPlus");
  TH1F* h1_stop_METUCMinus = ((TH1F*)stop->Get("METUCMinus/"+hist1) )->Clone("h1_stop_METUCMinus");
  TH1F* h1_stop_bTagPlus = ((TH1F*)stop->Get("bTagPlus/"+hist1) )->Clone("h1_stop_bTagPlus");
  TH1F* h1_stop_bTagMinus = ((TH1F*)stop->Get("bTagMinus/"+hist1) )->Clone("h1_stop_bTagMinus");

  TH1F* h2_stop_base = ((TH1F*)stop->Get("base/"+hist2) )->Clone("h2_stop_base");
  TH1F* h2_stop_JESPlus = ((TH1F*)stop->Get("JESPlus/"+hist2) )->Clone("h2_stop_JESPlus");
  TH1F* h2_stop_JESMinus = ((TH1F*)stop->Get("JESMinus/"+hist2) )->Clone("h2_stop_JESMinus");
  TH1F* h2_stop_JERPlus = ((TH1F*)stop->Get("JERPlus/"+hist2) )->Clone("h2_stop_JERPlus");
  TH1F* h2_stop_JERMinus = ((TH1F*)stop->Get("JERMinus/"+hist2) )->Clone("h2_stop_JERMinus");
  TH1F* h2_stop_METUCPlus = ((TH1F*)stop->Get("METUCPlus/"+hist2) )->Clone("h2_stop_METUCPlus");
  TH1F* h2_stop_METUCMinus = ((TH1F*)stop->Get("METUCMinus/"+hist2) )->Clone("h2_stop_METUCMinus");
  TH1F* h2_stop_bTagPlus = ((TH1F*)stop->Get("bTagPlus/"+hist2) )->Clone("h2_stop_bTagPlus");
  TH1F* h2_stop_bTagMinus = ((TH1F*)stop->Get("bTagMinus/"+hist2) )->Clone("h2_stop_bTagMinus");

  TFile *diboson = new TFile(inFile+"all_VV.root");
  TH1F* h1_diboson_base = ((TH1F*)diboson->Get("base/"+hist1) )->Clone("h1_diboson_base");
  TH1F* h1_diboson_JESPlus = ((TH1F*)diboson->Get("JESPlus/"+hist1) )->Clone("h1_diboson_JESPlus");
  TH1F* h1_diboson_JESMinus = ((TH1F*)diboson->Get("JESMinus/"+hist1) )->Clone("h1_diboson_JESMinus");
  TH1F* h1_diboson_JERPlus = ((TH1F*)diboson->Get("JERPlus/"+hist1) )->Clone("h1_diboson_JERPlus");
  TH1F* h1_diboson_JERMinus = ((TH1F*)diboson->Get("JERMinus/"+hist1) )->Clone("h1_diboson_JERMinus");
  TH1F* h1_diboson_METUCPlus = ((TH1F*)diboson->Get("METUCPlus/"+hist1) )->Clone("h1_diboson_METUCPlus");
  TH1F* h1_diboson_METUCMinus = ((TH1F*)diboson->Get("METUCMinus/"+hist1) )->Clone("h1_diboson_METUCMinus");
  TH1F* h1_diboson_bTagPlus = ((TH1F*)diboson->Get("bTagPlus/"+hist1) )->Clone("h1_diboson_bTagPlus");
  TH1F* h1_diboson_bTagMinus = ((TH1F*)diboson->Get("bTagMinus/"+hist1) )->Clone("h1_diboson_bTagMinus");

  TH1F* h2_diboson_base = ((TH1F*)diboson->Get("base/"+hist2) )->Clone("h2_diboson_base");
  TH1F* h2_diboson_JESPlus = ((TH1F*)diboson->Get("JESPlus/"+hist2) )->Clone("h2_diboson_JESPlus");
  TH1F* h2_diboson_JESMinus = ((TH1F*)diboson->Get("JESMinus/"+hist2) )->Clone("h2_diboson_JESMinus");
  TH1F* h2_diboson_JERPlus = ((TH1F*)diboson->Get("JERPlus/"+hist2) )->Clone("h2_diboson_JERPlus");
  TH1F* h2_diboson_JERMinus = ((TH1F*)diboson->Get("JERMinus/"+hist2) )->Clone("h2_diboson_JERMinus");
  TH1F* h2_diboson_METUCPlus = ((TH1F*)diboson->Get("METUCPlus/"+hist2) )->Clone("h2_diboson_METUCPlus");
  TH1F* h2_diboson_METUCMinus = ((TH1F*)diboson->Get("METUCMinus/"+hist2) )->Clone("h2_diboson_METUCMinus");
  TH1F* h2_diboson_bTagPlus = ((TH1F*)diboson->Get("bTagPlus/"+hist2) )->Clone("h2_diboson_bTagPlus");
  TH1F* h2_diboson_bTagMinus = ((TH1F*)diboson->Get("bTagMinus/"+hist2) )->Clone("h2_diboson_bTagMinus");


  double wh_br = 0.32; // 2x(1-x)

  h1_wh_base->Scale(1.0/h1_wh_toppt->GetBinContent(2));
  h1_wh_JESPlus->Scale(1.0/h1_wh_toppt->GetBinContent(2));
  h1_wh_JESMinus->Scale(1.0/h1_wh_toppt->GetBinContent(2));
  h1_wh_JERPlus->Scale(1.0/h1_wh_toppt->GetBinContent(2));
  h1_wh_JERMinus->Scale(1.0/h1_wh_toppt->GetBinContent(2));
  h1_wh_METUCPlus->Scale(1.0/h1_wh_toppt->GetBinContent(2));
  h1_wh_METUCMinus->Scale(1.0/h1_wh_toppt->GetBinContent(2));
  h1_wh_bTagPlus->Scale(1.0/h1_wh_toppt->GetBinContent(2));
  h1_wh_bTagMinus->Scale(1.0/h1_wh_toppt->GetBinContent(2));

  h2_wh_base->Scale(1.0/h2_wh_toppt->GetBinContent(2));
  h2_wh_JESPlus->Scale(1.0/h2_wh_toppt->GetBinContent(2));
  h2_wh_JESMinus->Scale(1.0/h2_wh_toppt->GetBinContent(2));
  h2_wh_JERPlus->Scale(1.0/h2_wh_toppt->GetBinContent(2));
  h2_wh_JERMinus->Scale(1.0/h2_wh_toppt->GetBinContent(2));
  h2_wh_METUCPlus->Scale(1.0/h2_wh_toppt->GetBinContent(2));
  h2_wh_METUCMinus->Scale(1.0/h2_wh_toppt->GetBinContent(2));
  h2_wh_bTagPlus->Scale(1.0/h2_wh_toppt->GetBinContent(2));
  h2_wh_bTagMinus->Scale(1.0/h2_wh_toppt->GetBinContent(2));

  h1_ttbar_base->Scale(1.0/h1_ttbar_toppt->GetBinContent(2));
  h1_ttbar_JESPlus->Scale(1.0/h1_ttbar_toppt->GetBinContent(2));
  h1_ttbar_JESMinus->Scale(1.0/h1_ttbar_toppt->GetBinContent(2));
  h1_ttbar_JERPlus->Scale(1.0/h1_ttbar_toppt->GetBinContent(2));
  h1_ttbar_JERMinus->Scale(1.0/h1_ttbar_toppt->GetBinContent(2));
  h1_ttbar_METUCPlus->Scale(1.0/h1_ttbar_toppt->GetBinContent(2));
  h1_ttbar_METUCMinus->Scale(1.0/h1_ttbar_toppt->GetBinContent(2));
  h1_ttbar_bTagPlus->Scale(1.0/h1_ttbar_toppt->GetBinContent(2));
  h1_ttbar_bTagMinus->Scale(1.0/h1_ttbar_toppt->GetBinContent(2));
  
  h2_ttbar_base->Scale(1.0/h2_ttbar_toppt->GetBinContent(2));
  h2_ttbar_JESPlus->Scale(1.0/h2_ttbar_toppt->GetBinContent(2));
  h2_ttbar_JESMinus->Scale(1.0/h2_ttbar_toppt->GetBinContent(2));
  h2_ttbar_JERPlus->Scale(1.0/h2_ttbar_toppt->GetBinContent(2));
  h2_ttbar_JERMinus->Scale(1.0/h2_ttbar_toppt->GetBinContent(2));
  h2_ttbar_METUCPlus->Scale(1.0/h2_ttbar_toppt->GetBinContent(2));
  h2_ttbar_METUCMinus->Scale(1.0/h2_ttbar_toppt->GetBinContent(2));
  h2_ttbar_bTagPlus->Scale(1.0/h2_ttbar_toppt->GetBinContent(2));
  h2_ttbar_bTagMinus->Scale(1.0/h2_ttbar_toppt->GetBinContent(2));

  
  TH1F* h1_TotalBkg = h1_wjet_base->Clone("h1_TotalBkg");
  h1_TotalBkg->Reset();
  h1_TotalBkg->Add(h1_ttbar_base);
  h1_TotalBkg->Add(h1_wjet_base);
  h1_TotalBkg->Add(h1_zjet_base);
  h1_TotalBkg->Add(h1_qcd_base);
  h1_TotalBkg->Add(h1_stop_base);
  h1_TotalBkg->Add(h1_diboson_base);
  
  TH1F* h2_TotalBkg = h2_wjet_base->Clone("h2_TotalBkg");
  h2_TotalBkg->Reset();
  h2_TotalBkg->Add(h2_ttbar_base);
  h2_TotalBkg->Add(h2_wjet_base);
  h2_TotalBkg->Add(h2_zjet_base);
  h2_TotalBkg->Add(h2_qcd_base);
  h2_TotalBkg->Add(h2_stop_base);
  h2_TotalBkg->Add(h2_diboson_base);
  
  TFile *data = new TFile(inFile+"data_selection_.root");
  TH1F* h1_data = ((TH1F*)data->Get("base/"+hist1) )->Clone("h1_data");
  TH1F* h2_data = ((TH1F*)data->Get("base/"+hist2) )->Clone("h2_data");
 
 
  double sError_wh = 0.0;
  double sError_ttbar = 0.0;
  double sError_wjets = 0.0;
  double sError_zjets = 0.0;
  double sError_qcd = 0.0;
  double sError_stop = 0.0;
  double sError_diboson = 0.0;
  double sError_totbkg = 0.0;

  ofstream outFile, outFile_incl, outFile_unc_incl, outFile_unc1, outFile_unc2;  
  outFile.open(hist1+"_svCat2_cutflow.tex");
  outFile_incl.open(hist1+"_incl_cutflow.tex");
  outFile_unc_incl.open("unc_csbar_systematics_tab_incl.tex");
  outFile_unc1.open("unc_csbar_systematics_tab_1.tex");
  outFile_unc2.open("unc_csbar_systematics_tab_2.tex");

  outFile<< fixed << showpoint <<setprecision(1);
  outFile_incl<< fixed << showpoint <<setprecision(1);
  outFile_unc_incl<< fixed << showpoint <<setprecision(1);
  outFile_unc1<< fixed << showpoint <<setprecision(1);
  outFile_unc2<< fixed << showpoint <<setprecision(1);

  if(NOTANPAS){
    outFile<<"\\documentclass[landscape,letterpaper]{article}"<<endl;   
    outFile<<"\\pagestyle{empty}"<<endl;   
    outFile<<"\\usepackage{epsfig}"<<endl;   
    outFile<<"\\usepackage{amsmath}"<<endl;   
    outFile<<"\\usepackage{array}"<<endl;   
    outFile<<"\\usepackage{multirow}"<<endl;   
    outFile<<"\\usepackage[cm]{fullpage}"<<endl;   
    outFile<<"\\textheight = 8.in"<<endl;   
    outFile<<"\\textwidth 7.0in"<<endl;   
    outFile<<"\\begin{document}"<<endl;   
    outFile<<"\\begin{center}"<<endl;   
    outFile<<"\\begin{LARGE}"<<endl;
  }   
  if(NOTANPAS){
    outFile_incl<<"\\documentclass[landscape,letterpaper]{article}"<<endl;
    outFile_incl<<"\\pagestyle{empty}"<<endl;
    outFile_incl<<"\\usepackage{epsfig}"<<endl;
    outFile_incl<<"\\usepackage{amsmath}"<<endl;
    outFile_incl<<"\\usepackage{array}"<<endl;
    outFile_incl<<"\\usepackage{multirow}"<<endl;
    outFile_incl<<"\\usepackage[cm]{fullpage}"<<endl;
    outFile_incl<<"\\textheight = 8.in"<<endl;
    outFile_incl<<"\\textwidth 7.0in"<<endl;
    outFile_incl<<"\\begin{document}"<<endl;
    outFile_incl<<"\\begin{center}"<<endl;
    outFile_incl<<"\\begin{LARGE}"<<endl;
  }
  // unc for inclusive table
  if(NOTANPAS){
    outFile_unc_incl<<"\\documentclass[landscape,letterpaper]{article}"<<endl;
    outFile_unc_incl<<"\\pagestyle{empty}"<<endl;
    outFile_unc_incl<<"\\usepackage{epsfig}"<<endl;
    outFile_unc_incl<<"\\usepackage{amsmath}"<<endl;
    outFile_unc_incl<<"\\usepackage{array}"<<endl;
    outFile_unc_incl<<"\\usepackage{multirow}"<<endl;
    outFile_unc_incl<<"\\usepackage[cm]{fullpage}"<<endl;
    outFile_unc_incl<<"\\textheight = 8.in"<<endl;
    outFile_unc_incl<<"\\textwidth 7.0in"<<endl;
    outFile_unc_incl<<"\\begin{document}"<<endl;
    outFile_unc_incl<<"\\begin{center}"<<endl;
    outFile_unc_incl<<"\\begin{LARGE}"<<endl;
  }

  
// for unc table1
  if(NOTANPAS){
    outFile_unc1<<"\\documentclass[landscape,letterpaper]{article}"<<endl;
    outFile_unc1<<"\\pagestyle{empty}"<<endl;
    outFile_unc1<<"\\usepackage{epsfig}"<<endl;
    outFile_unc1<<"\\usepackage{amsmath}"<<endl;
    outFile_unc1<<"\\usepackage{array}"<<endl;
    outFile_unc1<<"\\usepackage{multirow}"<<endl;
    outFile_unc1<<"\\usepackage[cm]{fullpage}"<<endl;
    outFile_unc1<<"\\textheight = 8.in"<<endl;
    outFile_unc1<<"\\textwidth 7.0in"<<endl;
    outFile_unc1<<"\\begin{document}"<<endl;
    outFile_unc1<<"\\begin{center}"<<endl;
    outFile_unc1<<"\\begin{LARGE}"<<endl;
  }
// for unc table2
  if(NOTANPAS){
    outFile_unc2<<"\\documentclass[landscape,letterpaper]{article}"<<endl;
    outFile_unc2<<"\\pagestyle{empty}"<<endl;
    outFile_unc2<<"\\usepackage{epsfig}"<<endl;
    outFile_unc2<<"\\usepackage{amsmath}"<<endl;
    outFile_unc2<<"\\usepackage{array}"<<endl;
    outFile_unc2<<"\\usepackage{multirow}"<<endl;
    outFile_unc2<<"\\usepackage[cm]{fullpage}"<<endl;
    outFile_unc2<<"\\textheight = 8.in"<<endl;
    outFile_unc2<<"\\textwidth 7.0in"<<endl;
    outFile_unc2<<"\\begin{document}"<<endl;
    outFile_unc2<<"\\begin{center}"<<endl;
    outFile_unc2<<"\\begin{LARGE}"<<endl;
  }
  // for unc incl table
  outFile_unc_incl<<"\\begin{tabular}{|l|c|c|c|c|c|c|c|}"<<endl;
  outFile_unc_incl<<"\\hline "<<endl;
  outFile_unc_incl<<" & WH & $t\\bar{t} _{\\mu + {\\rm jets} }$ & W+jets & Z+jets & Single top & Dibosons \\\\"<<endl;
  outFile_unc_incl<<"\\hline "<<endl;
  outFile_unc_incl<<"\\hline "<<endl;

  // for unc table1
  outFile_unc1<<"\\begin{tabular}{|l|c|c|c|c|c|c|c|}"<<endl;
  outFile_unc1<<"\\hline "<<endl;
  outFile_unc1<<" & WH & $t\\bar{t} _{\\mu + {\\rm jets} }$ & W+jets & Z+jets & Single top & Dibosons \\\\"<<endl;
  outFile_unc1<<"\\hline "<<endl;
  outFile_unc1<<"\\hline "<<endl;

  // for unc table2
  outFile_unc2<<"\\begin{tabular}{|l|c|c|c|c|c|c|c|}"<<endl;
  outFile_unc2<<"\\hline "<<endl;
  outFile_unc2<<" & WH & $t\\bar{t} _{\\mu + {\\rm jets} }$ & W+jets & Z+jets & Single top & Dibosons \\\\"<<endl;
  outFile_unc2<<"\\hline "<<endl;
  outFile_unc2<<"\\hline "<<endl;

  // for events Yeild table inclusive
  outFile_incl<<"\\begin{tabular}{ | c| c| }"<<endl;
  outFile_incl<<"\\multicolumn{2}{c}{ } \\\\"<<endl;
  outFile_incl<<"\\hline "<<endl;
  outFile_incl<<"\\multicolumn{1}{|c|}{} & \\multicolumn{1}{|c|}{Events Details} \\\\" << endl;
  outFile_incl<<"\\hline "<<endl;
  outFile_incl<<"\\multicolumn{1}{|c|}{Source} & \\multicolumn{1}{|c|}{N$_{\\rm events}$ $\\pm$ MC stat $\\pm$ syst} \\\\"<<endl;
  outFile_incl<<"\\hline "<<endl;
  outFile_incl<<"\\hline "<<endl;

  outFile_incl<<"HW, $M_{H}=120~GeV/c^{2}$"<<" & "<<wh_br*h1_wh_base->IntegralAndError(1,40,sError_wh)<<" $\\pm$ "<<wh_br*sError_wh<<" $\\pm$ "<< wh_br*sqrt(pow(TMath::Max(fabs(h1_wh_JESPlus->Integral() - h1_wh_base->Integral()), fabs(h1_wh_base->Integral() - h1_wh_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wh_JERPlus->Integral() - h1_wh_base->Integral()), fabs(h1_wh_base->Integral() - h1_wh_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wh_METUCPlus->Integral() - h1_wh_base->Integral()), fabs(h1_wh_base->Integral() - h1_wh_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wh_bTagPlus->Integral() - h1_wh_base->Integral()), fabs(h1_wh_base->Integral() - h1_wh_bTagMinus->Integral())), 2) + pow(h1_wh_base->Integral()*effUncWH, 2)) << " \\\\ "<< endl;


  // for events yeild table category
  outFile<<"\\begin{tabular}{ | c| c| c| }"<<endl;   
  outFile<<"\\multicolumn{3}{c}{ } \\\\"<<endl;   
  outFile<<"\\hline "<<endl;
  outFile<<"\\multicolumn{1}{|c|}{} & \\multicolumn{1}{|c|}{Category I} & \\multicolumn{1}{|c|}{Category II} \\\\" << endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\multicolumn{1}{|c|}{Source} & \\multicolumn{1}{|c|}{N$_{\\rm events}$ $\\pm$ MC stat $\\pm$ syst} & \\multicolumn{1}{|c|}{N$_{\\rm events}$ $\\pm$ MC stat $\\pm$ syst} \\\\"<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<"\\hline "<<endl; 
  outFile<<"HW, $M_{H}=120~GeV/c^{2}$"<<" & "<<wh_br*h1_wh_base->IntegralAndError(1,40,sError_wh)<<" $\\pm$ "<<wh_br*sError_wh<<" $\\pm$ "<< wh_br*sqrt(pow(TMath::Max(fabs(h1_wh_JESPlus->Integral() - h1_wh_base->Integral()), fabs(h1_wh_base->Integral() - h1_wh_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wh_JERPlus->Integral() - h1_wh_base->Integral()), fabs(h1_wh_base->Integral() - h1_wh_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wh_METUCPlus->Integral() - h1_wh_base->Integral()), fabs(h1_wh_base->Integral() - h1_wh_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wh_bTagPlus->Integral() - h1_wh_base->Integral()), fabs(h1_wh_base->Integral() - h1_wh_bTagMinus->Integral())), 2) + pow(h1_wh_base->Integral()*effUncWH, 2)) <<" & "<< wh_br*h2_wh_base->IntegralAndError(1,40,sError_wh)<<" $\\pm$ "<<wh_br*sError_wh<<" $\\pm$ "<< wh_br*sqrt(pow(TMath::Max(fabs(h2_wh_JESPlus->Integral() - h2_wh_base->Integral()), fabs(h2_wh_base->Integral() - h2_wh_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_wh_JERPlus->Integral() - h2_wh_base->Integral()), fabs(h2_wh_base->Integral() - h2_wh_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_wh_METUCPlus->Integral() - h2_wh_base->Integral()), fabs(h2_wh_base->Integral() - h2_wh_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_wh_bTagPlus->Integral() - h2_wh_base->Integral()), fabs(h2_wh_base->Integral() - h2_wh_bTagMinus->Integral())), 2) + pow(h2_wh_base->Integral()*effUncWH, 2)) <<" \\\\ "<<endl; 
  
  double unc_WH_cat1  =  sqrt(pow(TMath::Max(fabs(h1_wh_JESPlus->Integral() - h1_wh_base->Integral()), fabs(h1_wh_base->Integral() - h1_wh_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wh_JERPlus->Integral() - h1_wh_base->Integral()), fabs(h1_wh_base->Integral() - h1_wh_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wh_METUCPlus->Integral() - h1_wh_base->Integral()), fabs(h1_wh_base->Integral() - h1_wh_METUCMinus->Integral())), 2))/h1_wh_base->Integral() ;
  
  double unc_WH_cat2 =  sqrt(pow(TMath::Max(fabs(h2_wh_JESPlus->Integral() - h2_wh_base->Integral()), fabs(h2_wh_base->Integral() - h2_wh_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_wh_JERPlus->Integral() - h2_wh_base->Integral()), fabs(h2_wh_base->Integral() - h2_wh_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_wh_METUCPlus->Integral() - h2_wh_base->Integral()), fabs(h2_wh_base->Integral() - h2_wh_METUCMinus->Integral())), 2))/h2_wh_base->Integral() ;
  
  float bTagUnc_wh_cat1 = (h1_wh_base->Integral() > 0) ? (TMath::Max(fabs(h1_wh_bTagPlus->Integral() - h1_wh_base->Integral()), fabs(h1_wh_base->Integral() - h1_wh_bTagMinus->Integral()))/h1_wh_base->Integral()) : 0.00;

  float bTagUnc_wh_cat2 = (h2_wh_base->Integral() > 0) ? (TMath::Max(fabs(h2_wh_bTagPlus->Integral() - h2_wh_base->Integral()), fabs(h2_wh_base->Integral() - h2_wh_bTagMinus->Integral()))/h2_wh_base->Integral()) : 0.00;
  Double_t sError = 0;
  Double_t norm = h1_wh_base->IntegralAndError(1, h1_wh_base->GetNbinsX(), sError);
  float statUnc_wh_cat1 = (norm > 0) ? (fabs(sError)/norm) : 0.00;
  Double_t sError = 0;
  Double_t norm = h2_wh_base->IntegralAndError(1, h2_wh_base->GetNbinsX(), sError);
  float statUnc_wh_cat2 = (norm > 0) ? (fabs(sError)/norm) : 0.00;
  
  outFile_incl<<"\\hline "<<endl;
  outFile_incl<<"SM $t\\bar{t}$"<<" & "<<h1_ttbar_base->IntegralAndError(1,40,sError_ttbar)<<" $\\pm$ "<<sError_ttbar<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h1_ttbar_JESPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_ttbar_JERPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_ttbar_METUCPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_ttbar_bTagPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_bTagMinus->Integral())), 2) + pow(h1_ttbar_base->Integral()*effUncTTbar, 2)) << " \\\\ "<<endl;
  outFile_incl<<"\\hline "<<endl;

  outFile<<"\\hline "<<endl;  
  outFile<<"SM $t\\bar{t}$"<<" & "<<h1_ttbar_base->IntegralAndError(1,40,sError_ttbar)<<" $\\pm$ "<<sError_ttbar<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h1_ttbar_JESPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_ttbar_JERPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_ttbar_METUCPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_ttbar_bTagPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_bTagMinus->Integral())), 2) + pow(h1_ttbar_base->Integral()*effUncTTbar, 2)) <<" & " << h2_ttbar_base->IntegralAndError(1,40,sError_ttbar)<<" $\\pm$ "<<sError_ttbar<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h2_ttbar_JESPlus->Integral() - h2_ttbar_base->Integral()), fabs(h2_ttbar_base->Integral() - h2_ttbar_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_ttbar_JERPlus->Integral() - h2_ttbar_base->Integral()), fabs(h2_ttbar_base->Integral() - h2_ttbar_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_ttbar_METUCPlus->Integral() - h2_ttbar_base->Integral()), fabs(h2_ttbar_base->Integral() - h2_ttbar_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_ttbar_bTagPlus->Integral() - h2_ttbar_base->Integral()), fabs(h2_ttbar_base->Integral() - h2_ttbar_bTagMinus->Integral())), 2) + pow(h2_ttbar_base->Integral()*effUncTTbar, 2))  <<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;  

  double unc_ttbar_cat1 =  sqrt(pow(TMath::Max(fabs(h1_ttbar_JESPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_ttbar_JERPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_ttbar_METUCPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_METUCMinus->Integral())), 2))/h1_ttbar_base->Integral() ;

  double unc_ttbar_cat2 =  sqrt(pow(TMath::Max(fabs(h2_ttbar_JESPlus->Integral() - h2_ttbar_base->Integral()), fabs(h2_ttbar_base->Integral() - h2_ttbar_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_ttbar_JERPlus->Integral() - h2_ttbar_base->Integral()), fabs(h2_ttbar_base->Integral() - h2_ttbar_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_ttbar_METUCPlus->Integral() - h2_ttbar_base->Integral()), fabs(h2_ttbar_base->Integral() - h2_ttbar_METUCMinus->Integral())), 2))/h2_ttbar_base->Integral() ;

  float bTagUnc_ttbar_cat1 = (h1_ttbar_base->Integral() > 0) ? (TMath::Max(fabs(h1_ttbar_bTagPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_bTagMinus->Integral()))/h1_ttbar_base->Integral()) : 0.00;
  float bTagUnc_ttbar_cat2 = (h2_ttbar_base->Integral() > 0) ? (TMath::Max(fabs(h2_ttbar_bTagPlus->Integral() - h2_ttbar_base->Integral()), fabs(h2_ttbar_base->Integral() - h2_ttbar_bTagMinus->Integral()))/h2_ttbar_base->Integral()) : 0.00;

  Double_t sError = 0;
  Double_t norm = h1_ttbar_base->IntegralAndError(1, h1_ttbar_base->GetNbinsX(), sError);
  float statUnc_ttbar_cat1 = (norm > 0) ? (fabs(sError)/norm) : 0.00;

  Double_t sError = 0;
  Double_t norm = h2_ttbar_base->IntegralAndError(1, h2_ttbar_base->GetNbinsX(), sError);
  float statUnc_ttbar_cat2 = (norm > 0) ? (fabs(sError)/norm) : 0.00;
  
  outFile_incl<<"W+Jets"<<" & "<<h1_wjet_base->IntegralAndError(1,40,sError_wjets)<<" $\\pm$ "<<sError_wjets<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h1_wjet_JESPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wjet_METUCPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wjet_bTagPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_bTagMinus->Integral())), 2) + pow(h1_wjet_base->Integral()*effUncWJets, 2)) <<" \\\\ "<<endl;
  outFile_incl<<"\\hline "<<endl;

  outFile<<"W+Jets"<<" & "<<h1_wjet_base->IntegralAndError(1,40,sError_wjets)<<" $\\pm$ "<<sError_wjets<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h1_wjet_JESPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wjet_METUCPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wjet_bTagPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_bTagMinus->Integral())), 2) + pow(h1_wjet_base->Integral()*effUncWJets, 2)) <<" & "<<h2_wjet_base->IntegralAndError(1,40,sError_wjets)<<" $\\pm$ "<<sError_wjets<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h2_wjet_JESPlus->Integral() - h2_wjet_base->Integral()), fabs(h2_wjet_base->Integral() - h2_wjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_wjet_METUCPlus->Integral() - h2_wjet_base->Integral()), fabs(h2_wjet_base->Integral() - h2_wjet_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_wjet_bTagPlus->Integral() - h2_wjet_base->Integral()), fabs(h2_wjet_base->Integral() - h2_wjet_bTagMinus->Integral())), 2) + pow(h2_wjet_base->Integral()*effUncWJets, 2))  <<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;

  double unc_wjets_cat1 =  sqrt(pow(TMath::Max(fabs(h1_wjet_JESPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wjet_METUCPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_METUCMinus->Integral())), 2))/h1_wjet_base->Integral() ;

  // debug1
  double debug1 = pow(TMath::Max(fabs(h1_wjet_JESPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_JESMinus->Integral())), 2);
  cout << "debug1=  " << debug1 << endl;
  double debug2 = pow(TMath::Max(fabs(h1_wjet_JERPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_JERMinus->Integral())), 2);
  cout << "debug2=  " << debug2 << endl;
  double debug3 = pow(TMath::Max(fabs(h1_wjet_METUCPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_METUCMinus->Integral())), 2);
  cout << "debug3=  " << debug3 << endl;
  cout << "base_int=  " << h1_wjet_base->Integral() << endl;
  cout << "wjet_JERPlus=  "<< h1_wjet_JERPlus->Integral() << endl;
  cout << "wjet_JERMinus=  " << h1_wjet_JERMinus << endl;


  double unc_wjets_cat2 = sqrt(pow(TMath::Max(fabs(h2_wjet_JESPlus->Integral() - h2_wjet_base->Integral()), fabs(h2_wjet_base->Integral() - h2_wjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_wjet_METUCPlus->Integral() - h2_wjet_base->Integral()), fabs(h2_wjet_base->Integral() - h2_wjet_METUCMinus->Integral())), 2))/h2_wjet_base->Integral() ;

  float bTagUnc_wjet_cat1 = (h1_wjet_base->Integral() > 0) ? (TMath::Max(fabs(h1_wjet_bTagPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_bTagMinus->Integral()))/h1_wjet_base->Integral()) : 0.00;
  float bTagUnc_wjet_cat2 = (h2_wjet_base->Integral() > 0) ? (TMath::Max(fabs(h2_wjet_bTagPlus->Integral() - h2_wjet_base->Integral()), fabs(h2_wjet_base->Integral() - h2_wjet_bTagMinus->Integral()))/h2_wjet_base->Integral()) : 0.00;

  Double_t sError = 0;
  Double_t norm = h1_wjet_base->IntegralAndError(1, h1_wjet_base->GetNbinsX(), sError);
  float statUnc_wjet_cat1 = (norm > 0) ? (fabs(sError)/norm) : 0.00;

  Double_t sError = 0;
  Double_t norm = h2_wjet_base->IntegralAndError(1, h2_wjet_base->GetNbinsX(), sError);
  float statUnc_wjet_cat2 = (norm > 0) ? (fabs(sError)/norm) : 0.00;

  outFile_incl<<"Z+Jets"<<" & "<<h1_zjet_base->IntegralAndError(1,40,sError_zjets)<<" $\\pm$ "<<sError_zjets<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h1_zjet_JESPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_zjet_METUCPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_zjet_bTagPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_bTagMinus->Integral())), 2) + pow(h1_zjet_base->Integral()*effUncZJets, 2)) <<" \\\\ "<<endl;
  outFile_incl<<"\\hline "<< endl;

  outFile<<"Z+Jets"<<" & "<<h1_zjet_base->IntegralAndError(1,40,sError_zjets)<<" $\\pm$ "<<sError_zjets<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h1_zjet_JESPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_zjet_METUCPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_zjet_bTagPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_bTagMinus->Integral())), 2) + pow(h1_zjet_base->Integral()*effUncZJets, 2)) <<" & "<<h2_zjet_base->IntegralAndError(1,40,sError_zjets)<<" $\\pm$ "<<sError_zjets<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h2_zjet_JESPlus->Integral() - h2_zjet_base->Integral()), fabs(h2_zjet_base->Integral() - h2_zjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_zjet_METUCPlus->Integral() - h2_zjet_base->Integral()), fabs(h2_zjet_base->Integral() - h2_zjet_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_zjet_bTagPlus->Integral() - h2_zjet_base->Integral()), fabs(h2_zjet_base->Integral() - h2_zjet_bTagMinus->Integral())), 2) + pow(h2_zjet_base->Integral()*effUncZJets, 2))  <<" \\\\ "<<endl;

  outFile<<"\\hline "<<endl;
 double unc_zjets_cat1 = sqrt(pow(TMath::Max(fabs(h1_zjet_JESPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_zjet_METUCPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_METUCMinus->Integral())), 2))/h1_zjet_base->Integral() ;

 double unc_zjets_cat2 = sqrt(pow(TMath::Max(fabs(h2_zjet_JESPlus->Integral() - h2_zjet_base->Integral()), fabs(h2_zjet_base->Integral() - h2_zjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_zjet_METUCPlus->Integral() - h2_zjet_base->Integral()), fabs(h2_zjet_base->Integral() - h2_zjet_METUCMinus->Integral())), 2))/h2_zjet_base->Integral() ;

 float bTagUnc_zjet_cat1 = (h1_zjet_base->Integral() > 0) ? (TMath::Max(fabs(h1_zjet_bTagPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_bTagMinus->Integral()))/h1_zjet_base->Integral()) : 0.00;
 float bTagUnc_zjet_cat2 = (h2_zjet_base->Integral() > 0) ? (TMath::Max(fabs(h2_zjet_bTagPlus->Integral() - h2_zjet_base->Integral()), fabs(h2_zjet_base->Integral() - h2_zjet_bTagMinus->Integral()))/h2_zjet_base->Integral()) : 0.00;

 Double_t sError = 0;
 Double_t norm = h1_zjet_base->IntegralAndError(1, h1_zjet_base->GetNbinsX(), sError);
 float statUnc_zjet_cat1 = (norm > 0) ? (fabs(sError)/norm) : 0.00;

 Double_t sError = 0;
 Double_t norm = h2_zjet_base->IntegralAndError(1, h2_zjet_base->GetNbinsX(), sError);
 float statUnc_zjet_cat2 = (norm > 0) ? (fabs(sError)/norm) : 0.00;

  // outFile<<"QCD"<<" & "<<h1_qcd_base->IntegralAndError(1,40,sError_qcd)<<" $\\pm$ "<<sError_qcd<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h1_qcd_JESPlus->Integral() - h1_qcd_base->Integral()), fabs(h1_qcd_base->Integral() - h1_qcd_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_qcd_JERPlus->Integral() - h1_qcd_base->Integral()), fabs(h1_qcd_base->Integral() - h1_qcd_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_qcd_METUCPlus->Integral() - h1_qcd_base->Integral()), fabs(h1_qcd_base->Integral() - h1_qcd_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_qcd_bTagPlus->Integral() - h1_qcd_base->Integral()), fabs(h1_qcd_base->Integral() - h1_qcd_bTagMinus->Integral())), 2) + pow(h1_qcd_base->Integral()*effUncQcd, 2)) <<" \\\\ "<<endl;
 outFile_incl <<"SingleTop"<<" & "<<h1_stop_base->IntegralAndError(1,40,sError_stop)<<" $\\pm$ "<<sError_stop<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h1_stop_JESPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_stop_METUCPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_stop_bTagPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_bTagMinus->Integral())), 2) + pow(h1_stop_base->Integral()*effUncSingletop, 2)) <<" \\\\ " << endl;
   outFile_incl<<"\\hline "<<endl;

 outFile<<"SingleTop"<<" & "<<h1_stop_base->IntegralAndError(1,40,sError_stop)<<" $\\pm$ "<<sError_stop<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h1_stop_JESPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_stop_METUCPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_stop_bTagPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_bTagMinus->Integral())), 2) + pow(h1_stop_base->Integral()*effUncSingletop, 2))<<" & "<<h2_stop_base->IntegralAndError(1,40,sError_stop)<<" $\\pm$ "<<sError_stop<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h2_stop_JESPlus->Integral() - h2_stop_base->Integral()), fabs(h2_stop_base->Integral() - h2_stop_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_stop_METUCPlus->Integral() - h2_stop_base->Integral()), fabs(h2_stop_base->Integral() - h2_stop_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_stop_bTagPlus->Integral() - h2_stop_base->Integral()), fabs(h2_stop_base->Integral() - h2_stop_bTagMinus->Integral())), 2) + pow(h2_stop_base->Integral()*effUncSingletop, 2)) <<" \\\\ "<<endl;

 outFile<<"\\hline "<<endl;
 double unc_singletop_cat1 =  sqrt(pow(TMath::Max(fabs(h1_stop_JESPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_stop_METUCPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_METUCMinus->Integral())), 2))/h1_stop_base->Integral() ;

 double unc_singletop_cat2 =  sqrt(pow(TMath::Max(fabs(h2_stop_JESPlus->Integral() - h2_stop_base->Integral()), fabs(h2_stop_base->Integral() - h2_stop_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_stop_METUCPlus->Integral() - h2_stop_base->Integral()), fabs(h2_stop_base->Integral() - h2_stop_METUCMinus->Integral())), 2))/h2_stop_base->Integral() ;

 float bTagUnc_stop_cat1 = (h1_stop_base->Integral() > 0) ? (TMath::Max(fabs(h1_stop_bTagPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_bTagMinus->Integral()))/h1_stop_base->Integral()) : 0.00;
 float bTagUnc_stop_cat2 = (h2_stop_base->Integral() > 0) ? (TMath::Max(fabs(h2_stop_bTagPlus->Integral() - h2_stop_base->Integral()), fabs(h2_stop_base->Integral() - h2_stop_bTagMinus->Integral()))/h2_stop_base->Integral()) : 0.00;

 Double_t sError = 0;
 Double_t norm = h1_stop_base->IntegralAndError(1, h1_stop_base->GetNbinsX(), sError);
 float statUnc_stop_cat1 = (norm > 0) ? (fabs(sError)/norm) : 0.00;

 Double_t sError = 0;
 Double_t norm = h2_stop_base->IntegralAndError(1, h2_stop_base->GetNbinsX(), sError);
 float statUnc_stop_cat2 = (norm > 0) ? (fabs(sError)/norm) : 0.00;
 
 outFile_incl <<"Dibosons"<<" & "<<h1_diboson_base->IntegralAndError(1,40,sError_diboson)<<" $\\pm$ "<<sError_diboson<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h1_diboson_JESPlus->Integral() - h1_diboson_base->Integral()), fabs(h1_diboson_base->Integral() - h1_diboson_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_diboson_METUCPlus->Integral() - h1_diboson_base->Integral()), fabs(h1_diboson_base->Integral() - h1_diboson_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_diboson_bTagPlus->Integral() - h1_diboson_base->Integral()), fabs(h1_diboson_base->Integral() - h1_diboson_bTagMinus->Integral())), 2) + pow(h1_diboson_base->Integral()*effUncDiboson, 2)) << " \\\\ "<< endl;

 outFile_incl<<"\\hline "<<endl;

 outFile<<"Dibosons"<<" & "<<h1_diboson_base->IntegralAndError(1,40,sError_diboson)<<" $\\pm$ "<<sError_diboson<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h1_diboson_JESPlus->Integral() - h1_diboson_base->Integral()), fabs(h1_diboson_base->Integral() - h1_diboson_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_diboson_METUCPlus->Integral() - h1_diboson_base->Integral()), fabs(h1_diboson_base->Integral() - h1_diboson_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_diboson_bTagPlus->Integral() - h1_diboson_base->Integral()), fabs(h1_diboson_base->Integral() - h1_diboson_bTagMinus->Integral())), 2) + pow(h1_diboson_base->Integral()*effUncDiboson, 2)) <<" & "<< h2_diboson_base->IntegralAndError(1,40,sError_diboson)<<" $\\pm$ "<<sError_diboson<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h2_diboson_JESPlus->Integral() - h2_diboson_base->Integral()), fabs(h2_diboson_base->Integral() - h2_diboson_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_diboson_METUCPlus->Integral() - h2_diboson_base->Integral()), fabs(h2_diboson_base->Integral() - h2_diboson_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_diboson_bTagPlus->Integral() - h2_diboson_base->Integral()), fabs(h2_diboson_base->Integral() - h2_diboson_bTagMinus->Integral())), 2) + pow(h2_diboson_base->Integral()*effUncDiboson, 2)) <<" \\\\ "<<endl;

  outFile<<"\\hline "<<endl;  

  double unc_diboson_cat1 =  sqrt(pow(TMath::Max(fabs(h1_diboson_JESPlus->Integral() - h1_diboson_base->Integral()), fabs(h1_diboson_base->Integral() - h1_diboson_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_diboson_METUCPlus->Integral() - h1_diboson_base->Integral()), fabs(h1_diboson_base->Integral() - h1_diboson_METUCMinus->Integral())), 2))/h1_diboson_base->Integral() ;

  double unc_diboson_cat2 = sqrt(pow(TMath::Max(fabs(h2_diboson_JESPlus->Integral() - h2_diboson_base->Integral()), fabs(h2_diboson_base->Integral() - h2_diboson_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_diboson_METUCPlus->Integral() - h2_diboson_base->Integral()), fabs(h2_diboson_base->Integral() - h2_diboson_METUCMinus->Integral())), 2))/h2_diboson_base->Integral() ;

  float bTagUnc_diboson_cat1 = (h1_diboson_base->Integral() > 0) ? (TMath::Max(fabs(h1_diboson_bTagPlus->Integral() - h1_diboson_base->Integral()), fabs(h1_diboson_base->Integral() - h1_diboson_bTagMinus->Integral()))/h1_diboson_base->Integral()) : 0.00;
  float bTagUnc_diboson_cat2 = (h2_diboson_base->Integral() > 0) ? (TMath::Max(fabs(h2_diboson_bTagPlus->Integral() - h2_diboson_base->Integral()), fabs(h2_diboson_base->Integral() - h2_diboson_bTagMinus->Integral()))/h2_diboson_base->Integral()) : 0.00;

  Double_t sError = 0;
  Double_t norm = h1_diboson_base->IntegralAndError(1, h1_diboson_base->GetNbinsX(), sError);
  float statUnc_diboson_cat1 = (norm > 0) ? (fabs(sError)/norm) : 0.00;

  Double_t sError = 0;
  Double_t norm = h2_diboson_base->IntegralAndError(1, h2_diboson_base->GetNbinsX(), sError);
  float statUnc_diboson_cat2 = (norm > 0) ? (fabs(sError)/norm) : 0.00;

  double tmp_ttbar_tot_unc_cat1 = pow(TMath::Max(fabs(h1_ttbar_JESPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_ttbar_JERPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_ttbar_METUCPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_ttbar_bTagPlus->Integral() - h1_ttbar_base->Integral()), fabs(h1_ttbar_base->Integral() - h1_ttbar_bTagMinus->Integral())), 2) + pow(h1_ttbar_base->Integral()*effUncTTbar, 2) ;
   
  double tmp_ttbar_tot_unc_cat2 = pow(TMath::Max(fabs(h2_ttbar_JESPlus->Integral() - h2_ttbar_base->Integral()), fabs(h2_ttbar_base->Integral() - h2_ttbar_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_ttbar_JERPlus->Integral() - h2_ttbar_base->Integral()), fabs(h2_ttbar_base->Integral() - h2_ttbar_JERMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_ttbar_METUCPlus->Integral() - h2_ttbar_base->Integral()), fabs(h2_ttbar_base->Integral() - h2_ttbar_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_ttbar_bTagPlus->Integral() - h2_ttbar_base->Integral()), fabs(h2_ttbar_base->Integral() - h2_ttbar_bTagMinus->Integral())), 2) + pow(h2_ttbar_base->Integral()*effUncTTbar, 2) ;

  
  outFile_incl <<"Total Bkg"<<" & "<<h1_TotalBkg->IntegralAndError(1,40,sError_totbkg)<<" $\\pm$ "<<sError_totbkg<<" $\\pm$ "<<sqrt( tmp_ttbar_tot_unc_cat1 + pow(TMath::Max(fabs(h1_wjet_JESPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wjet_METUCPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wjet_bTagPlus->Integral() - h1_wjet_base->Integral()),fabs(h1_wjet_base->Integral() - h1_wjet_bTagMinus->Integral())), 2) + pow(h1_wjet_base->Integral()*effUncWJets, 2) + pow(TMath::Max(fabs(h1_zjet_JESPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_zjet_METUCPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_zjet_bTagPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_bTagMinus->Integral())), 2) + pow(h1_zjet_base->Integral()*effUncZJets, 2) + pow(TMath::Max(fabs(h1_stop_JESPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_stop_METUCPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_stop_bTagPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_bTagMinus->Integral())), 2) + pow(h1_stop_base->Integral()*effUncSingletop, 2) ) << " \\\\ "<<endl;

  outFile_incl<<"\\hline "<<endl;
  outFile_incl<<"\\hline "<<endl;


  outFile<<"Total Bkg"<<" & "<<h1_TotalBkg->IntegralAndError(1,40,sError_totbkg)<<" $\\pm$ "<<sError_totbkg<<" $\\pm$ "<<sqrt( tmp_ttbar_tot_unc_cat1 + pow(TMath::Max(fabs(h1_wjet_JESPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wjet_METUCPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_wjet_bTagPlus->Integral() - h1_wjet_base->Integral()), fabs(h1_wjet_base->Integral() - h1_wjet_bTagMinus->Integral())), 2) + pow(h1_wjet_base->Integral()*effUncWJets, 2) + pow(TMath::Max(fabs(h1_zjet_JESPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_zjet_METUCPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_zjet_bTagPlus->Integral() - h1_zjet_base->Integral()), fabs(h1_zjet_base->Integral() - h1_zjet_bTagMinus->Integral())), 2) + pow(h1_zjet_base->Integral()*effUncZJets, 2) + pow(TMath::Max(fabs(h1_stop_JESPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_stop_METUCPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h1_stop_bTagPlus->Integral() - h1_stop_base->Integral()), fabs(h1_stop_base->Integral() - h1_stop_bTagMinus->Integral())), 2) + pow(h1_stop_base->Integral()*effUncSingletop, 2) )  <<" & "<<h2_TotalBkg->IntegralAndError(1,40,sError_totbkg)<<" $\\pm$ "<<sError_totbkg<<" $\\pm$ "<<sqrt(tmp_ttbar_tot_unc_cat2 +  pow(TMath::Max(fabs(h2_wjet_JESPlus->Integral() - h2_wjet_base->Integral()), fabs(h2_wjet_base->Integral() - h2_wjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_wjet_METUCPlus->Integral() - h2_wjet_base->Integral()), fabs(h2_wjet_base->Integral() - h2_wjet_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_wjet_bTagPlus->Integral() - h2_wjet_base->Integral()), fabs(h2_wjet_base->Integral() - h2_wjet_bTagMinus->Integral())), 2) + pow(h2_wjet_base->Integral()*effUncWJets, 2) + pow(TMath::Max(fabs(h2_zjet_JESPlus->Integral() - h2_zjet_base->Integral()), fabs(h2_zjet_base->Integral() - h2_zjet_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_zjet_METUCPlus->Integral() - h2_zjet_base->Integral()), fabs(h2_zjet_base->Integral() - h2_zjet_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_zjet_bTagPlus->Integral() - h2_zjet_base->Integral()), fabs(h2_zjet_base->Integral() - h2_zjet_bTagMinus->Integral())), 2) + pow(h2_zjet_base->Integral()*effUncZJets, 2) + pow(TMath::Max(fabs(h2_stop_JESPlus->Integral() - h2_stop_base->Integral()), fabs(h2_stop_base->Integral() - h2_stop_JESMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_stop_METUCPlus->Integral() - h2_stop_base->Integral()), fabs(h2_stop_base->Integral() - h2_stop_METUCMinus->Integral())), 2) + pow(TMath::Max(fabs(h2_stop_bTagPlus->Integral() - h2_stop_base->Integral()), fabs(h2_stop_base->Integral() - h2_stop_bTagMinus->Integral())), 2) + pow(h2_stop_base->Integral()*effUncSingletop, 2) )   <<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;   
  outFile<<"\\hline "<<endl;  

  outFile_incl<<"Data "<<" & "<<h1_data->Integral()<<" \\\\ "<<endl;
  outFile<<"Data "<<" & "<<h1_data->Integral()<<" & "<< h2_data->Integral()<<" \\\\ "<<endl; 

  outFile_incl<<"\\hline "<<endl;
  outFile_incl<<"\\hline "<<endl;
  outFile_incl<<"\\end{tabular}"<<endl;

  outFile<<"\\hline "<<endl;    
  outFile<<"\\hline "<<endl;   
  outFile<<"\\end{tabular}"<<endl; 
  // for unc table1
  outFile_unc1 << " JES+JER+MET "<<"&"<< 100*unc_WH_cat1 << " & " << 100*unc_ttbar_cat1 << " & " << 100*unc_wjets_cat1 << "&" << 100*unc_zjets_cat1 <<"&" << 100*unc_singletop_cat1 << "&" << 100*unc_diboson_cat1 << " \\\\ " << endl;
  outFile_unc1<<"\\hline "<<endl;
  outFile_unc1 << " b-jet tagging "<<"&"<< 100*bTagUnc_wh_cat1 << "&"<< 100*bTagUnc_ttbar_cat1 << "&"<< "-" << "&"<< "-" << "&"<< 100*bTagUnc_stop_cat1 << "&"<< "-" << " \\\\ " << endl;
  outFile_unc1<<"\\hline "<<endl;
  outFile_unc1 << "Jet$\\rightarrow$b mis-id" <<"&"<< "-" << "&" <<" - " << "&" << 100*bTagUnc_wjet_cat1 << "&" << 100*bTagUnc_zjet_cat1 << "&" << " - "<<"&" << 100*bTagUnc_diboson_cat1 << " \\\\ " << endl;
  outFile_unc1<<"\\hline "<<endl;
  outFile_unc1<< "Lepton selections "<<"&"<< "1.0" << "&"<< "1.0" <<"&"<< "1.0" <<"&"<< "1.0" <<"&"<< "1.0" <<"&"<< "1.0" <<" \\\\ " << endl;
  outFile_unc1<<"\\hline "<<endl;
  outFile_unc1<< "Cross-section" <<"&" << "6.0" <<"&" << "6.0" <<"&" << "5.0" <<"&" << "4.0" <<"&" << "6.0" <<"&" << "10.0" << " \\\\ " << endl;  
  outFile_unc1<<"\\hline "<<endl;
  outFile_unc1<< "MC Statistics" << "&" << 100*statUnc_wh_cat1 << "&" << 100*statUnc_ttbar_cat1 <<"&" <<100*statUnc_wjet_cat1 <<"&"<<100*statUnc_zjet_cat1 << "&" <<100*statUnc_stop_cat1 <<"&"<<100*statUnc_diboson_cat1<<  " \\\\ " << endl;
  outFile_unc1<<"\\hline "<<endl;
  //  outFile_unc1<<"Category migration" <<"&"<<"10.0" << "&" << "12.0" << "&" << "15.0" << "&" << "15.0" << "&" << "12.0" << "&" << "18.0" << " \\\\" << endl;
  //  outFile_unc1<<"\\hline "<<endl;
  outFile_unc1<< "Luminosity "<<"&"<< "2.6" << "&"<< "2.6" <<"&"<< "2.6" <<"&"<< "2.6" <<"&"<< "2.6" <<"&"<< "2.6" <<" \\\\ " << endl;
  outFile_unc1<<"\\hline "<<endl;
  outFile_unc1<<"\\hline "<<endl;
  outFile_unc1<<"\\end{tabular}"<<endl;

// for unc table2
  outFile_unc2 << " JES+JER+MET "<<"&"<< 100*unc_WH_cat2 << " & " << 100*unc_ttbar_cat2 << " & " << 100*unc_wjets_cat2 << "&" << 100*unc_zjets_cat2 <<"&" << 100*unc_singletop_cat2 << "&" << 100*unc_diboson_cat2 << " \\\\ " << endl;
  outFile_unc2<<"\\hline "<<endl;
  outFile_unc2 << " b-jet tagging "<<"&"<< 100*bTagUnc_wh_cat2 << "&"<< 100*bTagUnc_ttbar_cat2 << "&"<< "-" << "&"<< "-" << "&"<< 100*bTagUnc_stop_cat2 << "&"<< "-" << " \\\\ " << endl;
  outFile_unc2<<"\\hline "<<endl;
  outFile_unc2 << "Jet$\\rightarrow$b mis-id" <<"&"<< "-" << "&" <<" - " << "&" << 100*bTagUnc_wjet_cat2 << "&" << 100*bTagUnc_zjet_cat2 << "&" << " - "<<"&" << 100*bTagUnc_diboson_cat2 << " \\\\ " << endl;
  outFile_unc2<<"\\hline "<<endl;
  outFile_unc2<< "Lepton selections "<<"&"<< "1.0" << "&"<< "1.0" <<"&"<< "1.0" <<"&"<< "1.0" <<"&"<< "1.0" <<"&"<< "1.0" <<" \\\\ " << endl;
  outFile_unc2<<"\\hline "<<endl;
  outFile_unc2<< "Cross-section" <<"&" << "6.0" <<"&" << "6.0" <<"&" << "5.0" <<"&" << "4.0" <<"&" << "6.0" <<"&" << "10.0" << " \\\\ " << endl;
  outFile_unc2<<"\\hline "<<endl;
  outFile_unc2<< "MC Statictics" << "&" << 100*statUnc_wh_cat2 << "&" << 100*statUnc_ttbar_cat2 <<"&" <<100*statUnc_wjet_cat2 <<"&"<<100*statUnc_zjet_cat2 << "&" <<100*statUnc_stop_cat2 <<"&"<<100*statUnc_diboson_cat2<<  " \\\\ " << endl;
  outFile_unc2<<"\\hline "<<endl;
  outFile_unc2<<"Category migration" <<"&"<<"2.5" << "&" << "1.5" << "&" << "4.0" << "&" << "4.0" << "&" << "2.5" << "&" << "4.0" << " \\\\" << endl;
  outFile_unc2<<"\\hline "<<endl;
  outFile_unc2<< "Luminosity "<<"&"<< "2.6" << "&"<< "2.6" <<"&"<< "2.6" <<"&"<< "2.6" <<"&"<< "2.6" <<"&"<< "2.6" <<" \\\\ " << endl;
  outFile_unc2<<"\\hline "<<endl;
  outFile_unc2<<"\\hline "<<endl;
  outFile_unc2<<"\\end{tabular}"<<endl;

  if(NOTANPAS){
    outFile<<"\\end{LARGE}"<<endl;  
    outFile<<"\\end{center}"<<endl;  
    outFile<<"\\end{document}"<<endl; 
  }

  if(NOTANPAS){
    outFile_incl<<"\\end{LARGE}"<<endl;
    outFile_incl<<"\\end{center}"<<endl;
    outFile_incl<<"\\end{document}"<<endl;
  }
  // for incl unc table 
  if(NOTANPAS){
    outFile_unc_incl<<"\\end{LARGE}"<<endl;
    outFile_unc_incl<<"\\end{center}"<<endl;
    outFile_unc_incl<<"\\end{document}"<<endl;
  }


  // for unc table1
  if(NOTANPAS){
    outFile_unc1<<"\\end{LARGE}"<<endl;
    outFile_unc1<<"\\end{center}"<<endl;
    outFile_unc1<<"\\end{document}"<<endl;
  }

  // for unc table2
  if(NOTANPAS){
    outFile_unc2<<"\\end{LARGE}"<<endl;
    outFile_unc2<<"\\end{center}"<<endl;
    outFile_unc2<<"\\end{document}"<<endl;
  }

  outFile.close(); 
  outFile_incl.close();
  outFile_unc_incl.close();
  outFile_unc1.close();
  outFile_unc2.close();
  cout << "unc_WH_cat1 = " << unc_WH_cat1 << "  unc_WH_cat2=  " << unc_WH_cat2 << endl;
  cout << "unc_ttbar_cat1 =  " << unc_ttbar_cat1 <<"  unc_ttbar_cat2=  "<<unc_ttbar_cat2<< endl;
  cout << "unc_wjets_cat1 =  " << unc_wjets_cat1 << "unc_wjets_cat2=  "<<unc_wjets_cat2<<endl;
  cout << "unc_zjets_cat1 =  " << unc_zjets_cat1 << "unc_zjets_cat2=  "<<unc_zjets_cat2<<endl;
  cout << "unc_singletop_cat1 =  "<< unc_singletop_cat1 << "unc_singletop_cat2=  "<<unc_singletop_cat2<<endl;
  cout << "unc_diboson_cat1 =  " << unc_diboson_cat1 << "unc_diboson_cat2=  " << unc_diboson_cat2 << endl;
} 

void bkgestimation(TString hist_muplus="mjj_kfit", TString hist_muminus="mjj_kfit")  
{  
  TString inFile("/afs/cern.ch/work/g/gkole/chargedHiggs/8TeV/JERStusy/uncertainty/ExampleAnalysis/8TeV/2M/kinfit11/v3/test5/");
  
  TFile *wh = new TFile(inFile+"all_Hplus.root");  
  TH1F* h_plus_wh_base = ((TH1F*)wh->Get("base/"+hist_muplus))->Clone("h_plus_wh_base");  
  TH1F* h_plus_wh_JESPlus = ((TH1F*)wh->Get("JESPlus/"+hist_muplus) )->Clone("h_plus_wh_JESPlus"); 
  TH1F* h_plus_wh_JESMinus = ((TH1F*)wh->Get("JESMinus/"+hist_muplus) )->Clone("h_plus_wh_JESMinus");
  TH1F* h_minus_wh_base = ((TH1F*)wh->Get("base/"+hist_muminus))->Clone("h_minus_wh_base");
  TH1F* h_minus_wh_JESPlus = ((TH1F*)wh->Get("JESPlus/"+hist_muminus) )->Clone("h_minus_wh_JESPlus");
  TH1F* h_minus_wh_JESMinus = ((TH1F*)wh->Get("JESMinus/"+hist_muminus) )->Clone("h_minus_wh_JESMinus");

  TFile *ttbar = new TFile(inFile+"all_TTJetsP.root");  
  TH1F* h_plus_ttbar_base = ((TH1F*)ttbar->Get("base/"+hist_muplus) )->Clone("h_plus_ttbar_base");   
  TH1F* h_plus_ttbar_JESPlus = ((TH1F*)ttbar->Get("JESPlus/"+hist_muplus) )->Clone("h_plus_ttbar_JESPlus");  
  TH1F* h_plus_ttbar_JESMinus = ((TH1F*)ttbar->Get("JESMinus/"+hist_muplus) )->Clone("h_plus_ttbar_JESMinus");
  TH1F* h_minus_ttbar_base = ((TH1F*)ttbar->Get("base/"+hist_muminus) )->Clone("h_minus_ttbar_base");
  TH1F* h_minus_ttbar_JESPlus = ((TH1F*)ttbar->Get("JESPlus/"+hist_muminus) )->Clone("h_minus_ttbar_JESPlus");
  TH1F* h_minus_ttbar_JESMinus = ((TH1F*)ttbar->Get("JESMinus/"+hist_muminus) )->Clone("h_minus_ttbar_JESMinus");


  TFile *wjet = new TFile(inFile+"all_WJets.root");
  TH1F* h_plus_wjet_base = ((TH1F*)wjet->Get("base/"+hist_muplus) )->Clone("h_plus_wjet_base");
  TH1F* h_plus_wjet_JESPlus = ((TH1F*)wjet->Get("JESPlus/"+hist_muplus) )->Clone("h_plus_wjet_JESPlus");
  TH1F* h_plus_wjet_JESMinus = ((TH1F*)wjet->Get("JESMinus/"+hist_muplus) )->Clone("h_plus_wjet_JESMinus");
  TH1F* h_minus_wjet_base = ((TH1F*)wjet->Get("base/"+hist_muminus) )->Clone("h_minus_wjet_base");
  TH1F* h_minus_wjet_JESPlus = ((TH1F*)wjet->Get("JESPlus/"+hist_muminus) )->Clone("h_minus_wjet_JESPlus");
  TH1F* h_minus_wjet_JESMinus = ((TH1F*)wjet->Get("JESMinus/"+hist_muminus) )->Clone("h_minus_wjet_JESMinus");


  TFile *zjet = new TFile(inFile+"all_DY.root");
  TH1F* h_plus_zjet_base = ((TH1F*)zjet->Get("base/"+hist_muplus) )->Clone("h_plus_zjet_base");
  TH1F* h_plus_zjet_JESPlus = ((TH1F*)zjet->Get("JESPlus/"+hist_muplus) )->Clone("h_plus_zjet_JESPlus");
  TH1F* h_plus_zjet_JESMinus = ((TH1F*)zjet->Get("JESMinus/"+hist_muplus) )->Clone("h_plus_zjet_JESMinus");
  TH1F* h_minus_zjet_base = ((TH1F*)zjet->Get("base/"+hist_muminus) )->Clone("h_minus_zjet_base");
  TH1F* h_minus_zjet_JESPlus = ((TH1F*)zjet->Get("JESPlus/"+hist_muminus) )->Clone("h_minus_zjet_JESPlus");
  TH1F* h_minus_zjet_JESMinus = ((TH1F*)zjet->Get("JESMinus/"+hist_muminus) )->Clone("h_minus_zjet_JESMinus");

  TFile *qcd = new TFile(inFile+"all_QCD.root");
  TH1F* h_plus_qcd_base = ((TH1F*)qcd->Get("base/"+hist_muplus) )->Clone("h_plus_qcd_base");
  TH1F* h_plus_qcd_JESPlus = ((TH1F*)qcd->Get("JESPlus/"+hist_muplus) )->Clone("h_plus_qcd_JESPlus");
  TH1F* h_plus_qcd_JESMinus = ((TH1F*)qcd->Get("JESMinus/"+hist_muplus) )->Clone("h_plus_qcd_JESMinus");
  TH1F* h_minus_qcd_base = ((TH1F*)qcd->Get("base/"+hist_muminus) )->Clone("h_minus_qcd_base");
  TH1F* h_minus_qcd_JESPlus = ((TH1F*)qcd->Get("JESPlus/"+hist_muminus) )->Clone("h_minus_qcd_JESPlus");
  TH1F* h_minus_qcd_JESMinus = ((TH1F*)qcd->Get("JESMinus/"+hist_muminus) )->Clone("h_minus_qcd_JESMinus");


  TFile *stop = new TFile(inFile+"all_ST.root");
  TH1F* h_plus_stop_base = ((TH1F*)stop->Get("base/"+hist_muplus) )->Clone("h_plus_stop_base");
  TH1F* h_plus_stop_JESPlus = ((TH1F*)stop->Get("JESPlus/"+hist_muplus) )->Clone("h_plus_stop_JESPlus");
  TH1F* h_plus_stop_JESMinus = ((TH1F*)stop->Get("JESMinus/"+hist_muplus) )->Clone("h_plus_stop_JESMinus");
  TH1F* h_minus_stop_base = ((TH1F*)stop->Get("base/"+hist_muminus) )->Clone("h_minus_stop_base");
  TH1F* h_minus_stop_JESPlus = ((TH1F*)stop->Get("JESPlus/"+hist_muminus) )->Clone("h_minus_stop_JESPlus");
  TH1F* h_minus_stop_JESMinus = ((TH1F*)stop->Get("JESMinus/"+hist_muminus) )->Clone("h_minus_stop_JESMinus");


  TFile *diboson = new TFile(inFile+"all_VV.root");
  TH1F* h_plus_diboson_base = ((TH1F*)diboson->Get("base/"+hist_muplus) )->Clone("h_plus_diboson_base");
  TH1F* h_plus_diboson_JESPlus = ((TH1F*)diboson->Get("JESPlus/"+hist_muplus) )->Clone("h_plus_diboson_JESPlus");
  TH1F* h_plus_diboson_JESMinus = ((TH1F*)diboson->Get("JESMinus/"+hist_muplus) )->Clone("h_plus_diboson_JESMinus");
  TH1F* h_minus_diboson_base = ((TH1F*)diboson->Get("base/"+hist_muminus) )->Clone("h_minus_diboson_base");
  TH1F* h_minus_diboson_JESPlus = ((TH1F*)diboson->Get("JESPlus/"+hist_muminus) )->Clone("h_minus_diboson_JESPlus");
  TH1F* h_minus_diboson_JESMinus = ((TH1F*)diboson->Get("JESMinus/"+hist_muminus) )->Clone("h_minus_diboson_JESMinus");

  TH1F* h_plus_TotalBkg = h_plus_wjet_base->Clone("h_plus_TotalBkg");
  h_plus_TotalBkg->Reset();
  h_plus_TotalBkg->Add(h_plus_ttbar_base);
  h_plus_TotalBkg->Add(h_plus_wjet_base);
  h_plus_TotalBkg->Add(h_plus_zjet_base);
  h_plus_TotalBkg->Add(h_plus_qcd_base);
  h_plus_TotalBkg->Add(h_plus_stop_base);
  h_plus_TotalBkg->Add(h_plus_diboson_base);
  TH1F* h_minus_TotalBkg = h_minus_wjet_base->Clone("h_minus_TotalBkg");
  h_minus_TotalBkg->Reset();
  h_minus_TotalBkg->Add(h_minus_ttbar_base);
  h_minus_TotalBkg->Add(h_minus_wjet_base);
  h_minus_TotalBkg->Add(h_minus_zjet_base);
  h_minus_TotalBkg->Add(h_minus_qcd_base);
  h_minus_TotalBkg->Add(h_minus_stop_base);
  h_minus_TotalBkg->Add(h_minus_diboson_base);

  TFile *data = new TFile(inFile+"data_selection_.root");
  TH1F* h_plus_data = ((TH1F*)data->Get("base/"+hist_muplus) )->Clone("h_plus_data");
  TH1F* h_minus_data = ((TH1F*)data->Get("base/"+hist_muminus) )->Clone("h_minus_data");


  double sError_wh = 0.0;
  double sError_ttbar = 0.0;
  double sError_wjets = 0.0;
  double sError_zjets = 0.0;
  double sError_qcd = 0.0;
  double sError_stop = 0.0;
  double sError_diboson = 0.0;
  double sError_totbkg = 0.0;

  ofstream outFile;
  outFile.open("bkgestimation_cutflow.tex");
  outFile<<"\\documentclass[landscape,letterpaper]{article}"<<endl;
  outFile<<"\\pagestyle{empty}"<<endl;
  outFile<<"\\usepackage{epsfig}"<<endl;
  outFile<<"\\usepackage{amsmath}"<<endl;
  outFile<<"\\usepackage{array}"<<endl;
  outFile<<"\\usepackage{multirow}"<<endl;   
  outFile<<"\\usepackage[cm]{fullpage}"<<endl;
  outFile<<"\\textheight = 8.in"<<endl; 
  outFile<<"\\textwidth 7.0in"<<endl;   
  outFile<<"\\begin{document}"<<endl;  
  outFile<<"\\begin{center}"<<endl;   
  outFile<<"\\begin{LARGE}"<<endl;   
  outFile<<"\\begin{tabular}{ | c| c| c| }"<<endl;
  outFile<<"\\multicolumn{3}{c}{ } \\\\"<<endl;
  outFile<<"\\hline "<<endl;  
  outFile<<"\\multicolumn{1}{|c|}{Source} & \\multicolumn{1}{|c|}{N$_{\\mu^{+}}$ events} & \\multicolumn{1}{|c|}{N$_{\\mu^{-}}$ events} \\\\"<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<"\\hline "<<endl; 

  outFile<<"SM TTbar"<<" & " << h_plus_ttbar_base->IntegralAndError(1,40,sError_ttbar) << "$\\pm$ " << sError_ttbar <<" & "<< h_minus_ttbar_base->IntegralAndError(1,40,sError_ttbar) << "$\\pm$ " << sError_ttbar<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"W + Jets"<<" & " << h_plus_wjet_base->IntegralAndError(1,40,sError_wjets) << "$\\pm$ " << sError_wjets  <<" & "<< h_minus_wjet_base->IntegralAndError(1,40,sError_wjets) << "$\\pm$ " << sError_wjets<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"Z + Jets"<<" & " << h_plus_zjet_base->IntegralAndError(1,40,sError_zjets) << "$\\pm$ " << sError_zjets <<" & "<< h_minus_zjet_base->IntegralAndError(1,40,sError_zjets) << "$\\pm$ " << sError_zjets<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"QCD"<<" & " << h_plus_qcd_base->IntegralAndError(1,40,sError_qcd) << "$\\pm$ " << sError_qcd <<" & "<< h_minus_qcd_base->IntegralAndError(1,40,sError_qcd) << "$\\pm$ " << sError_qcd<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"Single Top"<<" & " << h_plus_stop_base->IntegralAndError(1,40,sError_stop) << "$\\pm$ " << sError_stop <<" & "<< h_minus_stop_base->IntegralAndError(1,40,sError_stop) << "$\\pm$ " << sError_stop<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"Diboson"<<" & " << h_plus_diboson_base->IntegralAndError(1,40,sError_diboson) << "$\\pm$ " << sError_diboson <<" & "<< h_minus_diboson_base->IntegralAndError(1,40,sError_diboson) << "$\\pm$ " << sError_diboson<<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;

  outFile<<"Total Bkg"<<" & "<<h_plus_TotalBkg->IntegralAndError(1,40,sError_totbkg)<<" $\\pm$ "<<sError_totbkg<<" & "<<h_minus_TotalBkg->IntegralAndError(1,40,sError_totbkg)<<" $\\pm$ "<<sError_totbkg <<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"Data "<<" & "<<h_plus_data->Integral()<<" $\\pm$ "<< sqrt(h_plus_data->Integral())<<" & " << h_minus_data->Integral()<<" $\\pm$ "<<sqrt(h_minus_data->Integral())<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\end{tabular}"<<endl;
  outFile<<"\\begin{itemize}"<<endl;
  double wPdata = h_plus_data->Integral() - h_plus_stop_base->Integral() ;
  double wMdata = h_minus_data->Integral() - h_minus_stop_base->Integral() ;
  double scale_factor = (wPdata - wMdata)/(h_plus_wjet_base->Integral() - h_minus_wjet_base->Integral()); 
  outFile<<"\\item wPdata $(Data^{+} - singletop^{+})$ = " << wPdata << endl;
  outFile<<"\\item wMdata $(Data^{-} - singletop^{-})$ = " << wMdata << endl;
  outFile<<"\\item scale factor ($\\frac{wPdata - wMdata}{W^+_{mc} - W^-_{mc}}$) = " << scale_factor << endl;
  outFile<<"\\end{itemize}"<<endl;
  outFile<<"\\end{LARGE}"<<endl;
  outFile<<"\\end{center}"<<endl;
  outFile<<"\\end{document}"<<endl;
  
  outFile.close();
}

