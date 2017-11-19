#include <iostream>
#include <fstream>
#include <iomanip>

void makeCutFlowTable(bool NOTANPAS = false)
{
  TString inFile("$PWD/");

  cout << "inFile  " << inFile << endl;

  TFile *ttbar = new TFile(inFile+"all_TTJetsP.root");
  TH1F* h_ttbar = ((TH1F*)ttbar->Get("base/Iso/cutflow") )->Clone("h_ttbar");

  TH1F* httbar_toppt = (TH1F*)ttbar->Get("base/SF_topPtWeights"); 
  double avgTop = httbar_toppt->GetMean();
  double wh_scale = 0.32; // 2x(1-x) assuming x = 0.2  
  h_ttbar->Scale(1.0/avgTop);

  ofstream outFile;
  outFile.open("muon_cutflow.tex");
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
  //outFile<<"%\\end{LARGE}"<<endl; 
  outFile<<"\\end{center}"<<endl; 
  if(NOTANPAS){
    outFile<<"\\end{document}"<<endl;
  }
  outFile.close();

}
string printCutFlowHist(string histname, TH1F *h, int Nbin, double sf){
   string allBinValue = histname;
   for(int i =1; i<Nbin; i++){
     ostringstream convert;
     string Result("");
     //allBinValue += "h->GetBinContent(i)" + " & ";
     convert <<sf*h->GetBinContent(i);
     Result = convert.str(); 
     allBinValue += " & "+Result;
   }
   return allBinValue ;
}

void makeCutFlowTableAny(ofstream &outFile, string sys, int Nbin, 
		TH1F *h_wh,
		TH1F *h_ttbar,
		TH1F *h_stop,
		TH1F *h_wjets,
		TH1F *h_dyjets,
		TH1F *h_qcd,
		TH1F *h_st,
		TH1F *h_vv,
		TH1F *h_data){
  outFile<<""<<endl;
  outFile<<"%%%%%%%%%%%%%%%%%% NEXT TABLE: "+sys+" %%%%%%%%%%%%%%%%%"<<endl;
  outFile<<"\\begin{table}"<<endl; 
  outFile<<"\\begin{center}"<<endl; 
  //outFile<<"%\\begin{LARGE}"<<endl;
  outFile<<"\\footnotesize\\setlength{\\tabcolsep}{4.5pt}"<<endl;
  outFile<<"\\begin{tabular}{ | c| c| c| c| c| c| c| c| c| c| c| c| c|}"<<endl; 
  outFile<<"\\multicolumn{5}{c}{ } \\\\"<<endl; 
  outFile<<"\\hline "<<endl;
  outFile<<"\\multicolumn{1}{| c|}{"+ sys +" } & \\multicolumn{1}{ c|}{ $HLT\\_IsoMu24$ } & \\multicolumn{1}{ c|}{ $N_{muon}=1$ } & \\multicolumn{1}{ c|}{ $N_{ele}=0$} & \\multicolumn{1}{ c|}{ Muon SF } & \\multicolumn{1}{ c|}{ $\\mu^{Iso} < 0.15$ } &\\multicolumn{1}{ c|}{ $N_{jets}\\ge 4$ } & \\multicolumn{1}{ c|}{ $\\not\\!\\!E_T \\ge 20GeV$ }&  \\multicolumn{1}{ c|}{MT $\\ge$ 20 GeV} & \\multicolumn{1}{ c |}{ $\\ge$ 2btag }& \\multicolumn{1}{ c |}{ BTag SF }& \\multicolumn{1}{ c|}{Kin Fit } & \\multicolumn{1}{ c|}{KinFit cuts}  \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  TH1F* h_MC = (TH1F*)h_ttbar->Clone("Bkg MC");
  h_MC->Add(h_stop);
  h_MC->Add(h_wjets);
  h_MC->Add(h_dyjets);
  h_MC->Add(h_qcd);
  h_MC->Add(h_vv);
  
  TH1F* h_ratio = (TH1F*)h_data->Clone("ratio");
  h_ratio->Divide(h_MC);

  outFile<<printCutFlowHist("WH, M120", 	h_wh, 		Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<printCutFlowHist("TTJets", 	h_ttbar, 	Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<printCutFlowHist("STop", 	h_stop, 	Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<printCutFlowHist("WJets", 	h_wjets, 	Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<printCutFlowHist("DYJets", 	h_dyjets, 	Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<printCutFlowHist("QCD", 	h_qcd, 		Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<printCutFlowHist("VV", 	h_vv, 		Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<printCutFlowHist("Bkg", 	h_MC, 		Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<printCutFlowHist("Data", 	h_data, 	Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<printCutFlowHist("Data/Bkg", 	h_ratio, 	Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\end{tabular}"<<endl;
  //outFile<<"%\\end{LARGE}"<<endl; 
  outFile<<"\\end{center}"<<endl; 
  outFile<<"\\\caption{Number of evets after various cuts for sys: " +sys+ "}"<<endl; 
  outFile<<"\\end{table}"<<endl; 
}

void getSignalCutFlowTable(ofstream &outFile, string sys, int Nbin, 
		TH1F *h_wh80,
		TH1F *h_wh90,
		TH1F *h_wh100,
		TH1F *h_wh120,
		TH1F *h_wh140,
		TH1F *h_wh150,
		TH1F *h_wh155,
		TH1F *h_wh160
		){
  outFile<<""<<endl;
  outFile<<"%%%%%%%%%%%%%%%%%% NEXT TABLE: "+sys+" %%%%%%%%%%%%%%%%%"<<endl;
  outFile<<"\\begin{table}"<<endl; 
  outFile<<"\\begin{center}"<<endl; 
  //outFile<<"%\\begin{LARGE}"<<endl;
  outFile<<"\\footnotesize\\setlength{\\tabcolsep}{4.5pt}"<<endl;
  outFile<<"\\begin{tabular}{ | c| c| c| c| c| c| c| c| c| c| c| c| c|}"<<endl; 
  outFile<<"\\multicolumn{5}{c}{ } \\\\"<<endl; 
  outFile<<"\\hline "<<endl;
  outFile<<"\\multicolumn{1}{| c|}{"+ sys +" } & \\multicolumn{1}{ c|}{ $HLT\\_IsoMu24$ } & \\multicolumn{1}{ c|}{ $N_{muon}=1$ } & \\multicolumn{1}{ c|}{ $N_{ele}=0$} & \\multicolumn{1}{ c|}{ Muon SF } & \\multicolumn{1}{ c|}{ $\\mu^{Iso} < 0.15$ } &\\multicolumn{1}{ c|}{ $N_{jets}\\ge 4$ } & \\multicolumn{1}{ c|}{ $\\not\\!\\!E_T \\ge 20GeV$ }&  \\multicolumn{1}{ c|}{MT $\\ge$ 20 GeV} & \\multicolumn{1}{ c |}{ $\\ge$ 2btag }& \\multicolumn{1}{ c |}{ BTag SF }& \\multicolumn{1}{ c|}{Kin Fit } & \\multicolumn{1}{ c|}{KinFit cuts}  \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;

  outFile<<printCutFlowHist("WH, M80", 	h_wh80, 		Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<printCutFlowHist("WH, M90",  h_wh90, 		Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  outFile<<printCutFlowHist("WH, M100",  h_wh100, 		Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  outFile<<printCutFlowHist("WH, M120",  h_wh120, 		Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  outFile<<printCutFlowHist("WH, M140",  h_wh140, 		Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  outFile<<printCutFlowHist("WH, M150",  h_wh150, 		Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  outFile<<printCutFlowHist("WH, M155",  h_wh155, 		Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  outFile<<printCutFlowHist("WH, M160",  h_wh160, 		Nbin, 1)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\end{tabular}"<<endl;
  //outFile<<"%\\end{LARGE}"<<endl; 
  outFile<<"\\end{center}"<<endl; 
  outFile<<"\\\caption{Number of evets after various cuts for sys: " +sys+ "}"<<endl; 
  outFile<<"\\end{table}"<<endl; 
}

void makeBaseTable(TString histSubDir="Iso/", TString histName="cutflow"){  
  
  TString inFile("$PWD/");
  TString histPath(histSubDir+histName);
  TFile *wh 	  		= new TFile(inFile+"all_Hplus120.root"); 
  TFile *ttbar    		= new TFile(inFile+"all_TTJetsP.root"); 
  TFile *wjet  			= new TFile(inFile+"all_WJets.root"); 
  TFile *zjet  			= new TFile(inFile+"all_DY.root");
  TFile *qcd  			= new TFile(inFile+"all_QCD.root");
  TFile *stop  			= new TFile(inFile+"all_ST.root");
  TFile *diboson 		= new TFile(inFile+"all_VV.root");
  
  TH1F* h_wh_base  		= ((TH1F*)wh->Get("base/"+histPath))->Clone("h_wh_base"); 
  TH1F* h_ttbar_base  		= ((TH1F*)ttbar->Get("base/"+histPath) )->Clone("h_ttbar_base");  
  TH1F* h_wjet_base  		= ((TH1F*)wjet->Get("base/"+histPath) )->Clone("h_wjet_base");  
  TH1F* h_zjet_base  		= ((TH1F*)zjet->Get("base/"+histPath) )->Clone("h_zjet_base");
  TH1F* h_qcd_base 		= ((TH1F*)qcd->Get("base/"+histPath) )->Clone("h_qcd_base");
  TH1F* h_stop_base  		= ((TH1F*)stop->Get("base/"+histPath) )->Clone("h_stop_base");
  TH1F* h_diboson_base 		= ((TH1F*)diboson->Get("base/"+histPath) )->Clone("h_diboson_base");

  //get average of top pt-reweighting 
  TH1F* hsig_avgTop120 = (TH1F*)wh->Get("base/SF_topPtWeights"); 
  double sig_avgTop120 = hsig_avgTop120->GetMean();
  double sf120 = 1;
  h_wh_base->Scale(sf120/sig_avgTop120);
  
  TH1F* httbar_avgTop = (TH1F*)ttbar->Get("base/SF_topPtWeights"); 
  double ttbar_avgTop = httbar_avgTop->GetMean();
  double sf_ttbar = 1.0;
  h_ttbar_base->Scale(sf_ttbar/ttbar_avgTop);

  TH1F* h_TotalBkg = h_wjet_base->Clone("h_TotalBkg");
  h_TotalBkg->Reset();
  h_TotalBkg->Add(h_ttbar_base);
  h_TotalBkg->Add(h_wjet_base);
  h_TotalBkg->Add(h_zjet_base);
  h_TotalBkg->Add(h_qcd_base);
  h_TotalBkg->Add(h_stop_base);
  h_TotalBkg->Add(h_diboson_base);

  TFile *data = new TFile(inFile+"all_muData.root");
  TH1F* h_data = ((TH1F*)data->Get("base/"+histPath))->Clone("h_data");

  ofstream outFile; 
  outFile.open("baseCutFlowTable.tex"); 
   
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
  //Add another table with JESUP
  double Nbin = 13;
  makeCutFlowTableAny(outFile, "base", Nbin, h_wh_base,	h_ttbar_base, h_stop_base, h_wjet_base, h_zjet_base, h_qcd_base, h_stop_base, h_diboson_base, h_data);
  outFile<<"\\end{document}"<<endl;  
  outFile.close(); 
} 
  
void makeSignalTable(TString histSubDir="Iso/", TString histName="cutflow"){  
  
  TString inFile("$PWD/");
  TString histPath(histSubDir+histName);
  TFile *wh80   		= new TFile(inFile+"all_Hplus80.root"); 
  TFile *wh90   		= new TFile(inFile+"all_Hplus90.root"); 
  TFile *wh100  		= new TFile(inFile+"all_Hplus100.root"); 
  TFile *wh120  		= new TFile(inFile+"all_Hplus120.root"); 
  TFile *wh140  		= new TFile(inFile+"all_Hplus140.root"); 
  TFile *wh150  		= new TFile(inFile+"all_Hplus150.root"); 
  TFile *wh155  		= new TFile(inFile+"all_Hplus155.root"); 
  TFile *wh160  		= new TFile(inFile+"all_Hplus160.root"); 
  
  //Hplus  M80   signal 
  TH1F* h_wh80_base  			= ((TH1F*)wh80->Get("base/"+histPath))->Clone("h_wh_base"); 
  TH1F* h_wh80_JESPlus 			= ((TH1F*)wh80->Get("JESPlus/"+histPath) )->Clone("h_wh_JESPlus");
  TH1F* h_wh80_JESMinus 		= ((TH1F*)wh80->Get("JESMinus/"+histPath) )->Clone("h_wh_JESMinus");
  TH1F* h_wh80_JERPlus 			= ((TH1F*)wh80->Get("JERPlus/"+histPath) )->Clone("h_wh_JERPlus");
  TH1F* h_wh80_JERMinus 		= ((TH1F*)wh80->Get("JERMinus/"+histPath) )->Clone("h_wh_JERMinus");
  TH1F* h_wh80_METUCPlus 		= ((TH1F*)wh80->Get("METUCPlus/"+histPath) )->Clone("h_wh_METUCPlus");
  TH1F* h_wh80_METUCMinus 		= ((TH1F*)wh80->Get("METUCMinus/"+histPath) )->Clone("h_wh_METUCMinus");
  TH1F* h_wh80_bTagPlus 		= ((TH1F*)wh80->Get("bTagPlus/"+histPath) )->Clone("h_wh_bTagPlus");
  TH1F* h_wh80_bTagMinus 		= ((TH1F*)wh80->Get("bTagMinus/"+histPath) )->Clone("h_wh_bTagMinus");
  
  //Hplus  M90   signal 
  TH1F* h_wh90_base  			= ((TH1F*)wh90->Get("base/"+histPath))->Clone("h_wh_base"); 
  TH1F* h_wh90_JESPlus 			= ((TH1F*)wh90->Get("JESPlus/"+histPath) )->Clone("h_wh_JESPlus");
  TH1F* h_wh90_JESMinus 		= ((TH1F*)wh90->Get("JESMinus/"+histPath) )->Clone("h_wh_JESMinus");
  TH1F* h_wh90_JERPlus 			= ((TH1F*)wh90->Get("JERPlus/"+histPath) )->Clone("h_wh_JERPlus");
  TH1F* h_wh90_JERMinus 		= ((TH1F*)wh90->Get("JERMinus/"+histPath) )->Clone("h_wh_JERMinus");
  TH1F* h_wh90_METUCPlus 		= ((TH1F*)wh90->Get("METUCPlus/"+histPath) )->Clone("h_wh_METUCPlus");
  TH1F* h_wh90_METUCMinus 		= ((TH1F*)wh90->Get("METUCMinus/"+histPath) )->Clone("h_wh_METUCMinus");
  TH1F* h_wh90_bTagPlus 		= ((TH1F*)wh90->Get("bTagPlus/"+histPath) )->Clone("h_wh_bTagPlus");
  TH1F* h_wh90_bTagMinus 		= ((TH1F*)wh90->Get("bTagMinus/"+histPath) )->Clone("h_wh_bTagMinus");
  
  //Hplus  M100   signal 
  TH1F* h_wh100_base  			= ((TH1F*)wh100->Get("base/"+histPath))->Clone("h_wh_base"); 
  TH1F* h_wh100_JESPlus 		= ((TH1F*)wh100->Get("JESPlus/"+histPath) )->Clone("h_wh_JESPlus");
  TH1F* h_wh100_JESMinus 		= ((TH1F*)wh100->Get("JESMinus/"+histPath) )->Clone("h_wh_JESMinus");
  TH1F* h_wh100_JERPlus 		= ((TH1F*)wh100->Get("JERPlus/"+histPath) )->Clone("h_wh_JERPlus");
  TH1F* h_wh100_JERMinus 		= ((TH1F*)wh100->Get("JERMinus/"+histPath) )->Clone("h_wh_JERMinus");
  TH1F* h_wh100_METUCPlus 		= ((TH1F*)wh100->Get("METUCPlus/"+histPath) )->Clone("h_wh_METUCPlus");
  TH1F* h_wh100_METUCMinus 		= ((TH1F*)wh100->Get("METUCMinus/"+histPath) )->Clone("h_wh_METUCMinus");
  TH1F* h_wh100_bTagPlus 		= ((TH1F*)wh100->Get("bTagPlus/"+histPath) )->Clone("h_wh_bTagPlus");
  TH1F* h_wh100_bTagMinus 		= ((TH1F*)wh100->Get("bTagMinus/"+histPath) )->Clone("h_wh_bTagMinus");
  
  
  //Hplus  M120   signal 
  TH1F* h_wh120_base  			= ((TH1F*)wh120->Get("base/"+histPath))->Clone("h_wh_base"); 
  TH1F* h_wh120_JESPlus 		= ((TH1F*)wh120->Get("JESPlus/"+histPath) )->Clone("h_wh_JESPlus");
  TH1F* h_wh120_JESMinus 		= ((TH1F*)wh120->Get("JESMinus/"+histPath) )->Clone("h_wh_JESMinus");
  TH1F* h_wh120_JERPlus 		= ((TH1F*)wh120->Get("JERPlus/"+histPath) )->Clone("h_wh_JERPlus");
  TH1F* h_wh120_JERMinus 		= ((TH1F*)wh120->Get("JERMinus/"+histPath) )->Clone("h_wh_JERMinus");
  TH1F* h_wh120_METUCPlus 		= ((TH1F*)wh120->Get("METUCPlus/"+histPath) )->Clone("h_wh_METUCPlus");
  TH1F* h_wh120_METUCMinus 		= ((TH1F*)wh120->Get("METUCMinus/"+histPath) )->Clone("h_wh_METUCMinus");
  TH1F* h_wh120_bTagPlus 		= ((TH1F*)wh120->Get("bTagPlus/"+histPath) )->Clone("h_wh_bTagPlus");
  TH1F* h_wh120_bTagMinus 		= ((TH1F*)wh120->Get("bTagMinus/"+histPath) )->Clone("h_wh_bTagMinus");
  
  
  //Hplus  M140   signal 
  TH1F* h_wh140_base  			= ((TH1F*)wh140->Get("base/"+histPath))->Clone("h_wh_base"); 
  TH1F* h_wh140_JESPlus 		= ((TH1F*)wh140->Get("JESPlus/"+histPath) )->Clone("h_wh_JESPlus");
  TH1F* h_wh140_JESMinus 		= ((TH1F*)wh140->Get("JESMinus/"+histPath) )->Clone("h_wh_JESMinus");
  TH1F* h_wh140_JERPlus 		= ((TH1F*)wh140->Get("JERPlus/"+histPath) )->Clone("h_wh_JERPlus");
  TH1F* h_wh140_JERMinus 		= ((TH1F*)wh140->Get("JERMinus/"+histPath) )->Clone("h_wh_JERMinus");
  TH1F* h_wh140_METUCPlus 		= ((TH1F*)wh140->Get("METUCPlus/"+histPath) )->Clone("h_wh_METUCPlus");
  TH1F* h_wh140_METUCMinus 		= ((TH1F*)wh140->Get("METUCMinus/"+histPath) )->Clone("h_wh_METUCMinus");
  TH1F* h_wh140_bTagPlus 		= ((TH1F*)wh140->Get("bTagPlus/"+histPath) )->Clone("h_wh_bTagPlus");
  TH1F* h_wh140_bTagMinus 		= ((TH1F*)wh140->Get("bTagMinus/"+histPath) )->Clone("h_wh_bTagMinus");
  
  
  //Hplus  M150   signal 
  TH1F* h_wh150_base  			= ((TH1F*)wh150->Get("base/"+histPath))->Clone("h_wh_base"); 
  TH1F* h_wh150_JESPlus 		= ((TH1F*)wh150->Get("JESPlus/"+histPath) )->Clone("h_wh_JESPlus");
  TH1F* h_wh150_JESMinus 		= ((TH1F*)wh150->Get("JESMinus/"+histPath) )->Clone("h_wh_JESMinus");
  TH1F* h_wh150_JERPlus 		= ((TH1F*)wh150->Get("JERPlus/"+histPath) )->Clone("h_wh_JERPlus");
  TH1F* h_wh150_JERMinus 		= ((TH1F*)wh150->Get("JERMinus/"+histPath) )->Clone("h_wh_JERMinus");
  TH1F* h_wh150_METUCPlus 		= ((TH1F*)wh150->Get("METUCPlus/"+histPath) )->Clone("h_wh_METUCPlus");
  TH1F* h_wh150_METUCMinus 		= ((TH1F*)wh150->Get("METUCMinus/"+histPath) )->Clone("h_wh_METUCMinus");
  TH1F* h_wh150_bTagPlus 		= ((TH1F*)wh150->Get("bTagPlus/"+histPath) )->Clone("h_wh_bTagPlus");
  TH1F* h_wh150_bTagMinus 		= ((TH1F*)wh150->Get("bTagMinus/"+histPath) )->Clone("h_wh_bTagMinus");
  
  
  //Hplus  M155   signal 
  TH1F* h_wh155_base  			= ((TH1F*)wh155->Get("base/"+histPath))->Clone("h_wh_base"); 
  TH1F* h_wh155_JESPlus 		= ((TH1F*)wh155->Get("JESPlus/"+histPath) )->Clone("h_wh_JESPlus");
  TH1F* h_wh155_JESMinus 		= ((TH1F*)wh155->Get("JESMinus/"+histPath) )->Clone("h_wh_JESMinus");
  TH1F* h_wh155_JERPlus 		= ((TH1F*)wh155->Get("JERPlus/"+histPath) )->Clone("h_wh_JERPlus");
  TH1F* h_wh155_JERMinus 		= ((TH1F*)wh155->Get("JERMinus/"+histPath) )->Clone("h_wh_JERMinus");
  TH1F* h_wh155_METUCPlus 		= ((TH1F*)wh155->Get("METUCPlus/"+histPath) )->Clone("h_wh_METUCPlus");
  TH1F* h_wh155_METUCMinus 		= ((TH1F*)wh155->Get("METUCMinus/"+histPath) )->Clone("h_wh_METUCMinus");
  TH1F* h_wh155_bTagPlus 		= ((TH1F*)wh155->Get("bTagPlus/"+histPath) )->Clone("h_wh_bTagPlus");
  TH1F* h_wh155_bTagMinus 		= ((TH1F*)wh155->Get("bTagMinus/"+histPath) )->Clone("h_wh_bTagMinus");
  
  
  //Hplus  M160   signal 
  TH1F* h_wh160_base  			= ((TH1F*)wh160->Get("base/"+histPath))->Clone("h_wh_base"); 
  TH1F* h_wh160_JESPlus 		= ((TH1F*)wh160->Get("JESPlus/"+histPath) )->Clone("h_wh_JESPlus");
  TH1F* h_wh160_JESMinus 		= ((TH1F*)wh160->Get("JESMinus/"+histPath) )->Clone("h_wh_JESMinus");
  TH1F* h_wh160_JERPlus 		= ((TH1F*)wh160->Get("JERPlus/"+histPath) )->Clone("h_wh_JERPlus");
  TH1F* h_wh160_JERMinus 		= ((TH1F*)wh160->Get("JERMinus/"+histPath) )->Clone("h_wh_JERMinus");
  TH1F* h_wh160_METUCPlus 		= ((TH1F*)wh160->Get("METUCPlus/"+histPath) )->Clone("h_wh_METUCPlus");
  TH1F* h_wh160_METUCMinus 		= ((TH1F*)wh160->Get("METUCMinus/"+histPath) )->Clone("h_wh_METUCMinus");
  TH1F* h_wh160_bTagPlus 		= ((TH1F*)wh160->Get("bTagPlus/"+histPath) )->Clone("h_wh_bTagPlus");
  TH1F* h_wh160_bTagMinus 		= ((TH1F*)wh160->Get("bTagMinus/"+histPath) )->Clone("h_wh_bTagMinus");
  
  

  //get average of top pt-reweighting 
  TH1F* hsig_avgTop80 = (TH1F*)wh80->Get("base/SF_topPtWeights"); 
  double sig_avgTop80 = hsig_avgTop80->GetMean();
  double sf80 = 1.0209;
  h_wh80_base->Scale(sf80/sig_avgTop80);
  h_wh80_JESPlus->Scale(sf80/sig_avgTop80);
  h_wh80_JESMinus->Scale(sf80/sig_avgTop80);
  h_wh80_JERPlus->Scale(sf80/sig_avgTop80);
  h_wh80_JERMinus->Scale(sf80/sig_avgTop80);
  h_wh80_METUCPlus->Scale(sf80/sig_avgTop80);
  h_wh80_METUCMinus->Scale(sf80/sig_avgTop80);
  h_wh80_bTagPlus->Scale(sf80/sig_avgTop80);
  h_wh80_bTagMinus->Scale(sf80/sig_avgTop80);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop90 = (TH1F*)wh90->Get("base/SF_topPtWeights"); 
  double sig_avgTop90 = hsig_avgTop90->GetMean();
  double sf90 = 1.0871;
  h_wh90_base->Scale(sf90/sig_avgTop90);
  h_wh90_JESPlus->Scale(sf90/sig_avgTop90);
  h_wh90_JESMinus->Scale(sf90/sig_avgTop90);
  h_wh90_JERPlus->Scale(sf90/sig_avgTop90);
  h_wh90_JERMinus->Scale(sf90/sig_avgTop90);
  h_wh90_METUCPlus->Scale(sf90/sig_avgTop90);
  h_wh90_METUCMinus->Scale(sf90/sig_avgTop90);
  h_wh90_bTagPlus->Scale(sf90/sig_avgTop90);
  h_wh90_bTagMinus->Scale(sf90/sig_avgTop90);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop100 = (TH1F*)wh100->Get("base/SF_topPtWeights"); 
  double sig_avgTop100 = hsig_avgTop100->GetMean();
  double sf100 = 1;
  h_wh100_base->Scale(sf100/sig_avgTop100);
  h_wh100_JESPlus->Scale(sf100/sig_avgTop100);
  h_wh100_JESMinus->Scale(sf100/sig_avgTop100);
  h_wh100_JERPlus->Scale(sf100/sig_avgTop100);
  h_wh100_JERMinus->Scale(sf100/sig_avgTop100);
  h_wh100_METUCPlus->Scale(sf100/sig_avgTop100);
  h_wh100_METUCMinus->Scale(sf100/sig_avgTop100);
  h_wh100_bTagPlus->Scale(sf100/sig_avgTop100);
  h_wh100_bTagMinus->Scale(sf100/sig_avgTop100);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop120 = (TH1F*)wh120->Get("base/SF_topPtWeights"); 
  double sig_avgTop120 = hsig_avgTop120->GetMean();
  double sf120 = 1;
  h_wh120_base->Scale(sf120/sig_avgTop120);
  h_wh120_JESPlus->Scale(sf120/sig_avgTop120);
  h_wh120_JESMinus->Scale(sf120/sig_avgTop120);
  h_wh120_JERPlus->Scale(sf120/sig_avgTop120);
  h_wh120_JERMinus->Scale(sf120/sig_avgTop120);
  h_wh120_METUCPlus->Scale(sf120/sig_avgTop120);
  h_wh120_METUCMinus->Scale(sf120/sig_avgTop120);
  h_wh120_bTagPlus->Scale(sf120/sig_avgTop120);
  h_wh120_bTagMinus->Scale(sf120/sig_avgTop120);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop140 = (TH1F*)wh140->Get("base/SF_topPtWeights"); 
  double sig_avgTop140 = hsig_avgTop140->GetMean();
  double sf140 = 1;
  h_wh140_base->Scale(sf140/sig_avgTop140);
  h_wh140_JESPlus->Scale(sf140/sig_avgTop140);
  h_wh140_JESMinus->Scale(sf140/sig_avgTop140);
  h_wh140_JERPlus->Scale(sf140/sig_avgTop140);
  h_wh140_JERMinus->Scale(sf140/sig_avgTop140);
  h_wh140_METUCPlus->Scale(sf140/sig_avgTop140);
  h_wh140_METUCMinus->Scale(sf140/sig_avgTop140);
  h_wh140_bTagPlus->Scale(sf140/sig_avgTop140);
  h_wh140_bTagMinus->Scale(sf140/sig_avgTop140);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop150 = (TH1F*)wh150->Get("base/SF_topPtWeights"); 
  double sig_avgTop150 = hsig_avgTop150->GetMean();
  double sf150 = 1;
  h_wh150_base->Scale(sf150/sig_avgTop150);
  h_wh150_JESPlus->Scale(sf150/sig_avgTop150);
  h_wh150_JESMinus->Scale(sf150/sig_avgTop150);
  h_wh150_JERPlus->Scale(sf150/sig_avgTop150);
  h_wh150_JERMinus->Scale(sf150/sig_avgTop150);
  h_wh150_METUCPlus->Scale(sf150/sig_avgTop150);
  h_wh150_METUCMinus->Scale(sf150/sig_avgTop150);
  h_wh150_bTagPlus->Scale(sf150/sig_avgTop150);
  h_wh150_bTagMinus->Scale(sf150/sig_avgTop150);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop155 = (TH1F*)wh155->Get("base/SF_topPtWeights"); 
  double sig_avgTop155 = hsig_avgTop155->GetMean();
  double sf155 = 1;
  h_wh155_base->Scale(sf155/sig_avgTop155);
  h_wh155_JESPlus->Scale(sf155/sig_avgTop155);
  h_wh155_JESMinus->Scale(sf155/sig_avgTop155);
  h_wh155_JERPlus->Scale(sf155/sig_avgTop155);
  h_wh155_JERMinus->Scale(sf155/sig_avgTop155);
  h_wh155_METUCPlus->Scale(sf155/sig_avgTop155);
  h_wh155_METUCMinus->Scale(sf155/sig_avgTop155);
  h_wh155_bTagPlus->Scale(sf155/sig_avgTop155);
  h_wh155_bTagMinus->Scale(sf155/sig_avgTop155);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop160 = (TH1F*)wh160->Get("base/SF_topPtWeights"); 
  double sig_avgTop160 = hsig_avgTop160->GetMean();
  double sf160 = 1;
  h_wh160_base->Scale(sf160/sig_avgTop160);
  h_wh160_JESPlus->Scale(sf160/sig_avgTop160);
  h_wh160_JESMinus->Scale(sf160/sig_avgTop160);
  h_wh160_JERPlus->Scale(sf160/sig_avgTop160);
  h_wh160_JERMinus->Scale(sf160/sig_avgTop160);
  h_wh160_METUCPlus->Scale(sf160/sig_avgTop160);
  h_wh160_METUCMinus->Scale(sf160/sig_avgTop160);
  h_wh160_bTagPlus->Scale(sf160/sig_avgTop160);
  h_wh160_bTagMinus->Scale(sf160/sig_avgTop160);

  ofstream outFile; 
  outFile.open("signalTable.tex"); 
   
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
  
  //Add another table with JESUP
  double Nbin = 13;
  getSignalCutFlowTable(outFile, "base", Nbin, 		h_wh80_base,  		h_wh90_base,  		h_wh100_base,  	  	h_wh120_base,  	   h_wh140_base,      h_wh150_base,  	 h_wh155_base,       h_wh160_base      ); 
  getSignalCutFlowTable(outFile, "JESPlus", Nbin, 	h_wh80_JESPlus, 	h_wh90_JESPlus,         h_wh100_JESPlus,        h_wh120_JESPlus,   h_wh140_JESPlus,   h_wh150_JESPlus,   h_wh155_JESPlus,    h_wh160_JESPlus   ); 
  getSignalCutFlowTable(outFile, "JESMinus", Nbin, 	h_wh80_JESMinus, 	h_wh90_JESMinus,        h_wh100_JESMinus,       h_wh120_JESMinus,  h_wh140_JESMinus,  h_wh150_JESMinus,  h_wh155_JESMinus,   h_wh160_JESMinus  );
  getSignalCutFlowTable(outFile, "JERPlus", Nbin, 	h_wh80_JERPlus, 	h_wh90_JERPlus,         h_wh100_JERPlus,        h_wh120_JERPlus,   h_wh140_JERPlus,   h_wh150_JERPlus,   h_wh155_JERPlus,    h_wh160_JERPlus   );
  getSignalCutFlowTable(outFile, "JERMinus", Nbin, 	h_wh80_JERMinus, 	h_wh90_JERMinus,        h_wh100_JERMinus,       h_wh120_JERMinus,  h_wh140_JERMinus,  h_wh150_JERMinus,  h_wh155_JERMinus,   h_wh160_JERMinus  );
  getSignalCutFlowTable(outFile, "METUCPlus", Nbin, 	h_wh80_METUCPlus, 	h_wh90_METUCPlus,       h_wh100_METUCPlus,      h_wh120_METUCPlus, h_wh140_METUCPlus, h_wh150_METUCPlus, h_wh155_METUCPlus,  h_wh160_METUCPlus );
  getSignalCutFlowTable(outFile, "METUCMinus", Nbin, 	h_wh80_METUCMinus, 	h_wh90_METUCMinus,      h_wh100_METUCMinus,     h_wh120_METUCMinus,h_wh140_METUCMinus,h_wh150_METUCMinus,h_wh155_METUCMinus, h_wh160_METUCMinus);
  getSignalCutFlowTable(outFile, "bTagPlus", Nbin, 	h_wh80_bTagPlus, 	h_wh90_bTagPlus,        h_wh100_bTagPlus,       h_wh120_bTagPlus,  h_wh140_bTagPlus,  h_wh150_bTagPlus,  h_wh155_bTagPlus,   h_wh160_bTagPlus  );
  getSignalCutFlowTable(outFile, "bTagMinus", Nbin, 	h_wh80_bTagMinus, 	h_wh90_bTagMinus,       h_wh100_bTagMinus,      h_wh120_bTagMinus, h_wh140_bTagMinus, h_wh150_bTagMinus, h_wh155_bTagMinus,  h_wh160_bTagMinus );
  outFile<<"\\end{document}"<<endl;  
  outFile.close(); 
} 

void makeSummaryTableWithSys(TString histSubDir="Iso/", TString histName="cutflow"){  
  
  TString inFile("$PWD/");
  TString histPath(histSubDir+histName);
  TFile *wh 	  		= new TFile(inFile+"all_Hplus120.root"); 
  TFile *ttbar    		= new TFile(inFile+"all_TTJetsP.root"); 
  TFile *wjet  			= new TFile(inFile+"all_WJets.root"); 
  TFile *zjet  			= new TFile(inFile+"all_DY.root");
  TFile *qcd  			= new TFile(inFile+"all_QCD.root");
  TFile *stop  			= new TFile(inFile+"all_ST.root");
  TFile *diboson 		= new TFile(inFile+"all_VV.root");
  
  //Hplus  M120   signal 
  TH1F* h_wh_base  		= ((TH1F*)wh->Get("base/"+histPath))->Clone("h_wh_base"); 
  TH1F* h_wh_JESPlus 		= ((TH1F*)wh->Get("JESPlus/"+histPath) )->Clone("h_wh_JESPlus");
  TH1F* h_wh_JESMinus 		= ((TH1F*)wh->Get("JESMinus/"+histPath) )->Clone("h_wh_JESMinus");
  TH1F* h_wh_JERPlus 		= ((TH1F*)wh->Get("JERPlus/"+histPath) )->Clone("h_wh_JERPlus");
  TH1F* h_wh_JERMinus 		= ((TH1F*)wh->Get("JERMinus/"+histPath) )->Clone("h_wh_JERMinus");
  TH1F* h_wh_METUCPlus 		= ((TH1F*)wh->Get("METUCPlus/"+histPath) )->Clone("h_wh_METUCPlus");
  TH1F* h_wh_METUCMinus 	= ((TH1F*)wh->Get("METUCMinus/"+histPath) )->Clone("h_wh_METUCMinus");
  TH1F* h_wh_bTagPlus 		= ((TH1F*)wh->Get("bTagPlus/"+histPath) )->Clone("h_wh_bTagPlus");
  TH1F* h_wh_bTagMinus 		= ((TH1F*)wh->Get("bTagMinus/"+histPath) )->Clone("h_wh_bTagMinus");
  
  //ttbar+jets
  TH1F* h_ttbar_base  		= ((TH1F*)ttbar->Get("base/"+histPath) )->Clone("h_ttbar_base");  
  TH1F* h_ttbar_JESPlus 	= ((TH1F*)ttbar->Get("JESPlus/"+histPath) )->Clone("h_ttbar_JESPlus"); 
  TH1F* h_ttbar_JESMinus 	= ((TH1F*)ttbar->Get("JESMinus/"+histPath) )->Clone("h_ttbar_JESMinus"); 
  TH1F* h_ttbar_JERPlus 	= ((TH1F*)ttbar->Get("JERPlus/"+histPath) )->Clone("h_ttbar_JERPlus"); 
  TH1F* h_ttbar_JERMinus 	= ((TH1F*)ttbar->Get("JERMinus/"+histPath) )->Clone("h_ttbar_JERMinus"); 
  TH1F* h_ttbar_METUCPlus 	= ((TH1F*)ttbar->Get("METUCPlus/"+histPath) )->Clone("h_ttbar_METUCPlus"); 
  TH1F* h_ttbar_METUCMinus 	= ((TH1F*)ttbar->Get("METUCMinus/"+histPath) )->Clone("h_ttbar_METUCMinus"); 
  TH1F* h_ttbar_bTagPlus 	= ((TH1F*)ttbar->Get("bTagPlus/"+histPath) )->Clone("h_ttbar_bTagPlus"); 
  TH1F* h_ttbar_bTagMinus 	= ((TH1F*)ttbar->Get("bTagMinus/"+histPath) )->Clone("h_ttbar_bTagMinus"); 
  //w+jets
  TH1F* h_wjet_base  		= ((TH1F*)wjet->Get("base/"+histPath) )->Clone("h_wjet_base");  
  TH1F* h_wjet_JESPlus 		= ((TH1F*)wjet->Get("JESPlus/"+histPath) )->Clone("h_wjet_JESPlus"); 
  TH1F* h_wjet_JESMinus 	= ((TH1F*)wjet->Get("JESMinus/"+histPath) )->Clone("h_wjet_JESMinus"); 
  TH1F* h_wjet_JERPlus 		= ((TH1F*)wjet->Get("JERPlus/"+histPath) )->Clone("h_wjet_JERPlus"); 
  TH1F* h_wjet_JERMinus 	= ((TH1F*)wjet->Get("JERMinus/"+histPath) )->Clone("h_wjet_JERMinus"); 
  TH1F* h_wjet_METUCPlus 	= ((TH1F*)wjet->Get("METUCPlus/"+histPath) )->Clone("h_wjet_METUCPlus"); 
  TH1F* h_wjet_METUCMinus 	= ((TH1F*)wjet->Get("METUCMinus/"+histPath) )->Clone("h_wjet_METUCMinus"); 
  TH1F* h_wjet_bTagPlus 	= ((TH1F*)wjet->Get("bTagPlus/"+histPath) )->Clone("h_wjet_bTagPlus"); 
  TH1F* h_wjet_bTagMinus 	= ((TH1F*)wjet->Get("bTagMinus/"+histPath) )->Clone("h_wjet_bTagMinus"); 
  //dy+jets
  TH1F* h_zjet_base  		= ((TH1F*)zjet->Get("base/"+histPath) )->Clone("h_zjet_base");
  TH1F* h_zjet_JESPlus 		= ((TH1F*)zjet->Get("JESPlus/"+histPath) )->Clone("h_zjet_JESPlus");
  TH1F* h_zjet_JESMinus 	= ((TH1F*)zjet->Get("JESMinus/"+histPath) )->Clone("h_zjet_JESMinus");
  TH1F* h_zjet_JERPlus 		= ((TH1F*)zjet->Get("JERPlus/"+histPath) )->Clone("h_zjet_JERPlus");
  TH1F* h_zjet_JERMinus 	= ((TH1F*)zjet->Get("JERMinus/"+histPath) )->Clone("h_zjet_JERMinus");
  TH1F* h_zjet_METUCPlus 	= ((TH1F*)zjet->Get("METUCPlus/"+histPath) )->Clone("h_zjet_METUCPlus");
  TH1F* h_zjet_METUCMinus 	= ((TH1F*)zjet->Get("METUCMinus/"+histPath) )->Clone("h_zjet_METUCMinus");
  TH1F* h_zjet_bTagPlus 	= ((TH1F*)zjet->Get("bTagPlus/"+histPath) )->Clone("h_zjet_bTagPlus");
  TH1F* h_zjet_bTagMinus 	= ((TH1F*)zjet->Get("bTagMinus/"+histPath) )->Clone("h_zjet_bTagMinus");
  //qcd
  TH1F* h_qcd_base 		= ((TH1F*)qcd->Get("base/"+histPath) )->Clone("h_qcd_base");
  TH1F* h_qcd_JESPlus 		= ((TH1F*)qcd->Get("JESPlus/"+histPath) )->Clone("h_qcd_JESPlus");
  TH1F* h_qcd_JESMinus 		= ((TH1F*)qcd->Get("JESMinus/"+histPath) )->Clone("h_qcd_JESMinus");
  TH1F* h_qcd_JERPlus 		= ((TH1F*)qcd->Get("JERPlus/"+histPath) )->Clone("h_qcd_JERPlus");
  TH1F* h_qcd_JERMinus 		= ((TH1F*)qcd->Get("JERMinus/"+histPath) )->Clone("h_qcd_JERMinus");
  TH1F* h_qcd_METUCPlus 	= ((TH1F*)qcd->Get("METUCPlus/"+histPath) )->Clone("h_qcd_METUCPlus");
  TH1F* h_qcd_METUCMinus 	= ((TH1F*)qcd->Get("METUCMinus/"+histPath) )->Clone("h_qcd_METUCMinus");
  TH1F* h_qcd_bTagPlus 		= ((TH1F*)qcd->Get("bTagPlus/"+histPath) )->Clone("h_qcd_bTagPlus");
  TH1F* h_qcd_bTagMinus 	= ((TH1F*)qcd->Get("bTagMinus/"+histPath) )->Clone("h_qcd_bTagMinus");
  //single top
  TH1F* h_stop_base  		= ((TH1F*)stop->Get("base/"+histPath) )->Clone("h_stop_base");
  TH1F* h_stop_JESPlus 		= ((TH1F*)stop->Get("JESPlus/"+histPath) )->Clone("h_stop_JESPlus");
  TH1F* h_stop_JESMinus 	= ((TH1F*)stop->Get("JESMinus/"+histPath) )->Clone("h_stop_JESMinus");
  TH1F* h_stop_JERPlus 		= ((TH1F*)stop->Get("JERPlus/"+histPath) )->Clone("h_stop_JERPlus");
  TH1F* h_stop_JERMinus 	= ((TH1F*)stop->Get("JERMinus/"+histPath) )->Clone("h_stop_JERMinus");
  TH1F* h_stop_METUCPlus 	= ((TH1F*)stop->Get("METUCPlus/"+histPath) )->Clone("h_stop_METUCPlus");
  TH1F* h_stop_METUCMinus 	= ((TH1F*)stop->Get("METUCMinus/"+histPath) )->Clone("h_stop_METUCMinus");
  TH1F* h_stop_bTagPlus 	= ((TH1F*)stop->Get("bTagPlus/"+histPath) )->Clone("h_stop_bTagPlus");
  TH1F* h_stop_bTagMinus 	= ((TH1F*)stop->Get("bTagMinus/"+histPath) )->Clone("h_stop_bTagMinus");
  //vv
  TH1F* h_diboson_base 		= ((TH1F*)diboson->Get("base/"+histPath) )->Clone("h_diboson_base");
  TH1F* h_diboson_JESPlus 	= ((TH1F*)diboson->Get("JESPlus/"+histPath) )->Clone("h_diboson_JESPlus");
  TH1F* h_diboson_JESMinus 	= ((TH1F*)diboson->Get("JESMinus/"+histPath) )->Clone("h_diboson_JESMinus");
  TH1F* h_diboson_JERPlus 	= ((TH1F*)diboson->Get("JERPlus/"+histPath) )->Clone("h_diboson_JERPlus");
  TH1F* h_diboson_JERMinus 	= ((TH1F*)diboson->Get("JERMinus/"+histPath) )->Clone("h_diboson_JERMinus");
  TH1F* h_diboson_METUCPlus 	= ((TH1F*)diboson->Get("METUCPlus/"+histPath) )->Clone("h_diboson_METUCPlus");
  TH1F* h_diboson_METUCMinus 	= ((TH1F*)diboson->Get("METUCMinus/"+histPath) )->Clone("h_diboson_METUCMinus");
  TH1F* h_diboson_bTagPlus 	= ((TH1F*)diboson->Get("bTagPlus/"+histPath) )->Clone("h_diboson_bTagPlus");
  TH1F* h_diboson_bTagMinus 	= ((TH1F*)diboson->Get("bTagMinus/"+histPath) )->Clone("h_diboson_bTagMinus");

  //get average of top pt-reweighting 
  TH1F* hsig_avgTop120 = (TH1F*)wh->Get("base/SF_topPtWeights"); 
  double sig_avgTop120 = hsig_avgTop120->GetMean();
  double sf120 = 1;
  h_wh_base->Scale(sf120/sig_avgTop120);
  h_wh_JESPlus->Scale(sf120/sig_avgTop120);
  h_wh_JESMinus->Scale(sf120/sig_avgTop120);
  h_wh_JERPlus->Scale(sf120/sig_avgTop120);
  h_wh_JERMinus->Scale(sf120/sig_avgTop120);
  h_wh_METUCPlus->Scale(sf120/sig_avgTop120);
  h_wh_METUCMinus->Scale(sf120/sig_avgTop120);
  h_wh_bTagPlus->Scale(sf120/sig_avgTop120);
  h_wh_bTagMinus->Scale(sf120/sig_avgTop120);

  
  TH1F* httbar_avgTop = (TH1F*)ttbar->Get("base/SF_topPtWeights"); 
  double ttbar_avgTop = httbar_avgTop->GetMean();
  double sf_top = 1.05684;
  h_ttbar_base->Scale(sf_top/ttbar_avgTop);
  h_ttbar_JESPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_JESMinus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_JERPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_JERMinus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_METUCPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_METUCMinus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_bTagPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_bTagMinus->Scale(sf_top/ttbar_avgTop);

  TH1F* h_TotalBkg = h_wjet_base->Clone("h_TotalBkg");
  h_TotalBkg->Reset();
  h_TotalBkg->Add(h_ttbar_base);
  h_TotalBkg->Add(h_wjet_base);
  h_TotalBkg->Add(h_zjet_base);
  h_TotalBkg->Add(h_qcd_base);
  h_TotalBkg->Add(h_stop_base);
  h_TotalBkg->Add(h_diboson_base);


  TFile *data = new TFile(inFile+"all_muData.root");
  TH1F* h_data = ((TH1F*)data->Get("base/"+histPath))->Clone("h_data");


  ofstream outFile; 
  outFile.open("summaryWithSyst.tex"); 
   
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
  //outFile<<"\\begin{LARGE}"<<endl;  
  outFile<<"\\begin{tabular}{ | c| c| }"<<endl;  
  outFile<<"\\multicolumn{2}{c}{ } \\\\"<<endl;  
  outFile<<"\\hline "<<endl; 
  outFile<<"\\multicolumn{1}{|c|}{Source} & \\multicolumn{1}{|c|}{N$_{\rm events}$ $\\pm$ MC stat $\\pm$ JES/MET scale $\\pm$ bTag } \\\\"<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<"\\hline "<<endl;
  int b = 11;
  outFile<<"HW, $M_{H}=120~GeV/c^{2}$"<<" & "<<h_wh_base->GetBinContent(b)<<" $\\pm$ "<<h_wh_base->GetBinError(b)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_wh_JESPlus->GetBinContent(b) - h_wh_base->GetBinContent(b)), fabs(h_wh_base->GetBinContent(b) - h_wh_JESMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_wh_JERPlus->GetBinContent(b) - h_wh_base->GetBinContent(b)), fabs(h_wh_base->GetBinContent(b) - h_wh_JERMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_wh_METUCPlus->GetBinContent(b) - h_wh_base->GetBinContent(b)), fabs(h_wh_base->GetBinContent(b) - h_wh_METUCMinus->GetBinContent(b))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_wh_bTagPlus->GetBinContent(b) - h_wh_base->GetBinContent(b)), fabs(h_wh_base->GetBinContent(b) - h_wh_bTagMinus->GetBinContent(b))) <<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;  
  outFile<<"SM $t\\bar{t}$"<<" & "<<h_ttbar_base->GetBinContent(b)<<" $\\pm$ "<<h_ttbar_base->GetBinError(b)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_ttbar_JESPlus->GetBinContent(b) - h_ttbar_base->GetBinContent(b)), fabs(h_ttbar_base->GetBinContent(b) - h_ttbar_JESMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_ttbar_JERPlus->GetBinContent(b) - h_ttbar_base->GetBinContent(b)), fabs(h_ttbar_base->GetBinContent(b) - h_ttbar_JERMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_ttbar_METUCPlus->GetBinContent(b) - h_ttbar_base->GetBinContent(b)), fabs(h_ttbar_base->GetBinContent(b) - h_ttbar_METUCMinus->GetBinContent(b))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_ttbar_bTagPlus->GetBinContent(b) - h_ttbar_base->GetBinContent(b)), fabs(h_ttbar_base->GetBinContent(b) - h_ttbar_bTagMinus->GetBinContent(b))) <<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;

  outFile<<"W+Jets"<<" & "<<h_wjet_base->GetBinContent(b)<<" $\\pm$ "<<h_wjet_base->GetBinError(b)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_wjet_JESPlus->GetBinContent(b) - h_wjet_base->GetBinContent(b)), fabs(h_wjet_base->GetBinContent(b) - h_wjet_JESMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_wjet_JERPlus->GetBinContent(b) - h_wjet_base->GetBinContent(b)), fabs(h_wjet_base->GetBinContent(b) - h_wjet_JERMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_wjet_METUCPlus->GetBinContent(b) - h_wjet_base->GetBinContent(b)), fabs(h_wjet_base->GetBinContent(b) - h_wjet_METUCMinus->GetBinContent(b))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_wjet_bTagPlus->GetBinContent(b) - h_wjet_base->GetBinContent(b)), fabs(h_wjet_base->GetBinContent(b) - h_wjet_bTagMinus->GetBinContent(b))) <<" \\\\ "<<endl;  
  outFile<<"\\hline "<<endl;  
  
  outFile<<"Z+Jets"<<" & "<<h_zjet_base->GetBinContent(b)<<" $\\pm$ "<<h_zjet_base->GetBinError(b)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_zjet_JESPlus->GetBinContent(b) - h_zjet_base->GetBinContent(b)), fabs(h_zjet_base->GetBinContent(b) - h_zjet_JESMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_zjet_JERPlus->GetBinContent(b) - h_zjet_base->GetBinContent(b)), fabs(h_zjet_base->GetBinContent(b) - h_zjet_JERMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_zjet_METUCPlus->GetBinContent(b) - h_zjet_base->GetBinContent(b)), fabs(h_zjet_base->GetBinContent(b) - h_zjet_METUCMinus->GetBinContent(b))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_zjet_bTagPlus->GetBinContent(b) - h_zjet_base->GetBinContent(b)), fabs(h_zjet_base->GetBinContent(b) - h_zjet_bTagMinus->GetBinContent(b))) <<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;

  outFile<<"QCD"<<" & "<<h_qcd_base->GetBinContent(b)<<" $\\pm$ "<<h_qcd_base->GetBinError(b)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_qcd_JESPlus->GetBinContent(b) - h_qcd_base->GetBinContent(b)), fabs(h_qcd_base->GetBinContent(b) - h_qcd_JESMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_qcd_JERPlus->GetBinContent(b) - h_qcd_base->GetBinContent(b)), fabs(h_qcd_base->GetBinContent(b) - h_qcd_JERMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_qcd_METUCPlus->GetBinContent(b) - h_qcd_base->GetBinContent(b)), fabs(h_qcd_base->GetBinContent(b) - h_qcd_METUCMinus->GetBinContent(b))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_qcd_bTagPlus->GetBinContent(b) - h_qcd_base->GetBinContent(b)), fabs(h_qcd_base->GetBinContent(b) - h_qcd_bTagMinus->GetBinContent(b))) <<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;

  outFile<<"SingleTop"<<" & "<<h_stop_base->GetBinContent(b)<<" $\\pm$ "<<h_stop_base->GetBinError(b)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_stop_JESPlus->GetBinContent(b) - h_stop_base->GetBinContent(b)), fabs(h_stop_base->GetBinContent(b) - h_stop_JESMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_stop_JERPlus->GetBinContent(b) - h_stop_base->GetBinContent(b)), fabs(h_stop_base->GetBinContent(b) - h_stop_JERMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_stop_METUCPlus->GetBinContent(b) - h_stop_base->GetBinContent(b)), fabs(h_stop_base->GetBinContent(b) - h_stop_METUCMinus->GetBinContent(b))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_stop_bTagPlus->GetBinContent(b) - h_stop_base->GetBinContent(b)), fabs(h_stop_base->GetBinContent(b) - h_stop_bTagMinus->GetBinContent(b))) <<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;

  outFile<<"Dibosons"<<" & "<<h_diboson_base->GetBinContent(b)<<" $\\pm$ "<<h_diboson_base->GetBinError(b)<<" $\\pm$ "<< sqrt(pow(TMath::Max(fabs(h_diboson_JESPlus->GetBinContent(b) - h_diboson_base->GetBinContent(b)), fabs(h_diboson_base->GetBinContent(b) - h_diboson_JESMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_diboson_JERPlus->GetBinContent(b) - h_diboson_base->GetBinContent(b)), fabs(h_diboson_base->GetBinContent(b) - h_diboson_JERMinus->GetBinContent(b))), 2) + pow(TMath::Max(fabs(h_diboson_METUCPlus->GetBinContent(b) - h_diboson_base->GetBinContent(b)), fabs(h_diboson_base->GetBinContent(b) - h_diboson_METUCMinus->GetBinContent(b))), 2)) <<" $\\pm$ "<<TMath::Max(fabs(h_diboson_bTagPlus->GetBinContent(b) - h_diboson_base->GetBinContent(b)), fabs(h_diboson_base->GetBinContent(b) - h_diboson_bTagMinus->GetBinContent(b))) <<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"Total Bkg"<<" & "<<h_TotalBkg->GetBinContent(b)<<" $\\\pm$ "<<h_TotalBkg->GetBinError(b)<<" $\\\pm$ "<<" -- "<<" $\\\pm$ "<<" -- "<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;   
  outFile<<"\\hline "<<endl;  
  outFile<<" Data "<<" & "<<h_data->GetBinContent(b)<<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;    
  outFile<<"\\hline "<<endl;   
  outFile<<"\\end{tabular}"<<endl; 
  //outFile<<"\\end{LARGE}"<<endl;  
  outFile<<"\\end{center}"<<endl;
  outFile<<"\\\caption{First table}"<<endl;
  outFile<<"\\end{table}"<<endl;

  //Add another table with JESUP
  double Nbin = 13;
  makeCutFlowTableAny(outFile, "base", Nbin, 		h_wh_base,  	h_ttbar_base, h_stop_base, h_wjet_base, h_zjet_base, h_qcd_base, h_stop_base, h_diboson_base, h_data);
  makeCutFlowTableAny(outFile, "JESPlus", Nbin, 	h_wh_JESPlus, 	h_ttbar_JESPlus, h_stop_JESPlus, h_wjet_JESPlus, h_zjet_JESPlus, h_qcd_JESPlus, h_stop_JESPlus, h_diboson_JESPlus, h_data);
  makeCutFlowTableAny(outFile, "JESMinus", Nbin, 	h_wh_JESMinus, 	h_ttbar_JESMinus,h_stop_JERMinus, h_wjet_JESMinus, h_zjet_JESMinus, h_qcd_JESMinus, h_stop_JESMinus, h_diboson_JESMinus, h_data);
  makeCutFlowTableAny(outFile, "JERPlus", Nbin, 	h_wh_JERPlus, 	h_ttbar_JERPlus, h_stop_JERPlus, h_wjet_JERPlus, h_zjet_JERPlus, h_qcd_JERPlus, h_stop_JERPlus, h_diboson_JERPlus, h_data);
  makeCutFlowTableAny(outFile, "JERMinus", Nbin, 	h_wh_JERMinus, 	h_ttbar_JERMinus, h_stop_JERMinus, h_wjet_JERMinus, h_zjet_JERMinus, h_qcd_JERMinus, h_stop_JERMinus, h_diboson_JERMinus, h_data);
  makeCutFlowTableAny(outFile, "METUCPlus", Nbin, 	h_wh_METUCPlus, h_ttbar_METUCPlus, h_stop_METUCPlus, h_wjet_METUCPlus, h_zjet_METUCPlus, h_qcd_METUCPlus, h_stop_METUCPlus, h_diboson_METUCPlus, h_data);
  makeCutFlowTableAny(outFile, "METUCMinus", Nbin, 	h_wh_METUCMinus, h_ttbar_METUCMinus, h_stop_METUCMinus, h_wjet_METUCMinus, h_zjet_METUCMinus, h_qcd_METUCMinus, h_stop_METUCMinus, h_diboson_METUCMinus, h_data);
  makeCutFlowTableAny(outFile, "bTagPlus", Nbin, 	h_wh_bTagPlus, 	h_ttbar_bTagPlus, h_stop_bTagPlus, h_wjet_bTagPlus, h_zjet_bTagPlus, h_qcd_bTagPlus, h_stop_bTagPlus, h_diboson_bTagPlus, h_data);
  makeCutFlowTableAny(outFile, "bTagMinus", Nbin, 	h_wh_bTagMinus, h_ttbar_bTagMinus, h_stop_bTagMinus, h_wjet_bTagMinus, h_zjet_bTagMinus, h_qcd_bTagMinus, h_stop_bTagMinus, h_diboson_bTagMinus, h_data);
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
  
  TFile *data = new TFile(inFile+"all_muData.root");
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

  TFile *data = new TFile(inFile+"all_muData.root");
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

