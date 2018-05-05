#include <iostream>
#include <fstream>
#include <iomanip>

bool isMuChannel = true;
bool isEleChannel = false;

double sysUncJESTopPt( TH1F * h_JESPlus, TH1F * h_base, TH1F * h_JESMinus, TH1F * h_JERPlus, TH1F * h_JERMinus, TH1F * h_TopPtPlus, TH1F * h_TopPtMinus, int b){
  double uncJES = pow(TMath::Max(fabs(h_JESPlus->GetBinContent(b) - h_base->GetBinContent(b)), fabs(h_base->GetBinContent(b) - h_JESMinus->GetBinContent(b))), 2);
  double uncJER = pow(TMath::Max(fabs(h_JERPlus->GetBinContent(b) - h_base->GetBinContent(b)), fabs(h_base->GetBinContent(b) - h_JERMinus->GetBinContent(b))), 2);
  double uncTop = pow(TMath::Max(fabs(h_TopPtPlus->GetBinContent(b) - h_base->GetBinContent(b)), fabs(h_base->GetBinContent(b) - h_TopPtMinus->GetBinContent(b))), 2);
  double unc = sqrt(uncJES +uncJER  +uncTop);
  return unc ;
}

double sysUncBCTag (TH1F * h_bTagPlus, TH1F * h_base, TH1F * h_bTagMinus, int b){
	double uncTag = TMath::Max(fabs(h_bTagPlus->GetBinContent(b) - h_base->GetBinContent(b)), fabs(h_base->GetBinContent(b) - h_bTagMinus->GetBinContent(b)));
	return uncTag;
}	

//Relative sys unc
double relSysUncJetMET( TH1F * h_JESPlus, TH1F * h_base, TH1F * h_JESMinus, TH1F * h_JERPlus, TH1F * h_JERMinus, int b){
  double uncJES = pow(TMath::Max(fabs(h_JESPlus->GetBinContent(b) - h_base->GetBinContent(b)), fabs(h_base->GetBinContent(b) - h_JESMinus->GetBinContent(b))), 2);
  double uncJER = pow(TMath::Max(fabs(h_JERPlus->GetBinContent(b) - h_base->GetBinContent(b)), fabs(h_base->GetBinContent(b) - h_JERMinus->GetBinContent(b))), 2);
  double unc = 100*(sqrt(uncJES +uncJER)/h_base->GetBinContent(b));
  return unc ;
}

double relSysUncTopPt( TH1F * h_base, TH1F * h_TopPtPlus, TH1F * h_TopPtMinus, int b){
  double uncTop = pow(TMath::Max(fabs(h_TopPtPlus->GetBinContent(b) - h_base->GetBinContent(b)), fabs(h_base->GetBinContent(b) - h_TopPtMinus->GetBinContent(b))), 2);
  double unc = 100*(sqrt(uncTop)/h_base->GetBinContent(b));
  return unc ;
}
double relSysUncBCTag (TH1F * h_bTagPlus, TH1F * h_base, TH1F * h_bTagMinus, int b){
	double uncTag = TMath::Max(fabs(h_bTagPlus->GetBinContent(b) - h_base->GetBinContent(b)), fabs(h_base->GetBinContent(b) - h_bTagMinus->GetBinContent(b)));
	return 100*(uncTag/h_base->GetBinContent(b));
}

double relStatUnc (TH1F * h_base, int b){
	return 100*(h_base->GetBinError(b)/h_base->GetBinContent(b));
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
  outFile<<"\\footnotesize\\setlength{\\tabcolsep}{0.3pt}"<<endl;
  outFile<<"\\begin{tabular}{ | c| c| c| c| c| c| c| c| c| c| c| c| c|c|c|c|}"<<endl; 
  outFile<<"\\multicolumn{5}{c}{ } \\\\"<<endl; 
  outFile<<"\\hline "<<endl;
  if(isMuChannel) outFile<<"\\multicolumn{1}{| c|}{"+ sys +" } & \\multicolumn{1}{ c|}{ $HLT\\_IsoMu24$ } & \\multicolumn{1}{ c|}{ $N_{muon}=1$ } & \\multicolumn{1}{ c|}{ $N_{ele}=0$} & \\multicolumn{1}{ c|}{ Muon SF } & \\multicolumn{1}{ c|}{ $I_{rel}^\\mu < 0.15$ } &\\multicolumn{1}{ c|}{ $N_{jets}\\ge 4$ } & \\multicolumn{1}{ c|}{ $\\not\\!\\!E_T \\ge 20GeV$ }&  \\multicolumn{1}{ c |}{ $\\ge$ 2btag }& \\multicolumn{1}{ c |}{ BTag SF }& \\multicolumn{1}{ c|}{fit converges } & \\multicolumn{1}{ c|}{$Pt_{jets}^{kf}\\ge 25 GeV$ } & \\multicolumn{1}{ c|}{$\\Delta R_{jets}^{pf,kf}\\le 0.2$} & \\multicolumn{1}{ c|}{CTagL} & \\multicolumn{1}{ c|}{CTagLSF} \\\\ "<<endl;
  if(isEleChannel) outFile<<"\\multicolumn{1}{| c|}{"+ sys +" } & \\multicolumn{1}{ c|}{ $HLT\\_Ele27$ } & \\multicolumn{1}{ c|}{ $N_{ele}=1$ } & \\multicolumn{1}{ c|}{ $N_{muon}=0$} & \\multicolumn{1}{ c|}{ Ele SF } & \\multicolumn{1}{ c|}{ $I_{rel}^e < 0.08$ } &\\multicolumn{1}{ c|}{ $N_{jets}\\ge 4$ } & \\multicolumn{1}{ c|}{ $\\not\\!\\!E_T \\ge 20GeV$ }&  \\multicolumn{1}{ c |}{ $\\ge$ 2btag }& \\multicolumn{1}{ c |}{ BTag SF }& \\multicolumn{1}{ c|}{fit converges } & \\multicolumn{1}{ c|}{$Pt_{jets}^{kf}\\ge 25 GeV$ } & \\multicolumn{1}{ c|}{$\\Delta R_{jets}^{pf,kf}\\le 0.2$} & \\multicolumn{1}{ c|}{CTagL} & \\multicolumn{1}{ c|}{CTagLSF} \\\\ "<<endl;
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
  outFile<<"\\caption{Number of evets after various cuts for sys: " +sys+ "}"<<endl; 
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
  outFile<<"\\footnotesize\\setlength{\\tabcolsep}{0.3pt}"<<endl;
  outFile<<"\\begin{tabular}{ | c| c| c| c| c| c| c| c| c| c| c| c| c|c|c|c|}"<<endl; 
  outFile<<"\\multicolumn{5}{c}{ } \\\\"<<endl; 
  outFile<<"\\hline "<<endl;
  if(isMuChannel) outFile<<"\\multicolumn{1}{| c|}{"+ sys +" } & \\multicolumn{1}{ c|}{ $HLT\\_IsoMu24$ } & \\multicolumn{1}{ c|}{ $N_{muon}=1$ } & \\multicolumn{1}{ c|}{ $N_{ele}=0$} & \\multicolumn{1}{ c|}{ Muon SF } & \\multicolumn{1}{ c|}{ $I_{rel}^\\mu < 0.15$ } &\\multicolumn{1}{ c|}{ $N_{jets}\\ge 4$ } & \\multicolumn{1}{ c|}{ $\\not\\!\\!E_T \\ge 20GeV$ }&  \\multicolumn{1}{ c |}{ $\\ge$ 2btag }& \\multicolumn{1}{ c |}{ BTag SF }& \\multicolumn{1}{ c|}{fit converges } & \\multicolumn{1}{ c|}{$Pt_{jets}^{kf}\\ge 25 GeV$ } & \\multicolumn{1}{ c|}{$\\Delta R_{jets}^{pf,kf}\\le 0.2$} & \\multicolumn{1}{ c|}{CTagL} & \\multicolumn{1}{ c|}{CTagLSF} \\\\ "<<endl;
  if(isEleChannel) outFile<<"\\multicolumn{1}{| c|}{"+ sys +" } & \\multicolumn{1}{ c|}{ $HLT\\_Ele27$ } & \\multicolumn{1}{ c|}{ $N_{ele}=1$ } & \\multicolumn{1}{ c|}{ $N_{muon}=0$} & \\multicolumn{1}{ c|}{ Ele SF } & \\multicolumn{1}{ c|}{ $I_{rel}^e < 0.08$ } &\\multicolumn{1}{ c|}{ $N_{jets}\\ge 4$ } & \\multicolumn{1}{ c|}{ $\\not\\!\\!E_T \\ge 20GeV$ }&  \\multicolumn{1}{ c |}{ $\\ge$ 2btag }& \\multicolumn{1}{ c |}{ BTag SF }& \\multicolumn{1}{ c|}{fit converges } & \\multicolumn{1}{ c|}{$Pt_{jets}^{kf}\\ge 25 GeV$ } & \\multicolumn{1}{ c|}{$\\Delta R_{jets}^{pf,kf}\\le 0.2$} & \\multicolumn{1}{ c|}{CTagL} & \\multicolumn{1}{ c|}{CTagLSF} \\\\ "<<endl;
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
  outFile<<"\\caption{Number of evets after various cuts for sys: " +sys+ "}"<<endl; 
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
  
  TH1F* h_wh_base  		= (TH1F*)wh->Get("base/"+histPath)->Clone("h_wh_base"); 
  TH1F* h_ttbar_base  		= (TH1F*)ttbar->Get("base/"+histPath)->Clone("h_ttbar_base");  
  TH1F* h_wjet_base  		= (TH1F*)wjet->Get("base/"+histPath)->Clone("h_wjet_base");  
  TH1F* h_zjet_base  		= (TH1F*)zjet->Get("base/"+histPath)->Clone("h_zjet_base");
  TH1F* h_qcd_base 		= (TH1F*)qcd->Get("base/"+histPath)->Clone("h_qcd_base");
  TH1F* h_stop_base  		= (TH1F*)stop->Get("base/"+histPath)->Clone("h_stop_base");
  TH1F* h_diboson_base 		= (TH1F*)diboson->Get("base/"+histPath)->Clone("h_diboson_base");

  //get average of top pt-reweighting 
  TH1F* hsig_avgTop120 = (TH1F*)wh->Get("base/SF_topPtWeights"); 
  double sig_avgTop120 = hsig_avgTop120->GetMean();
  double sf120 = 1;
  h_wh_base->Scale(sf120/sig_avgTop120);
  
  TH1F* httbar_avgTop = (TH1F*)ttbar->Get("base/SF_topPtWeights"); 
  double ttbar_avgTop = httbar_avgTop->GetMean();
  double sf_ttbar = 1.0;
  h_ttbar_base->Scale(sf_ttbar/ttbar_avgTop);

  TH1F* h_TotalBkg = (TH1F*)h_wjet_base->Clone("h_TotalBkg");
  h_TotalBkg->Reset();
  h_TotalBkg->Add(h_ttbar_base);
  h_TotalBkg->Add(h_wjet_base);
  h_TotalBkg->Add(h_zjet_base);
  h_TotalBkg->Add(h_qcd_base);
  h_TotalBkg->Add(h_stop_base);
  h_TotalBkg->Add(h_diboson_base);

  TFile *data = new TFile(inFile+"all_muData.root");
  TH1F* h_data = (TH1F*)data->Get("base/"+histPath)->Clone("h_data");

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
  double Nbin = 15;
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
  TH1F* h_wh80_base  			= (TH1F*)wh80->Get("base/"+histPath)->Clone("h_wh_base"); 
  TH1F* h_wh80_JESPlus 			= (TH1F*)wh80->Get("JESPlus/"+histPath)->Clone("h_wh_JESPlus");
  TH1F* h_wh80_JESMinus 		= (TH1F*)wh80->Get("JESMinus/"+histPath)->Clone("h_wh_JESMinus");
  TH1F* h_wh80_JERPlus 			= (TH1F*)wh80->Get("JERPlus/"+histPath)->Clone("h_wh_JERPlus");
  TH1F* h_wh80_JERMinus 		= (TH1F*)wh80->Get("JERMinus/"+histPath)->Clone("h_wh_JERMinus");
  TH1F* h_wh80_TopPtPlus 		= (TH1F*)wh80->Get("TopPtPlus/"+histPath)->Clone("h_wh_TopPtPlus");
  TH1F* h_wh80_TopPtMinus 		= (TH1F*)wh80->Get("TopPtMinus/"+histPath)->Clone("h_wh_TopPtMinus");
  TH1F* h_wh80_bTagPlus 		= (TH1F*)wh80->Get("bTagPlus/"+histPath)->Clone("h_wh_bTagPlus");
  TH1F* h_wh80_bTagMinus 		= (TH1F*)wh80->Get("bTagMinus/"+histPath)->Clone("h_wh_bTagMinus");
  TH1F* h_wh80_cTagPlus 		= (TH1F*)wh80->Get("cTagPlus/"+histPath)->Clone("h_wh_cTagPlus");
  TH1F* h_wh80_cTagMinus 		= (TH1F*)wh80->Get("cTagMinus/"+histPath)->Clone("h_wh_cTagMinus");
  
  //Hplus  M90   signal 
  TH1F* h_wh90_base  			= (TH1F*)wh90->Get("base/"+histPath)->Clone("h_wh_base"); 
  TH1F* h_wh90_JESPlus 			= (TH1F*)wh90->Get("JESPlus/"+histPath)->Clone("h_wh_JESPlus");
  TH1F* h_wh90_JESMinus 		= (TH1F*)wh90->Get("JESMinus/"+histPath)->Clone("h_wh_JESMinus");
  TH1F* h_wh90_JERPlus 			= (TH1F*)wh90->Get("JERPlus/"+histPath)->Clone("h_wh_JERPlus");
  TH1F* h_wh90_JERMinus 		= (TH1F*)wh90->Get("JERMinus/"+histPath)->Clone("h_wh_JERMinus");
  TH1F* h_wh90_TopPtPlus 		= (TH1F*)wh90->Get("TopPtPlus/"+histPath)->Clone("h_wh_TopPtPlus");
  TH1F* h_wh90_TopPtMinus 		= (TH1F*)wh90->Get("TopPtMinus/"+histPath)->Clone("h_wh_TopPtMinus");
  TH1F* h_wh90_bTagPlus 		= (TH1F*)wh90->Get("bTagPlus/"+histPath)->Clone("h_wh_bTagPlus");
  TH1F* h_wh90_bTagMinus 		= (TH1F*)wh90->Get("bTagMinus/"+histPath)->Clone("h_wh_bTagMinus");
  TH1F* h_wh90_cTagPlus 		= (TH1F*)wh90->Get("cTagPlus/"+histPath)->Clone("h_wh_cTagPlus");
  TH1F* h_wh90_cTagMinus 		= (TH1F*)wh90->Get("cTagMinus/"+histPath)->Clone("h_wh_cTagMinus");
  
  //Hplus  M100   signal 
  TH1F* h_wh100_base  			= (TH1F*)wh100->Get("base/"+histPath)->Clone("h_wh_base"); 
  TH1F* h_wh100_JESPlus 		= (TH1F*)wh100->Get("JESPlus/"+histPath)->Clone("h_wh_JESPlus");
  TH1F* h_wh100_JESMinus 		= (TH1F*)wh100->Get("JESMinus/"+histPath)->Clone("h_wh_JESMinus");
  TH1F* h_wh100_JERPlus 		= (TH1F*)wh100->Get("JERPlus/"+histPath)->Clone("h_wh_JERPlus");
  TH1F* h_wh100_JERMinus 		= (TH1F*)wh100->Get("JERMinus/"+histPath)->Clone("h_wh_JERMinus");
  TH1F* h_wh100_TopPtPlus 		= (TH1F*)wh100->Get("TopPtPlus/"+histPath)->Clone("h_wh_TopPtPlus");
  TH1F* h_wh100_TopPtMinus 		= (TH1F*)wh100->Get("TopPtMinus/"+histPath)->Clone("h_wh_TopPtMinus");
  TH1F* h_wh100_bTagPlus 		= (TH1F*)wh100->Get("bTagPlus/"+histPath)->Clone("h_wh_bTagPlus");
  TH1F* h_wh100_bTagMinus 		= (TH1F*)wh100->Get("bTagMinus/"+histPath)->Clone("h_wh_bTagMinus");
  TH1F* h_wh100_cTagPlus 		= (TH1F*)wh100->Get("cTagPlus/"+histPath)->Clone("h_wh_cTagPlus");
  TH1F* h_wh100_cTagMinus 		= (TH1F*)wh100->Get("cTagMinus/"+histPath)->Clone("h_wh_cTagMinus");
  
  
  //Hplus  M120   signal 
  TH1F* h_wh120_base  			= (TH1F*)wh120->Get("base/"+histPath)->Clone("h_wh_base"); 
  TH1F* h_wh120_JESPlus 		= (TH1F*)wh120->Get("JESPlus/"+histPath)->Clone("h_wh_JESPlus");
  TH1F* h_wh120_JESMinus 		= (TH1F*)wh120->Get("JESMinus/"+histPath)->Clone("h_wh_JESMinus");
  TH1F* h_wh120_JERPlus 		= (TH1F*)wh120->Get("JERPlus/"+histPath)->Clone("h_wh_JERPlus");
  TH1F* h_wh120_JERMinus 		= (TH1F*)wh120->Get("JERMinus/"+histPath)->Clone("h_wh_JERMinus");
  TH1F* h_wh120_TopPtPlus 		= (TH1F*)wh120->Get("TopPtPlus/"+histPath)->Clone("h_wh_TopPtPlus");
  TH1F* h_wh120_TopPtMinus 		= (TH1F*)wh120->Get("TopPtMinus/"+histPath)->Clone("h_wh_TopPtMinus");
  TH1F* h_wh120_bTagPlus 		= (TH1F*)wh120->Get("bTagPlus/"+histPath)->Clone("h_wh_bTagPlus");
  TH1F* h_wh120_bTagMinus 		= (TH1F*)wh120->Get("bTagMinus/"+histPath)->Clone("h_wh_bTagMinus");
  TH1F* h_wh120_cTagPlus 		= (TH1F*)wh120->Get("cTagPlus/"+histPath)->Clone("h_wh_cTagPlus");
  TH1F* h_wh120_cTagMinus 		= (TH1F*)wh120->Get("cTagMinus/"+histPath)->Clone("h_wh_cTagMinus");
  
  
  //Hplus  M140   signal 
  TH1F* h_wh140_base  			= (TH1F*)wh140->Get("base/"+histPath)->Clone("h_wh_base"); 
  TH1F* h_wh140_JESPlus 		= (TH1F*)wh140->Get("JESPlus/"+histPath)->Clone("h_wh_JESPlus");
  TH1F* h_wh140_JESMinus 		= (TH1F*)wh140->Get("JESMinus/"+histPath)->Clone("h_wh_JESMinus");
  TH1F* h_wh140_JERPlus 		= (TH1F*)wh140->Get("JERPlus/"+histPath)->Clone("h_wh_JERPlus");
  TH1F* h_wh140_JERMinus 		= (TH1F*)wh140->Get("JERMinus/"+histPath)->Clone("h_wh_JERMinus");
  TH1F* h_wh140_TopPtPlus 		= (TH1F*)wh140->Get("TopPtPlus/"+histPath)->Clone("h_wh_TopPtPlus");
  TH1F* h_wh140_TopPtMinus 		= (TH1F*)wh140->Get("TopPtMinus/"+histPath)->Clone("h_wh_TopPtMinus");
  TH1F* h_wh140_bTagPlus 		= (TH1F*)wh140->Get("bTagPlus/"+histPath)->Clone("h_wh_bTagPlus");
  TH1F* h_wh140_bTagMinus 		= (TH1F*)wh140->Get("bTagMinus/"+histPath)->Clone("h_wh_bTagMinus");
  TH1F* h_wh140_cTagPlus 		= (TH1F*)wh140->Get("cTagPlus/"+histPath)->Clone("h_wh_cTagPlus");
  TH1F* h_wh140_cTagMinus 		= (TH1F*)wh140->Get("cTagMinus/"+histPath)->Clone("h_wh_cTagMinus");
  
  
  //Hplus  M150   signal 
  TH1F* h_wh150_base  			= (TH1F*)wh150->Get("base/"+histPath)->Clone("h_wh_base"); 
  TH1F* h_wh150_JESPlus 		= (TH1F*)wh150->Get("JESPlus/"+histPath)->Clone("h_wh_JESPlus");
  TH1F* h_wh150_JESMinus 		= (TH1F*)wh150->Get("JESMinus/"+histPath)->Clone("h_wh_JESMinus");
  TH1F* h_wh150_JERPlus 		= (TH1F*)wh150->Get("JERPlus/"+histPath)->Clone("h_wh_JERPlus");
  TH1F* h_wh150_JERMinus 		= (TH1F*)wh150->Get("JERMinus/"+histPath)->Clone("h_wh_JERMinus");
  TH1F* h_wh150_TopPtPlus 		= (TH1F*)wh150->Get("TopPtPlus/"+histPath)->Clone("h_wh_TopPtPlus");
  TH1F* h_wh150_TopPtMinus 		= (TH1F*)wh150->Get("TopPtMinus/"+histPath)->Clone("h_wh_TopPtMinus");
  TH1F* h_wh150_bTagPlus 		= (TH1F*)wh150->Get("bTagPlus/"+histPath)->Clone("h_wh_bTagPlus");
  TH1F* h_wh150_bTagMinus 		= (TH1F*)wh150->Get("bTagMinus/"+histPath)->Clone("h_wh_bTagMinus");
  TH1F* h_wh150_cTagPlus 		= (TH1F*)wh150->Get("cTagPlus/"+histPath)->Clone("h_wh_cTagPlus");
  TH1F* h_wh150_cTagMinus 		= (TH1F*)wh150->Get("cTagMinus/"+histPath)->Clone("h_wh_cTagMinus");
  
  
  //Hplus  M155   signal 
  TH1F* h_wh155_base  			= (TH1F*)wh155->Get("base/"+histPath)->Clone("h_wh_base"); 
  TH1F* h_wh155_JESPlus 		= (TH1F*)wh155->Get("JESPlus/"+histPath)->Clone("h_wh_JESPlus");
  TH1F* h_wh155_JESMinus 		= (TH1F*)wh155->Get("JESMinus/"+histPath)->Clone("h_wh_JESMinus");
  TH1F* h_wh155_JERPlus 		= (TH1F*)wh155->Get("JERPlus/"+histPath)->Clone("h_wh_JERPlus");
  TH1F* h_wh155_JERMinus 		= (TH1F*)wh155->Get("JERMinus/"+histPath)->Clone("h_wh_JERMinus");
  TH1F* h_wh155_TopPtPlus 		= (TH1F*)wh155->Get("TopPtPlus/"+histPath)->Clone("h_wh_TopPtPlus");
  TH1F* h_wh155_TopPtMinus 		= (TH1F*)wh155->Get("TopPtMinus/"+histPath)->Clone("h_wh_TopPtMinus");
  TH1F* h_wh155_bTagPlus 		= (TH1F*)wh155->Get("bTagPlus/"+histPath)->Clone("h_wh_bTagPlus");
  TH1F* h_wh155_bTagMinus 		= (TH1F*)wh155->Get("bTagMinus/"+histPath)->Clone("h_wh_bTagMinus");
  TH1F* h_wh155_cTagPlus 		= (TH1F*)wh155->Get("cTagPlus/"+histPath)->Clone("h_wh_cTagPlus");
  TH1F* h_wh155_cTagMinus 		= (TH1F*)wh155->Get("cTagMinus/"+histPath)->Clone("h_wh_cTagMinus");
  
  
  //Hplus  M160   signal 
  TH1F* h_wh160_base  			= (TH1F*)wh160->Get("base/"+histPath)->Clone("h_wh_base"); 
  TH1F* h_wh160_JESPlus 		= (TH1F*)wh160->Get("JESPlus/"+histPath)->Clone("h_wh_JESPlus");
  TH1F* h_wh160_JESMinus 		= (TH1F*)wh160->Get("JESMinus/"+histPath)->Clone("h_wh_JESMinus");
  TH1F* h_wh160_JERPlus 		= (TH1F*)wh160->Get("JERPlus/"+histPath)->Clone("h_wh_JERPlus");
  TH1F* h_wh160_JERMinus 		= (TH1F*)wh160->Get("JERMinus/"+histPath)->Clone("h_wh_JERMinus");
  TH1F* h_wh160_TopPtPlus 		= (TH1F*)wh160->Get("TopPtPlus/"+histPath)->Clone("h_wh_TopPtPlus");
  TH1F* h_wh160_TopPtMinus 		= (TH1F*)wh160->Get("TopPtMinus/"+histPath)->Clone("h_wh_TopPtMinus");
  TH1F* h_wh160_bTagPlus 		= (TH1F*)wh160->Get("bTagPlus/"+histPath)->Clone("h_wh_bTagPlus");
  TH1F* h_wh160_bTagMinus 		= (TH1F*)wh160->Get("bTagMinus/"+histPath)->Clone("h_wh_bTagMinus");
  TH1F* h_wh160_cTagPlus 		= (TH1F*)wh160->Get("cTagPlus/"+histPath)->Clone("h_wh_cTagPlus");
  TH1F* h_wh160_cTagMinus 		= (TH1F*)wh160->Get("cTagMinus/"+histPath)->Clone("h_wh_cTagMinus");
  
  

  //get average of top pt-reweighting 
  TH1F* hsig_avgTop80 = (TH1F*)wh80->Get("base/SF_topPtWeights"); 
  double sig_avgTop80 = hsig_avgTop80->GetMean();
  double sf80 = 1;
  h_wh80_base->Scale(sf80/sig_avgTop80);
  h_wh80_JESPlus->Scale(sf80/sig_avgTop80);
  h_wh80_JESMinus->Scale(sf80/sig_avgTop80);
  h_wh80_JERPlus->Scale(sf80/sig_avgTop80);
  h_wh80_JERMinus->Scale(sf80/sig_avgTop80);
  h_wh80_TopPtPlus->Scale(sf80/sig_avgTop80);
  h_wh80_TopPtMinus->Scale(sf80/sig_avgTop80);
  h_wh80_bTagPlus->Scale(sf80/sig_avgTop80);
  h_wh80_bTagMinus->Scale(sf80/sig_avgTop80);
  h_wh80_cTagPlus->Scale(sf80/sig_avgTop80);
  h_wh80_cTagMinus->Scale(sf80/sig_avgTop80);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop90 = (TH1F*)wh90->Get("base/SF_topPtWeights"); 
  double sig_avgTop90 = hsig_avgTop90->GetMean();
  double sf90 = 1.0;
  h_wh90_base->Scale(sf90/sig_avgTop90);
  h_wh90_JESPlus->Scale(sf90/sig_avgTop90);
  h_wh90_JESMinus->Scale(sf90/sig_avgTop90);
  h_wh90_JERPlus->Scale(sf90/sig_avgTop90);
  h_wh90_JERMinus->Scale(sf90/sig_avgTop90);
  h_wh90_TopPtPlus->Scale(sf90/sig_avgTop90);
  h_wh90_TopPtMinus->Scale(sf90/sig_avgTop90);
  h_wh90_bTagPlus->Scale(sf90/sig_avgTop90);
  h_wh90_bTagMinus->Scale(sf90/sig_avgTop90);
  h_wh90_cTagPlus->Scale(sf90/sig_avgTop90);
  h_wh90_cTagMinus->Scale(sf90/sig_avgTop90);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop100 = (TH1F*)wh100->Get("base/SF_topPtWeights"); 
  double sig_avgTop100 = hsig_avgTop100->GetMean();
  double sf100 = 1;
  h_wh100_base->Scale(sf100/sig_avgTop100);
  h_wh100_JESPlus->Scale(sf100/sig_avgTop100);
  h_wh100_JESMinus->Scale(sf100/sig_avgTop100);
  h_wh100_JERPlus->Scale(sf100/sig_avgTop100);
  h_wh100_JERMinus->Scale(sf100/sig_avgTop100);
  h_wh100_TopPtPlus->Scale(sf100/sig_avgTop100);
  h_wh100_TopPtMinus->Scale(sf100/sig_avgTop100);
  h_wh100_bTagPlus->Scale(sf100/sig_avgTop100);
  h_wh100_bTagMinus->Scale(sf100/sig_avgTop100);
  h_wh100_cTagPlus->Scale(sf100/sig_avgTop100);
  h_wh100_cTagMinus->Scale(sf100/sig_avgTop100);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop120 = (TH1F*)wh120->Get("base/SF_topPtWeights"); 
  double sig_avgTop120 = hsig_avgTop120->GetMean();
  double sf120 = 1;
  h_wh120_base->Scale(sf120/sig_avgTop120);
  h_wh120_JESPlus->Scale(sf120/sig_avgTop120);
  h_wh120_JESMinus->Scale(sf120/sig_avgTop120);
  h_wh120_JERPlus->Scale(sf120/sig_avgTop120);
  h_wh120_JERMinus->Scale(sf120/sig_avgTop120);
  h_wh120_TopPtPlus->Scale(sf120/sig_avgTop120);
  h_wh120_TopPtMinus->Scale(sf120/sig_avgTop120);
  h_wh120_bTagPlus->Scale(sf120/sig_avgTop120);
  h_wh120_bTagMinus->Scale(sf120/sig_avgTop120);
  h_wh120_cTagPlus->Scale(sf120/sig_avgTop120);
  h_wh120_cTagMinus->Scale(sf120/sig_avgTop120);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop140 = (TH1F*)wh140->Get("base/SF_topPtWeights"); 
  double sig_avgTop140 = hsig_avgTop140->GetMean();
  double sf140 = 1;
  h_wh140_base->Scale(sf140/sig_avgTop140);
  h_wh140_JESPlus->Scale(sf140/sig_avgTop140);
  h_wh140_JESMinus->Scale(sf140/sig_avgTop140);
  h_wh140_JERPlus->Scale(sf140/sig_avgTop140);
  h_wh140_JERMinus->Scale(sf140/sig_avgTop140);
  h_wh140_TopPtPlus->Scale(sf140/sig_avgTop140);
  h_wh140_TopPtMinus->Scale(sf140/sig_avgTop140);
  h_wh140_bTagPlus->Scale(sf140/sig_avgTop140);
  h_wh140_bTagMinus->Scale(sf140/sig_avgTop140);
  h_wh140_cTagPlus->Scale(sf140/sig_avgTop140);
  h_wh140_cTagMinus->Scale(sf140/sig_avgTop140);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop150 = (TH1F*)wh150->Get("base/SF_topPtWeights"); 
  double sig_avgTop150 = hsig_avgTop150->GetMean();
  double sf150 = 1;
  h_wh150_base->Scale(sf150/sig_avgTop150);
  h_wh150_JESPlus->Scale(sf150/sig_avgTop150);
  h_wh150_JESMinus->Scale(sf150/sig_avgTop150);
  h_wh150_JERPlus->Scale(sf150/sig_avgTop150);
  h_wh150_JERMinus->Scale(sf150/sig_avgTop150);
  h_wh150_TopPtPlus->Scale(sf150/sig_avgTop150);
  h_wh150_TopPtMinus->Scale(sf150/sig_avgTop150);
  h_wh150_bTagPlus->Scale(sf150/sig_avgTop150);
  h_wh150_bTagMinus->Scale(sf150/sig_avgTop150);
  h_wh150_cTagPlus->Scale(sf150/sig_avgTop150);
  h_wh150_cTagMinus->Scale(sf150/sig_avgTop150);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop155 = (TH1F*)wh155->Get("base/SF_topPtWeights"); 
  double sig_avgTop155 = hsig_avgTop155->GetMean();
  double sf155 = 1;
  h_wh155_base->Scale(sf155/sig_avgTop155);
  h_wh155_JESPlus->Scale(sf155/sig_avgTop155);
  h_wh155_JESMinus->Scale(sf155/sig_avgTop155);
  h_wh155_JERPlus->Scale(sf155/sig_avgTop155);
  h_wh155_JERMinus->Scale(sf155/sig_avgTop155);
  h_wh155_TopPtPlus->Scale(sf155/sig_avgTop155);
  h_wh155_TopPtMinus->Scale(sf155/sig_avgTop155);
  h_wh155_bTagPlus->Scale(sf155/sig_avgTop155);
  h_wh155_bTagMinus->Scale(sf155/sig_avgTop155);
  h_wh155_cTagPlus->Scale(sf155/sig_avgTop155);
  h_wh155_cTagMinus->Scale(sf155/sig_avgTop155);


  //get average of top pt-reweighting 
  TH1F* hsig_avgTop160 = (TH1F*)wh160->Get("base/SF_topPtWeights"); 
  double sig_avgTop160 = hsig_avgTop160->GetMean();
  double sf160 = 1;
  h_wh160_base->Scale(sf160/sig_avgTop160);
  h_wh160_JESPlus->Scale(sf160/sig_avgTop160);
  h_wh160_JESMinus->Scale(sf160/sig_avgTop160);
  h_wh160_JERPlus->Scale(sf160/sig_avgTop160);
  h_wh160_JERMinus->Scale(sf160/sig_avgTop160);
  h_wh160_TopPtPlus->Scale(sf160/sig_avgTop160);
  h_wh160_TopPtMinus->Scale(sf160/sig_avgTop160);
  h_wh160_bTagPlus->Scale(sf160/sig_avgTop160);
  h_wh160_bTagMinus->Scale(sf160/sig_avgTop160);
  h_wh160_cTagPlus->Scale(sf160/sig_avgTop160);
  h_wh160_cTagMinus->Scale(sf160/sig_avgTop160);

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
  double Nbin = 15;
  getSignalCutFlowTable(outFile, "base", Nbin, 		h_wh80_base,  		h_wh90_base,  		h_wh100_base,  	  	h_wh120_base,  	   h_wh140_base,      h_wh150_base,  	 h_wh155_base,       h_wh160_base      ); 
  getSignalCutFlowTable(outFile, "JESPlus", Nbin, 	h_wh80_JESPlus, 	h_wh90_JESPlus,         h_wh100_JESPlus,        h_wh120_JESPlus,   h_wh140_JESPlus,   h_wh150_JESPlus,   h_wh155_JESPlus,    h_wh160_JESPlus   ); 
  getSignalCutFlowTable(outFile, "JESMinus", Nbin, 	h_wh80_JESMinus, 	h_wh90_JESMinus,        h_wh100_JESMinus,       h_wh120_JESMinus,  h_wh140_JESMinus,  h_wh150_JESMinus,  h_wh155_JESMinus,   h_wh160_JESMinus  );
  getSignalCutFlowTable(outFile, "JERPlus", Nbin, 	h_wh80_JERPlus, 	h_wh90_JERPlus,         h_wh100_JERPlus,        h_wh120_JERPlus,   h_wh140_JERPlus,   h_wh150_JERPlus,   h_wh155_JERPlus,    h_wh160_JERPlus   );
  getSignalCutFlowTable(outFile, "JERMinus", Nbin, 	h_wh80_JERMinus, 	h_wh90_JERMinus,        h_wh100_JERMinus,       h_wh120_JERMinus,  h_wh140_JERMinus,  h_wh150_JERMinus,  h_wh155_JERMinus,   h_wh160_JERMinus  );
  getSignalCutFlowTable(outFile, "TopPtPlus", Nbin, 	h_wh80_TopPtPlus, 	h_wh90_TopPtPlus,       h_wh100_TopPtPlus,      h_wh120_TopPtPlus, h_wh140_TopPtPlus, h_wh150_TopPtPlus, h_wh155_TopPtPlus,  h_wh160_TopPtPlus );
  getSignalCutFlowTable(outFile, "TopPtMinus", Nbin, 	h_wh80_TopPtMinus, 	h_wh90_TopPtMinus,      h_wh100_TopPtMinus,     h_wh120_TopPtMinus,h_wh140_TopPtMinus,h_wh150_TopPtMinus,h_wh155_TopPtMinus, h_wh160_TopPtMinus);
  getSignalCutFlowTable(outFile, "bTagPlus", Nbin, 	h_wh80_bTagPlus, 	h_wh90_bTagPlus,        h_wh100_bTagPlus,       h_wh120_bTagPlus,  h_wh140_bTagPlus,  h_wh150_bTagPlus,  h_wh155_bTagPlus,   h_wh160_bTagPlus  );
  getSignalCutFlowTable(outFile, "bTagMinus", Nbin, 	h_wh80_bTagMinus, 	h_wh90_bTagMinus,       h_wh100_bTagMinus,      h_wh120_bTagMinus, h_wh140_bTagMinus, h_wh150_bTagMinus, h_wh155_bTagMinus,  h_wh160_bTagMinus );
  getSignalCutFlowTable(outFile, "cTagPlus", Nbin, 	h_wh80_cTagPlus, 	h_wh90_cTagPlus,        h_wh100_cTagPlus,       h_wh120_cTagPlus,  h_wh140_cTagPlus,  h_wh150_cTagPlus,  h_wh155_cTagPlus,   h_wh160_cTagPlus  );
  getSignalCutFlowTable(outFile, "cTagMinus", Nbin, 	h_wh80_cTagMinus, 	h_wh90_cTagMinus,       h_wh100_cTagMinus,      h_wh120_cTagMinus, h_wh140_cTagMinus, h_wh150_cTagMinus, h_wh155_cTagMinus,  h_wh160_cTagMinus );
  outFile<<"\\end{document}"<<endl;  
  outFile.close(); 
} 

void makeSummaryTableWithSys(TString histSubDir="Iso/", TString histName="cutflow"){  
  
  TString inFile("$PWD/");
  TString histPath(histSubDir+histName);
  TFile *wh80 	  		= new TFile(inFile+"all_Hplus80.root"); 
  TFile *wh90 	  		= new TFile(inFile+"all_Hplus90.root"); 
  TFile *wh100 	  		= new TFile(inFile+"all_Hplus100.root"); 
  TFile *wh120 	  		= new TFile(inFile+"all_Hplus120.root"); 
  TFile *wh140 	  		= new TFile(inFile+"all_Hplus140.root"); 
  TFile *wh150 	  		= new TFile(inFile+"all_Hplus150.root"); 
  TFile *wh155 	  		= new TFile(inFile+"all_Hplus155.root"); 
  TFile *wh160 	  		= new TFile(inFile+"all_Hplus160.root"); 
  
  TFile *ttbar    		= new TFile(inFile+"all_TTJetsP.root"); 
  TFile *wjet  			= new TFile(inFile+"all_WJets.root"); 
  TFile *zjet  			= new TFile(inFile+"all_DY.root");
  TFile *qcd  			= new TFile(inFile+"all_QCD.root");
  TFile *stop  			= new TFile(inFile+"all_ST.root");
  TFile *diboson 		= new TFile(inFile+"all_VV.root");
  
  //Hplus  wh80   signal 
  TH1F* h_wh80_base  		= (TH1F*)wh80->Get("base/"+histPath)->Clone("h_wh80_base"); 
  TH1F* h_wh80_JESPlus 		= (TH1F*)wh80->Get("JESPlus/"+histPath)->Clone("h_wh80_JESPlus");
  TH1F* h_wh80_JESMinus 		= (TH1F*)wh80->Get("JESMinus/"+histPath)->Clone("h_wh80_JESMinus");
  TH1F* h_wh80_PileupPlus 		= (TH1F*)wh80->Get("PileupPlus/"+histPath)->Clone("h_wh80_PileupPlus");
  TH1F* h_wh80_PileupMinus 		= (TH1F*)wh80->Get("PileupMinus/"+histPath)->Clone("h_wh80_PileupMinus");
  TH1F* h_wh80_JERPlus 		= (TH1F*)wh80->Get("JERPlus/"+histPath)->Clone("h_wh80_JERPlus");
  TH1F* h_wh80_JERMinus 		= (TH1F*)wh80->Get("JERMinus/"+histPath)->Clone("h_wh80_JERMinus");
  TH1F* h_wh80_TopPtPlus 		= (TH1F*)wh80->Get("TopPtPlus/"+histPath)->Clone("h_wh80_TopPtPlus");
  TH1F* h_wh80_TopPtMinus 	= (TH1F*)wh80->Get("TopPtMinus/"+histPath)->Clone("h_wh80_TopPtMinus");
  TH1F* h_wh80_bTagPlus 		= (TH1F*)wh80->Get("bTagPlus/"+histPath)->Clone("h_wh80_bTagPlus");
  TH1F* h_wh80_bTagMinus 		= (TH1F*)wh80->Get("bTagMinus/"+histPath)->Clone("h_wh80_bTagMinus");
  TH1F* h_wh80_cTagPlus 		= (TH1F*)wh80->Get("cTagPlus/"+histPath)->Clone("h_wh80_cTagPlus");
  TH1F* h_wh80_cTagMinus 		= (TH1F*)wh80->Get("cTagMinus/"+histPath)->Clone("h_wh80_cTagMinus");
  
  //Hplus  wh90   signal 
  TH1F* h_wh90_base  		= (TH1F*)wh90->Get("base/"+histPath)->Clone("h_wh90_base"); 
  TH1F* h_wh90_JESPlus 		= (TH1F*)wh90->Get("JESPlus/"+histPath)->Clone("h_wh90_JESPlus");
  TH1F* h_wh90_JESMinus 		= (TH1F*)wh90->Get("JESMinus/"+histPath)->Clone("h_wh90_JESMinus");
  TH1F* h_wh90_PileupPlus 		= (TH1F*)wh90->Get("PileupPlus/"+histPath)->Clone("h_wh90_PileupPlus");
  TH1F* h_wh90_PileupMinus 		= (TH1F*)wh90->Get("PileupMinus/"+histPath)->Clone("h_wh90_PileupMinus");
  TH1F* h_wh90_JERPlus 		= (TH1F*)wh90->Get("JERPlus/"+histPath)->Clone("h_wh90_JERPlus");
  TH1F* h_wh90_JERMinus 		= (TH1F*)wh90->Get("JERMinus/"+histPath)->Clone("h_wh90_JERMinus");
  TH1F* h_wh90_TopPtPlus 		= (TH1F*)wh90->Get("TopPtPlus/"+histPath)->Clone("h_wh90_TopPtPlus");
  TH1F* h_wh90_TopPtMinus 	= (TH1F*)wh90->Get("TopPtMinus/"+histPath)->Clone("h_wh90_TopPtMinus");
  TH1F* h_wh90_bTagPlus 		= (TH1F*)wh90->Get("bTagPlus/"+histPath)->Clone("h_wh90_bTagPlus");
  TH1F* h_wh90_bTagMinus 		= (TH1F*)wh90->Get("bTagMinus/"+histPath)->Clone("h_wh90_bTagMinus");
  TH1F* h_wh90_cTagPlus 		= (TH1F*)wh90->Get("cTagPlus/"+histPath)->Clone("h_wh90_cTagPlus");
  TH1F* h_wh90_cTagMinus 		= (TH1F*)wh90->Get("cTagMinus/"+histPath)->Clone("h_wh90_cTagMinus");
  
  
  //Hplus  wh100   signal 
  TH1F* h_wh100_base  		= (TH1F*)wh100->Get("base/"+histPath)->Clone("h_wh100_base"); 
  TH1F* h_wh100_JESPlus 		= (TH1F*)wh100->Get("JESPlus/"+histPath)->Clone("h_wh100_JESPlus");
  TH1F* h_wh100_JESMinus 		= (TH1F*)wh100->Get("JESMinus/"+histPath)->Clone("h_wh100_JESMinus");
  TH1F* h_wh100_PileupPlus 		= (TH1F*)wh100->Get("PileupPlus/"+histPath)->Clone("h_wh100_PileupPlus");
  TH1F* h_wh100_PileupMinus 		= (TH1F*)wh100->Get("PileupMinus/"+histPath)->Clone("h_wh100_PileupMinus");
  TH1F* h_wh100_JERPlus 		= (TH1F*)wh100->Get("JERPlus/"+histPath)->Clone("h_wh100_JERPlus");
  TH1F* h_wh100_JERMinus 		= (TH1F*)wh100->Get("JERMinus/"+histPath)->Clone("h_wh100_JERMinus");
  TH1F* h_wh100_TopPtPlus 		= (TH1F*)wh100->Get("TopPtPlus/"+histPath)->Clone("h_wh100_TopPtPlus");
  TH1F* h_wh100_TopPtMinus 	= (TH1F*)wh100->Get("TopPtMinus/"+histPath)->Clone("h_wh100_TopPtMinus");
  TH1F* h_wh100_bTagPlus 		= (TH1F*)wh100->Get("bTagPlus/"+histPath)->Clone("h_wh100_bTagPlus");
  TH1F* h_wh100_bTagMinus 		= (TH1F*)wh100->Get("bTagMinus/"+histPath)->Clone("h_wh100_bTagMinus");
  TH1F* h_wh100_cTagPlus 		= (TH1F*)wh100->Get("cTagPlus/"+histPath)->Clone("h_wh100_cTagPlus");
  TH1F* h_wh100_cTagMinus 		= (TH1F*)wh100->Get("cTagMinus/"+histPath)->Clone("h_wh100_cTagMinus");
  
  //Hplus  wh120   signal 
  TH1F* h_wh120_base  		= (TH1F*)wh120->Get("base/"+histPath)->Clone("h_wh120_base"); 
  TH1F* h_wh120_JESPlus 		= (TH1F*)wh120->Get("JESPlus/"+histPath)->Clone("h_wh120_JESPlus");
  TH1F* h_wh120_JESMinus 		= (TH1F*)wh120->Get("JESMinus/"+histPath)->Clone("h_wh120_JESMinus");
  TH1F* h_wh120_PileupPlus 		= (TH1F*)wh120->Get("PileupPlus/"+histPath)->Clone("h_wh120_PileupPlus");
  TH1F* h_wh120_PileupMinus 		= (TH1F*)wh120->Get("PileupMinus/"+histPath)->Clone("h_wh120_PileupMinus");
  TH1F* h_wh120_JERPlus 		= (TH1F*)wh120->Get("JERPlus/"+histPath)->Clone("h_wh120_JERPlus");
  TH1F* h_wh120_JERMinus 		= (TH1F*)wh120->Get("JERMinus/"+histPath)->Clone("h_wh120_JERMinus");
  TH1F* h_wh120_TopPtPlus 		= (TH1F*)wh120->Get("TopPtPlus/"+histPath)->Clone("h_wh120_TopPtPlus");
  TH1F* h_wh120_TopPtMinus 	= (TH1F*)wh120->Get("TopPtMinus/"+histPath)->Clone("h_wh120_TopPtMinus");
  TH1F* h_wh120_bTagPlus 		= (TH1F*)wh120->Get("bTagPlus/"+histPath)->Clone("h_wh120_bTagPlus");
  TH1F* h_wh120_bTagMinus 		= (TH1F*)wh120->Get("bTagMinus/"+histPath)->Clone("h_wh120_bTagMinus");
  TH1F* h_wh120_cTagPlus 		= (TH1F*)wh120->Get("cTagPlus/"+histPath)->Clone("h_wh120_cTagPlus");
  TH1F* h_wh120_cTagMinus 		= (TH1F*)wh120->Get("cTagMinus/"+histPath)->Clone("h_wh120_cTagMinus");
  
  //Hplus  wh140   signal 
  TH1F* h_wh140_base  		= (TH1F*)wh140->Get("base/"+histPath)->Clone("h_wh140_base"); 
  TH1F* h_wh140_JESPlus 		= (TH1F*)wh140->Get("JESPlus/"+histPath)->Clone("h_wh140_JESPlus");
  TH1F* h_wh140_JESMinus 		= (TH1F*)wh140->Get("JESMinus/"+histPath)->Clone("h_wh140_JESMinus");
  TH1F* h_wh140_PileupPlus 		= (TH1F*)wh140->Get("PileupPlus/"+histPath)->Clone("h_wh140_PileupPlus");
  TH1F* h_wh140_PileupMinus 		= (TH1F*)wh140->Get("PileupMinus/"+histPath)->Clone("h_wh140_PileupMinus");
  TH1F* h_wh140_JERPlus 		= (TH1F*)wh140->Get("JERPlus/"+histPath)->Clone("h_wh140_JERPlus");
  TH1F* h_wh140_JERMinus 		= (TH1F*)wh140->Get("JERMinus/"+histPath)->Clone("h_wh140_JERMinus");
  TH1F* h_wh140_TopPtPlus 		= (TH1F*)wh140->Get("TopPtPlus/"+histPath)->Clone("h_wh140_TopPtPlus");
  TH1F* h_wh140_TopPtMinus 	= (TH1F*)wh140->Get("TopPtMinus/"+histPath)->Clone("h_wh140_TopPtMinus");
  TH1F* h_wh140_bTagPlus 		= (TH1F*)wh140->Get("bTagPlus/"+histPath)->Clone("h_wh140_bTagPlus");
  TH1F* h_wh140_bTagMinus 		= (TH1F*)wh140->Get("bTagMinus/"+histPath)->Clone("h_wh140_bTagMinus");
  TH1F* h_wh140_cTagPlus 		= (TH1F*)wh140->Get("cTagPlus/"+histPath)->Clone("h_wh140_cTagPlus");
  TH1F* h_wh140_cTagMinus 		= (TH1F*)wh140->Get("cTagMinus/"+histPath)->Clone("h_wh140_cTagMinus");
  
  
  //Hplus  wh150   signal 
  TH1F* h_wh150_base  		= (TH1F*)wh150->Get("base/"+histPath)->Clone("h_wh150_base"); 
  TH1F* h_wh150_JESPlus 		= (TH1F*)wh150->Get("JESPlus/"+histPath)->Clone("h_wh150_JESPlus");
  TH1F* h_wh150_JESMinus 		= (TH1F*)wh150->Get("JESMinus/"+histPath)->Clone("h_wh150_JESMinus");
  TH1F* h_wh150_PileupPlus 		= (TH1F*)wh150->Get("PileupPlus/"+histPath)->Clone("h_wh150_PileupPlus");
  TH1F* h_wh150_PileupMinus 		= (TH1F*)wh150->Get("PileupMinus/"+histPath)->Clone("h_wh150_PileupMinus");
  TH1F* h_wh150_JERPlus 		= (TH1F*)wh150->Get("JERPlus/"+histPath)->Clone("h_wh150_JERPlus");
  TH1F* h_wh150_JERMinus 		= (TH1F*)wh150->Get("JERMinus/"+histPath)->Clone("h_wh150_JERMinus");
  TH1F* h_wh150_TopPtPlus 		= (TH1F*)wh150->Get("TopPtPlus/"+histPath)->Clone("h_wh150_TopPtPlus");
  TH1F* h_wh150_TopPtMinus 	= (TH1F*)wh150->Get("TopPtMinus/"+histPath)->Clone("h_wh150_TopPtMinus");
  TH1F* h_wh150_bTagPlus 		= (TH1F*)wh150->Get("bTagPlus/"+histPath)->Clone("h_wh150_bTagPlus");
  TH1F* h_wh150_bTagMinus 		= (TH1F*)wh150->Get("bTagMinus/"+histPath)->Clone("h_wh150_bTagMinus");
  TH1F* h_wh150_cTagPlus 		= (TH1F*)wh150->Get("cTagPlus/"+histPath)->Clone("h_wh150_cTagPlus");
  TH1F* h_wh150_cTagMinus 		= (TH1F*)wh150->Get("cTagMinus/"+histPath)->Clone("h_wh150_cTagMinus");
  
  //Hplus  wh155   signal 
  TH1F* h_wh155_base  		= (TH1F*)wh155->Get("base/"+histPath)->Clone("h_wh155_base"); 
  TH1F* h_wh155_JESPlus 		= (TH1F*)wh155->Get("JESPlus/"+histPath)->Clone("h_wh155_JESPlus");
  TH1F* h_wh155_JESMinus 		= (TH1F*)wh155->Get("JESMinus/"+histPath)->Clone("h_wh155_JESMinus");
  TH1F* h_wh155_PileupPlus 		= (TH1F*)wh155->Get("PileupPlus/"+histPath)->Clone("h_wh155_PileupPlus");
  TH1F* h_wh155_PileupMinus 		= (TH1F*)wh155->Get("PileupMinus/"+histPath)->Clone("h_wh155_PileupMinus");
  TH1F* h_wh155_JERPlus 		= (TH1F*)wh155->Get("JERPlus/"+histPath)->Clone("h_wh155_JERPlus");
  TH1F* h_wh155_JERMinus 		= (TH1F*)wh155->Get("JERMinus/"+histPath)->Clone("h_wh155_JERMinus");
  TH1F* h_wh155_TopPtPlus 		= (TH1F*)wh155->Get("TopPtPlus/"+histPath)->Clone("h_wh155_TopPtPlus");
  TH1F* h_wh155_TopPtMinus 	= (TH1F*)wh155->Get("TopPtMinus/"+histPath)->Clone("h_wh155_TopPtMinus");
  TH1F* h_wh155_bTagPlus 		= (TH1F*)wh155->Get("bTagPlus/"+histPath)->Clone("h_wh155_bTagPlus");
  TH1F* h_wh155_bTagMinus 		= (TH1F*)wh155->Get("bTagMinus/"+histPath)->Clone("h_wh155_bTagMinus");
  TH1F* h_wh155_cTagPlus 		= (TH1F*)wh155->Get("cTagPlus/"+histPath)->Clone("h_wh155_cTagPlus");
  TH1F* h_wh155_cTagMinus 		= (TH1F*)wh155->Get("cTagMinus/"+histPath)->Clone("h_wh155_cTagMinus");
  
  //Hplus  wh160   signal 
  TH1F* h_wh160_base  		= (TH1F*)wh160->Get("base/"+histPath)->Clone("h_wh160_base"); 
  TH1F* h_wh160_JESPlus 		= (TH1F*)wh160->Get("JESPlus/"+histPath)->Clone("h_wh160_JESPlus");
  TH1F* h_wh160_JESMinus 		= (TH1F*)wh160->Get("JESMinus/"+histPath)->Clone("h_wh160_JESMinus");
  TH1F* h_wh160_PileupPlus 		= (TH1F*)wh160->Get("PileupPlus/"+histPath)->Clone("h_wh160_PileupPlus");
  TH1F* h_wh160_PileupMinus 		= (TH1F*)wh160->Get("PileupMinus/"+histPath)->Clone("h_wh160_PileupMinus");
  TH1F* h_wh160_JERPlus 		= (TH1F*)wh160->Get("JERPlus/"+histPath)->Clone("h_wh160_JERPlus");
  TH1F* h_wh160_JERMinus 		= (TH1F*)wh160->Get("JERMinus/"+histPath)->Clone("h_wh160_JERMinus");
  TH1F* h_wh160_TopPtPlus 		= (TH1F*)wh160->Get("TopPtPlus/"+histPath)->Clone("h_wh160_TopPtPlus");
  TH1F* h_wh160_TopPtMinus 	= (TH1F*)wh160->Get("TopPtMinus/"+histPath)->Clone("h_wh160_TopPtMinus");
  TH1F* h_wh160_bTagPlus 		= (TH1F*)wh160->Get("bTagPlus/"+histPath)->Clone("h_wh160_bTagPlus");
  TH1F* h_wh160_bTagMinus 		= (TH1F*)wh160->Get("bTagMinus/"+histPath)->Clone("h_wh160_bTagMinus");
  TH1F* h_wh160_cTagPlus 		= (TH1F*)wh160->Get("cTagPlus/"+histPath)->Clone("h_wh160_cTagPlus");
  TH1F* h_wh160_cTagMinus 		= (TH1F*)wh160->Get("cTagMinus/"+histPath)->Clone("h_wh160_cTagMinus");
  
  
  //ttbar+jets
  TH1F* h_ttbar_base  		= (TH1F*)ttbar->Get("base/"+histPath)->Clone("h_ttbar_base");  
  TH1F* h_ttbar_JESPlus 	= (TH1F*)ttbar->Get("JESPlus/"+histPath)->Clone("h_ttbar_JESPlus"); 
  TH1F* h_ttbar_JESMinus 	= (TH1F*)ttbar->Get("JESMinus/"+histPath)->Clone("h_ttbar_JESMinus"); 
  TH1F* h_ttbar_PileupPlus 	= (TH1F*)ttbar->Get("PileupPlus/"+histPath)->Clone("h_ttbar_PileupPlus"); 
  TH1F* h_ttbar_PileupMinus 	= (TH1F*)ttbar->Get("PileupMinus/"+histPath)->Clone("h_ttbar_PileupMinus"); 
  TH1F* h_ttbar_JERPlus 	= (TH1F*)ttbar->Get("JERPlus/"+histPath)->Clone("h_ttbar_JERPlus"); 
  TH1F* h_ttbar_JERMinus 	= (TH1F*)ttbar->Get("JERMinus/"+histPath)->Clone("h_ttbar_JERMinus"); 
  TH1F* h_ttbar_TopPtPlus 	= (TH1F*)ttbar->Get("TopPtPlus/"+histPath)->Clone("h_ttbar_TopPtPlus"); 
  TH1F* h_ttbar_TopPtMinus 	= (TH1F*)ttbar->Get("TopPtMinus/"+histPath)->Clone("h_ttbar_TopPtMinus"); 
  TH1F* h_ttbar_bTagPlus 	= (TH1F*)ttbar->Get("bTagPlus/"+histPath)->Clone("h_ttbar_bTagPlus"); 
  TH1F* h_ttbar_bTagMinus 	= (TH1F*)ttbar->Get("bTagMinus/"+histPath)->Clone("h_ttbar_bTagMinus"); 
  TH1F* h_ttbar_cTagPlus 	= (TH1F*)ttbar->Get("cTagPlus/"+histPath)->Clone("h_ttbar_cTagPlus"); 
  TH1F* h_ttbar_cTagMinus 	= (TH1F*)ttbar->Get("cTagMinus/"+histPath)->Clone("h_ttbar_cTagMinus"); 
  //w+jets
  TH1F* h_wjet_base  		= (TH1F*)wjet->Get("base/"+histPath)->Clone("h_wjet_base");  
  TH1F* h_wjet_JESPlus 		= (TH1F*)wjet->Get("JESPlus/"+histPath)->Clone("h_wjet_JESPlus"); 
  TH1F* h_wjet_JESMinus 	= (TH1F*)wjet->Get("JESMinus/"+histPath)->Clone("h_wjet_JESMinus"); 
  TH1F* h_wjet_PileupPlus 		= (TH1F*)wjet->Get("PileupPlus/"+histPath)->Clone("h_wjet_PileupPlus"); 
  TH1F* h_wjet_PileupMinus 	= (TH1F*)wjet->Get("PileupMinus/"+histPath)->Clone("h_wjet_PileupMinus"); 
  TH1F* h_wjet_JERPlus 		= (TH1F*)wjet->Get("JERPlus/"+histPath)->Clone("h_wjet_JERPlus"); 
  TH1F* h_wjet_JERMinus 	= (TH1F*)wjet->Get("JERMinus/"+histPath)->Clone("h_wjet_JERMinus"); 
  TH1F* h_wjet_TopPtPlus 	= (TH1F*)wjet->Get("TopPtPlus/"+histPath)->Clone("h_wjet_TopPtPlus"); 
  TH1F* h_wjet_TopPtMinus 	= (TH1F*)wjet->Get("TopPtMinus/"+histPath)->Clone("h_wjet_TopPtMinus"); 
  TH1F* h_wjet_bTagPlus 	= (TH1F*)wjet->Get("bTagPlus/"+histPath)->Clone("h_wjet_bTagPlus"); 
  TH1F* h_wjet_bTagMinus 	= (TH1F*)wjet->Get("bTagMinus/"+histPath)->Clone("h_wjet_bTagMinus"); 
  TH1F* h_wjet_cTagPlus 	= (TH1F*)wjet->Get("cTagPlus/"+histPath)->Clone("h_wjet_cTagPlus"); 
  TH1F* h_wjet_cTagMinus 	= (TH1F*)wjet->Get("cTagMinus/"+histPath)->Clone("h_wjet_cTagMinus"); 
  //dy+jets
  TH1F* h_zjet_base  		= (TH1F*)zjet->Get("base/"+histPath)->Clone("h_zjet_base");
  TH1F* h_zjet_JESPlus 		= (TH1F*)zjet->Get("JESPlus/"+histPath)->Clone("h_zjet_JESPlus");
  TH1F* h_zjet_JESMinus 	= (TH1F*)zjet->Get("JESMinus/"+histPath)->Clone("h_zjet_JESMinus");
  TH1F* h_zjet_PileupPlus 		= (TH1F*)zjet->Get("PileupPlus/"+histPath)->Clone("h_zjet_PileupPlus");
  TH1F* h_zjet_PileupMinus 	= (TH1F*)zjet->Get("PileupMinus/"+histPath)->Clone("h_zjet_PileupMinus");
  TH1F* h_zjet_JERPlus 		= (TH1F*)zjet->Get("JERPlus/"+histPath)->Clone("h_zjet_JERPlus");
  TH1F* h_zjet_JERMinus 	= (TH1F*)zjet->Get("JERMinus/"+histPath)->Clone("h_zjet_JERMinus");
  TH1F* h_zjet_TopPtPlus 	= (TH1F*)zjet->Get("TopPtPlus/"+histPath)->Clone("h_zjet_TopPtPlus");
  TH1F* h_zjet_TopPtMinus 	= (TH1F*)zjet->Get("TopPtMinus/"+histPath)->Clone("h_zjet_TopPtMinus");
  TH1F* h_zjet_bTagPlus 	= (TH1F*)zjet->Get("bTagPlus/"+histPath)->Clone("h_zjet_bTagPlus");
  TH1F* h_zjet_bTagMinus 	= (TH1F*)zjet->Get("bTagMinus/"+histPath)->Clone("h_zjet_bTagMinus");
  TH1F* h_zjet_cTagPlus 	= (TH1F*)zjet->Get("cTagPlus/"+histPath)->Clone("h_zjet_cTagPlus");
  TH1F* h_zjet_cTagMinus 	= (TH1F*)zjet->Get("cTagMinus/"+histPath)->Clone("h_zjet_cTagMinus");
  //qcd
  TH1F* h_qcd_base 		= (TH1F*)qcd->Get("base/"+histPath)->Clone("h_qcd_base");
  TH1F* h_qcd_JESPlus 		= (TH1F*)qcd->Get("JESPlus/"+histPath)->Clone("h_qcd_JESPlus");
  TH1F* h_qcd_JESMinus 		= (TH1F*)qcd->Get("JESMinus/"+histPath)->Clone("h_qcd_JESMinus");
  TH1F* h_qcd_PileupPlus 		= (TH1F*)qcd->Get("PileupPlus/"+histPath)->Clone("h_qcd_PileupPlus");
  TH1F* h_qcd_PileupMinus 		= (TH1F*)qcd->Get("PileupMinus/"+histPath)->Clone("h_qcd_PileupMinus");
  TH1F* h_qcd_JERPlus 		= (TH1F*)qcd->Get("JERPlus/"+histPath)->Clone("h_qcd_JERPlus");
  TH1F* h_qcd_JERMinus 		= (TH1F*)qcd->Get("JERMinus/"+histPath)->Clone("h_qcd_JERMinus");
  TH1F* h_qcd_TopPtPlus 	= (TH1F*)qcd->Get("TopPtPlus/"+histPath)->Clone("h_qcd_TopPtPlus");
  TH1F* h_qcd_TopPtMinus 	= (TH1F*)qcd->Get("TopPtMinus/"+histPath)->Clone("h_qcd_TopPtMinus");
  TH1F* h_qcd_bTagPlus 		= (TH1F*)qcd->Get("bTagPlus/"+histPath)->Clone("h_qcd_bTagPlus");
  TH1F* h_qcd_bTagMinus 	= (TH1F*)qcd->Get("bTagMinus/"+histPath)->Clone("h_qcd_bTagMinus");
  TH1F* h_qcd_cTagPlus 		= (TH1F*)qcd->Get("cTagPlus/"+histPath)->Clone("h_qcd_cTagPlus");
  TH1F* h_qcd_cTagMinus 	= (TH1F*)qcd->Get("cTagMinus/"+histPath)->Clone("h_qcd_cTagMinus");
  //single top
  TH1F* h_stop_base  		= (TH1F*)stop->Get("base/"+histPath)->Clone("h_stop_base");
  TH1F* h_stop_JESPlus 		= (TH1F*)stop->Get("JESPlus/"+histPath)->Clone("h_stop_JESPlus");
  TH1F* h_stop_JESMinus 	= (TH1F*)stop->Get("JESMinus/"+histPath)->Clone("h_stop_JESMinus");
  TH1F* h_stop_PileupPlus 		= (TH1F*)stop->Get("PileupPlus/"+histPath)->Clone("h_stop_PileupPlus");
  TH1F* h_stop_PileupMinus 	= (TH1F*)stop->Get("PileupMinus/"+histPath)->Clone("h_stop_PileupMinus");
  TH1F* h_stop_JERPlus 		= (TH1F*)stop->Get("JERPlus/"+histPath)->Clone("h_stop_JERPlus");
  TH1F* h_stop_JERMinus 	= (TH1F*)stop->Get("JERMinus/"+histPath)->Clone("h_stop_JERMinus");
  TH1F* h_stop_TopPtPlus 	= (TH1F*)stop->Get("TopPtPlus/"+histPath)->Clone("h_stop_TopPtPlus");
  TH1F* h_stop_TopPtMinus 	= (TH1F*)stop->Get("TopPtMinus/"+histPath)->Clone("h_stop_TopPtMinus");
  TH1F* h_stop_bTagPlus 	= (TH1F*)stop->Get("bTagPlus/"+histPath)->Clone("h_stop_bTagPlus");
  TH1F* h_stop_bTagMinus 	= (TH1F*)stop->Get("bTagMinus/"+histPath)->Clone("h_stop_bTagMinus");
  TH1F* h_stop_cTagPlus 	= (TH1F*)stop->Get("cTagPlus/"+histPath)->Clone("h_stop_cTagPlus");
  TH1F* h_stop_cTagMinus 	= (TH1F*)stop->Get("cTagMinus/"+histPath)->Clone("h_stop_cTagMinus");
  //vv
  TH1F* h_diboson_base 		= (TH1F*)diboson->Get("base/"+histPath)->Clone("h_diboson_base");
  TH1F* h_diboson_JESPlus 	= (TH1F*)diboson->Get("JESPlus/"+histPath)->Clone("h_diboson_JESPlus");
  TH1F* h_diboson_JESMinus 	= (TH1F*)diboson->Get("JESMinus/"+histPath)->Clone("h_diboson_JESMinus");
  TH1F* h_diboson_PileupPlus 	= (TH1F*)diboson->Get("PileupPlus/"+histPath)->Clone("h_diboson_PileupPlus");
  TH1F* h_diboson_PileupMinus 	= (TH1F*)diboson->Get("PileupMinus/"+histPath)->Clone("h_diboson_PileupMinus");
  TH1F* h_diboson_JERPlus 	= (TH1F*)diboson->Get("JERPlus/"+histPath)->Clone("h_diboson_JERPlus");
  TH1F* h_diboson_JERMinus 	= (TH1F*)diboson->Get("JERMinus/"+histPath)->Clone("h_diboson_JERMinus");
  TH1F* h_diboson_TopPtPlus 	= (TH1F*)diboson->Get("TopPtPlus/"+histPath)->Clone("h_diboson_TopPtPlus");
  TH1F* h_diboson_TopPtMinus 	= (TH1F*)diboson->Get("TopPtMinus/"+histPath)->Clone("h_diboson_TopPtMinus");
  TH1F* h_diboson_bTagPlus 	= (TH1F*)diboson->Get("bTagPlus/"+histPath)->Clone("h_diboson_bTagPlus");
  TH1F* h_diboson_bTagMinus 	= (TH1F*)diboson->Get("bTagMinus/"+histPath)->Clone("h_diboson_bTagMinus");
  TH1F* h_diboson_cTagPlus 	= (TH1F*)diboson->Get("cTagPlus/"+histPath)->Clone("h_diboson_cTagPlus");
  TH1F* h_diboson_cTagMinus 	= (TH1F*)diboson->Get("cTagMinus/"+histPath)->Clone("h_diboson_cTagMinus");

  //get average of top pt-reweighting 
  TH1F* hsig_avgTopwh80 = (TH1F*)wh80->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh80 = hsig_avgTopwh80->GetMean();
  double sfwh80 = 1;
  h_wh80_base->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_JESPlus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_JESMinus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_PileupPlus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_PileupMinus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_JERPlus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_JERMinus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_TopPtPlus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_TopPtMinus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_bTagPlus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_bTagMinus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_cTagPlus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_cTagMinus->Scale(sfwh80/sig_avgTopwh80);

  TH1F* hsig_avgTopwh90 = (TH1F*)wh90->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh90 = hsig_avgTopwh90->GetMean();
  double sfwh90 = 1;
  h_wh90_base->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_JESPlus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_JESMinus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_PileupPlus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_PileupMinus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_JERPlus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_JERMinus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_TopPtPlus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_TopPtMinus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_bTagPlus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_bTagMinus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_cTagPlus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_cTagMinus->Scale(sfwh90/sig_avgTopwh90);

  TH1F* hsig_avgTopwh100 = (TH1F*)wh100->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh100 = hsig_avgTopwh100->GetMean();
  double sfwh100 = 1;
  h_wh100_base->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_JESPlus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_JESMinus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_PileupPlus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_PileupMinus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_JERPlus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_JERMinus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_TopPtPlus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_TopPtMinus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_bTagPlus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_bTagMinus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_cTagPlus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_cTagMinus->Scale(sfwh100/sig_avgTopwh100);

  TH1F* hsig_avgTopwh120 = (TH1F*)wh120->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh120 = hsig_avgTopwh120->GetMean();
  double sfwh120 = 1;
  h_wh120_base->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_JESPlus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_JESMinus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_PileupPlus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_PileupMinus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_JERPlus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_JERMinus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_TopPtPlus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_TopPtMinus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_bTagPlus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_bTagMinus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_cTagPlus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_cTagMinus->Scale(sfwh120/sig_avgTopwh120);

  TH1F* hsig_avgTopwh140 = (TH1F*)wh140->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh140 = hsig_avgTopwh140->GetMean();
  double sfwh140 = 1;
  h_wh140_base->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_JESPlus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_JESMinus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_PileupPlus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_PileupMinus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_JERPlus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_JERMinus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_TopPtPlus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_TopPtMinus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_bTagPlus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_bTagMinus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_cTagPlus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_cTagMinus->Scale(sfwh140/sig_avgTopwh140);

  TH1F* hsig_avgTopwh150 = (TH1F*)wh150->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh150 = hsig_avgTopwh150->GetMean();
  double sfwh150 = 1;
  h_wh150_base->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_JESPlus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_JESMinus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_PileupPlus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_PileupMinus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_JERPlus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_JERMinus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_TopPtPlus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_TopPtMinus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_bTagPlus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_bTagMinus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_cTagPlus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_cTagMinus->Scale(sfwh150/sig_avgTopwh150);

  TH1F* hsig_avgTopwh155 = (TH1F*)wh155->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh155 = hsig_avgTopwh155->GetMean();
  double sfwh155 = 1;
  h_wh155_base->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_JESPlus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_JESMinus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_PileupPlus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_PileupMinus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_JERPlus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_JERMinus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_TopPtPlus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_TopPtMinus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_bTagPlus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_bTagMinus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_cTagPlus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_cTagMinus->Scale(sfwh155/sig_avgTopwh155);

  TH1F* hsig_avgTopwh160 = (TH1F*)wh160->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh160 = hsig_avgTopwh160->GetMean();
  double sfwh160 = 1;
  h_wh160_base->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_JESPlus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_JESMinus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_PileupPlus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_PileupMinus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_JERPlus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_JERMinus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_TopPtPlus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_TopPtMinus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_bTagPlus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_bTagMinus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_cTagPlus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_cTagMinus->Scale(sfwh160/sig_avgTopwh160);

  
  TH1F* httbar_avgTop = (TH1F*)ttbar->Get("base/SF_topPtWeights"); 
  double ttbar_avgTop = httbar_avgTop->GetMean();
  double sf_top = 1;
  h_ttbar_base->Scale(sf_top/ttbar_avgTop);
  h_ttbar_JESPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_JESMinus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_PileupPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_PileupMinus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_JERPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_JERMinus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_TopPtPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_TopPtMinus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_bTagPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_bTagMinus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_cTagPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_cTagMinus->Scale(sf_top/ttbar_avgTop);

  TH1F* h_TotalBkg = (TH1F*)h_wjet_base->Clone("h_TotalBkg");
  h_TotalBkg->Reset();
  h_TotalBkg->Add(h_ttbar_base);
  h_TotalBkg->Add(h_wjet_base);
  h_TotalBkg->Add(h_zjet_base);
  h_TotalBkg->Add(h_qcd_base);
  h_TotalBkg->Add(h_stop_base);
  h_TotalBkg->Add(h_diboson_base);


  TFile *data = new TFile(inFile+"all_muData.root");
  TH1F* h_data = (TH1F*)data->Get("base/"+histPath)->Clone("h_data");


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
  outFile<<"\\multicolumn{1}{|c|}{Source} & \\multicolumn{1}{|c|}{N$_{\r events}$ $\\pm$ MC stat $\\pm$ JEC/MET/Top $\\pm$ bTag $\\pm$ cTag $\\pm$ Pileup} \\\\"<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<"\\hline "<<endl;
  int b = 14;

  
  outFile<<"HW, $M_{H}=80~GeV/c^{2}$"<< " & "<<h_wh80_base->GetBinContent(b)<<" $\\pm$ "<<h_wh80_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_wh80_JESPlus, h_wh80_base, h_wh80_JESMinus, h_wh80_JERPlus, h_wh80_JERMinus, h_wh80_TopPtPlus, h_wh80_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_wh80_bTagPlus, h_wh80_base, h_wh80_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_wh80_cTagPlus, h_wh80_base, h_wh80_cTagMinus, b)<<" $\\pm$ "<< sysUncBCTag (h_wh80_PileupPlus, h_wh80_base, h_wh80_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  
  outFile<<"HW, $M_{H}=90~GeV/c^{2}$"<< " & "<<h_wh90_base->GetBinContent(b)<<" $\\pm$ "<<h_wh90_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_wh90_JESPlus, h_wh90_base, h_wh90_JESMinus, h_wh90_JERPlus, h_wh90_JERMinus, h_wh90_TopPtPlus, h_wh90_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_wh90_bTagPlus, h_wh90_base, h_wh90_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_wh90_cTagPlus, h_wh90_base, h_wh90_cTagMinus, b)<<" $\\pm$ "<< sysUncBCTag (h_wh90_PileupPlus, h_wh90_base, h_wh90_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  
  outFile<<"HW, $M_{H}=100~GeV/c^{2}$"<< " & "<<h_wh100_base->GetBinContent(b)<<" $\\pm$ "<<h_wh100_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_wh100_JESPlus, h_wh100_base, h_wh100_JESMinus, h_wh100_JERPlus, h_wh100_JERMinus, h_wh100_TopPtPlus, h_wh100_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_wh100_bTagPlus, h_wh100_base, h_wh100_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_wh100_cTagPlus, h_wh100_base, h_wh100_cTagMinus, b)<<" $\\pm$ "<< sysUncBCTag (h_wh100_PileupPlus, h_wh100_base, h_wh100_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  
  outFile<<"HW, $M_{H}=120~GeV/c^{2}$"<< " & "<<h_wh120_base->GetBinContent(b)<<" $\\pm$ "<<h_wh120_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_wh120_JESPlus, h_wh120_base, h_wh120_JESMinus, h_wh120_JERPlus, h_wh120_JERMinus, h_wh120_TopPtPlus, h_wh120_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_wh120_bTagPlus, h_wh120_base, h_wh120_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_wh120_cTagPlus, h_wh120_base, h_wh120_cTagMinus, b)<<" $\\pm$ "<< sysUncBCTag (h_wh120_PileupPlus, h_wh120_base, h_wh120_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  
  outFile<<"HW, $M_{H}=140~GeV/c^{2}$"<< " & "<<h_wh140_base->GetBinContent(b)<<" $\\pm$ "<<h_wh140_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_wh140_JESPlus, h_wh140_base, h_wh140_JESMinus, h_wh140_JERPlus, h_wh140_JERMinus, h_wh140_TopPtPlus, h_wh140_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_wh140_bTagPlus, h_wh140_base, h_wh140_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_wh140_cTagPlus, h_wh140_base, h_wh140_cTagMinus, b)<<" $\\pm$ "<< sysUncBCTag (h_wh140_PileupPlus, h_wh140_base, h_wh140_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  
  outFile<<"HW, $M_{H}=150~GeV/c^{2}$"<< " & "<<h_wh150_base->GetBinContent(b)<<" $\\pm$ "<<h_wh150_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_wh150_JESPlus, h_wh150_base, h_wh150_JESMinus, h_wh150_JERPlus, h_wh150_JERMinus, h_wh150_TopPtPlus, h_wh150_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_wh150_bTagPlus, h_wh150_base, h_wh150_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_wh150_cTagPlus, h_wh150_base, h_wh150_cTagMinus, b)<<" $\\pm$ "<< sysUncBCTag (h_wh150_PileupPlus, h_wh150_base, h_wh150_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  
  outFile<<"HW, $M_{H}=155~GeV/c^{2}$"<< " & "<<h_wh155_base->GetBinContent(b)<<" $\\pm$ "<<h_wh155_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_wh155_JESPlus, h_wh155_base, h_wh155_JESMinus, h_wh155_JERPlus, h_wh155_JERMinus, h_wh155_TopPtPlus, h_wh155_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_wh155_bTagPlus, h_wh155_base, h_wh155_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_wh155_cTagPlus, h_wh155_base, h_wh155_cTagMinus, b)<<" $\\pm$ "<< sysUncBCTag (h_wh155_PileupPlus, h_wh155_base, h_wh155_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  
  outFile<<"HW, $M_{H}=160~GeV/c^{2}$"<< " & "<<h_wh160_base->GetBinContent(b)<<" $\\pm$ "<<h_wh160_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_wh160_JESPlus, h_wh160_base, h_wh160_JESMinus, h_wh160_JERPlus, h_wh160_JERMinus, h_wh160_TopPtPlus, h_wh160_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_wh160_bTagPlus, h_wh160_base, h_wh160_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_wh160_cTagPlus, h_wh160_base, h_wh160_cTagMinus, b)<<" $\\pm$ "<< sysUncBCTag (h_wh160_PileupPlus, h_wh160_base, h_wh160_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  
  
  //===================================
  outFile<<"SM $t\\bar{t}$"<<" & "<<h_ttbar_base->GetBinContent(b)<<" $\\pm$ "<<h_ttbar_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_ttbar_JESPlus, h_ttbar_base, h_ttbar_JESMinus, h_ttbar_JERPlus, h_ttbar_JERMinus, h_ttbar_TopPtPlus, h_ttbar_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_ttbar_bTagPlus, h_ttbar_base, h_ttbar_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_ttbar_cTagPlus, h_ttbar_base, h_ttbar_cTagMinus, b)<<" $\\pm$ "<<sysUncBCTag (h_ttbar_PileupPlus, h_ttbar_base, h_ttbar_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"W+Jets"<< " & "<<h_wjet_base->GetBinContent(b)<<" $\\pm$ "<<h_wjet_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_wjet_JESPlus, h_wjet_base, h_wjet_JESMinus, h_wjet_JERPlus, h_wjet_JERMinus, h_wjet_TopPtPlus, h_wjet_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_wjet_bTagPlus, h_wjet_base, h_wjet_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_wjet_cTagPlus, h_wjet_base, h_wjet_cTagMinus, b)<<" $\\pm$ "<<sysUncBCTag (h_wjet_PileupPlus, h_wjet_base, h_wjet_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  outFile<<"Z+Jets"<< " & "<<h_zjet_base->GetBinContent(b)<<" $\\pm$ "<<h_zjet_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_zjet_JESPlus, h_zjet_base, h_zjet_JESMinus, h_zjet_JERPlus, h_zjet_JERMinus, h_zjet_TopPtPlus, h_zjet_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_zjet_bTagPlus, h_zjet_base, h_zjet_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_zjet_cTagPlus, h_zjet_base, h_zjet_cTagMinus, b)<<" $\\pm$ "<<sysUncBCTag (h_zjet_PileupPlus, h_zjet_base, h_zjet_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"QCD"<< " & "<<h_qcd_base->GetBinContent(b)<<" $\\pm$ "<<h_qcd_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_qcd_JESPlus, h_qcd_base, h_qcd_JESMinus, h_qcd_JERPlus, h_qcd_JERMinus, h_qcd_TopPtPlus, h_qcd_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_qcd_bTagPlus, h_qcd_base, h_qcd_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_qcd_cTagPlus, h_qcd_base, h_qcd_cTagMinus, b)<<" $\\pm$ "<<sysUncBCTag (h_qcd_PileupPlus, h_qcd_base, h_qcd_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"SingleTop"<< " & "<<h_stop_base->GetBinContent(b)<<" $\\pm$ "<<h_stop_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_stop_JESPlus, h_stop_base, h_stop_JESMinus, h_stop_JERPlus, h_stop_JERMinus, h_stop_TopPtPlus, h_stop_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_stop_bTagPlus, h_stop_base, h_stop_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_stop_cTagPlus, h_stop_base, h_stop_cTagMinus, b)<<" $\\pm$ "<<sysUncBCTag (h_stop_PileupPlus, h_stop_base, h_stop_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"VV"<< " & "<<h_diboson_base->GetBinContent(b)<<" $\\pm$ "<<h_diboson_base->GetBinError(b)<<" $\\pm$ "<< sysUncJESTopPt( h_diboson_JESPlus, h_diboson_base, h_diboson_JESMinus, h_diboson_JERPlus, h_diboson_JERMinus, h_diboson_TopPtPlus, h_diboson_TopPtMinus, b) <<" $\\pm$ "<< sysUncBCTag (h_diboson_bTagPlus, h_diboson_base, h_diboson_bTagMinus, b) << " $\\pm$ "<< sysUncBCTag (h_diboson_cTagPlus, h_diboson_base, h_diboson_cTagMinus, b)<<" $\\pm$ "<<sysUncBCTag (h_diboson_PileupPlus, h_diboson_base, h_diboson_PileupMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"Total Bkg"<<" & "<<h_TotalBkg->GetBinContent(b)<<" $\\pm$ "<<h_TotalBkg->GetBinError(b)<<" $\\pm$ "<<" -- "<<" $\\pm$ "<<" -- "<<"$\\pm$ "<<" -- "<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;   
  outFile<<"\\hline "<<endl;  
  outFile<<" Data "<<" & "<<h_data->GetBinContent(b)<<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;    
  outFile<<"\\hline "<<endl;   
  outFile<<"\\end{tabular}"<<endl; 
  //outFile<<"\\end{LARGE}"<<endl;  
  outFile<<"\\end{center}"<<endl;
  outFile<<"\\caption{First table}"<<endl;
  outFile<<"\\end{table}"<<endl;

  //Add another table with JESUP
  double Nbin = 15;
 
  makeCutFlowTableAny(outFile, "base", Nbin, 		h_wh120_base,  	h_ttbar_base, h_stop_base, h_wjet_base, h_zjet_base, h_qcd_base, h_stop_base, h_diboson_base, h_data);
  makeCutFlowTableAny(outFile, "JESPlus", Nbin, 	h_wh120_JESPlus, 	h_ttbar_JESPlus, h_stop_JESPlus, h_wjet_JESPlus, h_zjet_JESPlus, h_qcd_JESPlus, h_stop_JESPlus, h_diboson_JESPlus, h_data);
  makeCutFlowTableAny(outFile, "JESMinus", Nbin, 	h_wh120_JESMinus, 	h_ttbar_JESMinus,h_stop_JERMinus, h_wjet_JESMinus, h_zjet_JESMinus, h_qcd_JESMinus, h_stop_JESMinus, h_diboson_JESMinus, h_data);
  makeCutFlowTableAny(outFile, "JERPlus", Nbin, 	h_wh120_JERPlus, 	h_ttbar_JERPlus, h_stop_JERPlus, h_wjet_JERPlus, h_zjet_JERPlus, h_qcd_JERPlus, h_stop_JERPlus, h_diboson_JERPlus, h_data);
  makeCutFlowTableAny(outFile, "JERMinus", Nbin, 	h_wh120_JERMinus, 	h_ttbar_JERMinus, h_stop_JERMinus, h_wjet_JERMinus, h_zjet_JERMinus, h_qcd_JERMinus, h_stop_JERMinus, h_diboson_JERMinus, h_data);
  makeCutFlowTableAny(outFile, "TopPtPlus", Nbin, 	h_wh120_TopPtPlus, h_ttbar_TopPtPlus, h_stop_TopPtPlus, h_wjet_TopPtPlus, h_zjet_TopPtPlus, h_qcd_TopPtPlus, h_stop_TopPtPlus, h_diboson_TopPtPlus, h_data);
  makeCutFlowTableAny(outFile, "TopPtMinus", Nbin, 	h_wh120_TopPtMinus, h_ttbar_TopPtMinus, h_stop_TopPtMinus, h_wjet_TopPtMinus, h_zjet_TopPtMinus, h_qcd_TopPtMinus, h_stop_TopPtMinus, h_diboson_TopPtMinus, h_data);
  makeCutFlowTableAny(outFile, "bTagPlus", Nbin, 	h_wh120_bTagPlus, 	h_ttbar_bTagPlus, h_stop_bTagPlus, h_wjet_bTagPlus, h_zjet_bTagPlus, h_qcd_bTagPlus, h_stop_bTagPlus, h_diboson_bTagPlus, h_data);
  makeCutFlowTableAny(outFile, "bTagMinus", Nbin, 	h_wh120_bTagMinus, h_ttbar_bTagMinus, h_stop_bTagMinus, h_wjet_bTagMinus, h_zjet_bTagMinus, h_qcd_bTagMinus, h_stop_bTagMinus, h_diboson_bTagMinus, h_data);
  makeCutFlowTableAny(outFile, "cTagPlus", Nbin, 	h_wh120_cTagPlus, 	h_ttbar_cTagPlus, h_stop_cTagPlus, h_wjet_cTagPlus, h_zjet_cTagPlus, h_qcd_cTagPlus, h_stop_cTagPlus, h_diboson_cTagPlus, h_data);
  makeCutFlowTableAny(outFile, "cTagMinus", Nbin, 	h_wh120_cTagMinus, h_ttbar_cTagMinus, h_stop_cTagMinus, h_wjet_cTagMinus, h_zjet_cTagMinus, h_qcd_cTagMinus, h_stop_cTagMinus, h_diboson_cTagMinus, h_data);
  makeCutFlowTableAny(outFile, "PileupPlus", Nbin, 	h_wh120_PileupPlus, 	h_ttbar_PileupPlus, h_stop_PileupPlus, h_wjet_PileupPlus, h_zjet_PileupPlus, h_qcd_PileupPlus, h_stop_PileupPlus, h_diboson_PileupPlus, h_data);
  makeCutFlowTableAny(outFile, "PileupMinus", Nbin, 	h_wh120_PileupMinus, h_ttbar_PileupMinus, h_stop_PileupMinus, h_wjet_PileupMinus, h_zjet_PileupMinus, h_qcd_PileupMinus, h_stop_PileupMinus, h_diboson_PileupMinus, h_data);
  outFile<<"\\end{document}"<<endl;  
  outFile.close(); 
} 

void makeFinalSysUncTable(TString histSubDir="Iso/", TString histName="cutflow"){  
  
  TString inFile("$PWD/");
  TString histPath(histSubDir+histName);
  TFile *wh80 	  		= new TFile(inFile+"all_Hplus80.root"); 
  TFile *wh90 	  		= new TFile(inFile+"all_Hplus90.root"); 
  TFile *wh100 	  		= new TFile(inFile+"all_Hplus100.root"); 
  TFile *wh120 	  		= new TFile(inFile+"all_Hplus120.root"); 
  TFile *wh140 	  		= new TFile(inFile+"all_Hplus140.root"); 
  TFile *wh150 	  		= new TFile(inFile+"all_Hplus150.root"); 
  TFile *wh155 	  		= new TFile(inFile+"all_Hplus155.root"); 
  TFile *wh160 	  		= new TFile(inFile+"all_Hplus160.root"); 
  
  TFile *ttbar    		= new TFile(inFile+"all_TTJetsP.root"); 
  TFile *wjet  			= new TFile(inFile+"all_WJets.root"); 
  TFile *zjet  			= new TFile(inFile+"all_DY.root");
  TFile *qcd  			= new TFile(inFile+"all_QCD.root");
  TFile *stop  			= new TFile(inFile+"all_ST.root");
  TFile *diboson 		= new TFile(inFile+"all_VV.root");
  
  //Hplus  wh80   signal 
  TH1F* h_wh80_base  		= (TH1F*)wh80->Get("base/"+histPath)->Clone("h_wh80_base"); 
  TH1F* h_wh80_JESPlus 		= (TH1F*)wh80->Get("JESPlus/"+histPath)->Clone("h_wh80_JESPlus");
  TH1F* h_wh80_JESMinus 		= (TH1F*)wh80->Get("JESMinus/"+histPath)->Clone("h_wh80_JESMinus");
  TH1F* h_wh80_PileupPlus 		= (TH1F*)wh80->Get("PileupPlus/"+histPath)->Clone("h_wh80_PileupPlus");
  TH1F* h_wh80_PileupMinus 		= (TH1F*)wh80->Get("PileupMinus/"+histPath)->Clone("h_wh80_PileupMinus");
  TH1F* h_wh80_JERPlus 		= (TH1F*)wh80->Get("JERPlus/"+histPath)->Clone("h_wh80_JERPlus");
  TH1F* h_wh80_JERMinus 		= (TH1F*)wh80->Get("JERMinus/"+histPath)->Clone("h_wh80_JERMinus");
  TH1F* h_wh80_TopPtPlus 		= (TH1F*)wh80->Get("TopPtPlus/"+histPath)->Clone("h_wh80_TopPtPlus");
  TH1F* h_wh80_TopPtMinus 	= (TH1F*)wh80->Get("TopPtMinus/"+histPath)->Clone("h_wh80_TopPtMinus");
  TH1F* h_wh80_bTagPlus 		= (TH1F*)wh80->Get("bTagPlus/"+histPath)->Clone("h_wh80_bTagPlus");
  TH1F* h_wh80_bTagMinus 		= (TH1F*)wh80->Get("bTagMinus/"+histPath)->Clone("h_wh80_bTagMinus");
  TH1F* h_wh80_cTagPlus 		= (TH1F*)wh80->Get("cTagPlus/"+histPath)->Clone("h_wh80_cTagPlus");
  TH1F* h_wh80_cTagMinus 		= (TH1F*)wh80->Get("cTagMinus/"+histPath)->Clone("h_wh80_cTagMinus");
  
  //Hplus  wh90   signal 
  TH1F* h_wh90_base  		= (TH1F*)wh90->Get("base/"+histPath)->Clone("h_wh90_base"); 
  TH1F* h_wh90_JESPlus 		= (TH1F*)wh90->Get("JESPlus/"+histPath)->Clone("h_wh90_JESPlus");
  TH1F* h_wh90_JESMinus 		= (TH1F*)wh90->Get("JESMinus/"+histPath)->Clone("h_wh90_JESMinus");
  TH1F* h_wh90_PileupPlus 		= (TH1F*)wh90->Get("PileupPlus/"+histPath)->Clone("h_wh90_PileupPlus");
  TH1F* h_wh90_PileupMinus 		= (TH1F*)wh90->Get("PileupMinus/"+histPath)->Clone("h_wh90_PileupMinus");
  TH1F* h_wh90_JERPlus 		= (TH1F*)wh90->Get("JERPlus/"+histPath)->Clone("h_wh90_JERPlus");
  TH1F* h_wh90_JERMinus 		= (TH1F*)wh90->Get("JERMinus/"+histPath)->Clone("h_wh90_JERMinus");
  TH1F* h_wh90_TopPtPlus 		= (TH1F*)wh90->Get("TopPtPlus/"+histPath)->Clone("h_wh90_TopPtPlus");
  TH1F* h_wh90_TopPtMinus 	= (TH1F*)wh90->Get("TopPtMinus/"+histPath)->Clone("h_wh90_TopPtMinus");
  TH1F* h_wh90_bTagPlus 		= (TH1F*)wh90->Get("bTagPlus/"+histPath)->Clone("h_wh90_bTagPlus");
  TH1F* h_wh90_bTagMinus 		= (TH1F*)wh90->Get("bTagMinus/"+histPath)->Clone("h_wh90_bTagMinus");
  TH1F* h_wh90_cTagPlus 		= (TH1F*)wh90->Get("cTagPlus/"+histPath)->Clone("h_wh90_cTagPlus");
  TH1F* h_wh90_cTagMinus 		= (TH1F*)wh90->Get("cTagMinus/"+histPath)->Clone("h_wh90_cTagMinus");
  
  
  //Hplus  wh100   signal 
  TH1F* h_wh100_base  		= (TH1F*)wh100->Get("base/"+histPath)->Clone("h_wh100_base"); 
  TH1F* h_wh100_JESPlus 		= (TH1F*)wh100->Get("JESPlus/"+histPath)->Clone("h_wh100_JESPlus");
  TH1F* h_wh100_JESMinus 		= (TH1F*)wh100->Get("JESMinus/"+histPath)->Clone("h_wh100_JESMinus");
  TH1F* h_wh100_PileupPlus 		= (TH1F*)wh100->Get("PileupPlus/"+histPath)->Clone("h_wh100_PileupPlus");
  TH1F* h_wh100_PileupMinus 		= (TH1F*)wh100->Get("PileupMinus/"+histPath)->Clone("h_wh100_PileupMinus");
  TH1F* h_wh100_JERPlus 		= (TH1F*)wh100->Get("JERPlus/"+histPath)->Clone("h_wh100_JERPlus");
  TH1F* h_wh100_JERMinus 		= (TH1F*)wh100->Get("JERMinus/"+histPath)->Clone("h_wh100_JERMinus");
  TH1F* h_wh100_TopPtPlus 		= (TH1F*)wh100->Get("TopPtPlus/"+histPath)->Clone("h_wh100_TopPtPlus");
  TH1F* h_wh100_TopPtMinus 	= (TH1F*)wh100->Get("TopPtMinus/"+histPath)->Clone("h_wh100_TopPtMinus");
  TH1F* h_wh100_bTagPlus 		= (TH1F*)wh100->Get("bTagPlus/"+histPath)->Clone("h_wh100_bTagPlus");
  TH1F* h_wh100_bTagMinus 		= (TH1F*)wh100->Get("bTagMinus/"+histPath)->Clone("h_wh100_bTagMinus");
  TH1F* h_wh100_cTagPlus 		= (TH1F*)wh100->Get("cTagPlus/"+histPath)->Clone("h_wh100_cTagPlus");
  TH1F* h_wh100_cTagMinus 		= (TH1F*)wh100->Get("cTagMinus/"+histPath)->Clone("h_wh100_cTagMinus");
  
  //Hplus  wh120   signal 
  TH1F* h_wh120_base  		= (TH1F*)wh120->Get("base/"+histPath)->Clone("h_wh120_base"); 
  TH1F* h_wh120_JESPlus 		= (TH1F*)wh120->Get("JESPlus/"+histPath)->Clone("h_wh120_JESPlus");
  TH1F* h_wh120_JESMinus 		= (TH1F*)wh120->Get("JESMinus/"+histPath)->Clone("h_wh120_JESMinus");
  TH1F* h_wh120_PileupPlus 		= (TH1F*)wh120->Get("PileupPlus/"+histPath)->Clone("h_wh120_PileupPlus");
  TH1F* h_wh120_PileupMinus 		= (TH1F*)wh120->Get("PileupMinus/"+histPath)->Clone("h_wh120_PileupMinus");
  TH1F* h_wh120_JERPlus 		= (TH1F*)wh120->Get("JERPlus/"+histPath)->Clone("h_wh120_JERPlus");
  TH1F* h_wh120_JERMinus 		= (TH1F*)wh120->Get("JERMinus/"+histPath)->Clone("h_wh120_JERMinus");
  TH1F* h_wh120_TopPtPlus 		= (TH1F*)wh120->Get("TopPtPlus/"+histPath)->Clone("h_wh120_TopPtPlus");
  TH1F* h_wh120_TopPtMinus 	= (TH1F*)wh120->Get("TopPtMinus/"+histPath)->Clone("h_wh120_TopPtMinus");
  TH1F* h_wh120_bTagPlus 		= (TH1F*)wh120->Get("bTagPlus/"+histPath)->Clone("h_wh120_bTagPlus");
  TH1F* h_wh120_bTagMinus 		= (TH1F*)wh120->Get("bTagMinus/"+histPath)->Clone("h_wh120_bTagMinus");
  TH1F* h_wh120_cTagPlus 		= (TH1F*)wh120->Get("cTagPlus/"+histPath)->Clone("h_wh120_cTagPlus");
  TH1F* h_wh120_cTagMinus 		= (TH1F*)wh120->Get("cTagMinus/"+histPath)->Clone("h_wh120_cTagMinus");
  
  //Hplus  wh140   signal 
  TH1F* h_wh140_base  		= (TH1F*)wh140->Get("base/"+histPath)->Clone("h_wh140_base"); 
  TH1F* h_wh140_JESPlus 		= (TH1F*)wh140->Get("JESPlus/"+histPath)->Clone("h_wh140_JESPlus");
  TH1F* h_wh140_JESMinus 		= (TH1F*)wh140->Get("JESMinus/"+histPath)->Clone("h_wh140_JESMinus");
  TH1F* h_wh140_PileupPlus 		= (TH1F*)wh140->Get("PileupPlus/"+histPath)->Clone("h_wh140_PileupPlus");
  TH1F* h_wh140_PileupMinus 		= (TH1F*)wh140->Get("PileupMinus/"+histPath)->Clone("h_wh140_PileupMinus");
  TH1F* h_wh140_JERPlus 		= (TH1F*)wh140->Get("JERPlus/"+histPath)->Clone("h_wh140_JERPlus");
  TH1F* h_wh140_JERMinus 		= (TH1F*)wh140->Get("JERMinus/"+histPath)->Clone("h_wh140_JERMinus");
  TH1F* h_wh140_TopPtPlus 		= (TH1F*)wh140->Get("TopPtPlus/"+histPath)->Clone("h_wh140_TopPtPlus");
  TH1F* h_wh140_TopPtMinus 	= (TH1F*)wh140->Get("TopPtMinus/"+histPath)->Clone("h_wh140_TopPtMinus");
  TH1F* h_wh140_bTagPlus 		= (TH1F*)wh140->Get("bTagPlus/"+histPath)->Clone("h_wh140_bTagPlus");
  TH1F* h_wh140_bTagMinus 		= (TH1F*)wh140->Get("bTagMinus/"+histPath)->Clone("h_wh140_bTagMinus");
  TH1F* h_wh140_cTagPlus 		= (TH1F*)wh140->Get("cTagPlus/"+histPath)->Clone("h_wh140_cTagPlus");
  TH1F* h_wh140_cTagMinus 		= (TH1F*)wh140->Get("cTagMinus/"+histPath)->Clone("h_wh140_cTagMinus");
  
  
  //Hplus  wh150   signal 
  TH1F* h_wh150_base  		= (TH1F*)wh150->Get("base/"+histPath)->Clone("h_wh150_base"); 
  TH1F* h_wh150_JESPlus 		= (TH1F*)wh150->Get("JESPlus/"+histPath)->Clone("h_wh150_JESPlus");
  TH1F* h_wh150_JESMinus 		= (TH1F*)wh150->Get("JESMinus/"+histPath)->Clone("h_wh150_JESMinus");
  TH1F* h_wh150_PileupPlus 		= (TH1F*)wh150->Get("PileupPlus/"+histPath)->Clone("h_wh150_PileupPlus");
  TH1F* h_wh150_PileupMinus 		= (TH1F*)wh150->Get("PileupMinus/"+histPath)->Clone("h_wh150_PileupMinus");
  TH1F* h_wh150_JERPlus 		= (TH1F*)wh150->Get("JERPlus/"+histPath)->Clone("h_wh150_JERPlus");
  TH1F* h_wh150_JERMinus 		= (TH1F*)wh150->Get("JERMinus/"+histPath)->Clone("h_wh150_JERMinus");
  TH1F* h_wh150_TopPtPlus 		= (TH1F*)wh150->Get("TopPtPlus/"+histPath)->Clone("h_wh150_TopPtPlus");
  TH1F* h_wh150_TopPtMinus 	= (TH1F*)wh150->Get("TopPtMinus/"+histPath)->Clone("h_wh150_TopPtMinus");
  TH1F* h_wh150_bTagPlus 		= (TH1F*)wh150->Get("bTagPlus/"+histPath)->Clone("h_wh150_bTagPlus");
  TH1F* h_wh150_bTagMinus 		= (TH1F*)wh150->Get("bTagMinus/"+histPath)->Clone("h_wh150_bTagMinus");
  TH1F* h_wh150_cTagPlus 		= (TH1F*)wh150->Get("cTagPlus/"+histPath)->Clone("h_wh150_cTagPlus");
  TH1F* h_wh150_cTagMinus 		= (TH1F*)wh150->Get("cTagMinus/"+histPath)->Clone("h_wh150_cTagMinus");
  
  //Hplus  wh155   signal 
  TH1F* h_wh155_base  		= (TH1F*)wh155->Get("base/"+histPath)->Clone("h_wh155_base"); 
  TH1F* h_wh155_JESPlus 		= (TH1F*)wh155->Get("JESPlus/"+histPath)->Clone("h_wh155_JESPlus");
  TH1F* h_wh155_JESMinus 		= (TH1F*)wh155->Get("JESMinus/"+histPath)->Clone("h_wh155_JESMinus");
  TH1F* h_wh155_PileupPlus 		= (TH1F*)wh155->Get("PileupPlus/"+histPath)->Clone("h_wh155_PileupPlus");
  TH1F* h_wh155_PileupMinus 		= (TH1F*)wh155->Get("PileupMinus/"+histPath)->Clone("h_wh155_PileupMinus");
  TH1F* h_wh155_JERPlus 		= (TH1F*)wh155->Get("JERPlus/"+histPath)->Clone("h_wh155_JERPlus");
  TH1F* h_wh155_JERMinus 		= (TH1F*)wh155->Get("JERMinus/"+histPath)->Clone("h_wh155_JERMinus");
  TH1F* h_wh155_TopPtPlus 		= (TH1F*)wh155->Get("TopPtPlus/"+histPath)->Clone("h_wh155_TopPtPlus");
  TH1F* h_wh155_TopPtMinus 	= (TH1F*)wh155->Get("TopPtMinus/"+histPath)->Clone("h_wh155_TopPtMinus");
  TH1F* h_wh155_bTagPlus 		= (TH1F*)wh155->Get("bTagPlus/"+histPath)->Clone("h_wh155_bTagPlus");
  TH1F* h_wh155_bTagMinus 		= (TH1F*)wh155->Get("bTagMinus/"+histPath)->Clone("h_wh155_bTagMinus");
  TH1F* h_wh155_cTagPlus 		= (TH1F*)wh155->Get("cTagPlus/"+histPath)->Clone("h_wh155_cTagPlus");
  TH1F* h_wh155_cTagMinus 		= (TH1F*)wh155->Get("cTagMinus/"+histPath)->Clone("h_wh155_cTagMinus");
  
  //Hplus  wh160   signal 
  TH1F* h_wh160_base  		= (TH1F*)wh160->Get("base/"+histPath)->Clone("h_wh160_base"); 
  TH1F* h_wh160_JESPlus 		= (TH1F*)wh160->Get("JESPlus/"+histPath)->Clone("h_wh160_JESPlus");
  TH1F* h_wh160_JESMinus 		= (TH1F*)wh160->Get("JESMinus/"+histPath)->Clone("h_wh160_JESMinus");
  TH1F* h_wh160_PileupPlus 		= (TH1F*)wh160->Get("PileupPlus/"+histPath)->Clone("h_wh160_PileupPlus");
  TH1F* h_wh160_PileupMinus 		= (TH1F*)wh160->Get("PileupMinus/"+histPath)->Clone("h_wh160_PileupMinus");
  TH1F* h_wh160_JERPlus 		= (TH1F*)wh160->Get("JERPlus/"+histPath)->Clone("h_wh160_JERPlus");
  TH1F* h_wh160_JERMinus 		= (TH1F*)wh160->Get("JERMinus/"+histPath)->Clone("h_wh160_JERMinus");
  TH1F* h_wh160_TopPtPlus 		= (TH1F*)wh160->Get("TopPtPlus/"+histPath)->Clone("h_wh160_TopPtPlus");
  TH1F* h_wh160_TopPtMinus 	= (TH1F*)wh160->Get("TopPtMinus/"+histPath)->Clone("h_wh160_TopPtMinus");
  TH1F* h_wh160_bTagPlus 		= (TH1F*)wh160->Get("bTagPlus/"+histPath)->Clone("h_wh160_bTagPlus");
  TH1F* h_wh160_bTagMinus 		= (TH1F*)wh160->Get("bTagMinus/"+histPath)->Clone("h_wh160_bTagMinus");
  TH1F* h_wh160_cTagPlus 		= (TH1F*)wh160->Get("cTagPlus/"+histPath)->Clone("h_wh160_cTagPlus");
  TH1F* h_wh160_cTagMinus 		= (TH1F*)wh160->Get("cTagMinus/"+histPath)->Clone("h_wh160_cTagMinus");
  
  
  //ttbar+jets
  TH1F* h_ttbar_base  		= (TH1F*)ttbar->Get("base/"+histPath)->Clone("h_ttbar_base");  
  TH1F* h_ttbar_JESPlus 	= (TH1F*)ttbar->Get("JESPlus/"+histPath)->Clone("h_ttbar_JESPlus"); 
  TH1F* h_ttbar_JESMinus 	= (TH1F*)ttbar->Get("JESMinus/"+histPath)->Clone("h_ttbar_JESMinus"); 
  TH1F* h_ttbar_PileupPlus 	= (TH1F*)ttbar->Get("PileupPlus/"+histPath)->Clone("h_ttbar_PileupPlus"); 
  TH1F* h_ttbar_PileupMinus 	= (TH1F*)ttbar->Get("PileupMinus/"+histPath)->Clone("h_ttbar_PileupMinus"); 
  TH1F* h_ttbar_JERPlus 	= (TH1F*)ttbar->Get("JERPlus/"+histPath)->Clone("h_ttbar_JERPlus"); 
  TH1F* h_ttbar_JERMinus 	= (TH1F*)ttbar->Get("JERMinus/"+histPath)->Clone("h_ttbar_JERMinus"); 
  TH1F* h_ttbar_TopPtPlus 	= (TH1F*)ttbar->Get("TopPtPlus/"+histPath)->Clone("h_ttbar_TopPtPlus"); 
  TH1F* h_ttbar_TopPtMinus 	= (TH1F*)ttbar->Get("TopPtMinus/"+histPath)->Clone("h_ttbar_TopPtMinus"); 
  TH1F* h_ttbar_bTagPlus 	= (TH1F*)ttbar->Get("bTagPlus/"+histPath)->Clone("h_ttbar_bTagPlus"); 
  TH1F* h_ttbar_bTagMinus 	= (TH1F*)ttbar->Get("bTagMinus/"+histPath)->Clone("h_ttbar_bTagMinus"); 
  TH1F* h_ttbar_cTagPlus 	= (TH1F*)ttbar->Get("cTagPlus/"+histPath)->Clone("h_ttbar_cTagPlus"); 
  TH1F* h_ttbar_cTagMinus 	= (TH1F*)ttbar->Get("cTagMinus/"+histPath)->Clone("h_ttbar_cTagMinus"); 
  //w+jets
  TH1F* h_wjet_base  		= (TH1F*)wjet->Get("base/"+histPath)->Clone("h_wjet_base");  
  TH1F* h_wjet_JESPlus 		= (TH1F*)wjet->Get("JESPlus/"+histPath)->Clone("h_wjet_JESPlus"); 
  TH1F* h_wjet_JESMinus 	= (TH1F*)wjet->Get("JESMinus/"+histPath)->Clone("h_wjet_JESMinus"); 
  TH1F* h_wjet_PileupPlus 		= (TH1F*)wjet->Get("PileupPlus/"+histPath)->Clone("h_wjet_PileupPlus"); 
  TH1F* h_wjet_PileupMinus 	= (TH1F*)wjet->Get("PileupMinus/"+histPath)->Clone("h_wjet_PileupMinus"); 
  TH1F* h_wjet_JERPlus 		= (TH1F*)wjet->Get("JERPlus/"+histPath)->Clone("h_wjet_JERPlus"); 
  TH1F* h_wjet_JERMinus 	= (TH1F*)wjet->Get("JERMinus/"+histPath)->Clone("h_wjet_JERMinus"); 
  TH1F* h_wjet_TopPtPlus 	= (TH1F*)wjet->Get("TopPtPlus/"+histPath)->Clone("h_wjet_TopPtPlus"); 
  TH1F* h_wjet_TopPtMinus 	= (TH1F*)wjet->Get("TopPtMinus/"+histPath)->Clone("h_wjet_TopPtMinus"); 
  TH1F* h_wjet_bTagPlus 	= (TH1F*)wjet->Get("bTagPlus/"+histPath)->Clone("h_wjet_bTagPlus"); 
  TH1F* h_wjet_bTagMinus 	= (TH1F*)wjet->Get("bTagMinus/"+histPath)->Clone("h_wjet_bTagMinus"); 
  TH1F* h_wjet_cTagPlus 	= (TH1F*)wjet->Get("cTagPlus/"+histPath)->Clone("h_wjet_cTagPlus"); 
  TH1F* h_wjet_cTagMinus 	= (TH1F*)wjet->Get("cTagMinus/"+histPath)->Clone("h_wjet_cTagMinus"); 
  //dy+jets
  TH1F* h_zjet_base  		= (TH1F*)zjet->Get("base/"+histPath)->Clone("h_zjet_base");
  TH1F* h_zjet_JESPlus 		= (TH1F*)zjet->Get("JESPlus/"+histPath)->Clone("h_zjet_JESPlus");
  TH1F* h_zjet_JESMinus 	= (TH1F*)zjet->Get("JESMinus/"+histPath)->Clone("h_zjet_JESMinus");
  TH1F* h_zjet_PileupPlus 		= (TH1F*)zjet->Get("PileupPlus/"+histPath)->Clone("h_zjet_PileupPlus");
  TH1F* h_zjet_PileupMinus 	= (TH1F*)zjet->Get("PileupMinus/"+histPath)->Clone("h_zjet_PileupMinus");
  TH1F* h_zjet_JERPlus 		= (TH1F*)zjet->Get("JERPlus/"+histPath)->Clone("h_zjet_JERPlus");
  TH1F* h_zjet_JERMinus 	= (TH1F*)zjet->Get("JERMinus/"+histPath)->Clone("h_zjet_JERMinus");
  TH1F* h_zjet_TopPtPlus 	= (TH1F*)zjet->Get("TopPtPlus/"+histPath)->Clone("h_zjet_TopPtPlus");
  TH1F* h_zjet_TopPtMinus 	= (TH1F*)zjet->Get("TopPtMinus/"+histPath)->Clone("h_zjet_TopPtMinus");
  TH1F* h_zjet_bTagPlus 	= (TH1F*)zjet->Get("bTagPlus/"+histPath)->Clone("h_zjet_bTagPlus");
  TH1F* h_zjet_bTagMinus 	= (TH1F*)zjet->Get("bTagMinus/"+histPath)->Clone("h_zjet_bTagMinus");
  TH1F* h_zjet_cTagPlus 	= (TH1F*)zjet->Get("cTagPlus/"+histPath)->Clone("h_zjet_cTagPlus");
  TH1F* h_zjet_cTagMinus 	= (TH1F*)zjet->Get("cTagMinus/"+histPath)->Clone("h_zjet_cTagMinus");
  //qcd
  TH1F* h_qcd_base 		= (TH1F*)qcd->Get("base/"+histPath)->Clone("h_qcd_base");
  TH1F* h_qcd_JESPlus 		= (TH1F*)qcd->Get("JESPlus/"+histPath)->Clone("h_qcd_JESPlus");
  TH1F* h_qcd_JESMinus 		= (TH1F*)qcd->Get("JESMinus/"+histPath)->Clone("h_qcd_JESMinus");
  TH1F* h_qcd_PileupPlus 		= (TH1F*)qcd->Get("PileupPlus/"+histPath)->Clone("h_qcd_PileupPlus");
  TH1F* h_qcd_PileupMinus 		= (TH1F*)qcd->Get("PileupMinus/"+histPath)->Clone("h_qcd_PileupMinus");
  TH1F* h_qcd_JERPlus 		= (TH1F*)qcd->Get("JERPlus/"+histPath)->Clone("h_qcd_JERPlus");
  TH1F* h_qcd_JERMinus 		= (TH1F*)qcd->Get("JERMinus/"+histPath)->Clone("h_qcd_JERMinus");
  TH1F* h_qcd_TopPtPlus 	= (TH1F*)qcd->Get("TopPtPlus/"+histPath)->Clone("h_qcd_TopPtPlus");
  TH1F* h_qcd_TopPtMinus 	= (TH1F*)qcd->Get("TopPtMinus/"+histPath)->Clone("h_qcd_TopPtMinus");
  TH1F* h_qcd_bTagPlus 		= (TH1F*)qcd->Get("bTagPlus/"+histPath)->Clone("h_qcd_bTagPlus");
  TH1F* h_qcd_bTagMinus 	= (TH1F*)qcd->Get("bTagMinus/"+histPath)->Clone("h_qcd_bTagMinus");
  TH1F* h_qcd_cTagPlus 		= (TH1F*)qcd->Get("cTagPlus/"+histPath)->Clone("h_qcd_cTagPlus");
  TH1F* h_qcd_cTagMinus 	= (TH1F*)qcd->Get("cTagMinus/"+histPath)->Clone("h_qcd_cTagMinus");
  //single top
  TH1F* h_stop_base  		= (TH1F*)stop->Get("base/"+histPath)->Clone("h_stop_base");
  TH1F* h_stop_JESPlus 		= (TH1F*)stop->Get("JESPlus/"+histPath)->Clone("h_stop_JESPlus");
  TH1F* h_stop_JESMinus 	= (TH1F*)stop->Get("JESMinus/"+histPath)->Clone("h_stop_JESMinus");
  TH1F* h_stop_PileupPlus 		= (TH1F*)stop->Get("PileupPlus/"+histPath)->Clone("h_stop_PileupPlus");
  TH1F* h_stop_PileupMinus 	= (TH1F*)stop->Get("PileupMinus/"+histPath)->Clone("h_stop_PileupMinus");
  TH1F* h_stop_JERPlus 		= (TH1F*)stop->Get("JERPlus/"+histPath)->Clone("h_stop_JERPlus");
  TH1F* h_stop_JERMinus 	= (TH1F*)stop->Get("JERMinus/"+histPath)->Clone("h_stop_JERMinus");
  TH1F* h_stop_TopPtPlus 	= (TH1F*)stop->Get("TopPtPlus/"+histPath)->Clone("h_stop_TopPtPlus");
  TH1F* h_stop_TopPtMinus 	= (TH1F*)stop->Get("TopPtMinus/"+histPath)->Clone("h_stop_TopPtMinus");
  TH1F* h_stop_bTagPlus 	= (TH1F*)stop->Get("bTagPlus/"+histPath)->Clone("h_stop_bTagPlus");
  TH1F* h_stop_bTagMinus 	= (TH1F*)stop->Get("bTagMinus/"+histPath)->Clone("h_stop_bTagMinus");
  TH1F* h_stop_cTagPlus 	= (TH1F*)stop->Get("cTagPlus/"+histPath)->Clone("h_stop_cTagPlus");
  TH1F* h_stop_cTagMinus 	= (TH1F*)stop->Get("cTagMinus/"+histPath)->Clone("h_stop_cTagMinus");
  //vv
  TH1F* h_diboson_base 		= (TH1F*)diboson->Get("base/"+histPath)->Clone("h_diboson_base");
  TH1F* h_diboson_JESPlus 	= (TH1F*)diboson->Get("JESPlus/"+histPath)->Clone("h_diboson_JESPlus");
  TH1F* h_diboson_JESMinus 	= (TH1F*)diboson->Get("JESMinus/"+histPath)->Clone("h_diboson_JESMinus");
  TH1F* h_diboson_PileupPlus 	= (TH1F*)diboson->Get("PileupPlus/"+histPath)->Clone("h_diboson_PileupPlus");
  TH1F* h_diboson_PileupMinus 	= (TH1F*)diboson->Get("PileupMinus/"+histPath)->Clone("h_diboson_PileupMinus");
  TH1F* h_diboson_JERPlus 	= (TH1F*)diboson->Get("JERPlus/"+histPath)->Clone("h_diboson_JERPlus");
  TH1F* h_diboson_JERMinus 	= (TH1F*)diboson->Get("JERMinus/"+histPath)->Clone("h_diboson_JERMinus");
  TH1F* h_diboson_TopPtPlus 	= (TH1F*)diboson->Get("TopPtPlus/"+histPath)->Clone("h_diboson_TopPtPlus");
  TH1F* h_diboson_TopPtMinus 	= (TH1F*)diboson->Get("TopPtMinus/"+histPath)->Clone("h_diboson_TopPtMinus");
  TH1F* h_diboson_bTagPlus 	= (TH1F*)diboson->Get("bTagPlus/"+histPath)->Clone("h_diboson_bTagPlus");
  TH1F* h_diboson_bTagMinus 	= (TH1F*)diboson->Get("bTagMinus/"+histPath)->Clone("h_diboson_bTagMinus");
  TH1F* h_diboson_cTagPlus 	= (TH1F*)diboson->Get("cTagPlus/"+histPath)->Clone("h_diboson_cTagPlus");
  TH1F* h_diboson_cTagMinus 	= (TH1F*)diboson->Get("cTagMinus/"+histPath)->Clone("h_diboson_cTagMinus");

  //get average of top pt-reweighting 
  TH1F* hsig_avgTopwh80 = (TH1F*)wh80->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh80 = hsig_avgTopwh80->GetMean();
  double sfwh80 = 1;
  h_wh80_base->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_JESPlus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_JESMinus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_PileupPlus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_PileupMinus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_JERPlus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_JERMinus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_TopPtPlus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_TopPtMinus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_bTagPlus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_bTagMinus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_cTagPlus->Scale(sfwh80/sig_avgTopwh80);
  h_wh80_cTagMinus->Scale(sfwh80/sig_avgTopwh80);

  TH1F* hsig_avgTopwh90 = (TH1F*)wh90->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh90 = hsig_avgTopwh90->GetMean();
  double sfwh90 = 1;
  h_wh90_base->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_JESPlus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_JESMinus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_PileupPlus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_PileupMinus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_JERPlus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_JERMinus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_TopPtPlus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_TopPtMinus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_bTagPlus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_bTagMinus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_cTagPlus->Scale(sfwh90/sig_avgTopwh90);
  h_wh90_cTagMinus->Scale(sfwh90/sig_avgTopwh90);

  TH1F* hsig_avgTopwh100 = (TH1F*)wh100->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh100 = hsig_avgTopwh100->GetMean();
  double sfwh100 = 1;
  h_wh100_base->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_JESPlus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_JESMinus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_PileupPlus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_PileupMinus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_JERPlus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_JERMinus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_TopPtPlus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_TopPtMinus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_bTagPlus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_bTagMinus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_cTagPlus->Scale(sfwh100/sig_avgTopwh100);
  h_wh100_cTagMinus->Scale(sfwh100/sig_avgTopwh100);

  TH1F* hsig_avgTopwh120 = (TH1F*)wh120->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh120 = hsig_avgTopwh120->GetMean();
  double sfwh120 = 1;
  h_wh120_base->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_JESPlus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_JESMinus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_PileupPlus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_PileupMinus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_JERPlus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_JERMinus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_TopPtPlus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_TopPtMinus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_bTagPlus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_bTagMinus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_cTagPlus->Scale(sfwh120/sig_avgTopwh120);
  h_wh120_cTagMinus->Scale(sfwh120/sig_avgTopwh120);

  TH1F* hsig_avgTopwh140 = (TH1F*)wh140->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh140 = hsig_avgTopwh140->GetMean();
  double sfwh140 = 1;
  h_wh140_base->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_JESPlus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_JESMinus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_PileupPlus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_PileupMinus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_JERPlus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_JERMinus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_TopPtPlus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_TopPtMinus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_bTagPlus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_bTagMinus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_cTagPlus->Scale(sfwh140/sig_avgTopwh140);
  h_wh140_cTagMinus->Scale(sfwh140/sig_avgTopwh140);

  TH1F* hsig_avgTopwh150 = (TH1F*)wh150->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh150 = hsig_avgTopwh150->GetMean();
  double sfwh150 = 1;
  h_wh150_base->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_JESPlus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_JESMinus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_PileupPlus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_PileupMinus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_JERPlus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_JERMinus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_TopPtPlus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_TopPtMinus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_bTagPlus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_bTagMinus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_cTagPlus->Scale(sfwh150/sig_avgTopwh150);
  h_wh150_cTagMinus->Scale(sfwh150/sig_avgTopwh150);

  TH1F* hsig_avgTopwh155 = (TH1F*)wh155->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh155 = hsig_avgTopwh155->GetMean();
  double sfwh155 = 1;
  h_wh155_base->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_JESPlus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_JESMinus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_PileupPlus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_PileupMinus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_JERPlus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_JERMinus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_TopPtPlus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_TopPtMinus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_bTagPlus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_bTagMinus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_cTagPlus->Scale(sfwh155/sig_avgTopwh155);
  h_wh155_cTagMinus->Scale(sfwh155/sig_avgTopwh155);

  TH1F* hsig_avgTopwh160 = (TH1F*)wh160->Get("base/SF_topPtWeights"); 
  double sig_avgTopwh160 = hsig_avgTopwh160->GetMean();
  double sfwh160 = 1;
  h_wh160_base->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_JESPlus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_JESMinus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_PileupPlus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_PileupMinus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_JERPlus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_JERMinus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_TopPtPlus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_TopPtMinus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_bTagPlus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_bTagMinus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_cTagPlus->Scale(sfwh160/sig_avgTopwh160);
  h_wh160_cTagMinus->Scale(sfwh160/sig_avgTopwh160);

  
  TH1F* httbar_avgTop = (TH1F*)ttbar->Get("base/SF_topPtWeights"); 
  double ttbar_avgTop = httbar_avgTop->GetMean();
  double sf_top = 1;
  h_ttbar_base->Scale(sf_top/ttbar_avgTop);
  h_ttbar_PileupPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_PileupMinus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_JERPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_JERMinus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_TopPtPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_TopPtMinus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_bTagPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_bTagMinus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_cTagPlus->Scale(sf_top/ttbar_avgTop);
  h_ttbar_cTagMinus->Scale(sf_top/ttbar_avgTop);

  TH1F* h_TotalBkg = (TH1F*)h_wjet_base->Clone("h_TotalBkg");
  h_TotalBkg->Reset();
  h_TotalBkg->Add(h_ttbar_base);
  h_TotalBkg->Add(h_wjet_base);
  h_TotalBkg->Add(h_zjet_base);
  h_TotalBkg->Add(h_qcd_base);
  h_TotalBkg->Add(h_stop_base);
  h_TotalBkg->Add(h_diboson_base);


  TFile *data = new TFile(inFile+"all_muData.root");
  TH1F* h_data = (TH1F*)data->Get("base/"+histPath)->Clone("h_data");


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
  //outFile<<"%\\begin{LARGE}"<<endl;
  outFile<<"\\footnotesize\\setlength{\\tabcolsep}{0.3pt}"<<endl;
  outFile<<"\\begin{tabular}{ | c| c| c| c| c| c| c| c| c| c| c| c| c|c|c|c|}"<<endl; 
  outFile<<"\\multicolumn{5}{c}{ } \\\\"<<endl; 
  outFile<<"\\hline "<<endl;
  if(isMuChannel) outFile<<"\\multicolumn{1}{| p{2cm}| }{Process } & \\multicolumn{1}{ p{0.8cm}|}{\\rotatebox{90}{Luminosity}} & \\multicolumn{1}{ p{0.8cm}|}{\\rotatebox{90}{Pileup reweighting} } & \\multicolumn{1}{ p{0.8cm}|}{\\rotatebox{90}{Lepton selections}} & \\multicolumn{1}{ p{0.8cm}|}{\\rotatebox{90}{JES+JER+MET}} & \\multicolumn{1}{ p{0.8cm}|}{ \\rotatebox{90}{b-jet tagging} } &\\multicolumn{1}{ p{0.8cm}|}{ \\rotatebox{90}{Jet$\\rightarrow bjet ~missID$} } & \\multicolumn{1}{ p{0.8cm}|}{ \\rotatebox{90}{c-jet tagging} }&  \\multicolumn{1}{ p{0.8cm}|}{ \\rotatebox{90}{Jet$\\rightarrow $cjet ~missID} }& \\multicolumn{1}{ p{0.8cm}|}{ \\rotatebox{90}{Normalization}  }& \\multicolumn{1}{ p{0.8cm}|}{\\rotatebox{90}{MC Statistics}  } & \\multicolumn{1}{ p{0.8cm}|}{\\rotatebox{90}{Top-Pt reweighting} }  \\\\ "<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<"\\hline "<<endl;
  int b = 14;

  outFile<<"WH80"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_wh80_base, h_wh80_PileupPlus, h_wh80_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_wh80_JESPlus, h_wh80_base, 
	            h_wh80_JESMinus, h_wh80_JERPlus, h_wh80_JERMinus, b) <<
	  " &  "<< relSysUncBCTag(h_wh80_bTagPlus, h_wh80_base, h_wh80_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_wh80_cTagPlus, h_wh80_base, h_wh80_cTagMinus, b)<<
	  " &  "<<"-"<<
	  " &  "<<"6.1"<<" & "<<relStatUnc(h_wh80_base, b)<<
	  " &  "<<relSysUncTopPt(h_wh80_base, h_wh80_TopPtPlus, h_wh80_TopPtMinus, b)<<
	  " \\\\[0.1cm] \\hline"<<endl;


  outFile<<"WH90"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_wh90_base, h_wh90_PileupPlus, h_wh90_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_wh90_JESPlus, h_wh90_base, 
	            h_wh90_JESMinus, h_wh90_JERPlus, h_wh90_JERMinus, b) <<
	  " &  "<< relSysUncBCTag(h_wh90_bTagPlus, h_wh90_base, h_wh90_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_wh90_cTagPlus, h_wh90_base, h_wh90_cTagMinus, b)<<
	  " &  "<<"-"<<
	  " &  "<<"6.1"<<" & "<<relStatUnc(h_wh90_base, b)<<
	  " &  "<<relSysUncTopPt(h_wh90_base, h_wh90_TopPtPlus, h_wh90_TopPtMinus, b)<<
	  " \\\\[0.1cm] \\hline"<<endl;

  outFile<<"WH100"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_wh100_base, h_wh100_PileupPlus, h_wh100_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_wh100_JESPlus, h_wh100_base, 
	            h_wh100_JESMinus, h_wh100_JERPlus, h_wh100_JERMinus, b) <<
	  " &  "<< relSysUncBCTag(h_wh100_bTagPlus, h_wh100_base, h_wh100_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_wh100_cTagPlus, h_wh100_base, h_wh100_cTagMinus, b)<<
	  " &  "<<"-"<<
	  " &  "<<"6.1"<<" & "<<relStatUnc(h_wh100_base, b)<<
	  " &  "<<relSysUncTopPt(h_wh100_base, h_wh100_TopPtPlus, h_wh100_TopPtMinus, b)<<
	  " \\\\[0.1cm] \\hline"<<endl;

  outFile<<"WH120"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_wh120_base, h_wh120_PileupPlus, h_wh120_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_wh120_JESPlus, h_wh120_base, 
	            h_wh120_JESMinus, h_wh120_JERPlus, h_wh120_JERMinus, b) <<
	  " &  "<< relSysUncBCTag(h_wh120_bTagPlus, h_wh120_base, h_wh120_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_wh120_cTagPlus, h_wh120_base, h_wh120_cTagMinus, b)<<
	  " &  "<<"-"<<
	  " &  "<<"6.1"<<" & "<<relStatUnc(h_wh120_base, b)<<
	  " &  "<<relSysUncTopPt(h_wh120_base, h_wh120_TopPtPlus, h_wh120_TopPtMinus, b)<<
	  " \\\\[0.1cm] \\hline"<<endl;

  outFile<<"WH140"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_wh140_base, h_wh140_PileupPlus, h_wh140_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_wh140_JESPlus, h_wh140_base, 
	            h_wh140_JESMinus, h_wh140_JERPlus, h_wh140_JERMinus, b) <<
	  " &  "<< relSysUncBCTag(h_wh140_bTagPlus, h_wh140_base, h_wh140_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_wh140_cTagPlus, h_wh140_base, h_wh140_cTagMinus, b)<<
	  " &  "<<"-"<<
	  " &  "<<"6.1"<<" & "<<relStatUnc(h_wh140_base, b)<<
	  " &  "<<relSysUncTopPt(h_wh140_base, h_wh140_TopPtPlus, h_wh140_TopPtMinus, b)<<
	  " \\\\[0.1cm] \\hline"<<endl;

  outFile<<"WH150"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_wh150_base, h_wh150_PileupPlus, h_wh150_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_wh150_JESPlus, h_wh150_base, 
	            h_wh150_JESMinus, h_wh150_JERPlus, h_wh150_JERMinus, b) <<
	  " &  "<< relSysUncBCTag(h_wh150_bTagPlus, h_wh150_base, h_wh150_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_wh150_cTagPlus, h_wh150_base, h_wh150_cTagMinus, b)<<
	  " &  "<<"-"<<
	  " &  "<<"6.1"<<" & "<<relStatUnc(h_wh150_base, b)<<
	  " &  "<<relSysUncTopPt(h_wh150_base, h_wh150_TopPtPlus, h_wh150_TopPtMinus, b)<<
	  " \\\\[0.1cm] \\hline"<<endl;

  outFile<<"WH155"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_wh155_base, h_wh155_PileupPlus, h_wh155_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_wh155_JESPlus, h_wh155_base, 
	            h_wh155_JESMinus, h_wh155_JERPlus, h_wh155_JERMinus, b) <<
	  " &  "<< relSysUncBCTag(h_wh155_bTagPlus, h_wh155_base, h_wh155_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_wh155_cTagPlus, h_wh155_base, h_wh155_cTagMinus, b)<<
	  " &  "<<"-"<<
	  " &  "<<"6.1"<<" & "<<relStatUnc(h_wh155_base, b)<<
	  " &  "<<relSysUncTopPt(h_wh155_base, h_wh155_TopPtPlus, h_wh155_TopPtMinus, b)<<
	  " \\\\[0.1cm] \\hline"<<endl;

  outFile<<"WH160"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_wh160_base, h_wh160_PileupPlus, h_wh160_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_wh160_JESPlus, h_wh160_base, 
	            h_wh160_JESMinus, h_wh160_JERPlus, h_wh160_JERMinus, b) <<
	  " &  "<< relSysUncBCTag(h_wh160_bTagPlus, h_wh160_base, h_wh160_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_wh160_cTagPlus, h_wh160_base, h_wh160_cTagMinus, b)<<
	  " &  "<<"-"<<
	  " &  "<<"6.1"<<" & "<<relStatUnc(h_wh160_base, b)<<
	  " &  "<<relSysUncTopPt(h_wh160_base, h_wh160_TopPtPlus, h_wh160_TopPtMinus, b)<<
	  " \\\\[0.1cm] \\hline"<<endl;


  outFile<<"$t\\bar{t}$ + jets"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_ttbar_base, h_ttbar_PileupPlus, h_ttbar_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_ttbar_JESPlus, h_ttbar_base, 
	            h_ttbar_JESMinus, h_ttbar_JERPlus, h_ttbar_JERMinus, b) <<
	  " &  "<< relSysUncBCTag(h_ttbar_bTagPlus, h_ttbar_base, h_ttbar_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_ttbar_cTagPlus, h_ttbar_base, h_ttbar_cTagMinus, b)<<
	  " &  "<<"-"<<
	  " &  "<<"6.1"<<" & "<<relStatUnc(h_ttbar_base, b)<<
	  " &  "<<relSysUncTopPt(h_ttbar_base, h_ttbar_TopPtPlus, h_ttbar_TopPtMinus, b)<<
	  " \\\\[0.1cm] \\hline"<<endl;

  outFile<<"Single ~t"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_stop_base, h_stop_PileupPlus, h_stop_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_stop_JESPlus, h_stop_base, 
	            h_stop_JESMinus, h_stop_JERPlus, h_stop_JERMinus, b) <<
	  " &  "<< relSysUncBCTag(h_stop_bTagPlus, h_stop_base, h_stop_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_stop_cTagPlus, h_stop_base, h_stop_cTagMinus, b)<<
	  " &  "<<"-"<<
	  " &  "<<"5.0"<<" & "<<relStatUnc(h_stop_base, b)<<
	  " &  "<<"-"<<
	  " \\\\[0.1cm] \\hline"<<endl;

  outFile<<"W + jets"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_wjet_base, h_wjet_PileupPlus, h_wjet_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_wjet_JESPlus, h_wjet_base, 
	            h_wjet_JESMinus, h_wjet_JERPlus, h_wjet_JERMinus, b) <<
	  " &  "<<"-"<<
	  " &  "<< relSysUncBCTag(h_wjet_bTagPlus, h_wjet_base, h_wjet_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_wjet_cTagPlus, h_wjet_base, h_wjet_cTagMinus, b)<<
	  " &  "<<"5.0"<<" & "<<relStatUnc(h_wjet_base, b)<<
	  " &  "<<"-"<<
	  " \\\\[0.1cm] \\hline"<<endl;

  outFile<<"$Z/\\gamma$ + jets"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_zjet_base, h_zjet_PileupPlus, h_zjet_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_zjet_JESPlus, h_zjet_base, 
	            h_zjet_JESMinus, h_zjet_JERPlus, h_zjet_JERMinus, b) <<
	  " &  "<<"-"<<
	  " &  "<< relSysUncBCTag(h_zjet_bTagPlus, h_zjet_base, h_zjet_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_zjet_cTagPlus, h_zjet_base, h_zjet_cTagMinus, b)<<
	  " &  "<<"4.5"<<" & "<<relStatUnc(h_zjet_base, b)<<
	  " &  "<<"-"<<
	  " \\\\[0.1cm] \\hline"<<endl;

  outFile<<"VV"<< " & "<<2.5<< std::setprecision(2)<< 
	  " &  "<< relSysUncBCTag(h_diboson_base, h_diboson_PileupPlus, h_diboson_PileupMinus, b)<<
	  " &  "<<3.3<< " & "<<relSysUncJetMET( h_diboson_JESPlus, h_diboson_base, 
	            h_diboson_JESMinus, h_diboson_JERPlus, h_diboson_JERMinus, b) <<
	  " &  "<<"-"<<
	  " &  "<< relSysUncBCTag(h_diboson_bTagPlus, h_diboson_base, h_diboson_bTagMinus, b) << 
	  " &  "<<"-"<<
	  " &  "<<relSysUncBCTag(h_diboson_cTagPlus, h_diboson_base, h_diboson_cTagMinus, b)<<
	  " &  "<<"4.0"<<" & "<<relStatUnc(h_diboson_base, b)<<
	  " &  "<<"-"<<
	  " \\\\[0.1cm] \\hline"<<endl;

  outFile<<"QCD"<< " & -"<<
	  " &  "<<"-"<<
	  " &  "<<"-"<<
	  " &  "<<"-"<<
	  " &  "<<"-"<<
	  " &  "<<"-"<<
	  " &  "<<"-"<<
	  " &  "<<"-"<<
	  " &  "<<"60"<<" & "<<relStatUnc(h_wh80_base, b)<<
	  " &  "<<"-"<<
	  " \\\\[0.1cm] \\hline"<<endl;
/*
  outFile<<"$t\\bar{t}$"<< "&"<<std::setprecision(2)<<h_ttbar_base->GetBinContent(b)<<"$\\pm$ "<< 
  relStatUnc(h_ttbar_base, b)<<" $\\pm$ "<< 
  relSysUncJetMET( h_ttbar_JESPlus, h_ttbar_base, h_ttbar_JESMinus, h_ttbar_JERPlus, 
		  h_ttbar_JERMinus, b) <<" $\\pm$ "<< 
  relSysUncBCTag(h_ttbar_base, h_ttbar_PileupPlus, h_ttbar_PileupMinus, b)<<" $\\pm$ "<<
  relSysUncTopPt(h_ttbar_base, h_ttbar_TopPtPlus, h_ttbar_TopPtMinus, b)<<" $\\pm$ "<<
  relSysUncBCTag(h_ttbar_bTagPlus, h_ttbar_base, h_ttbar_bTagMinus, b) << " $\\pm$ "<< 
  relSysUncBCTag(h_ttbar_cTagPlus, h_ttbar_base, h_ttbar_cTagMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  


  outFile<<"W + Jets"<< " & "<<std::setprecision(2)<<h_wjet_base->GetBinContent(b)<<" $\\pm$ "<< 
  relStatUnc(h_wjet_base, b)<<" $\\pm$ "<< 
  relSysUncJetMET( h_wjet_JESPlus, h_wjet_base, h_wjet_JESMinus, h_wjet_JERPlus, 
		  h_wjet_JERMinus, b) <<" $\\pm$ - "<< 
  relSysUncBCTag(h_wjet_bTagPlus, h_wjet_base, h_wjet_bTagMinus, b) << " $\\pm$ "<< 
  relSysUncBCTag(h_wjet_PileupPlus, h_wjet_base, h_wjet_PileupMinus, b) << " $\\pm$ "<< 
  relSysUncBCTag(h_wjet_cTagPlus, h_wjet_base, h_wjet_cTagMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  

  outFile<<"Z + Jets"<< " & "<<std::setprecision(2)<<h_zjet_base->GetBinContent(b)<<" $\\pm$ "<< 
  relStatUnc(h_zjet_base, b)<<" $\\pm$ "<< 
  relSysUncJetMET( h_zjet_JESPlus, h_zjet_base, h_zjet_JESMinus, h_zjet_JERPlus, 
		  h_zjet_JERMinus, b) <<" $\\pm$ - "<< 
  relSysUncBCTag(h_zjet_PileupPlus, h_zjet_base, h_zjet_PileupMinus, b) << " $\\pm$ "<< 
  relSysUncBCTag(h_zjet_bTagPlus, h_zjet_base, h_zjet_bTagMinus, b) << " $\\pm$ "<< 
  relSysUncBCTag(h_zjet_cTagPlus, h_zjet_base, h_zjet_cTagMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  
  outFile<<"Single t"<< " & "<<std::setprecision(2)<<h_stop_base->GetBinContent(b)<<" $\\pm$ "<< 
  relStatUnc(h_stop_base, b)<<" $\\pm$ "<< 
  relSysUncJetMET( h_stop_JESPlus, h_stop_base, h_stop_JESMinus, h_stop_JERPlus, 
		  h_stop_JERMinus, b) <<" $\\pm$ - "<< 
  relSysUncBCTag(h_stop_PileupPlus, h_stop_base, h_stop_PileupMinus, b) << " $\\pm$ "<< 
  relSysUncBCTag(h_stop_bTagPlus, h_stop_base, h_stop_bTagMinus, b) << " $\\pm$ "<< 
  relSysUncBCTag(h_stop_cTagPlus, h_stop_base, h_stop_cTagMinus, b)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  

  outFile<<"VV"<< " & "<<std::setprecision(2)<<h_diboson_base->GetBinContent(b)<<" $\\pm$ "<< 
  relStatUnc(h_diboson_base, b)<<" $\\pm$ "<< 
  relSysUncJetMET( h_diboson_JESPlus, h_diboson_base, h_diboson_JESMinus, h_diboson_JERPlus, 
		  h_diboson_JERMinus, b) <<" $\\pm$ - "<< 
  relSysUncBCTag(h_diboson_PileupPlus, h_diboson_base, h_diboson_PileupMinus, b) << " $\\pm$ "<< 
  relSysUncBCTag(h_diboson_bTagPlus, h_diboson_base, h_diboson_bTagMinus, b) << " $\\pm$ "<< 
  relSysUncBCTag(h_diboson_cTagPlus, h_diboson_base, h_diboson_cTagMinus, b)<<" \\\\ "<<endl;
  */
  outFile<<"\\hline "<<endl;  
  outFile<<"\\end{tabular}"<<endl; 
  //outFile<<"\\end{LARGE}"<<endl;  
  outFile<<"\\end{center}"<<endl;
  outFile<<"\\caption{First table}"<<endl;
  outFile<<"\\end{table}"<<endl;
  outFile<<"\\end{document}"<<endl;  
  outFile.close(); 
} 
