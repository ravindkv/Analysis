#include "MyStackHisto.h"

///////////////////////////////////////////  
//CHANNEL
bool isMuChannel = true;
bool isEleChannel = false;

bool isDDqcd = true;
TFile *f_QCD_dd = new TFile("all_QCD_dd.root","RECREATE");
  
//SAVE HISTOS ON DISK
bool isSaveHisto = true;
///////////////////////////////////////////  

void plotStackedHisto(TString baseDir, TString isoDir, TString histDir, TString histName, TString xTitle,bool isData=false, bool isSig=false, double xmin=0, double xmax=10, double unc=false){
  MyStackHisto MyS;
  string hist_name (histName);

  //user's input for data driven qcd
  double qcd_sf_btag     =  0.0;
  double qcd_sf_kfit     =  0.0;
  double qcd_sf_btag_err =  0.0;
  double qcd_sf_kfit_err =  0.0;
  //muon channel
  if(isMuChannel){
    qcd_sf_btag     =  2.652 ;
    qcd_sf_kfit     =  1.504 ;
    qcd_sf_btag_err =  0.1504;
    qcd_sf_kfit_err =  0.2206;
  }
  //electron channel
  if(isEleChannel){
    qcd_sf_btag     =  2.428 ;
    qcd_sf_kfit     =  1.996 ;
    qcd_sf_btag_err =  0.1331;
    qcd_sf_kfit_err =  0.2419;
  }
  //qcd scale factors for data-driven QCD
  //vector<double> sfAndErr = MyS.getQcdSF(fData, fTT, fST, fWJ, fDY, fVV, histDir, histName);
  //double qcdSF = sfAndErr[0];
  //double qcdErr = sfAndErr[1];
  //Pad
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(3);
  TCanvas *c1 = new TCanvas();
  //TCanvas *c1 = new TCanvas("c1", "Data_MC", 400, 600);
  const float xpad[2] = {0.,1};
  const float ypad[4] = {0.,0.30,0.30,0.98};
  if(isData){
    c1->Divide(1, 2, 0, 0); c1->cd(1);
    gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);
    if(isData) gPad->SetLogy(true);
    if(hist_name.find("mjj") != string::npos) gPad->SetLogy(false);
  }

  //-------------------------------
  //Legends
  //-------------------------------
  TLegend* leg = new TLegend(0.7518792,0.3061504,0.9512081,0.8798861,NULL,"brNDC");
  //TLegend* leg = new TLegend(0.7618792,0.3061504,0.9712081,0.8798861,NULL,"brNDC");
  if(hist_name.find("pt") != string::npos || hist_name.find("mt") != string::npos || hist_name.find("Fit") != string::npos ||hist_name.find("RelIso") != string::npos){
    leg = new TLegend(0.6018792,0.6061504,0.9912081,0.8898861,NULL,"brNDC");
    leg->SetNColumns(2);
  }
  leg->SetFillStyle(0); leg->SetBorderSize(0);
  leg->SetFillColor(10); leg->SetTextSize(0.06);
  
  //-------------------------------
  // stack MC Bkg histo
  //-------------------------------
  THStack* hStack = new THStack("hStack","");
  TH1F* hVV = MyS.getHisto(fVV, baseDir, isoDir, histDir, histName);
  TH1F* hMC = (TH1F*)hVV->Clone("hMC");
  int col_depth =2;
  if(baseDir=="baseLowMET/") col_depth = 1;
  hVV->SetFillColor(11 + col_depth);
  leg->AddEntry(hVV,"VV","F");
  hStack->Add(hVV);
  MyS.stackHisto(fDY, "Z/#gamma+jets", baseDir, isoDir, histDir, histName, kViolet+col_depth, 1,   hStack, hMC, leg);
  MyS.stackHisto(fWJ, "W+ jets", baseDir, isoDir, histDir, histName, kPink +col_depth , 1,   hStack, hMC, leg);
  
  // trim the histDir string
  std::string histDir_str;
  std::string histDir_orig(histDir);
  std::remove_copy(histDir_orig.begin(), histDir_orig.end(), std::back_inserter(histDir_str), '/');
  TString histDir_(histDir_str);
  // QCD from Data
  ///if(histDir=="") isDDqcd = false;
  TH1F * hQCD_dd = MyS.getHisto(fQCD, baseDir, isoDir, histDir, histName);
  hQCD_dd->Reset(); // initialize empty hist
  if(isDDqcd){
    //hQCD_dd = MyS.getDDqcd(baseDir, "NonIso/", histDir, histName,  qcdSF,  qcdErr);
    if(baseDir=="baseLowMET/")hQCD_dd = MyS.getDDqcd(baseDir, "NonIso/", histDir, histName); // don't apply QCD sf in lowMET region.
    else if(histDir=="BTag/") hQCD_dd = MyS.getDDqcd(baseDir, "NonIso/", histDir, histName, qcd_sf_btag, qcd_sf_btag_err);
    else if(histDir=="KinFit/") hQCD_dd = MyS.getDDqcd(baseDir, "NonIso/", histDir, histName,  qcd_sf_kfit,  qcd_sf_kfit_err);
    else hQCD_dd = MyS.getDDqcd(baseDir, "NonIso/", histDir, histName,  qcd_sf_kfit,  qcd_sf_kfit_err);
    hQCD_dd->SetFillColor(kGreen+2);
    hQCD_dd->GetXaxis()->SetRangeUser(xmin,xmax);
    //create same dir to the data driven qcd file
    std::string histPath = std::string(baseDir+isoDir+histDir_);
    TDirectory *d = f_QCD_dd->GetDirectory(histPath.c_str());
    if(!d) f_QCD_dd->mkdir(histPath.c_str());
    f_QCD_dd->cd(histPath.c_str());
    //hQCD->Draw();
    hQCD_dd->Write();
    leg->AddEntry(hQCD_dd,"QCD","F");
    hStack->Add(hQCD_dd);
    hMC->Add(hQCD_dd);
  }
  else MyS.stackHisto(fQCD, "QCD", baseDir, isoDir, histDir, histName, kGreen +col_depth, 1,   hStack, hMC, leg);
  MyS.stackHisto(fST, "Single t", baseDir, isoDir, histDir, histName, kYellow+col_depth , 1,   hStack, hMC, leg);
  MyS.stackHisto(fTT,"t#bar{t} + jets", baseDir, isoDir, histDir, histName, kCyan+col_depth, 1,   hStack, hMC, leg);

  gPad->SetTopMargin(0.10);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  hStack->Draw("HIST");
  hStack->SetMinimum(1.0);
  hStack ->GetXaxis()->SetRangeUser(xmin, xmax);
  //cout<<hStack->GetMaximum()<<endl;
  if(isData){
    hStack->GetYaxis()->SetTitleOffset(0.70);
    hStack->GetYaxis()->SetTitleSize(0.10);   
    hStack->GetYaxis()->SetLabelSize(0.09);   
    hStack->GetYaxis()->SetTickLength(0.04); 
    hStack->GetYaxis()->SetTitle("Events");
    hStack->GetXaxis()->SetTitleOffset(1.20);
  }
  else{
  hStack->GetYaxis()->SetTitle("Events");
  hStack->GetXaxis()->SetTitle(xTitle);
  hStack->GetXaxis()->SetTitleSize(0.07);
  hStack->GetXaxis()->SetLabelSize(0.07);   
  hStack->GetXaxis()->SetTickLength(0.05); 
  if(histDir_str.find("PtbJet") != string::npos)
    hStack->GetXaxis()->SetNdivisions(5);
  hStack->GetYaxis()->SetNdivisions(5);
  hStack->GetYaxis()->SetTitleSize(0.10);   
  hStack->GetYaxis()->SetTitleOffset(1.15);
  hStack->GetYaxis()->SetLabelSize(0.07);   
  hStack->GetYaxis()->SetTickLength(0.04); 
  gPad->SetLeftMargin(0.22);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.20);
  hStack->GetXaxis()->SetTitleOffset(1.20);
  }

  //-------------------------------------///
  //unc band
  //-------------------------------------///
  if(unc){
  TGraphAsymmErrors *UncBand;
  UncBand = MyS.UNCGRAPH(
            MyS.addHistoForUnc("base/", 	 isoDir, histDir, histName,     isDDqcd),
      	    MyS.addHistoForUnc("JESPlus/",      isoDir, histDir, histName,     isDDqcd),
      	    MyS.addHistoForUnc("JESMinus/",     isoDir, histDir, histName,     isDDqcd),
      	    MyS.addHistoForUnc("JERPlus/",      isoDir, histDir, histName,     isDDqcd),
      	    MyS.addHistoForUnc("JERMinus/",     isoDir, histDir, histName,     isDDqcd),
      	    MyS.addHistoForUnc("bcTagPlus1/",     isoDir, histDir, histName,     isDDqcd),
      	    MyS.addHistoForUnc("bcTagMinus1/",    isoDir, histDir, histName,     isDDqcd),
      	    MyS.addHistoForUnc("bcTagPlus2/",     isoDir, histDir, histName,     isDDqcd),
      	    MyS.addHistoForUnc("bcTagMinus2/",    isoDir, histDir, histName,     isDDqcd),
      	    MyS.addHistoForUnc("bcTagPlus3/",     isoDir, histDir, histName,     isDDqcd),
      	    MyS.addHistoForUnc("bcTagMinus3/",    isoDir, histDir, histName,     isDDqcd),
      	    MyS.addHistoForUnc("PileupPlus/",   isoDir, histDir, histName,     isDDqcd),
      	    MyS.addHistoForUnc("PileupMinus/",  isoDir, histDir, histName, 	isDDqcd),
	    hQCD_dd, true, false);
  UncBand->SetFillColor(kSpring +9);
  UncBand->SetFillStyle(3008);
  UncBand->Draw(" E2 same");
  leg->AddEntry(UncBand, "Unc","F"); 
  }

  //-------------------------------
  //Data
  //-------------------------------
  TH1F* hData = MyS.getHisto(fData, baseDir, isoDir, histDir, histName);;
  ///MyS.decorateHisto(hData, "", xTitle, "Events");
  hData->SetFillColor(kBlack);
  hData->SetMarkerStyle(20); hData->SetMarkerSize(1.2);
  if(isData)hData->Draw("Esame"); 
  if(isData)leg->AddEntry(hData,"Data","PE"); 

  //-------------------------------
  //Signal 
  //-------------------------------
  TH1F* hSig = MyS.getHisto(fSig, baseDir, isoDir, histDir, histName);
  if(histDir_str.find("PtbJet") != string::npos)
    hSig = MyS.getHisto(fSig, baseDir, isoDir, histDir, histName);
  ///MyS.decorateHisto(hSig, "", xTitle, "Events");
  hSig->SetLineColor(kRed); hSig->SetLineStyle(2);
  hSig->SetLineWidth(3); hSig->SetFillColor(0);
  if(isSig)hSig->Draw("HISTSAME"); 
  leg->AddEntry(hSig, "Signal","L"); 
  //leg->AddEntry(hSig, "#splitline{Signal}{M_{H^{+}} = 120 GeV}","L"); 
  leg->Draw();
  double yMax = 0;
  if(hData->GetMaximum() > hSig->GetMaximum()) yMax = hData->GetMaximum();
  else yMax = hSig->GetMaximum();
  if(yMax < hMC->GetMaximum()) yMax = hMC->GetMaximum();

  if(isData) hStack->SetMaximum(4.0*hStack->GetMaximum());
  else hStack->SetMaximum(1.1*yMax);
  if(baseDir=="baseLowMET/" && hist_name.find("mjj") != string::npos)
	  hStack->SetMaximum(1.1*yMax);

  //TPaveText *cct = MyS.paveText(0.70,0.8554,0.80,0.8562, 0, 19, 1, 0, 132);
  TPaveText *cct = MyS.paveText(0.40,0.8454,0.50,0.8462, 0, 19, 1, 0, 132);
  cct->SetTextSize(0.09);
  if(histDir_str.find("PtbJet") != string::npos)
    cct->AddText("M_{H^{+}} = 90 GeV");
  else cct->AddText("M_{H^{+}} = 120 GeV");

  //-------------------------------------///
  //  Draw Pave Text 
  //-------------------------------------///
  //hist name
  TPaveText *hLable = MyS.paveText(0.6513423,0.7754898,0.6010067,0.8962187, 0, 19, 1, 0, 132);
  hLable->SetTextSize(0.080);
  hLable->AddText(xTitle);
  
  //channel
  TPaveText *ch = MyS.paveText(0.823,0.9154898,0.9210067,0.9762187, 0, 19, 1, 0, 132);
  ch->SetTextSize(0.12);
  if(isMuChannel) ch->AddText("#mu + jets");
  if(isEleChannel) ch->AddText("e + jets");
  //CMS prili
  TPaveText *pt = MyS.paveText(0.01,0.9554,0.82,0.9562, 0, 19, 1, 0, 132);
  if(isData) pt->SetTextSize(0.085);
  else pt->SetTextSize(0.059);
  if(baseDir=="baseLowMET/") pt->AddText(histDir_+"(E_{T}^{miss} < 20 GeV): #sqrt{s} = 13 TeV, 35.9 fb^{-1}");
  else pt->AddText(histDir_+": #sqrt{s} = 13 TeV, 35.9 fb^{-1}");
  //TText *text = pt->AddText(histDir+": CMS Preliminary, #sqrt{s} = 13 TeV, 35.9 fb^{-1}");
  pt->Draw();
  if(isSig) cct->Draw();
  ch->Draw();
  //hLable->Draw();
  gPad->RedrawAxis();
  c1->Update();
  
  //-------------------------------------///
  // Ratio = DATA/Bkg
  //-------------------------------------///
  if(isData){
    c1->cd(2);
    gPad->SetTopMargin(0); gPad->SetBottomMargin(0.5); //gPad->SetGridy();
    if(histDir=="") gPad->SetBottomMargin(0.55);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05);
    gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
    TH1F *hRatio = (TH1F*)hData->Clone("hRatio");
    hRatio->Reset();
    hRatio->Add(hData);
    hRatio->Divide(hMC); 
    MyS.decorateHisto(hRatio, "", xTitle, "#frac{Data}{Bkg}");
    hRatio->SetFillColor(kBlack);
    if(baseDir=="baseLowMET/") hRatio->GetYaxis()->SetRangeUser(0, 2);
    else hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatio->GetXaxis()->SetRangeUser(xmin, xmax);
    hRatio->GetYaxis()->SetTitleOffset(0.27);
    hRatio->GetXaxis()->SetTitleOffset(1.10);
    hRatio->SetMarkerStyle(20); hRatio->SetMarkerSize(1.2);
    hRatio->GetYaxis()->SetTitleSize(0.22); 
    hRatio->GetXaxis()->SetTitleSize(0.23);
    hRatio->GetXaxis()->SetLabelSize(0.22); 
    hRatio->GetYaxis()->SetLabelSize(0.12); 
    if(hist_name.find("mjj") != string::npos){
      hRatio->GetXaxis()->SetTitleSize(0.15); 
      hRatio->GetXaxis()->SetTitleOffset(1.40);
    }
    //lable x-axis, for cutflow
    if(histName=="cutflow"){
      vector<string >cut_label;
      if(isEleChannel){
        cut_label.push_back("Ele trigger");
        cut_label.push_back("No. of ele = 1");
      }
      if(isMuChannel){
        cut_label.push_back("Mu trigger");
        cut_label.push_back("No. of mu = 1");
      }
      cut_label.push_back("No. of jets #geq 4");
      cut_label.push_back("MET #geq 20 GeV");
      cut_label.push_back("No. of b-jets #geq 2");
      cut_label.push_back("KinFit Sel");
      cut_label.push_back("No. of c-jet =1");
      for(int istep=0; istep<cut_label.size(); istep++ ){
       hRatio->GetXaxis()->SetBinLabel(istep+1, cut_label[istep].c_str());
      }
      hRatio->GetXaxis()->LabelsOption("u");
      hRatio->GetXaxis()->SetTickLength(0.08); 
      hRatio->GetXaxis()->SetLabelOffset(0.08);
      hRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
      /*
      gPad->SetBottomMargin(0.6); //gPad->SetGridy();
      hRatio->GetXaxis()->SetLabelOffset(0.04);
      hRatio->GetXaxis()->SetLabelSize(0.12);
      hRatio->GetXaxis()->SetTitleOffset(1.20);
      hRatio->GetYaxis()->SetTitleOffset(1.00);
      */
    }
    //unc band
    hRatio->Draw("E"); // use "P" or "AP"
    if(unc){
    TGraphAsymmErrors *UncBand_Ratio;
    UncBand_Ratio = MyS.UNCGRAPH(
	    MyS.addHistoForUnc("base/", 	 isoDir, histDir, histName,  isDDqcd),
      	    MyS.addHistoForUnc("JESPlus/",     isoDir, histDir, histName,  isDDqcd),
      	    MyS.addHistoForUnc("JESMinus/",    isoDir, histDir, histName,  isDDqcd),
      	    MyS.addHistoForUnc("JERPlus/",     isoDir, histDir, histName,  isDDqcd),
      	    MyS.addHistoForUnc("JERMinus/",    isoDir, histDir, histName,  isDDqcd),
      	    MyS.addHistoForUnc("bcTagPlus1/",  isoDir, histDir, histName,  isDDqcd),
      	    MyS.addHistoForUnc("bcTagMinus1/", isoDir, histDir, histName,  isDDqcd),
      	    MyS.addHistoForUnc("bcTagPlus2/",  isoDir, histDir, histName,  isDDqcd),
      	    MyS.addHistoForUnc("bcTagMinus2/", isoDir, histDir, histName,  isDDqcd),
      	    MyS.addHistoForUnc("bcTagPlus3/",  isoDir, histDir, histName,  isDDqcd),
      	    MyS.addHistoForUnc("bcTagMinus3/", isoDir, histDir, histName,  isDDqcd),
      	    MyS.addHistoForUnc("PileupPlus/",  isoDir, histDir, histName,  isDDqcd),
      	    MyS.addHistoForUnc("PileupMinus/", isoDir, histDir, histName,  isDDqcd),
	    hQCD_dd, false, true);
    UncBand_Ratio->SetFillColor(kSpring +9);
    UncBand_Ratio->SetFillStyle(3008);
    UncBand_Ratio->Draw("E2 same");
    }
    hRatio->Draw("E same"); // use "P" or "AP"
    //base line at 1
    TF1 *baseLine = new TF1("baseLine","1", -100, 2000); 
    baseLine->SetLineColor(kBlack);
    baseLine->Draw("SAME");
    c1->Update();
  }
  if(isSaveHisto){
    mkdir(histDir_, S_IRWXU);
    TString outFile("$PWD/");
    outFile += histDir_+"/"+histName;
    if(isMuChannel) outFile += "_mu"+histDir_+".pdf";
    if(isEleChannel) outFile += "_ele"+histDir_+".pdf";
    c1->SaveAs(outFile);
    //c1->Close();
  }
}

void MyStackHisto(){

  TString histDir="KinFit/"; // BTag/, KinFit/, PtbJetBin/;
  TString baseDir = "base/"; // base/, baseLowMET/;

  TString isoDir = "Iso/";
  bool isDataMjj=true;
  //flags
  bool isData = true;
  bool isSig = true;
  bool isUnc = true;
  bool isMCqcd = false;
  if(isMuChannel && isMCqcd){
   plotStackedHisto(baseDir, "", "", "RelIso_1Mu","I_{rel}^{e}", isData,  isSig,  0.0, 1.0,  false);
   plotStackedHisto(baseDir, "", "", "pt_met_1Mu","E_{T}^{miss}[GeV]", isData,  isSig,  0.0, 500.0, false);
  }
  if(isEleChannel && isMCqcd){
   plotStackedHisto(baseDir, "", "", "RelIso_1Ele","I_{rel}^{e}", isData,  isSig,  0.0, 1.0,  false);
   plotStackedHisto(baseDir, "", "", "pt_met_1Ele","E_{T}^{miss}[GeV]", isData,  isSig,  0.0, 500.0, false);
  }
  if(isMCqcd) plotStackedHisto(baseDir, isoDir, "", "cutflow","cutflow", isData,  isSig,  0.0, 10.0, false);

 if(histDir=="BTag/"){
   plotStackedHisto(baseDir, isoDir, histDir, "mjj", "M_{jj}[GeV]", isDataMjj, isSig,  0, 400, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "pfCISV", "pfCISV2BJetTags", isData, isSig,  -0.1, 1.5, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "chi2OfKinFit", "#chi^{2} of KF", isData, isSig,  0, 100, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "probOfKinFit", "prob of KF", isData, isSig,  0, 1, isUnc);
 }
 if(histDir=="KinFit/"){
   plotStackedHisto(baseDir, isoDir, histDir, "mjj_kfit", "M_{jj}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "mjj_kfit_CTagIncL", "M_{jj}^{Inc_CTagL}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "mjj_kfit_CTagIncM", "M_{jj}^{Inc_CTagM}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "mjj_kfit_CTagIncT", "M_{jj}^{Inc_CTagT}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "mjj_kfit_CTagExO", "M_{jj}^{Ex_CTagO}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "mjj_kfit_CTagExL", "M_{jj}^{Ex_CTagL}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "mjj_kfit_CTagExM", "M_{jj}^{Ex_CTagM}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "mjj_kfit_CTagExT", "M_{jj}^{Ex_CTagT}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "pt_bjetH", "Pt_{bjet}^{Had} [GeV]", isData, isSig,  0, 500, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "pfCCvsL", "pfCombinedCvsLJetTags", isData, isSig,  -1.1, 2.0, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "pfCCvsB", "pfCombinedCvsBJetTags", isData, isSig,  -1.1, 2.0, isUnc);
 }
 if(histDir=="PtbJetBin/"){
   plotStackedHisto(baseDir, isoDir, "KinFit/", "mjj_kfit", "M_{jj}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, "KinFit/", "mjj_kfit_CTagIncL", "M_{jj}^{Inc_CTagL}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, "KinFit/", "mjj_kfit_CTagIncM", "M_{jj}^{Inc_CTagM}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, "KinFit/", "mjj_kfit_CTagIncT", "M_{jj}^{Inc_CTagT}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, "KinFit/", "mjj_kfit_CTagExO", "M_{jj}^{Ex_CTagO}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, "KinFit/", "mjj_kfit_CTagExL", "M_{jj}^{Ex_CTagL}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, "KinFit/", "mjj_kfit_CTagExM", "M_{jj}^{Ex_CTagM}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, "KinFit/", "mjj_kfit_CTagExT", "M_{jj}^{Ex_CTagT}[GeV]", isDataMjj, isSig,  0, 250, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetInc/", "mjj_kfit_25To42", "M_{jj}^{Inc}(25 < Pt_{bjet}^{Had} #leq 42)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetInc/", "mjj_kfit_42To57", "M_{jj}^{Inc}(42 < Pt_{bjet}^{Had} #leq 57)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetInc/", "mjj_kfit_57To74", "M_{jj}^{Inc}(57 < Pt_{bjet}^{Had} #leq 74)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetInc/", "mjj_kfit_74To99", "M_{jj}^{Inc}(74 < Pt_{bjet}^{Had} #leq 99)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetInc/", "mjj_kfit_99To500","M_{jj}^{Inc}(99 < Pt_{bjet}^{Had} #leq 500)[GeV]", isDataMjj, isSig,  0, 200, isUnc);

   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExO/", "mjj_kfit_25To42", "M_{jj}^{Inc}(25 < Pt_{bjet}^{Had} #leq 42)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExO/", "mjj_kfit_42To57", "M_{jj}^{Inc}(42 < Pt_{bjet}^{Had} #leq 57)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExO/", "mjj_kfit_57To74", "M_{jj}^{Inc}(57 < Pt_{bjet}^{Had} #leq 74)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExO/", "mjj_kfit_74To99", "M_{jj}^{Inc}(74 < Pt_{bjet}^{Had} #leq 99)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExO/", "mjj_kfit_99To500","M_{jj}^{Inc}(99 < Pt_{bjet}^{Had} #leq 500)[GeV]", isDataMjj, isSig,  0, 200, isUnc);

   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExL/", "mjj_kfit_25To42", "M_{jj}^{Inc}(25 < Pt_{bjet}^{Had} #leq 42)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExL/", "mjj_kfit_42To57", "M_{jj}^{Inc}(42 < Pt_{bjet}^{Had} #leq 57)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExL/", "mjj_kfit_57To74", "M_{jj}^{Inc}(57 < Pt_{bjet}^{Had} #leq 74)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExL/", "mjj_kfit_74To99", "M_{jj}^{Inc}(74 < Pt_{bjet}^{Had} #leq 99)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExL/", "mjj_kfit_99To500","M_{jj}^{Inc}(99 < Pt_{bjet}^{Had} #leq 500)[GeV]", isDataMjj, isSig,  0, 200, isUnc);

   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExM/", "mjj_kfit_25To42", "M_{jj}^{Inc}(25 < Pt_{bjet}^{Had} #leq 42)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExM/", "mjj_kfit_42To57", "M_{jj}^{Inc}(42 < Pt_{bjet}^{Had} #leq 57)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExM/", "mjj_kfit_57To74", "M_{jj}^{Inc}(57 < Pt_{bjet}^{Had} #leq 74)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExM/", "mjj_kfit_74To99", "M_{jj}^{Inc}(74 < Pt_{bjet}^{Had} #leq 99)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExM/", "mjj_kfit_99To500","M_{jj}^{Inc}(99 < Pt_{bjet}^{Had} #leq 500)[GeV]", isDataMjj, isSig,  0, 200, isUnc);

   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExT/", "mjj_kfit_25To42", "M_{jj}^{Inc}(25 < Pt_{bjet}^{Had} #leq 42)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExT/", "mjj_kfit_42To57", "M_{jj}^{Inc}(42 < Pt_{bjet}^{Had} #leq 57)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExT/", "mjj_kfit_57To74", "M_{jj}^{Inc}(57 < Pt_{bjet}^{Had} #leq 74)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExT/", "mjj_kfit_74To99", "M_{jj}^{Inc}(74 < Pt_{bjet}^{Had} #leq 99)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
   plotStackedHisto(baseDir, isoDir, "PtbJetCTagExT/", "mjj_kfit_99To500","M_{jj}^{Inc}(99 < Pt_{bjet}^{Had} #leq 500)[GeV]", isDataMjj, isSig,  0, 200, isUnc);
  }
 if(histDir=="BTag/" ||histDir=="KinFit/"){
   plotStackedHisto(baseDir, isoDir, histDir, "eta_jet", "#eta^{jets}", isData, isSig,  -3.0, 5.0, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "pt_jet", "Pt^{jets} [GeV]", isData, isSig,  0.0, 700.0, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "nvtx", "N^{vertex}", isData, isSig,  0, 70, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "rhoAll", "#rho", isData, isSig,  0, 70, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "final_multi_jet", "N^{jets} [GeV]", isData, isSig,  3, 15, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "CSVL_count", "N^{b-jets}", isData, isSig,  1, 10, isUnc);
   if(baseDir=="baseLowMET/")plotStackedHisto(baseDir, isoDir, histDir, "final_pt_met", "E_{T}^{miss}[GeV]", isData, isSig,  0.0, 40.0, isUnc);
   else plotStackedHisto(baseDir, isoDir, histDir, "final_pt_met", "E_{T}^{miss}[GeV]", isData, isSig,  0.0, 500.0, isUnc);
   plotStackedHisto(baseDir, isoDir, histDir, "wmt", "MT[GeV]", isData, isSig,  0.0, 250.0, isUnc);
   if(isMuChannel){
    plotStackedHisto(baseDir, isoDir, histDir, "pt_mu", "Pt^{#mu} [GeV]", isData, isSig,  0.0, 500.0, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "eta_mu", "#eta^{#mu} [GeV]", isData, isSig,  -3.0, 5.0, isUnc);
   }
   if(isEleChannel){
     plotStackedHisto(baseDir, isoDir, histDir, "pt_ele", "Pt^{e} [GeV]", isData, isSig,  0.0, 500.0, isUnc);
     plotStackedHisto(baseDir, isoDir, histDir, "eta_ele", "#eta^{e} [GeV]", isData, isSig,  -3.0, 5.0, isUnc);
   }
  }
} 
