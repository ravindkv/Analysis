
//TFile *fData     = TFile::Open("all_EleData.root");
TFile *fData     = TFile::Open("all_muData.root");
//TFile *fSig     = TFile::Open("all_Hplus150.root");
TFile* fSig	= TFile::Open("all_TTJetsP.root");
TFile* fTT	= TFile::Open("all_TTJetsP.root");

//--------------------------------------------
//function to make 2d-histo from two array
//--------------------------------------------
TH2F* get2dHisto(TFile *histFile, TString histName, TString xaxis_title){
  TH2F* hist = (TH2F*)(histFile->Get("base//"+histName))->Clone(histName);
  //gStyle->SetPalette(55);
  hist->GetXaxis()->SetTitle(xaxis_title);
  return hist;
}
TH2F* decorate2dHisto(TH2F* hist, TString myTit, TString xTit, TString yTit, TString zTit, int color){
  hist->SetTitle(myTit);
  hist->SetFillColor(color);
  hist->GetXaxis()->SetTitle(xTit);
  hist->GetYaxis()->SetTitle(yTit);
  hist->GetZaxis()->SetTitle(zTit);
  //hist->GetYaxis()->SetRangeUser(1, 3000);
  hist->GetYaxis()->SetRangeUser(1, 1.4* hist->GetMaximum());
  hist->GetYaxis()->SetTitleOffset(1.00);
  hist->GetXaxis()->SetTitleOffset(1.00);
  hist->GetYaxis()->SetTitleSize(0.07);   hist->GetXaxis()->SetTitleSize(0.07);
  hist->GetXaxis()->SetLabelSize(0.07);   hist->GetXaxis()->LabelsOption("u"); // extra
  hist->GetYaxis()->SetLabelSize(0.07);   hist->GetXaxis()->LabelsOption("u"); // extra
  hist->GetXaxis()->SetTickLength(0.05); 
  hist->GetXaxis()->SetNdivisions(5); 
  hist->GetYaxis()->SetTickLength(0.04); 
  //hist->SetMinimum(pMin);
  return hist;
}

void getAll2dHisto(){
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas();
  c1->cd();
  gPad->SetBottomMargin(0.20);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  TH2F *hist = (TH2F*) get2dHisto(fData, "RelIso_MET" , "RelIso vs MET");
  //decorate2dHisto(hist, "#mu + jets", "I_{rel}^{#mu}", "E_{T}^{miss}", "", 1);
  decorate2dHisto(hist, "e + jets", "I_{rel}^{e}", "E_{T}^{miss}", "", 1);
  hist->Draw("colz");
  //c1->SaveAs("RelIso_MET_mu.pdf");
  c1->SaveAs("RelIso_MET_ele.pdf");
}
void getAndDeco2dHisto(TFile *histFile, TString dir, TString histName, TString xaxis_title){
  TH2F* hist = (TH2F*)(histFile->Get("base/Iso/"+dir+"/"+histName))->Clone(histName);
  TCanvas *c1 = new TCanvas();
  //gStyle->SetPalette(55);
  gStyle->SetOptStat(111);
  gStyle->SetFrameLineWidth(4);
  gPad->SetBottomMargin(0.20);
  gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  hist->GetXaxis()->SetTitle(xaxis_title);
  hist->SetTitle("Binning: 10%, 10%, 10% : "+dir+": #mu+jets: t#bar{t}");
  //hist->SetFillColor(color);
  //hist->GetXaxis()->SetTitle(xTit);
  //hist->GetYaxis()->SetTitle(yTit);
  //hist->GetZaxis()->SetTitle(zTit);
  //hist->GetYaxis()->SetRangeUser(1, 3000);
  hist->GetYaxis()->SetTitle("M_{csb_{Lep}} [GeV]");
  //hist->GetYaxis()->SetRangeUser(1, 1.4* hist->GetMaximum());
  hist->GetYaxis()->SetTitleOffset(1.00);
  hist->GetXaxis()->SetRangeUser(0, 300);
  hist->GetXaxis()->SetTitleOffset(1.20);
  hist->GetYaxis()->SetTitleSize(0.07);   hist->GetXaxis()->SetTitleSize(0.07);
  hist->GetXaxis()->SetLabelSize(0.07);   hist->GetXaxis()->LabelsOption("u"); // extra
  hist->GetYaxis()->SetLabelSize(0.07);   hist->GetXaxis()->LabelsOption("u"); // extra
  hist->GetXaxis()->SetTickLength(0.05); 
  hist->GetXaxis()->SetNdivisions(5); 
  hist->GetYaxis()->SetNdivisions(5); 
  hist->GetYaxis()->SetTickLength(0.04); 
  hist->Draw("colz");
  c1->SaveAs(dir+"_"+histName+".pdf");
}
void getTProfile(TFile *histFile, TString dir, TString histName, TString xaxis_title){
  TProfile* hist = (TProfile*)(histFile->Get("base/Iso/"+dir+"/"+histName))->Clone(histName);
  TCanvas *c1 = new TCanvas();
  //gStyle->SetPalette(55);
  gStyle->SetOptStat(111);
  gStyle->SetFrameLineWidth(4);
  gPad->SetBottomMargin(0.20);
  gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  hist->GetXaxis()->SetTitle(xaxis_title);
  hist->SetTitle("Binning: 10%, 10%, 10% : "+dir+": #mu+jets: t#bar{t}");
  //hist->GetYaxis()->SetRangeUser(1, 3000);
  hist->GetYaxis()->SetTitle("M_{jj} [GeV]");
  //hist->GetYaxis()->SetRangeUser(1, 1.4* hist->GetMaximum());
  hist->GetYaxis()->SetTitleOffset(1.00);
  //hist->GetXaxis()->SetRangeUser(0, 300);
  hist->GetXaxis()->SetTitleOffset(1.20);
  hist->GetYaxis()->SetTitleSize(0.07);   hist->GetXaxis()->SetTitleSize(0.07);
  hist->GetXaxis()->SetLabelSize(0.07);   hist->GetXaxis()->LabelsOption("u"); // extra
  hist->GetYaxis()->SetLabelSize(0.07);   hist->GetXaxis()->LabelsOption("u"); // extra
  hist->GetXaxis()->SetTickLength(0.05); 
  //hist->GetXaxis()->SetNdivisions(5); 
  hist->SetLineColor(3);
  hist->SetLineWidth(3);
  hist->SetMarkerSize(1);
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(4); 
  hist->GetYaxis()->SetNdivisions(5); 
  hist->GetYaxis()->SetTickLength(0.04); 
  hist->Draw();
  c1->SaveAs(dir+"_"+histName+".pdf");
}

void getHisto(TFile *histFile, TString dir, TString histName, TString xaxis_title){
  TH1F* hist; 
  if(!(histFile->Get("base/Iso/"+dir+"/"+histName))){
    hist = (TH1F*)(fTT->Get("base/Iso/"+dir+"/"+histName))->Clone(histName);
    hist->Add(hist, -1);
  }else hist = (TH1F*)(histFile->Get("base/Iso/"+dir+"/"+histName))->Clone(histName);
  
  TCanvas *c1 = new TCanvas();
  gStyle->SetOptStat(11111111);
  gStyle->SetFrameLineWidth(4);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  hist->GetYaxis()->SetTitle("Events");
  //hist->GetYaxis()->SetRangeUser(1, 3000);
  hist->SetLineWidth(3);
  //hist->SetTitle(dir+": #mu+jets: WH150");
  //hist->SetTitle(dir+": #mu+jets: t#bar{t}");
  //hist->SetTitle("Binning: 10g, 10g, 10g : "+dir+": #mu+jets: t#bar{t}");
  //hist->SetTitle("Binning: 5%, 5%, 10% : "+dir+": #mu+jets: t#bar{t}");
  //hist->SetTitle("Binning: 5%, 5%, 20% : "+dir+": #mu+jets: t#bar{t}");
  hist->SetTitle("Binning: 10%, 10%, 10% : "+dir+": #mu+jets: t#bar{t}");
  hist->GetYaxis()->SetRangeUser(1, 1.2* hist->GetMaximum());
  hist->GetXaxis()->SetRangeUser(0, 200);
  hist->GetXaxis()->SetTitle(xaxis_title);
  hist->GetXaxis()->SetTitleOffset(1.15);
  hist->GetYaxis()->SetTitleSize(0.05);   hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);   hist->GetXaxis()->LabelsOption("u"); // extra
  hist->GetYaxis()->SetLabelSize(0.05);   hist->GetXaxis()->LabelsOption("u"); // extra
  hist->GetXaxis()->SetTickLength(0.05); 
  hist->GetYaxis()->SetTickLength(0.04); 
  hist->GetYaxis()->SetTitleSize(0.07);   
  hist->GetYaxis()->SetTitleOffset(1.05);
  cout<<dir<<":"<<histName<<": "<<setw(10)<<hist->Integral()<<setw(10)<<hist->GetMean()<<setw(10)<<hist->GetRMS()<<endl;
  hist->Draw();
  c1->SaveAs(dir+"_"+histName+".pdf");
  //c1->Close();
}


void getAllHistoBin10p10p10p(){
  /*
  getHisto(fSig, "KinFit", "mjj_kfit" , "M_{jj}^{Inc}[GeV]");
  getHisto(fSig, "PtbJetInc", "mjj_kfit_25To35" , "M_{jj}^{Inc}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetInc", "mjj_kfit_35To42" , "M_{jj}^{Inc}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetInc", "mjj_kfit_42To50" , "M_{jj}^{Inc}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetInc", "mjj_kfit_50To57" , "M_{jj}^{Inc}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetInc", "mjj_kfit_57To65" , "M_{jj}^{Inc}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetInc", "mjj_kfit_65To74" , "M_{jj}^{Inc}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetInc", "mjj_kfit_74To84" , "M_{jj}^{Inc}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetInc", "mjj_kfit_84To99" , "M_{jj}^{Inc}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetInc", "mjj_kfit_99To124", "M_{jj}^{Inc}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetInc", "mjj_kfit_124To500","M_{jj}^{Inc}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");
  */
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_25To35" , "M_{jj}^{CatL}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_35To42" , "M_{jj}^{CatL}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_42To50" , "M_{jj}^{CatL}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_50To57" , "M_{jj}^{CatL}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_57To65" , "M_{jj}^{CatL}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_65To74" , "M_{jj}^{CatL}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_74To84" , "M_{jj}^{CatL}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_84To99" , "M_{jj}^{CatL}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_99To124", "M_{jj}^{CatL}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_124To500","M_{jj}^{CatL}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatM", "mjj_kfit_25To35" , "M_{jj}^{CatM}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_35To42" , "M_{jj}^{CatM}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_42To50" , "M_{jj}^{CatM}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_50To57" , "M_{jj}^{CatM}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_57To65" , "M_{jj}^{CatM}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_65To74" , "M_{jj}^{CatM}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_74To84" , "M_{jj}^{CatM}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_84To99" , "M_{jj}^{CatM}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_99To124", "M_{jj}^{CatM}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_124To500","M_{jj}^{CatM}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatT", "mjj_kfit_25To35" , "M_{jj}^{CatT}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_35To42" , "M_{jj}^{CatT}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_42To50" , "M_{jj}^{CatT}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_50To57" , "M_{jj}^{CatT}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_57To65" , "M_{jj}^{CatT}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_65To74" , "M_{jj}^{CatT}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_74To84" , "M_{jj}^{CatT}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_84To99" , "M_{jj}^{CatT}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_99To124", "M_{jj}^{CatT}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_124To500","M_{jj}^{CatT}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");
 }

void getAllHistoBin5p5p20p(){
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_25To32",  "M_{jj}^{CatL}(25 < Pt_{bjet}^{Had} #leq 32)[GeV]"); 
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_32To37",  "M_{jj}^{CatL}(32 < Pt_{bjet}^{Had} #leq 37)[GeV]"); 
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_37To41",  "M_{jj}^{CatL}(37 < Pt_{bjet}^{Had} #leq 41)[GeV]");  
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_41To44",  "M_{jj}^{CatL}(41 < Pt_{bjet}^{Had} #leq 44)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_44To48",  "M_{jj}^{CatL}(44 < Pt_{bjet}^{Had} #leq 48)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_48To52",  "M_{jj}^{CatL}(48 < Pt_{bjet}^{Had} #leq 52)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_52To55",  "M_{jj}^{CatL}(52 < Pt_{bjet}^{Had} #leq 55)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_55To59",  "M_{jj}^{CatL}(55 < Pt_{bjet}^{Had} #leq 59)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_59To63",  "M_{jj}^{CatL}(59 < Pt_{bjet}^{Had} #leq 63)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_63To67",  "M_{jj}^{CatL}(63 < Pt_{bjet}^{Had} #leq 67)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_67To72",  "M_{jj}^{CatL}(67 < Pt_{bjet}^{Had} #leq 72)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_72To76",  "M_{jj}^{CatL}(72 < Pt_{bjet}^{Had} #leq 76)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_76To81",  "M_{jj}^{CatL}(76 < Pt_{bjet}^{Had} #leq 81)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_81To87",  "M_{jj}^{CatL}(81 < Pt_{bjet}^{Had} #leq 87)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_87To93",  "M_{jj}^{CatL}(87 < Pt_{bjet}^{Had} #leq 93)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_93To102", "M_{jj}^{CatL}(93 < Pt_{bjet}^{Had} #leq 102)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_102To112","M_{jj}^{CatL}(102 < Pt_{bjet}^{Had} #leq 112)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_112To127","M_{jj}^{CatL}(112 < Pt_{bjet}^{Had} #leq 127)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_127To152","M_{jj}^{CatL}(127 < Pt_{bjet}^{Had} #leq 152)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_152To500","M_{jj}^{CatL}(152 < Pt_{bjet}^{Had} #leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatM", "mjj_kfit_25To32",  "M_{jj}^{CatM}(25 < Pt_{bjet}^{Had} #leq 32)[GeV]"); 
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_32To37",  "M_{jj}^{CatM}(32 < Pt_{bjet}^{Had} #leq 37)[GeV]"); 
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_37To41",  "M_{jj}^{CatM}(37 < Pt_{bjet}^{Had} #leq 41)[GeV]");  
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_41To44",  "M_{jj}^{CatM}(41 < Pt_{bjet}^{Had} #leq 44)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_44To48",  "M_{jj}^{CatM}(44 < Pt_{bjet}^{Had} #leq 48)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_48To52",  "M_{jj}^{CatM}(48 < Pt_{bjet}^{Had} #leq 52)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_52To55",  "M_{jj}^{CatM}(52 < Pt_{bjet}^{Had} #leq 55)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_55To59",  "M_{jj}^{CatM}(55 < Pt_{bjet}^{Had} #leq 59)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_59To63",  "M_{jj}^{CatM}(59 < Pt_{bjet}^{Had} #leq 63)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_63To67",  "M_{jj}^{CatM}(63 < Pt_{bjet}^{Had} #leq 67)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_67To72",  "M_{jj}^{CatM}(67 < Pt_{bjet}^{Had} #leq 72)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_72To76",  "M_{jj}^{CatM}(72 < Pt_{bjet}^{Had} #leq 76)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_76To81",  "M_{jj}^{CatM}(76 < Pt_{bjet}^{Had} #leq 81)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_81To87",  "M_{jj}^{CatM}(81 < Pt_{bjet}^{Had} #leq 87)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_87To93",  "M_{jj}^{CatM}(87 < Pt_{bjet}^{Had} #leq 93)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_93To102", "M_{jj}^{CatM}(93 < Pt_{bjet}^{Had} #leq 102)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_102To112","M_{jj}^{CatM}(102 < Pt_{bjet}^{Had} #leq 112)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_112To127","M_{jj}^{CatM}(112 < Pt_{bjet}^{Had} #leq 127)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_127To152","M_{jj}^{CatM}(127 < Pt_{bjet}^{Had} #leq 152)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_152To500","M_{jj}^{CatM}(152 < Pt_{bjet}^{Had} #leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatT", "mjj_kfit_25To43"  , "M_{jj}^{CatT}(25 < Pt_{bjet}^{Had} #leq 43)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_43To57"  , "M_{jj}^{CatT}(43 < Pt_{bjet}^{Had} #leq 57)[GeV$]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_57To74"  , "M_{jj}^{CatT}(57 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_74To99"  , "M_{jj}^{CatT}(74 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_99To500" , "M_{jj}^{CatT}(99 < Pt_{bjet}^{Had} #leq 500)[GeV]");
 }
void getAllHistoBin5p5p10p(){
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_25To32",  "M_{jj}^{CatL}(25 < Pt_{bjet}^{Had} #leq 32)[GeV]"); 
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_32To37",  "M_{jj}^{CatL}(32 < Pt_{bjet}^{Had} #leq 37)[GeV]"); 
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_37To41",  "M_{jj}^{CatL}(37 < Pt_{bjet}^{Had} #leq 41)[GeV]");  
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_41To44",  "M_{jj}^{CatL}(41 < Pt_{bjet}^{Had} #leq 44)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_44To48",  "M_{jj}^{CatL}(44 < Pt_{bjet}^{Had} #leq 48)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_48To52",  "M_{jj}^{CatL}(48 < Pt_{bjet}^{Had} #leq 52)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_52To55",  "M_{jj}^{CatL}(52 < Pt_{bjet}^{Had} #leq 55)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_55To59",  "M_{jj}^{CatL}(55 < Pt_{bjet}^{Had} #leq 59)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_59To63",  "M_{jj}^{CatL}(59 < Pt_{bjet}^{Had} #leq 63)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_63To67",  "M_{jj}^{CatL}(63 < Pt_{bjet}^{Had} #leq 67)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_67To72",  "M_{jj}^{CatL}(67 < Pt_{bjet}^{Had} #leq 72)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_72To76",  "M_{jj}^{CatL}(72 < Pt_{bjet}^{Had} #leq 76)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_76To81",  "M_{jj}^{CatL}(76 < Pt_{bjet}^{Had} #leq 81)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_81To87",  "M_{jj}^{CatL}(81 < Pt_{bjet}^{Had} #leq 87)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_87To93",  "M_{jj}^{CatL}(87 < Pt_{bjet}^{Had} #leq 93)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_93To102", "M_{jj}^{CatL}(93 < Pt_{bjet}^{Had} #leq 102)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_102To112","M_{jj}^{CatL}(102 < Pt_{bjet}^{Had} #leq 112)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_112To127","M_{jj}^{CatL}(112 < Pt_{bjet}^{Had} #leq 127)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_127To152","M_{jj}^{CatL}(127 < Pt_{bjet}^{Had} #leq 152)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_152To500","M_{jj}^{CatL}(152 < Pt_{bjet}^{Had} #leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatM", "mjj_kfit_25To32",  "M_{jj}^{CatM}(25 < Pt_{bjet}^{Had} #leq 32)[GeV]"); 
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_32To37",  "M_{jj}^{CatM}(32 < Pt_{bjet}^{Had} #leq 37)[GeV]"); 
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_37To41",  "M_{jj}^{CatM}(37 < Pt_{bjet}^{Had} #leq 41)[GeV]");  
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_41To44",  "M_{jj}^{CatM}(41 < Pt_{bjet}^{Had} #leq 44)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_44To48",  "M_{jj}^{CatM}(44 < Pt_{bjet}^{Had} #leq 48)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_48To52",  "M_{jj}^{CatM}(48 < Pt_{bjet}^{Had} #leq 52)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_52To55",  "M_{jj}^{CatM}(52 < Pt_{bjet}^{Had} #leq 55)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_55To59",  "M_{jj}^{CatM}(55 < Pt_{bjet}^{Had} #leq 59)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_59To63",  "M_{jj}^{CatM}(59 < Pt_{bjet}^{Had} #leq 63)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_63To67",  "M_{jj}^{CatM}(63 < Pt_{bjet}^{Had} #leq 67)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_67To72",  "M_{jj}^{CatM}(67 < Pt_{bjet}^{Had} #leq 72)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_72To76",  "M_{jj}^{CatM}(72 < Pt_{bjet}^{Had} #leq 76)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_76To81",  "M_{jj}^{CatM}(76 < Pt_{bjet}^{Had} #leq 81)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_81To87",  "M_{jj}^{CatM}(81 < Pt_{bjet}^{Had} #leq 87)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_87To93",  "M_{jj}^{CatM}(87 < Pt_{bjet}^{Had} #leq 93)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_93To102", "M_{jj}^{CatM}(93 < Pt_{bjet}^{Had} #leq 102)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_102To112","M_{jj}^{CatM}(102 < Pt_{bjet}^{Had} #leq 112)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_112To127","M_{jj}^{CatM}(112 < Pt_{bjet}^{Had} #leq 127)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_127To152","M_{jj}^{CatM}(127 < Pt_{bjet}^{Had} #leq 152)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_152To500","M_{jj}^{CatM}(152 < Pt_{bjet}^{Had} #leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatT", "mjj_kfit_25To35" , "M_{jj}^{CatT}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_35To42" , "M_{jj}^{CatT}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_42To50" , "M_{jj}^{CatT}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_50To57" , "M_{jj}^{CatT}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_57To65" , "M_{jj}^{CatT}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_65To74" , "M_{jj}^{CatT}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_74To84" , "M_{jj}^{CatT}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_84To99" , "M_{jj}^{CatT}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_99To124", "M_{jj}^{CatT}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_124To500","M_{jj}^{CatT}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");
 }

void getAllHistoBin10g10g10g(){

  getHisto(fSig, "PtbJetCatL", "mjj_kfit_25To35" , "M_{jj}^{CatL}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_35To45" , "M_{jj}^{CatL}(35 < Pt_{bjet}^{Had} #leq 45)[GeV$]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_45To55" , "M_{jj}^{CatL}(45 < Pt_{bjet}^{Had} #leq 55)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_55To65" , "M_{jj}^{CatL}(55 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_65To75" , "M_{jj}^{CatL}(65 < Pt_{bjet}^{Had} #leq 75)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_75To85" , "M_{jj}^{CatL}(75 < Pt_{bjet}^{Had} #leq 85)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_85To95" , "M_{jj}^{CatL}(85 < Pt_{bjet}^{Had} #leq 95)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_95To105" , "M_{jj}^{CatL}(95 < Pt_{bjet}^{Had} #leq 105)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_105To150", "M_{jj}^{CatL}(105 < Pt_{bjet}^{Had} #leq 150)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_150To500","M_{jj}^{CatL}(150 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatM", "mjj_kfit_25To35" , "M_{jj}^{CatM}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_35To45" , "M_{jj}^{CatM}(35 < Pt_{bjet}^{Had} #leq 45)[GeV$]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_45To55" , "M_{jj}^{CatM}(45 < Pt_{bjet}^{Had} #leq 55)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_55To65" , "M_{jj}^{CatM}(55 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_65To75" , "M_{jj}^{CatM}(65 < Pt_{bjet}^{Had} #leq 75)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_75To85" , "M_{jj}^{CatM}(75 < Pt_{bjet}^{Had} #leq 85)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_85To95" , "M_{jj}^{CatM}(85 < Pt_{bjet}^{Had} #leq 95)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_95To105" , "M_{jj}^{CatM}(95 < Pt_{bjet}^{Had} #leq 105)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_105To150", "M_{jj}^{CatM}(105 < Pt_{bjet}^{Had} #leq 150)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_150To500","M_{jj}^{CatM}(150 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatT", "mjj_kfit_25To35" , "M_{jj}^{CatT}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_35To45" , "M_{jj}^{CatT}(35 < Pt_{bjet}^{Had} #leq 45)[GeV$]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_45To55" , "M_{jj}^{CatT}(45 < Pt_{bjet}^{Had} #leq 55)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_55To65" , "M_{jj}^{CatT}(55 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_65To75" , "M_{jj}^{CatT}(65 < Pt_{bjet}^{Had} #leq 75)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_75To85" , "M_{jj}^{CatT}(75 < Pt_{bjet}^{Had} #leq 85)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_85To95" , "M_{jj}^{CatT}(85 < Pt_{bjet}^{Had} #leq 95)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_95To105" , "M_{jj}^{CatT}(95 < Pt_{bjet}^{Had} #leq 105)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_105To150", "M_{jj}^{CatT}(105 < Pt_{bjet}^{Had} #leq 150)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_150To500","M_{jj}^{CatT}(150 < Pt_{bjet}^{Had}#leq 500)[GeV]");
 }
void getAllHistoBin10p10p10p_triJetcsbHad(){
  getHisto(fSig, "PtbJetCatL", "triJet_csbHad_25To35" , "M_{csb_{Had}}^{CatL}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbHad_35To42" , "M_{csb_{Had}}^{CatL}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbHad_42To50" , "M_{csb_{Had}}^{CatL}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbHad_50To57" , "M_{csb_{Had}}^{CatL}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbHad_57To65" , "M_{csb_{Had}}^{CatL}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbHad_65To74" , "M_{csb_{Had}}^{CatL}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbHad_74To84" , "M_{csb_{Had}}^{CatL}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbHad_84To99" , "M_{csb_{Had}}^{CatL}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbHad_99To124", "M_{csb_{Had}}^{CatL}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbHad_124To500","M_{csb_{Had}}^{CatL}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatM", "triJet_csbHad_25To35" , "M_{csb_{Had}}^{CatM}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbHad_35To42" , "M_{csb_{Had}}^{CatM}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbHad_42To50" , "M_{csb_{Had}}^{CatM}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbHad_50To57" , "M_{csb_{Had}}^{CatM}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbHad_57To65" , "M_{csb_{Had}}^{CatM}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbHad_65To74" , "M_{csb_{Had}}^{CatM}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbHad_74To84" , "M_{csb_{Had}}^{CatM}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbHad_84To99" , "M_{csb_{Had}}^{CatM}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbHad_99To124", "M_{csb_{Had}}^{CatM}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbHad_124To500","M_{csb_{Had}}^{CatM}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatT", "triJet_csbHad_25To35" , "M_{csb_{Had}}^{CatT}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbHad_35To42" , "M_{csb_{Had}}^{CatT}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbHad_42To50" , "M_{csb_{Had}}^{CatT}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbHad_50To57" , "M_{csb_{Had}}^{CatT}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbHad_57To65" , "M_{csb_{Had}}^{CatT}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbHad_65To74" , "M_{csb_{Had}}^{CatT}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbHad_74To84" , "M_{csb_{Had}}^{CatT}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbHad_84To99" , "M_{csb_{Had}}^{CatT}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbHad_99To124", "M_{csb_{Had}}^{CatT}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbHad_124To500","M_{csb_{Had}}^{CatT}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");
 }
void getAllHistoBin10p10p10p_triJetcsbLep(){
  getHisto(fSig, "PtbJetCatL", "triJet_csbLep_25To35" , "M_{csb_{Lep}}^{CatL}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbLep_35To42" , "M_{csb_{Lep}}^{CatL}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbLep_42To50" , "M_{csb_{Lep}}^{CatL}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbLep_50To57" , "M_{csb_{Lep}}^{CatL}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbLep_57To65" , "M_{csb_{Lep}}^{CatL}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbLep_65To74" , "M_{csb_{Lep}}^{CatL}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbLep_74To84" , "M_{csb_{Lep}}^{CatL}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbLep_84To99" , "M_{csb_{Lep}}^{CatL}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbLep_99To124", "M_{csb_{Lep}}^{CatL}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatL", "triJet_csbLep_124To500","M_{csb_{Lep}}^{CatL}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatM", "triJet_csbLep_25To35" , "M_{csb_{Lep}}^{CatM}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbLep_35To42" , "M_{csb_{Lep}}^{CatM}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbLep_42To50" , "M_{csb_{Lep}}^{CatM}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbLep_50To57" , "M_{csb_{Lep}}^{CatM}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbLep_57To65" , "M_{csb_{Lep}}^{CatM}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbLep_65To74" , "M_{csb_{Lep}}^{CatM}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbLep_74To84" , "M_{csb_{Lep}}^{CatM}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbLep_84To99" , "M_{csb_{Lep}}^{CatM}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbLep_99To124", "M_{csb_{Lep}}^{CatM}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatM", "triJet_csbLep_124To500","M_{csb_{Lep}}^{CatM}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatT", "triJet_csbLep_25To35" , "M_{csb_{Lep}}^{CatT}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbLep_35To42" , "M_{csb_{Lep}}^{CatT}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbLep_42To50" , "M_{csb_{Lep}}^{CatT}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbLep_50To57" , "M_{csb_{Lep}}^{CatT}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbLep_57To65" , "M_{csb_{Lep}}^{CatT}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbLep_65To74" , "M_{csb_{Lep}}^{CatT}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbLep_74To84" , "M_{csb_{Lep}}^{CatT}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbLep_84To99" , "M_{csb_{Lep}}^{CatT}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbLep_99To124", "M_{csb_{Lep}}^{CatT}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatT", "triJet_csbLep_124To500","M_{csb_{Lep}}^{CatT}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");
 }
//TH2F* getAndDeco2dHisto(TFile *histFile, TString histName, TString xaxis_title){
void getAllHistoBin10p10p10p_triJet_csbHad_csbLep_50To57(){
  getAndDeco2dHisto(fSig, "PtbJetCatL", "triJet_csbHad_csbLep_25To35" , "M_{csb_{Had}}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatL", "triJet_csbHad_csbLep_35To42" , "M_{csb_{Had}}(35 < Pt_{bjet}^{Had} #leq 42)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatL", "triJet_csbHad_42To50_csbLep" , "M_{csb_{Had}}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatL", "triJet_csbHad_csbLep_50To57" , "M_{csb_{Had}}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatL", "triJet_csbHad_csbLep_57To65" , "M_{csb_{Had}}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatL", "triJet_csbHad_csbLep_65To74" , "M_{csb_{Had}}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatL", "triJet_csbHad_csbLep_74To84" , "M_{csb_{Had}}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatL", "triJet_csbHad_csbLep_84To99" , "M_{csb_{Had}}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatL", "triJet_csbHad_csbLep_99To124", "M_{csb_{Had}}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatL", "triJet_csbHad_csbLep_124To500","M_{csb_{Had}}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getAndDeco2dHisto(fSig, "PtbJetCatM", "triJet_csbHad_csbLep_25To35" , "M_{csb_{Had}}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatM", "triJet_csbHad_csbLep_35To42" , "M_{csb_{Had}}(35 < Pt_{bjet}^{Had} #leq 42)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatM", "triJet_csbHad_42To50_csbLep" , "M_{csb_{Had}}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatM", "triJet_csbHad_csbLep_50To57" , "M_{csb_{Had}}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatM", "triJet_csbHad_csbLep_57To65" , "M_{csb_{Had}}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatM", "triJet_csbHad_csbLep_65To74" , "M_{csb_{Had}}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatM", "triJet_csbHad_csbLep_74To84" , "M_{csb_{Had}}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatM", "triJet_csbHad_csbLep_84To99" , "M_{csb_{Had}}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatM", "triJet_csbHad_csbLep_99To124", "M_{csb_{Had}}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatM", "triJet_csbHad_csbLep_124To500","M_{csb_{Had}}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getAndDeco2dHisto(fSig, "PtbJetCatT", "triJet_csbHad_csbLep_25To35" , "M_{csb_{Had}}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatT", "triJet_csbHad_csbLep_35To42" , "M_{csb_{Had}}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getAndDeco2dHisto(fSig, "PtbJetCatT", "triJet_csbHad_42To50_csbLep" , "M_{csb_{Had}}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatT", "triJet_csbHad_csbLep_50To57" , "M_{csb_{Had}}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatT", "triJet_csbHad_csbLep_57To65" , "M_{csb_{Had}}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatT", "triJet_csbHad_csbLep_65To74" , "M_{csb_{Had}}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatT", "triJet_csbHad_csbLep_74To84" , "M_{csb_{Had}}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatT", "triJet_csbHad_csbLep_84To99" , "M_{csb_{Had}}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatT", "triJet_csbHad_csbLep_99To124", "M_{csb_{Had}}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getAndDeco2dHisto(fSig, "PtbJetCatT", "triJet_csbHad_csbLep_124To500","M_{csb_{Had}}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");
 }
void getAllHistoBin10p10p10p_ClassA(){
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_25To35_ClassA" , "M_{jj}^{ClassA}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_35To42_ClassA" , "M_{jj}^{ClassA}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_42To50_ClassA" , "M_{jj}^{ClassA}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_50To57_ClassA" , "M_{jj}^{ClassA}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_57To65_ClassA" , "M_{jj}^{ClassA}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_65To74_ClassA" , "M_{jj}^{ClassA}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_74To84_ClassA" , "M_{jj}^{ClassA}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_84To99_ClassA" , "M_{jj}^{ClassA}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_99To124_ClassA", "M_{jj}^{ClassA}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_124To500_ClassA","M_{jj}^{ClassA}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatM", "mjj_kfit_25To35_ClassA" , "M_{jj}^{ClassA}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_35To42_ClassA" , "M_{jj}^{ClassA}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_42To50_ClassA" , "M_{jj}^{ClassA}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_50To57_ClassA" , "M_{jj}^{ClassA}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_57To65_ClassA" , "M_{jj}^{ClassA}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_65To74_ClassA" , "M_{jj}^{ClassA}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_74To84_ClassA" , "M_{jj}^{ClassA}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_84To99_ClassA" , "M_{jj}^{ClassA}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_99To124_ClassA", "M_{jj}^{ClassA}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_124To500_ClassA","M_{jj}^{ClassA}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatT", "mjj_kfit_25To35_ClassA" , "M_{jj}^{ClassA}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_35To42_ClassA" , "M_{jj}^{ClassA}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_42To50_ClassA" , "M_{jj}^{ClassA}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_50To57_ClassA" , "M_{jj}^{ClassA}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_57To65_ClassA" , "M_{jj}^{ClassA}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_65To74_ClassA" , "M_{jj}^{ClassA}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_74To84_ClassA" , "M_{jj}^{ClassA}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_84To99_ClassA" , "M_{jj}^{ClassA}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_99To124_ClassA", "M_{jj}^{ClassA}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_124To500_ClassA","M_{jj}^{ClassA}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");
 }
void getAllHistoBin10p10p10p_ClassB(){
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_25To35_ClassB" , "M_{jj}^{ClassB}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_35To42_ClassB" , "M_{jj}^{ClassB}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_42To50_ClassB" , "M_{jj}^{ClassB}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_50To57_ClassB" , "M_{jj}^{ClassB}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_57To65_ClassB" , "M_{jj}^{ClassB}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_65To74_ClassB" , "M_{jj}^{ClassB}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_74To84_ClassB" , "M_{jj}^{ClassB}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_84To99_ClassB" , "M_{jj}^{ClassB}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_99To124_ClassB", "M_{jj}^{ClassB}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_124To500_ClassB","M_{jj}^{ClassB}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatM", "mjj_kfit_25To35_ClassB" , "M_{jj}^{ClassB}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_35To42_ClassB" , "M_{jj}^{ClassB}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_42To50_ClassB" , "M_{jj}^{ClassB}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_50To57_ClassB" , "M_{jj}^{ClassB}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_57To65_ClassB" , "M_{jj}^{ClassB}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_65To74_ClassB" , "M_{jj}^{ClassB}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_74To84_ClassB" , "M_{jj}^{ClassB}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_84To99_ClassB" , "M_{jj}^{ClassB}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_99To124_ClassB", "M_{jj}^{ClassB}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_124To500_ClassB","M_{jj}^{ClassB}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatT", "mjj_kfit_25To35_ClassB" , "M_{jj}^{ClassB}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_35To42_ClassB" , "M_{jj}^{ClassB}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_42To50_ClassB" , "M_{jj}^{ClassB}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_50To57_ClassB" , "M_{jj}^{ClassB}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_57To65_ClassB" , "M_{jj}^{ClassB}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_65To74_ClassB" , "M_{jj}^{ClassB}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_74To84_ClassB" , "M_{jj}^{ClassB}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_84To99_ClassB" , "M_{jj}^{ClassB}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_99To124_ClassB", "M_{jj}^{ClassB}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_124To500_ClassB","M_{jj}^{ClassB}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");
 }

void getAllHistoBin10p10p10p_ClassAB(){
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_25To35" , "M_{jj}^{ClassAB}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_35To42" , "M_{jj}^{ClassAB}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_42To50" , "M_{jj}^{ClassAB}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_50To57" , "M_{jj}^{ClassAB}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_57To65" , "M_{jj}^{ClassAB}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_65To74" , "M_{jj}^{ClassAB}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_74To84" , "M_{jj}^{ClassAB}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_84To99" , "M_{jj}^{ClassAB}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_99To124", "M_{jj}^{ClassAB}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatL", "mjj_kfit_124To500","M_{jj}^{ClassAB}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatM", "mjj_kfit_25To35" , "M_{jj}^{ClassAB}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_35To42" , "M_{jj}^{ClassAB}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_42To50" , "M_{jj}^{ClassAB}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_50To57" , "M_{jj}^{ClassAB}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_57To65" , "M_{jj}^{ClassAB}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_65To74" , "M_{jj}^{ClassAB}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_74To84" , "M_{jj}^{ClassAB}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_84To99" , "M_{jj}^{ClassAB}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_99To124", "M_{jj}^{ClassAB}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatM", "mjj_kfit_124To500","M_{jj}^{ClassAB}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");

  getHisto(fSig, "PtbJetCatT", "mjj_kfit_25To35" , "M_{jj}^{ClassAB}(25 < Pt_{bjet}^{Had} #leq 35)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_35To42" , "M_{jj}^{ClassAB}(35 < Pt_{bjet}^{Had} #leq 42)[GeV$]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_42To50" , "M_{jj}^{ClassAB}(42 < Pt_{bjet}^{Had} #leq 50)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_50To57" , "M_{jj}^{ClassAB}(50 < Pt_{bjet}^{Had} #leq 57)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_57To65" , "M_{jj}^{ClassAB}(57 < Pt_{bjet}^{Had} #leq 65)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_65To74" , "M_{jj}^{ClassAB}(65 < Pt_{bjet}^{Had} #leq 74)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_74To84" , "M_{jj}^{ClassAB}(74 < Pt_{bjet}^{Had} #leq 84)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_84To99" , "M_{jj}^{ClassAB}(84 < Pt_{bjet}^{Had} #leq 99)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_99To124", "M_{jj}^{ClassAB}(99 < Pt_{bjet}^{Had} #leq 124)[GeV]");
  getHisto(fSig, "PtbJetCatT", "mjj_kfit_124To500","M_{jj}^{ClassAB}(124 < Pt_{bjet}^{Had}#leq 500)[GeV]");
 }
void getAllTProfileBin10p10p10p_mjj_kfit_eta_bjetH(){
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_bjetH_25To35" , "#eta_{bjet}^{Had}(25 < Pt_{bjet}^{Had} #leq 35 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_bjetH_35To42" , "#eta_{bjet}^{Had}(35 < Pt_{bjet}^{Had} #leq 42 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_bjetH_42To50" , "#eta_{bjet}^{Had}(42 < Pt_{bjet}^{Had} #leq 50 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_bjetH_50To57" , "#eta_{bjet}^{Had}(50 < Pt_{bjet}^{Had} #leq 57 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_bjetH_57To65" , "#eta_{bjet}^{Had}(57 < Pt_{bjet}^{Had} #leq 65 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_bjetH_65To74" , "#eta_{bjet}^{Had}(65 < Pt_{bjet}^{Had} #leq 74 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_bjetH_74To84" , "#eta_{bjet}^{Had}(74 < Pt_{bjet}^{Had} #leq 84 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_bjetH_84To99" , "#eta_{bjet}^{Had}(84 < Pt_{bjet}^{Had} #leq 99 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_bjetH_99To124", "#eta_{bjet}^{Had}(99 < Pt_{bjet}^{Had} #leq 124 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_bjetH_124To500","#eta_{bjet}^{Had}(124 < Pt_{bjet}^{Had}#leq 500 GeV)");

  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_bjetH_25To35" , "#eta_{bjet}^{Had}(25 < Pt_{bjet}^{Had} #leq 35 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_bjetH_35To42" , "#eta_{bjet}^{Had}(35 < Pt_{bjet}^{Had} #leq 42 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_bjetH_42To50" , "#eta_{bjet}^{Had}(42 < Pt_{bjet}^{Had} #leq 50 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_bjetH_50To57" , "#eta_{bjet}^{Had}(50 < Pt_{bjet}^{Had} #leq 57 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_bjetH_57To65" , "#eta_{bjet}^{Had}(57 < Pt_{bjet}^{Had} #leq 65 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_bjetH_65To74" , "#eta_{bjet}^{Had}(65 < Pt_{bjet}^{Had} #leq 74 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_bjetH_74To84" , "#eta_{bjet}^{Had}(74 < Pt_{bjet}^{Had} #leq 84 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_bjetH_84To99" , "#eta_{bjet}^{Had}(84 < Pt_{bjet}^{Had} #leq 99 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_bjetH_99To124", "#eta_{bjet}^{Had}(99 < Pt_{bjet}^{Had} #leq 124 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_bjetH_124To500","#eta_{bjet}^{Had}(124 < Pt_{bjet}^{Had}#leq 500 GeV)");

  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_bjetH_25To35" , "#eta_{bjet}^{Had}(25 < Pt_{bjet}^{Had} #leq 35 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_bjetH_35To42" , "#eta_{bjet}^{Had}(35 < Pt_{bjet}^{Had} #leq 42 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_bjetH_42To50" , "#eta_{bjet}^{Had}(42 < Pt_{bjet}^{Had} #leq 50 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_bjetH_50To57" , "#eta_{bjet}^{Had}(50 < Pt_{bjet}^{Had} #leq 57 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_bjetH_57To65" , "#eta_{bjet}^{Had}(57 < Pt_{bjet}^{Had} #leq 65 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_bjetH_65To74" , "#eta_{bjet}^{Had}(65 < Pt_{bjet}^{Had} #leq 74 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_bjetH_74To84" , "#eta_{bjet}^{Had}(74 < Pt_{bjet}^{Had} #leq 84 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_bjetH_84To99" , "#eta_{bjet}^{Had}(84 < Pt_{bjet}^{Had} #leq 99 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_bjetH_99To124", "#eta_{bjet}^{Had}(99 < Pt_{bjet}^{Had} #leq 124 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_bjetH_124To500","#eta_{bjet}^{Had}(124 < Pt_{bjet}^{Had}#leq 500 GeV)");
 }
void getAllTProfileBin10p10p10p_mjj_kfit_eta_diJet(){
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_diJet_25To35" , "#eta_{jj}(25 < Pt_{bjet}^{Had} #leq 35 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_diJet_35To42" , "#eta_{jj}(35 < Pt_{bjet}^{Had} #leq 42 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_diJet_42To50" , "#eta_{jj}(42 < Pt_{bjet}^{Had} #leq 50 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_diJet_50To57" , "#eta_{jj}(50 < Pt_{bjet}^{Had} #leq 57 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_diJet_57To65" , "#eta_{jj}(57 < Pt_{bjet}^{Had} #leq 65 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_diJet_65To74" , "#eta_{jj}(65 < Pt_{bjet}^{Had} #leq 74 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_diJet_74To84" , "#eta_{jj}(74 < Pt_{bjet}^{Had} #leq 84 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_diJet_84To99" , "#eta_{jj}(84 < Pt_{bjet}^{Had} #leq 99 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_diJet_99To124", "#eta_{jj}(99 < Pt_{bjet}^{Had} #leq 124 GeV)");
  getTProfile(fSig, "PtbJetCatL", "mjj_kfit_eta_diJet_124To500","#eta_{jj}(124 < Pt_{bjet}^{Had}#leq 500 GeV)");

  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_diJet_25To35" , "#eta_{jj}(25 < Pt_{bjet}^{Had} #leq 35 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_diJet_35To42" , "#eta_{jj}(35 < Pt_{bjet}^{Had} #leq 42 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_diJet_42To50" , "#eta_{jj}(42 < Pt_{bjet}^{Had} #leq 50 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_diJet_50To57" , "#eta_{jj}(50 < Pt_{bjet}^{Had} #leq 57 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_diJet_57To65" , "#eta_{jj}(57 < Pt_{bjet}^{Had} #leq 65 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_diJet_65To74" , "#eta_{jj}(65 < Pt_{bjet}^{Had} #leq 74 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_diJet_74To84" , "#eta_{jj}(74 < Pt_{bjet}^{Had} #leq 84 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_diJet_84To99" , "#eta_{jj}(84 < Pt_{bjet}^{Had} #leq 99 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_diJet_99To124", "#eta_{jj}(99 < Pt_{bjet}^{Had} #leq 124 GeV)");
  getTProfile(fSig, "PtbJetCatM", "mjj_kfit_eta_diJet_124To500","#eta_{jj}(124 < Pt_{bjet}^{Had}#leq 500 GeV)");

  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_diJet_25To35" , "#eta_{jj}(25 < Pt_{bjet}^{Had} #leq 35 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_diJet_35To42" , "#eta_{jj}(35 < Pt_{bjet}^{Had} #leq 42 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_diJet_42To50" , "#eta_{jj}(42 < Pt_{bjet}^{Had} #leq 50 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_diJet_50To57" , "#eta_{jj}(50 < Pt_{bjet}^{Had} #leq 57 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_diJet_57To65" , "#eta_{jj}(57 < Pt_{bjet}^{Had} #leq 65 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_diJet_65To74" , "#eta_{jj}(65 < Pt_{bjet}^{Had} #leq 74 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_diJet_74To84" , "#eta_{jj}(74 < Pt_{bjet}^{Had} #leq 84 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_diJet_84To99" , "#eta_{jj}(84 < Pt_{bjet}^{Had} #leq 99 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_diJet_99To124", "#eta_{jj}(99 < Pt_{bjet}^{Had} #leq 124 GeV)");
  getTProfile(fSig, "PtbJetCatT", "mjj_kfit_eta_diJet_124To500","#eta_{jj}(124 < Pt_{bjet}^{Had}#leq 500 GeV)");
 }
