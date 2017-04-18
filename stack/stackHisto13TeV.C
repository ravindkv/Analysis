#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

using namespace std;
//TGraphAsymmErrors *FULLGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus,  TH1F *hMETUCPlus, TH1F *hMETUCMinus, TH1F *bTagPlus, TH1F *bTagMinus, TH1F *h_QCD_iso) {
TGraphAsymmErrors *FULLGRAPH(TH1F *hCentral, TH1F *h_QCD_iso) {
  
  TGraphAsymmErrors *gr;
  //int n1 = hJESPlus->GetNbinsX(); int n2 = hJESMinus->GetNbinsX();
  int n1 = 10; int n2 = 10;
  if(n1 != n2) cout<<"ATTENTION!!!!!  BinMismatch"<<endl;
  double *Yval, *errorU, *errorD, *XerrorU, *XerrorD, *Xval ;
  
  Yval = new double[n1]; errorU = new double[n1]; errorD = new double[n1];
  XerrorU=new double[n1]; XerrorD=new double[n1]; Xval=new double[n1];
  cout << "No. of bins= " << n1 << endl;
  for(int i=0; i<n1; i++)
    {
      Yval[i]   = hCentral->GetBinContent(i+1);
      
      //errorU[i] = sqrt(pow(fabs(hJESPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(fabs(hJERPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(fabs(hMETUCPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(fabs(bTagPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(0.02*hCentral->GetBinContent(i+1),2) + pow(0.10*hCentral->GetBinContent(i+1),2) + pow(hCentral->GetBinError(i+1),2) + pow(0.03*hCentral->GetBinContent(i+1),2) + pow(0.2*h_QCD_iso->GetBinContent(i+1),2)); 
      //errorD[i] = sqrt(pow(fabs(hCentral->GetBinContent(i+1) - hJESMinus->GetBinContent(i+1)),2) + pow(fabs(hCentral->GetBinContent(i+1) - hJERMinus->GetBinContent(i+1)),2) + pow(fabs(hCentral->GetBinContent(i+1) - hMETUCMinus->GetBinContent(i+1)),2) + pow(fabs(hCentral->GetBinContent(i+1) - bTagMinus->GetBinContent(i+1)),2) + pow(0.02*hCentral->GetBinContent(i+1),2) + pow(0.10*hCentral->GetBinContent(i+1),2) + pow(hCentral->GetBinError(i+1),2) + pow(0.03*hCentral->GetBinContent(i+1),2) + pow(0.2*h_QCD_iso->GetBinContent(i+1),2));
      errorU[i] = hCentral->GetBinCenter(i+1);
      errorD[i] = hCentral->GetBinCenter(i+1);
      Xval[i]   = hCentral->GetBinCenter(i+1);
      XerrorU[i]= hCentral->GetBinWidth(i+1)/2.;
      XerrorD[i]= hCentral->GetBinWidth(i+1)/2.;
    } // for(int i=0; i<n1; i++)
  
  gr = new TGraphAsymmErrors(n1, Xval, Yval, XerrorD, XerrorU, errorD, errorU);
  
  return gr;
  
  delete [] Yval; delete [] errorU; delete [] errorD; delete [] XerrorU; delete [] XerrorD; delete [] Xval;

} // TGraphAsymmErrors *FULLGRAPH(TH1F *hCentral, TH1F *hUpper, TH1F *hLower)                                                                                                  

//TGraphAsymmErrors *RATIOGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus,  TH1F *hMETUCPlus, TH1F *hMETUCMinus, TH1F *bTagPlus, TH1F *bTagMinus, TH1F *h_QCD_iso) {
TGraphAsymmErrors *RATIOGRAPH(TH1F *hCentral, TH1F *h_QCD_iso) {
  
  TGraphAsymmErrors *gr;
  //int n1 = hJESPlus->GetNbinsX(); int n2 = hJESMinus->GetNbinsX();
  int n1 = 10; int n2 = 10;
  if(n1 != n2) cout<<"ATTENTION!!!!!  BinMismatch"<<endl;
  double *Yval, *errorU, *errorD, *XerrorU, *XerrorD, *Xval ;
  Yval = new double[n1]; errorU = new double[n1]; errorD = new double[n1];
  XerrorU=new double[n1]; XerrorD=new double[n1]; Xval=new double[n1];
  cout << "No. of bins= " << n1 << endl;
  
  for(int i=0; i<n1; i++)
    {
      Yval[i]   = 0.0;
      //errorU[i] = sqrt(pow(fabs(hJESPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(fabs(hJERPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(fabs(hMETUCPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(fabs(bTagPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(0.02*hCentral->GetBinContent(i+1),2) + pow(0.10*hCentral->GetBinContent(i+1),2) + pow(hCentral->GetBinError(i+1),2) + pow(0.03*hCentral->GetBinContent(i+1),2) + pow(0.2*h_QCD_iso->GetBinContent(i+1),2));
      //errorD[i] = sqrt(pow(fabs(hCentral->GetBinContent(i+1) - hJESMinus->GetBinContent(i+1)),2) + pow(fabs(hCentral->GetBinContent(i+1) - hJERMinus->GetBinContent(i+1)),2) + pow(fabs(hCentral->GetBinContent(i+1) - hMETUCMinus->GetBinContent(i+1)),2) + pow(fabs(hCentral->GetBinContent(i+1) - bTagMinus->GetBinContent(i+1)),2) + pow(0.02*hCentral->GetBinContent(i+1),2) + pow(0.10*hCentral->GetBinContent(i+1),2) + pow(hCentral->GetBinError(i+1),2) + pow(0.03*hCentral->GetBinContent(i+1),2) + pow(0.2*h_QCD_iso->GetBinContent(i+1),2));
      errorU[i] = errorU[i]/hCentral->GetBinContent(i+1);
      errorD[i] = errorD[i]/hCentral->GetBinContent(i+1);
      Xval[i]   = hCentral->GetBinCenter(i+1);
      XerrorU[i]= hCentral->GetBinWidth(i+1)/2.;
      XerrorD[i]= hCentral->GetBinWidth(i+1)/2.;
      
    } // for(int i=0; i<n1; i++)
  
  gr = new TGraphAsymmErrors(n1, Xval, Yval, XerrorD, XerrorU, errorD, errorU);
  return gr;
  delete [] Yval; delete [] errorU; delete [] errorD; delete [] XerrorU; delete [] XerrorD; delete [] Xval;
} // TGraphAsymmErrors *RATIOGRAPH(TH1F *hCentral, TH1F *hUpper, TH1F *hLower)


//example_stack(dir+"/pt_ele", "P_{T}[GeV]", 5, true, true, true, false, true, 25.0, 100.0);
void example_stack(TString histname, TString xaxis_title, int bin, bool log=false, bool drawdata=true, bool ratio=false, bool drawsignal=false, bool axisrange=false, double xmin=0, double xmax=10,bool label=false )
{
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();
  gStyle->SetOptStat(0);
  
  TString inFile("/afs/cern.ch/work/r/rverma/private/analysis/CMSSW_7_2_3/src/Analysis/stack/13TeV/outputDir/");

  c1 = new TCanvas();
  const float xpad[2] = {0.,1};
  //  const float ypad[4] = {0.,0.3,0.3,1.0};
  const float ypad[4] = {0.,0.2351916,0.2351916,0.98};
  if(ratio){
    c1->Divide(1,2);
    c1->cd(1);
    //gPad->SetBottomMargin(0);
    gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);
    //p->SetGridx();
    //p->SetGridy();
    gPad->SetLogy(log);
  }
  
  double scale_factor = 18.866686;
  
  TLegend* leg = new TLegend(0.7818792,0.6261504,0.9312081,0.9198861,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  //  TPaveText *pt = new TPaveText(0.15,0.92,0.9,0.99, "brNDC");
  TPaveText *pt = new TPaveText(0.15,0.94,0.93,0.97, "brNDC"); // good_v1
  pt->SetBorderSize(1);
  pt->SetFillColor(19);
  pt->SetFillStyle(0);
  pt->SetLineColor(0);
  pt->SetTextFont(132);

  THStack* MuptStack = new THStack("MuptStack","");
  TFile* f1 = TFile::Open(inFile+"HplusM120_MuMC_20170409_Ntuple_Merged_Anal.root");
  //TFile* f1 = TFile::Open(inFile+"diboson_ele_selection.root");
  if(f1 == 0) return;
  if(f1->IsZombie()){f1->Close(); return;}

  THStack* MuptStack = new THStack("MuptStack","");
  TH1F* h1_base = (f1->Get("base/"+histname))->Clone("h1_base");

  h1_base->Scale(scale_factor);  
  h1_base->SetFillColor(kYellow);
  if(axisrange){
    h1_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  MuptStack->Add(h1_base);
  
  TH1F* hMC = (TH1F*)h1_base->Clone("hMC");
  hMC->Reset();  
  
  TH1F* hSig = (TH1F*)h1_base->Clone("hSig"); 
  hSig->Reset(); 
  hSig->SetLineColor(kCyan);
  hSig->Add(h1_base);

  TFile* f2 = TFile::Open(inFile+"HplusM120_MuMC_20170409_Ntuple_Merged_Anal.root");
  //TFile* f2 = TFile::Open(inFile+"qcd_ele_selection.root");
  if(f2 == 0) return;
  if(f2->IsZombie()){f2->Close(); return;}
  TH1F* h2_base = (f2->Get("base/"+histname))->Clone("h2_base");

  h2_base->Scale(scale_factor);  
  h2_base->SetFillColor(kBlue);
  if(axisrange){
    h2_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }

  TString inFile_noniso("/afs/cern.ch/work/r/rverma/private/analysis/CMSSW_7_2_3/src/Analysis/stack/13TeV/outputDir/");
  //  double qcd_SF = 1.16388; double err_qcd_SF = 0.0757889; // after 1 electron selection:     1.16388 +/- 0.0757889 
  double qcd_SF = 1.32294; double err_qcd_SF = 0.00323672; // after 1 electron selection:     1.32294 +/- 0.00323672 in lowMET region 

  //TFile *ttbar_noniso = new TFile(inFile_noniso+"ttbar_ele_selection.root");
  TFile *ttbar_noniso = new TFile(inFile_noniso+"TTJets_MuMC_20170409_Ntuple_Merged_Anal.root");
  TH1F* h_ttbar_noniso = ((TH1F*)ttbar_noniso->Get("base/"+histname) )->Clone("h_ttbar_noniso");
  TH1F* h_ttbar_noniso_toppt =((TH1F*)ttbar_noniso->Get("base/AvTopPtWeight") )->Clone("h_ttbar_noniso_toppt");
  
  //TFile *wjet_noniso = new TFile(inFile_noniso+"wjet_ele_selection.root");
  TFile *ttbar_noniso = new TFile(inFile_noniso+"TTJets_MuMC_20170409_Ntuple_Merged_Anal.root");
  TH1F* h_wjet_noniso = ((TH1F*)wjet_noniso->Get("base/"+histname) )->Clone("h_wjet_noniso");

  //TFile *zjet_noniso = new TFile(inFile_noniso+"zjet_ele_selection.root");
  TFile *ttbar_noniso = new TFile(inFile_noniso+"TTJets_MuMC_20170409_Ntuple_Merged_Anal.root");
  TH1F* h_zjet_noniso = ((TH1F*)zjet_noniso->Get("base/"+histname) )->Clone("h_zjet_noniso");

  //TFile *stop_noniso = new TFile(inFile_noniso+"singletop_ele_selection.root");
  TFile *stop_noniso = new TFile(inFile_noniso+"TTJets_MuMC_20170409_Ntuple_Merged_Anal.root");
  TH1F* h_stop_noniso = ((TH1F*)stop_noniso->Get("base/"+histname) )->Clone("h_stop_noniso");
  
  //TFile *diboson_noniso = new TFile(inFile_noniso+"diboson_ele_selection.root");
  TFile *diboson_noniso = new TFile(inFile_noniso+"TTJets_MuMC_20170409_Ntuple_Merged_Anal.root");
  TH1F* h_diboson_noniso = ((TH1F*)diboson_noniso->Get("base/"+histname) )->Clone("h_diboson_noniso");

  //TFile *data_noniso = new TFile(inFile_noniso+"data_ele_selection.root");
  TFile *data_noniso = new TFile(inFile_noniso+"MuRunBv1_MuData_20170409_Ntuple_Merged_Anal.root");
  TH1F* h_data_noniso = ((TH1F*)data_noniso->Get("base/"+histname) )->Clone("h_data_noniso");
  
  h_ttbar_noniso->Scale(scale_factor/h_ttbar_noniso_toppt->GetBinContent(2));
  h_wjet_noniso->Scale(scale_factor);
  h_zjet_noniso->Scale(scale_factor);
  h_stop_noniso->Scale(scale_factor);
  h_diboson_noniso->Scale(scale_factor);

  // calculate total bkgs
  
  TH1F* h_noniso_TotalBkg = h_wjet_noniso->Clone("h_noniso_TotalBkg");
  h_noniso_TotalBkg->Reset();
  h_noniso_TotalBkg->Add(h_ttbar_noniso);
  h_noniso_TotalBkg->Add(h_wjet_noniso);
  h_noniso_TotalBkg->Add(h_zjet_noniso);
  h_noniso_TotalBkg->Add(h_stop_noniso);
  h_noniso_TotalBkg->Add(h_diboson_noniso);

  // define new noniso-qcd 
  TH1F* h_QCD_noniso = h_wjet_noniso->Clone("h_QCD_noniso");
  h_QCD_noniso->Reset();
  h_QCD_noniso->Add(h_data_noniso); 
  h_QCD_noniso->Add(h_noniso_TotalBkg,-1); // done (data - all other bkgs) in non-iso region
  
  // define new iso qcd
  TH1F* h_QCD_iso = h_wjet_noniso->Clone("h_QCD_iso");
  h_QCD_iso->Reset();
  h_QCD_iso->Add(h_QCD_noniso);
  h_QCD_iso->Scale(qcd_SF); // qcd_SF is the scale factor from non-iso region to iso region
  
  for(int ibin = 1; ibin < h_QCD_iso->GetNbinsX(); ibin++){
    //  cout << "BinContent:  "<< h_QCD_iso_new->GetBinContent(ibin) << endl;
    if(h_QCD_iso->GetBinContent(ibin) < 0.0){
      h_QCD_iso->SetBinContent(ibin,1.0e-5);
    }
  }
  

  h_QCD_iso->Rebin(bin);
  h_QCD_iso->SetFillColor(kBlue);
  if(axisrange){
    h_QCD_iso->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  MuptStack->Add(h_QCD_iso);
  hMC->Add(h_QCD_iso);
  hSig->Add(h_QCD_iso); 

  
  //TFile* f3 = TFile::Open(inFile+"zjet_ele_selection.root");
  TFile* f3 = TFile::Open(inFile+"HplusM120_MuMC_20170409_Ntuple_Merged_Anal.root");
  if(f3 == 0) return;
  if(f3->IsZombie()){f3->Close(); return;}
  TH1F* h3_base = (f3->Get("base/"+histname))->Clone("h3_base");

  h3_base->Scale(scale_factor);  
  h3_base->SetFillColor(kRed);
  if(axisrange){
    h3_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  MuptStack->Add(h3_base);
  hMC->Add(h3_base);
  hSig->Add(h3_base);


  //TFile* f4 = TFile::Open(inFile+"singletop_ele_selection.root");
  TFile* f4 = TFile::Open(inFile+"HplusM120_MuMC_20170409_Ntuple_Merged_Anal.root");
  if(f4 == 0) return;
  if(f4->IsZombie()){f4->Close(); return;}
  TH1F* h4_base = (f4->Get("base/"+histname))->Clone("h4_base");

  h4_base->Scale(scale_factor); 
  h4_base->SetFillColor(kGreen);
  if(axisrange){
    h4_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  MuptStack->Add(h4_base);
  hMC->Add(h4_base);
  hSig->Add(h4_base);


  //TFile* f5 = TFile::Open(inFile+"wjet_ele_selection.root");
  TFile* f5 = TFile::Open(inFile+"HplusM120_MuMC_20170409_Ntuple_Merged_Anal.root");
  if(f5 == 0) return;
  if(f5->IsZombie()){f5->Close(); return;}
  TH1F* h5_base = (f5->Get("base/"+histname))->Clone("h5_base");

  h5_base->Scale(scale_factor);  
  h5_base->SetFillColor(kMagenta+3);
  if(axisrange){
    h5_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  MuptStack->Add(h5_base);
  hMC->Add(h5_base);
  hSig->Add(h5_base);

  //TFile* f6 = TFile::Open(inFile+"ttbar_ele_selection.root");
  TFile* f6 = TFile::Open(inFile+"TTJets_MuMC_20170409_Ntuple_Merged_Anal.root");
  if(f6 == 0) return;
  if(f6->IsZombie()){f6->Close(); return;}
  TH1F* h_ttbar_toppt =((TH1F*)f6->Get("base/AvTopPtWeight") )->Clone("h_ttbar_toppt");
  TH1F* h6_base = (f6->Get("base/"+histname))->Clone("h6_base");

  h6_base->Scale(scale_factor/h_ttbar_toppt->GetBinContent(2)); 
  h6_base->Rebin(bin); 
  h6_base->SetFillColor(kCyan);
  if(axisrange){
    h6_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }

  if(histname.Contains("test")){
    cout << "producing plot for cutflow" << endl;
    h6_base->SetBinContent(2,h_ttbar_toppt->GetBinContent(2)*h6_base->GetBinContent(2));
    cout << "1 lepton in ttbar=  " << h6_base->GetBinContent(2) <<endl;
  }
  MuptStack->Add(h6_base);
  hMC->Add(h6_base);
 
  double br = 0.1; 
  hSig->Add(h6_base, (1.0-br)*(1.0-br));


  //TFile* f7 = TFile::Open(inFile+"wh_M_120_ele_selection.root");
  TFile* f7 = TFile::Open(inFile+"HplusM120_MuMC_20170409_Ntuple_Merged_Anal.root");
  if(f7 == 0) return;
  if(f7->IsZombie()){f7->Close(); return;}
  TH1F* h_wh120_toppt =((TH1F*)f7->Get("base/AvTopPtWeight") )->Clone("h_wh120_toppt");
  TH1F* h7_base = (f7->Get("base/"+histname))->Clone("h7_base");

  h7_base->Scale(scale_factor/h_wh120_toppt->GetBinContent(2));  
  h7_base->Rebin(bin); 
  
  if(axisrange){
    h7_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  
  if(histname.Contains("test")){
    cout << "producing plot for cutflow for signal" << endl;
    h7_base->SetBinContent(2,h_wh120_toppt->GetBinContent(2)*h7_base->GetBinContent(2));
  }

  hSig->Add(h7_base, 2*br*(1.0-br));
  h7_base->SetLineStyle(2);
  h7_base->SetLineWidth(2);

  hSig->SetLineColor(kMagenta+0);
  hSig->SetLineStyle(2);
  hSig->SetLineWidth(3);
  hSig->SetFillColor(0);
  
  ///////////////////// DATA ////////////////////////

  //TFile* data_ = TFile::Open(inFile+"data_ele_selection.root");
  TFile* data_ = TFile::Open(inFile+"MuRunBv1_MuData_20170409_Ntuple_Merged_Anal.root");
  if(data_ == 0) return;
  if(data_->IsZombie()){data_->Close(); return;}
  TH1F* data = (data_->Get("base/"+histname))->Clone("data");
  data->SetBinErrorOption(TH1::kPoisson); // added this line for Asymmetric Error Bars for Poisson Event Counts
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.9);
  data->Rebin(bin);
  if(axisrange){
    data->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  //  chargedHiggs signal
  data->SetFillColor(kBlack);

  if(label){
    TString steps[15] = {"","= 1 electron","#geq 4 jets","#slash{E}_{T} #geq 20GeV", "#geq 2 b-tagged jets","After KinFit","","","","","","","","",""};
    const size_t nsteps = sizeof(steps)/sizeof(TString);
    for(int istep=0; istep<nsteps; istep++ ){
      data->GetXaxis()->SetBinLabel(istep+1, steps[istep]);
      data->GetXaxis()->SetLabelSize(0.06);
      data->GetXaxis()->LabelsOption("u");
    }
  }
    
  ///////////////////// Combine ////////////////////////

  TGraphAsymmErrors *UncBand_;
  UncBand_ = FULLGRAPH(hMC, h_QCD_iso);
  //UncBand_ = FULLGRAPH(hMC, hMC_JESPlus, hMC_JESMinus, hMC_JERPlus, hMC_JERMinus, hMC_METUCPlus, hMC_METUCMinus, hMC_bTagPlus, hMC_bTagMinus, h_QCD_iso);
  UncBand_->SetFillColor(1);
  UncBand_->SetFillStyle(3001);
  
  // legend adding
  if(drawdata)leg->AddEntry(data,"Data","P"); // not for dijet mass
  leg->AddEntry(h1_base,"VV","F");
  leg->AddEntry(h2_base,"QCD","F");
  leg->AddEntry(h3_base,"Z+jets","F");
  leg->AddEntry(h4_base,"Single t","F");
  leg->AddEntry(h5_base,"W+jets","F");
  leg->AddEntry(h6_base,"t#bar{t}","F");
  leg->AddEntry(UncBand_,"Unc","F");
  if(drawsignal)leg->AddEntry(hSig,"Sig+bkgs","L"); // only for hSig histogram

  
  TPaveText *cct = new TPaveText(0.2013423,0.7754898,0.4010067,0.8962187,"brNDC");
  cct->SetFillColor(19);
  cct->SetFillStyle(0);
  cct->SetLineColor(0);
  cct->SetBorderSize(1);
  cct->AddText("M_{H^{+}} = 120 GeV");
  cct->AddText("Br(t#rightarrow H^{+}b) = 0.1");

  if(histname.Contains("mjj_kfit_Id_probfit1")) {
   TPaveText *ch = new TPaveText(0.201,0.70548,0.4010067,0.7562187,"brNDC");
  }else{
   TPaveText *ch = new TPaveText(0.453,0.7754898,0.6510067,0.8962187,"brNDC");
  }
  ch->SetFillColor(19);
  ch->SetFillStyle(0);
  ch->SetLineColor(0);
  ch->SetBorderSize(1);
  ch->AddText("e + jets");

  if(histname.Contains("mjj_kfit_Id_probfit1")) {
   pt->SetTextSize(0.049);
   TText *text = pt->AddText("CMS Preliminary,    #sqrt{s} = 8 TeV,    19.7 fb^{-1}");
   text->SetTextAlign(11);
  }else{
   pt->SetTextSize(0.059);
   TText *text = pt->AddText("CMS Preliminary,    #sqrt{s} = 8 TeV,    19.7 fb^{-1}");
   text->SetTextAlign(11);
  }

  //  data->SetAxisRange(1.0, (data->GetMaximum())*1000 ,"Y"); // good_v1
  if(histname.Contains("Pre_RelIso") || histname.Contains("cutflow")) {
   data->SetAxisRange(1.0, 5.0e9 ,"Y");
  }else{
   data->SetAxisRange(1.0, 1.0e6 ,"Y");
  }

  if(histname.Contains("mjj_kfit_Id_probfit1") && drawdata){
   data->SetAxisRange(0.0, 4000 ,"Y");
  }
  data->GetYaxis()->SetTitleOffset(1.45);
  data->GetYaxis()->SetTitle("Events");
  data->GetXaxis()->SetTitle(xaxis_title);
  //  UncBand_test->SetLineWidth(8);
  //  bool uncflag = true;
  if(drawdata){
    data->Draw("E"); // not for dijet mass
    MuptStack->Draw("HISTSAME"); // no histsame for dijet calculation
    MuptStack->GetXaxis()->SetTitle(xaxis_title);
    MuptStack->Draw("HISTSAME");
    bool uncflag = true;
    if(histname.Contains("pt_mu") || histname.Contains("eta_mu") || histname.Contains("pt_jet") || histname.Contains("eta_jet") || histname.Contains("final_multi_jet") || histname.Contains("final_pt_met") || histname.Contains("btag_jet") || histname.Contains("Final_RelIso") || histname.Contains("svMass_jID") || histname.Contains("nvtx") || histname.Contains("test1") || histname.Contains("multi_jet") || histname.Contains("pt_met") || histname.Contains("CSVM_count") )uncflag = false;
    string histpart = histname;
    if(histpart.find("BTag") != string::npos || histpart.find("KinFit") != string::npos ) uncflag = true;
    if(uncflag)    UncBand_->Draw("E2 same");
  }
  else{
    MuptStack->SetMinimum(1.0);
    MuptStack->SetMaximum(1.1*(MuptStack->GetMaximum()));
    //    stack->GetXaxis()->SetLimits(xmin, xmax);
    //    stack->Draw("H"); // make sure it's redrawn
    MuptStack->Draw("HIST");
    MuptStack->GetXaxis()->SetTitle(xaxis_title);
    MuptStack->GetYaxis()->SetTitle("Events");
    MuptStack->GetYaxis()->SetTitleOffset(1.45);
    MuptStack->Draw("HISTSAME");
    bool uncflag = true;
    if(histname.Contains("pt_mu") || histname.Contains("eta_mu") || histname.Contains("pt_jet") || histname.Contains("eta_jet") || histname.Contains("final_multi_jet") || histname.Contains("final_pt_met") || histname.Contains("btag_jet") || histname.Contains("Final_RelIso") || histname.Contains("svMass_jID") || histname.Contains("nvtx") || histname.Contains("test1") || histname.Contains("multi_jet") || histname.Contains("pt_met") || histname.Contains("CSVM_count"))uncflag = false;
    string histpart = histname;
    if(histpart.find("BTag") != string::npos || histpart.find("KinFit") != string::npos ) uncflag = true;
    if(uncflag)UncBand_->Draw("E2 same");
  }
  if(drawsignal)hSig->Draw("HISTSAME"); // only for hSig histogram 
  //  if(drawsignal)h7_base->Draw("HISTSAME"); // for sv_mass plot only

  if(drawdata)data->Draw("ESAME"); // not for dijet mass 
  MuptStack->SetMinimum(1.0);
  leg->Draw();
  pt->Draw();
  if(drawsignal) cct->Draw();
  ch->Draw();
  gPad->RedrawAxis();
  c1->Update();

  if(ratio){
    c1->cd(2);
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.3);
    gPad->SetGridy();
    gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
    //    gPad->SetPad(xpad[0],0.05,xpad[1],ypad[2]);

    TH1F *hRatio = (TH1F*)data->Clone("hRatio");
    hRatio->Reset();
    hRatio->Add(data);
    hRatio->Add(hMC, -1);
    hRatio->Divide(hMC);
    TGraphAsymmErrors *UncBand_Ratio;
    //UncBand_Ratio = RATIOGRAPH(hMC, hMC_JESPlus, hMC_JESMinus, hMC_JERPlus, hMC_JERMinus, hMC_METUCPlus, hMC_METUCMinus, hMC_bTagPlus, hMC_bTagMinus, h_QCD_iso);
    UncBand_Ratio = RATIOGRAPH(hMC, h_QCD_iso);
    UncBand_Ratio->SetFillColor(30);
    UncBand_Ratio->SetFillStyle(3001);

    // upto this ok then need to define one tgraph
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(1.0);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetRangeUser(-0.5, 0.5);
    //    hRatio->GetXaxis()->SetRangeUser(xmin, xmax);
    hRatio->GetXaxis()->SetTitle(xaxis_title);
    hRatio->GetYaxis()->SetTitleOffset(0.5);
    hRatio->GetYaxis()->SetTitle("#frac{Data-Bkg}{Bkg}");
    hRatio->GetYaxis()->CenterTitle();
    hRatio->GetYaxis()->SetTitleSize(0.1);
    hRatio->GetXaxis()->SetTitleSize(0.1);
    hRatio->GetXaxis()->SetLabelSize(0.17); // 0.1
    hRatio->GetXaxis()->LabelsOption("u"); // extra
    hRatio->GetYaxis()->SetLabelSize(0.06);
    hRatio->Draw("E"); // use "P" or "AP"
    bool uncflag = true;
    if(histname.Contains("pt_mu") || histname.Contains("eta_mu") || histname.Contains("pt_jet") || histname.Contains("eta_jet") || histname.Contains("final_multi_jet") || histname.Contains("final_pt_met") || histname.Contains("btag_jet") || histname.Contains("Final_RelIso") || histname.Contains("svMass_jID") || histname.Contains("nvtx") || histname.Contains("test1") || histname.Contains("multi_jet") || histname.Contains("pt_met") || histname.Contains("CSVM_count"))uncflag = false;    // "" this will be CSVM_count
    cout << "uncflag:   " << uncflag << endl;
    string histpart = histname; 
    if(histpart.find("BTag") != string::npos || histpart.find("KinFit") != string::npos ) uncflag = true;
    if(uncflag)UncBand_Ratio->Draw(" E2 same");
    hRatio->Draw("E same");
    c1->Update();
  }

  //TString outFile("plots/KineFit_ele_Sel_v1_test1_datadriven_from_data/base/");
  TString outFile(" ");
  outFile += histname;
  outFile += "_ele.png";
  c1->SaveAs(outFile);
  c1->Close();

}

void example_stack_all(){

   // final required plots:
  example_stack("mjj_kfit_Id_probfit1","M_{jj}[GeV]",1,false,true,true,true); // unblinded
  example_stack("cutflow","",1,true,true,true,true,true,1,5,true); // log, drawdata, ratio, drawsignal, axisrange, xmin, xmax, label // good_v1
  example_stack("Pre_RelIso","Relative Isolation",1,true,true, true,false,true,0,0.29);
  
}
void example_stack_kfit(TString dir="BTag"){
  example_stack(dir+"/pt_ele","P_{T}[GeV]",5,true, true,true, false, true,25.0,100.0); 
  example_stack(dir+"/eta_ele","#eta^{ele}",1,true,true,true,false,true,-2.5,2.5);
  example_stack(dir+"/pt_jet","P_{T}[GeV]",5,true,true,true,false, true,20.0,100.0);
  example_stack(dir+"/eta_jet","#eta^{jets}",1,true,true,true,false,true,-3.0,3.0);

  if(dir == "KinFit"){
    example_stack(dir+"/kfJet1_pt","P_{T}[GeV]",1,true,true,true,true, true,20.0,200.0);
    example_stack(dir+"/kfJet2_pt","P_{T}[GeV]",1,true,true,true,true, true,20.0,200.0);
    example_stack(dir+"/kfJet1_eta","#eta^{jets}",1,true,true,true,true,true,-3.0,3.0);
    example_stack(dir+"/kfJet2_eta","#eta^{jets}",1,true,true,true,true,true,-3.0,3.0);
  }

  example_stack(dir+"/final_multi_jet","N_{jets}",1,true,true,true,false,true,3,11);
  example_stack(dir+"/final_pt_met","E_{T}^{miss}[GeV]",10,true,true,true,false,true,15.0,100.0);
  example_stack(dir+"/wmt","m_{T}[GeV]",2,true,true,true,true,true,0,175);
  example_stack(dir+"/CSVM_count","btagged jet multiplicity",1,true,true, true,false,true,1,8);
  example_stack(dir+"/nvtx","N_{vertex}",1,true,true,true);
}
