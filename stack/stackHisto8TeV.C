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
TGraphAsymmErrors *FULLGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus,  TH1F *hMETUCPlus, TH1F *hMETUCMinus, TH1F *bTagPlus, TH1F *bTagMinus, TH1F *h_QCD_iso) {
  
  TGraphAsymmErrors *gr;
  int n1 = hJESPlus->GetNbinsX(); int n2 = hJESMinus->GetNbinsX();
  if(n1 != n2) cout<<"ATTENTION!!!!!  BinMismatch"<<endl;
  double *Yval, *errorU, *errorD, *XerrorU, *XerrorD, *Xval ;
  
  Yval = new double[n1]; errorU = new double[n1]; errorD = new double[n1];
  XerrorU=new double[n1]; XerrorD=new double[n1]; Xval=new double[n1];
  cout << "No. of bins= " << n1 << endl;
  for(int i=0; i<n1; i++)
    {
      Yval[i]   = hCentral->GetBinContent(i+1);
      
      errorU[i] = sqrt(pow(fabs(hJESPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(fabs(hJERPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(fabs(hMETUCPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(fabs(bTagPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(0.02*hCentral->GetBinContent(i+1),2) + pow(0.10*hCentral->GetBinContent(i+1),2) + pow(hCentral->GetBinError(i+1),2) + pow(0.03*hCentral->GetBinContent(i+1),2) + pow(0.2*h_QCD_iso->GetBinContent(i+1),2)); 
      errorD[i] = sqrt(pow(fabs(hCentral->GetBinContent(i+1) - hJESMinus->GetBinContent(i+1)),2) + pow(fabs(hCentral->GetBinContent(i+1) - hJERMinus->GetBinContent(i+1)),2) + pow(fabs(hCentral->GetBinContent(i+1) - hMETUCMinus->GetBinContent(i+1)),2) + pow(fabs(hCentral->GetBinContent(i+1) - bTagMinus->GetBinContent(i+1)),2) + pow(0.02*hCentral->GetBinContent(i+1),2) + pow(0.10*hCentral->GetBinContent(i+1),2) + pow(hCentral->GetBinError(i+1),2) + pow(0.03*hCentral->GetBinContent(i+1),2) + pow(0.2*h_QCD_iso->GetBinContent(i+1),2));

      Xval[i]   = hCentral->GetBinCenter(i+1);
      XerrorU[i]= hCentral->GetBinWidth(i+1)/2.;
      XerrorD[i]= hCentral->GetBinWidth(i+1)/2.;
      /*
      cout <<"x value= " << Xval[i]<<"Y central value=  "<< Yval[i] << "errorU[i]= "<< errorU[i] <<" & "<< "errorD[i]= "<< errorD[i] << endl;
      cout << "bin=  "<< i+1 << "central= "<< hCentral->GetBinContent(i+1) << "jesup=  " << hJESPlus->GetBinContent(i+1)<< "jesdown= " << hJESMinus->GetBinContent(i+1) << endl;
      cout << "bin=  "<< i+1 << "central= "<< hCentral->GetBinContent(i+1) << "metup=  " << hMETUCPlus->GetBinContent(i+1)<< "metdown= " << hMETUCMinus->GetBinContent(i+1) << endl;
      cout << "bin=  "<< i+1 << "central= "<< hCentral->GetBinContent(i+1) << "btagup=  " << bTagPlus->GetBinContent(i+1)<< "btagdown= " << bTagMinus->GetBinContent(i+1) << endl;
      */
    } // for(int i=0; i<n1; i++)
  
  gr = new TGraphAsymmErrors(n1, Xval, Yval, XerrorD, XerrorU, errorD, errorU);
  
  return gr;
  
  delete [] Yval; delete [] errorU; delete [] errorD; delete [] XerrorU; delete [] XerrorD; delete [] Xval;

} // TGraphAsymmErrors *FULLGRAPH(TH1F *hCentral, TH1F *hUpper, TH1F *hLower)                                                                                                  

TGraphAsymmErrors *RATIOGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus,  TH1F *hMETUCPlus, TH1F *hMETUCMinus, TH1F *bTagPlus, TH1F *bTagMinus, TH1F *h_QCD_iso) {
  
  TGraphAsymmErrors *gr;
  
  int n1 = hJESPlus->GetNbinsX(); int n2 = hJESMinus->GetNbinsX();
  
  if(n1 != n2) cout<<"ATTENTION!!!!!  BinMismatch"<<endl;
  
  double *Yval, *errorU, *errorD, *XerrorU, *XerrorD, *Xval ;
  
  Yval = new double[n1]; errorU = new double[n1]; errorD = new double[n1];
  XerrorU=new double[n1]; XerrorD=new double[n1]; Xval=new double[n1];
  cout << "No. of bins= " << n1 << endl;
  for(int i=0; i<n1; i++)
    {
      Yval[i]   = 0.0;
      errorU[i] = sqrt(pow(fabs(hJESPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(fabs(hJERPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(fabs(hMETUCPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(fabs(bTagPlus->GetBinContent(i+1) - hCentral->GetBinContent(i+1)),2) + pow(0.02*hCentral->GetBinContent(i+1),2) + pow(0.10*hCentral->GetBinContent(i+1),2) + pow(hCentral->GetBinError(i+1),2) + pow(0.03*hCentral->GetBinContent(i+1),2) + pow(0.2*h_QCD_iso->GetBinContent(i+1),2));
      errorD[i] = sqrt(pow(fabs(hCentral->GetBinContent(i+1) - hJESMinus->GetBinContent(i+1)),2) + pow(fabs(hCentral->GetBinContent(i+1) - hJERMinus->GetBinContent(i+1)),2) + pow(fabs(hCentral->GetBinContent(i+1) - hMETUCMinus->GetBinContent(i+1)),2) + pow(fabs(hCentral->GetBinContent(i+1) - bTagMinus->GetBinContent(i+1)),2) + pow(0.02*hCentral->GetBinContent(i+1),2) + pow(0.10*hCentral->GetBinContent(i+1),2) + pow(hCentral->GetBinError(i+1),2) + pow(0.03*hCentral->GetBinContent(i+1),2) + pow(0.2*h_QCD_iso->GetBinContent(i+1),2));
      errorU[i] = errorU[i]/hCentral->GetBinContent(i+1);
      errorD[i] = errorD[i]/hCentral->GetBinContent(i+1);
      Xval[i]   = hCentral->GetBinCenter(i+1);
      XerrorU[i]= hCentral->GetBinWidth(i+1)/2.;
      XerrorD[i]= hCentral->GetBinWidth(i+1)/2.;
      
      //      cout << " errorU[i]=  "<<  errorU[i] << "   errorD[i]=  " << errorD[i] << endl; 
      /*
	cout <<"x value= " << Xval[i]<<"Y central value=  "<< Yval[i] << "errorU[i]= "<< errorU[i] <<" & "<< "errorD[i]= "<< errorD[i] << endl;
	cout << "bin=  "<< i+1 << "central= "<< hCentral->GetBinContent(i+1) << "jesup=  " << hJESPlus->GetBinContent(i+1)<< "jesdown= " << hJESMinus->GetBinContent(i+1) << endl;
	cout << "bin=  "<< i+1 << "central= "<< hCentral->GetBinContent(i+1) << "metup=  " << hMETUCPlus->GetBinContent(i+1)<< "metdown= " << hMETUCMinus->GetBinContent(i+1) << endl;
	cout << "bin=  "<< i+1 << "central= "<< hCentral->GetBinContent(i+1) << "btagup=  " << bTagPlus->GetBinContent(i+1)<< "btagdown= " << bTagMinus->GetBinContent(i+1) << endl;
      */
    } // for(int i=0; i<n1; i++)
  
  gr = new TGraphAsymmErrors(n1, Xval, Yval, XerrorD, XerrorU, errorD, errorU);
  return gr;
  delete [] Yval; delete [] errorU; delete [] errorD; delete [] XerrorU; delete [] XerrorD; delete [] Xval;
} // TGraphAsymmErrors *RATIOGRAPH(TH1F *hCentral, TH1F *hUpper, TH1F *hLower)



void example_stack(TString histname, TString xaxis_title, int bin, bool log=false, bool drawdata=true, bool ratio=false, bool drawsignal=false, bool axisrange=false, double xmin=0, double xmax=10,bool label=false )
{
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();
  gStyle->SetOptStat(0);
  
  TString inFile("/afs/cern.ch/work/g/gkole/chargedHiggs/8TeV/ele_ch/SLC6/AfterApproved/8TeV/2M/MET_corr/Sel_app_v2/test1/merged/");

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
  //pt->SetTextSize(0.059);
  //TText *text = pt->AddText("CMS Preliminary,    #sqrt{s} = 8 TeV,    19.7 fb^{-1}");
  //text->SetTextAlign(11);
  //  pt->Draw();
  //  leg->SetHeader("#splitline{CMS Preliminary}{   @ #sqrt{s} = 7 TeV}");

  THStack* MuptStack = new THStack("MuptStack","");

  TFile* f1 = TFile::Open(inFile+"diboson_ele_selection.root");
  if(f1 == 0) return;
  if(f1->IsZombie()){f1->Close(); return;}

  THStack* MuptStack = new THStack("MuptStack","");
  TH1F* h1_base = (f1->Get("base/"+histname))->Clone("h1_base");
  TH1F* h1_JESPlus = ((TH1F*)f1->Get("JESPlus/"+histname) )->Clone("h1_JESPlus");
  TH1F* h1_JESMinus = ((TH1F*)f1->Get("JESMinus/"+histname) )->Clone("h1_JESMinus");
  TH1F* h1_JERPlus = ((TH1F*)f1->Get("JERPlus/"+histname) )->Clone("h1_JERPlus");
  TH1F* h1_JERMinus = ((TH1F*)f1->Get("JERMinus/"+histname) )->Clone("h1_JERMinus");
  TH1F* h1_METUCPlus = ((TH1F*)f1->Get("METUCPlus/"+histname) )->Clone("h1_METUCPlus");
  TH1F* h1_METUCMinus = ((TH1F*)f1->Get("METUCMinus/"+histname) )->Clone("h1_METUCMinus");
  TH1F* h1_bTagPlus = ((TH1F*)f1->Get("bTagPlus/"+histname) )->Clone("h1_bTagPlus");
  TH1F* h1_bTagMinus = ((TH1F*)f1->Get("bTagMinus/"+histname) )->Clone("h1_bTagMinus");

  h1_base->Scale(scale_factor);  h1_JESPlus->Scale(scale_factor);  h1_JESMinus->Scale(scale_factor);h1_JERPlus->Scale(scale_factor);  h1_JERMinus->Scale(scale_factor);
  h1_METUCPlus->Scale(scale_factor); h1_METUCMinus->Scale(scale_factor); h1_bTagPlus->Scale(scale_factor); h1_bTagMinus->Scale(scale_factor);
  h1_base->Rebin(bin); h1_JESPlus->Rebin(bin);  h1_JESMinus->Rebin(bin);  h1_JERPlus->Rebin(bin);  h1_JERMinus->Rebin(bin);
  h1_METUCPlus->Rebin(bin); h1_METUCMinus->Rebin(bin); h1_bTagPlus->Rebin(bin); h1_bTagMinus->Rebin(bin);

  h1_base->SetFillColor(kYellow);
  if(axisrange){
    h1_base->GetXaxis()->SetRangeUser(xmin,xmax);
    h1_JESPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h1_JESMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h1_JERPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h1_JERMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h1_METUCPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h1_METUCMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h1_bTagPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h1_bTagMinus->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  MuptStack->Add(h1_base);
  

  TH1F* hMC = (TH1F*)h1_base->Clone("hMC");
  TH1F* hMC_JESPlus = (TH1F*)h1_base->Clone("hMC_JESPlus");
  TH1F* hMC_JESMinus = (TH1F*)h1_base->Clone("hMC_JESMinus");
  TH1F* hMC_JERPlus = (TH1F*)h1_base->Clone("hMC_JERPlus");
  TH1F* hMC_JERMinus = (TH1F*)h1_base->Clone("hMC_JERMinus");
  TH1F* hMC_METUCPlus = (TH1F*)h1_base->Clone("hMC_METUCPlus");
  TH1F* hMC_METUCMinus = (TH1F*)h1_base->Clone("hMC_METUCMinus");
  TH1F* hMC_bTagPlus = (TH1F*)h1_base->Clone("hMC_bTagPlus");
  TH1F* hMC_bTagMinus = (TH1F*)h1_base->Clone("hMC_bTagMinus");

  hMC->Reset();  hMC_JESPlus->Reset();  hMC_JESMinus->Reset(); hMC_JERPlus->Reset();  hMC_JERMinus->Reset(); hMC_METUCPlus->Reset(); hMC_METUCMinus->Reset(); hMC_bTagPlus->Reset(); hMC_bTagMinus->Reset();
  hMC->Add(h1_base);
  hMC_JESPlus->Add(h1_JESPlus); hMC_JERPlus->Add(h1_JERPlus); hMC_METUCPlus->Add(h1_METUCPlus); hMC_bTagPlus->Add(h1_bTagPlus); // only + line
  hMC_JESMinus->Add(h1_JESMinus); hMC_JERMinus->Add(h1_JERMinus); hMC_METUCMinus->Add(h1_METUCMinus); hMC_bTagMinus->Add(h1_bTagMinus); // only - line

  TH1F* hSig = (TH1F*)h1_base->Clone("hSig"); 
  hSig->Reset(); 
  hSig->SetLineColor(kCyan);
  hSig->Add(h1_base);

  TFile* f2 = TFile::Open(inFile+"qcd_ele_selection.root");
  if(f2 == 0) return;
  if(f2->IsZombie()){f2->Close(); return;}
  TH1F* h2_base = (f2->Get("base/"+histname))->Clone("h2_base");
  TH1F* h2_JESPlus = ((TH1F*)f2->Get("JESPlus/"+histname) )->Clone("h2_JESPlus");
  TH1F* h2_JESMinus = ((TH1F*)f2->Get("JESMinus/"+histname) )->Clone("h2_JESMinus");
  TH1F* h2_JERPlus = ((TH1F*)f2->Get("JERPlus/"+histname) )->Clone("h2_JERPlus");
  TH1F* h2_JERMinus = ((TH1F*)f2->Get("JERMinus/"+histname) )->Clone("h2_JERMinus");
  TH1F* h2_METUCPlus = ((TH1F*)f2->Get("METUCPlus/"+histname) )->Clone("h2_METUCPlus");
  TH1F* h2_METUCMinus = ((TH1F*)f2->Get("METUCMinus/"+histname) )->Clone("h2_METUCMinus");
  TH1F* h2_bTagPlus = ((TH1F*)f2->Get("bTagPlus/"+histname) )->Clone("h2_bTagPlus");
  TH1F* h2_bTagMinus = ((TH1F*)f2->Get("bTagMinus/"+histname) )->Clone("h2_bTagMinus");

  h2_base->Scale(scale_factor);  h2_JESPlus->Scale(scale_factor);  h2_JESMinus->Scale(scale_factor);h2_JERPlus->Scale(scale_factor);  h2_JERMinus->Scale(scale_factor);
  h2_METUCPlus->Scale(scale_factor); h2_METUCMinus->Scale(scale_factor); h2_bTagPlus->Scale(scale_factor); h2_bTagMinus->Scale(scale_factor);
  h2_base->Rebin(bin); h2_JESPlus->Rebin(bin);  h2_JESMinus->Rebin(bin); h2_JERPlus->Rebin(bin);  h2_JERMinus->Rebin(bin);
  h2_METUCPlus->Rebin(bin); h2_METUCMinus->Rebin(bin); h2_bTagPlus->Rebin(bin); h2_bTagMinus->Rebin(bin);

  h2_base->SetFillColor(kBlue);
  if(axisrange){
    h2_base->GetXaxis()->SetRangeUser(xmin,xmax);
    h2_JESPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h2_JESMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h2_JERPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h2_JERMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h2_METUCPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h2_METUCMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h2_bTagPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h2_bTagMinus->GetXaxis()->SetRangeUser(xmin,xmax);

  }

  // Here qcd is data-driven
  
  // MuptStack->Add(h2_base);
  // hMC->Add(h2_base);
  // hMC_JESPlus->Add(h2_JESPlus); hMC_JERPlus->Add(h2_JERPlus); hMC_METUCPlus->Add(h2_METUCPlus); hMC_bTagPlus->Add(h2_bTagPlus); // only + line                                                               
  // hMC_JESMinus->Add(h2_JESMinus); hMC_JERMinus->Add(h2_JERMinus); hMC_METUCMinus->Add(h2_METUCMinus); hMC_bTagMinus->Add(h2_bTagMinus); // only - line                                                         
  // hSig->Add(h2_base);
  
  TString inFile_noniso("/afs/cern.ch/work/g/gkole/chargedHiggs/8TeV/ele_ch/SLC6/AfterApproved/8TeV/2M/MET_corr/Sel_app_v2/test1_antiiso/merged/");
  //  double qcd_SF = 1.16388; double err_qcd_SF = 0.0757889; // after 1 electron selection:     1.16388 +/- 0.0757889 
  double qcd_SF = 1.32294; double err_qcd_SF = 0.00323672; // after 1 electron selection:     1.32294 +/- 0.00323672 in lowMET region 

  TFile *ttbar_noniso = new TFile(inFile_noniso+"ttbar_ele_selection.root");
  TH1F* h_ttbar_noniso = ((TH1F*)ttbar_noniso->Get("base/"+histname) )->Clone("h_ttbar_noniso");
  TH1F* h_ttbar_noniso_toppt =((TH1F*)ttbar_noniso->Get("base/AvTopPtWeight") )->Clone("h_ttbar_noniso_toppt");
  
  TFile *wjet_noniso = new TFile(inFile_noniso+"wjet_ele_selection.root");
  TH1F* h_wjet_noniso = ((TH1F*)wjet_noniso->Get("base/"+histname) )->Clone("h_wjet_noniso");

  TFile *zjet_noniso = new TFile(inFile_noniso+"zjet_ele_selection.root");
  TH1F* h_zjet_noniso = ((TH1F*)zjet_noniso->Get("base/"+histname) )->Clone("h_zjet_noniso");

  TFile *stop_noniso = new TFile(inFile_noniso+"singletop_ele_selection.root");
  TH1F* h_stop_noniso = ((TH1F*)stop_noniso->Get("base/"+histname) )->Clone("h_stop_noniso");
  
  TFile *diboson_noniso = new TFile(inFile_noniso+"diboson_ele_selection.root");
  TH1F* h_diboson_noniso = ((TH1F*)diboson_noniso->Get("base/"+histname) )->Clone("h_diboson_noniso");

  TFile *data_noniso = new TFile(inFile_noniso+"data_ele_selection.root");
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
  hMC_JESPlus->Add(h_QCD_iso); hMC_JERPlus->Add(h_QCD_iso); hMC_METUCPlus->Add(h_QCD_iso); hMC_bTagPlus->Add(h_QCD_iso); // only + line
  hMC_JESMinus->Add(h_QCD_iso); hMC_JERMinus->Add(h_QCD_iso); hMC_METUCMinus->Add(h_QCD_iso); hMC_bTagMinus->Add(h_QCD_iso); // only - line
  
  hSig->Add(h_QCD_iso); 

  
  
  TFile* f3 = TFile::Open(inFile+"zjet_ele_selection.root");
  if(f3 == 0) return;
  if(f3->IsZombie()){f3->Close(); return;}
  TH1F* h3_base = (f3->Get("base/"+histname))->Clone("h3_base");
  TH1F* h3_JESPlus = ((TH1F*)f3->Get("JESPlus/"+histname) )->Clone("h3_JESPlus");
  TH1F* h3_JESMinus = ((TH1F*)f3->Get("JESMinus/"+histname) )->Clone("h3_JESMinus");
  TH1F* h3_JERPlus = ((TH1F*)f3->Get("JERPlus/"+histname) )->Clone("h3_JERPlus");
  TH1F* h3_JERMinus = ((TH1F*)f3->Get("JERMinus/"+histname) )->Clone("h3_JERMinus");
  TH1F* h3_METUCPlus = ((TH1F*)f3->Get("METUCPlus/"+histname) )->Clone("h3_METUCPlus");
  TH1F* h3_METUCMinus = ((TH1F*)f3->Get("METUCMinus/"+histname) )->Clone("h3_METUCMinus");
  TH1F* h3_bTagPlus = ((TH1F*)f3->Get("bTagPlus/"+histname) )->Clone("h3_bTagPlus");
  TH1F* h3_bTagMinus = ((TH1F*)f3->Get("bTagMinus/"+histname) )->Clone("h3_bTagMinus");

  h3_base->Scale(scale_factor);  h3_JESPlus->Scale(scale_factor);  h3_JESMinus->Scale(scale_factor); h3_JERPlus->Scale(scale_factor);  h3_JERMinus->Scale(scale_factor);
  h3_METUCPlus->Scale(scale_factor); h3_METUCMinus->Scale(scale_factor); h3_bTagPlus->Scale(scale_factor); h3_bTagMinus->Scale(scale_factor);
  h3_base->Rebin(bin); h3_JESPlus->Rebin(bin);  h3_JESMinus->Rebin(bin); h3_JERPlus->Rebin(bin);  h3_JERMinus->Rebin(bin);
  h3_METUCPlus->Rebin(bin); h3_METUCMinus->Rebin(bin); h3_bTagPlus->Rebin(bin); h3_bTagMinus->Rebin(bin);

  h3_base->SetFillColor(kRed);
  if(axisrange){
    h3_base->GetXaxis()->SetRangeUser(xmin,xmax);
    h3_JESPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h3_JESMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h3_JERPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h3_JERMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h3_METUCPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h3_METUCMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h3_bTagPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h3_bTagMinus->GetXaxis()->SetRangeUser(xmin,xmax);

  }
  MuptStack->Add(h3_base);
  hMC->Add(h3_base);
  hMC_JESPlus->Add(h3_JESPlus); hMC_JERPlus->Add(h3_JERPlus); hMC_METUCPlus->Add(h3_METUCPlus); hMC_bTagPlus->Add(h3_bTagPlus); // only + line
  hMC_JESMinus->Add(h3_JESMinus); hMC_JERMinus->Add(h3_JERMinus); hMC_METUCMinus->Add(h3_METUCMinus); hMC_bTagMinus->Add(h3_bTagMinus); // only - line
  hSig->Add(h3_base);


  TFile* f4 = TFile::Open(inFile+"singletop_ele_selection.root");
  if(f4 == 0) return;
  if(f4->IsZombie()){f4->Close(); return;}
  TH1F* h4_base = (f4->Get("base/"+histname))->Clone("h4_base");
  TH1F* h4_JESPlus = ((TH1F*)f4->Get("JESPlus/"+histname) )->Clone("h4_JESPlus");
  TH1F* h4_JESMinus = ((TH1F*)f4->Get("JESMinus/"+histname) )->Clone("h4_JESMinus");
  TH1F* h4_JERPlus = ((TH1F*)f4->Get("JERPlus/"+histname) )->Clone("h4_JERPlus");
  TH1F* h4_JERMinus = ((TH1F*)f4->Get("JERMinus/"+histname) )->Clone("h4_JERMinus");
  TH1F* h4_METUCPlus = ((TH1F*)f4->Get("METUCPlus/"+histname) )->Clone("h4_METUCPlus");
  TH1F* h4_METUCMinus = ((TH1F*)f4->Get("METUCMinus/"+histname) )->Clone("h4_METUCMinus");
  TH1F* h4_bTagPlus = ((TH1F*)f4->Get("bTagPlus/"+histname) )->Clone("h4_bTagPlus");
  TH1F* h4_bTagMinus = ((TH1F*)f4->Get("bTagMinus/"+histname) )->Clone("h4_bTagMinus");

  h4_base->Scale(scale_factor);  h4_JESPlus->Scale(scale_factor);  h4_JESMinus->Scale(scale_factor); h4_JERPlus->Scale(scale_factor);  h4_JERMinus->Scale(scale_factor);
  h4_METUCPlus->Scale(scale_factor); h4_METUCMinus->Scale(scale_factor); h4_bTagPlus->Scale(scale_factor); h4_bTagMinus->Scale(scale_factor);
  h4_base->Rebin(bin); h4_JESPlus->Rebin(bin);  h4_JESMinus->Rebin(bin); h4_JERPlus->Rebin(bin);  h4_JERMinus->Rebin(bin);
  h4_METUCPlus->Rebin(bin); h4_METUCMinus->Rebin(bin); h4_bTagPlus->Rebin(bin); h4_bTagMinus->Rebin(bin);

  h4_base->SetFillColor(kGreen);
  if(axisrange){
    h4_base->GetXaxis()->SetRangeUser(xmin,xmax);
    h4_JESPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h4_JESMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h4_JERPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h4_JERMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h4_METUCPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h4_METUCMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h4_bTagPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h4_bTagMinus->GetXaxis()->SetRangeUser(xmin,xmax);

  }
  MuptStack->Add(h4_base);
  hMC->Add(h4_base);
  hMC_JESPlus->Add(h4_JESPlus); hMC_JERPlus->Add(h4_JERPlus); hMC_METUCPlus->Add(h4_METUCPlus); hMC_bTagPlus->Add(h4_bTagPlus); // only + line
  hMC_JESMinus->Add(h4_JESMinus); hMC_JERMinus->Add(h4_JERMinus); hMC_METUCMinus->Add(h4_METUCMinus); hMC_bTagMinus->Add(h4_bTagMinus); // only - line
  
  hSig->Add(h4_base);


  TFile* f5 = TFile::Open(inFile+"wjet_ele_selection.root");
  if(f5 == 0) return;
  if(f5->IsZombie()){f5->Close(); return;}
  TH1F* h5_base = (f5->Get("base/"+histname))->Clone("h5_base");
  TH1F* h5_JESPlus = ((TH1F*)f5->Get("JESPlus/"+histname) )->Clone("h5_JESPlus");
  TH1F* h5_JESMinus = ((TH1F*)f5->Get("JESMinus/"+histname) )->Clone("h5_JESMinus");
  TH1F* h5_JERPlus = ((TH1F*)f5->Get("JERPlus/"+histname) )->Clone("h5_JERPlus");
  TH1F* h5_JERMinus = ((TH1F*)f5->Get("JERMinus/"+histname) )->Clone("h5_JERMinus");
  TH1F* h5_METUCPlus = ((TH1F*)f5->Get("METUCPlus/"+histname) )->Clone("h5_METUCPlus");
  TH1F* h5_METUCMinus = ((TH1F*)f5->Get("METUCMinus/"+histname) )->Clone("h5_METUCMinus");
  TH1F* h5_bTagPlus = ((TH1F*)f5->Get("bTagPlus/"+histname) )->Clone("h5_bTagPlus");
  TH1F* h5_bTagMinus = ((TH1F*)f5->Get("bTagMinus/"+histname) )->Clone("h5_bTagMinus");

  h5_base->Scale(scale_factor);  h5_JESPlus->Scale(scale_factor);  h5_JESMinus->Scale(scale_factor); h5_JERPlus->Scale(scale_factor);  h5_JERMinus->Scale(scale_factor);
  h5_METUCPlus->Scale(scale_factor); h5_METUCMinus->Scale(scale_factor); h5_bTagPlus->Scale(scale_factor); h5_bTagMinus->Scale(scale_factor);
  h5_base->Rebin(bin); h5_JESPlus->Rebin(bin);  h5_JESMinus->Rebin(bin); h5_JERPlus->Rebin(bin);  h5_JERMinus->Rebin(bin);
  h5_METUCPlus->Rebin(bin); h5_METUCMinus->Rebin(bin); h5_bTagPlus->Rebin(bin); h5_bTagMinus->Rebin(bin);

  h5_base->SetFillColor(kMagenta+3);
  if(axisrange){
    h5_base->GetXaxis()->SetRangeUser(xmin,xmax);
    h5_JESPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h5_JESMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h5_JERPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h5_JERMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h5_METUCPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h5_METUCMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h5_bTagPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h5_bTagMinus->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  MuptStack->Add(h5_base);
  hMC->Add(h5_base);
  hMC_JESPlus->Add(h5_JESPlus); hMC_JERPlus->Add(h5_JERPlus);  hMC_METUCPlus->Add(h5_METUCPlus); hMC_bTagPlus->Add(h5_bTagPlus); // only + line
  hMC_JESMinus->Add(h5_JESMinus); hMC_JERMinus->Add(h5_JERMinus); hMC_METUCMinus->Add(h5_METUCMinus); hMC_bTagMinus->Add(h5_bTagMinus); // only - line
  hSig->Add(h5_base);
  //  cout << "hMC_JESPlus get bin after 1 mu just only wjet=  "<< h6_JESPlus->GetBinContent(2) << "same for hMC  " << h6_base->GetBinContent(2) << "same for hMC_JESMinus=   "<<h6_JESMinus->GetBinContent(2) <<  endl;

  //  cout << "hMC_JESPlus get bin after 1 mu after full addition=  "<< hMC_JESPlus->GetBinContent(2) << "same for hMC  " << hMC->GetBinContent(2) << "same for hMC_JESMinus=    "  <<hMC_JESMinus->GetBinContent(2) <<  endl;


  TFile* f6 = TFile::Open(inFile+"ttbar_ele_selection.root");
  if(f6 == 0) return;
  if(f6->IsZombie()){f6->Close(); return;}
  TH1F* h_ttbar_toppt =((TH1F*)f6->Get("base/AvTopPtWeight") )->Clone("h_ttbar_toppt");
  TH1F* h6_base = (f6->Get("base/"+histname))->Clone("h6_base");
  TH1F* h6_JESPlus = ((TH1F*)f6->Get("JESPlus/"+histname) )->Clone("h6_JESPlus");
  TH1F* h6_JESMinus = ((TH1F*)f6->Get("JESMinus/"+histname) )->Clone("h6_JESMinus");
  TH1F* h6_JERPlus = ((TH1F*)f6->Get("JERPlus/"+histname) )->Clone("h6_JERPlus");
  TH1F* h6_JERMinus = ((TH1F*)f6->Get("JERMinus/"+histname) )->Clone("h6_JERMinus");
  TH1F* h6_METUCPlus = ((TH1F*)f6->Get("METUCPlus/"+histname) )->Clone("h6_METUCPlus");
  TH1F* h6_METUCMinus = ((TH1F*)f6->Get("METUCMinus/"+histname) )->Clone("h6_METUCMinus");
  TH1F* h6_bTagPlus = ((TH1F*)f6->Get("bTagPlus/"+histname) )->Clone("h6_bTagPlus");
  TH1F* h6_bTagMinus = ((TH1F*)f6->Get("bTagMinus/"+histname) )->Clone("h6_bTagMinus");

  h6_base->Scale(scale_factor/h_ttbar_toppt->GetBinContent(2));  h6_JESPlus->Scale(scale_factor/h_ttbar_toppt->GetBinContent(2));  h6_JESMinus->Scale(scale_factor/h_ttbar_toppt->GetBinContent(2));
  h6_JERPlus->Scale(scale_factor/h_ttbar_toppt->GetBinContent(2));  h6_JERMinus->Scale(scale_factor/h_ttbar_toppt->GetBinContent(2));
  h6_METUCPlus->Scale(scale_factor/h_ttbar_toppt->GetBinContent(2)); h6_METUCMinus->Scale(scale_factor/h_ttbar_toppt->GetBinContent(2)); h6_bTagPlus->Scale(scale_factor/h_ttbar_toppt->GetBinContent(2)); h6_bTagMinus->Scale(scale_factor/h_ttbar_toppt->GetBinContent(2));
  h6_base->Rebin(bin); h6_JESPlus->Rebin(bin);  h6_JESMinus->Rebin(bin); h6_JERPlus->Rebin(bin);  h6_JERMinus->Rebin(bin);
  h6_METUCPlus->Rebin(bin); h6_METUCMinus->Rebin(bin); h6_bTagPlus->Rebin(bin); h6_bTagMinus->Rebin(bin);
  
  h6_base->SetFillColor(kCyan);
  if(axisrange){
    h6_base->GetXaxis()->SetRangeUser(xmin,xmax);
    h6_JESPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h6_JESMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h6_JERPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h6_JERMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h6_METUCPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h6_METUCMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h6_bTagPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h6_bTagMinus->GetXaxis()->SetRangeUser(xmin,xmax);
  
  }
  if(histname.Contains("test")){
    cout << "producing plot for cutflow" << endl;
    h6_base->SetBinContent(2,h_ttbar_toppt->GetBinContent(2)*h6_base->GetBinContent(2));
    h6_JESPlus->SetBinContent(2,h_ttbar_toppt->GetBinContent(2)*h6_base->GetBinContent(2));
    h6_JESMinus->SetBinContent(2,h_ttbar_toppt->GetBinContent(2)*h6_JESMinus->GetBinContent(2));
    h6_METUCPlus->SetBinContent(2,h_ttbar_toppt->GetBinContent(2)*h6_METUCPlus->GetBinContent(2));
    h6_METUCMinus->SetBinContent(2,h_ttbar_toppt->GetBinContent(2)*h6_METUCMinus->GetBinContent(2));
    h6_bTagPlus->SetBinContent(2,h_ttbar_toppt->GetBinContent(2)*h6_bTagPlus->GetBinContent(2));
    h6_bTagMinus->SetBinContent(2,h_ttbar_toppt->GetBinContent(2)*h6_bTagMinus->GetBinContent(2));
    
    cout << "1 lepton in ttbar=  " << h6_base->GetBinContent(2) <<endl;
  }
  MuptStack->Add(h6_base);
  hMC->Add(h6_base);
  hMC_JESPlus->Add(h6_JESPlus);  hMC_JERPlus->Add(h6_JERPlus); hMC_METUCPlus->Add(h6_METUCPlus); hMC_bTagPlus->Add(h6_bTagPlus); // only + line
  hMC_JESMinus->Add(h6_JESMinus); hMC_JERMinus->Add(h6_JERMinus); hMC_METUCMinus->Add(h6_METUCMinus); hMC_bTagMinus->Add(h6_bTagMinus); // only - line
 
  double br = 0.1; 
  hSig->Add(h6_base, (1.0-br)*(1.0-br));


  TFile* f7 = TFile::Open(inFile+"wh_M_120_ele_selection.root");
  if(f7 == 0) return;
  if(f7->IsZombie()){f7->Close(); return;}
  TH1F* h_wh120_toppt =((TH1F*)f7->Get("base/AvTopPtWeight") )->Clone("h_wh120_toppt");
  TH1F* h7_base = (f7->Get("base/"+histname))->Clone("h7_base");
  TH1F* h7_JESPlus = ((TH1F*)f7->Get("JESPlus/"+histname) )->Clone("h7_JESPlus");
  TH1F* h7_JESMinus = ((TH1F*)f7->Get("JESMinus/"+histname) )->Clone("h7_JESMinus");
  TH1F* h7_JERPlus = ((TH1F*)f7->Get("JERPlus/"+histname) )->Clone("h7_JERPlus");
  TH1F* h7_JERMinus = ((TH1F*)f7->Get("JERMinus/"+histname) )->Clone("h7_JERMinus");
  TH1F* h7_METUCPlus = ((TH1F*)f7->Get("METUCPlus/"+histname) )->Clone("h7_METUCPlus");
  TH1F* h7_METUCMinus = ((TH1F*)f7->Get("METUCMinus/"+histname) )->Clone("h7_METUCMinus");
  TH1F* h7_bTagPlus = ((TH1F*)f7->Get("bTagPlus/"+histname) )->Clone("h7_bTagPlus");
  TH1F* h7_bTagMinus = ((TH1F*)f7->Get("bTagMinus/"+histname) )->Clone("h7_bTagMinus");

  h7_base->Scale(scale_factor/h_wh120_toppt->GetBinContent(2));  h7_JESPlus->Scale(scale_factor/h_wh120_toppt->GetBinContent(2));  h7_JESMinus->Scale(scale_factor/h_wh120_toppt->GetBinContent(2));
  h7_JERPlus->Scale(scale_factor/h_wh120_toppt->GetBinContent(2));  h7_JERMinus->Scale(scale_factor/h_wh120_toppt->GetBinContent(2));
  h7_METUCPlus->Scale(scale_factor/h_wh120_toppt->GetBinContent(2)); h7_METUCMinus->Scale(scale_factor/h_wh120_toppt->GetBinContent(2)); h7_bTagPlus->Scale(scale_factor/h_wh120_toppt->GetBinContent(2)); h7_bTagMinus->Scale(scale_factor/h_wh120_toppt->GetBinContent(2));
  h7_base->Rebin(bin); h7_JESPlus->Rebin(bin);  h7_JESMinus->Rebin(bin); h7_JERPlus->Rebin(bin);  h7_JERMinus->Rebin(bin);
  h7_METUCPlus->Rebin(bin); h7_METUCMinus->Rebin(bin); h7_bTagPlus->Rebin(bin); h7_bTagMinus->Rebin(bin);

  
  if(axisrange){
    h7_base->GetXaxis()->SetRangeUser(xmin,xmax);
    h7_JESPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h7_JESMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h7_JERPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h7_JERMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h7_METUCPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h7_METUCMinus->GetXaxis()->SetRangeUser(xmin,xmax);
    h7_bTagPlus->GetXaxis()->SetRangeUser(xmin,xmax);
    h7_bTagMinus->GetXaxis()->SetRangeUser(xmin,xmax);

  }
  
  if(histname.Contains("test")){
    cout << "producing plot for cutflow for signal" << endl;
    h7_base->SetBinContent(2,h_wh120_toppt->GetBinContent(2)*h7_base->GetBinContent(2));
    h7_JESPlus->SetBinContent(2,h_wh120_toppt->GetBinContent(2)*h7_JESPlus->GetBinContent(2));
    h7_JESMinus->SetBinContent(2,h_wh120_toppt->GetBinContent(2)*h7_JESMinus->GetBinContent(2));
    h7_JERPlus->SetBinContent(2,h_wh120_toppt->GetBinContent(2)*h7_JERPlus->GetBinContent(2));
    h7_JERMinus->SetBinContent(2,h_wh120_toppt->GetBinContent(2)*h7_JERMinus->GetBinContent(2));
    h7_METUCPlus->SetBinContent(2,h_wh120_toppt->GetBinContent(2)*h7_METUCPlus->GetBinContent(2));
    h7_METUCMinus->SetBinContent(2,h_wh120_toppt->GetBinContent(2)*h7_METUCMinus->GetBinContent(2));
    h7_bTagPlus->SetBinContent(2,h_wh120_toppt->GetBinContent(2)*h7_bTagPlus->GetBinContent(2));
    h7_bTagMinus->SetBinContent(2,h_wh120_toppt->GetBinContent(2)*h7_bTagMinus->GetBinContent(2));
  }

  hSig->Add(h7_base, 2*br*(1.0-br));
  h7_base->SetLineStyle(2);
  h7_base->SetLineWidth(2);

  hSig->SetLineColor(kMagenta+0);
  hSig->SetLineStyle(2);
  hSig->SetLineWidth(3);
  hSig->SetFillColor(0);
  

  TFile* data_ = TFile::Open(inFile+"data_ele_selection.root");
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
    
  TGraphAsymmErrors *UncBand_;
  UncBand_ = FULLGRAPH(hMC, hMC_JESPlus, hMC_JESMinus, hMC_JERPlus, hMC_JERMinus, hMC_METUCPlus, hMC_METUCMinus, hMC_bTagPlus, hMC_bTagMinus, h_QCD_iso);
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
    UncBand_Ratio = RATIOGRAPH(hMC, hMC_JESPlus, hMC_JESMinus, hMC_JERPlus, hMC_JERMinus, hMC_METUCPlus, hMC_METUCMinus, hMC_bTagPlus, hMC_bTagMinus, h_QCD_iso);
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

  TString outFile("plots/KineFit_ele_Sel_v1_test1_datadriven_from_data/base/");
  outFile += histname;
  outFile += "_ele.png";
  c1->SaveAs(outFile);
  c1->Close();

}

void example_stack_all(){

//   example_stack("pt_mu","P_{T}[GeV/c]",5,true, true,true, false, true,20.0,100.0); // log, drawdata, ratio, drawsignal, axisrange,xmin,xmax
//   example_stack("eta_mu","#eta^{#mu}",1,true,true,true,false,true,-2.5,3.5);
//   example_stack("Pre_RelIso","Muon Rel Iso",1,true,true, true,false,true,0,0.28);
//   example_stack("pt_jet","P_{T}[GeV/c]",6,true,true,true,false, true,20.0,100.0);
//   example_stack("eta_jet","#eta^{jets}",1,true,true,true,false,true,-3.0,3.0);
//   example_stack("final_multi_jet","N_{jets}",1,true,true,true,false, true,3.0,11.0);
//   example_stack("final_pt_met","E_{T}^{miss}[GeV/c]",10,true,true,true);
//   example_stack("btag_jet","CSVM",1,true,true, true,false,true,-1.2,1.2);
//   example_stack("Final_RelIso","Muon Rel Iso",1,true,true, true,false,true,0,0.15);
//   //  example_stack("svMass_jID","M_{sv}[GeV/c]",1,true,true,true,true, true, 0.25, 5.0);
//   example_stack("nvtx","N_{vertex}",1,true,true,true,false,true, 0, 45);
//   example_stack("CSVM_count","btag",1,true,true, true,false,true,1,8);
//   example_stack("Pre_RelIso","Muon Rel Iso",1,true,true, true,false,true,0,0.28);
//   example_stack("multi_jet","N_{jets}",1,true,true,true,false, true,0.0,9.0);
//   example_stack("pt_met","E_{T}^{miss}[GeV/c]",10,true,true,true);

  
   // final required plots:
//  example_stack("mjj_kfit_Id_probfit1","M_{jj}[GeV]",1,false,false,false,true); // blinded
  example_stack("mjj_kfit_Id_probfit1","M_{jj}[GeV]",1,false,true,true,true); // unblinded

  
  example_stack("cutflow","",1,true,true,true,true,true,1,5,true); // log, drawdata, ratio, drawsignal, axisrange, xmin, xmax, label // good_v1
//  example_stack("cutflow","",1,true,true,true,true,true,1,6,true);


  //  example_stack("cutflow","cutflow",1,true,true,true,false,0,10,true);
  // example_stack("diJet_Mass",1,false,false,false);
  // example_stack("mjj_kinfit_btag",1,false,false,false,true);
  //  example_stack("base/mjj_kfit",1,false,false,false,true);
  //  example_stack("base/mjj_kfit_sv",1,false,false,false,true);
  //  example_stack("base/svmass",1,false,false,false,true);
  // example_stack("pri_vtxs",1,true,true,true);
  // example_stack("nvtx_alljet",1,true,true,true);
  // example_stack("w_mt",4,true,true,true);
  // //  example_stack("angle",1,true,true, true);
  // // example_stack("btag_value_pre",1,true,true, true);
  // example_stack("btag_value_final",1,true,true, true);
  // example_stack("pre_btag_jet_mult",1,true,true, true);

  // example_stack("btag_jet_mult",1,true,true, true);
  // example_stack("w_mt_kinfit",1,true,true, true);
  // example_stack("w_mt",1,true,true, true);
  // example_stack("Muon_mult_final",1,true,true, true);
  example_stack("Pre_RelIso","Relative Isolation",1,true,true, true,false,true,0,0.29);
  //  example_stack("Final_RelIso","Relative Isolation",1,true,true, true,false,true,0,0.09);
  // example_stack("deltaEta_mu_bjet",1,true,true, true);
  // example_stack("deltaPhi_mu_bjet",1,true,true, true);
  // example_stack("deltaR_mu_bjet",1,true,true, true);
  // example_stack("dijetMass_mu_bjet",1,true,true, true);
  // example_stack("pujetid",1,true,true,true);
  
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
