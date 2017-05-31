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
void example_stack(TString histname, TString xaxis_title, int bin, bool log=false, bool drawdata=true, bool ratio=false, bool drawsignal=false, bool axisrange=false, double xmin=0, double xmax=10,bool label=false )
{
  ////////////////////////////////////////////////////////////
  // 
  // Basic decorations for stacked plots
  //
  ////////////////////////////////////////////////////////////
  
  gStyle->SetOptStat(0);
  c1 = new TCanvas();
  const float xpad[2] = {0.,1};
  //const float ypad[4] = {0.,0.3,0.3,1.0};
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
  
  TLegend* leg = new TLegend(0.7818792,0.6261504,0.9312081,0.9198861,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  //TPaveText *pt = new TPaveText(0.15,0.92,0.9,0.99, "brNDC");
  TPaveText *pt = new TPaveText(0.15,0.94,0.93,0.97, "brNDC"); // good_v1
  pt->SetBorderSize(1);
  pt->SetFillColor(19);
  pt->SetFillStyle(0);
  pt->SetLineColor(0);
  pt->SetTextFont(132);

  ////////////////////////////////////////////////////////////
  // 
  // Stack all the other MC histos on top of the base histo
  //
  ////////////////////////////////////////////////////////////

  TString inFile("13TeV/outputDir/");
  
  //VV is the base histo
  //TFile* f1 = TFile::Open(inFile+"WW_Merged.root");
  TFile* f1 = TFile::Open(inFile+"all_VV.root");
  if(f1 == 0) return;
  if(f1->IsZombie()){f1->Close(); return;}
  TH1F* h1_base = (f1->Get("base/"+histname))->Clone("h1_base");
  double scale_factor = 1.0;
  h1_base->Scale(scale_factor);  
  h1_base->SetFillColor(kYellow);
  if(axisrange){
    h1_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  leg->AddEntry(h1_base,"VV","F");
  //Define stacked histo
  THStack* MuptStack = new THStack("MuptStack","");
  MuptStack->Add(h1_base);
  
  //hMC = all Bkg MC samples
  TH1F* hMC = (TH1F*)h1_base->Clone("hMC");
  hMC->Reset();  
  
  //hSig = Bkg MC+ signal MC
  TH1F* hSig = (TH1F*)h1_base->Clone("hSig"); 
  hSig->Reset(); 
  hSig->SetLineColor(kCyan);
  hSig->Add(h1_base);
  hSig->SetLineColor(kMagenta+0);
  hSig->SetLineStyle(2);
  hSig->SetLineWidth(3);
  hSig->SetFillColor(0);
  
  //stack all the MC samples
  stackHisto(inFile, "all_QCD.root", "all_QCD", histname, kBlue, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  //stackHisto(inFile, "all_DY.root", "all_DY", histname, kRed, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  stackHisto(inFile, "DYJetsToLL_Merged.root", "all_DY", histname, kRed, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  stackHisto(inFile, "all_ST.root", "all_ST", histname, kGreen , 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  //stackHisto(inFile, "all_WJets.root", "all_WJets", histname, kMagenta+3, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  stackHisto(inFile, "WJetsToLNu_Merged.root", "all_WJets", histname, kMagenta+3, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  stackHisto(inFile, "all_TTJets.root","ttbar", histname, kCyan, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  //stackHisto(inFile, "all_H120.root","Hplus", histname,  6, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig, leg);
  
  
  ////////////////////////////////////////////////////////////
  //
  // Plot DATA points on the stacked MC histos
  //
  ////////////////////////////////////////////////////////////

  TFile* data_ = TFile::Open(inFile+"all_MuData.root");
  if(data_ == 0) return;
  if(data_->IsZombie()){data_->Close(); return;}
  //TH1F* data = (data_->Get("base/BTag/"+histname))->Clone("data");
  TH1F* data = (data_->Get("base/"+histname))->Clone("data");
  ///data->SetBinErrorOption(TH1::kPoisson); // added this line for Asymmetric Error Bars for Poisson Event Counts
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.9);
  data->Rebin(bin);
  if(axisrange){
    data->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  data->SetFillColor(kBlack);
  
  //lable x-axis, for cutflow
  if(label){
    TString steps[15] = {"","= 1 muon","#geq 4 jets","#slash{E}_{T} #geq 0GeV", "#geq 2 b-tagged jets","After KinFit","","","","","","","","",""};
    const size_t nsteps = sizeof(steps)/sizeof(TString);
    for(int istep=0; istep<nsteps; istep++ ){
      data->GetXaxis()->SetBinLabel(istep+1, steps[istep]);
      data->GetXaxis()->SetLabelSize(0.06);
      data->GetXaxis()->LabelsOption("u");
    }
  }
    
  
  if(drawdata)leg->AddEntry(data,"Data","P"); 
  if(drawsignal)leg->AddEntry(hSig,"Sig+bkgs","L"); // only for hSig histogram
  TPaveText *cct = new TPaveText(0.2013423,0.7754898,0.4010067,0.8962187,"brNDC");
  cct->SetFillColor(19);
  cct->SetFillStyle(0);
  cct->SetLineColor(0);
  cct->SetBorderSize(1);
  cct->AddText("M_{H^{+}} = 120 GeV");
  ///cct->AddText("Br(t#rightarrow H^{+}b) = 0.1");
  TPaveText *ch = new TPaveText(0.753,0.7754898,0.6510067,0.8962187,"brNDC");
  ch->SetFillColor(19);
  ch->SetFillStyle(0);
  ch->SetLineColor(0);
  ch->SetBorderSize(1);
  ch->AddText("#mu + jets");
  pt->SetTextSize(0.059);
  TText *text = pt->AddText("CMS Preliminary,    #sqrt{s} = 13 TeV,    35.5 fb^{-1}");
  //TText *text = pt->AddText("#sqrt{s} = 13 TeV");
  text->SetTextAlign(11);
  if(histname.Contains("Pre_RelIso") || histname.Contains("cutflow")) {
   data->SetAxisRange(1.0, 5.0e9 ,"Y");
  }else{
   data->SetAxisRange(1.0, 1.0e6 ,"Y");
  }

  data->GetYaxis()->SetTitleOffset(1.45);
  data->GetYaxis()->SetTitle("Events");
  data->GetXaxis()->SetTitle(xaxis_title);
  if(drawdata){
    data->Draw("E"); // not for dijet mass
    MuptStack->Draw("HISTSAME"); // no histsame for dijet calculation
    MuptStack->GetXaxis()->SetTitle(xaxis_title);
    MuptStack->SetMinimum(1.0);
    MuptStack->SetMaximum(1.1*(MuptStack->GetMaximum()));
    //stack->GetXaxis()->SetLimits(xmin, xmax);
    //stack->Draw("H"); // make sure it's redrawn
    MuptStack->Draw("HIST");
    MuptStack->GetXaxis()->SetTitle(xaxis_title);
    MuptStack->GetYaxis()->SetTitle("Events");
    MuptStack->GetYaxis()->SetTitleOffset(1.45);
    MuptStack->Draw("HISTSAME");
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
  
  ////////////////////////////////////////////////////////////
  //
  // DATA - Bkg
  //
  ////////////////////////////////////////////////////////////

  if(ratio){
    c1->cd(2);
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.3);
    gPad->SetGridy();
    gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
    //gPad->SetPad(xpad[0],0.05,xpad[1],ypad[2]);
    TH1F *hRatio = (TH1F*)data->Clone("hRatio");
    //hRatio->Reset();
    //hRatio->Add(data);
    //hRatio->Add(hMC, -1);
    //hRatio->Divide(hMC);
    int nbin= data->GetSize();
    cout<<"bin = "<<nbin<<endl;
    double div=0.0; 
    double sigma=0.0;
    double a1 =0.0;
    double binC_data= 0.0;
    double binC_mc = 0.0;

    for(int i=1; i<nbin; i++){
      binC_data = data->GetBinContent(i);
      binC_mc = hMC->GetBinContent(i);
      cout<<i<<endl;
      cout<<binC_mc<<endl;
      cout<<binC_data<<endl;
      cout<<endl;
      if(binC_mc !=0&& binC_data!=0){
        div = binC_data/binC_mc;
        a1 = sqrt(1.0/binC_data + 1.0/binC_mc);
        sigma = div*a1;
        hRatio->SetBinContent(i, div);
        hRatio->SetBinError(i, sigma);
      }
    }
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(1.0);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetRangeUser(-2.5, 2.5);
    //hRatio->GetXaxis()->SetRangeUser(xmin, xmax);
    hRatio->GetXaxis()->SetTitle(xaxis_title);
    hRatio->GetYaxis()->SetTitleOffset(0.5);
    hRatio->GetYaxis()->SetTitle("#frac{Data}{MC}");
    hRatio->GetYaxis()->CenterTitle();
    hRatio->GetYaxis()->SetTitleSize(0.1);
    hRatio->GetXaxis()->SetTitleSize(0.1);
    hRatio->GetXaxis()->SetLabelSize(0.17); // 0.1
    hRatio->GetXaxis()->LabelsOption("u"); // extra
    hRatio->GetYaxis()->SetLabelSize(0.06);
    hRatio->Draw("e1"); // use "P" or "AP"
    //hRatio->Draw("E same");
    c1->Update();
  }

  //TString outFile("plots/KineFit_ele_Sel_v1_test1_datadriven_from_data/base/");
  TString outFile("/home/ravindra/lxplus/StackHisto/Btag");
  outFile += histname;
  outFile += "_ele.png";
  c1->SaveAs(outFile);
  //c1->Close();

}
void stackHisto(TString path, TString filename, TString lable, TString histname, int color= 0, double scale=1, bool axisrange=false, double xmin=0, double xmax=10, THStack* MuptStack, TH1F* hMC, TH1F* hSig, TLegend* leg){
  TFile* f2 = TFile::Open(path+filename);
  if(f2 == 0) return;
  if(f2->IsZombie()){f2->Close(); return;}
  TH1F* h2_base = (f2->Get("base/"+histname))->Clone("h2_base");
  h2_base->Scale(scale);  
  h2_base->SetFillColor(color);
  if(axisrange){
    h2_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  leg->AddEntry(h2_base,lable,"F");
  MuptStack->Add(h2_base);
  hMC->Add(h2_base);
  hSig->Add(h2_base);
}

void example_stack_all(){

   // final required plots:
  //example_stack("pt_mu","",1,true,true,true,true,true,1,5, true); // log, drawdata, ratio, drawsignal, axisrange, xmin, xmax, label // good_v1
  //example_stack("pt_mu","", 1,true,true, true,false,true,0, 5, true);
  // example_stack("pt_mu","", 1,true,true, true,false,true,0, 5, false);
  //example_stack("pt_mu", "pt_mu", 1, true, true, true, true, false, 0, 5, false);
  //example_stack("pt_mu", "pt_mu", 1, true, true, true, false, false, 0, 5, false);
  //example_stack("eta_mu", "eta_mu", 1, true, true, true, false, false, 0, 5, false);
  example_stack("pt_jet", "pt_jet", 1, true, true, true, false, false, 0, 5, false);
  example_stack("eta_jet", "pt_jet", 1, true, true, true, false, false, 0, 5, false);
  example_stack("mjj", "mjj", 1, true, false, false, false, false, 0, 5, false);
  //example_stack("multi_jet", "multi_jet", 1, true, true, true, false, false, 0, 5, false);
  example_stack("nvtx", "nvtx", 1, true, true, true, false, false, 0, 5, false);
  //example_stack("nvtx_nocut", "nvtx_nocut", 1, true, true, true, false, false, 0, 5, false);
  //example_stack("nvtx_1mu", "nvtx_1mu", 1, true, true, true, false, false, 0, 5, false);
 // example_stack("/wmt","m_{T}[GeV]", 2, true, true, true, false, false,0,175);
  //example_stack("cutflow", "cutflow", 1, true, true, true, false, false, 0, 5, true);
  //example_stack("Pre_RelIso", "Pre_RelIso", 1, true, true, true, false, false, 0, 5, false);

  //example_stack("pt_mu","", 1,true,true, false,false,false,0, 5, false);
  //example_stack("pt_mu","", 1,true,false, false,false,false,0, 5, false);
  //example_stack("pt_mu","", 1,false,true, true,false,true,0, 5, true);
  
}
void example_stack_kfit(TString dir="BTag"){
  //example_stack(dir+"/pt_jet","jet P_{T}[GeV]", 1, true,true,true,false, false,20.0,100.0);
  //example_stack(dir+"/eta_jet","#eta^{jets}", 1, true,true,true,false,false,-3.0,3.0);
  //example_stack(dir+"/pt_mu","#mu P_{T}[GeV]", 1, true,true,true,false, false,20.0,100.0);
  //example_stack(dir+"/eta_mu","#eta^{#mu}", 1, true,true,true,false,false,-3.0,3.0);
  ////example_stack(dir+"/bDiscr_Loose","bDiscr", 1, true,true,true,false,false,0.0,1.0);
  //example_stack("pfCISV","pfCombinedInclusiveSecondaryVertexV2BJetTags", 1, true,true,true,false,false,0.0,1.0);
  //example_stack("pfCMVA","pfCombinedMVAV2BJetTags", 1, true,true,true,false,false,0.0,1.0);
  //example_stack("pfCCvsL","pfCombinedCvsLJetTags", 1, true,true,true,false,false,0.0,1.0);
  //example_stack("pfCCvs","pfCombinedCvsBJetTags", 1, true,true,true,false,false,0.0,1.0);
  //example_stack("mjj","m_{jj}",1,true,true,true);
  //////example_stack(dir+"/final_pt_mu","#pt_mu", 1, true,true,true,false, false,20.0,100.0);
  //example_stack(dir+"/final_pt_met","#pt_met", 1, true,true,true,false, false,20.0,100.0);
 // example_stack(dir+"/final_RelIso_mu", "final_RelIso", 1, true, true, true, false, false, 0, 5);

  if(dir == "KinFit"){
    example_stack(dir+"/kfJet1_pt","P_{T}[GeV]",1,true,true,true,true, true,20.0,200.0);
    example_stack(dir+"/kfJet2_pt","P_{T}[GeV]",1,true,true,true,true, true,20.0,200.0);
    example_stack(dir+"/kfJet1_eta","#eta^{jets}",1,true,true,true,true,true,-3.0,3.0);
    example_stack(dir+"/kfJet2_eta","#eta^{jets}",1,true,true,true,true,true,-3.0,3.0);
  }

  example_stack(dir+"/final_multi_jet","N_{jets}",1,true,true,true,false,true,3,15);
  //example_stack(dir+"/wmt","m_{T}[GeV]",2,true,true,true,false,false,0,175);
  /////example_stack(dir+"/CSVM_count","btagged jet multiplicity",1,true,true, true,false,true,1,8);
  example_stack(dir+"/nvtx","N_{vertex}",1,true,true,true);
  //example_stack(dir+"/rho","rho",1,true,true,true);
  //example_stack(dir+"/rho_0","rho_0",1,true,true,true);
  example_stack(dir+"/chi2","chi2",1,true,true,true);
  example_stack(dir+"/ndof","ndof",1,true,true,true);
} 
  
