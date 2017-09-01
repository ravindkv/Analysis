#include "interface/HistogramPlotter.hh" 
#include <iostream> 
#include <iomanip> 
#include <math.h>
 
ClassImp(HistogramPlotter) 

void HistogramPlotter::CreateAnalHistos(TString cutflowType, TFile* outFile_)
{

  //Define Histograms 
  //base histo
  InitHist(cutflowType, "", outFile_); 
  addHisto("cutflow", cutflowType, 15 , 0., 15.); 
  addHisto("totalEvents", cutflowType, 10 , 0., 10000000000.); 
  addHisto("intimepu", cutflowType, 6000, 0., 1000.);
  addHisto("outoftimepu", cutflowType, 6000, 0., 1000.);
  addHisto("totalpu", cutflowType, 6000, 0., 1000.);
  addHisto("trueintimepu", cutflowType, 6000, 0., 1000.);
  addHisto("trueoutoftimepu", cutflowType, 6000, 0., 1000.);
  addHisto("truetotalpu", cutflowType, 6000, 0., 1000.);
  addHisto("RelIso_mu",cutflowType, 40, 0, 1.0);
  addHisto("hepNUP", cutflowType, 100, 1., 20.);
  //Scale factors
  addHisto("SF_hepNUP_WJets",cutflowType, 1000, 0, 1000);
  addHisto("SF_hepNUP_DYJets",cutflowType, 1000, 0, 1000);
  addHisto("SF_sampleWeight",cutflowType, 1000, 0, 1000);
  addHisto("SF_weightPU",cutflowType, 1000, 0, 1000);
  addHisto("SF_topPtWeights",cutflowType, 1000, 0, 10);
  addHisto("SF_muonSF",cutflowType, 1000, 0, 10);
 
  ////////////////// Isolation ///////////////// 
  //base/Iso histo
  InitHist("Iso", cutflowType, outFile_);
  addHisto("cutflow", cutflowType+"/Iso", 15 , 0., 15.); 
  addHisto("pre_RelIso_mu",cutflowType+"/Iso", 40,0,0.5);
  addHisto("final_RelIso_mu",cutflowType+"/Iso", 40,0,0.5);
  addHisto("nvtx", cutflowType+"/Iso", 100, 0., 100.);
  addHisto("nvtx_6Kbins", cutflowType+"/Iso", 6000, 0., 1000.);
  addHisto("rhoAll", cutflowType+"/Iso", 100, 0., 100.);
  addHisto("chi2", cutflowType+"/Iso", 100, 0., 500.);
  addHisto("ndof", cutflowType+"/Iso", 100, 0., 500.);
  
  //LOOSE: base/Iso/Btag histo
  InitHist("Iso/BTag", cutflowType, outFile_);
  addHisto("pfCISV", cutflowType+"/Iso/BTag", 50, -2., 2.);
  addHisto("pfCMVA", cutflowType+"/Iso/BTag", 50, -2., 2.);
  addHisto("pfCCvsL",cutflowType+"/Iso/BTag", 50, -2., 2.);
  addHisto("pfCCvsB", cutflowType+"/Iso/BTag", 50, -2., 2.);
  addHisto("CSVL_count", cutflowType+"/Iso/BTag", 50,0,10);
  addHisto("pfCCvsL_0",cutflowType+"/Iso/BTag", 50, -2., 2.);
  addHisto("pfCCvsL_1",cutflowType+"/Iso/BTag", 50, -2., 2.);
  addHisto("pfCCvsB_0", cutflowType+"/Iso/BTag", 50, -2., 2.);
  addHisto("pfCCvsB_1", cutflowType+"/Iso/BTag", 50, -2., 2.);
  addHisto("mjj",cutflowType+"/Iso/BTag", 400, 0, 2000);
  addHisto("final_RelIso_mu",cutflowType+"/Iso/BTag", 40,0,0.5);
  addHisto("final_multi_jet", cutflowType+"/Iso/BTag", 10,0,10);
  addHisto("nvtx", cutflowType+"/Iso/BTag", 100, 0., 100.);
  addHisto("nvtx_6Kbins", cutflowType+"/Iso/BTag", 6000, 0., 1000.);
  addHisto("rhoAll", cutflowType+"/Iso/BTag", 100, 0., 100.);
  addHisto("chi2", cutflowType+"/Iso/BTag", 100, 0., 500.);
  addHisto("ndof", cutflowType+"/Iso/BTag", 100, 0., 500.);
  addHisto("wmt", cutflowType+"/Iso/BTag", 50, 0., 500.);
  addHisto("pt_bjet", cutflowType+"/Iso/BTag", 50, 0., 500.);
  addHisto("eta_bjet", cutflowType+"/Iso/BTag", 50, -5.0, 5.0);
  
  //base/Iso/KinFit histo
  InitHist("Iso/KinFit", cutflowType, outFile_);
  addHisto("final_RelIso_mu",cutflowType+"/Iso/KinFit", 40,0,0.5);
  addHisto("final_multi_jet", cutflowType+"/Iso/KinFit", 10,0,10);
  addHisto("CSVL_count", cutflowType+"/Iso/KinFit", 50,0,10);
  addHisto("CSVM_count", cutflowType+"/Iso/KinFit", 10,0,10);
  addHisto("wmt", cutflowType+"/Iso/KinFit", 50, 0., 500.);
  addHisto("rhoAll", cutflowType+"/Iso/KinFit", 100, 0., 100.);
  addHisto("chi2", cutflowType+"/Iso/KinFit", 100, 0., 500.);
  addHisto("ndof", cutflowType+"/Iso/KinFit", 100, 0., 500.);
  addHisto("nvtx", cutflowType+"/Iso/KinFit", 100, 0., 100.);
  addHisto("nvtx_6Kbins", cutflowType+"/Iso/KinFit", 6000, 0., 1000.);
  addHisto("kfJet1_pt", cutflowType+"/Iso/KinFit", 40, 0, 200);
  addHisto("kfJet2_pt", cutflowType+"/Iso/KinFit", 40, 0, 200);
  addHisto("kfJet1_eta", cutflowType+"/Iso/KinFit", 60, -3.0, 3.0);
  addHisto("kfJet2_eta", cutflowType+"/Iso/KinFit", 60, -3.0, 3.0);
  addHisto("kfJet1_phi", cutflowType+"/Iso/KinFit", 63, -M_PI, M_PI);
  addHisto("kfJet2_phi", cutflowType+"/Iso/KinFit", 63, -M_PI, M_PI);
  addHisto("mjj_kfit",cutflowType+"/Iso/KinFit", 400, 0, 2000);
  addHisto("mjj_kfit_CTag",cutflowType+"/Iso/KinFit", 400, 0, 2000);
  addHisto("mjj_kfit_noCTag",cutflowType+"/Iso/KinFit", 400, 0, 2000);
  addHisto("mjj_kfit_Id",cutflowType+"/Iso/KinFit", 400, 0, 2000);
  addHisto("mjj_kfit_Id_probfit1",cutflowType+"/Iso/KinFit", 400, 0, 2000);
  addHisto("mjj_kfit_Id_probfit2",cutflowType+"/Iso/KinFit", 400, 0, 2000);
  
  ////////////////// NonIsolation ///////////////// 
  //base/NonIso histo
  InitHist("NonIso", cutflowType, outFile_);
  addHisto("cutflow", cutflowType+"/NonIso", 15 , 0., 15.); 
  addHisto("final_RelIso_mu",cutflowType+"/NonIso", 40,0,0.5);
  addHisto("nvtx", cutflowType+"/NonIso", 100, 0., 100.);
  addHisto("nvtx_6Kbins", cutflowType+"/NonIso", 6000, 0., 1000.);
  addHisto("rhoAll", cutflowType+"/NonIso", 100, 0., 100.);
  addHisto("chi2", cutflowType+"/NonIso", 100, 0., 500.);
  addHisto("ndof", cutflowType+"/NonIso", 100, 0., 500.);
  addHisto("pfCISV", cutflowType+"/NonIso", 50, -2., 2.);
  addHisto("pfCMVA", cutflowType+"/NonIso", 50, -2., 2.);
  addHisto("pfCCvsL",cutflowType+"/NonIso", 50, -2., 2.);
  addHisto("pfCCvsB", cutflowType+"/NonIso", 50, -2., 2.);
  addHisto("pfCCvsL_0",cutflowType+"/NonIso", 50, -2., 2.);
  addHisto("pfCCvsL_1",cutflowType+"/NonIso", 50, -2., 2.);
  addHisto("pfCCvsB_0", cutflowType+"/NonIso", 50, -2., 2.);
  addHisto("pfCCvsB_1", cutflowType+"/NonIso", 50, -2., 2.);
  addHisto("CSVL_count", cutflowType+"/NonIso", 50,0,10);
  addHisto("pt_metJESJER", cutflowType+"/NonIso", 50, 0., 500.);
  addHisto("mjj",cutflowType+"/NonIso", 400, 0, 2000);
  
  //base/NonIso/Btag histo
  InitHist("NonIso/BTag", cutflowType, outFile_);
  addHisto("final_RelIso_mu",cutflowType+"/NonIso/BTag", 40,0,0.5);
  addHisto("final_multi_jet", cutflowType+"/NonIso/BTag", 10,0,10);
  addHisto("nvtx", cutflowType+"/NonIso/BTag", 100, 0., 100.);
  addHisto("nvtx_6Kbins", cutflowType+"/NonIso/BTag", 6000, 0., 1000.);
  addHisto("rhoAll", cutflowType+"/NonIso/BTag", 100, 0., 100.);
  addHisto("chi2", cutflowType+"/NonIso/BTag", 100, 0., 500.);
  addHisto("ndof", cutflowType+"/NonIso/BTag", 100, 0., 500.);
  addHisto("wmt", cutflowType+"/NonIso/BTag", 50, 0., 500.);
  addHisto("mjj",cutflowType+"/NonIso/BTag", 400, 0, 2000);
  addHisto("pt_bjet", cutflowType+"/NonIso/BTag", 50, 0., 500.);
  addHisto("eta_bjet", cutflowType+"/NonIso/BTag", 50, -5.0, 5.0);
  
  //base/NonIso/KinFit histo
  InitHist("NonIso/KinFit", cutflowType, outFile_);
  addHisto("final_RelIso_mu",cutflowType+"/NonIso/KinFit", 40,0,0.5);
  addHisto("final_multi_jet", cutflowType+"/NonIso/KinFit", 10,0,10);
  addHisto("CSVL_count", cutflowType+"/NonIso/KinFit", 50,0,10);
  addHisto("wmt", cutflowType+"/NonIso/KinFit", 50, 0., 500.);
  addHisto("rhoAll", cutflowType+"/NonIso/KinFit", 100, 0., 100.);
  addHisto("chi2", cutflowType+"/NonIso/KinFit", 100, 0., 500.);
  addHisto("ndof", cutflowType+"/NonIso/KinFit", 100, 0., 500.);
  addHisto("nvtx", cutflowType+"/NonIso/KinFit", 100, 0., 100.);
  addHisto("nvtx_6Kbins", cutflowType+"/NonIso/KinFit", 6000, 0., 1000.);
  addHisto("kfJet1_pt", cutflowType+"/NonIso/KinFit", 40, 0, 200);
  addHisto("kfJet2_pt", cutflowType+"/NonIso/KinFit", 40, 0, 200);
  addHisto("kfJet1_eta", cutflowType+"/NonIso/KinFit", 60, -3.0, 3.0);
  addHisto("kfJet2_eta", cutflowType+"/NonIso/KinFit", 60, -3.0, 3.0);
  addHisto("kfJet1_phi", cutflowType+"/NonIso/KinFit", 63, -M_PI, M_PI);
  addHisto("kfJet2_phi", cutflowType+"/NonIso/KinFit", 63, -M_PI, M_PI);
  addHisto("mjj_kfit",cutflowType+"/NonIso/KinFit", 400, 0, 2000);
  addHisto("mjj_kfit_Id",cutflowType+"/NonIso/KinFit", 400, 0, 2000);
  addHisto("mjj_kfit_Id_probfit1",cutflowType+"/NonIso/KinFit", 400, 0, 2000);
  addHisto("mjj_kfit_Id_probfit2",cutflowType+"/NonIso/KinFit", 400, 0, 2000);
}


  //void HistogramPlotter::InitHist(TString dirname, TString parentDir, TFile *file)
  void HistogramPlotter::InitHist(TString dirname, TString parentDir, TFile *file)
{
  std::string name(dirname);
  std::string fullname;
  if(parentDir.Length() != 0){
    TDirectory *d = file->GetDirectory(parentDir.Data());
    d->mkdir(name.c_str());
    fullname = std::string(parentDir+"/"+dirname);
  }
  else{file->mkdir(name.c_str()); fullname = name;}
  TDirectory *d = file->GetDirectory(fullname.c_str());
  file->cd(fullname.c_str());
  addHisto("pt_jet", fullname, 50, 0., 500.);
  TH1 *h1 = getHisto("pt_jet", fullname);
  h1->SetDirectory(d);
  /*
  addHisto("test_hist","fullname",10, 0., 10.);
  h1=getHisto("test_hist", fullname);
  h1->SetDirectory(d);
  */
  addHisto("eta_jet", fullname, 50, -5.0, 5.0);
  h1 = getHisto("eta_jet", fullname);
  h1->SetDirectory(d);
  addHisto("phi_jet", fullname, 63, -M_PI, M_PI);
  h1 = getHisto("phi_jet", fullname);
  h1->SetDirectory(d);
  addHisto("multi_jet", fullname, 100, 0., 20.);
  h1 = getHisto("multi_jet", fullname);
  h1->SetDirectory(d);
/*
  addHisto("pt_ele", fullname, 50, 0., 500.);
  h1 = getHisto("pt_ele", fullname);
  h1->SetDirectory(d);
  addHisto("eta_ele", fullname, 50, -5.0, 5.0);
  h1 = getHisto("eta_ele", fullname);
  h1->SetDirectory(d);
  addHisto("phi_ele", fullname, 63, -M_PI, M_PI);
  h1 = getHisto("phi_ele", fullname);
  h1->SetDirectory(d);
*/
  addHisto("pt_mu", fullname, 50, 0., 500.);
  h1 = getHisto("pt_mu", fullname);
  h1->SetDirectory(d);
  addHisto("eta_mu", fullname, 50, -5.0, 5.0);
  h1 = getHisto("eta_mu", fullname);
  h1->SetDirectory(d);
  addHisto("phi_mu", fullname, 63, -M_PI, M_PI);
  h1 = getHisto("phi_mu", fullname);
  h1->SetDirectory(d);

  addHisto("pt_met", fullname, 50, 0., 500.);
  h1 = getHisto("pt_met", fullname);
  h1->SetDirectory(d);
  addHisto("final_pt_met", fullname, 50, 0., 500.);
  h1 = getHisto("final_pt_met", fullname);
  h1->SetDirectory(d);

  addHisto("phi_met", fullname, 63, -M_PI, M_PI);
  h1 = getHisto("phi_met", fullname);
  h1->SetDirectory(d);
  addHisto("final_phi_met", fullname, 63, -M_PI, M_PI);
  h1 = getHisto("final_phi_met", fullname);
  h1->SetDirectory(d);
}  

void HistogramPlotter::addHisto(TString name, TString dirname, int range, double min, double max)
{
  //TString fullname = name+"_"+dirname; 
  TString fullname = dirname+"/"+name;
  std::string hname(fullname);
  histos1_[fullname] = new TH1F(name.Data(), hname.c_str(), range, min, max); 
  //histos1_[fullname] = new TH1D(name.Data(), hname.c_str(), range, min, max); 
  //histos1_[fullname]->Sumw2();
}

void HistogramPlotter::add2DHisto(TString name, TString dirname, int range1, double min1, double max1, int range2, double min2, double max2)
{
  //TString fullname = name+"_"+dirname;
  TString fullname = dirname+"/"+name;
  std::string hname(fullname); 
  histos2_[fullname] = new TH2D(name.Data(), hname.c_str(), range1, min1, max1, range2, min2, max2);
  //histos2_[fullname]->Sumw2();
}

void HistogramPlotter::fillHistoPU(TString name, TString dirname, int bin, double binCont)
{
  TH1* h = getHisto(name, dirname);
  if(h != 0) h->SetBinContent(bin, binCont);
}

void HistogramPlotter::fillHisto(TString name, TString dirname, double value, double weight)
{
  //TString fullname = name+"_"+dirname;
  //TString fullname = dirname+"/"+name;
  TH1* h = getHisto(name, dirname);
  if(h != 0) h->Fill(value, weight);
}

TH1* HistogramPlotter::getHisto(TString name, TString dirname)
{
  //TString fullname = name+"_"+dirname;
  TString fullname = dirname+"/"+name;
  TH1 * h = 0;
  if(histos1_.find(fullname) != histos1_.end())h = histos1_[fullname];
  else if(histos2_.find(fullname) != histos2_.end())h = histos2_[fullname];
  return h;
}


