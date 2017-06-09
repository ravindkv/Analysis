#include "interface/HistogramPlotter.hh" 
#include <iostream> 
#include <iomanip> 
#include <math.h>
 
ClassImp(HistogramPlotter) 

void HistogramPlotter::CreateAnalHistos(TString cutflowType, TFile* outFile_)
{

  //Define Histograms 
  InitHist(cutflowType, "", outFile_); 
  addHisto("cutflow", cutflowType, 50 , 0., 10.); 
  addHisto("totalEvents", cutflowType, 50 , 0., 100000000000.); 
  addHisto("pre_RelIso_mu",cutflowType, 40,0,0.5);
  addHisto("final_RelIso_mu",cutflowType, 40,0,0.5);
  addHisto("final_multi_mu",cutflowType, 10,0,10);
  addHisto("final_multi_jet", cutflowType, 100,0,10);
  addHisto("pfCISV", cutflowType, 200, -5., 5.);
  addHisto("pfCMVA", cutflowType, 200, -5., 5.);
  addHisto("pfCCvsL",cutflowType, 200, -5., 5.);
  addHisto("pfCCvsB", cutflowType, 200, -5., 5.);
  
  addHisto("pfCCvsL_0",cutflowType, 200, -5., 5.);
  addHisto("pfCCvsL_1",cutflowType, 200, -5., 5.);
  addHisto("pfCCvsB_0", cutflowType, 200, -5., 5.);
  addHisto("pfCCvsB_1", cutflowType, 200, -5., 5.);

  addHisto("CSVL_count", cutflowType, 50,0,10);
  addHisto("CSVM_count", cutflowType, 50,0,10);
  addHisto("wmt", cutflowType, 50, 0., 200.);
  addHisto("nvtx", cutflowType, 100, 0., 100.);
  addHisto("nvtx_6Kbins", cutflowType, 6000, 0., 1000.);
  addHisto("intimepu", cutflowType, 6000, 0., 1000.);
  addHisto("outoftimepu", cutflowType, 6000, 0., 1000.);
  addHisto("totalpu", cutflowType, 6000, 0., 1000.);
  addHisto("trueintimepu", cutflowType, 6000, 0., 1000.);
  addHisto("trueoutoftimepu", cutflowType, 6000, 0., 1000.);
  addHisto("truetotalpu", cutflowType, 6000, 0., 1000.);
  addHisto("rhoAll", cutflowType, 100, 0., 100.);
  addHisto("rhoAll0", cutflowType, 100, 0., 100.);
  addHisto("chi2", cutflowType, 100, 0., 500.);
  addHisto("ndof", cutflowType, 100, 0., 500.);
  addHisto("mjj", cutflowType, 50, 20., 250.);
  addHisto("hepNUP", cutflowType, 100, 1., 20.);
  
/*  
  addHisto("mjj_kfit", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_deltaR", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_deltaR_le70", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_deltaR_ge70", cutflowType, 40, 0.,200.);
  
  addHisto("mjj_kfit_Id", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt50", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt75", cutflowType, 40, 0.,200.);

  addHisto("mjj_kfit_drCut1", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_drCut2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_drCut3", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_drCut4", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_drCut5", cutflowType, 40, 0.,200.);

  addHisto("svMass_jID",cutflowType, 20,0.0,5.0);

  addHisto("mjj_kfit_Id_svCat1", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut1_svCat1", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut2_svCat1", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut3_svCat1", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut3_svCat1_muplus", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut3_svCat1_muminus", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut4_svCat1", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut5_svCat1", cutflowType, 40, 0.,200.);

  addHisto("mjj_kfit_Id_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut1_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut2_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut3_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut3_svCat2_muplus", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut3_svCat2_muminus", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut4_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut5_svCat2", cutflowType, 40, 0.,200.);

  addHisto("mjj_kfit_Id_pt50_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt50_drCut1_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt50_drCut2_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt50_drCut3_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt50_drCut3_svCat2_muplus", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt50_drCut3_svCat2_muminus", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt50_drCut4_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt50_drCut5_svCat2", cutflowType, 40, 0.,200.);

  addHisto("mjj_kfit_Id_pt75_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt75_drCut1_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt75_drCut2_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt75_drCut3_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt75_drCut3_svCat2_muplus", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt75_drCut3_svCat2_muminus", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt75_drCut4_svCat2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt75_drCut5_svCat2", cutflowType, 40, 0.,200.);

  addHisto("mjj_kfit_Id_drCut1", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut3", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut3_muplus", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut3_muminus", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut4", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_drCut5", cutflowType, 40, 0.,200.);

  addHisto("mjj_kfit_Id_pt50_drCut1", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt50_drCut2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt50_drCut3", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt50_drCut4", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt50_drCut5", cutflowType, 40, 0.,200.);

  addHisto("mjj_kfit_Id_pt75_drCut1", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt75_drCut2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt75_drCut3", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt75_drCut4", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_Id_pt75_drCut5", cutflowType, 40, 0.,200.);

  // match jet
  addHisto("mjj_kfit_match", cutflowType, 40, 0.,200.);
  addHisto("match_dijet_gamma", cutflowType, 50, 1, 10);
  addHisto("match_dijet_gamma_deltaR", cutflowType, 50, 1, 10);
  addHisto("deltaEta_match_jet12", cutflowType, 50,0,5.0);
  addHisto("deltaPhi_match_jet12", cutflowType, 50,0,7.0);
  addHisto("deltaR_match_jet12", cutflowType, 50, 0.0, 10.);
  // mis match
  addHisto("mjj_kfit_mismatch", cutflowType, 40, 0.,200.);
  addHisto("mismatch_dijet_gamma", cutflowType, 50, 1, 10);
  addHisto("mismatch_dijet_gamma_deltaR", cutflowType, 50, 1, 10);
  addHisto("deltaEta_mismatch_jet12", cutflowType, 50,0,5.0);
  addHisto("deltaPhi_mismatch_jet12", cutflowType, 50,0,7.0);
  addHisto("deltaR_mismatch_jet12", cutflowType, 50, 0.0, 10.);
  addHisto("gen_1st_quark_pt", cutflowType, 40, 0, 200);
  addHisto("gen_2nd_quark_pt", cutflowType, 40, 0, 200);
  addHisto("gen_1st_quark_eta", cutflowType, 60, -3.0, 3.0);
  addHisto("gen_2nd_quark_eta", cutflowType, 60, -3.0, 3.0);
  
  addHisto("mjj_kfit_chi2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_chi2_p2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_chi2_p4", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_chi2_p6", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_chi2_p8", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_chi2_1p0", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_chi2_1p2", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_chi2_1p4", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_chi2_1p6", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_chi2_1p8", cutflowType, 40, 0.,200.);
  addHisto("mjj_kfit_chi2_2p0", cutflowType, 40, 0.,200.);

  addHisto("mjj_kfit_sv", cutflowType, 40, 0.,200.);
  addHisto("svmass",cutflowType, 20,0.0,5.0);
  addHisto("chi2_fit", cutflowType, 100, 0.0, 9.0);
  addHisto("prob_fit",cutflowType, 40, -0.5, 1.5);
  addHisto("prob_fit_p2",cutflowType, 40, -0.5, 1.5);
  addHisto("prob_fit_p4",cutflowType, 40, -0.5, 1.5);
  addHisto("prob_fit_p6",cutflowType, 40, -0.5, 1.5);
  addHisto("prob_fit_p8",cutflowType, 40, -0.5, 1.5);
  addHisto("prob_fit_1p0",cutflowType, 40, -0.5, 1.5);
  addHisto("prob_fit_1p2",cutflowType, 40, -0.5, 1.5);
  addHisto("prob_fit_1p4",cutflowType, 40, -0.5, 1.5);
  addHisto("prob_fit_1p6",cutflowType, 40, -0.5, 1.5);
  addHisto("prob_fit_1p8",cutflowType, 40, -0.5, 1.5);
  addHisto("prob_fit_2p0",cutflowType, 40, -0.5, 1.5);
  add2DHisto("dijet_prob",cutflowType, 40, -0.5, 1.5, 40, 0.0, 200.);
  add2DHisto("dijet_chi2",cutflowType, 90, 0.0, 9.0, 40, 0.0, 200.);

  addHisto("AvTopPtWeight",cutflowType, 2, 0., 2.); //Average weight to be divided for each histogram
  addHisto("SVEffUncert", cutflowType, 3, 0., 3.); //To store SV Eff uncert.
 */
  InitHist("BTag", cutflowType, outFile_); 
  addHisto("pre_RelIso_mu",cutflowType+"/BTag", 40,0,0.5);
  addHisto("final_RelIso_mu",cutflowType+"/BTag", 40,0,0.5);
  addHisto("final_multi_mu",cutflowType+"/BTag", 10,0,10);
  addHisto("final_multi_jet", cutflowType+"/BTag", 100,0,10);
  addHisto("CSVL_count", cutflowType+"/BTag", 10,0,10);
  addHisto("CSVM_count", cutflowType+"/BTag", 10,0,10);
  addHisto("wmt", cutflowType+"/BTag", 50, 0., 200.);
  addHisto("rhoAll", cutflowType+"/BTag", 100, 0., 100.);
  addHisto("rhoAll0", cutflowType+"/BTag", 100, 0., 100.);
  addHisto("chi2", cutflowType+"/BTag", 100, 0., 500.);
  addHisto("ndof", cutflowType+"/BTag", 100, 0., 500.);
  addHisto("nvtx", cutflowType+"/BTag", 100, 0., 100.);
  addHisto("nvtx_6Kbins", cutflowType+"/BTag", 6000, 0., 1000.);
  addHisto("hepNUP", cutflowType+"/BTag", 100, 1., 20.);
/*
  InitHist("KinFit", cutflowType, outFile_);
  addHisto("pre_RelIso_mu",cutflowType+"/KinFit", 40,0,0.5);
  addHisto("final_RelIso_mu",cutflowType+"/KinFit", 40,0,0.5);
  addHisto("final_multi_mu",cutflowType+"/KinFit", 10,0,10);
  addHisto("final_multi_jet", cutflowType+"/KinFit", 10,0,10);
  addHisto("CSVL_count", cutflowType+"/KinFit", 10,0,10);
  addHisto("CSVM_count", cutflowType+"/KinFit", 10,0,10);
  addHisto("wmt", cutflowType+"/KinFit", 100, 0., 200.);
  addHisto("nvtx", cutflowType+"/KinFit", 50, 0., 50.);
  addHisto("kfJet1_pt", cutflowType+"/KinFit", 40, 0, 200);
  addHisto("kfJet2_pt", cutflowType+"/KinFit", 40, 0, 200);
  addHisto("kfJet1_eta", cutflowType+"/KinFit", 60, -3.0, 3.0);
  addHisto("kfJet2_eta", cutflowType+"/KinFit", 60, -3.0, 3.0);
  addHisto("kfJet1_phi", cutflowType+"/KinFit", 63, -M_PI, M_PI);
  addHisto("kfJet2_phi", cutflowType+"/KinFit", 63, -M_PI, M_PI);
  addHisto("recoJet1_pt", cutflowType+"/KinFit", 40, 0, 200);
  addHisto("recoJet2_pt", cutflowType+"/KinFit", 40, 0, 200);
  addHisto("recoJet1_eta", cutflowType+"/KinFit", 60, -3.0, 3.0);
  addHisto("recoJet2_eta", cutflowType+"/KinFit", 60, -3.0, 3.0);
  addHisto("recoJet1_phi", cutflowType+"/KinFit", 63, -M_PI, M_PI);
  addHisto("recoJet2_phi", cutflowType+"/KinFit", 63, -M_PI, M_PI);
*/
  
  //InitHist("OneLep2J", cutflowType, outFile_); 
  //InitHist("OneLep2J1B", cutflowType, outFile_); 

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
  addHisto("bDiscr_Loose", fullname, 200, -5., 5.);
  h1 = getHisto("bDiscr_Loose", fullname);
  h1->SetDirectory(d);

  addHisto("pt_ele", fullname, 50, 0., 500.);
  h1 = getHisto("pt_ele", fullname);
  h1->SetDirectory(d);
  addHisto("eta_ele", fullname, 50, -5.0, 5.0);
  h1 = getHisto("eta_ele", fullname);
  h1->SetDirectory(d);
  addHisto("phi_ele", fullname, 63, -M_PI, M_PI);
  h1 = getHisto("phi_ele", fullname);
  h1->SetDirectory(d);

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
  histos1_[fullname] = new TH1D(name.Data(), hname.c_str(), range, min, max); 
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

