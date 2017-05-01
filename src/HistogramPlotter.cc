#include "interface/HistogramPlotter.hh" 
#include <iostream> 
#include <iomanip> 
#include <math.h>
 
ClassImp(HistogramPlotter) 

void HistogramPlotter::CreateAnalHistos(TString cutflowType, TFile* outFile_)
{

  //Define Histograms 
  InitHist(cutflowType, "", outFile_, false); 
  addHisto("cutflow", cutflowType, 50 , 0., 10.); 
  addHisto("Pre_RelIso",cutflowType, 40,0,0.5);
  addHisto("Final_RelIso",cutflowType, 40,0,0.5);
  addHisto("Muon_mult_final",cutflowType, 10,0,10);
  addHisto("final_multi_jet", cutflowType, 100,0,10);
  addHisto("CSVL_count", cutflowType, 50,0,10);
  addHisto("CSVM_count", cutflowType, 50,0,10);
  addHisto("wmt", cutflowType, 50, 0., 200.);
  addHisto("nvtx", cutflowType, 50, 0., 50.);
  addHisto("nvtx_nocut", cutflowType, 50, 0., 50.);
  addHisto("nvtx_1mu", cutflowType, 50, 0., 50.);
  addHisto("nvtx_1mu_4jet", cutflowType, 50, 0., 50.);
  addHisto("nvtx_1mu_4jet_btag", cutflowType, 50, 0., 50.);
  addHisto("nvtx_1mu_4jet_btag_kinfit", cutflowType, 50, 0., 50.);
  addHisto("diJet_Mass", cutflowType, 50, 20., 250.);
  
  //NVTX FROM 0 To 10
  InitHist("nvtx_10", cutflowType, outFile_, true); 
  //befor 1 mu cut
  addHisto("multi_mu_10_b1mu", cutflowType+"/nvtx_10", 100, 0., 20.);
  addHisto("pt_mu_10_b1mu", cutflowType+"/nvtx_10", 50, 0., 500.);
  addHisto("eta_mu_10_b1mu", cutflowType+"/nvtx_10", 50, -5.0, 5.0);
  addHisto("pt_met_10_b1mu", cutflowType+"/nvtx_10", 50, 0., 500.);
  //after 1 mu cut
  addHisto("pt_mu_10_a1mu", cutflowType+"/nvtx_10", 50, 0., 500.);
  addHisto("eta_mu_10_a1mu", cutflowType+"/nvtx_10", 50, -5.0, 5.0);
  addHisto("pt_met_10_a1mu", cutflowType+"/nvtx_10", 50, 0., 500.);
  //before 4 jet cut
  addHisto("multi_jet_10_b4jet", cutflowType+"/nvtx_10", 100, 0., 20.);
  addHisto("pt_jet_10_b4jet", cutflowType+"/nvtx_10", 50, 0., 500.);
  addHisto("eta_jet_10_b4jet", cutflowType+"/nvtx_10", 50, -5.0, 5.0);
  addHisto("pt_met_10_b4jet", cutflowType+"/nvtx_10", 50, 0., 500.);
  addHisto("bdiscr_10_b4jet", cutflowType+"/nvtx_10", 50, 0, 2.0);
  //after 4 jet cut
  addHisto("multi_jet_10_a4jet", cutflowType+"/nvtx_10", 100, 0., 20.);
  addHisto("pt_jet_10_a4jet", cutflowType+"/nvtx_10", 50, 0., 500.);
  addHisto("eta_jet_10_a4jet", cutflowType+"/nvtx_10", 50, -5.0, 5.0);
  addHisto("pt_met_10_a4jet", cutflowType+"/nvtx_10", 50, 0., 500.);
  addHisto("bdiscr_10_a4jet", cutflowType+"/nvtx_10", 50, 0, 2.0);
  addHisto("multi_bjet_10_a4jet", cutflowType+"/nvtx_10", 100, 0., 20.);
  //after btag
  addHisto("pt_bjet_10_a4jet", cutflowType+"/nvtx_10", 50, 0., 500.);
  addHisto("eta_bjet_10_a4jet", cutflowType+"/nvtx_10", 50, -5.0, 5.0);
  addHisto("pt_met_10_a4jet_btag", cutflowType+"/nvtx_10", 50, 0., 500.);
  addHisto("bdiscr_10_a4jet_btag", cutflowType+"/nvtx_10", 50, 0, 2.0);
  
  //NVTX FROM 10 To 20
  InitHist("nvtx_20", cutflowType, outFile_, true); 
  //befor 1 mu cut
  addHisto("multi_mu_20_b1mu", cutflowType+"/nvtx_20", 100, 0., 20.);
  addHisto("pt_mu_20_b1mu", cutflowType+"/nvtx_20", 50, 0., 500.);
  addHisto("eta_mu_20_b1mu", cutflowType+"/nvtx_20", 50, -5.0, 5.0);
  addHisto("pt_met_20_b1mu", cutflowType+"/nvtx_20", 50, 0., 500.);
  //after 1 mu cut
  addHisto("pt_mu_20_a1mu", cutflowType+"/nvtx_20", 50, 0., 500.);
  addHisto("eta_mu_20_a1mu", cutflowType+"/nvtx_20", 50, -5.0, 5.0);
  addHisto("pt_met_20_a1mu", cutflowType+"/nvtx_20", 50, 0., 500.);
  //before 4 jet cut
  addHisto("multi_jet_20_b4jet", cutflowType+"/nvtx_20", 100, 0., 20.);
  addHisto("pt_jet_20_b4jet", cutflowType+"/nvtx_20", 50, 0., 500.);
  addHisto("eta_jet_20_b4jet", cutflowType+"/nvtx_20", 50, -5.0, 5.0);
  addHisto("pt_met_20_b4jet", cutflowType+"/nvtx_20", 50, 0., 500.);
  addHisto("bdiscr_20_b4jet", cutflowType+"/nvtx_20", 50, 0, 2.0);
  //after 4 jet cut
  addHisto("multi_jet_20_a4jet", cutflowType+"/nvtx_20", 100, 0., 20.);
  addHisto("pt_jet_20_a4jet", cutflowType+"/nvtx_20", 50, 0., 500.);
  addHisto("eta_jet_20_a4jet", cutflowType+"/nvtx_20", 50, -5.0, 5.0);
  addHisto("pt_met_20_a4jet", cutflowType+"/nvtx_20", 50, 0., 500.);
  addHisto("bdiscr_20_a4jet", cutflowType+"/nvtx_20", 50, 0, 2.0);
  addHisto("multi_bjet_20_a4jet", cutflowType+"/nvtx_20", 100, 0., 20.);
  //after btag
  addHisto("pt_bjet_20_a4jet", cutflowType+"/nvtx_20", 50, 0., 500.);
  addHisto("eta_bjet_20_a4jet", cutflowType+"/nvtx_20", 50, -5.0, 5.0);
  addHisto("pt_met_20_a4jet_btag", cutflowType+"/nvtx_20", 50, 0., 500.);
  addHisto("bdiscr_20_a4jet_btag", cutflowType+"/nvtx_20", 50, 0, 2.0);
  
  //NVTX FROM 20 To 30
  InitHist("nvtx_30", cutflowType, outFile_, true); 
  //befor 1 mu cut
  addHisto("multi_mu_30_b1mu", cutflowType+"/nvtx_30", 100, 0., 20.);
  addHisto("pt_mu_30_b1mu", cutflowType+"/nvtx_30", 50, 0., 500.);
  addHisto("pt_met_30_b1mu", cutflowType+"/nvtx_30", 50, 0., 500.);
  addHisto("eta_mu_30_b1mu", cutflowType+"/nvtx_30", 50, -5.0, 5.0);
  //after 1 mu cut
  addHisto("pt_mu_30_a1mu", cutflowType+"/nvtx_30", 50, 0., 500.);
  addHisto("pt_met_30_a1mu", cutflowType+"/nvtx_30", 50, 0., 500.);
  addHisto("eta_mu_30_a1mu", cutflowType+"/nvtx_30", 50, -5.0, 5.0);
  //before 4 jet cut
  addHisto("multi_jet_30_b4jet", cutflowType+"/nvtx_30", 100, 0., 20.);
  addHisto("pt_jet_30_b4jet", cutflowType+"/nvtx_30", 50, 0., 500.);
  addHisto("eta_jet_30_b4jet", cutflowType+"/nvtx_30", 50, -5.0, 5.0);
  addHisto("pt_met_30_b4jet", cutflowType+"/nvtx_30", 50, 0., 500.);
  addHisto("bdiscr_30_b4jet", cutflowType+"/nvtx_30", 50, 0, 2.0);
  //after 4 jet cut
  addHisto("multi_jet_30_a4jet", cutflowType+"/nvtx_30", 100, 0., 20.);
  addHisto("pt_jet_30_a4jet", cutflowType+"/nvtx_30", 50, 0., 500.);
  addHisto("pt_met_30_a4jet", cutflowType+"/nvtx_30", 50, 0., 500.);
  addHisto("eta_jet_30_a4jet", cutflowType+"/nvtx_30", 50, -5.0, 5.0);
  addHisto("bdiscr_30_a4jet", cutflowType+"/nvtx_30", 50, 0, 2.0);
  addHisto("multi_bjet_30_a4jet", cutflowType+"/nvtx_30", 100, 0., 20.);
  //after btag
  addHisto("pt_bjet_30_a4jet", cutflowType+"/nvtx_30", 50, 0., 500.);
  addHisto("eta_bjet_30_a4jet", cutflowType+"/nvtx_30", 50, -5.0, 5.0);
  addHisto("pt_met_30_a4jet_btag", cutflowType+"/nvtx_30", 50, 0., 500.);
  addHisto("bdiscr_30_a4jet_btag", cutflowType+"/nvtx_30", 50, 0, 2.0);
  
  //NVTX FROM 30 To 40
  InitHist("nvtx_40", cutflowType, outFile_, true); 
  //befor 1 mu cut
  addHisto("multi_mu_40_b1mu", cutflowType+"/nvtx_40", 100, 0., 20.);
  addHisto("pt_mu_40_b1mu", cutflowType+"/nvtx_40", 50, 0., 500.);
  addHisto("pt_met_40_b1mu", cutflowType+"/nvtx_40", 50, 0., 500.);
  addHisto("eta_mu_40_b1mu", cutflowType+"/nvtx_40", 50, -5.0, 5.0);
  //after 1 mu cut
  addHisto("pt_mu_40_a1mu", cutflowType+"/nvtx_40", 50, 0., 500.);
  addHisto("pt_met_40_a1mu", cutflowType+"/nvtx_40", 50, 0., 500.);
  addHisto("eta_mu_40_a1mu", cutflowType+"/nvtx_40", 50, -5.0, 5.0);
  //before 4 jet cut
  addHisto("multi_jet_40_b4jet", cutflowType+"/nvtx_40", 100, 0., 20.);
  addHisto("pt_jet_40_b4jet", cutflowType+"/nvtx_40", 50, 0., 500.);
  addHisto("pt_met_40_b4jet", cutflowType+"/nvtx_40", 50, 0., 500.);
  addHisto("eta_jet_40_b4jet", cutflowType+"/nvtx_40", 50, -5.0, 5.0);
  addHisto("bdiscr_40_b4jet", cutflowType+"/nvtx_40", 50, 0, 2.0);
  //after 4 jet cut
  addHisto("multi_jet_40_a4jet", cutflowType+"/nvtx_40", 100, 0., 20.);
  addHisto("pt_jet_40_a4jet", cutflowType+"/nvtx_40", 50, 0., 500.);
  addHisto("eta_jet_40_a4jet", cutflowType+"/nvtx_40", 50, -5.0, 5.0);
  addHisto("pt_met_40_a4jet", cutflowType+"/nvtx_40", 50, 0., 500.);
  addHisto("bdiscr_40_a4jet", cutflowType+"/nvtx_40", 50, 0, 2.0);
  addHisto("multi_bjet_40_a4jet", cutflowType+"/nvtx_40", 100, 0., 20.);
  //after btag
  addHisto("pt_bjet_40_a4jet", cutflowType+"/nvtx_40", 50, 0., 500.);
  addHisto("eta_bjet_40_a4jet", cutflowType+"/nvtx_40", 50, -5.0, 5.0);
  addHisto("pt_met_40_a4jet_btag", cutflowType+"/nvtx_40", 50, 0., 500.);
  addHisto("bdiscr_40_a4jet_btag", cutflowType+"/nvtx_40", 50, 0, 2.0);
  
  //NVTX FROM 40 To 50
  InitHist("nvtx_50", cutflowType, outFile_, true); 
  //befor 1 mu cut
  addHisto("multi_mu_50_b1mu", cutflowType+"/nvtx_50", 100, 0., 20.);
  addHisto("pt_mu_50_b1mu", cutflowType+"/nvtx_50", 50, 0., 500.);
  addHisto("pt_met_50_b1mu", cutflowType+"/nvtx_50", 50, 0., 500.);
  addHisto("eta_mu_50_b1mu", cutflowType+"/nvtx_50", 50, -5.0, 5.0);
  //after 1 mu cut
  addHisto("pt_mu_50_a1mu", cutflowType+"/nvtx_50", 50, 0., 500.);
  addHisto("pt_met_50_a1mu", cutflowType+"/nvtx_50", 50, 0., 500.);
  addHisto("eta_mu_50_a1mu", cutflowType+"/nvtx_50", 50, -5.0, 5.0);
  //before 4 jet cut
  addHisto("multi_jet_50_b4jet", cutflowType+"/nvtx_50", 100, 0., 20.);
  addHisto("pt_jet_50_b4jet", cutflowType+"/nvtx_50", 50, 0., 500.);
  addHisto("pt_met_50_b4jet", cutflowType+"/nvtx_50", 50, 0., 500.);
  addHisto("eta_jet_50_b4jet", cutflowType+"/nvtx_50", 50, -5.0, 5.0);
  addHisto("bdiscr_50_b4jet", cutflowType+"/nvtx_50", 50, 0, 2.0);
  //after 4 jet cut
  addHisto("multi_jet_50_a4jet", cutflowType+"/nvtx_50", 100, 0., 20.);
  addHisto("pt_jet_50_a4jet", cutflowType+"/nvtx_50", 50, 0., 500.);
  addHisto("pt_met_50_a4jet", cutflowType+"/nvtx_50", 50, 0., 500.);
  addHisto("eta_jet_50_a4jet", cutflowType+"/nvtx_50", 50, -5.0, 5.0);
  addHisto("bdiscr_50_a4jet", cutflowType+"/nvtx_50", 50, 0, 2.0);
  addHisto("multi_bjet_50_a4jet", cutflowType+"/nvtx_50", 100, 0., 20.);
  //after btag
  addHisto("pt_bjet_50_a4jet", cutflowType+"/nvtx_50", 50, 0., 500.);
  addHisto("eta_bjet_50_a4jet", cutflowType+"/nvtx_50", 50, -5.0, 5.0);
  addHisto("pt_met_50_a4jet_btag", cutflowType+"/nvtx_50", 50, 0., 500.);
  addHisto("bdiscr_50_a4jet_btag", cutflowType+"/nvtx_50", 50, 0, 2.0);
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
  InitHist("BTag", cutflowType, outFile_, false); 
  addHisto("Pre_RelIso",cutflowType+"/BTag", 40,0,0.5);
  addHisto("Final_RelIso",cutflowType+"/BTag", 40,0,0.5);
  addHisto("Muon_mult_final",cutflowType+"/BTag", 10,0,10);
  addHisto("final_multi_jet", cutflowType+"/BTag", 100,0,10);
  addHisto("CSVL_count", cutflowType+"/BTag", 10,0,10);
  addHisto("CSVM_count", cutflowType+"/BTag", 10,0,10);
  addHisto("wmt", cutflowType+"/BTag", 50, 0., 200.);
  addHisto("nvtx", cutflowType+"/BTag", 50, 0., 50.);
/*
  InitHist("KinFit", cutflowType, outFile_);
  addHisto("Pre_RelIso",cutflowType+"/KinFit", 40,0,0.5);
  addHisto("Final_RelIso",cutflowType+"/KinFit", 40,0,0.5);
  addHisto("Muon_mult_final",cutflowType+"/KinFit", 10,0,10);
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
  void HistogramPlotter::InitHist(TString dirname, TString parentDir, TFile *file, bool isNvtx=false)
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
  if(!isNvtx){
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
  addHisto("btag_jet", fullname, 50, -10., 10.);
  h1 = getHisto("btag_jet", fullname);
  h1->SetDirectory(d);
  addHisto("btagmulti_jet", fullname, 100, 0., 10.);
  h1 = getHisto("btagmulti_jet", fullname);
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
  }
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

