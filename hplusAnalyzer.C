// wjet_inc, zjet_inc, qcd, stop-all, diboson (cat migration unc)
#include <iomanip>
#include <iostream>
#include <fstream>

#include "TRandom2.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TH2.h"
#include "TTimeStamp.h"
#include "Math/VectorUtil.h"

#include "interface/Reader.h"
#include "interface/ObjectSelector.hh"
#include "interface/MomentumVec.h"
#include "interface/LumiReweighting.h"
#include "interface/UncertaintyComputer.hh"
#include "interface/HistogramPlotter.hh"

using namespace std;

class hplusAnalyzer : public ObjectSelector, HistogramPlotter
{
public :
  hplusAnalyzer() : ObjectSelector(), HistogramPlotter()
  {
    DRMIN_JET = 0.5;
    DRMIN_ELE = 0.5;
    METCUT_   = 30.0;
    //LumiWeights_ = reweight::LumiReWeighting("MC_Pileup_Summer2012_600bins.root","Data_Pileup_2012B_600bins.root", "pileup", "pileup");
    LumiWeights_ = reweight::LumiReWeighting("mcPileup.root","dataPileup.root", "pileup", "pileup");
    PShiftDown_ = reweight::PoissonMeanShifter(-0.5);
    PShiftUp_ = reweight::PoissonMeanShifter(0.5);
    //MC cross sections at 13 TeV: 
    //https://github.com/BristolTopGroup/AnalysisSoftware/blob/master/python/DataSetInfo_13TeV.py
    //https://github.com/BristolTopGroup/AnalysisSoftware/blob/master/python/DataSetInfo_13TeV_25ns.py
    //https://indico.cern.ch/event/617002/contributions/2490586/attachments/1419016/2173704/update_27022017.pdf
    xss["TTJets"]            =  831.76;
    xss["ST_tW"]             =  35.6; 
    xss["ST_t"]              =  136.2; 
    xss["ST_s"]              =  7.3;
    xss["WJetsToLNu"]        =  61526.7;  
    xss["W1JetsToLNu"]       =  9493;
    xss["W2JetsToLNu"]       =  3120;
    xss["W3JetsToLNu"]       =  942.3;
    xss["W4JetsToLNu"]       =  524.2;
    xss["DYJetsToLL"]        =  6025.2;
    xss["DY1JetsToLL"]       =  1016;
    xss["DY2JetsToLL"]       =  331.3;
    xss["DY3JetsToLL"]       =  96.6;
    xss["DY4JetsToLL"]       =  51.4;
    xss["WW"]                =  63.21;
    xss["WZ"]                =  22.82;
    xss["ZZ"]                =  10.32;
    xss["QCD_Pt-15to20_Mu"]  =  3819570;
    xss["QCD_Pt-20to30_EM"]  =  5352960;
    xss["QCD_Pt-120to170_EM"]=  62964;
    xss["QCD_Pt-170to300_EM"]=  18810;
    xss["HplusM80"]          =  1;
    xss["HplusM90"]          =  1;
    xss["HplusM100"]         =  1;
    xss["HplusM120"]         =  1; 
    xss["HplusM140"]         =  1;
    xss["HplusM150"]         =  1;
    xss["HplusM155"]         =  1;
    xss["HplusM160"]         =  1;
    
 xss["WJETS"] = 36257.0;  
    //xss["TTBAR"] = 245.8; 
    xss["ZJETS"] = 3504.0; 
    xss["QCD"] = 134680;
    //stop sample 
    xss["TOPS"]  = 3.79; 
    xss["TOPT"]  = 56.4; 
    xss["TOPW"]  = 11.1;
    //sbar sample 
    xss["TBARS"]  = 1.76; 
    xss["TBART"]  = 30.7; 
    xss["TBARW"]  = 11.1;
    //di-boson samples
    xss["WW"] = 33.61; 
    xss["WZ"] = 12.63; 
    xss["ZZ"] = 5.196;
    //signal sample 
    xss["WH"] = 245.8; xss["HH"] = 245.8;

    //muon Trigger/ID/ISo SFs, in bins of eta (from muon POG)
    //SFs for different lumi period are weighted by lumi fraction.
                                                                                           
    double lumiA = 808.411; double lumiB = 4044.0;
    double lumiC = 495.003+6432.0; double lumiD = 7274.0;
    double lumiTotal = lumiA+lumiB+lumiC+lumiD;
    //Trigger SF for HLT_IsoMu24_eta2p1
    double sfEta1 = (lumiA*0.956+lumiB*0.9798+lumiC*0.9841+lumiD*0.98151)/lumiTotal; // 0<|eta|<0.9 
    double sfEta2 = (lumiA*0.9528+lumiB*0.9618+lumiC*0.9688+lumiD*0.96156)/lumiTotal; // 0.9<|eta|<1.2
    double sfEta3 = (lumiA*0.9809+lumiB*0.9814+lumiC*1.0021+lumiD*0.99721)/lumiTotal; // 1.2<|eta|<2.1
    //multiply mu ID/Iso SFs
    sfEta1 = sfEta1*0.9939*1.0004;
    sfEta2 = sfEta2*0.9902*1.0031;
    sfEta3 = sfEta3*0.9970*1.0050;
    muSF["sfEta1"] = sfEta1;
    muSF["sfEta2"] = sfEta2;
    muSF["sfEta3"] = sfEta3;
  };
  ~hplusAnalyzer() {
    delete evR;
  };
  
  void CutFlowAnalysis(TString url,  string myKey="PFlow", bool isData = true, string evtType="data");
  void CutFlowProcessor(TString url,  string myKey="PFlow", TString cutflowType="base", bool isData = true, TFile *outFile_=0);
  //void CreateAnalHistos(TString flowType, TFile* outFile_);
  void processEvents();
  float reweightHEPNUPWJets(int hepNUP);
  float reweightHEPNUPDYJets(int hepNUP);

private :
  double DRMIN_JET, DRMIN_ELE, METCUT_;
  Reader *evR;
  
  reweight::LumiReWeighting LumiWeights_;
  reweight::PoissonMeanShifter PShiftUp_;   //pileup syst up
  reweight::PoissonMeanShifter PShiftDown_; //pileup syst down 
  std::map<string, double> xss;
  std::map<string, double> muSF;
  ofstream outfile_;
};
void hplusAnalyzer::CutFlowAnalysis(TString url, string myKey, bool isData, string evtType){
  
  TString outFile("13TeV/outputDir/");
  //TString outFile("");
  TString Filename_ = outFile+evtType+"_Anal.root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel( 9 );
  
  TString debug_Filename_ = Filename_+"_debug.txt";
  string debug_file(debug_Filename_);
  outfile_.open(debug_file.c_str());

  CutFlowProcessor(url, myKey, "base", isData, outFile_);
  if(not isData){
    
//     CutFlowProcessor(url, myKey, "JESPlus", isData, outFile_, Cat);
//     CutFlowProcessor(url, myKey, "JESMinus", isData, outFile_, Cat);
//     CutFlowProcessor(url, myKey, "JERPlus", isData, outFile_, Cat);
//     CutFlowProcessor(url, myKey, "JERMinus", isData, outFile_, Cat);
//     CutFlowProcessor(url, myKey, "METUCPlus", isData, outFile_, Cat);
//     CutFlowProcessor(url, myKey, "METUCMinus", isData, outFile_, Cat);
//     CutFlowProcessor(url, myKey, "bTagPlus", isData, outFile_, Cat);
//     CutFlowProcessor(url, myKey, "bTagMinus", isData, outFile_, Cat);
//     CutFlowProcessor(url, myKey, "TopPtPlus", isData, outFile_, Cat);
//     CutFlowProcessor(url, myKey, "TopPtMinus", isData, outFile_, Cat);

  }
  outfile_.close();
  outFile_->Write(); 
  outFile_->Close(); 
}
  
void hplusAnalyzer::CutFlowProcessor(TString url,  string myKey, TString cutflowType, bool isData, TFile *outFile_){
    int input_count = 0;
  
  outfile_<<"///// Begin processing "<<cutflowType<<"  selection ///////"<<endl;
  cout <<"///// Begin processing "<<cutflowType<<"  selection ///////"<<endl;
  
  bool isPFlow = (myKey.find("PFlow") != string::npos) ? true : false;

  string eAlgo("Electrons"), mAlgo("Muons"), jAlgo("Jets"), metAlgo("mets");
  //  if(myKey.find("PFlow") != string::npos)jAlgo = "PFlow";
  //  else if(myKey.find("HPS") != string::npos)jAlgo = "JetsAK5PF";
  
  //Uncertainty variations, JES, JER, MET unclustered, bTag
  ///int jes = 0, jer = 0, metuc = 0, bscale = 0;
  int jes = 1, jer = 1, metuc = 1, bscale = 0;
  if(cutflowType.Contains("JESPlus"))jes = 1;
  else if (cutflowType.Contains("JESMinus"))jes = -1;
  else if (cutflowType.Contains("JERPlus"))jer = 1;
  else if (cutflowType.Contains("JERMinus"))jer = -1;
  else if (cutflowType.Contains("METUCPlus"))metuc = 1;
  else if (cutflowType.Contains("METUCMinus"))metuc = -1;
  else if (cutflowType.Contains("bTagPlus"))bscale = 1;
  else if (cutflowType.Contains("bTagMinus"))bscale = -1; 

  evR = new Reader();
  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }

  int nEntries = evR->AssignEventTreeFrom(f);
  if( nEntries == 0) {return; }
  //get initial number of events
  TH1F* inputcf = (TH1F*)(f->Get("allEventsFilter/totalEvents"))->Clone("inputcf");
  double initialEvents = inputcf->GetBinContent(1);

  cout<<"Input file : "<<url<<endl;
  outfile_<<"Input file : "<<url<<endl;
  outfile_<<"Available input sample : "<<initialEvents<<endl;
  //define histograms 
  CreateAnalHistos(cutflowType, outFile_);
  
  double sampleWeight(1.0);
  //cross-sections (pb)
  //double Lumi = 34390.0;
  double Lumi = 2000.0;
  //double sigma_TTJets = 831.76;
  //double sigma_HplusM120 = 1;
  //sampleWeight = sigma_TTJets *Lumi /1000000;
  //sampleWeight = sigma_TTJets *Lumi / initialEvents;
  ///sampleWeight = sigma_HplusM120 *Lumi / initialEvents;
  outfile_<<"sampleWeight = "<<sampleWeight<<endl;

  MyEvent *ev;
  int nTriggEvent = 0, nSelEvents = 0, matchjetcount= 0, threepairjet = 0;
  double nVerticesFailCount = 0.0;
  double corr_jet_pair_cs = 0.0, corr_jet_pair_ud = 0.0;
  double matchedJet_q = 0.0, matchedJet_b = 0.0, matched_quark_eta_pt = 0.0, not_matchedJet_q = 0.0;
  double TotalTopPtWeights = 0, TotalLplusJEvents = 0; 
  double TotalSVEff = 0, TotalEvtsWithSV = 0.;
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;
    ev = evR->GetNewEvent(i);
    if(ev==0) continue;
    
    //print event number
    if(i%1000==0) cout<<"\033[01;32mEvent number = \033[00m"<< i << endl;
 
    //cout<<"sampleWeight basic = "<<sampleWeight<<endl;
    // apply PU re-weighting
    double evtWeight = 1.0;
    if(!isData){
      //get sample information
      string sampleName = ev->sampleInfo.sampleName;
      //cout<<"sampleName = "<<sampleName<<endl;
      if(sampleName.find("WJETS") != string::npos || sampleName.find("W1JETS") != string::npos || sampleName.find("W2JETS") != string::npos || sampleName.find("W3JETS") != string::npos || sampleName.find("W4JETS") != string::npos){
	int hepNUP = ev->sampleInfo.hepNUP;
        sampleWeight = reweightHEPNUPWJets(hepNUP) * (Lumi/1000.0);
      }
      else if(sampleName.find("ZJETS") != string::npos || sampleName.find("Z1JETS") != string::npos || sampleName.find("Z2JETS") != string::npos || sampleName.find("Z3JETS") != string::npos || sampleName.find("Z4JETS") != string::npos){
        int hepNUP = ev->sampleInfo.hepNUP;
        sampleWeight = reweightHEPNUPDYJets(hepNUP) * (Lumi/1000.0);
      }
      if(i < 1){
        sampleWeight = xss[sampleName] * Lumi / initialEvents;
        outfile_<<"Scale factor for lumi "<<Lumi<<" pb is "<< sampleWeight<<endl;
      }
      
      evtWeight *= sampleWeight; // upto this only sigma*lumi weight is applied
      //vector<double>puweights = ev->sampleInfo.puWeights;  // this line and below this line applies the pile up corrections
      //if(puweights.size() > 0) evtWeight *= puweights[0];

      vector<double>pu = ev->sampleInfo.truepileup;
      if(pu.size() > 0) {
	float npu = pu[0];
	double weight = LumiWeights_.weight(npu);
	//        if (cutflowType.Contains("PUPlus"))weight = weight*PShiftUp_.ShiftWeight( npu );
	//        else if (cutflowType.Contains("PUMinus"))weight = weight*PShiftDown_.ShiftWeight( npu );
	//cout<<"pu weight = "<<weight<<endl;
    evtWeight *= weight;  // pile up weight is also applied
      }
    }

    //Apply Top re-weighting weights here (after pile pu weights).                       
    double topPtWeights_offline = 1.0;
    if(!isData){
      //get sample information
      string sampleName = ev->sampleInfo.sampleName;

      if(sampleName.find("TTJets") != string::npos || sampleName.find("WH") != string::npos || sampleName.find("HplusM120") != string::npos){

        vector<double>topptweights = ev->sampleInfo.topPtWeights;
        if(topptweights.size() > 0){
          topPtWeights_offline = topptweights[0];
          if(cutflowType.Contains("TopPtPlus"))
            topPtWeights_offline = topPtWeights_offline*topPtWeights_offline;
          else if(cutflowType.Contains("TopPtMinus"))
            topPtWeights_offline = 1.0;
        }
      }
    }
    //    outfile_<< "topPtWeights_offline" << topPtWeights_offline << endl;

    TotalTopPtWeights += topPtWeights_offline;
    //cout<<"topPtWeights_offline = "<<topPtWeights_offline<<endl;
    TotalLplusJEvents++;
    evtWeight *= topPtWeights_offline; //Multiply to the total weights
    //cout<<"evtWeight = "<<evtWeight<<endl;
    //cout<<endl;
    //double evtWeight = 1.0;
    /// for data, evtWeight = 1
    //evtWeight *= sampleWeight;
    
    double nCutPass = 0.0;// double nCutPass_plus = 0.0; double nCutPass_minus = 0.0;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
    //    fillHisto("cutflow_plus", cutflowType, nCutPass_plus, evtWeight);
    //    fillHisto("cutflow_minus", cutflowType, nCutPass_minus, evtWeight);
    
    bool passTrig = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("Mu") != string::npos) passTrig = true;
    }
    if(!passTrig){
            //cout << "not satisfying trigger" << endl;
      continue;
    }
    nTriggEvent++;
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
   
    //get objects //
    vector<MyVertex> Vertices = ev->PrimaryVtxs;
    if(Vertices.size() <= 0){
      nVerticesFailCount+=evtWeight;
      cout<<" no vertexes , exit"<<endl;
      continue;
    }
    
    vector<MyMuon> pfMuons = evR->getMuons(ev, mAlgo);
    vector<MyMuon> pfMuons_noiso = evR->getMuons(ev, mAlgo);
    vector<MyElectron> pfElectrons = evR->getElectrons(ev, eAlgo);
    vector<MyJet> pfJets = evR->getJets(ev, jAlgo);
    MyMET met = evR->getMET(ev, metAlgo);
    ///vector<MyTau>pfTaus; pfTaus.clear(); 

    // preselect objects 
    vector<int> m_init; m_init.clear();
    preSelectMuons(&m_init, pfMuons, Vertices[0], isPFlow);
    vector<int> m_init_noiso; m_init_noiso.clear();
    preSelectMuonsNoIso(&m_init_noiso, pfMuons_noiso, Vertices[0], isPFlow);
    vector<int> e_init; e_init.clear();
    preSelectElectrons(&e_init, pfElectrons, Vertices[0], isPFlow);
    vector<int> j_init; j_init.clear();
    preSelectJets(jAlgo, &j_init, pfJets, jes, jer);
    
    // clean objects //
    vector<int> e_final; e_final.clear();
    ElectronCleaning( pfElectrons, pfMuons, &e_init, &e_final, &m_init, DRMIN_ELE);
    vector<int> j_final; j_final.clear();
    vector<int> t_final; t_final.clear();
    JetCleaning(pfJets, pfMuons, pfElectrons,  &j_init, &j_final, &m_init, &e_final, DRMIN_JET);
    
    //Get MC partons
    vector<MyLorentzVector> bquarks; bquarks.clear();
    vector<MyLorentzVector> lquarks; lquarks.clear();
    MyLorentzVector mcTop, mcAntiTop;
    mcTop.SetPxPyPzE(0., 0., 0., 0.); mcAntiTop.SetPxPyPzE(0., 0., 0., 0.);
    if(!isData){
      vector<MyMCParticle>allMCParticles = ev->mcParticles;
      for(size_t imc=0; imc < allMCParticles.size(); ++imc){
	if(abs(allMCParticles[imc].pid) == 5 && allMCParticles[imc].mother.size() > 0 && (abs(allMCParticles[imc].mother[0])==6) )
	  bquarks.push_back(allMCParticles[imc].p4Gen);
	else if(abs(allMCParticles[imc].pid) <= 4 && allMCParticles[imc].mother.size() > 0 && (abs(allMCParticles[imc].mother[0])==24 || abs(allMCParticles[imc].mother[0])==37) )
	  lquarks.push_back(allMCParticles[imc].p4Gen); 
	else if(allMCParticles[imc].pid == 6 && allMCParticles[imc].status == 3)
	  mcTop = allMCParticles[imc].p4Gen;
	else if(allMCParticles[imc].pid == -6 && allMCParticles[imc].status == 3)
	  mcAntiTop = allMCParticles[imc].p4Gen;
      }
    }
    
    ///
    //Apply Lepton selection//////////////////////////////
    if(m_init_noiso.size() > 0){
      int m_i_noiso = m_init_noiso[0];
      double mRelIso_no_iso = pfMuons_noiso[m_i_noiso].pfRelIso;
      fillHisto("Pre_RelIso",cutflowType, mRelIso_no_iso, evtWeight);
    }
    double pri_vtxs = Vertices.size();
    fillHisto("nvtx_nocut", cutflowType, pri_vtxs, evtWeight);
    int nLepton = m_init.size();  // this condition proof that only muon + jet events
    //Variation with nvtx, before 1 muon ------------------------
    for(int v=0; v<6; v++){
      if(Vertices.size() >= 10*v && Vertices.size() <= 10*(v+1)){
        TString nvtx_dir = Form("nvtx_%d", 10*(v+1)); 
        fillHisto(Form("multi_mu_%d_b1mu", 10*(v+1)) , cutflowType+"/"+nvtx_dir, nLepton, evtWeight);
        fillHisto(Form("pt_met_%d_b1mu", 10*(v+1)), cutflowType+"/"+nvtx_dir, met.p4.pt(), evtWeight);
        for(int m=0; m<m_init.size(); m++){
          fillHisto(Form("pt_mu_%d_b1mu", 10*(v+1)), cutflowType+"/"+nvtx_dir, pfMuons[m].p4.pt(), evtWeight); 
          fillHisto(Form("eta_mu_%d_b1mu", 10*(v+1)) , cutflowType+"/"+nvtx_dir, pfMuons[m].p4.eta(), evtWeight);
        }
      }
    }//-----------------------
    if(nLepton != 1)continue;
    ///if( looseMuonVeto( m_init[0],pfMuons, isPFlow) ) continue;
    ///if( looseElectronVeto(-1,pfElectrons, isPFlow) ) continue;
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
    
    //Variation with nvtx, after 1 muon ------------------------
    for(int v=0; v<6; v++){
      if(Vertices.size() >= 10*v && Vertices.size() <= 10*(v+1)){
        TString nvtx_dir = Form("nvtx_%d", 10*(v+1)); 
        fillHisto(Form("pt_met_%d_a1mu", 10*(v+1)), cutflowType+"/"+nvtx_dir, met.p4.pt(), evtWeight);
        for(int m=0; m<m_init.size(); m++){
          fillHisto(Form("pt_mu_%d_a1mu", 10*(v+1)), cutflowType+"/"+nvtx_dir, pfMuons[m].p4.pt(), evtWeight); 
          fillHisto(Form("eta_mu_%d_a1mu", 10*(v+1)) , cutflowType+"/"+nvtx_dir, pfMuons[m].p4.eta(), evtWeight);
        }
      }
    }//-----------------------
    //apply muon SF to eventWeights 
    int m_i = m_init[0];
/* 
    double musfWeight = 1.0;
    if(fabs(pfMuons[m_i].p4.eta()) < 0.9)musfWeight = muSF["sfEta1"];
    else if(fabs(pfMuons[m_i].p4.eta()) > 0.9 && fabs(pfMuons[m_i].p4.eta()) < 1.2)musfWeight = muSF["sfEta2"];
    else musfWeight = muSF["sfEta3"];
    evtWeight *= musfWeight;
    //cout<<"evtWeight musf = "<<evtWeight<<endl;
*/   
    int count_muon = m_init.size();
    int muCharge = pfMuons[m_i].charge;

    ///double nCutPass = 0.0;// double nCutPass_plus = 0.0; double nCutPass_minus = 0.0;
    double mRelIso = pfMuons[m_i].pfRelIso;
    fillHisto("pt_mu", cutflowType, pfMuons[m_i].p4.pt(), evtWeight);
    fillHisto("Final_RelIso",cutflowType, mRelIso, evtWeight);
    fillHisto("Muon_mult_final",cutflowType, count_muon, evtWeight);

    //Apply Top re-weighting weights here (after lepton selection, such that we can use coefficients of l+jets events).
    ///double topPtWeights_offline = 1.0;
    if(!isData){
      //get sample information
      string sampleName = ev->sampleInfo.sampleName;
      if(sampleName.find("TTJets") != string::npos || sampleName.find("WH") != string::npos || sampleName.find("HplusM120") != string::npos){
	topPtWeights_offline = sqrt(exp(0.159 - 0.00141*mcTop.pt())*exp(0.159 - 0.00141*mcAntiTop.pt()));
	if(cutflowType.Contains("TopPtPlus"))
	  topPtWeights_offline = topPtWeights_offline*topPtWeights_offline;
	else if(cutflowType.Contains("TopPtMinus"))
	  topPtWeights_offline = 1.0;
      }
    }
    //    cout << "topPtWeights_offline" << topPtWeights_offline << endl;
    TotalTopPtWeights += topPtWeights_offline; 
    TotalLplusJEvents++;
    evtWeight *= topPtWeights_offline; //Multiply to the total weights

    // Fill histogram after trigger and one offline isolated muon and applied 2nd lepton veto
    int nJet = j_final.size();
    fillHisto("eta_mu", cutflowType, pfMuons[m_i].p4.eta(), evtWeight);
    fillHisto("phi_mu", cutflowType, pfMuons[m_i].p4.phi(), evtWeight);
    // vertex just after one lepton selection
    //double pri_vtxs = Vertices.size();
    fillHisto("nvtx", cutflowType, pri_vtxs, evtWeight);
    fillHisto("nvtx_1mu", cutflowType, pri_vtxs, evtWeight);

    // Apply Jet Selection
    //Variation with nvtx, before 4 jet ------------------------
    ///int threejet = 0;
    ///float highestJetPt = 0;
    for(int v=0; v<6; v++){
      if(Vertices.size() >= 10*v && Vertices.size() <= 10*(v+1)){
        TString nvtx_dir = Form("nvtx_%d", 10*(v+1)); 
        fillHisto(Form("multi_jet_%d_b4jet", 10*(v+1)), cutflowType+"/"+nvtx_dir, j_final.size(), evtWeight);
        fillHisto(Form("pt_met_%d_b4jet", 10*(v+1)), cutflowType+"/"+nvtx_dir, met.p4.pt(), evtWeight);
        double bDiscr = 0.0;
        for(size_t ijet = 0; ijet < j_final.size(); ijet++){
          int ind_jet = j_final[ijet];
          bDiscr = pfJets[ind_jet].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
          //double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer); 
          double jetPt = pfJets[ind_jet].p4.pt(); 
          ///fillHisto("pt_jet", cutflowType, jetPt, evtWeight);
          ///fillHisto("eta_jet", cutflowType, pfJets[ind_jet].p4.eta(), evtWeight);
          ///fillHisto("phi_jet", cutflowType, pfJets[ind_jet].p4.phi(), evtWeight);
          ///if(jetPt > 30 )threejet++;
          ///if(jetPt > highestJetPt)highestJetPt=jetPt;
          fillHisto(Form("pt_jet_%d_b4jet", 10*(v+1)), cutflowType+"/"+nvtx_dir, jetPt, evtWeight); 
          fillHisto(Form("eta_jet_%d_b4jet", 10*(v+1)) , cutflowType+"/"+nvtx_dir, pfJets[ind_jet].p4.eta(), evtWeight);
          fillHisto(Form("bdiscr_%d_b4jet", 10*(v+1)), cutflowType+"/"+nvtx_dir, bDiscr, evtWeight); 
        }
      }
    }//-----------------------

    fillHisto("multi_jet", cutflowType, nJet, evtWeight);
    if(nJet < 4)continue;  // this condition implies event should contain at least 4 jets
    fillHisto("final_multi_jet", cutflowType, nJet, evtWeight);
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
    fillHisto("nvtx_1mu_4jet", cutflowType, pri_vtxs, evtWeight);

    //Variation with nvtx, after 4 jet ------------------------
    int threejet = 0;
    float highestJetPt = 0;
    for(int v=0; v<6; v++){
      if(Vertices.size() >= 10*v && Vertices.size() <= 10*(v+1)){
        TString nvtx_dir = Form("nvtx_%d", 10*(v+1)); 
        fillHisto(Form("multi_jet_%d_a4jet", 10*(v+1)), cutflowType+"/"+nvtx_dir, j_final.size(), evtWeight);
        fillHisto(Form("pt_met_%d_a4jet", 10*(v+1)), cutflowType+"/"+nvtx_dir, met.p4.pt(), evtWeight);
        double bDiscr = 0.0;
        for(size_t ijet = 0; ijet < j_final.size(); ijet++){
          int ind_jet = j_final[ijet];
          bDiscr = pfJets[ind_jet].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
          //double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer); 
          double jetPt = pfJets[ind_jet].p4.pt(); 
          fillHisto("pt_jet", cutflowType, jetPt, evtWeight);
          fillHisto("eta_jet", cutflowType, pfJets[ind_jet].p4.eta(), evtWeight);
          fillHisto("phi_jet", cutflowType, pfJets[ind_jet].p4.phi(), evtWeight);
          if(jetPt > 30 )threejet++;
          if(jetPt > highestJetPt)highestJetPt=jetPt;
          fillHisto(Form("pt_jet_%d_a4jet", 10*(v+1)), cutflowType+"/"+nvtx_dir, jetPt, evtWeight); 
          fillHisto(Form("eta_jet_%d_a4jet", 10*(v+1)) , cutflowType+"/"+nvtx_dir, pfJets[ind_jet].p4.eta(), evtWeight);
          fillHisto(Form("bdiscr_%d_a4jet", 10*(v+1)), cutflowType+"/"+nvtx_dir, bDiscr, evtWeight); 
        }
      }
    }//-----------------------
    
    
    //    if(threejet < 3) continue; // three jets should have pt > 30 GeV do not need for ATLAS selection
    // Met distribution   
    double   leptonPt(0), deltaPhi(0);
    double metPt = 0; //met.p4.pt();
    if(!metuc)metPt = metWithJESJER(pfJets, &j_final, met, jes, jer);
    else metPt = metWithUncl(pfJets, &j_final, pfMuons, &m_init, pfElectrons, &e_final, met, metuc);
    fillHisto("pt_met", cutflowType, metPt, evtWeight);
//-----------------    
    ////if(metPt < 20) continue;  // Missing transverse energy cut 30 GeV(CMS) for ATLAS 20 GeV 
    fillHisto("final_pt_met", cutflowType, metPt, evtWeight);
    ///nCutPass++;
    ///fillHisto("cutflow", cutflowType, nCutPass, evtWeight);

    int m_j = m_init[0];
    leptonPt = TMath::Abs(pfMuons[m_j].p4.pt());
    deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(pfMuons[m_j].p4, met.p4);
    double mt = sqrt (  2*leptonPt*metPt*(1 - cos(deltaPhi) ) ) ;
    fillHisto("wmt", cutflowType, mt, evtWeight);
    /*
    if((mt + metPt ) < 30 ) continue;  // extra condition to reduce multijet bkgs
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
    */
    
    // Here putting the criteria that event should contain at least one b-tagged jet  
    int count_CSVL = 0; int count_CSVM = 0; //int count_CSVT = 0; 
    vector<int> j_final_nob; j_final_nob.clear();
    double bDiscr = 0.0 ;
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      bDiscr = pfJets[ind_jet].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
      //cout<<"bDiscr = "<<bDiscr<<endl;
      if(bDiscr > 0.5426){
        count_CSVL++;
        fillHisto("btag_jet", cutflowType, bDiscr, evtWeight); 
      }
      else j_final_nob.push_back(ind_jet);  
    }
    fillHisto("CSVL_count", cutflowType, count_CSVL, evtWeight);
    
    //Get the invariant mass of csbar pair
    if(j_final_nob.size() >= 2){
      int first_index = j_final_nob[0];
      int sec_index = j_final_nob[1];
      fillHisto("pt_cjet", cutflowType, pfJets[first_index].p4.pt(), evtWeight);
      fillHisto("pt_sjet", cutflowType, pfJets[sec_index].p4.pt(), evtWeight);
      MyLorentzVector diJet = pfJets[first_index].p4 + pfJets[sec_index].p4;
      //cout<<"diJet.M() = "<<diJet.M()<<endl;
      fillHisto("diJet_Mass", cutflowType, diJet.M(), evtWeight);
    }
    /*
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
            if(pfJets[ind_jet].bDiscriminator["combinedSecondaryVertexBJetTags"] > 0.679)count_CSVM++;
      bool isBtag = getBtagWithSF(pfJets[ind_jet], isData, bscale, true); 
      if(isBtag)count_CSVM++; 
    }// 13TeV, it is 0.8484, the tight one is 0.9535 

    //    if(count_CSVL <= 2) continue; // Atfirst demanding for ATLAS 3L and then 0M
    */
    if(count_CSVL <= 1) continue; // Demanding for 2M b-tagged jets
    fillHisto("CSVM_count", cutflowType, count_CSVM, evtWeight);
    //    fillHisto("btagmulti_jet", cutflowType, count, evtWeight);
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight); 
    nSelEvents++;  // this is the counter for counting final number of events  
    fillHisto("nvtx_1mu_4jet_btag", cutflowType, pri_vtxs, evtWeight);

    //Variation with nvtx, after b tagging ------------------------
    ///int threejet = 0;
    ///float highestJetPt = 0;
    for(int v=0; v<6; v++){
      if(Vertices.size() >= 10*v && Vertices.size() <= 10*(v+1)){
        TString nvtx_dir = Form("nvtx_%d", 10*(v+1)); 
        fillHisto(Form("multi_bjet_%d_a4jet", 10*(v+1)), cutflowType+"/"+nvtx_dir, j_final.size(), evtWeight);
        fillHisto(Form("pt_met_%d_a4jet_btag", 10*(v+1)), cutflowType+"/"+nvtx_dir, met.p4.pt(), evtWeight);
        double bDiscr = 0.0 ;
        for(size_t ijet = 0; ijet < j_final.size(); ijet++){
            int ind_jet = j_final[ijet];
            //double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer); 
            double jetPt = pfJets[ind_jet].p4.pt(); 
            fillHisto(Form("pt_bjet_%d_a4jet", 10*(v+1)), cutflowType+"/"+nvtx_dir, jetPt, evtWeight); 
            fillHisto(Form("eta_bjet_%d_a4jet", 10*(v+1)) , cutflowType+"/"+nvtx_dir, pfJets[ind_jet].p4.eta(), evtWeight);
            bDiscr = pfJets[ind_jet].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
            //cout<<"bDiscr = "<<bDiscr<<endl;
            if(bDiscr > 0.5426){
            fillHisto(Form("bdiscr_%d_a4jet_btag", 10*(v+1)), cutflowType+"/"+nvtx_dir, bDiscr, evtWeight); 
            }
        }
      }
    }//-----------------------
   
    // add set of plots after BTag:
    fillHisto("pt_mu", cutflowType+"/BTag", pfMuons[m_i].p4.pt(), evtWeight);
    fillHisto("eta_mu", cutflowType+"/BTag", pfMuons[m_i].p4.eta(), evtWeight);
    fillHisto("phi_mu", cutflowType+"/BTag", pfMuons[m_i].p4.phi(), evtWeight);
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      //double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      double jetPt = pfJets[ind_jet].p4.pt(); 
      fillHisto("pt_jet", cutflowType+"/BTag", jetPt, evtWeight);
      fillHisto("eta_jet", cutflowType+"/BTag", pfJets[ind_jet].p4.eta(), evtWeight);
      fillHisto("phi_jet", cutflowType+"/BTag", pfJets[ind_jet].p4.phi(), evtWeight);
    }
    fillHisto("final_multi_jet", cutflowType+"/BTag", nJet, evtWeight);
    fillHisto("final_pt_met", cutflowType+"/BTag", metPt, evtWeight);
    fillHisto("nvtx", cutflowType+"/BTag", pri_vtxs, evtWeight);
    fillHisto("wmt", cutflowType+"/BTag", mt, evtWeight);

        input_count++;
        if(input_count%10==0)
        cout << "input count: "<< input_count << endl;
        //if(input_count > 100000) break;
        if(i > 5000) break;
      }//event loop
      
      outfile_ << "Number of times HadP and HadQ matches with jets" << matchjetcount << endl;
      outfile_ << "total number of selected events " << nSelEvents <<endl; 
      outfile_ << "No of times three pair jet matched:    " << threepairjet << endl;
      outfile_ << "No of correct jet pair assign to W:   " << matchedJet_q << endl; 
      outfile_ << "No of not correct jet pair assign to W:   " << not_matchedJet_q << endl;
  outfile_ << "No of correct jet pair not assign to W but |eta| < 2.5 and pt > 25:   " << matched_quark_eta_pt << endl;
  outfile_ << "No of correct jet matched to b (from top) :   " << matchedJet_b << endl;
  outfile_ << "Average weights for top-reweighting : "<<TotalTopPtWeights<<"/"<<TotalLplusJEvents<<" = "<<TotalTopPtWeights/TotalLplusJEvents<<endl; 
  outfile_ << "Efficiency Uncertainty on the SV finding : "<<TotalSVEff<<"/"<<TotalEvtsWithSV<<" = "<<TotalSVEff/TotalEvtsWithSV<<endl;
  fillHisto("AvTopPtWeight", cutflowType, 1.0, TotalTopPtWeights/TotalLplusJEvents);
  fillHisto("SVEffUncert", cutflowType, 2.0, TotalSVEff/TotalEvtsWithSV);
  fillHisto("SVEffUncert", cutflowType, 1.0, TotalEvtsWithSV);
  fillHisto("SVEffUncert", cutflowType, 0.0, TotalSVEff);
  f->Close(); 
  delete f;
}

void hplusAnalyzer::processEvents(){ 

  CutFlowAnalysis("TTJets_ntuple_MuChannel.root", "PF", true, "TTJets_MuMC_check"); 
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", false, "outputFile");
  
  //MC samples
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/multicrab_29April17/MuMC_20170429/TTJets_MuMC_20170429/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TTJets_MuMC_20170429/170429_110042/0000/TTJets_MuMC_20170429_Ntuple_9.root", "PF", true, "TTJets_MuMC");
  
  //DATA samples
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/multicrab_29April17/MuData_20170429/MuRunHv3_MuData_20170429/SingleMuon/MuRunHv3_MuData_20170429/170429_110948/0000/MuRunHv3_MuData_20170429_Ntuple_9.root", "PF", true, "MuRunHv3");
  
  //for condor submission
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", true, "outputFile");
} 

float hplusAnalyzer::reweightHEPNUPWJets(int hepNUP) {

  int nJets = hepNUP-5;
  if(nJets==0)      return 0.515637279;
  else if(nJets==1) return 0.185061618;
  else if(nJets==2) return 0.062474528;
  else if(nJets==3) return 0.038296287;
  else if(nJets>=4) return 0.019044067;
  else return 1 ;
}

float hplusAnalyzer::reweightHEPNUPDYJets(int hepNUP){

  int nJets = hepNUP-5;
  if(nJets==0)      return 0.117971671;
  else if(nJets==1) return 0.022925801;
  else if(nJets==2) return 0.009107319;
  else if(nJets==3) return 0.005462804;
  else if(nJets>=4) return 0.004252079;
  else return 1 ;
}

