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
    /////////////// Pileup reweigting //////////////////
    //Pileup informations for Data: 
    //
    LumiWeights_ = reweight::LumiReWeighting("trueInTimePU_mcDY.root","trueMinBiasPU_dataMu.root", "pileup", "pileup");
    PShiftDown_ = reweight::PoissonMeanShifter(-0.5);
    PShiftUp_ = reweight::PoissonMeanShifter(0.5);
    
    /////////////// MC cross sections at 13 TeV: /////// 
    //https://github.com/BristolTopGroup/AnalysisSoftware/blob/master/python/DataSetInfo_13TeV.py
    //https://github.com/BristolTopGroup/AnalysisSoftware/blob/master/python/DataSetInfo_13TeV_25ns.py
    //https://indico.cern.ch/event/617002/contributions/2490586/attachments/1419016/2173704/update_27022017.pdf
    //evtDBS= event at Data Base Server i.e in DAS (https://cmsweb.cern.ch/das/).
    xss["DY1JetsToLL"]       =  1016;          evtDBS["DY1JetsToLL"]       =  62627174;
    xss["DY2JetsToLL"]       =  331.3;         evtDBS["DY2JetsToLL"]       =  19970551;
    xss["DY3JetsToLL"]       =  96.6;          evtDBS["DY3JetsToLL"]       =  5856110;
    xss["DY4JetsToLL"]       =  51.4;          evtDBS["DY4JetsToLL"]       =  4197868;
    xss["DYJetsToLL"]        =  6025.2;        evtDBS["DYJetsToLL"]        =  49144274;
    xss["HplusM100"]         =  1;             evtDBS["HplusM100"]         =  996170; 
    xss["HplusM120"]         =  1;             evtDBS["HplusM120"]         =  994498; 
    xss["HplusM140"]         =  1;             evtDBS["HplusM140"]         =  987730; 
    xss["HplusM150"]         =  1;             evtDBS["HplusM150"]         =  990645;
    xss["HplusM155"]         =  1;             evtDBS["HplusM155"]         =  952984;
    xss["HplusM160"]         =  1;             evtDBS["HplusM160"]         =  992264;
    xss["HplusM80"]          =  1;             evtDBS["HplusM80"]          =  976710;
    xss["HplusM90"]          =  1;             evtDBS["HplusM90"]          =  988480;
    xss["QCD_Pt-15to20_Mu"]  =  3819570;       evtDBS["QCD_Pt-15to20_Mu"]  =  4141251;
    xss["QCD_Pt-20to30_Mu"]  =  2960198;       evtDBS["QCD_Pt-20to30_Mu"]  =  31475157;
    xss["QCD_Pt-30to50_Mu"]  =  1652471;       evtDBS["QCD_Pt-30to50_Mu"]  =  29954815;
    xss["QCD_Pt-50to80_Mu"]  =  437504;        evtDBS["QCD_Pt-50to80_Mu"]  =  19806915;
    xss["QCD_Pt-80to120_Mu"] =  106033;        evtDBS["QCD_Pt-80to120_Mu"] =  13786971;
    xss["QCD_Pt-120to170_Mu"]=  25190;         evtDBS["QCD_Pt-120to170_Mu"]=  8042721;
    xss["QCD_Pt-170to300_Mu"]=  8654;          evtDBS["QCD_Pt-170to300_Mu"]=  7947159;
    xss["ST_s"]              =  7.3;           evtDBS["ST_s"]              =  2989199;
    xss["ST_t"]              =  136.2;         evtDBS["ST_t"]              =  38811017;
    xss["ST_tW"]             =  35.6;          evtDBS["ST_tW"]             =  6933094;
    xss["TTJets"]            =  831.76;        evtDBS["TTJets"]            =  10139950;   
    xss["W1JetsToLNu"]       =  9493;          evtDBS["W1JetsToLNu"]       =  45367044;
    xss["W2JetsToLNu"]       =  3120;          evtDBS["W2JetsToLNu"]       =  29878415;
    xss["W3JetsToLNu"]       =  942.3;         evtDBS["W3JetsToLNu"]       =  19798117;
    xss["W4JetsToLNu"]       =  524.2;         evtDBS["W4JetsToLNu"]       =  9170576;
    xss["WJetsToLNu"]        =  61526.7;       evtDBS["WJetsToLNu"]        =  29705748;
    xss["WW"]                =  63.21;         evtDBS["WW"]                =  994012;
    xss["WZ"]                =  22.82;         evtDBS["WZ"]                =  1000000;
    xss["ZZ"]                =  10.32;         evtDBS["ZZ"]                =  990064; 
    
    //Lumis(inverse pb) of single muon DATA at 13TeV
    //https://docs.google.com/spreadsheets/d/1lQyfcY0gnG_IgFrtnbBES_HV1M7ARQM9qCw01vsnxSk/edit?usp=sharing
    double lumiB = 5403; 
    double lumiC = 2395;
    double lumiD = 4255; 
    double lumiE = 4053;
    double lumiF = 3105;
    double lumiG = 7544;
    double lumiH = 8529+216;
    double lumiTotal = lumiB+ lumiC+ lumiD+ lumiE+ lumiF+ lumiG+ lumiH;
    
    //muon Trigger/ID/ISo SFs, in bins of eta (from muon POG)
    //SFs for different lumi period are weighted by lumi fraction.
    //Trigger SF for HLT_IsoMu24_eta2p1
    double sfEta1 = (lumiE*0.956+lumiB*0.9798+lumiC*0.9841+lumiD*0.98151)/lumiTotal; // 0<|eta|<0.9 
    double sfEta2 = (lumiE*0.9528+lumiB*0.9618+lumiC*0.9688+lumiD*0.96156)/lumiTotal; // 0.9<|eta|<1.2
    double sfEta3 = (lumiE*0.9809+lumiB*0.9814+lumiC*1.0021+lumiD*0.99721)/lumiTotal; // 1.2<|eta|<2.1
    //multiply mu ID/Iso SFs
    sfEta1 = sfEta1*0.9939*1.0004;
    sfEta2 = sfEta2*0.9902*1.0031;
    sfEta3 = sfEta3*0.9970*1.0050;
    muSF["sfEta1"] = 1; //sfEta1;
    muSF["sfEta2"] = 1; //sfEta2;
    muSF["sfEta3"] = 1; //sfEta3;
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
  std::map<string, double> evtDBS;
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
  outfile_.close();
  outFile_->Write(); 
  outFile_->Close(); 
}
  
void hplusAnalyzer::CutFlowProcessor(TString url,  string myKey, TString cutflowType, bool isData, TFile *outFile_){
    int input_count = 0;
  
  outfile_<<"///// Begin processing "<<cutflowType<<"  selection ///////"<<endl;
  cout <<"///// Begin processing "<<cutflowType<<"  selection ///////"<<endl;
  
  bool isPFlow = (myKey.find("PFlow") != string::npos) ? true : false;

  string eAlgo("Electrons"), mAlgo("Muons"), jAlgo("Jets"), metAlgo("METs");
  //  if(myKey.find("PFlow") != string::npos)jAlgo = "PFlow";
  //  else if(myKey.find("HPS") != string::npos)jAlgo = "JetsAK5PF";
  
  //Uncertainty variations, JES, JER, MET unclustered, bTag
  ///int jes = 0, jer = 0, metuc = 0, bscale = 0;
  int jes = 1, jer = 1, metuc = 1; 
  if(cutflowType.Contains("JESPlus"))jes = 1;
  else if (cutflowType.Contains("JESMinus"))jes = -1;
  else if (cutflowType.Contains("JERPlus"))jer = 1;
  else if (cutflowType.Contains("JERMinus"))jer = -1;
  else if (cutflowType.Contains("METUCPlus"))metuc = 1;
  else if (cutflowType.Contains("METUCMinus"))metuc = -1;
  //else if (cutflowType.Contains("bTagPlus"))bscale = 1;
  //else if (cutflowType.Contains("bTagMinus"))bscale = -1; 

  evR = new Reader();
  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }
 
  //---------------------------------------------------//
  //get initial number of events, from ntuples
  //store initial informations, in a txt file
  //---------------------------------------------------//
  double lumiTotal = 35500;
  int nEntries = evR->AssignEventTreeFrom(f);
  if( nEntries == 0) {return; }
  TH1F* inputcf = (TH1F*)(f->Get("allEventsFilter/totalEvents"))->Clone("inputcf");
  double initialEvents = inputcf->GetBinContent(1);
  cout<<"input file: "<<url<<endl;
  outfile_<<"input file: "<<url<<endl;
  outfile_<<"totalEvents: "<<initialEvents<<endl;
  CreateAnalHistos(cutflowType, outFile_);
  fillHisto("totalEvents", cutflowType, initialEvents, 1);
  MyEvent *ev;
  int nTriggEvent = 0, nSelEvents = 0, matchjetcount= 0, threepairjet = 0;
  double nVerticesFailCount = 0.0;
  ///double corr_jet_pair_cs = 0.0, corr_jet_pair_ud = 0.0;
  double matchedJet_q = 0.0, matchedJet_b = 0.0, matched_quark_eta_pt = 0.0, not_matchedJet_q = 0.0;
  double TotalTopPtWeights = 0, TotalLplusJEvents = 0; 
  double TotalSVEff = 0, TotalEvtsWithSV = 0.;
  double sampleWeight(1.0);
  
  //---------------------------------------------------//
  //loop over each event, of the ntuple
  //---------------------------------------------------//
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;
    ev = evR->GetNewEvent(i);
    if(ev==0) continue;
    if(i%1000==0) cout<<"\033[01;32mEvent number = \033[00m"<< i << endl;
 
    //---------------------------------------------------//
    //apply lumi, k factor and pileup weight
    //---------------------------------------------------//
    double evtWeight = 1.0;
    if(!isData){
      //lumi weight
      string sampleName = ev->sampleInfo.sampleName;
      sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
      if(i < 1){
        outfile_<<"Sample weight for "<<sampleName<<" is= "<< sampleWeight<<endl;
      }
      evtWeight *= sampleWeight; 
      //k factor weights for WJets samples (we have to update it ==========)
      /*
      if(sampleName.find("WJETS") != string::npos || sampleName.find("W1JETS") != string::npos || sampleName.find("W2JETS") != string::npos || sampleName.find("W3JETS") != string::npos || sampleName.find("W4JETS") != string::npos){
	int hepNUP = ev->sampleInfo.hepNUP;
        sampleWeight = reweightHEPNUPWJets(hepNUP) * (Lumi/1000.0);
      }
      else if(sampleName.find("ZJETS") != string::npos || sampleName.find("Z1JETS") != string::npos || sampleName.find("Z2JETS") != string::npos || sampleName.find("Z3JETS") != string::npos || sampleName.find("Z4JETS") != string::npos){
        int hepNUP = ev->sampleInfo.hepNUP;
        sampleWeight = reweightHEPNUPDYJets(hepNUP) * (Lumi/1000.0);
      }
      */
      //pileup weight 
      vector<double>pu = ev->sampleInfo.truepileup;
      if(pu.size() > 0) {
	  float npu = pu[0];
	  double weight = LumiWeights_.weight(npu);
      evtWeight *= weight;  
      }
    } 

    //---------------------------------------------------//
    // apply top re-weighting weight
    //---------------------------------------------------//
    double topPtWeights_offline = 1.0;
    if(!isData){
      string sampleName = ev->sampleInfo.sampleName;
      //(we have to update it ==========)
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
    TotalTopPtWeights += topPtWeights_offline;
    //cout<<"topPtWeights_offline = "<<topPtWeights_offline<<endl;
    TotalLplusJEvents++;
    evtWeight *= topPtWeights_offline; //Multiply to the total weights
    
    //---------------------------------------------------//
    //apply muon triggers
    //---------------------------------------------------//
    double nCutPass = 0.0;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
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
   
    //---------------------------------------------------//
    //get all objets e.g. leptons, jets, vertices etc.
    //---------------------------------------------------//
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
    
    //---------------------------------------------------//
    //apply selection cuts on leptons
    //---------------------------------------------------//
    if(m_init_noiso.size() > 0){
      int m_i_noiso = m_init_noiso[0];
      double mRelIso_no_iso = pfMuons_noiso[m_i_noiso].pfRelIso;
      fillHisto("pre_RelIso_mu",cutflowType, mRelIso_no_iso, evtWeight);
    }
    double pri_vtxs = Vertices.size();
    double rho_vtxs = 0.0;
    int nLepton = m_init.size();  // this condition proof that only muon + jet events
    if(nLepton != 1)continue;
    if( looseMuonVeto( m_init[0],pfMuons, isPFlow) ) continue;
    if( looseElectronVeto(-1,pfElectrons, isPFlow) ) continue;
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
    
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
    ///int muCharge = pfMuons[m_i].charge;

    ///double nCutPass = 0.0;// double nCutPass_plus = 0.0; double nCutPass_minus = 0.0;
    double mRelIso = pfMuons[m_i].pfRelIso;
    fillHisto("pt_mu", cutflowType, pfMuons[m_i].p4.pt(), evtWeight);
    fillHisto("final_RelIso_mu",cutflowType, mRelIso, evtWeight);
    fillHisto("final_multi_mu",cutflowType, count_muon, evtWeight);

    // Fill histogram after trigger and one offline isolated muon and applied 2nd lepton veto
    int nJet = j_final.size();
    fillHisto("eta_mu", cutflowType, pfMuons[m_i].p4.eta(), evtWeight);
    fillHisto("phi_mu", cutflowType, pfMuons[m_i].p4.phi(), evtWeight);
    // vertex just after one lepton selection
    //double pri_vtxs = Vertices.size();
    fillHisto("nvtx", cutflowType, pri_vtxs, evtWeight);
    fillHisto("rho_0", cutflowType, Vertices[0].rho, evtWeight);
    for( std::size_t n=0; n<Vertices.size(); n++){
      rho_vtxs = Vertices[n].rho;
      fillHisto("rho", cutflowType, rho_vtxs, evtWeight);
      }
    //---------------------------------------------------//
    // Apply Jet Selection
    //---------------------------------------------------//

    fillHisto("multi_jet", cutflowType, nJet, evtWeight);
    if(nJet < 4)continue;  // this condition implies event should contain at least 4 jets
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      fillHisto("pt_jet", cutflowType, pfJets[ind_jet].p4.pt(), evtWeight);
      fillHisto("eta_jet", cutflowType, pfJets[ind_jet].p4.eta(), evtWeight);
      fillHisto("phi_jet", cutflowType, pfJets[ind_jet].p4.phi(), evtWeight);
    }
    fillHisto("final_multi_jet", cutflowType, nJet, evtWeight);
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);
    
    //---------------------------------------------------//
    //apply MET selection   
    //---------------------------------------------------//
    double   leptonPt(0), deltaPhi(0);
    double metPt = 0; 
    metPt= met.p4.pt();
    ///if(!metuc)metPt = metWithJESJER(pfJets, &j_final, met, jes, jer);
    ///else metPt = metWithUncl(pfJets, &j_final, pfMuons, &m_init, pfElectrons, &e_final, met, metuc);
    fillHisto("pt_met", cutflowType, metPt, evtWeight);
    fillHisto("phi_met", cutflowType, met.p4.phi(), evtWeight);
    if(metPt < 20) continue;  // Missing transverse energy cut 30 GeV(CMS) for ATLAS 20 GeV 
    fillHisto("final_pt_met", cutflowType, metPt, evtWeight);
    fillHisto("final_phi_met", cutflowType, met.p4.phi(), evtWeight);
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);

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
    
    //---------------------------------------------------//
    //apply B-tagging, C-tagging
    //---------------------------------------------------//
    int count_CSVL = 0; int count_CSVM = 0; //int count_CSVT = 0; 
    vector<int> j_final_nob; j_final_nob.clear();
    double pfCISV = 0.0; //pfCombinedInclusiveSecondaryVertexV2BJetTags
    double pfCMVA = 0.0; //pfCombinedMVAV2BJetTags
    double pfCCvsL = 0.0; //pfCombinedCvsLJetTags
    double pfCCvs = 0.0; //pfCombinedCvsBJetTags
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      pfCISV = pfJets[ind_jet].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
      pfCMVA = pfJets[ind_jet].bDiscriminator["pfCombinedMVAV2BJetTags"];
      pfCCvsL= pfJets[ind_jet].bDiscriminator["pfCombinedCvsLJetTags"];
      pfCCvs = pfJets[ind_jet].bDiscriminator["pfCombinedCvsBJetTags"];
      fillHisto("pfCISV", cutflowType, pfCISV , evtWeight); 
      fillHisto("pfCMVA", cutflowType, pfCMVA , evtWeight); 
      fillHisto("pfCCvsL", cutflowType, pfCCvsL, evtWeight); 
      fillHisto("pfCCvs", cutflowType, pfCCvs , evtWeight); 
      if(pfCISV > 0.5426){
        count_CSVL++;
        fillHisto("bDiscr_Loose", cutflowType+"/BTag", pfCISV, evtWeight); 
      }
      else j_final_nob.push_back(ind_jet);  
    }
    fillHisto("CSVL_count", cutflowType, count_CSVL, evtWeight);
    
    //Get the invariant mass of csbar pair
    if(j_final_nob.size() >= 2){
      int first_index = j_final_nob[0];
      int sec_index = j_final_nob[1];
      fillHisto("pt_sjet", cutflowType, pfJets[sec_index].p4.pt(), evtWeight);
      MyLorentzVector diJet = pfJets[first_index].p4 + pfJets[sec_index].p4;
      fillHisto("mjj", cutflowType, diJet.M(), evtWeight);
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
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight); 
    nSelEvents++;  // this is the counter for counting final number of events  

    //---------------------------------------------------//
    // add set of plots after BTag:
    //---------------------------------------------------//
    fillHisto("pt_mu", cutflowType+"/BTag", pfMuons[m_i].p4.pt(), evtWeight);
    fillHisto("eta_mu", cutflowType+"/BTag", pfMuons[m_i].p4.eta(), evtWeight);
    fillHisto("phi_mu", cutflowType+"/BTag", pfMuons[m_i].p4.phi(), evtWeight);
    fillHisto("final_RelIso_mu",cutflowType+"/BTag", mRelIso, evtWeight);
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
    fillHisto("rho_0", cutflowType+"/BTag", Vertices[0].rho, evtWeight);
    for(std::size_t n=0; n<Vertices.size(); n++){
      rho_vtxs = Vertices[n].rho;
      fillHisto("rho", cutflowType+"/BTag", rho_vtxs, evtWeight);
      }
    fillHisto("wmt", cutflowType+"/BTag", mt, evtWeight);

    input_count++;
    if(input_count%10==0)
    cout << "input count: "<< input_count << endl;
    //if(input_count > 100000) break;
    //if(i > 5000) break;
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
  
  //Data, MC sample from lxplus and T2
  //CutFlowAnalysis("TTJets_ntuple_MuChannel.root", "PF", true, "TTJets_MuMC_check"); 
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/multicrab_02May17/MuData_20170502/MuRunCv1_MuData_20170502/SingleMuon/MuRunCv1_MuData_20170502/170502_125714/0000/MuRunCv1_MuData_20170502_Ntuple_2.root", "PF", true, "MuRunCv1_tot");
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/multicrab_29April17/MuMC_20170429/TTJets_MuMC_20170429/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TTJets_MuMC_20170429/170429_110042/0000/TTJets_MuMC_20170429_Ntuple_3.root", "PF", false, "TTJets_MuMC");
  
  //condor submission
  //for Data 
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", true, "outputFile");
  //for MC
  CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", false, "outputFile");
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

