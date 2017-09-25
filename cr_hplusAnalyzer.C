
#include "hplusAnalyzer.h"
#include <map>

using namespace std;
void hplusAnalyzer::CutFlowAnalysis(TString url, string myKey, string evtType){
  
  TString outFile("13TeV/outputDir/");
  TString Filename_ = outFile+evtType+"_Anal.root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel( 9 );
  
  TString debug_Filename_ = Filename_+"_debug.txt";
  string debug_file(debug_Filename_);
  outfile_.open(debug_file.c_str());

  CutFlowProcessor(url, myKey, "base", outFile_);
  outfile_.close();
  outFile_->Write(); 
  outFile_->Close(); 
}

//---------------------------------------------------//
//Process the cuts, event by event
//---------------------------------------------------//  
void hplusAnalyzer::CutFlowProcessor(TString url,  string myKey, TString cutflowType, TFile *outFile_){
  int input_count = 0;
  int input_count_NonIso = 0;
  bool isPFlow = (myKey.find("PFlow") != string::npos) ? true : false;
  string eAlgo("Electrons"), mAlgo("Muons"), jAlgo("Jets"), metAlgo("METs");
  
  //Uncertainty variations, JES, JER, MET unclustered, bTag
  ///int jes = 0, jer = 0, metuc = 0, bscale = 0;
  int jes = 0, jer = 0; 
  if(cutflowType.Contains("JESPlus"))jes = 1;
  else if (cutflowType.Contains("JESMinus"))jes = -1;
  else if (cutflowType.Contains("JERPlus"))jer = 1;
  else if (cutflowType.Contains("JERMinus"))jer = -1;
  //else if (cutflowType.Contains("METUCPlus"))metuc = 1;
  //else if (cutflowType.Contains("METUCMinus"))metuc = -1;
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
  double lumiTotal = 35136;
  int nEntries = evR->AssignEventTreeFrom(f);
  if(nEntries == 0) {return; }
  
  TH1F* inputcf = (TH1F*)(f->Get("allEventsFilter/totalEvents"));
  TH1F* intimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/intimepu"));
  TH1F* outoftimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/outoftimepu"));
  TH1F* totalpu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/totalpu"));
  TH1F* trueintimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/trueintimepu"));
  TH1F* trueoutoftimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/trueoutoftimepu"));
  TH1F* truetotalpu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/truetotalpu"));
  double initialEvents = inputcf->GetBinContent(1);

  cout<<"\033[01;32m input file: \033[00m"<<url<<"\n"<<endl;
  outfile_<<"input file: "<<url<<endl;
  outfile_<<"totalEvents: "<<initialEvents<<endl;
  CreateAnalHistos(cutflowType, outFile_);
  for(int nbin=1; nbin<6001; nbin++){
    fillHistoPU("intimepu", cutflowType, nbin, intimepu->GetBinContent(nbin));  
    fillHistoPU("outoftimepu", cutflowType, nbin, outoftimepu->GetBinContent(nbin));  
    fillHistoPU("totalpu", cutflowType, nbin, totalpu->GetBinContent(nbin));  
    fillHistoPU("trueintimepu", cutflowType, nbin, trueintimepu->GetBinContent(nbin));  
    fillHistoPU("trueoutoftimepu", cutflowType, nbin, trueoutoftimepu->GetBinContent(nbin));  
    fillHistoPU("truetotalpu", cutflowType, nbin, truetotalpu->GetBinContent(nbin));  
  }
  fillHisto("totalEvents", cutflowType, initialEvents, 1);
  MyEvent *ev;
  int nTriggEvent = 0, nSelEvents = 0, matchjetcount= 0, threepairjet = 0;
  double nVerticesFailCount = 0.0;
  ///double corr_jet_pair_cs = 0.0, corr_jet_pair_ud = 0.0;
  double matchedJet_q = 0.0, matchedJet_b = 0.0, matched_quark_eta_pt = 0.0, not_matchedJet_q = 0.0;
  double TotalTopPtWeights = 0, TotalLplusJEvents = 0; 
  double TotalSVEff = 0, TotalEvtsWithSV = 0.;
  
  //---------------------------------------------------//
  //BTag SF: read CSV file for SF, 2D histos for eff 
  //---------------------------------------------------//      
  //https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
  const std::string & filename 		= "stack/CSVv2_Moriond17_B_H.csv";
  const std::string & tagger 		= "CSVv2";
  const std::string & measurementType 	= "comb";
  const std::string & sysType 		= "central"; 
  const std::vector<std::string> & otherSysTypes = {"up", "down"};
  //b-quark
  BTagCalibrationReader readCSVbL= readCSVfile(filename, tagger, BTagEntry::OP_LOOSE,
    	      measurementType, sysType, otherSysTypes, BTagEntry::FLAV_B);
  //c-quark
  BTagCalibrationReader readCSVcL= readCSVfile(filename, tagger, BTagEntry::OP_LOOSE,
    	      measurementType, sysType, otherSysTypes, BTagEntry::FLAV_C);
  //other(light) quarks and gluon
  BTagCalibrationReader readCSVlL= readCSVfile(filename, tagger, BTagEntry::OP_LOOSE,
    	      "incl", sysType, otherSysTypes, BTagEntry::FLAV_UDSG);
  
  //getBTagEffHistos(f);
  TString histPath("myMiniTreeProducer/Jets/");
  TH2D* h2_BTaggingEff_Denom_b 		= (TH2D*)(f->Get(histPath+"h2_BTaggingEff_Denom_b"));
  TH2D* h2_BTaggingEff_Denom_c 		= (TH2D*)(f->Get(histPath+"h2_BTaggingEff_Denom_c"));
  TH2D* h2_BTaggingEff_Denom_udsg 	= (TH2D*)(f->Get(histPath+"h2_BTaggingEff_Denom_udsg")); 
  TH2D* h2_BTaggingEff_Num_bL 		= (TH2D*)(f->Get(histPath+"h2_BTaggingEff_Num_bL"));
  TH2D* h2_BTaggingEff_Num_cL 		= (TH2D*)(f->Get(histPath+"h2_BTaggingEff_Num_cL"));
  TH2D* h2_BTaggingEff_Num_udsgL 	= (TH2D*)(f->Get(histPath+"h2_BTaggingEff_Num_udsgL")); 
  
  //---------------------------------------------------//
  //muon scale factors from 2D histograms 
  //---------------------------------------------------//      
  //https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
  //https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2 
  //Trigger SF
  TFile *f_trigSF_BCDEF 	= new TFile("stack/muonSF/triggreSF_BCDEF.root");
  TFile *f_trigSF_GH 		= new TFile("stack/muonSF/triggreSF_GH.root");
  TH2D *h2_trigSF_BCDEF 	= (TH2D*)f_trigSF_BCDEF->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
  TH2D *h2_trigSF_GH 		= (TH2D*)f_trigSF_GH->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
  //Identification SF
  TFile *f_idSF_BCDEF 		= new TFile("stack/muonSF/idSF_BCDEF.root");
  TFile *f_idSF_GH 		= new TFile("stack/muonSF/idSF_GH.root");
  TH2D *h2_idSF_BCDEF 		= (TH2D*)f_idSF_BCDEF->Get("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
  TH2D *h2_idSF_GH 		= (TH2D*)f_idSF_GH->Get("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
  //Isolation SF
  TFile *f_isoSF_BCDEF 		= new TFile("stack/muonSF/isoSF_BCDEF.root");
  TFile *f_isoSF_GH 		= new TFile("stack/muonSF/isoSF_GH.root");
  TH2D *h2_isoSF_BCDEF 		= (TH2D*)f_isoSF_BCDEF->Get("TightISO_MediumID_pt_eta/abseta_pt_ratio");
  TH2D *h2_isoSF_GH 		= (TH2D*)f_isoSF_GH->Get("TightISO_MediumID_pt_eta/abseta_pt_ratio");
  
  //---------------------------------------------------//
  //Electron scale factors from 2D histograms 
  //---------------------------------------------------//      
  //https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Efficiencies_and_scale_factors
  //Reconstruction SF
  TFile *f_ele_recoSF 	  	= new TFile("stack/eleSF/ele_recoSF.root");
  TH2D *h2_ele_recoSF 		= (TH2D*)f_ele_recoSF->Get("EGamma_SF2D");
  //Identification SF
  TFile *f_ele_veto_idSF 	= new TFile("stack/eleSF/ele_veto_idSF.root");
  TH2D *h2_ele_veto_idSF 	= (TH2D*)f_ele_veto_idSF->Get("EGamma_SF2D");

  TFile *f_ele_loose_idSF 	= new TFile("stack/eleSF/ele_loose_idSF.root");
  TH2D *h2_ele_loose_idSF 	= (TH2D*)f_ele_loose_idSF->Get("EGamma_SF2D");
  
  TFile *f_ele_medium_idSF 	= new TFile("stack/eleSF/ele_medium_idSF.root");
  TH2D *h2_ele_medium_idSF 	= (TH2D*)f_ele_medium_idSF->Get("EGamma_SF2D");
  
  TFile *f_ele_tight_idSF 	= new TFile("stack/eleSF/ele_tight_idSF.root");
  TH2D *h2_ele_tight_idSF 	= (TH2D*)f_ele_tight_idSF->Get("EGamma_SF2D");
  
  //---------------------------------------------------//
  //get sigma of Jet Pt-resolution
  //---------------------------------------------------//
  TH1F* jetPtReso = (TH1F*)(f->Get("myMiniTreeProducer/Jets/JER_Jets"));
  double sigma_jetPtreso = jetPtReso->GetStdDev();
  cout<<"sigma_jetPtreso = "<<sigma_jetPtreso<<endl;
  
  //---------------------------------------------------//
  //loop over each event, of the ntuple
  //---------------------------------------------------//
  double kfCount = 0;
  double bTagCount = 0;
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
    if(!ev->isData){
      string sampleName = ev->sampleInfo.sampleName;
      //k-factor weight (along with lumi weight) 
      if(sampleName.find("WJetsToLNu") != string::npos || sampleName.find("W1JetsToLNu") != string::npos || sampleName.find("W2JetsToLNu") != string::npos || sampleName.find("W3JetsToLNu") != string::npos || sampleName.find("W4JetsToLNu") != string::npos){
        int hepNUP = ev->sampleInfo.hepNUP;
        double weightK = reweightHEPNUPWJets(hepNUP) * (lumiTotal/1000.0);
        if(i < 1){
          outfile_<<"Sample weight for "<<sampleName<<" is= "<< weightK<<endl;
        }
        evtWeight *= weightK;  
        fillHisto("hepNUP", cutflowType, hepNUP, evtWeight);
        fillHisto("SF_hepNUP_WJets", cutflowType, weightK, 1);
      }
      else if(sampleName.find("DYJetsToLL") != string::npos || sampleName.find("DY1JetsToLL") != string::npos || sampleName.find("DY2JetsToLL") != string::npos || sampleName.find("DY3JetsToLL") != string::npos || sampleName.find("DY4JetsToLL") != string::npos){
        int hepNUP = ev->sampleInfo.hepNUP;
        std::vector<int> hepIDUP = ev->sampleInfo.hepIDUP;
        std::vector<int> hepISTUP = ev->sampleInfo.hepISTUP;
        int countZ = 0;
        for(size_t p=0; p<hepIDUP.size(); p++){
          if(hepIDUP[p]==23 && hepISTUP[p]==2)
            countZ = countZ + 1;
        }
        if(countZ==0) hepNUP = hepNUP+1;
        double weightK = reweightHEPNUPDYJets(hepNUP) * (lumiTotal/1000.0);
        evtWeight *= weightK;  
        if(i < 1){
          outfile_<<"Sample weight for "<<sampleName<<" is= "<< weightK<<endl;
        }
        fillHisto("hepNUP", cutflowType, hepNUP, evtWeight);
        fillHisto("SF_hepNUP_DYJets", cutflowType, weightK, 1);
      }
      //lumi weight
      else {
      double sampleWeight(1.0);
      sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
      evtWeight *= sampleWeight; 
      fillHisto("SF_sampleWeight", cutflowType, sampleWeight, 1);
      if(i < 1){
        outfile_<<"Sample weight for "<<sampleName<<" is= "<< sampleWeight<<endl;
        }
      }
      //pileup weight
      vector<double>pu = ev->sampleInfo.truepileup;
      if(pu.size() > 0) {
      float npu = pu[0];
      double weightPU = LumiWeights_.weight(npu);
      evtWeight *= weightPU;  
      fillHisto("SF_weightPU", cutflowType, weightPU, 1);
      }
    } 
    //No cut
    double nCutPass = 0.0;
    fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight);
    
    //---------------------------------------------------//
    //apply muon triggers
    //---------------------------------------------------//
    double nCutPass_NonIso = 0.0;
    fillHisto("cutflow", cutflowType+"/NonIso", nCutPass_NonIso, evtWeight);
    
    bool passTrig = false;
    bool passTrig_Mu = false;
    bool passTrig_Ele = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("Ele") != string::npos) passTrig_Ele = true;
      if(trig[it].find("Mu") != string::npos) passTrig_Mu = true;
    }
    if(passTrig_Mu && passTrig_Ele) passTrig = true; 
    if(!passTrig_Mu){
    //cout << "not satisfying trigger" << endl;
      continue;
    }
    nTriggEvent++;
    nCutPass++;
    fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight);
    
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
    preSelectJets(jAlgo, &j_init, pfJets, jes, jer, sigma_jetPtreso);
    
    // clean objects //
    vector<int> e_final; e_final.clear();
    ElectronCleaning( pfElectrons, pfMuons, &e_init, &e_final, &m_init, DRMIN_ELE);
    vector<int> j_final; j_final.clear();
    vector<int> t_final; t_final.clear();
    JetCleaning(pfJets, pfMuons, pfElectrons,  &j_init, &j_final, &m_init, &e_final, DRMIN_JET);
    
    //---------------------------------------------------//
    //apply selection cuts on leptons
    //---------------------------------------------------//
    if(m_init_noiso.size() > 0){
      int m_i_noiso = m_init_noiso[0];
      double mRelIso_no_iso = pfMuons_noiso[m_i_noiso].pfRelIso;
      fillHisto("pre_RelIso_mu",cutflowType+"/Iso", mRelIso_no_iso, evtWeight);
      fillHisto("pre_RelIso_mu",cutflowType+"/NonIso", mRelIso_no_iso, evtWeight);
    }
    int nMuon = m_init_noiso.size();
    int nEle = e_final.size();
    double pri_vtxs = Vertices[0].totVtx;
    
    if(nMuon != 1)continue;
    //veto 0th muon, if other muons are stroger than the 0th.
    //we veto 0th if it is a fake muon
    //check if 0th muon has mediumMuon ID
    int m_i = m_init_noiso[0];
    if(!isMediumMuon(&pfMuons_noiso[m_i], isPFlow)) continue;
    if(looseMuonVeto(m_i, pfMuons_noiso, isPFlow) ) continue;
    nCutPass++; 
    fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight); // one lepton 
    
    if(nEle != 2)continue;
    int e_i = e_final[0];
    if(looseElectronVetoTemp(e_final[0], e_final[1], pfElectrons, Vertices[0], isPFlow)) continue;
    nCutPass++; 
    fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight); // one lepton 
    //events should not have only one electron, in muon+jets channel
    //Transverse mass b/w lepton and MET
    double metPt = 0; 
    double   leptonPt(0), deltaPhi(0);
    metPt = metWithJESJER(pfJets, &j_final, met, jes, jer, sigma_jetPtreso);
    leptonPt = TMath::Abs(pfMuons_noiso[m_i].p4.pt());
    deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(pfMuons_noiso[m_i].p4, met.p4);
    double mt = sqrt (  2*leptonPt*metPt*(1 - cos(deltaPhi) ) ) ;
    int count_muon = m_init.size();
    ///int muCharge = pfMuons[m_i].charge;
    
    //---------------------------------------------------//
    //apply muon SF to eventWeights 
    //---------------------------------------------------//
    double lumi_BCDEF = 19.04; double lumi_GH = 16.09 ;	
    double lumi = lumi_BCDEF + lumi_GH;
    //trigger 	
    double muSFtrig_BCDEF 	= getMuonTrigSF(h2_trigSF_BCDEF, pfMuons[m_i].p4.eta(), pfMuons[m_i].p4.pt());
    double muSFtrig_GH 	= getMuonTrigSF(h2_trigSF_GH, pfMuons[m_i].p4.eta(), pfMuons[m_i].p4.pt());
    double muSFtrig 	= (muSFtrig_BCDEF*lumi_BCDEF + muSFtrig_GH*lumi_GH)/lumi; 
    //identification: Please note, the 2D histo has Pt<120. So for Pt>120, we have taken muSFid = 1.0;
    double muSFid_BCDEF 	= getMuonSF(h2_idSF_BCDEF, pfMuons[m_i].p4.eta(), pfMuons[m_i].p4.pt());
    double muSFid_GH 	= getMuonSF(h2_idSF_GH, pfMuons[m_i].p4.eta(), pfMuons[m_i].p4.pt());
    double muSFid 		= (muSFid_BCDEF*lumi_BCDEF + muSFid_GH*lumi_GH)/lumi; 
    //isolation: Please note, the 2D histo has Pt<120. So for Pt>120, we have taken muSFiso = 1.0;
    double muSFiso_BCDEF 	= getMuonSF(h2_isoSF_BCDEF, pfMuons[m_i].p4.eta(), pfMuons[m_i].p4.pt());
    double muSFiso_GH 	= getMuonSF(h2_isoSF_GH, pfMuons[m_i].p4.eta(), pfMuons[m_i].p4.pt());
    double muSFiso 		= (muSFiso_BCDEF*lumi_BCDEF + muSFiso_GH*lumi_GH)/lumi; 
    //combined SF
    double muSF =1.0;
    if(!ev->isData){
      muSF = muSFtrig*muSFid*muSFiso;	
    }
    fillHisto("SF_muonSF", cutflowType, muSF, 1);
    evtWeight *= muSF;
    nCutPass++; 
    fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight); // one lepton 
    
    //---------------------------------------------------//
    //apply Electron SF to eventWeights 
    //---------------------------------------------------//
    //Reco	
    double ele_recoSF 		= getEleSF(h2_ele_recoSF, pfElectrons[e_i].eleSCEta, pfElectrons[e_i].p4.pt());
    //ID
    //double ele_veto_idSF  	= getEleSF(h2_ele_veto_idSF, pfElectrons[e_i].eleSCEta, pfElectrons[e_i].p4.pt());
    //double ele_loose_idSF  	= getEleSF(h2_ele_loose_idSF, pfElectrons[e_i].eleSCEta, pfElectrons[e_i].p4.pt());
    double ele_medium_idSF  	= getEleSF(h2_ele_medium_idSF, pfElectrons[e_i].eleSCEta, pfElectrons[e_i].p4.pt());
    //double ele_tight_idSF  	= getEleSF(h2_ele_tight_idSF, pfElectrons[e_i].eleSCEta, pfElectrons[e_i].p4.pt());
    double eleSF =1.0;
    if(!ev->isData){
      eleSF = ele_recoSF*ele_medium_idSF;	
    }
    fillHisto("SF_muonSF", cutflowType, muSF, 1);
    evtWeight *= eleSF;
    nCutPass++; 
    fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight); // one lepton 
    //---------------------------------------------------//
    // Iso(<0.12) and Non-iso(>0.12) region 
    //---------------------------------------------------//
    bool noisofound = false;
    bool isofound = false;
    double tmp_iso = pfMuons_noiso[m_i].pfRelIso;
    if(tmp_iso < 0.15 ) isofound = true;
    fillHisto("RelIso_mu",cutflowType, tmp_iso, evtWeight);
    if(tmp_iso > 0.15 && tmp_iso < 0.40) noisofound = true;
  
    ////////////////////////////////////////////////////////////////  
    //  Isolation region
    ////////////////////////////////////////////////////////////////  
    if(isofound){
      nCutPass++;
      fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight);
      ///double nCutPass = 0.0;// double nCutPass_plus = 0.0; double nCutPass_minus = 0.0;
      double mRelIso = pfMuons[m_i].pfRelIso;
      fillHisto("pt_mu", cutflowType+"/Iso", pfMuons[m_i].p4.pt(), evtWeight);
      fillHisto("final_RelIso_mu",cutflowType+"/Iso", mRelIso, evtWeight);

      // Fill histogram after trigger and one offline isolated muon and applied 2nd lepton veto
      fillHisto("eta_mu", cutflowType+"/Iso", pfMuons[m_i].p4.eta(), evtWeight);
      fillHisto("phi_mu", cutflowType+"/Iso", pfMuons[m_i].p4.phi(), evtWeight);
      // vertex just after one lepton selection
      //double pri_vtxs = Vertices.size();
      fillHisto("nvtx", cutflowType+"/Iso", pri_vtxs, evtWeight);
      fillHisto("nvtx_6Kbins", cutflowType+"/Iso", pri_vtxs, evtWeight);
      for( std::size_t n=0; n<Vertices.size(); n++){
        fillHisto("rhoAll", cutflowType+"/Iso", Vertices[n].rhoAll, evtWeight);
        fillHisto("chi2", cutflowType+"/Iso", Vertices[n].chi2, evtWeight);
        fillHisto("ndof", cutflowType+"/Iso", Vertices[n].ndof, evtWeight);
        }
      
      //---------------------------------------------------//
      // Apply Jet Selection
      //---------------------------------------------------//
      int count_jets = j_final.size();
      fillHisto("multi_jet", cutflowType+"/Iso", count_jets, evtWeight);
      //if(count_jets < 4)continue;  // events should have 4 or more jets
      if(count_jets < 2)continue;  // events should have 2 or more jets
      nCutPass++;
      fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight);
      for(size_t ijet = 0; ijet < j_final.size(); ijet++){
        int ind_jet = j_final[ijet];
        double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer, sigma_jetPtreso);
	fillHisto("pt_jet", cutflowType+"/Iso", jetPt, evtWeight);
        fillHisto("eta_jet", cutflowType+"/Iso", pfJets[ind_jet].p4.eta(), evtWeight);
        fillHisto("phi_jet", cutflowType+"/Iso", pfJets[ind_jet].p4.phi(), evtWeight);
      }
      fillHisto("final_multi_jet", cutflowType+"/Iso", count_jets, evtWeight);
      
      //---------------------------------------------------//
      //apply MET selection   
      //---------------------------------------------------//
      ///if(!metuc)metPt = metWithJESJER(pfJets, &j_final, met, jes, jer);
      ///else metPt = metWithUncl(pfJets, &j_final, pfMuons, &m_init, pfElectrons, &e_final, met, metuc);
      metPt = metWithJESJER(pfJets, &j_final, met, jes, jer, sigma_jetPtreso);
      fillHisto("pt_met", cutflowType+"/Iso", metPt, evtWeight);
      fillHisto("phi_met", cutflowType+"/Iso", met.p4.phi(), evtWeight);
      if(metPt < 20) continue;  // Missing transverse energy cut 30 GeV(CMS) for ATLAS 20 GeV 
      nCutPass++;
      fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight);
      fillHisto("final_pt_met", cutflowType+"/Iso", metPt, evtWeight);
      fillHisto("final_phi_met", cutflowType+"/Iso", met.p4.phi(), evtWeight);
      
      //---------------------------------------------------//
      //apply B-tagging, C-tagging
      //---------------------------------------------------//
      vector<int> j_final_nob; j_final_nob.clear();
      vector<int> j_final_b; j_final_b.clear();
      vector<double> bdiscr; bdiscr.clear();
      double pfCISV = 0.0; //pfCombinedInclusiveSecondaryVertexV2BJetTags
      double pfCMVA = 0.0; //pfCombinedMVAV2BJetTags
      double pfCCvsL = 0.0;//pfCombinedCvsLJetTags
      double pfCCvsB = 0.0; //pfCombinedCvsBJetTags
      double pfCCvsL_0 = 0.0;
      double pfCCvsL_1 = 0.0;
      double pfCCvsB_0 = 0.0;
      double pfCCvsB_1 = 0.0;
      int bscale = 0;

      //LOOSE BTAG
      int count_CSVL_SF = 0; 
      int count_CSVL = 0;
      ///bool isBtagL = true; bool isBtagM = true; isBtagT = true; 
      bool isBtag = false;
      for(size_t ijet = 0; ijet < j_final.size(); ijet++){
        int ind_jet = j_final[ijet];
        pfCISV = pfJets[ind_jet].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
        pfCMVA = pfJets[ind_jet].bDiscriminator["pfCombinedMVAV2BJetTags"];
        pfCCvsL= pfJets[ind_jet].bDiscriminator["pfCombinedCvsLJetTags"];
        pfCCvsB = pfJets[ind_jet].bDiscriminator["pfCombinedCvsBJetTags"];
        fillHisto("pfCISV", cutflowType+"/Iso/BTag", pfCISV , evtWeight); 
        fillHisto("pfCMVA", cutflowType+"/Iso/BTag", pfCMVA , evtWeight); 
        fillHisto("pfCCvsL", cutflowType+"/Iso/BTag", pfCCvsL, evtWeight); 
        fillHisto("pfCCvsB", cutflowType+"/Iso/BTag", pfCCvsB , evtWeight); 
        //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
        //b-quark
        if(abs(pfJets[ind_jet].partonFlavour) ==5)
          isBtag = getBtagWithSF(readCSVbL, h2_BTaggingEff_Num_bL, h2_BTaggingEff_Denom_b, pfJets[ind_jet], ev->isData, bscale); 
        //c-quark
        else if(abs(pfJets[ind_jet].partonFlavour) ==4) 
          isBtag = getBtagWithSF(readCSVcL, h2_BTaggingEff_Num_cL, h2_BTaggingEff_Denom_c, pfJets[ind_jet], ev->isData, bscale); 
        //other quarks and gluon
        else isBtag = getBtagWithSF(readCSVlL, h2_BTaggingEff_Num_udsgL, h2_BTaggingEff_Denom_udsg, pfJets[ind_jet], ev->isData, bscale); 
        
        if(isBtag){
          count_CSVL_SF++; 
          double jetPt = jetPtWithJESJER(pfJets[ijet], 0, 0, sigma_jetPtreso);
          fillHisto("pt_bjet", cutflowType+"/Iso"+"/BTag", jetPt, evtWeight);
          fillHisto("eta_bjet", cutflowType+"/Iso"+"/BTag", pfJets[ijet].p4.eta(), evtWeight);
          j_final_b.push_back(ind_jet);
          bdiscr.push_back(pfCISV);
        }
        else j_final_nob.push_back(ind_jet); 
        if(pfCISV >0.5426) count_CSVL ++;
      }
      if(count_CSVL ==2 ){
      nCutPass = 9;
      fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight); 
      }
      fillHisto("CSVL_count", cutflowType+"/Iso/BTag", count_CSVL_SF, evtWeight);
      if(count_CSVL_SF != 2) continue; // Demanding for 2L b-tagged jets
      //if(count_CSVL_SF <= 1) continue; // Demanding for 2L b-tagged jets
      nCutPass = 10;
      fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight); 
      
      ///if(mt < 30) continue; // mt cuts
      nCutPass = 11;
      fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight); 
      //---------------------------------------------------//
      // apply top re-weighting weight
      //---------------------------------------------------//
      double topPtWeights_offline = 1.0;
      if(!ev->isData){
        string sampleName = ev->sampleInfo.sampleName;
        if(sampleName.find("TTJets") != string::npos || sampleName.find("HplusM120") != string::npos)
	{
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
      fillHisto("SF_topPtWeights", cutflowType, topPtWeights_offline, 1);
      //cout<<"topPtWeights_offline = "<<topPtWeights_offline<<endl;
      TotalLplusJEvents++;
      evtWeight *= topPtWeights_offline; //Multiply to the total weights
      nCutPass = 12;
      fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight); 
      
      //---------------------------------------------------//
      //invariant mass of c sbar
      //---------------------------------------------------//

      //sort j_final_b w.r.t b-discriminator value(ascending order)
      std::map<double, int> bdiscr_sorted_bjets;
      for(unsigned long k=0; k<j_final_b.size(); k++){
        bdiscr_sorted_bjets.insert(pair <double, int> (bdiscr[k],j_final_b[k])); 
      }
      map <double, int> :: iterator bdiscr_itr;
      int index_of_2nd_bjet;
      int index_of_1st_bjet;
      vector<int> index_of_other_bjets;
      int total_bjets = j_final_b.size();
      for(bdiscr_itr = bdiscr_sorted_bjets.begin(); bdiscr_itr != bdiscr_sorted_bjets.end(); ++bdiscr_itr){
	 total_bjets --;
	 if(total_bjets==1) index_of_2nd_bjet = bdiscr_itr->second;  
	 else if(total_bjets==0) index_of_1st_bjet = bdiscr_itr->second;  
	 else index_of_other_bjets.push_back(bdiscr_itr->second);
      }
      //mjj will involve 2 non-bjet, highest pt jets
      if(j_final_b.size()==2){
        if(j_final_nob.size() >= 2){
          int index_of_1st_mjj = j_final_nob[0];
          int index_of_2nd_mjj = j_final_nob[1];
          MyLorentzVector diJet = pfJets[index_of_1st_mjj].p4 + pfJets[index_of_2nd_mjj].p4;
          fillHisto("mjj", cutflowType+"/Iso/BTag", diJet.M(), evtWeight);
        }
      }
      //Arrange other bjets and non-bjets in pt order in a list
      //mjj will involve 2 highest pt jets in this list
      else{ 
        std::map<double, int> pt_sorted_jets;
	for(unsigned long k=0; k<j_final_nob.size(); k++){
          double pt_of_nobjet = pfJets[j_final_nob[k]].p4.pt(); 
          pt_sorted_jets.insert(pair <double, int> (pt_of_nobjet, j_final_nob[k])); 
	}
        for(unsigned long k=0; k<index_of_other_bjets.size(); k++){
          double pt_of_other = pfJets[index_of_other_bjets[k]].p4.pt();
          pt_sorted_jets.insert(pair <double, int> (pt_of_other, index_of_other_bjets[k])); 
        }
        //select two highest pt jets
        int index_of_1st_mjj = 0;
        int index_of_2nd_mjj = 0;
        int total_jets_for_mjj = pt_sorted_jets.size();
        map <double, int> :: iterator itr_pt;
        for(itr_pt = pt_sorted_jets.begin(); itr_pt != pt_sorted_jets.end(); ++itr_pt){
           total_jets_for_mjj --;
           if(total_jets_for_mjj==1) index_of_2nd_mjj = itr_pt->second;  
           if(total_jets_for_mjj==0) index_of_1st_mjj = itr_pt->second;  
        }
	MyLorentzVector diJet = pfJets[index_of_1st_mjj].p4 + pfJets[index_of_2nd_mjj].p4;
        fillHisto("mjj", cutflowType+"/Iso/BTag", diJet.M(), evtWeight);
      }	

      //---------------------------------------------------//
      // add set of plots after BTag:
      //---------------------------------------------------//
      fillHisto("pt_mu", cutflowType+"/Iso"+"/BTag", pfMuons[m_i].p4.pt(), evtWeight);
      fillHisto("eta_mu", cutflowType+"/Iso"+"/BTag", pfMuons[m_i].p4.eta(), evtWeight);
      fillHisto("phi_mu", cutflowType+"/Iso"+"/BTag", pfMuons[m_i].p4.phi(), evtWeight);
      fillHisto("final_RelIso_mu",cutflowType+"/Iso"+"/BTag", mRelIso, evtWeight);
      for(size_t ijet = 0; ijet < j_final.size(); ijet++){
        int ind_jet = j_final[ijet];
        double jetPt = jetPtWithJESJER(pfJets[ind_jet], 0, 0, sigma_jetPtreso);
        fillHisto("pt_jet", cutflowType+"/Iso"+"/BTag", jetPt, evtWeight);
        fillHisto("eta_jet", cutflowType+"/Iso"+"/BTag", pfJets[ind_jet].p4.eta(), evtWeight);
        fillHisto("phi_jet", cutflowType+"/Iso"+"/BTag", pfJets[ind_jet].p4.phi(), evtWeight);
      }
      fillHisto("final_multi_jet", cutflowType+"/Iso"+"/BTag", count_jets, evtWeight);
      fillHisto("final_pt_met", cutflowType+"/Iso"+"/BTag", metPt, evtWeight);
      fillHisto("nvtx", cutflowType+"/Iso"+"/BTag", pri_vtxs, evtWeight);
      fillHisto("nvtx_6Kbins", cutflowType+"/Iso"+"/BTag", pri_vtxs, evtWeight);
      for(std::size_t n=0; n<Vertices.size(); n++){
        fillHisto("rhoAll", cutflowType+"/Iso"+"/BTag", Vertices[n].rhoAll, evtWeight);
        fillHisto("chi2", cutflowType+"/Iso"+"/BTag", Vertices[n].chi2, evtWeight);
        fillHisto("ndof", cutflowType+"/Iso"+"/BTag", Vertices[n].ndof, evtWeight);
        }
      fillHisto("wmt", cutflowType+"/Iso"+"/BTag", mt, evtWeight);
      
      bTagCount++;
      input_count++;
      if(input_count%10==0)
      cout << "input count iso: "<< input_count << endl;
      //if(i > 30000) break;
    }//isomuon
  }//event loop
  
  cout<<"bTagCount = "<<bTagCount<<endl;
  cout<<"kfCount = "<<kfCount<<endl;
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
  //CutFlowAnalysis("outFile_.root", "PF", "TTJets_MuMC_check"); 
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/", "PF", "");
  
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/ntuple_MuMC_kfitL_20170915/MuMC_20170915/TTJets_MuMC_20170915/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TTJets_MuMC_20170915/170914_234059/0000/TTJets_MuMC_20170915_Ntuple_1.root", "PF", "");
  
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/ntuple_MuData_kfitL_20170915/MuData_20170915/MuRunB2v2_MuData_20170915/SingleMuon/MuRunB2v2_MuData_20170915/170914_234448/0000/MuRunB2v2_MuData_20170915_Ntuple_1.root", "PF", "");

  //====================================
  //condor submission
  CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", "outputFile");
  //====================================
} 
