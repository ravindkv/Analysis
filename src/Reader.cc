#include "interface/Reader.h"

ClassImp(Reader)

using namespace std;

unsigned int Reader::AssignEventTreeFrom(TFile *f, TString path)
{
  
  if(f==0) return 0;
  
  myTree = (TTree *) f->Get( path + "/" + "MyTree" );
  TTree *t = myTree;
  unsigned int nEntries = t->GetEntriesFast();
  if(nEntries==0) return 0;
  cout<<"============================"<<endl;
  cout<<" Total events = "<<nEntries<<endl;
  cout<<"============================"<<endl;
  myEvent = new MyEvent;
  
  myEvent = 0; t->GetBranch("MyEvent")->SetAddress(&myEvent);
  
  
  return nEntries;
}

unsigned int Reader::AssignEventTreeFromList(const char *file_list_name, TString path)
{
  int nfiles;
  char tex[200], list_name[100];
  FILE *filelist;

  // get filelist name
  if( file_list_name!=NULL && strcmp(file_list_name, "") )
    strcpy(list_name, file_list_name);
  else strcpy(list_name, "file.list");
  printf("Reading file.list = %s\n", list_name);
  filelist = fopen(list_name, "r");
  nfiles = 0;
  TChain *chain = new TChain( path + "/" + "MyTree");
  while( !feof(filelist) )
    {
      fscanf(filelist, "%s\n", tex);
      if( strcmp(tex,"//")==0 ) break;
      //printf("%s\n", tex);
      chain->Add(tex);
      nfiles++;
      //if( nfiles==nmax ) break;
    }
  printf("%d files\n", nfiles);
  
  if(!chain) return 0;
  myTree = chain;
  myEvent = new MyEvent;
  
  //myEvent = 0; chain->GetBranch("MyEvent")->SetAddress(&myEvent);
  myTree->SetBranchAddress("MyEvent", &myEvent);
  
  unsigned int nEntries = myTree->GetEntries();
  if(nEntries==0) return 0;
  cout<<" nEntries "<<nEntries<<endl;
  return nEntries;
}

MyEvent *Reader::GetNewEvent(unsigned int ientry)
{
  
  TTree *t = myTree;
  
  MyEvent *ev = myEvent;
  if(t==0 || ev==0) return 0;
  
  //check number of entries in tree
  if(ientry >= t->GetEntriesFast()) return 0;
  
  //fill the mini event
  t->GetEntry(ientry);
  return ev;
}

MyEvent *Reader::GetNewEventFromList(unsigned int ientry)
{
  
  TTree *t = myTree;
  
  MyEvent *ev = myEvent;
  if(t==0 || ev==0) return 0;
  
  //check number of entries in tree
  if(ientry >= t->GetEntries()) return 0;
  
  //fill the mini event
  t->GetEvent(ientry);
  return ev;
}

Long64_t Reader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!myTree) return -5;
   Long64_t centry = myTree->LoadTree(entry);
   if (centry < 0) return centry;
   return centry;
}

vector<MyElectron> Reader::getElectrons(MyEvent* ev, string algo)
{
  vector<MyElectron> selElectrons;
  selElectrons.clear();
  vector<MyElectron> allElectrons = ev->Electrons;
  for(size_t iele=0; iele < allElectrons.size(); ++iele){
      selElectrons.push_back(allElectrons[iele]);
  }
  return selElectrons;
}

vector<MyMuon> Reader::getMuons(MyEvent* ev, string algo)
{
  vector<MyMuon> selMuons;
  selMuons.clear();
  vector<MyMuon> allMuons = ev->Muons;
  for(size_t imu=0; imu < allMuons.size(); ++imu){
      selMuons.push_back(allMuons[imu]);
  }
  return selMuons;
}

vector<MyJet> Reader::getJets(MyEvent* ev, string algo, int jes, int jer, bool isData)
{
  vector<MyJet> selJets;
  selJets.clear();
  vector<MyJet> allJets = ev->Jets;
  for(size_t ijet=0; ijet < allJets.size(); ++ijet){
    MyJet jet = allJets[ijet];
    double jesPt = 1.0;
    //apply JES uncert scaling 
    jesPt = (1+(jet.JECUncertainty*double(jes)));
    //calculate JER sf
    double gen_pt = jet.Genp4.pt();  
    double jet_pt = jet.p4.pt();
    double sigmaJER = jet.resolution ;
    //apply JER uncert, scaling
    double delR = DeltaR(jet.Genp4, jet.p4);
    double rCone = 0.4;
    double jerP4 = 1.0;
    if(!isData && delR<rCone/2 && abs(jet_pt -gen_pt)<3*sigmaJER*jet_pt ){
      double SF = getJERSF(jet.p4.eta(), jer);
      jerP4   = max(0.0, 1.0 + (SF - 1)*(jet_pt - gen_pt)/ jet_pt);
    }
    double newPx = jerP4*jet.p4.px()*jesPt;
    double newPy = jerP4*jet.p4.py()*jesPt;
    double newPz = jerP4*jet.p4.pz();
    double newE  = jerP4*jet.p4.E();
    jet.p4.SetCoordinates(newPx, newPy, newPz, newE);
    selJets.push_back(jet);
  }
  return selJets;
}
vector<MyJet> Reader::getJetsNoCorr(MyEvent* ev, string algo)
{
  vector<MyJet> selJets;
  selJets.clear();
  vector<MyJet> allJets = ev->Jets;
  for(size_t ijet=0; ijet < allJets.size(); ++ijet){
    selJets.push_back(allJets[ijet]);
  }
  return selJets;
}

MyMET Reader::getMET(MyEvent* ev, string algo) 
{ 
  MyMET selmet;
  vector<MyMET> allMETs = ev->mets; 
  for(size_t imet=0; imet < allMETs.size(); ++imet){ 
      selmet = allMETs[imet];
  }
  return selmet; 
} 

double Reader::getJERSF(double eta, int jer){
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
  //New SF with stat + sys unc
  //https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Spring16_25nsV10a_MC/Spring16_25nsV10a_MC_SF_AK4PFchs.txt
  std::vector<double> JEREtaMap;
  std::vector<double> JERSF;
  std::vector<double> JERSFUp;
  std::vector<double> JERSFDown;
  //eta values
  JEREtaMap.push_back(0.0); 
  JEREtaMap.push_back(0.5); 
  JEREtaMap.push_back(0.8); 
  JEREtaMap.push_back(1.1); 
  JEREtaMap.push_back(1.3); 
  JEREtaMap.push_back(1.7); 
  JEREtaMap.push_back(1.9); 
  JEREtaMap.push_back(2.1); 
  JEREtaMap.push_back(2.3); 
  JEREtaMap.push_back(2.5); 
  JEREtaMap.push_back(2.8); 
  JEREtaMap.push_back(3.0); 
  JEREtaMap.push_back(3.2); 
  JEREtaMap.push_back(4.7);
  //SF values
  JERSF.push_back(1.109); JERSFDown.push_back(1.044); JERSFUp.push_back(1.174);
  JERSF.push_back(1.138); JERSFDown.push_back(1.072); JERSFUp.push_back(1.204);
  JERSF.push_back(1.114); JERSFDown.push_back(1.050); JERSFUp.push_back(1.178);
  JERSF.push_back(1.123); JERSFDown.push_back(1.022); JERSFUp.push_back(1.224);
  JERSF.push_back(1.084); JERSFDown.push_back(0.985); JERSFUp.push_back(1.183);
  JERSF.push_back(1.082); JERSFDown.push_back(0.973); JERSFUp.push_back(1.191);
  JERSF.push_back(1.140); JERSFDown.push_back(1.020); JERSFUp.push_back(1.260);
  JERSF.push_back(1.067); JERSFDown.push_back(0.953); JERSFUp.push_back(1.181);
  JERSF.push_back(1.177); JERSFDown.push_back(0.967); JERSFUp.push_back(1.387);
  JERSF.push_back(1.364); JERSFDown.push_back(1.203); JERSFUp.push_back(1.525);
  JERSF.push_back(1.857); JERSFDown.push_back(1.654); JERSFUp.push_back(2.060);
  JERSF.push_back(1.328); JERSFDown.push_back(1.203); JERSFUp.push_back(1.453);
  JERSF.push_back(1.160); JERSFDown.push_back(1.013); JERSFUp.push_back(1.307);
  double SF = 1.0;
  for(size_t i = 0; i < 13; i++){
    if(TMath::Abs(eta) >= JEREtaMap[i] && TMath::Abs(eta) < JEREtaMap[i+1]){
      if(jer == 0)SF = JERSF[i];
      else if (jer == 1) SF = JERSFUp[i];
      else if(jer == -1) SF = JERSFDown[i];  
    }
  }
  return SF;
}
double Reader::DeltaR(MyLorentzVector aV, MyLorentzVector bV){
  double deta = TMath::Abs(aV.eta() - bV.eta());
  double dphi = TMath::Abs(aV.phi() - bV.phi());
  if(dphi > M_PI) dphi = 2*M_PI - dphi;
  double delR = sqrt(deta*deta + dphi*dphi);
  return delR;
}
