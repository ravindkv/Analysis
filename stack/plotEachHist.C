
void plotHist(TFile *file, TString chName, TString catName){
    TH1F * hData = (TH1F*)file->Get("data_obs");
    TH1F * hTT = (TH1F*)file->Get("ttbar");
    TH1F * hST = (TH1F*)file->Get("stop");
    TH1F * hSig = (TH1F*)file->Get("WH100");
    TH1F * hQCD = (TH1F*)file->Get("qcd");
    TH1F * hWJet = (TH1F*)file->Get("wjet");
    TH1F * hDYJet = (TH1F*)file->Get("zjet");
    TH1F * hVV = (TH1F*)file->Get("vv");

    TCanvas *can1 = new TCanvas();
    can1->Divide(2,2);
    can1->cd(1);
    hData->Draw("ALP");
    can1->cd(2);
    hTT->Draw("ALP");
    can1->cd(3);
    hST->Draw("ALP");
    can1->cd(4);
    hSig->Draw("ALP");
    can1->SaveAs(chName+"_"+catName+"_part1.pdf");

    TCanvas *can2 = new TCanvas();
    can2->Divide(2,2);
    can2->cd(1);
    hQCD->Draw("ALP");
    can2->cd(2);
    hDYJet->Draw("ALP");
    can2->cd(3);
    hWJet->Draw("ALP");
    can2->cd(4);
    hVV->Draw("ALP");
    can2->SaveAs(chName+"_"+catName+"_part2.pdf");
}

void plotEachHist(){
TFile * fEleL = new TFile("Shapes_hcs_13TeV_ele_KinFit_mjj_kfit_CTagExL_WH100.root");
TFile * fEleM = new TFile("Shapes_hcs_13TeV_ele_KinFit_mjj_kfit_CTagExM_WH100.root");
TFile * fEleT = new TFile("Shapes_hcs_13TeV_ele_KinFit_mjj_kfit_CTagExT_WH100.root");
TFile * fMuL  = new TFile("Shapes_hcs_13TeV_mu_KinFit_mjj_kfit_CTagExL_WH100.root");
TFile * fMuM  = new TFile("Shapes_hcs_13TeV_mu_KinFit_mjj_kfit_CTagExM_WH100.root");
TFile * fMuT  = new TFile("Shapes_hcs_13TeV_mu_KinFit_mjj_kfit_CTagExT_WH100.root");

plotHist(fEleL, "electron", "loose");
//plotHist(fEleM, "electron", "medium");
//plotHist(fEleT, "electron", "tight");
//p
//plotHist(fMuL, "muon", "loose");
//plotHist(fMuM, "muon", "medium");
//plotHist(fMuT, "muon", "tight");
}
