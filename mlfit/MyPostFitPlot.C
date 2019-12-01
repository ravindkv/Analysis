#include "MyPostFitPlot.h"

using namespace std;

//SAVE HISTOS ON DISK
bool isSaveHisto = true;
bool isData = true;
int iPeriod = 4;
int iPosX = 10;
double xMin = 1;
double xMax = 30;
double yMin = 0.8;
double yMax = 1.2;


//--------------------------------------------//
//stack histos
//--------------------------------------------//
void example_stack(TString chName, TString chLable, TString charmCat){
  //Pad
  setTDRStyle();
  TCanvas *c1 = myCanvas(chName, iPeriod, 1, isData);
  const float xpad[2] = {0.,1.0};
  const float ypad[4] = {0.,0.30,0.30,1.0};
  c1->Divide(1, 2); c1->cd(1);
  gPad->SetRightMargin(0.03);
  gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);

  //hMC = all Bkg MC samples
  TH1F* hTT   = getHisto(fitDiagOut, chName, "ttbar");
  TH1F* hST   = getHisto(fitDiagOut, chName, "stop");
  TH1F* hQCD  = getHisto(fitDiagOut, chName, "qcd");
  TH1F* hWjet = getHisto(fitDiagOut, chName, "wjet");
  TH1F* hDY   = getHisto(fitDiagOut, chName, "zjet");
  TH1F* hVV = getHisto(fitDiagOut, chName, "vv");
  TH1F* hMC   = getHisto(fitDiagOut, chName, "total_background");
  //Define stacked histo
  THStack* hStack = new THStack("hStack","");
  //VV is the base histo
  hVV->SetFillColor(kGreen +1);
  hStack->Add(hVV);
  stackHisto(hDY , hStack,  kOrange +1);
  stackHisto(hWjet , hStack,  kYellow +1);
  stackHisto(hQCD  , hStack,  kRed +1);
  stackHisto(hST , hStack,  kViolet +1);
  stackHisto(hTT, hStack,  kCyan + 1);

  gPad->SetTopMargin(0.09);
  gPad->SetBottomMargin(0.0);
  hStack->Draw("HIST");
  hStack->SetMinimum(1.0);
  //hStack ->GetXaxis()->SetRangeUser(xMin, xMax);
  TGaxis::SetMaxDigits(3);
  hStack->GetYaxis()->SetTitle("Events / 5 GeV");
  hStack->GetXaxis()->SetTitle("xTitle");

  //-------------------------------------///
  //unc band
  //-------------------------------------///
  TGraphAsymmErrors *UncBand;
  //UncBand = UNCGRAPH(addHistoForUnc(fitDiagOut, chName), true, false);
  UncBand = UNCGRAPH(hMC, true, false);
  UncBand->SetFillColor(kSpring +9);
  UncBand->SetFillStyle(3008);
  UncBand->Draw(" E2 same");

  //Data
  TGraphAsymmErrors *dataGraph = (TGraphAsymmErrors*)fitDiagOut->Get("shapes_fit_s/"+chName+"/data");
  TH1F *hData = convertGraphToHist(dataGraph);
  hData->SetFillColor(kBlack);
  hData->SetMarkerStyle(20); 
  hData->SetMarkerSize(1.2);
  c1->Update();
  if(isData)hData->Draw("Esame");

  //-------------------------------
  //Signal
  //-------------------------------
  TH1F* hSig = getHistoSig(fitDiagOut, chName, "WH100");
  hSig->SetLineColor(kRed); hSig->SetLineStyle(2);
  hSig->SetLineWidth(3); hSig->SetFillColor(0);
  hSig->Draw("HISTSAME");
  hStack->SetMaximum(1.1*hStack->GetMaximum());

  //-------------------------------------///
  //  Draw Pave Text
  //-------------------------------------///
  //channel
  TPaveText *ch = paveText(0.20,0.5061504,0.35,0.60, 0, 19, 1, 0, 132);
  ch->SetTextSize(0.05);
  if(chName.Contains("ch1") || chName.Contains("ch2") ||chName.Contains("ch3") ) 
      ch->AddText("#splitline{#splitline{#mu + jets}{After KF}}{"+charmCat+"}");
  if(chName.Contains("ch4") || chName.Contains("ch5") ||chName.Contains("ch6") ) 
      ch->AddText("#splitline{#splitline{e + jets}{After KF}}{"+charmCat+"}");
  gPad->RedrawAxis();
  c1->Update();
  

  //-------------------------------
  //Legends
  //-------------------------------
  TLegend* leg = new TLegend(0.6818792,0.3061504,0.8512081,0.908861,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(kBlack);
  leg->SetTextFont(42);
  leg->SetTextAngle(0);
  leg->SetTextSize(0.04);
  leg->SetTextAlign(12);
  if(isData)leg->AddEntry(hData,"Data","PE");
  leg->AddEntry(hTT,"t#bar{t} + jets","F");
  leg->AddEntry(hST,"Single t","F");
  leg->AddEntry(hQCD,"QCD","F");
  leg->AddEntry(hWjet,"W + jets","F");
  leg->AddEntry(hDY,"Z/#gamma + jets","F");
  leg->AddEntry(hVV,"VV","F");
  leg->AddEntry(UncBand, "Post-fit unc.","F");
  leg->AddEntry(hSig, "m_{H^{+}} = 100 GeV","L");
  leg->AddEntry((TObject*)0, "B(t #rightarrow H^{+} b) = 0.065","");
  c1->Update();
  ch->Draw();
  leg->Draw();

  //-------------------------------------///
  // Ratio = DATA/Bkg
  //-------------------------------------///
  c1->cd(2);
  c1->Update();
  c1->RedrawAxis();
  c1->GetFrame()->Draw();
  gPad->SetTopMargin(0); 
  gPad->SetBottomMargin(0.30); //gPad->SetGridy();
  gPad->SetRightMargin(0.03);
  gPad->SetTickx(0);
  gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
  TH1F *hRatio = (TH1F*)hData->Clone("hRatio");
  hRatio->Reset();
  hRatio->Add(hData);
  hRatio->Divide(hMC);
  hRatio->SetFillColor(kBlack);
  hRatio->GetYaxis()->SetRangeUser(yMin, yMax);
  hRatio->GetXaxis()->SetTitleSize(0.13);
  hRatio->GetXaxis()->SetLabelSize(0.15);
  hRatio->GetXaxis()->SetLabelFont(42);
  hRatio->GetXaxis()->SetLabelColor(kBlack);
  hRatio->GetXaxis()->SetAxisColor(kBlack);
  hRatio->GetXaxis()->SetTickLength(0.06);
  //hRatio->GetXaxis()->SetRangeUser(xMin, xMax);
  hRatio->GetXaxis()->SetTitleOffset(1);
  hRatio->GetXaxis()->SetLabelOffset(0.01);
  hRatio->SetMarkerStyle(20); 
  hRatio->SetMarkerSize(1.2);
  hRatio->GetYaxis()->SetTitleSize(0.13);
  hRatio->GetYaxis()->SetLabelSize(0.11);
  hRatio->GetYaxis()->SetLabelFont(42);
  hRatio->GetYaxis()->SetAxisColor(1);
  hRatio->GetYaxis()->SetTickLength(0.03);
  hRatio->GetYaxis()->SetNdivisions(6,5,0);
  hRatio->GetYaxis()->SetTitleOffset(0.5);
  hRatio->GetYaxis()->SetLabelOffset(0.01);
  hRatio->GetYaxis()->CenterTitle();
  hRatio->GetXaxis()->SetTitle("m_{jj} (GeV)");
  for(i = 0; i<30; i++){
    string lable = "    "+intToStr(i*5 + 20);
    if(i%3==!0)hRatio->GetXaxis()->SetBinLabel(i, lable.c_str());
    hRatio->GetXaxis()->LabelsOption("h");
    //hRatio->GetXaxis()->SetTickLength(0.08);
    hRatio->GetXaxis()->SetLabelOffset(0.02);
  }
  hRatio->GetYaxis()->SetTitle("#frac{Data}{Bkg}");
  
  //unc band
  hRatio->Draw("E"); // use "P" or "AP"
  TGraphAsymmErrors *UncBand_Ratio;
  UncBand_Ratio = UNCGRAPH(hMC, false, true);
  UncBand_Ratio->SetFillColor(kSpring +9);
  //UncBand_Ratio->SetFillStyle(3001);
  UncBand_Ratio->Draw("E2 same");
  hRatio->Draw("E same"); // use "P" or "AP"
  //base line at 1
  TF1 *baseLine = new TF1("baseLine","1", -100, 2000); 
  baseLine->SetLineColor(kBlack);
  baseLine->Draw("SAME");
  CMS_lumi(c1, iPeriod, iPosX);
  c1->Update();
  if(isSaveHisto){
    TString outFile("$PWD/output/mjj_postfit_"+chName+".pdf");
    c1->SaveAs(outFile);
    c1->Close();
  }
}


void MyPostFitPlot(){
  example_stack("ch1", "#mu + jets", "Loose");
  example_stack("ch2", "#mu + jets", "Medium");
  example_stack("ch3", "#mu + jets", "Tight");
  example_stack("ch4", "e + jets", "Loose");
  example_stack("ch5", "e + jets", "Medium");
  example_stack("ch6", "e + jets", "Tight");
} 

