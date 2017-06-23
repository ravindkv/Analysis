23c23,24
<   
---
>  
>  /* 
30d30
<   addHisto("final_multi_mu",cutflowType+"/Iso", 50,0,10);
36,44c36,86
<   addHisto("pfCISV", cutflowType+"/Iso", 50, -2., 2.);
<   addHisto("pfCMVA", cutflowType+"/Iso", 50, -2., 2.);
<   addHisto("pfCCvsL",cutflowType+"/Iso", 50, -2., 2.);
<   addHisto("pfCCvsB", cutflowType+"/Iso", 50, -2., 2.);
<   addHisto("pfCCvsL_0",cutflowType+"/Iso", 50, -2., 2.);
<   addHisto("pfCCvsL_1",cutflowType+"/Iso", 50, -2., 2.);
<   addHisto("pfCCvsB_0", cutflowType+"/Iso", 50, -2., 2.);
<   addHisto("pfCCvsB_1", cutflowType+"/Iso", 50, -2., 2.);
<   addHisto("CSVL_count", cutflowType+"/Iso", 50,0,10);
---
>   addHisto("wmt_1mu", cutflowType+"/Iso", 500, 0., 500.);
>   addHisto("wmt_1mu_4jet", cutflowType+"/Iso", 500, 0., 500.);
>   addHisto("wmt_1mu_4jet_btag", cutflowType+"/Iso", 500, 0., 500.);
>   
>   //TIGHT: base/Iso/Btag histo
>   InitHist("Iso/BTagT", cutflowType, outFile_);
>   addHisto("pfCISV", cutflowType+"/Iso/BTagT", 50, -2., 2.);
>   addHisto("pfCMVA", cutflowType+"/Iso/BTagT", 50, -2., 2.);
>   addHisto("pfCCvsL",cutflowType+"/Iso/BTagT", 50, -2., 2.);
>   addHisto("pfCCvsB", cutflowType+"/Iso/BTagT", 50, -2., 2.);
>   addHisto("CSVL_count", cutflowType+"/Iso/BTagT", 50,0,10);
>   addHisto("pfCCvsL_0",cutflowType+"/Iso/BTagT", 50, -2., 2.);
>   addHisto("pfCCvsL_1",cutflowType+"/Iso/BTagT", 50, -2., 2.);
>   addHisto("pfCCvsB_0", cutflowType+"/Iso/BTagT", 50, -2., 2.);
>   addHisto("pfCCvsB_1", cutflowType+"/Iso/BTagT", 50, -2., 2.);
>   addHisto("mjj",cutflowType+"/Iso/BTagT", 500, 0, 500);
>   addHisto("cutflow", cutflowType+"/Iso/BTagT", 50 , 0., 10.); 
>   addHisto("final_RelIso_mu",cutflowType+"/Iso/BTagT", 40,0,0.5);
>   addHisto("final_multi_jet", cutflowType+"/Iso/BTagT", 100,0,10);
>   addHisto("nvtx", cutflowType+"/Iso/BTagT", 100, 0., 100.);
>   addHisto("nvtx_6Kbins", cutflowType+"/Iso/BTagT", 6000, 0., 1000.);
>   addHisto("rhoAll", cutflowType+"/Iso/BTagT", 100, 0., 100.);
>   addHisto("chi2", cutflowType+"/Iso/BTagT", 100, 0., 500.);
>   addHisto("ndof", cutflowType+"/Iso/BTagT", 100, 0., 500.);
>   addHisto("wmt", cutflowType+"/Iso/BTagT", 500, 0., 500.);
>   addHisto("CSVT_count", cutflowType+"/Iso/BTagT", 50,0,10);
>   addHisto("pt_jetJESJER", cutflowType+"/Iso/BTagT", 50, 0., 500.);
>   
>   //MEDIUM: base/Iso/Btag histo
>   InitHist("Iso/BTagM", cutflowType, outFile_);
>   addHisto("pfCISV", cutflowType+"/Iso/BTagM", 50, -2., 2.);
>   addHisto("pfCMVA", cutflowType+"/Iso/BTagM", 50, -2., 2.);
>   addHisto("pfCCvsL",cutflowType+"/Iso/BTagM", 50, -2., 2.);
>   addHisto("pfCCvsB", cutflowType+"/Iso/BTagM", 50, -2., 2.);
>   addHisto("CSVL_count", cutflowType+"/Iso/BTagM", 50,0,10);
>   addHisto("pfCCvsL_0",cutflowType+"/Iso/BTagM", 50, -2., 2.);
>   addHisto("pfCCvsL_1",cutflowType+"/Iso/BTagM", 50, -2., 2.);
>   addHisto("pfCCvsB_0", cutflowType+"/Iso/BTagM", 50, -2., 2.);
>   addHisto("pfCCvsB_1", cutflowType+"/Iso/BTagM", 50, -2., 2.);
>   addHisto("mjj",cutflowType+"/Iso/BTagM", 500, 0, 500);
>   addHisto("cutflow", cutflowType+"/Iso/BTagM", 50 , 0., 10.); 
>   addHisto("final_RelIso_mu",cutflowType+"/Iso/BTagM", 40,0,0.5);
>   addHisto("final_multi_jet", cutflowType+"/Iso/BTagM", 100,0,10);
>   addHisto("nvtx", cutflowType+"/Iso/BTagM", 100, 0., 100.);
>   addHisto("nvtx_6Kbins", cutflowType+"/Iso/BTagM", 6000, 0., 1000.);
>   addHisto("rhoAll", cutflowType+"/Iso/BTagM", 100, 0., 100.);
>   addHisto("chi2", cutflowType+"/Iso/BTagM", 100, 0., 500.);
>   addHisto("ndof", cutflowType+"/Iso/BTagM", 100, 0., 500.);
>   addHisto("wmt", cutflowType+"/Iso/BTagM", 500, 0., 500.);
>   addHisto("CSVM_count", cutflowType+"/Iso/BTagM", 50,0,10);
>   addHisto("pt_jetJESJER", cutflowType+"/Iso/BTagM", 50, 0., 500.);
46c88
<   //base/Iso/Btag histo
---
>   //LOOSE: base/Iso/Btag histo
47a90,99
>   addHisto("pfCISV", cutflowType+"/Iso/BTag", 50, -2., 2.);
>   addHisto("pfCMVA", cutflowType+"/Iso/BTag", 50, -2., 2.);
>   addHisto("pfCCvsL",cutflowType+"/Iso/BTag", 50, -2., 2.);
>   addHisto("pfCCvsB", cutflowType+"/Iso/BTag", 50, -2., 2.);
>   addHisto("CSVL_count", cutflowType+"/Iso/BTag", 50,0,10);
>   addHisto("pfCCvsL_0",cutflowType+"/Iso/BTag", 50, -2., 2.);
>   addHisto("pfCCvsL_1",cutflowType+"/Iso/BTag", 50, -2., 2.);
>   addHisto("pfCCvsB_0", cutflowType+"/Iso/BTag", 50, -2., 2.);
>   addHisto("pfCCvsB_1", cutflowType+"/Iso/BTag", 50, -2., 2.);
>   addHisto("mjj",cutflowType+"/Iso/BTag", 500, 0, 500);
49d100
<   addHisto("final_multi_mu",cutflowType+"/Iso/BTag", 50,0,10);
56c107,108
<   addHisto("wmt", cutflowType+"/Iso/BTag", 50, 0., 200.);
---
>   addHisto("wmt", cutflowType+"/Iso/BTag", 500, 0., 500.);
>   addHisto("pt_jetJESJER", cutflowType+"/Iso/BTag", 50, 0., 500.);
61d112
<   addHisto("final_multi_mu",cutflowType+"/Iso/KinFit", 50,0,10);
78,81c129,132
<   addHisto("mjj_kfit",cutflowType+"/Iso/KinFit", 5000, 0, 500);
<   addHisto("mjj_kfit_Id",cutflowType+"/Iso/KinFit", 5000, 0, 500);
<   addHisto("mjj_kfit_Id_probfit1",cutflowType+"/Iso/KinFit", 5000, 0, 500);
<   addHisto("mjj_kfit_Id_probfit2",cutflowType+"/Iso/KinFit", 5000, 0, 500);
---
>   addHisto("mjj_kfit",cutflowType+"/Iso/KinFit", 500, 0, 500);
>   addHisto("mjj_kfit_Id",cutflowType+"/Iso/KinFit", 500, 0, 500);
>   addHisto("mjj_kfit_Id_probfit1",cutflowType+"/Iso/KinFit", 500, 0, 500);
>   addHisto("mjj_kfit_Id_probfit2",cutflowType+"/Iso/KinFit", 500, 0, 500);
88d138
<   addHisto("final_multi_mu",cutflowType+"/NonIso/BTag", 50,0,10);
102a153,154
>   addHisto("pt_metJESJER", cutflowType+"/NonIso", 50, 0., 500.);
>   addHisto("mjj",cutflowType+"/NonIso", 500, 0, 500);
107d158
<   addHisto("final_multi_mu",cutflowType+"/NonIso/BTag", 50,0,10);
114a166
>   addHisto("pt_jetJESJER", cutflowType+"/NonIso/BTag", 50, 0., 500.);
119d170
<   addHisto("final_multi_mu",cutflowType+"/NonIso/KinFit", 50,0,10);
135,139c186,190
<   addHisto("mjj_kfit",cutflowType+"/NonIso/KinFit", 5000, 0, 500);
<   addHisto("mjj_kfit_Id",cutflowType+"/NonIso/KinFit", 5000, 0, 500);
<   addHisto("mjj_kfit_Id_probfit1",cutflowType+"/NonIso/KinFit", 5000, 0, 500);
<   addHisto("mjj_kfit_Id_probfit2",cutflowType+"/NonIso/KinFit", 5000, 0, 500);
< 
---
>   addHisto("mjj_kfit",cutflowType+"/NonIso/KinFit", 500, 0, 500);
>   addHisto("mjj_kfit_Id",cutflowType+"/NonIso/KinFit", 500, 0, 500);
>   addHisto("mjj_kfit_Id_probfit1",cutflowType+"/NonIso/KinFit", 500, 0, 500);
>   addHisto("mjj_kfit_Id_probfit2",cutflowType+"/NonIso/KinFit", 500, 0, 500);
> */
214c265,266
<   histos1_[fullname] = new TH1D(name.Data(), hname.c_str(), range, min, max); 
---
>   histos1_[fullname] = new TH1F(name.Data(), hname.c_str(), range, min, max); 
>   //histos1_[fullname] = new TH1D(name.Data(), hname.c_str(), range, min, max); 
