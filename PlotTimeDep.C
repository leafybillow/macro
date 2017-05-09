void PlotTimeDep(int runnum){
  gROOT->SetBatch(kTRUE);
  TFile *rootfile = TFile::Open(Form("./ROOTfiles/parity16_%d_standard.root",runnum));
      
  int comb_first[6]={1,2,3,1,2,1};
  int comb_second[6]={2,3,4,3,4,4};
  for(int isam =0; isam < 8 ; isam++){
    for(int icomb =0 ; icomb <6 ;icomb++){
      R->SetAlias(Form("Ab%d_%d",isam+1,icomb),Form("1e6*(vqwk5_%d_b%d_cal-vqwk5_%d_b%d_cal)/(vqwk5_%d_b%d_cal+vqwk5_%d_b%d_cal)",isam,comb_first[icomb],isam,comb_second[icomb],isam,comb_first[icomb],isam,comb_second[icomb]));
    }
    R->SetAlias(Form("Ab%d_%d",isam+1,6),Form("1e6*(vqwk5_%d_b1_cal +vqwk5_%d_b2_cal - vqwk5_%d_b3_cal - vqwk5_%d_b4_cal)/(vqwk5_%d_b1_cal + vqwk5_%d_b2_cal+ vqwk5_%d_b3_cal + vqwk5_%d_b4_cal)",isam,isam,isam,isam,isam,isam,isam,isam));
    P->SetAlias(Form("Ab%d_%d",isam+1,7),Form("asym_blumi%d",isam+1));
  }
  for(int ibcm=0; ibcm<4; ibcm++){  // ibcm = vqwk index
    // Aq_<n> : <n> bcm number
    for(int icomb =0; icomb<6;icomb++){
      R->SetAlias(Form("Aq%d_%d",ibcm+1,icomb),Form("1e6*(vqwk0_%d_b%d_cal-vqwk0_%d_b%d_cal)/(vqwk0_%d_b%d_cal+vqwk0_%d_b%d_cal)",ibcm,comb_first[icomb],ibcm,comb_second[icomb],ibcm,comb_first[icomb],ibcm,comb_second[icomb]));
    }
    R->SetAlias(Form("Aq%d_%d",ibcm+1,6),Form("1e6*(vqwk0_%d_b1_cal +vqwk0_%d_b2_cal - vqwk0_%d_b3_cal - vqwk0_%d_b4_cal)/(vqwk0_%d_b1_cal + vqwk0_%d_b2_cal+ vqwk0_%d_b3_cal + vqwk0_%d_b4_cal)",ibcm,ibcm,ibcm,ibcm,ibcm,ibcm,ibcm,ibcm));
    P->SetAlias(Form("Aq%d_%d",ibcm+1,7),Form("asym_bcm%d",ibcm+1));
  }

  TString histo_title[8]={"(1-2)/(1+2)","(2-3)/(2+3)","(3-4)/(3+4)","(1-3)/(1+3)","(2-4)/(2+4)","(1-4)/(1+4)","(1+2-3-4)/(1+2+3+4)","Asym"};
  
  TCanvas *c1 = new TCanvas("c1","c1",1080,1920);
  c1->Divide(1,8);

  TCanvas *c_bcm = new TCanvas("c_bcm","c_bcm",1080,960);
  c_bcm->Divide(1,4);
  
  TGraph **g8 = new TGraph[8];

  // Plots order
  // Plots Raw data of
  //// bcms
  //// sams
  //// bpms
  
  // Plots Asym data of
  //// bcms
  //// sams

  /// cuts for bpm
  // run2347
  // TString cut_text ="(ev_num > 30e3 && ev_num<40e3)||(ev_num>70e3 && ev_num< 80e3)||(ev_num > 105e3 && ev_num<130e3)||(ev_num>145e3 && ev_num< 160e3)||(ev_num > 175e3 && ev_num<190e3)";

  //run2503
  //TString cut_text ="(ev_num > 10e3 && ev_num<14e3)||(ev_num>22e3 && ev_num< 28e3)||(ev_num > 32e3 && ev_num<40e3)||(ev_num>44e3 && ev_num< 50e3)";
  
  
  for(int ibcm = 0; ibcm<4; ibcm++){ // vqwk bcm index start from 0
    c_bcm->cd(ibcm+1);
    int n = R->Draw(Form("vqwk0_%d_cal:ev_num",ibcm),"ev_num>10","goff");
    g8[ibcm] = new TGraph(n,R->GetV2(),R->GetV1());
    g8[ibcm]->SetTitle("");
    g8[ibcm]->Draw("APL");

  }
  c_bcm->SaveAs(Form("run%d_bcm_raw.png",runnum));
  
  for(int isam = 0 ; isam < 8; isam++){
    c1->cd(isam+1);
    int n = R->Draw(Form("vqwk5_%d_cal:ev_num",isam),"ev_num>10","goff");
    g8[isam] = new TGraph(n,R->GetV2(),R->GetV1());
    g8[isam]->SetTitle("");
    g8[isam]->Draw("APL");
  }
  c1->SaveAs(Form("run%d_sam_raw.png",runnum));

  TString bpm_name[10] = {"bpm4bx","bpm4by","bpm4ax","bpm4ay","bpm12x","bpm12y","bpm8x","bpm8y","bpm14x","bpm14y"};
  TCanvas *c_bpm = new TCanvas("c_bpm","c_bpm",1080,3000);
  c_bpm->Divide(1,10);
  TGraph *g_bpm[10];

  
  for(int ibpm =0; ibpm<10;ibpm++){
    c_bpm->cd(ibpm+1);
    int n = R->Draw(Form("%s:ev_num",bpm_name[ibpm].Data()),"","goff");
    g_bpm[ibpm] = new TGraph(n,R->GetV2(),R->GetV1());
    g_bpm[ibpm]->SetTitle("");
    g_bpm[ibpm]->Draw("AP");
    TText *t1 = new TText(0,0.5,bpm_name[ibpm]);
    t1->SetNDC(kTRUE);
    t1->Draw("same");
  }
  c_bpm->SaveAs(Form("run%d_bpm_raw.png",runnum));


  
  for(int ibcm =1; ibcm<=4;ibcm++){ // real bcm index start from 1
    for(int icomb =0; icomb<7; icomb++){
      int n = R->Draw(Form("Aq%d_%d:ev_num",ibcm,icomb),"ev_num>10","goff");
      g8[icomb] = new TGraph(n,R->GetV2(),R->GetV1());
      c1->cd(icomb+1);
      g8[icomb]->SetTitle("");
      g8[icomb]->Draw("ALP");
      TText *t1 = new TText(0,0.5,histo_title[icomb]);
      t1->SetNDC(kTRUE);
      t1->Draw("same");
    }
    icomb = 7;
    int n = P->Draw(Form("Aq%d_%d:m_ev_num",ibcm,icomb),"m_ev_num>10","goff");
    g8[icomb] = new TGraph(n,P->GetV2(),P->GetV1());
    c1->cd(icomb+1);
    g8[icomb]->SetTitle("");
    g8[icomb]->Draw("ALP");
    TText *t1 = new TText(0,0.5,histo_title[icomb]);
    t1->SetNDC(kTRUE);
    t1->Draw("same");
    c1->SaveAs(Form("run%d_bcm%d_asym.png",runnum,ibcm));
  }
    
  for(int isam = 1; isam <=8;isam++){
    for(int icomb =0; icomb<7; icomb++){
      int n = R->Draw(Form("Ab%d_%d:ev_num",isam,icomb),"ev_num>10","goff");
      g8[icomb] = new TGraph(n,R->GetV2(),R->GetV1());
      c1->cd(icomb+1);
      g8[icomb]->SetTitle("");
      g8[icomb]->Draw("ALP");
      TText *t1 = new TText(0,0.5,histo_title[icomb]);
      t1->SetNDC(kTRUE);
      t1->Draw("same");
    }
    icomb = 7;
    int n = P->Draw(Form("Ab%d_%d:m_ev_num",isam,icomb),"m_ev_num>10","goff");
    g8[icomb] = new TGraph(n,P->GetV2(),P->GetV1());
    c1->cd(icomb+1);
    g8[icomb]->SetTitle("");
    g8[icomb]->Draw("ALP");
    TText *t1 = new TText(0,0.5,histo_title[icomb]);
    t1->SetNDC(kTRUE);
    t1->Draw("same");
    c1->SaveAs(Form("run%d_sam%d_asym.png",runnum,isam));
  }

  
}
