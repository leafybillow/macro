void GetPedestalNoise(int runnum){
  
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(1);
  TFile *panfile = TFile::Open(Form("./ROOTfiles/parity16_%d_standard.root",runnum));
  //run2503
  TString parity_cut = "(m_ev_num>10 && m_ev_num<4e3)||(m_ev_num>6e3 && m_ev_num<8e3)||(m_ev_num>17e3 && m_ev_num<18e3)||(m_ev_num>29100 && m_ev_num<29800)||(m_ev_num>41e3 &&m_ev_num<41500)";
  //run2721
  //  TString parity_cut = "(m_ev_num>10e3 && m_ev_num<11e3)||(m_ev_num>15e3 && m_ev_num<16e3)||(m_ev_num>28e3 && m_ev_num<30e3)";

  P->Draw(">>cutlist",parity_cut.Data());
  TEventList *cutlist = (TEventList*)gDirectory->Get("cutlist");
  P->SetEventList(cutlist);
  
  TH1D *hblumi[8];
  double evt_0, evt_1;
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->cd();
  
  P->Draw("(evt_bcm4[0]-evt_bcm4[1])-(evt_bcm3[0]-evt_bcm3[1])>>htbcm","","");
  // TH1D *hbcm = (TH1D*)gDirectory->FindObject("htbcm");
  // FitGaussian(hbcm);
  c1->SaveAs("bcm_pedestal_corre.pdf");
  
  for(int iblumi = 0; iblumi<8; iblumi++){
    P->Draw(Form("(evt_blumi%d[0]-evt_blumi%d[1])>>h%d",iblumi+1,iblumi+1,iblumi),"","goff");
    hblumi[iblumi] = (TH1D*)gDirectory->FindObject(Form("h%d",iblumi));
    c1->cd();
    hblumi[iblumi]->Draw();
    FitGaussian(hblumi[iblumi]);
    c1->SaveAs(Form("h%d.pdf",iblumi));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_blumi_pedestal_noise.pdf",runnum));
  gSystem->Exec("rm h*.pdf");

  for(int iblumi = 0; iblumi<8; iblumi++){
    c1->cd();
    P->Draw(Form("evt_blumi%d[0]-evt_blumi%d[1]:m_ev_num",iblumi+1,iblumi+1),"","COL");
    c1->SaveAs(Form("h%d.pdf",iblumi));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_blumi_pedestal_noise_2d.pdf",runnum));
  gSystem->Exec("rm h*.pdf");

  panfile->Close();
}

void FitGaussian(TH1D *h1){
  TF1 *gaus = new TF1("gaus","gaus",-10e3,10e3);
  double par[3]={0,0,0}; // fitting parameter
  par[1] = h1->GetMean();
  double rms = h1->GetRMS();
  par[2] = TMath::Sqrt(rms);
  par[0] = h1->GetBinContent(h1->GetMaximumBin());
  
  gaus->SetParameters(par);
  h1->Fit("gaus","QR+","",par[1]-2*rms,par[1]+2*rms);
}
