void PlotBCMBlumiCorrelation(int runnum){
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(1);

  TFile *panfile = TFile::Open(Form("ROOTfiles/parity16_%d_standard.root",runnum));

  //run 2721
  //TString cut_parity = "(ok_cut)&&((m_ev_num > 10 && m_ev_num < 5e3) ||(m_ev_num>20e3 && m_ev_num < 25e3))";
  //run2503 cut
  TString cut_parity="(m_ev_num > 10e3 && m_ev_num<14e3)||(m_ev_num>22e3 && m_ev_num< 28e3)||(m_ev_num > 32e3 && m_ev_num<40e3)||(m_ev_num>44e3 && m_ev_num< 50e3)";
  TTree *ptree = panfile->Get("P");
  ptree->Draw(">>cutlist",cut_parity);
  TEventList *cutlist = (TEventList*)gDirectory->Get("cutlist");
  ptree->SetEventList(cutlist);

  TCanvas *c2 = new TCanvas("c2","c2",1600,800);
  c2->Divide(2,1);

  // BCMs Correlations
  c2->cd(1);
  ptree->Draw("asym_bcm1:asym_bcm2>>hbcmcol12","","COLZ");

  c2->cd(2);
  ptree->Draw("asym_bcm1:asym_bcm2>>hbcmprof12","","PROFS");
  TH2D *h2 = (TH2D*)gDirectory->FindObject("hbcmprof12");
  FitLinear(h2);
  c2->SaveAs("h1.pdf");

  c2->cd(1);
  ptree->Draw("asym_bcm3:asym_bcm4>>hbcmcol34","","COLZ");

  c2->cd(2);
  ptree->Draw("asym_bcm3:asym_bcm4>>hbcmprof34","","PROFS");
  TH2D *h2 = (TH2D*)gDirectory->FindObject("hbcmprof34");
  FitLinear(h2);
  c2->SaveAs("h2.pdf");

  gSystem->Exec(Form("pdfunite h*.pdf run%d_bcm_correlation.pdf",runnum));
  gSystem->Exec("rm h*.pdf");
    
  for(int iblumi=1;iblumi<=8;iblumi++){
    for(int ibcm =1; ibcm<=4;ibcm++){
      c2->cd(1);
      ptree->Draw(Form("asym_blumi%d:asym_bcm%d>>hblumicol%d%d",iblumi,ibcm,iblumi,ibcm),"","COLZ");
      c2->cd(2);
      TString prof_name = Form("hblumprof%d%d",iblumi,ibcm);
      ptree->Draw(Form("asym_blumi%d:asym_bcm%d>>%s",iblumi,ibcm,prof_name.Data()),"","prof");
      TH2D *h2 = (TH2D*)gDirectory->FindObject(prof_name.Data());
      FitLinear(h2);
      c2->SaveAs(Form("%s.pdf",prof_name.Data()));
    }
    gSystem->Exec(Form("pdfunite h*.pdf run%d_blumi%d_correlation.pdf",runnum,iblumi));
    gSystem->Exec("rm h*.pdf");
  }
}

double FitLinear(TH2D *h1){
  TF1 *l1 = new TF1("l1","[0]*x+[1]",-10e5,10e5);
  l1->SetParNames("slope","intercept");
  l1->SetParameters(0,1);
  double rms = h1->GetRMS(1);
  double mean = h1->GetMean(1);
  h1->Fit("l1","QR+","",mean-2*rms,mean+2*rms);
  return l1->GetParameter(0);
}

