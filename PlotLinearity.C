void PlotLinearity(int runnum){
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(1);
  
  Plot4Peak(runnum,custom_cut(runnum).Data());  // Plot 4 peaks in Asym_bcms

  PlotCorrelation(runnum);  //Plot asym_sams vs asym_bcm
  
  Regression4peak(runnum);//Regression for each 4 peak
  
  Plot4PeakScan(runnum);

}

void Plot4PeakScan(int runnum){
  
  TString leg_text[4]={"101","001","110","010"};
  TString histo_name;
  TCanvas *c3 = new TCanvas("c3","c3",1200,400);
  c3->Divide(3,1);
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  TH2D *hsambcmcolz;
  TH1D *hsam;
  TH1D *hbcm;
  double reg_asym_sam[8][4]; // 8 regressed SAMs with 4 peaks mean value;
  double reg_asym_sam_err[8][4]; // 8 regressed SAMs with 4 peaks mean value;

  double reg_asym_n_sam[8][4][4]; // [isam][ibcm][ipeak]
  double reg_asym_n_sam_err[8][4][4];

  double asym_bcm[4][4];
  double asym_bcm_err[4][4]; 
  
  for(int ipeak=0;ipeak<4;ipeak++){
    TFile *reg_file = TFile::Open(Form("ROOTfiles/parity16_%d_regress_%s.root",runnum,leg_text[ipeak].Data()));
    TTree *reg_tree = reg_file->Get("reg");
    
    for(int isam =1; isam<=8;isam++){
	histo_name = Form("hregblumi%d",isam);
	reg_tree->Draw(Form("reg_asym_blumi%d>>%s",isam,histo_name.Data()),"ok_cut","goff");
	hsam = (TH1D*)gDirectory->FindObject(histo_name);
	c1->cd();
	hsam->Draw();
	FitGaussian(hsam);
	c1->SaveAs(Form("h%d.pdf",isam));
	reg_asym_sam[isam-1][ipeak] = hsam->GetMean();
	reg_asym_sam_err[isam-1][ipeak] = hsam->GetRMS()/TMath::Sqrt(hsam->GetEntries());
    }
    gSystem->Exec(Form("pdfunite h*.pdf run%d_regress_sam_1d_%s.pdf",runnum,leg_text[ipeak].Data()));
    gSystem->Exec("rm h*.pdf");
    
    for(int ibcm =1;ibcm<=4;ibcm++){
      for(int isam =1; isam<=8;isam++){
	histo_name = Form("hregblumi%dnbcm%d",isam,ibcm);
	reg_tree->Draw(Form("reg_asym_blumi%d-asym_bcm%d>>%s",isam,ibcm,histo_name.Data()),"ok_cut","goff");
	hsam = (TH1D*)gDirectory->FindObject(histo_name);
	c1->cd();
	hsam->Draw();
	FitGaussian(hsam);
	c1->SaveAs(Form("h%d.pdf",isam));
	reg_asym_n_sam[isam-1][ibcm-1][ipeak] = hsam->GetMean();
	reg_asym_n_sam_err[isam-1][ibcm-1][ipeak] = hsam->GetRMS()/TMath::Sqrt(hsam->GetEntries());
      }
      gSystem->Exec(Form("pdfunite h*.pdf run%d_regress_asym_sam-asym_bcm%d_1d_%s.pdf",runnum,ibcm,leg_text[ipeak].Data()));
      gSystem->Exec("rm h*.pdf");
    }
    
    for(int ibcm =1; ibcm<=4;ibcm++){
      histo_name = Form("hbcm%d",ibcm);
      reg_tree->Draw(Form("asym_bcm%d>>%s",ibcm,histo_name.Data()),"ok_cut","goff");
      hbcm = (TH1D*)gDirectory->FindObject(histo_name);
      c1->cd();
      hbcm->Draw();
      FitGaussian(hbcm);
      c1->SaveAs(Form("h%d.pdf",ibcm));
      asym_bcm[ibcm-1][ipeak] = hbcm->GetMean();
      asym_bcm_err[ibcm-1][ipeak] = hbcm->GetRMS()/TMath::Sqrt(hbcm->GetEntries());
    }
    gSystem->Exec(Form("pdfunite h*.pdf run%d_asym_bcm_1d_%s.pdf",runnum,leg_text[ipeak].Data()));
    gSystem->Exec("rm h*.pdf");

    reg_file->Close();
  }

  //Plot Linear Scan and fitting;
  TGraphErrors *g_lin[8];
  TCanvas *c42 = new TCanvas("c42","c42",1600,1200);
  c42->Divide(2,4);
  for(int ibcm =0;ibcm<4;ibcm++){
    for(int isam = 0;isam<8;isam++){
      g_lin[isam] = new TGraphErrors(4,asym_bcm[ibcm],reg_asym_sam[isam],asym_bcm_err[ibcm],reg_asym_sam_err[isam]);
      g_lin[isam]->SetTitle(Form("SAM %d",isam+1));
      g_lin[isam]->GetYaxis()->SetTitle(Form("reg_asym_blumi%d/ppm",isam+1));
      g_lin[isam]->GetXaxis()->SetTitle(Form("asym_bcm%d/ppm",ibcm+1));
      c42->cd(isam+1);
      g_lin[isam]->Draw("AP");
      g_lin[isam]->Fit("pol1","Q");
    }
    c42->SaveAs(Form("run%d_4peak_bcm%d.pdf",runnum,ibcm+1));
  }				   


  //Plot Residue and Fitting;
  TGraphErrors *g_res[8];

  for(int ibcm =0;ibcm<4;ibcm++){
    for(int isam = 0;isam<8;isam++){
      g_res[isam] = new TGraphErrors(4,asym_bcm[ibcm],reg_asym_n_sam[isam][ibcm],asym_bcm_err[ibcm],reg_asym_n_sam_err[isam][ibcm]);
      g_res[isam]->SetTitle(Form("SAM %d",isam+1));
      g_res[isam]->GetYaxis()->SetTitle(Form("reg_asym_blumi%d-asym_bcm%d/ppm",isam+1,ibcm+1));
      g_res[isam]->GetXaxis()->SetTitle(Form("asym_bcm%d/ppm",ibcm+1));
      c42->cd(isam+1);
      g_res[isam]->Draw("AP");
      g_res[isam]->Fit("pol1","Q");
    }
    c42->SaveAs(Form("run%d_4peak_bcm%d_residue.pdf",runnum,ibcm+1));
  }				   

}

void Regression4peak(int runnum){
  TString parity_cut[4] = {"&&evt_pairsynch[0]==0&&prev_hel==1",
			   "&&evt_pairsynch[0]==0&&prev_hel==0",
			   "&&evt_pairsynch[1]==0&&prev_hel==1",
			   "&&evt_pairsynch[1]==0&&prev_hel==0"};
  TString leg_text[4]={"101","001","110","010"};
  TString postpan_dir = "/home/yetao/workarea/jlab_2016/postpan";
  
  for(int ipeak=0;ipeak<4;ipeak++){
    gSystem->Exec(Form("cp /home/yetao/workarea/jlab_2016/postpan/control.conf_reg%d /home/yetao/workarea/jlab_2016/postpan/control.conf_reg%d_%s",runnum,runnum,leg_text[ipeak].Data()));
    FILE* conf_file;
    conf_file = fopen(Form("%s/control.conf_reg%d_%s",postpan_dir.Data(),runnum,leg_text[ipeak].Data()),"a");
    fprintf(conf_file,"%s\n",parity_cut[ipeak].Data());
    fclose(conf_file);
    gSystem->Exec(Form("%s/redana -r %d -C %s/control.conf_reg%d_%s",postpan_dir.Data(),runnum,postpan_dir.Data(),runnum,leg_text[ipeak].Data()));
    gSystem->Exec(Form("mv ./ROOTfiles/parity16_%d_regress.root ./ROOTfiles/parity16_%d_regress_%s.root",runnum,runnum,leg_text[ipeak].Data()));
  }
}

void PlotCorrelation(int runnum){
  TFile *reg_file = TFile::Open(Form("ROOTfiles/parity16_%d_regress.root",runnum));
  TTree *reg_tree = reg_file->Get("reg");

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  TString histo_name;
  TH2D *h2fit;
  
  for(int ibcm =1; ibcm<=4;ibcm++){
    for(int isam =1; isam<=8;isam++){
      c2->cd(1);
      histo_name = Form("hblumi%dbcm%dcolz",isam,ibcm);
      reg_tree->Draw(Form("asym_blumi%d:asym_bcm%d>>%s",isam,ibcm,histo_name.Data()),"ok_cut","COLZ");
      c2->cd(2);
      histo_name = Form("hblumi%dbcm%dprof",isam,ibcm);
      reg_tree->Draw(Form("asym_blumi%d:asym_bcm%d>>%s",isam,ibcm,histo_name.Data()),"ok_cut","prof");
      h2fit = (TH2D*)gDirectory->FindObject(histo_name);
      FitLinear(h2fit);
      c2->SaveAs(Form("h%d.pdf",isam));
    }
    gSystem->Exec(Form("pdfunite h*.pdf run%d_raw_asym_sams_vs_bcm%d.pdf",runnum,ibcm));
    gSystem->Exec("rm h*.pdf");
  }

  //residue correlation
  for(int ibcm =1; ibcm<=4;ibcm++){
    for(int isam =1; isam<=8;isam++){
      c2->cd(1);
      histo_name = Form("hnblumi%dbcm%dcolz",isam,ibcm);
      reg_tree->Draw(Form("asym_blumi%d-asym_bcm%d:asym_bcm%d>>%s",isam,ibcm,ibcm,histo_name.Data()),"ok_cut","COLZ");
      c2->cd(2);
      histo_name = Form("hnblumi%dbcm%dprof",isam,ibcm);
      reg_tree->Draw(Form("asym_blumi%d-asym_bcm%d:asym_bcm%d>>%s",isam,ibcm,ibcm,histo_name.Data()),"ok_cut","prof");
      h2fit = (TH2D*)gDirectory->FindObject(histo_name);
      FitLinear(h2fit);
      c2->SaveAs(Form("h%d.pdf",isam));
    }
    gSystem->Exec(Form("pdfunite h*.pdf run%d_raw_asym_n_sam_vs_bcm%d.pdf",runnum,ibcm));
    gSystem->Exec("rm h*.pdf");
  }

    //residue correlation after regression
  for(int ibcm =1; ibcm<=4;ibcm++){
    for(int isam =1; isam<=8;isam++){
      c2->cd(1);
      histo_name = Form("hregblumi%dbcm%dcolz",isam,ibcm);
      reg_tree->Draw(Form("reg_asym_blumi%d-asym_bcm%d:asym_bcm%d>>%s",isam,ibcm,ibcm,histo_name.Data()),"ok_cut","COLZ");
      c2->cd(2);
      histo_name = Form("hregblumi%dbcm%dprof",isam,ibcm);
      reg_tree->Draw(Form("reg_asym_blumi%d-asym_bcm%d:asym_bcm%d>>%s",isam,ibcm,ibcm,histo_name.Data()),"ok_cut","prof");
      h2fit = (TH2D*)gDirectory->FindObject(histo_name);
      FitLinear(h2fit);
      c2->SaveAs(Form("h%d.pdf",isam));
    }
    gSystem->Exec(Form("pdfunite h*.pdf run%d_reg_asym_n_sam_vs_bcm%d.pdf",runnum,ibcm));
    gSystem->Exec("rm h*.pdf");
  }
  
  reg_file->Close();
}

void Plot4Peak(int runnum, TString custom_cut){
  TFile *panfile = TFile::Open(Form("ROOTfiles/parity16_%d_standard.root",runnum));
  TTree *parity_tree = panfile->Get("P");
  
  TString parity_cut[4] = {"ok_cut&&evt_pairsynch[0]==0&&prev_hel==1&&",
			   "ok_cut&&evt_pairsynch[0]==0&&prev_hel==0&&",
			   "ok_cut&&evt_pairsynch[1]==0&&prev_hel==1&&",
			   "ok_cut&&evt_pairsynch[1]==0&&prev_hel==0&&"};

  for(int i=0;i<4;i++){
    parity_cut[i] = parity_cut[i] +custom_cut;
  }
  TString leg_text[4]={"101","001","110","010"};
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TH1D *h4peak[4];
  TString histo_name;
  // BCM
  for(int ibcm = 1; ibcm<=4; ibcm++){
    TLegend *leg = new TLegend(0.1,0.9,0.3,0.7);
    for(int ipeak =0;ipeak<4;ipeak++){
      histo_name = Form("h4peak%d",ipeak);
      parity_tree->Draw(Form("asym_bcm%d>>%s",ibcm,histo_name.Data()),parity_cut[ipeak].Data(),"goff");
      h4peak[ipeak] = (TH1D*)gDirectory->FindObject(histo_name);
      h4peak[ipeak]->SetTitle(Form("asym_bcm%d",ibcm));
      h4peak[ipeak]->SetLineColor(ipeak+1);
      leg->AddEntry(h4peak[ipeak],leg_text[ipeak],"l");
    }
    c1->cd();
    h4peak[0]->Draw("");
    h4peak[1]->Draw("same");h4peak[2]->Draw("same");h4peak[3]->Draw("same");
    leg->Draw();
    c1->SaveAs(Form("h%d.pdf",ibcm));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_bcm_4peak.pdf",runnum));
  gSystem->Exec("rm h*.pdf");
  // SAM
  for(int isam = 1; isam<=8; isam++){
    TLegend *leg = new TLegend(0.1,0.9,0.3,0.7);
    for(int ipeak =0;ipeak<4;ipeak++){
      histo_name = Form("h4peak%d%d",isam,ipeak);
      parity_tree->Draw(Form("asym_blumi%d>>%s",isam,histo_name.Data()),parity_cut[ipeak].Data(),"goff");
      h4peak[ipeak] = (TH1D*)gDirectory->FindObject(histo_name);
      h4peak[ipeak]->SetTitle(Form("asym_sam%d",isam));
      h4peak[ipeak]->SetLineColor(ipeak+1);
      leg->AddEntry(h4peak[ipeak],leg_text[ipeak],"l");
    }
    c1->cd();
    h4peak[0]->Draw("");
    h4peak[1]->Draw("same");h4peak[2]->Draw("same");h4peak[3]->Draw("same");
    leg->Draw();
    c1->SaveAs(Form("h%d.pdf",isam));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_sam_4peak.pdf",runnum));
  gSystem->Exec("rm h*.pdf");

  panfile->Close();
}

TString custom_cut(int runnum){
  TString text_cut;
  if(runnum ==2503)
    text_cut ="((m_ev_num > 10e3 && m_ev_num<14e3)||(m_ev_num>22e3 && m_ev_num< 28e3)||(m_ev_num > 32e3 && m_ev_num<40e3)||(m_ev_num>44e3 && m_ev_num< 50e3))";
  if(runnum == 2721)
    text_cut =" ((m_ev_num>10 && m_ev_num< 6000)||(m_ev_num>19e3 && m_ev_num<27000))";
  if(runnum ==2347)
    text_cut ="((m_ev_num > 30e3 && m_ev_num < 40e3) ||(m_ev_num>70e3 && m_ev_num < 80e3) ||(m_ev_num>105e3 && m_ev_num < 130e3) ||(m_ev_num>145e3 && m_ev_num < 160e3) ||(m_ev_num>175e3 && m_ev_num < 190e3))";
  else
    cout << "No cuts found" << endl;

  return text_cut;
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

void FitGaussian(TH1 *h1){
  TF1 *g1 = new TF1("g1","gaus",-10e5,10e5);
  double par[3]={0,0,0}; // fitting parameter
  par[1] = h1->GetMean();
  double rms = h1->GetRMS();
  par[2] = TMath::Sqrt(rms);
  par[0] = h1->GetBinContent(h1->GetMaximumBin());

  g1->SetParameters(par);
  h1->Fit("g1","QR+","",par[1]-3*rms,par[1]+3*rms);
}
