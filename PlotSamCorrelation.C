//Global Variables

TString peak_name[4]={"101","001","110","010"};
TString bpm_title[10] = {"4bx","4by","4ax","4ay","14x","14y","12x","12y","8x","8y"};

void PlotSamCorrelation(int runnum){

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(1);

  Plot_Sam_vs_Sam(runnum);
  Plot_Sam_vs_Bpm(runnum);
  
  for(int ipeak =0; ipeak<4;ipeak++){
    TString file_name = Form("./ROOTfiles/parity16_%d_regress_%s.root",runnum,peak_name[ipeak].Data());
    TFile *reg_file = TFile::Open(file_name.Data());
    TTree *reg_tree = reg_file->Get("reg");
  
    TString sam_name[8];
    TString draw_text[2];  
    // (3+7) vs (1+5) vs (4+8)
    sam_name[0] = Form("(reg_asym_n_blumi%d + reg_asym_n_blumi%d)",1,5);
    sam_name[1] = Form("(reg_asym_n_blumi%d + reg_asym_n_blumi%d)",3,7);
    sam_name[2] = Form("(reg_asym_n_blumi%d + reg_asym_n_blumi%d)",4,8);

    draw_text[0] = sam_name[0];
    draw_text[1] = sam_name[1];
    Plot4Can(reg_tree,draw_text,"comb",1);

    draw_text[0] = sam_name[0];
    draw_text[1] = sam_name[2];
    Plot4Can(reg_tree,draw_text,"comb",2);

    draw_text[0] = sam_name[1];
    draw_text[1] = sam_name[2];
    Plot4Can(reg_tree,draw_text,"comb",3);

    gSystem->Exec(Form("pdfunite h*.pdf run%d_sam_resolution_%s.pdf",runnum,peak_name[ipeak].Data()));
    gSystem->Exec("rm h*.pdf");
		
    // 1-5-like plots
  
    for(int isam = 0; isam<8;isam++){
      sam_name[isam] = Form("reg_asym_n_blumi%d",isam+1);
    }

    draw_text[0] = sam_name[0];
    draw_text[1] = sam_name[4];
    Plot2Can(reg_tree,draw_text,"sub",1);

    draw_text[0] = sam_name[2];
    draw_text[1] = sam_name[6];
    Plot2Can(reg_tree,draw_text,"sub",2);

    draw_text[0] = sam_name[3];
    draw_text[1] = sam_name[7];
    Plot2Can(reg_tree,draw_text,"sub",3);

    gSystem->Exec(Form("pdfunite h*.pdf run%d_sam_double_diff_%s.pdf",runnum,peak_name[ipeak].Data()));
    gSystem->Exec("rm h*.pdf");
  
    // Plot BCM double Difference
    draw_text[0] = "asym_bcm3";
    draw_text[1] = "asym_bcm4";
    Plot4Can(reg_tree,draw_text,"bcm",1);

    draw_text[0] = "asym_bcm1";
    draw_text[1] = "asym_bcm2";
    Plot4Can(reg_tree,draw_text,"bcm",2);

    gSystem->Exec(Form("pdfunite h*.pdf run%d_bcm_double_diff_%s.pdf",runnum,peak_name[ipeak].Data()));
    gSystem->Exec("rm h*.pdf");
  
    // 1+5 vs 1-5 like
  
    for(int isam = 0; isam<8;isam++){
      sam_name[isam] = Form("reg_asym_n_blumi%d",isam+1);
    }
  
    for(int iplot = 0; iplot<4; iplot++){
      draw_text[0] = sam_name[iplot];
      draw_text[1] = sam_name[iplot+4];
      Plot4Can(reg_tree,draw_text,"pm",iplot);
    }

    gSystem->Exec(Form("pdfunite h*.pdf run%d_cross_check_%s.pdf",runnum,peak_name[ipeak].Data()));
    gSystem->Exec("rm h*pdf");

    reg_file->Close();
    
  }


}

void Plot_Sam_vs_Sam(int runnum){
  for(int ipeak =0;ipeak<4;ipeak++){
    TFile *reg_file = TFile::Open(Form("./ROOTfiles/parity16_%d_regress_%s.root",
				       runnum,peak_name[ipeak].Data()));
    TTree *reg_tree = (TTree*)reg_file ->Get("reg");
    TCanvas *c99 = new TCanvas("c99","c99",900,900);
    c99->Divide(9,9);

    TString sam_text1, sam_text2;
    for(int isam =0;isam<8;isam++){
      for(int jsam=0;jsam<8;jsam++){
	c99->cd((isam+1)*9+jsam+2);
	sam_text1 = Form("reg_asym_n_blumi%d",isam+1);
	sam_text2 = Form("reg_asym_n_blumi%d",jsam+1);
	reg_tree->Draw(Form("%s:%s",sam_text1.Data(),sam_text2.Data()),"ok_cut","COLZ");
      }
    }
    c99->SaveAs(Form("h%d.pdf",ipeak));
    reg_file->Close();
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_sam_vs_sam_colz.pdf",runnum));
  gSystem->Exec("rm h*.pdf");
  
  for(int ipeak =0;ipeak<4;ipeak++){
    TFile *reg_file = TFile::Open(Form("./ROOTfiles/parity16_%d_regress_%s.root",
				       runnum,peak_name[ipeak].Data()));
    TTree *reg_tree = (TTree*)reg_file->Get("reg");
    TCanvas *c99 = new TCanvas("c99","c99",900,900);
    c99->Divide(9,9);

    TString sam_text1, sam_text2;
    for(int isam =0;isam<8;isam++){
      for(int jsam=0;jsam<8;jsam++){
	c99->cd((isam+1)*9+jsam+2);
	sam_text1 = Form("reg_asym_n_blumi%d",isam+1);
	sam_text2 = Form("reg_asym_n_blumi%d",jsam+1);
	reg_tree->Draw(Form("%s:%s",sam_text1.Data(),sam_text2.Data()),"ok_cut","prof");
      }
    }
    c99->SaveAs(Form("h%d.pdf",ipeak));
    reg_file->Close();
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_sam_vs_sam_prof.pdf",runnum));
  gSystem->Exec("rm h*.pdf");
  
}

void Plot_Sam_vs_Bpm(int runnum){
  
  TString sam_title;
  for(int ipeak =0; ipeak<4;ipeak++){
    TFile *reg_file = TFile::Open(Form("./ROOTfiles/parity16_%d_regress_%s.root",
				       runnum,peak_name[ipeak].Data()));
    TTree *reg_tree = reg_file->Get("reg");

    TCanvas *c911 = new TCanvas("c911","c911",1100,900);
    c911->Divide(11,9);
    for(int iblumi =1; iblumi<=8;iblumi++){
      c911->cd(iblumi*11+1);
      TText *t_blumi = new TText(0.4,0.4,Form("%d",iblumi));
      t_blumi->SetNDC(kTRUE);
      t_blumi->SetTextSize(0.5);
      t_blumi->Draw();
    }
    for(int ibpm =0;ibpm<10;ibpm++){
      c911->cd(ibpm+2);
      TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
      t_bpm->SetNDC(kTRUE);
      t_bpm->SetTextSize(0.4);
      t_bpm->Draw();
    }

    c911->cd(1);
    TText *t_peak = new TText(0.4,0.4,Form("%s",peak_name[ipeak].Data()));
    t_peak->SetNDC(kTRUE);
    t_peak->SetTextSize(0.2);
    t_peak->Draw();

    for(int isam =0;isam<8;isam++){
      for(int ibpm =0; ibpm<10;ibpm++){
	c911->cd((isam+1)*11+ibpm+2);
	sam_title = Form("reg_asym_n_blumi%d",isam+1);
	reg_tree->Draw(Form("%s:diff_bpm%s",sam_title.Data(),bpm_title[ibpm].Data()),
		       "ok_cut","COLZ");
      }
    }
    c911->SaveAs(Form("h%d.pdf",ipeak));
    reg_file->Close();
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_sam_vs_bpm_colz_4peak.pdf",runnum));
  gSystem->Exec("rm h*.pdf");
 
}

void Plot2Can(TTree *reg_tree,TString *draw_text,TString histo_pre,int plot_label){
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  TString histo_name;
  
  c2->cd(1);
  reg_tree->Draw(Form("(%s):(%s)",draw_text[0].Data(),draw_text[1].Data()),"ok_cut","COLZ");

  c2->cd(2);
  histo_name = Form("h%s%d",histo_pre.Data(),plot_label);
  reg_tree->Draw(Form("(%s)-(%s)>>%s",draw_text[0].Data(),draw_text[1].Data(),histo_name.Data()),"ok_cut");
  TH1D *hfit = (TH1D*)gDirectory->FindObject(histo_name);
  FitGaussian(hfit);

  c2->SaveAs(Form("h%d.pdf",plot_label));
  
}

void Plot4Can(TTree *reg_tree,TString *draw_text,TString histo_pre,int plot_label){
  TString histo_name;
  TString draw_title;
  TH1D *hfit;
  TCanvas *c4 = new TCanvas("c4","c4",1200,1200);
  c4->Divide(2,2);

  c4->cd(1);
  histo_name = Form("hcol%s%d",histo_pre.Data(),plot_label);
  draw_title = Form("(%s):(%s)>>%s",draw_text[0].Data(),draw_text[1].Data(),histo_name.Data());
  reg_tree->Draw(draw_title,"ok_cut","COLZ");
    
  c4->cd(2);
  histo_name = Form("h2%s%d",histo_pre.Data(),plot_label);
  reg_tree->Draw(Form("%s>>%s",draw_text[0].Data(),histo_name.Data()),"ok_cut");
  hfit = (TH1D*)gDirectory->FindObject(histo_name);
  FitGaussian(hfit);
  
  c4->cd(3);
  histo_name = Form("h3%s%d",histo_pre.Data(),plot_label);
  reg_tree->Draw(Form("%s>>%s",draw_text[1].Data(),histo_name.Data()),"ok_cut");
  hfit = (TH1D*)gDirectory->FindObject(histo_name);
  FitGaussian(hfit);
  
  c4->cd(4);
  histo_name = Form("h4%s%d",histo_pre.Data(),plot_label);
  reg_tree->Draw(Form("(%s)-(%s)>>%s",draw_text[0].Data(),draw_text[1].Data(),histo_name.Data()),"ok_cut");
  hfit = (TH1D*)gDirectory->FindObject(histo_name);
  FitGaussian(hfit);
  
  c4->SaveAs(Form("h%d.pdf",plot_label));
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

double FitLinear(TH2D *h1){
  TF1 *l1 = new TF1("l1","[0]*x+[1]",-10e5,10e5);
  l1->SetParNames("slope","intercept");
  l1->SetParameters(0,1);
  double rms = h1->GetRMS(1);
  double mean = h1->GetMean(1);
  h1->Fit("l1","QR+","",mean-2*rms,mean+2*rms);
  return l1->GetParameter(0);
}
