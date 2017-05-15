void PlotSamCorrelation(int runnum){

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(1);

  TString peak_name[4]={"101","001","110","010"};

  for(int ipeak =0; ipeak<4;ipeak++){
    TString file_name = Form("./ROOTfiles/parity16_%d_regress_%s.root",runnum,peak_name[ipeak].Data());
    TFile *reg_file = TFile::Open(file_name.Data());
    TTree *reg_tree = reg_file->Get("reg");
  
    TString sam_name[8];
    TString draw_text[2];  
    // (3+7) vs (1+5) vs (4+8)
    sam_name[0] = Form("(reg_asym_blumi%d + reg_asym_blumi%d)*0.5 -asym_bcm4",1,5);
    sam_name[1] = Form("(reg_asym_blumi%d + reg_asym_blumi%d)*0.5 -asym_bcm4",3,7);
    sam_name[2] = Form("(reg_asym_blumi%d + reg_asym_blumi%d)*0.5 -asym_bcm4",4,8);

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
      sam_name[isam] = Form("reg_asym_blumi%d -asym_bcm4",isam+1);
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
      sam_name[isam] = Form("reg_asym_blumi%d -asym_bcm4",isam+1);
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
  reg_tree->Draw(Form("(%s) -(%s)>>%s",draw_text[0].Data(),draw_text[1].Data(),histo_name.Data()),"ok_cut");
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
