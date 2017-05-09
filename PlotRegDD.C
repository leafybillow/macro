void PlotRegDD(int runnum){
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(1);
  /// run .x CombineRawReg
  TFile *regfile = TFile::Open(Form("./ROOTfiles/parity16_%d_combine.root",runnum));

  // plot normalized blumis correlation
  int blumi_a[4]={1,2,3,4};
  int blumi_b[4]={5,6,7,8};
  TString histo_name;
  TString tree_name;
  TString cut_text="ok_cut";
  
  // Plot Correlation between single SAM with Position Difference
  TCanvas *c810 = new TCanvas("c810","c810",1100,900);
  c810->Divide(11,9);
  TString bpm_name[10] = {"diff_bpm4bx","diff_bpm4by","diff_bpm4ax","diff_bpm4ay","diff_bpm14x","diff_bpm14y","diff_bpm12x","diff_bpm12y","diff_bpm8x","diff_bpm8y"};
  TString bpm_title[10] = {"4bx","4by","4ax","4ay","14x","14y","12x","12y","8x","8y"};
  TString bcm_title[4] = {"bcm1","bcm2","bcm3","bcm4"};

  TString bpm_old[10];
  TString suffix_text = "_old";
  for(int ibpm =0; ibpm<10;ibpm++){
    bpm_old[ibpm] = bpm_name[ibpm]+suffix_text;
  }

  //vector for printing table
  //  vector< vector<double> > vec2d_reg_blumi_bpm_coeff;
  vector< vector<double> > vec2d_reg_blumi_bpm_slope;
  vector< vector<double> > vec2d_raw_blumi_bpm_slope;
  vector< vector<double> > vec2d_reg_blumi_blumi_slope;
  vector< vector<double> > vec2d_raw_blumi_blumi_slope;
  vector< vector<double> > vec2d_bpm_bpm_slope;

  //  vector<double> vec_reg_blumi_bpm_coeff;
  vector<double> vec_reg_blumi_bpm_slope;
  vector<double> vec_raw_blumi_bpm_slope;
  vector<double> vec_reg_blumi_blumi_slope;
  vector<double> vec_raw_blumi_blumi_slope;
  vector<double> vec_bpm_bpm_slope;

  vector<double>  vec_bpm_rms;
  vector<double>  vec_raw_blumi_rms;
  vector<double>  vec_reg_blumi_rms;
  
  //////COLZ
  for(int iblumi = 1; iblumi<=8;iblumi++){
    for(int ibpm =0 ;ibpm<10;ibpm++){
      c810->cd((iblumi)*11+ibpm+2);
      TString blumi_name = Form("reg_asym_n_blumi%d",iblumi);
      TString tree_name = Form("%s:%s",blumi_name.Data(),bpm_old[ibpm].Data());
      reg->Draw(tree_name.Data(),cut_text,"COL");
    }
  }

  for(int iblumi =1; iblumi<=8;iblumi++){
    c810->cd(iblumi*11+1);
    TText *t_blumi = new TText(0.4,0.4,Form("%d",iblumi));
    t_blumi->SetNDC(kTRUE);
    t_blumi->SetTextSize(0.5);
    t_blumi->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    c810->cd(ibpm+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }
 
  c810->SaveAs(Form("run%d_reg_blumi_bpm_colz.pdf",runnum));

  ///// PROFILE
  for(int iblumi = 1; iblumi<=8;iblumi++){
    for(int ibpm =0 ;ibpm<10;ibpm++){
      c810->cd((iblumi)*11+ibpm+2);
      TString histo_name = Form("hblumivsbpm%d%d",iblumi,ibpm);
      TString blumi_name = Form("reg_asym_n_blumi%d",iblumi);
      TString tree_name = Form("%s:%s>>%s",blumi_name.Data(),bpm_old[ibpm].Data(),histo_name.Data());
      reg->Draw(tree_name.Data(),cut_text,"prof");
      TH2D *hfit = (TH2D*)gDirectory->FindObject(histo_name);
      double slope = FitLinear(hfit);
      vec_reg_blumi_bpm_slope.push_back(slope);
    }
    vec2d_reg_blumi_bpm_slope.push_back(vec_reg_blumi_bpm_slope);
    vec_reg_blumi_bpm_slope.clear();
  }

  for(int iblumi =1; iblumi<=8;iblumi++){
    c810->cd(iblumi*11+1);
    TText *t_blumi = new TText(0.4,0.4,Form("%d",iblumi));
    t_blumi->SetNDC(kTRUE);
    t_blumi->SetTextSize(0.5);
    t_blumi->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    c810->cd(ibpm+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }
 
  c810->SaveAs(Form("run%d_reg_blumi_bpm_prof.pdf",runnum));


  // Plot Correlation between combined SAMs with Position Difference
  
  TCanvas *c810 = new TCanvas("c810","c810",1100,900);
  c810->Divide(11,9);

  //////COLZ
  for(int ibpm = 0; ibpm<10;ibpm++){
    for(int iblumi = 0; iblumi<4;iblumi++){
      c810->cd((iblumi+1)*11+ibpm+2);
      TString blumi_name = Form("reg_asym_n_blumi%d+reg_asym_n_blumi%d",blumi_a[iblumi],blumi_b[iblumi]);
      reg->Draw(Form("%s:%s",blumi_name.Data(),bpm_old[ibpm].Data()),cut_text,"COL");

      c810->cd((iblumi+5)*11+ibpm+2);
      TString blumi_name = Form("reg_asym_n_blumi%d-reg_asym_n_blumi%d",blumi_a[iblumi],blumi_b[iblumi]);
      reg->Draw(Form("%s:%s",blumi_name.Data(),bpm_old[ibpm].Data()),cut_text,"COL");

    }
  }

  for(int iblumi =1; iblumi<=4;iblumi++){
    c810->cd(iblumi*11+1);
    TText *t_blumi = new TText(0.2,0.4,Form("%d + %d ",blumi_a[iblumi-1],blumi_b[iblumi-1]));
    t_blumi->SetNDC(kTRUE);
    t_blumi->SetTextSize(0.3);
    t_blumi->Draw();

    c810->cd((iblumi+4)*11+1);
    TText *t_blumi = new TText(0.2,0.4,Form("%d - %d ",blumi_a[iblumi-1],blumi_b[iblumi-1]));
    t_blumi->SetNDC(kTRUE);
    t_blumi->SetTextSize(0.3);
    t_blumi->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    c810->cd(ibpm+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }

  c810->SaveAs(Form("run%d_reg_combined_blumi_bpm_colz.pdf",runnum));


  //////PROF
  for(int ibpm = 0; ibpm<10;ibpm++){
    for(int iblumi = 0; iblumi<4;iblumi++){
      c810->cd((iblumi+1)*11+ibpm+2);
      
      TString histo_name = Form("hcombineblumivsbpmreg1%d%d",ibpm,iblumi);
      TString blumi_name = Form("reg_asym_n_blumi%d+reg_asym_n_blumi%d",blumi_a[iblumi],blumi_b[iblumi]);
      reg->Draw(Form("%s:%s>>%s",blumi_name.Data(),bpm_old[ibpm].Data(),histo_name.Data()),cut_text,"prof");
      TH2D *h2dnew = (TH2D*)gDirectory->FindObject(histo_name);
      FitLinear(h2dnew);


      c810->cd((iblumi+5)*11+ibpm+2);

      TString histo_name = Form("hcombineblumivsbpmraw2%d%d",ibpm,iblumi);
      TString blumi_name = Form("reg_asym_n_blumi%d-reg_asym_n_blumi%d",blumi_a[iblumi],blumi_b[iblumi]);
      reg->Draw(Form("%s:%s>>%s",blumi_name.Data(),bpm_old[ibpm].Data(),histo_name.Data()),cut_text,"prof");
      TH2D *h22dnew = (TH2D*)gDirectory->FindObject(histo_name);
      FitLinear(h22dnew);
    }
  }

  for(int iblumi =1; iblumi<=4;iblumi++){
    c810->cd(iblumi*11+1);
    TText *t_blumi = new TText(0.2,0.4,Form("%d + %d ",blumi_a[iblumi-1],blumi_b[iblumi-1]));
    t_blumi->SetNDC(kTRUE);
    t_blumi->SetTextSize(0.3);
    t_blumi->Draw();

    c810->cd((iblumi+4)*11+1);
    TText *t_blumi = new TText(0.2,0.4,Form("%d - %d ",blumi_a[iblumi-1],blumi_b[iblumi-1]));
    t_blumi->SetNDC(kTRUE);
    t_blumi->SetTextSize(0.3);
    t_blumi->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    c810->cd(ibpm+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }

  c810->SaveAs(Form("run%d_reg_combined_blumi_bpm_prof.pdf",runnum));

  // Plot Correlation between SAMs vs SAMs
  TCanvas *c88 = new TCanvas("c88","c88",900,900);
  c88->Divide(9,9);

  /////COLZ
  for(int iblumi =1; iblumi<=8;iblumi++){
    for(int jblumi = 1; jblumi<=8;jblumi++){
      c88->cd(iblumi*9+jblumi+1);
      reg->Draw(Form("reg_asym_n_blumi%d :reg_asym_n_blumi%d",iblumi,jblumi),cut_text,"COL");
    }
  }
  for(int icanvas = 1; icanvas<=8;icanvas++){
    c88->cd(icanvas+1);
    TText *t1 =new TText(0.4,0.4,Form("%d",icanvas));
    t1->SetTextSize(0.5);
    t1->SetNDC(kTRUE);
    t1->Draw();

    c88->cd(1+(icanvas)*9);
    TText *t2 = new TText(0.5,0.4,Form("%d",icanvas));
    t2->SetTextSize(0.5);
    t2->SetNDC(kTRUE);
    t2->Draw();
  }
  c88->SaveAs(Form("run%d_reg_blumivsblumi_colz.pdf",runnum));

  /////PROFILE
  for(int iblumi =1; iblumi<=8;iblumi++){
    for(int jblumi = 1; jblumi<=8;jblumi++){
      c88->cd(iblumi*9+jblumi+1);
      TString histo_name = Form("hblumivsblumireg%d%d",iblumi,jblumi);
      reg->Draw(Form("reg_asym_n_blumi%d :reg_asym_n_blumi%d>>%s",iblumi,jblumi,histo_name.Data()),cut_text,"prof");
      TH2D *h2d = (TH2D*)gDirectory->FindObject(histo_name);
      double slope = FitLinear(h2d);
      vec_reg_blumi_blumi_slope.push_back(slope);
    }
    vec2d_reg_blumi_blumi_slope.push_back(vec_reg_blumi_blumi_slope);
    vec_reg_blumi_blumi_slope.clear();
  }
  for(int icanvas = 1; icanvas<=8;icanvas++){
    c88->cd(icanvas+1);
    TText *t1 =new TText(0.4,0.4,Form("%d",icanvas));
    t1->SetTextSize(0.5);
    t1->SetNDC(kTRUE);
    t1->Draw();

    c88->cd(1+(icanvas)*9);
    TText *t2 = new TText(0.5,0.4,Form("%d",icanvas));
    t2->SetTextSize(0.5);
    t2->SetNDC(kTRUE);
    t2->Draw();
  }
  c88->SaveAs(Form("run%d_reg_blumivsblumi_prof.pdf",runnum));

  // Plot Normalized Regressed Combined Blumis Correlation
  TCanvas *c4 = new TCanvas("c4","c4",700,700);
  c4->Divide(2,2);
  for(int iplot = 0; iplot <3; iplot++){
    c4->cd(1);
    tree_name = Form("(reg_asym_n_blumi%d+reg_asym_n_blumi%d)*0.5 :(reg_asym_n_blumi4+reg_asym_n_blumi8)*0.5",blumi_a[iplot],blumi_b[iplot]);
    histo_name = Form("h01%d",iplot);
    reg->Draw(Form("%s>>%s",tree_name.Data(),histo_name.Data()),cut_text,"goff");
    
    h1 = (TH2D*)gDirectory->FindObject(histo_name);
    h1->SetTitle("Asym Regressed");
    h1->Draw("COLZ");
  
    c4->cd(2);
    tree_name = Form("(reg_asym_n_blumi%d+reg_asym_n_blumi%d)*0.5",blumi_a[iplot],blumi_b[iplot]);
    histo_name = Form("h02%d",iplot);
    reg->Draw(Form("%s>>%s",tree_name.Data(),histo_name.Data()),cut_text,"goff");
    h2 = (TH1D*)gDirectory->FindObject(histo_name);
    h2->SetTitle(Form("regressed normalized blumi(%d+%d)/2",blumi_a[iplot],blumi_b[iplot]));
    h2->Draw();
    FitGaussian(h2);
  
    c4->cd(3);
    reg->Draw("(reg_asym_n_blumi4+reg_asym_n_blumi8)*0.5 >>h3",cut_text);
    h3->SetTitle("regressed normalized blumi(4+8)/2");
    FitGaussian(h3);

    c4->cd(4);
    tree_name = Form("(reg_asym_n_blumi%d+reg_asym_n_blumi%d)*0.5 -(reg_asym_n_blumi4+reg_asym_n_blumi8)*0.5",blumi_a[iplot],blumi_b[iplot]);
    histo_name = Form("h04%d",iplot);
    reg->Draw(Form("%s>>%s",tree_name.Data(),histo_name.Data()),cut_text,"goff");
    
    h4 = (TH1D*)gDirectory->FindObject(histo_name);
    h4->SetTitle("Double Difference");
    h4->Draw();
    FitGaussian(h4);

    c4->SaveAs(Form("h%d.pdf",iplot));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_reg_blumi_correlation.pdf",runnum));
  gSystem->Exec("rm h*.pdf");
  
  // Blumi a+b and blumi a-b double difference
  int blumi_a[4]={1,2,3,4};
  int blumi_b[4]={5,6,7,8};
  TCanvas *c2 = new TCanvas("c2","c2",1800,600);
  c2->Divide(3,1);
  for(int iplot =0 ;iplot <4;iplot++){
    c2->cd(1);
    histo_name = Form("hcolz%d",iplot);
    tree_name = Form("reg_asym_n_blumi%d:reg_asym_n_blumi%d",blumi_a[iplot],blumi_b[iplot]);
    reg->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text,"goff");
    hcolznew = (TH2D*)gDirectory->FindObject(histo_name);
    hcolznew->SetTitle(tree_name);
    hcolznew->Draw("COLZ");
        
    c2->cd(2);
    histo_name = Form("h11%d",iplot);
    tree_name = Form("reg_asym_n_blumi%d + reg_asym_n_blumi%d",blumi_a[iplot],blumi_b[iplot]);
    reg->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text,"goff");

    h2new =(TH1D*)gDirectory->FindObject(histo_name);
    h2new->SetTitle(tree_name);
    h2new->Draw();
    FitGaussian(h2new);
    
    c2->cd(3);
    histo_name = Form("h12%d",iplot);
    tree_name = Form("reg_asym_n_blumi%d - reg_asym_n_blumi%d",blumi_a[iplot],blumi_b[iplot]);
    reg->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text,"goff");
    
    h3new =(TH1D*)gDirectory->FindObject(histo_name);
    h3new->SetTitle(tree_name);
    h3new->Draw();
    FitGaussian(h3new);

    c2->SaveAs(Form("h%d.pdf",iplot));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_reg_blumi_DD.pdf",runnum));
  gSystem->Exec("rm h*.pdf");


  
  // Plot bcm vs blumi
  
  TCanvas *c59 = new TCanvas("c59","c59",900,500);
  c59->Divide(9,5);
  ///// COLZ 
  for(int ibcm=1;ibcm<=4;ibcm++){
    for(int iblumi=0; iblumi<8;iblumi++){
      c59->cd(ibcm*9+iblumi+2);
      tree_name = Form("reg_asym_n_blumi%d:asym_bcm%d",iblumi+1,ibcm);
      reg->Draw(tree_name.Data(),cut_text,"COLZ");
    }
  }

  for(int ibcm =1; ibcm<=4;ibcm++){
    c59->cd(ibcm*9+1);
    TText *t_bcm = new TText(0.4,0.4,bcm_title[ibcm-1]);
    t_bcm->SetNDC(kTRUE);
    t_bcm->SetTextSize(0.2);
    t_bcm->Draw();
  }
  
  for(int iblumi =0;iblumi<8;iblumi++){
    c59->cd(iblumi+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%d",iblumi+1));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }
 
  c59->SaveAs(Form("run%d_reg_bcm_vs_blumi_colz.pdf",runnum));
  ///// PROFILE

  for(int ibcm=1;ibcm<=4;ibcm++){
    for(int iblumi=0; iblumi<8;iblumi++){
      c59->cd(ibcm*9+iblumi+2);
      tree_name = Form("reg_asym_n_blumi%d:asym_bcm%d",iblumi+1,ibcm);
      reg->Draw(tree_name.Data(),cut_text,"PROF");
    }
  }

  // for(int ibcm =1; ibcm<=4;ibcm++){
  //   c59->cd(ibcm*9+1);
  //   TText *t_bcm = new TText(0.4,0.4,bcm_title[ibcm-1]);
  //   t_bcm->SetNDC(kTRUE);
  //   t_bcm->SetTextSize(0.2);
  //   t_bcm->Draw();
  // }
  
  // for(int iblumi =0;iblumi<8;iblumi++){
  //   c59->cd(iblumi+2);
  //   TText *t_bpm = new TText(0.2,0.2,Form("%d",iblumi+1));
  //   t_bpm->SetNDC(kTRUE);
  //   t_bpm->SetTextSize(0.4);
  //   t_bpm->Draw();
  // }
 
  c59->SaveAs(Form("run%d_reg_bcm_vs_blumi_prof.pdf",runnum));

  /////////////////////////////////////////////////////////////////////////////////////////RAW PLOTS
  /////////////////////////////////////////////////////////////////////////////////////////RAW PLOTS
  /////////////////////////////////////////////////////////////////////////////////////////RAW PLOTS
  /////////////////////////////////////////////////////////////////////////////////////////RAW PLOTS
  TCanvas *c810 = new TCanvas("c810","c810",1100,900);
  c810->Divide(11,9);

  /// Raw blumi vs bpms
  //////COLZ
  for(int iblumi = 1; iblumi<=8;iblumi++){
    for(int ibpm =0 ;ibpm<10;ibpm++){
      c810->cd((iblumi)*11+ibpm+2);
      TString blumi_name = Form("asym_n_blumi%d",iblumi);
      TString tree_name = Form("%s:%s",blumi_name.Data(),bpm_old[ibpm].Data());
      reg->Draw(tree_name.Data(),cut_text,"COL");
    }
  }

  for(int iblumi =1; iblumi<=8;iblumi++){
    c810->cd(iblumi*11+1);
    TText *t_blumi = new TText(0.4,0.4,Form("%d",iblumi));
    t_blumi->SetNDC(kTRUE);
    t_blumi->SetTextSize(0.5);
    t_blumi->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    c810->cd(ibpm+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }
 
  c810->SaveAs(Form("run%d_raw_blumi_bpm_colz.pdf",runnum));

  ///// PROFILE
    for(int iblumi = 1; iblumi<=8;iblumi++){
    for(int ibpm =0 ;ibpm<10;ibpm++){
      c810->cd((iblumi)*11+ibpm+2);
      TString histo_name = Form("hblumivsbpmraw%d%d",iblumi,ibpm);
      TString blumi_name = Form("asym_n_blumi%d",iblumi);
      TString tree_name = Form("%s:%s>>%s",blumi_name.Data(),bpm_old[ibpm].Data(),histo_name.Data());
      reg->Draw(tree_name.Data(),cut_text,"prof");
      TH2D *hfit = (TH2D*)gDirectory->FindObject(histo_name);
      double slope = FitLinear(hfit);
      vec_raw_blumi_bpm_slope.push_back(slope);
    }
    vec2d_raw_blumi_bpm_slope.push_back(vec_raw_blumi_bpm_slope);
    vec_raw_blumi_bpm_slope.clear();
  }

  for(int iblumi =1; iblumi<=8;iblumi++){
    c810->cd(iblumi*11+1);
    TText *t_blumi = new TText(0.4,0.4,Form("%d",iblumi));
    t_blumi->SetNDC(kTRUE);
    t_blumi->SetTextSize(0.5);
    t_blumi->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    c810->cd(ibpm+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }
 
  c810->SaveAs(Form("run%d_raw_blumi_bpm_prof.pdf",runnum));


  // Plot Correlation between combined SAMs with Position Difference
  
  TCanvas *c810 = new TCanvas("c810","c810",1100,900);
  c810->Divide(11,9);

  //////COLZ
  for(int ibpm = 0; ibpm<10;ibpm++){
    for(int iblumi = 0; iblumi<4;iblumi++){
      c810->cd((iblumi+1)*11+ibpm+2);
      TString blumi_name = Form("asym_n_blumi%d+asym_n_blumi%d",blumi_a[iblumi],blumi_b[iblumi]);
      reg->Draw(Form("%s:%s",blumi_name.Data(),bpm_old[ibpm].Data()),cut_text,"COL");

      c810->cd((iblumi+5)*11+ibpm+2);
      TString blumi_name = Form("asym_n_blumi%d-asym_n_blumi%d",blumi_a[iblumi],blumi_b[iblumi]);
      reg->Draw(Form("%s:%s",blumi_name.Data(),bpm_old[ibpm].Data()),cut_text,"COL");

    }
  }

  for(int iblumi =1; iblumi<=4;iblumi++){
    c810->cd(iblumi*11+1);
    TText *t_blumi = new TText(0.2,0.4,Form("%d + %d ",blumi_a[iblumi-1],blumi_b[iblumi-1]));
    t_blumi->SetNDC(kTRUE);
    t_blumi->SetTextSize(0.3);
    t_blumi->Draw();

    c810->cd((iblumi+4)*11+1);
    TText *t_blumi = new TText(0.2,0.4,Form("%d - %d ",blumi_a[iblumi-1],blumi_b[iblumi-1]));
    t_blumi->SetNDC(kTRUE);
    t_blumi->SetTextSize(0.3);
    t_blumi->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    c810->cd(ibpm+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }

  c810->SaveAs(Form("run%d_raw_combined_blumi_bpm_colz.pdf",runnum));


  //////PROF
    for(int ibpm = 0; ibpm<10;ibpm++){
    for(int iblumi = 0; iblumi<4;iblumi++){
      c810->cd((iblumi+1)*11+ibpm+2);
      TString histo_name = Form("hcombineblumivsbpmraw1%d%d",ibpm,iblumi);
      TString blumi_name = Form("asym_n_blumi%d+asym_n_blumi%d",blumi_a[iblumi],blumi_b[iblumi]);
      reg->Draw(Form("%s:%s>>%s",blumi_name.Data(),bpm_old[ibpm].Data(),histo_name.Data()),cut_text,"prof");
      TH2D *h2dnew = (TH2D*)gDirectory->FindObject(histo_name);
      FitLinear(h2dnew);
      
      c810->cd((iblumi+5)*11+ibpm+2);
      TString histo_name = Form("hcombineblumivsbpmraw2%d%d",ibpm,iblumi);
      TString blumi_name = Form("asym_n_blumi%d-asym_n_blumi%d",blumi_a[iblumi],blumi_b[iblumi]);
      reg->Draw(Form("%s:%s>>histo_name",blumi_name.Data(),bpm_old[ibpm].Data(),histo_name.Data()),cut_text,"prof");

      TH2D *h22dnew =(TH2D*)gDirectory->FindObject(histo_name);
      FitLinear(h22dnew);

    }
  }

  for(int iblumi =1; iblumi<=4;iblumi++){
    c810->cd(iblumi*11+1);
    TText *t_blumi = new TText(0.2,0.4,Form("%d + %d ",blumi_a[iblumi-1],blumi_b[iblumi-1]));
    t_blumi->SetNDC(kTRUE);
    t_blumi->SetTextSize(0.3);
    t_blumi->Draw();

    c810->cd((iblumi+4)*11+1);
    TText *t_blumi = new TText(0.2,0.4,Form("%d - %d ",blumi_a[iblumi-1],blumi_b[iblumi-1]));
    t_blumi->SetNDC(kTRUE);
    t_blumi->SetTextSize(0.3);
    t_blumi->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    c810->cd(ibpm+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }

  c810->SaveAs(Form("run%d_raw_combined_blumi_bpm_prof.pdf",runnum));

  // Plot Correlation between SAMs vs SAMs
  TCanvas *c88 = new TCanvas("c88","c88",900,900);
  c88->Divide(9,9);

  /////COLZ
  for(int iblumi =1; iblumi<=8;iblumi++){
    for(int jblumi = 1; jblumi<=8;jblumi++){
      c88->cd(iblumi*9+jblumi+1);
      reg->Draw(Form("asym_n_blumi%d :asym_n_blumi%d",iblumi,jblumi),cut_text,"COL");
    }
  }
  for(int icanvas = 1; icanvas<=8;icanvas++){
    c88->cd(icanvas+1);
    TText *t1 =new TText(0.4,0.4,Form("%d",icanvas));
    t1->SetTextSize(0.5);
    t1->SetNDC(kTRUE);
    t1->Draw();

    c88->cd(1+(icanvas)*9);
    TText *t2 = new TText(0.5,0.4,Form("%d",icanvas));
    t2->SetTextSize(0.5);
    t2->SetNDC(kTRUE);
    t2->Draw();
  }
  c88->SaveAs(Form("run%d_raw_blumivsblumi_colz.pdf",runnum));

  /////PROFILE
  for(int iblumi =1; iblumi<=8;iblumi++){
    for(int jblumi = 1; jblumi<=8;jblumi++){
      c88->cd(iblumi*9+jblumi+1);
      histo_name = Form("hblumiblumiraw%d%d",iblumi,jblumi);
      reg->Draw(Form("asym_n_blumi%d :asym_n_blumi%d>>%s",iblumi,jblumi,histo_name.Data()),cut_text,"prof");
      h2d = (TH2D*)gDirectory->FindObject(histo_name);
      double slope = FitLinear(h2d);
      vec_raw_blumi_blumi_slope.push_back(slope);
    }
    vec2d_raw_blumi_blumi_slope.push_back(vec_raw_blumi_blumi_slope);
    vec_raw_blumi_blumi_slope.clear();
  }
  for(int icanvas = 1; icanvas<=8;icanvas++){
    c88->cd(icanvas+1);
    TText *t1 =new TText(0.4,0.4,Form("%d",icanvas));
    t1->SetTextSize(0.5);
    t1->SetNDC(kTRUE);
    t1->Draw();

    c88->cd(1+(icanvas)*9);
    TText *t2 = new TText(0.5,0.4,Form("%d",icanvas));
    t2->SetTextSize(0.5);
    t2->SetNDC(kTRUE);
    t2->Draw();
  }
  c88->SaveAs(Form("run%d_raw_blumivsblumi_prof.pdf",runnum));

  // Plot Normalized Raw Blumi Correlation
  TCanvas *c4 = new TCanvas("c4","c4",700,700);
  c4->Divide(2,2);
  for(int iplot = 0; iplot <3; iplot++){
    c4->cd(1);
    tree_name = Form("(asym_n_blumi%d+asym_n_blumi%d)*0.5 :(asym_n_blumi4+asym_n_blumi8)*0.5",blumi_a[iplot],blumi_b[iplot]);
    histo_name = Form("h01%d",iplot);
    reg->Draw(Form("%s>>%s",tree_name.Data(),histo_name.Data()),cut_text,"goff");
    
    h1 = (TH2D*)gDirectory->FindObject(histo_name);
    h1->SetTitle("Asym UnRegressed");
    h1->Draw("COLZ");
  
    c4->cd(2);
    tree_name = Form("(asym_n_blumi%d+asym_n_blumi%d)*0.5",blumi_a[iplot],blumi_b[iplot]);
    histo_name = Form("h02%d",iplot);
    reg->Draw(Form("%s>>%s",tree_name.Data(),histo_name.Data()),cut_text,"goff");
    h2 = (TH1D*)gDirectory->FindObject(histo_name);
    h2->SetTitle(Form("Unregressed normalized blumi(%d+%d)/2",blumi_a[iplot],blumi_b[iplot]));
    h2->Draw();
    FitGaussian(h2);
  
    c4->cd(3);
    reg->Draw("(asym_n_blumi4+asym_n_blumi8)*0.5 >>h3",cut_text);
    h3->SetTitle("Unregressed normalized blumi(4+8)/2");
    FitGaussian(h3);

    c4->cd(4);
    tree_name = Form("(asym_n_blumi%d+asym_n_blumi%d)*0.5 -(asym_n_blumi4+asym_n_blumi8)*0.5",blumi_a[iplot],blumi_b[iplot]);
    histo_name = Form("h04%d",iplot);
    reg->Draw(Form("%s>>%s",tree_name.Data(),histo_name.Data()),cut_text,"goff");
    
    h4 = (TH1D*)gDirectory->FindObject(histo_name);
    h4->SetTitle("Double Difference");
    h4->Draw();
    FitGaussian(h4);

    c4->SaveAs(Form("h%d.pdf",iplot));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_raw_blumi_correlation.pdf",runnum));
  gSystem->Exec("rm h*.pdf");
  
  // Blumi a+b and blumi a-b double difference
  int blumi_a[4]={1,2,3,4};
  int blumi_b[4]={5,6,7,8};
  TCanvas *c2 = new TCanvas("c2","c2",1800,600);
  c2->Divide(3,1);
  for(int iplot =0 ;iplot <4;iplot++){
    c2->cd(1);
    histo_name = Form("hcolz%d",iplot);
    tree_name = Form("asym_n_blumi%d:asym_n_blumi%d",blumi_a[iplot],blumi_b[iplot]);
    reg->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text,"goff");
    hcolznew = (TH2D*)gDirectory->FindObject(histo_name);
    hcolznew->SetTitle(tree_name);
    hcolznew->Draw("COLZ");
        
    c2->cd(2);
    histo_name = Form("h11%d",iplot);
    tree_name = Form("asym_n_blumi%d + asym_n_blumi%d",blumi_a[iplot],blumi_b[iplot]);
    reg->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text,"goff");

    h2new =(TH1D*)gDirectory->FindObject(histo_name);
    h2new->SetTitle(tree_name);
    h2new->Draw();
    FitGaussian(h2new);
    
    c2->cd(3);
    histo_name = Form("h12%d",iplot);
    tree_name = Form("asym_n_blumi%d - asym_n_blumi%d",blumi_a[iplot],blumi_b[iplot]);
    reg->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text,"goff");
    
    h3new =(TH1D*)gDirectory->FindObject(histo_name);
    h3new->SetTitle(tree_name);
    h3new->Draw();
    FitGaussian(h3new);

    c2->SaveAs(Form("h%d.pdf",iplot));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_raw_blumi_DD.pdf",runnum));
  gSystem->Exec("rm h*.pdf");


  
  // Plot bcm vs blumi
  
  TCanvas *c59 = new TCanvas("c59","c59",900,500);
  c59->Divide(9,5);
  ///// COLZ 
  for(int ibcm=1;ibcm<=4;ibcm++){
    for(int iblumi=0; iblumi<8;iblumi++){
      c59->cd(ibcm*9+iblumi+2);
      tree_name = Form("asym_n_blumi%d:asym_bcm%d",iblumi+1,ibcm);
      reg->Draw(tree_name.Data(),cut_text,"COLZ");
    }
  }

  for(int ibcm =1; ibcm<=4;ibcm++){
    c59->cd(ibcm*9+1);
    TText *t_bcm = new TText(0.4,0.4,bcm_title[ibcm-1]);
    t_bcm->SetNDC(kTRUE);
    t_bcm->SetTextSize(0.2);
    t_bcm->Draw();
  }
  
  for(int iblumi =0;iblumi<8;iblumi++){
    c59->cd(iblumi+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%d",iblumi+1));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }
 
  c59->SaveAs(Form("run%d_raw_bcm_vs_blumi_colz.pdf",runnum));
  
  ///// PROFILE
  for(int ibcm=1;ibcm<=4;ibcm++){
    for(int iblumi=0; iblumi<8;iblumi++){
      c59->cd(ibcm*9+iblumi+2);
      tree_name = Form("asym_n_blumi%d:asym_bcm%d",iblumi+1,ibcm);
      reg->Draw(tree_name.Data(),cut_text,"prof");
    }
  }

  for(int ibcm =1; ibcm<=4;ibcm++){
    c59->cd(ibcm*9+1);
    TText *t_bcm = new TText(0.4,0.4,bcm_title[ibcm-1]);
    t_bcm->SetNDC(kTRUE);
    t_bcm->SetTextSize(0.2);
    t_bcm->Draw();
  }
  
  for(int iblumi =0;iblumi<8;iblumi++){
    c59->cd(iblumi+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%d",iblumi+1));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }
 
  c59->SaveAs(Form("run%d_raw_bcm_vs_blumi_prof.pdf",runnum));

    //Plot Correlation between bpms 
  TCanvas *c1111 = new TCanvas("c1111","c1111",1100,1100);
  c1111->Divide(11,11);

  for( int ibpm = 0; ibpm<10; ibpm++){
    for(int jbpm =0; jbpm<10; jbpm++){
      c1111->cd((ibpm+1)*11+jbpm+2);
      reg->Draw(Form("%s:%s",bpm_old[ibpm].Data(),bpm_old[jbpm].Data()),cut_text,"COL");
    }
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    c1111->cd(ibpm+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    c1111->cd((ibpm+1)*11+1);
    TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }

  
  c1111->SaveAs(Form("run%d_raw_bpmvsbpm_colz.pdf",runnum));

  //Plot Correlation between bpms PROFILE
  TCanvas *c1111 = new TCanvas("c1111","c1111",1100,1100);
  c1111->Divide(11,11);

  for( int ibpm = 0; ibpm<10; ibpm++){
    for(int jbpm =0; jbpm<10; jbpm++){
      c1111->cd((ibpm+1)*11+jbpm+2);
      TString histo_name = Form("hprofbpm%dbpm%d",ibpm,jbpm);
      reg->Draw(Form("%s:%s>>%s",bpm_old[ibpm].Data(),bpm_old[jbpm].Data(),histo_name.Data()),cut_text,"PROF");
      TH2D *h2dfit = (TH2D*)gDirectory->FindObject(histo_name);
      double slope = FitLinear(h2dfit);
      vec_bpm_bpm_slope.push_back(slope);
    }
    vec2d_bpm_bpm_slope.push_back(vec_bpm_bpm_slope);
    vec_bpm_bpm_slope.clear();
  }
  
  for(int ibpm =0;ibpm<10;ibpm++){
    c1111->cd(ibpm+2);
    TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    c1111->cd((ibpm+1)*11+1);
    TText *t_bpm = new TText(0.2,0.2,Form("%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.4);
    t_bpm->Draw();
  }
  c1111->SaveAs(Form("run%d_raw_bpmvsbpm_prof.pdf",runnum));
  

  //////////////////////////////////// 1D profile
  //////////////////////////////////// 1D profile
  //////////////////////////////////// 1D profile

  TCanvas *c1 = new TCanvas("c1","c1",800,800);

  for(int iblumi = 0;iblumi<8;iblumi++){
    c1->cd();
    TString histo_name = Form("hblumiregdfit%d",iblumi);
    TString tree_name = Form("reg_asym_n_blumi%d>>%s",iblumi+1,histo_name.Data());
    reg->Draw(tree_name,cut_text);
    TH1D *hfit1dblumi =(TH1D*)gDirectory->FindObject(histo_name);
    double rms = hfit1dblumi->GetRMS();
    vec_reg_blumi_rms.push_back(rms);
    
    FitGaussian(hfit1dblumi);
    c1->SaveAs(Form("h%d.pdf",iblumi));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_reg_blumi_1D.pdf",runnum));
  gSystem->Exec("rm h*.pdf");

  for(int iblumi = 0;iblumi<8;iblumi++){
    c1->cd();
    TString histo_name = Form("hblumirawdfit%d",iblumi);
    TString tree_name = Form("asym_n_blumi%d>>%s",iblumi+1,histo_name.Data());
    reg->Draw(tree_name,cut_text);
    TH1D *hfit1dblumi =(TH1D*)gDirectory->FindObject(histo_name);
    
    double rms = hfit1dblumi->GetRMS();
    vec_raw_blumi_rms.push_back(rms);

    FitGaussian(hfit1dblumi);
    c1->SaveAs(Form("h%d.pdf",iblumi));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_raw_blumi_1D.pdf",runnum));
  gSystem->Exec("rm h*.pdf");

  for(int ibcm = 0;ibcm<4;ibcm++){
    c1->cd();
    TString histo_name = Form("hbcm1dfit%d",ibcm+1);
    TString tree_name = Form("asym_bcm%d>>%s",ibcm+1,histo_name.Data());
    reg->Draw(tree_name,cut_text);
    TH1D *hfit1dblumi =(TH1D*)gDirectory->FindObject(histo_name);
    FitGaussian(hfit1dblumi);
    c1->SaveAs(Form("h%d.pdf",ibcm));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_bcms_1D.pdf",runnum));
  gSystem->Exec("rm h*.pdf");


  for(int ibpm=0;ibpm<10;ibpm++){
    c1->cd();
    TString histo_name = Form("hbpmfit%d",ibpm);
    TString tree_name = Form("%s_old>>%s",bpm_name[ibpm].Data(),histo_name.Data());
    reg->Draw(tree_name,cut_text);
    TH1D *hfit1dbpm = (TH1D*)gDirectory->FindObject(histo_name);
    
    double rms = hfit1dbpm->GetRMS();
    vec_bpm_rms.push_back(rms);
    
    FitGaussian(hfit1dbpm);
    c1->SaveAs(Form("h%d.pdf",ibpm));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_bpms_1D.pdf",runnum));
  gSystem->Exec("rm h*.pdf");


  //PRINT tables
  FILE* res_table;
  res_table = fopen(Form("run%d_table.txt",runnum),"w");

  int size_bpm = 10;
  int size_blumi = 8;
  
  //////// reg_blumi vs bpm
  for(int ibpm=0;ibpm<size_bpm;ibpm++){
    if(ibpm==0)
      fprintf(res_table," ");
    fprintf(res_table,"\& %s",bpm_title[ibpm].Data());
  }
  fprintf(res_table," \\\\ \n \\hline \n");

  for(int iblumi =0; iblumi<size_blumi;iblumi++){
    fprintf(res_table,"reg\\_blumi%d",iblumi+1);
    for(int ibpm =0; ibpm<size_bpm; ibpm++){
      fprintf(res_table,"\& %1.3lf",vec2d_reg_blumi_bpm_slope[iblumi][ibpm]);
    }
    fprintf(res_table," \\\\ \n");
  }
  fprintf(res_table,"\n");
  
  //////// raw_blumi vs bpm
  for(int ibpm=0;ibpm<size_bpm;ibpm++){
    if(ibpm==0)
      fprintf(res_table," ");
    fprintf(res_table,"\& %s",bpm_title[ibpm].Data());
  }
  fprintf(res_table," \\\\ \n \\hline \n");

  for(int iblumi =0; iblumi<size_blumi;iblumi++){
    fprintf(res_table,"raw\\_blumi%d",iblumi+1);
    for(int ibpm =0; ibpm<size_bpm; ibpm++){
      fprintf(res_table,"\& %1.3lf",vec2d_raw_blumi_bpm_slope[iblumi][ibpm]);
    }
    fprintf(res_table," \\\\ \n");
  }
  fprintf(res_table,"\n");

  ///// reg_blumi vs reg_blumi
  for(int iblumi=0;iblumi<size_blumi;iblumi++){
    if(iblumi==0)
      fprintf(res_table," ");
    fprintf(res_table,"\& reg\\_blumi%d",iblumi+1);
  }
  fprintf(res_table," \\\\ \n \\hline \n");

  for(int iblumi=0; iblumi<size_blumi;iblumi++){
    fprintf(res_table,"reg\\_blumi%d",iblumi+1);
    for(int jblumi=0;jblumi<size_blumi;jblumi++){
      fprintf(res_table,"\& %1.3lf",vec2d_reg_blumi_blumi_slope[iblumi][jblumi]);
    }
    fprintf(res_table," \\\\ \n");
  }
  fprintf(res_table,"\n");
    
  ///// raw_blumi vs raw_blumi
  for(int iblumi=0;iblumi<size_blumi;iblumi++){
    if(iblumi==0)
      fprintf(res_table," ");
    fprintf(res_table,"\& raw\\_blumi%d",iblumi+1);
  }
  fprintf(res_table," \\\\ \n \\hline \n");

  for(int iblumi=0; iblumi<size_blumi;iblumi++){
    fprintf(res_table,"raw\\_blumi%d",iblumi+1);
    for(int jblumi=0;jblumi<size_blumi;jblumi++){
      fprintf(res_table,"\& %1.3lf",vec2d_raw_blumi_blumi_slope[iblumi][jblumi]);
    }
    fprintf(res_table," \\\\ \n");
  }
  fprintf(res_table,"\n");


  ///// bpm vs bpm
  for(int ibpm=0;ibpm<size_bpm;ibpm++){
    if(ibpm==0)
      fprintf(res_table," ");
    fprintf(res_table,"\& %s",bpm_title[ibpm].Data());
  }
  fprintf(res_table," \\\\ \n \\hline \n");

  for(int ibpm =0; ibpm<size_bpm;ibpm++){
    fprintf(res_table,"%s",bpm_title[ibpm].Data());
    for(int jbpm =0;jbpm<size_bpm;jbpm++){
      fprintf(res_table,"\& %1.3lf",vec2d_bpm_bpm_slope[ibpm][jbpm]);
    }
    fprintf(res_table," \\\\ \n");
  }
  fprintf(res_table,"\n");

  //// bpm rms
  for(int ibpm=0;ibpm<size_bpm;ibpm++){
    if(ibpm==0)
      fprintf(res_table," ");
    fprintf(res_table,"\& %s",bpm_title[ibpm].Data());
  }
  fprintf(res_table," \\\\ \n \\hline \n");
  
  for(int ibpm =0;ibpm<size_bpm;ibpm++){
    if(ibpm ==0)
      fprintf(res_table,"rms");
    fprintf(res_table,"\& %1.2lf",vec_bpm_rms[ibpm]);
  }
  fprintf(res_table,"\\\\ \n");
  fprintf(res_table,"\n");
  
  //// reg_blumi rms
  for(int iblumi=0;iblumi<size_blumi;iblumi++){
    if(iblumi==0)
      fprintf(res_table," ");
    fprintf(res_table,"\& reg\\_blumi%d",iblumi+1);
  }
  fprintf(res_table," \\\\ \n \\hline \n");

  for(int iblumi =0;iblumi<size_blumi;iblumi++){
    if(iblumi==0)
      fprintf(res_table,"rms");
    fprintf(res_table,"\& %1.2lf",vec_reg_blumi_rms[iblumi]);
  }
  fprintf(res_table,"\\\\ \n");
  fprintf(res_table,"\n");

  //// raw_blumi rms
  for(int iblumi=0;iblumi<size_blumi;iblumi++){
    if(iblumi==0)
      fprintf(res_table," ");
    fprintf(res_table,"\& raw\\_blumi%d",iblumi+1);
  }
  fprintf(res_table," \\\\ \n \\hline \n");

  for(int iblumi =0;iblumi<size_blumi;iblumi++){
    if(iblumi==0)
      fprintf(res_table,"rms");
    fprintf(res_table,"\& %1.2lf",vec_raw_blumi_rms[iblumi]);
  }
  fprintf(res_table,"\\\\ \n");
  fprintf(res_table,"\n");

  
  fclose(res_table);
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
