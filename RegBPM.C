void RegBPM(int runnum){
  
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(1);

  TString bpm_name[10] = {"_bpm4bx","_bpm4by","_bpm4ax","_bpm4ay","_bpm14x","_bpm14y","_bpm12x","_bpm12y","_bpm8x","_bpm8y"};
  TString bpm_title[10] = {"4bx","4by","4ax","4ay","14x","14y","12x","12y","8x","8y"};
  TString bpm_prefix = "diff";
  FILE *bpm_conffile;
  TString postpan_dir = "/home/yetao/workarea/jlab_2016/postpan/";

  vector< vector<double> > vec2d_reg_reg_slope;
  vector< vector<double> > vec2d_reg_raw_slope;
  vector< vector<double> > vec2d_raw_raw_slope;
  vector< vector<double> > vec2d_bpm_coeff;
    
  vector< double> vec_reg_reg_slope;
  vector< double> vec_reg_raw_slope;
  vector< double> vec_raw_raw_slope;
  vector< double> vec_bpm_coeff;
  
  vector< double> vec_raw_bpm_rms;
  vector< double> vec_reg_bpm_rms;

  for(int ireg =0; ireg <10; ireg++){
    
    //Prepare  configuration file
    gSystem->Exec(Form("cp %scontrol.conf_reg%d_bpm %scontrol.conf_reg%d%s",postpan_dir.Data(),runnum,postpan_dir.Data(),runnum,bpm_name[ireg].Data()));
    bpm_conffile = fopen(Form("%scontrol.conf_reg%d%s",postpan_dir.Data(),runnum,bpm_name[ireg].Data()),"a");

    for(int ibpm =0; ibpm<10;ibpm++){
      if(ibpm == ireg)
	fprintf(bpm_conffile,"dv diff%s\n",bpm_name[ibpm].Data());
      else
	fprintf(bpm_conffile,"iv diff%s\n",bpm_name[ibpm].Data());
    }
    fclose(bpm_conffile);

    // run postpan
    TString commandline = Form("%sredana -r %d -C %scontrol.conf_reg%d%s",postpan_dir.Data(),runnum,postpan_dir.Data(),runnum,bpm_name[ireg].Data());
    gSystem->Exec(commandline);
    gSystem->Exec(Form("mv ROOTfiles/parity16_%d_regress.root ROOTfiles/parity16_%d_regress%s.root",runnum,runnum,bpm_name[ireg].Data()));

  }
  
  // Extract regressed bpm and Combine to parity16_%runnum_regress_bpm.root
  TFile *combine_file = TFile::Open(Form("./ROOTfiles/parity16_%d_regress_bpm.root",runnum),"RECREATE");
  TTree *new_tree = new TTree("reg","combined_reg");
  TBranch *reg_bpm_branch[10];
  TBranch *raw_bpm_branch[10];
  TBranch *cut_branch;
  double reg_bpm_val[10];
  double raw_bpm_val[10];
  int ok_cut;
  //Print Slope in Latex format: title line
  
  

  for(int ibpm = 0;ibpm<10;ibpm++){
    TString raw_name = "diff"+bpm_name[ibpm];
    TString reg_name = "reg_"+raw_name;
    reg_bpm_branch[ibpm] = new_tree->Branch(reg_name,&reg_bpm_val[ibpm],Form("%s/D",reg_name.Data()));
    raw_bpm_branch[ibpm] = new_tree->Branch(raw_name,&raw_bpm_val[ibpm],Form("%s/D",raw_name.Data()));
  }
  cut_branch = new_tree->Branch("ok_cut",&ok_cut,"ok_cut/I");
    
    
  for(int ibpm = 0; ibpm<10;ibpm++){
    TFile *reg_file = TFile::Open(Form("./ROOTfiles/parity16_%d_regress%s.root",runnum,bpm_name[ibpm].Data()));
    TTree *reg_tree = reg_file->Get("reg");
    int nentries = reg_tree->GetEntries();
    if (ibpm ==0)
      new_tree->SetEntries(nentries);

    for(int ientry =0; ientry<nentries;ientry++){
      reg_tree->GetEntry(ientry);
      TString branch_name_temp = Form("reg_diff%s",bpm_name[ibpm].Data());
      reg_bpm_val[ibpm] = reg_tree->GetLeaf(branch_name_temp)->GetValue(0);
      reg_bpm_branch[ibpm]->Fill();
      
      TString branch_name_temp = Form("diff%s",bpm_name[9-ibpm].Data());
      raw_bpm_val[9-ibpm] = reg_tree->GetLeaf(branch_name_temp)->GetValue(0);
      raw_bpm_branch[9-ibpm]->Fill();

      if(ibpm ==0){
	ok_cut = reg->GetLeaf("ok_cut")->GetValue(0);
	cut_branch->Fill();
      }
    }
    
    // push back Regression Slope
    TTree *coeff_tree = reg_file->Get("regcoeffs");
    int nminirun = coeff_tree->GetEntries();  // number of array, every array has 9 components
    for(int ientry = 0; ientry<nminirun;ientry++){
      coeff_tree->GetEntry(ientry);
      for(int i =0; i<9;i++){
	double coeffs = coeff_tree->GetLeaf("coeff")->GetValue(i); // an array
	vec_bpm_coeff.push_back(coeffs);
      }
    }
    vec2d_bpm_coeff.push_back(vec_bpm_coeff);
    vec_bpm_coeff.clear();
    reg_file->Close();
  }
  combine_file->Write();
  
  //Make Plots
  TCanvas *c1111_col = new TCanvas("c1111_col","c1111_col",1100,1100);
  c1111_col->Divide(11,11);
  TCanvas *c1111_prof = new TCanvas("c1111_prof","c1111_prof",1100,1100);
  c1111_prof->Divide(11,11);

  // REG VS REG
  for(int ibpm = 0; ibpm<10; ibpm++){
    for(int jbpm=0; jbpm<10;jbpm++){
      
      c1111_prof->cd((ibpm+1)*11+jbpm+2);
      new_tree->Draw(Form("reg_diff%s:reg_diff%s>>h2prof%d%d",bpm_name[ibpm].Data(),bpm_name[jbpm].Data(),ibpm,jbpm),"ok_cut","PROF");
      TH2D *h2dfit = (TH2D*)gDirectory->FindObject(Form("h2prof%d%d",ibpm,jbpm));
      double slope = FitLinear(h2dfit);
      vec_reg_reg_slope.push_back(slope);
      
      c1111_col->cd((ibpm+1)*11+jbpm+2);
      new_tree->Draw(Form("reg_diff%s:reg_diff%s>>h2col%d%d",bpm_name[ibpm].Data(),bpm_name[jbpm].Data(),ibpm,jbpm),"ok_cut","COL");
    }
    vec2d_reg_reg_slope.push_back(vec_reg_reg_slope);
    vec_reg_reg_slope.clear();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    TText *t_bpm = new TText(0.2,0.2,Form("reg_%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.2);
    c1111_prof->cd(ibpm+2);
    t_bpm->Draw();
    c1111_col->cd(ibpm+2);
    t_bpm->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    TText *t_bpm = new TText(0.2,0.4,Form("reg_%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.2);
    c1111_prof->cd((ibpm+1)*11+1);
    t_bpm->Draw();
    c1111_col->cd((ibpm+1)*11+1);
    t_bpm->Draw();
  }
  c1111_col->SaveAs(Form("run%d_bpm_reg_vs_reg_colz.pdf",runnum));
  c1111_prof->SaveAs(Form("run%d_bpm_reg_vs_reg_prof.pdf",runnum));
  
  // REG VS RAW
  TCanvas *c1111_col = new TCanvas("c1111_col","c1111_col",1100,1100);
  c1111_col->Divide(11,11);
  TCanvas *c1111_prof = new TCanvas("c1111_prof","c1111_prof",1100,1100);
  c1111_prof->Divide(11,11);

  for(int ibpm = 0; ibpm<10; ibpm++){
    for(int jbpm=0; jbpm<10;jbpm++){

      c1111_prof->cd((ibpm+1)*11+jbpm+2);
      new_tree->Draw(Form("reg_diff%s:diff%s>>h2regrawprof%d%d",bpm_name[ibpm].Data(),bpm_name[jbpm].Data(),ibpm,jbpm),"ok_cut","PROF");
      TH2D *h2dfit = (TH2D*)gDirectory->FindObject(Form("h2regrawprof%d%d",ibpm,jbpm));
      double slope = FitLinear(h2dfit);
      vec_reg_raw_slope.push_back(slope);
      
      c1111_col->cd((ibpm+1)*11+jbpm+2);
      new_tree->Draw(Form("reg_diff%s:diff%s>>h2regrawcol%d%d",bpm_name[ibpm].Data(),bpm_name[jbpm].Data(),ibpm,jbpm),"ok_cut","COL");
    }
    vec2d_reg_raw_slope.push_back(vec_reg_raw_slope);
    vec_reg_raw_slope.clear();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    TText *t_bpm = new TText(0.2,0.2,Form("raw_%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.2);
    c1111_prof->cd(ibpm+2);
    t_bpm->Draw();
    c1111_col->cd(ibpm+2);
    t_bpm->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    TText *t_bpm = new TText(0.2,0.4,Form("reg_%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.2);
    c1111_prof->cd((ibpm+1)*11+1);
    t_bpm->Draw();
    c1111_col->cd((ibpm+1)*11+1);
    t_bpm->Draw();
  }
  c1111_col->SaveAs(Form("run%d_bpm_reg_vs_raw_colz.pdf",runnum));
  c1111_prof->SaveAs(Form("run%d_bpm_reg_vs_raw_prof.pdf",runnum));


  
  // RAW VS RAW
  TCanvas *c1111_col = new TCanvas("c1111_col","c1111_col",1100,1100);
  c1111_col->Divide(11,11);
  TCanvas *c1111_prof = new TCanvas("c1111_prof","c1111_prof",1100,1100);
  c1111_prof->Divide(11,11);

  for(int ibpm = 0; ibpm<10; ibpm++){
    for(int jbpm=0; jbpm<10;jbpm++){
      c1111_prof->cd((ibpm+1)*11+jbpm+2);
      new_tree->Draw(Form("diff%s:diff%s>>h2rawrawprof%d%d",bpm_name[ibpm].Data(),bpm_name[jbpm].Data(),ibpm,jbpm),"ok_cut","PROF");
      TH2D *h2dfit = (TH2D*)gDirectory->FindObject(Form("h2rawrawprof%d%d",ibpm,jbpm));
      double slope = FitLinear(h2dfit);
      vec_raw_raw_slope.push_back(slope);
      
      c1111_col->cd((ibpm+1)*11+jbpm+2);
      new_tree->Draw(Form("diff%s:diff%s>>h2rawrawcol%d%d",bpm_name[ibpm].Data(),bpm_name[jbpm].Data(),ibpm,jbpm),"ok_cut","COL");
    }
    vec2d_raw_raw_slope.push_back(vec_raw_raw_slope);
    vec_raw_raw_slope.clear();
  }
  
  for(int ibpm =0;ibpm<10;ibpm++){
    TText *t_bpm = new TText(0.2,0.2,Form("raw_%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.2);
    c1111_prof->cd(ibpm+2);
    t_bpm->Draw();
    c1111_col->cd(ibpm+2);
    t_bpm->Draw();
  }
  for(int ibpm =0;ibpm<10;ibpm++){
    TText *t_bpm = new TText(0.2,0.4,Form("raw_%s",bpm_title[ibpm].Data()));
    t_bpm->SetNDC(kTRUE);
    t_bpm->SetTextSize(0.2);
    c1111_prof->cd((ibpm+1)*11+1);
    t_bpm->Draw();
    c1111_col->cd((ibpm+1)*11+1);
    t_bpm->Draw();
  }
  
  c1111_col->SaveAs(Form("run%d_bpm_raw_vs_raw_colz.pdf",runnum));
  c1111_prof->SaveAs(Form("run%d_bpm_raw_vs_raw_prof.pdf",runnum));

  // 1D plot
  gStyle->SetOptFit(1);
  
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  
  for(int ibpm =0; ibpm<10; ibpm++){
    c1->cd();
    new_tree->Draw(Form("reg_diff%s>>hregbpm%d",bpm_name[ibpm].Data(),ibpm),"ok_cut","goff");
    TH1D *hfit = (TH1D*)gDirectory->FindObject(Form("hregbpm%d",ibpm));
    hfit->Draw();
    FitGaussian(hfit);
    double rms = hfit->GetRMS();
    vec_reg_bpm_rms.push_back(rms);
    c1->SaveAs(Form("h%d.pdf",ibpm));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_reg_diff_bpm_1d.pdf",runnum));
  gSystem->Exec("rm h*.pdf");

  for(int ibpm =0; ibpm<10; ibpm++){
    c1->cd();
    new_tree->Draw(Form("diff%s>>hrawbpm%d",bpm_name[ibpm].Data(),ibpm),"ok_cut","goff");
    TH1D *hfit = (TH1D*)gDirectory->FindObject(Form("hrawbpm%d",ibpm));
    hfit->Draw();
    FitGaussian(hfit);
    double rms = hfit->GetRMS();
    vec_raw_bpm_rms.push_back(rms);
    c1->SaveAs(Form("h%d.pdf",ibpm));
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_raw_diff_bpm_1d.pdf",runnum));
  gSystem->Exec("rm h*.pdf");


  //Print table
  FILE *res_text;
  res_text = fopen(Form("run%d_bpmreg_res.txt",runnum),"w");

  //Print slopes
  fprintf(res_text," ");
  for(int ibpm =0; ibpm<10;ibpm++){
    fprintf(res_text,"\& reg\\_%s ",bpm_title[ibpm].Data());
  }
  fprintf(res_text,"\\\\ \n");
  fprintf(res_text," \\hline \n");
  for(int ibpm = 0; ibpm<10; ibpm++){
    for(int jbpm =0; jbpm<10; jbpm++){
      if(jbpm ==0)
	fprintf(res_text,"reg\\_%s",bpm_title[ibpm].Data());
      fprintf(res_text,"\& %1.4lf",vec2d_reg_reg_slope[ibpm][jbpm]);
    }
    fprintf(res_text,"\\\\ \n");
  }
  fprintf(res_text,"\n");

  fprintf(res_text," ");
  for(int ibpm =0; ibpm<10;ibpm++){
    fprintf(res_text,"\& raw\\_%s ",bpm_title[ibpm].Data());
  }
  fprintf(res_text,"\\\\ \n");
  fprintf(res_text," \\hline \n");
  for(int ibpm = 0; ibpm<10; ibpm++){
    for(int jbpm =0; jbpm<10; jbpm++){
      if(jbpm ==0)
	fprintf(res_text,"reg\\_%s",bpm_title[ibpm].Data());
      fprintf(res_text,"\& %1.4lf",vec2d_reg_raw_slope[ibpm][jbpm]);
    }
    fprintf(res_text,"\\\\ \n");
  }
  fprintf(res_text,"\n");
  
  fprintf(res_text," ");
  for(int ibpm =0; ibpm<10;ibpm++){
    fprintf(res_text,"\& raw\\_%s ",bpm_title[ibpm].Data());
  }
  fprintf(res_text,"\\\\ \n");
  fprintf(res_text," \\hline \n");
  for(int ibpm = 0; ibpm<10; ibpm++){
    for(int jbpm =0; jbpm<10; jbpm++){
      if(jbpm ==0)
	fprintf(res_text,"raw\\_%s",bpm_title[ibpm].Data());
      fprintf(res_text,"\& %1.4lf",vec2d_raw_raw_slope[ibpm][jbpm]);
    }
    fprintf(res_text,"\\\\ \n");
  }
  fprintf(res_text,"\n");
  
  // Print Reg Coeff
  for(int iminirun =0 ; iminirun<nminirun;iminirun++){
    fprintf(res_text,"minirun-%d",iminirun);
    for(int ibpm =0; ibpm<10;ibpm++){
      fprintf(res_text,"\& %s ",bpm_title[ibpm].Data());
    }
    fprintf(res_text,"\\\\ \n");
    fprintf(res_text," \\hline \n");

    for(int ibpm =0; ibpm<10;ibpm++){
      for(int j=0;j<9;j++){
	if(j==0)
	  fprintf(res_text,"%s",bpm_title[ibpm].Data());
	if(j == ibpm)
 	  fprintf(res_text,"\& ");
	fprintf(res_text,"\& %1.4lf",vec2d_bpm_coeff[ibpm+9*iminirun][j]);
      }
      fprintf(res_text,"\\\\ \n");
    }
    fprintf(res_text,"\n");
  }

  //Print RMS 
  for(int ibpm =0; ibpm<10;ibpm++){
    fprintf(res_text,"\& %s ",bpm_title[ibpm].Data());
  }
  fprintf(res_text,"\\\\ \n");
  fprintf(res_text," \\hline \n");

  fprintf(res_text," raw rms");
  for(int ibpm = 0; ibpm<10;ibpm++){
    fprintf(res_text,"\& %1.2lf",vec_raw_bpm_rms[ibpm]);
  }
  fprintf(res_text,"\\\\ \n");

  fprintf(res_text," reg rms");
  for(int ibpm = 0; ibpm<10;ibpm++){
    fprintf(res_text,"\& %1.2lf",vec_reg_bpm_rms[ibpm]);
  }
  fprintf(res_text,"\\\\ \n");

  combine_file->Close();
  fclose(res_text);
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
