// Global Variables

TString postpan_dir = "/home/yetao/workarea/jlab_2016/postpan/";
TString peak_name[4]={"101","001","110","010"};
TString bpm_title[10] = {"4bx","4by","4ax","4ay","14x","14y","12x","12y","8x","8y"};

void RegBPM_4Peak(int runnum){
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(1);
  for(int ipeak =0; ipeak<4;ipeak++){
    for(int ibpm =0; ibpm<10; ibpm++){
      // WriteConfig(runnum,ibpm,ipeak);
      // RegBPM(runnum,ibpm,ipeak);
    }
    //CombineRegBPM(runnum,ipeak);
  }
  //Then make some other plots and print results
  PrintRMS(runnum);
  Plot_BPM_Profile(runnum);

}
void Plot_BPM_Profile(int runnum){
  TString histo_name;
  TH1D *hfit;
  TCanvas *c1  = new TCanvas("c1","c1",800,800);
  
  for(int ipeak =0; ipeak<4;ipeak++){
    TFile *reg_file = TFile::Open(Form("./ROOTfiles/parity16_%d_bpmreg_%s.root",
				       runnum,peak_name[ipeak].Data()));
    TTree *reg_tree = reg_file->Get("reg");
    
    for(int ibpm =0; ibpm<10;ibpm++){
      histo_name = Form("hregbpm%d%d",ibpm,ipeak);
      c1->cd();
      reg_tree->Draw(Form("reg_diff_bpm%s>>%s",
			  bpm_title[ibpm].Data(),histo_name.Data()),"ok_cut");
      hfit = (TH1D*)gDirectory->FindObject(histo_name);
      FitGaussian(hfit);
      hfit->SetTitle(Form("reg_diff_bpm%s - %s",bpm_title[ibpm].Data(),peak_name[ipeak].Data()));
      c1->SaveAs(Form("h%d%d.pdf",ipeak,ibpm));
    }
    reg_file->Close();
  }
  gSystem->Exec(Form("pdfunite h*.pdf run%d_bpm_1d_profile_4peak.pdf",runnum));
  gSystem->Exec("rm h*.pdf");
}

void PrintRMS(int runnum){
  FILE *table_txt;
  table_txt = fopen(Form("run%d_table_bpm_rms_4peak.txt",runnum),"w");

  for(int ipeak =0; ipeak<4;ipeak++){
    TFile *reg_file = TFile::Open(Form("./ROOTfiles/parity16_%d_bpmreg_%s.root",
				       runnum,peak_name[ipeak].Data()));
    TTree *reg_tree = reg_file->Get("reg");
    double bpm_rms[4][10];
    TString histo_name;
    TString leaf_name;
    TH1D *hrms;
    for(int ibpm =0; ibpm<10;ibpm++){
      histo_name = Form("hbpm%d",ibpm);
      leaf_name = Form("reg_diff_bpm%s",bpm_title[ibpm].Data());
      reg_tree->Draw(Form("%s>>%s",leaf_name.Data(),histo_name.Data()),"ok_cut","goff");
      hrms = (TH1D*)gDirectory->FindObject(histo_name);
      bpm_rms[ipeak][ibpm] = hrms->GetRMS();
    }
    reg_file->Close();
  }
  // Print table in Latex Format
  for(int ibpm =0; ibpm<10;ibpm++){
    fprintf(table_txt," \& %s",bpm_title[ibpm].Data());
  }
  fprintf(table_txt,"\\\\ \n");
  for(int ipeak =0;ipeak<4;ipeak++){
    fprintf(table_txt,"%s",peak_name[ipeak].Data());
    for(int ibpm =0; ibpm<10; ibpm++){
      fprintf(table_txt,"\& %1.3lf",bpm_rms[ipeak][ibpm]);
    }
    fprintf(table_txt,"\\\\ \n");
  }
  fclose(table_txt);
}

void CombineRegBPM(int runnum, int ipeak){
  TFile *combine_file = TFile::Open(Form("./ROOTfiles/parity16_%d_bpmreg_%s.root",
					 runnum,peak_name[ipeak].Data()),"RECREATE");
  TTree *new_tree = new TTree("reg","reg");

  TBranch *reg_bpm_branch[10];
  TBranch *raw_bpm_branch[10];
  TBranch *cut_branch;
  double reg_bpm_val[10];
  double raw_bpm_val[10];
  int ok_cut;
  
  for(int ibpm =0; ibpm<10; ibpm++){
    TString raw_name = "diff_bpm"+bpm_title[ibpm];
    TString reg_name = "reg_diff_bpm"+bpm_title[ibpm];

    reg_bpm_branch[ibpm] = new_tree->Branch(reg_name,&reg_bpm_val[ibpm],Form("%s/D",reg_name.Data()));
    raw_bpm_branch[ibpm] = new_tree->Branch(raw_name,&raw_bpm_val[ibpm],Form("%s/D",raw_name.Data()));
  }
  cut_branch = new_tree->Branch("ok_cut",&ok_cut,"ok_cut/I");

  TString branch_name_temp;
  for(int ibpm =0; ibpm<10; ibpm++){
    TFile *reg_file = TFile::Open(Form("./ROOTfiles/parity16_%d_regress_%s_%s.root",
				       runnum,bpm_title[ibpm].Data(),peak_name[ipeak].Data()));
    TTree *reg_tree = reg_file->Get("reg");
    int nentries = reg_tree->GetEntries();
    if(ibpm ==0)
      new_tree->SetEntries(nentries);

    for(int ientry =0;ientry<nentries;ientry++){
      reg_tree->GetEntry(ientry);
      branch_name_temp = Form("reg_diff_bpm%s",bpm_title[ibpm].Data());
      reg_bpm_val[ibpm] = reg_tree->GetLeaf(branch_name_temp)->GetValue(0);
      reg_bpm_branch[ibpm]->Fill();

      branch_name_temp = Form("diff_bpm%s",bpm_title[9-ibpm].Data());
      raw_bpm_val[9-ibpm] = reg_tree->GetLeaf(branch_name_temp)->GetValue(0);
      raw_bpm_branch[ibpm]->Fill();

      if(ibpm ==0){
	ok_cut = reg->GetLeaf("ok_cut")->GetValue(0);
	cut_branch->Fill();
      }
    }
    reg_file->Close();
  }
  combine_file->Write();
  combine_file->Close();
}
  

void RegBPM(int runnum,int ibpm,int ipeak){
  TString config_name = Form("%scontrol.conf_reg%d_bpm%s_%s",
			     postpan_dir.Data(),runnum,bpm_title[ibpm].Data(),peak_name[ipeak].Data());
  TString command_line = Form("%sredana -r %d -C %s",
			      postpan_dir.Data(),runnum,config_name.Data());
  
  gSystem->Exec(command_line.Data());
  gSystem->Exec(Form("mv ./ROOTfiles/parity16_%d_regress.root ./ROOTfiles/parity16_%d_regress_%s_%s.root",
		     runnum,
		     runnum,bpm_title[ibpm].Data(),peak_name[ipeak].Data()));
}

void WriteConfig(int runnum,int ibpm,int ipeak){
  TString new_config_file_name = Form("%scontrol.conf_reg%d_bpm%s_%s",
				      postpan_dir.Data(),runnum,bpm_title[ibpm].Data(),peak_name[ipeak].Data());
  gSystem->Exec(Form("cp %scontrol.conf_reg%d_bpm %s",
		     postpan_dir.Data(),runnum,new_config_file_name.Data()));

  FILE *new_config_file;
  new_config_file = fopen(new_config_file_name.Data(),"a"); //append

  TString parity_cut[4] = {"&&evt_pairsynch[0]==0&&prev_hel==1",
			   "&&evt_pairsynch[0]==0&&prev_hel==0",
			   "&&evt_pairsynch[1]==0&&prev_hel==1",
			   "&&evt_pairsynch[1]==0&&prev_hel==0"};
  
  fprintf(new_config_file,"%s\n",parity_cut[ipeak].Data());

  for(int jbpm =0; jbpm<10; jbpm++){
    if(jbpm == ibpm)
      fprintf(new_config_file,"dv diff_bpm%s\n",bpm_title[jbpm].Data());
    else
      fprintf(new_config_file,"iv diff_bpm%s\n",bpm_title[jbpm].Data());
  }

  fclose(new_config_file);
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
