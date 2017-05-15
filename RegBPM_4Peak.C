// Global Variables

TString postpan_dir = "/home/yetao/workarea/jlab_2016/postpan/";
TString peak_name[4]={"101","001","110","010"};
TString bpm_title[10] = {"4bx","4by","4ax","4ay","14x","14y","12x","12y","8x","8y"};

void RegBPM_4Peak(int runnum){

  for(int ipeak =0; ipeak<4;ipeak++){
    for(int ibpm =0; ibpm<10; ibpm++){
      // WriteConfig(runnum,ibpm,ipeak);
      // RegBPM(runnum,ibpm,ipeak);
    }
    //CombineRegBPM(runnum,ipeak);
  }

  //Then make some other plots and print results
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

