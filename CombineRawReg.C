void CombineRawReg(int runnum){

  TFile *regfile = TFile::Open(Form("./ROOTfiles/parity16_%d_regress.root",runnum));
  TTree *reg_tree = (TTree*)regfile->Get("reg");

  TFile *combinefile = TFile::Open(Form("./ROOTfiles/parity16_%d_combine.root",runnum),"RECREATE");
  TTree *new_tree  = reg_tree->CloneTree();

  TString bpm_name[10] = {"diff_bpm4bx","diff_bpm4by","diff_bpm4ax","diff_bpm4ay","diff_bpm14x","diff_bpm14y","diff_bpm12x","diff_bpm12y","diff_bpm8x","diff_bpm8y"};
  TString suffix_text = "_old";
  double diff_bpm_readout[10];

  TBranch *bpm_branch[10];
  for(int ibpm =0;ibpm<10;ibpm++){
    TString branch_name = bpm_name[ibpm]+suffix_text; 
    bpm_branch[ibpm] = new_tree->Branch(branch_name,&diff_bpm_readout[ibpm],Form("%s/D",branch_name.Data()));
  }

  TFile *panfile = TFile::Open(Form("./ROOTfiles/parity16_%d_standard.root",runnum));
  TTree *pan_tree = (TTree*)panfile->Get("P");
  int nentries = pan_tree->GetEntries();
  cout << nentries << endl;
  for(int ientry =0 ;ientry<nentries;ientry++){
    pan_tree->GetEntry(ientry);
    for(int ibpm =0; ibpm<10;ibpm++){
      diff_bpm_readout[ibpm] = pan_tree->GetLeaf(bpm_name[ibpm])->GetValue(0);
      bpm_branch[ibpm]->Fill();
    }
  }
  cout<< new_tree->GetEntries() << endl;
  combinefile->Write("",TObject::kOverwrite);
  
  panfile->Close();
  regfile->Close();
  combinefile->Close();
}
