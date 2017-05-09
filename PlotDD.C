void PlotDD(int runnum){
  gROOT->SetBatch(kTRUE);
  TFile *rootfile = TFile::Open(Form("./ROOTfiles/parity16_%d_standard.root",runnum));
  gStyle->SetOptFit(1);
  // Set Alias for SAMs and BCMs
  // Ab<n>_<m> : SAM-<n> - <m> combination
  // Aq<n>_<m> : BCM<n> -<m> combination
  double sam_cut = 5e3;
  double bcm_cut = 5e3;
  double sam_cut_asym = 3e3;
  double bcm_cut_asym = 3e3;
  int comb_first[6]={1,2,3,1,2,1};
  int comb_second[6]={2,3,4,3,4,4};
  for(int isam =0; isam < 8 ; isam++){
    for(int icomb =0 ; icomb <6 ;icomb++){
      R->SetAlias(Form("Ab%d_%d",isam+1,icomb),Form("1e6*(vqwk5_%d_b%d_cal-vqwk5_%d_b%d_cal)/(vqwk5_%d_b%d_cal+vqwk5_%d_b%d_cal)",isam,comb_first[icomb],isam,comb_second[icomb],isam,comb_first[icomb],isam,comb_second[icomb]));
    }
    R->SetAlias(Form("Ab%d_%d",isam+1,6),Form("1e6*(vqwk5_%d_b1_cal +vqwk5_%d_b2_cal - vqwk5_%d_b3_cal - vqwk5_%d_b4_cal)/(vqwk5_%d_b1_cal + vqwk5_%d_b2_cal+ vqwk5_%d_b3_cal + vqwk5_%d_b4_cal)",isam,isam,isam,isam,isam,isam,isam,isam));
    P->SetAlias(Form("Ab%d_%d",isam+1,7),Form("asym_blumi%d",isam+1));
  }
  for(int ibcm=0; ibcm<4; ibcm++){  // ibcm = vqwk index
    // Aq_<n> : <n> bcm number
    for(int icomb =0; icomb<6;icomb++){
      R->SetAlias(Form("Aq%d_%d",ibcm+1,icomb),Form("1e6*(vqwk0_%d_b%d_cal-vqwk0_%d_b%d_cal)/(vqwk0_%d_b%d_cal+vqwk0_%d_b%d_cal)",ibcm,comb_first[icomb],ibcm,comb_second[icomb],ibcm,comb_first[icomb],ibcm,comb_second[icomb]));
    }
    R->SetAlias(Form("Aq%d_%d",ibcm+1,6),Form("1e6*(vqwk0_%d_b1_cal +vqwk0_%d_b2_cal - vqwk0_%d_b3_cal - vqwk0_%d_b4_cal)/(vqwk0_%d_b1_cal + vqwk0_%d_b2_cal+ vqwk0_%d_b3_cal + vqwk0_%d_b4_cal)",ibcm,ibcm,ibcm,ibcm,ibcm,ibcm,ibcm,ibcm));
    P->SetAlias(Form("Aq%d_%d",ibcm+1,7),Form("asym_bcm%d",ibcm+1));
  }
  
  // Draw Plots
  TCanvas *c1 = new TCanvas("c1","c1",700,700);
  c1->Divide(2,2);
  int sam_a[2]={3,4};
  int sam_b[2]={7,8};

  TString histo_title[8]={"(1-2)/(1+2)","(2-3)/(2+3)","(3-4)/(3+4)","(1-3)/(1+3)","(2-4)/(2+4)","(1-4)/(1+4)","(1+2-3-4)/(1+2+3+4)","Asym"};

  TString x_title[4] = {"Digital Upstream BCM","Digital Downstream BCM","Analog Upstream BCM","Analog Downstream BCM"};

  TTree *tree[8] ={R,R,R,R,R,R,R,P};

  // run2721 cut
  TString cut_raw ="(ev_num > 10 && ev_num<5e3)||(ev_num>20e3 && ev_num< 25e3)" ;
  TString cut_parity = "(m_ev_num > 10 && m_ev_num < 5e3) ||(m_ev_num>20e3 && m_ev_num < 25e3)";

  //run2856 cut
  // TString cut_raw ="(ev_num > 35e3 && ev_num<70e3)" ;
  // TString cut_parity = "(m_ev_num > 35e3 && m_ev_num < 70e3)";

    //run2725 cut
  // TString cut_raw ="(ev_num > 2e4 && ev_num<6e4)" ;
  // TString cut_parity = "(m_ev_num > 2e4 && m_ev_num < 6e4)";
  
  //run2851 cut
  // TString cut_raw ="(ev_num > 10 && ev_num<103)||(ev_num>22e3 && ev_num< 44e3)" ;
  // TString cut_parity = "(m_ev_num > 10 && m_ev_num < 10e3) ||(m_ev_num>22e3 && m_ev_num < 44e3)";
  //run2347 cut
  // TString cut_raw ="(ev_num > 30e3 && ev_num<40e3)||(ev_num>70e3 && ev_num< 80e3)||(ev_num > 105e3 && ev_num<130e3)||(ev_num>145e3 && ev_num< 160e3)||(ev_num > 175e3 && ev_num<190e3)";
  // TString cut_parity = "(m_ev_num > 30e3 && m_ev_num < 40e3) ||(m_ev_num>70e3 && m_ev_num < 80e3) ||(m_ev_num>105e3 && m_ev_num < 130e3) ||(m_ev_num>145e3 && m_ev_num < 160e3) ||(m_ev_num>175e3 && m_ev_num < 190e3)";

  
  TString cut_text[8];
  for(int icomb=0; icomb<7;icomb++){
    cut_text[icomb] = cut_raw;
  }
  cut_text[7] = cut_parity;
  
  TCanvas *c_all = new TCanvas("c_all","c_all",1200,600);
  TH2D *h[8];
  c_all ->Divide(4,2);
  for(int icomb =0; icomb<8; icomb++){
    // cut_raw = Form("ev_num>10 && abs(Aq3_%d)< %f && abs(Aq4_%d)<%f && bcm%d > 500",icomb,bcm_cut,icomb,bcm_cut,ibcm+1);
    c_all->cd(icomb+1);
    tree[icomb]->Draw(Form("Aq3_%d : Aq4_%d>>hist%d",icomb,icomb,icomb),cut_text[icomb],"goff");
    h[icomb] = (TH2D*)gDirectory->FindObject(Form("hist%d",icomb));
    h[icomb]->SetTitle(histo_title[icomb]);
    h[icomb]->GetXaxis()->SetTitle("Analog Downstream BCM ");
    h[icomb]->GetYaxis()->SetTitleOffset(1.5);
    h[icomb]->GetYaxis()->SetTitle("Analog Upstream BCM");
    h[icomb]->Draw("COLZ");
  }
  c_all ->SaveAs(Form("run%d_call.pdf",runnum));

  TCanvas *c6 = new TCanvas("c6","c6",1200,800);
  c6->Divide(3,2);
  for(int icomb =0; icomb<3; icomb++){
    // cut_raw = Form("ev_num>10 && abs(Aq3_%d)< %f && abs(Aq4_%d)<%f && bcm%d > 500",icomb,bcm_cut,icomb,bcm_cut,ibcm+1);
    c6->cd(icomb+1);
    R->Draw(Form("0.25*(Ab3_%d+Ab7_%d+Ab4_%d+Ab8_%d):Aq3_%d>>hist",icomb,icomb,icomb,icomb,icomb),cut_raw,"COLZ");
    h[icomb] = (TH2D*)gDirectory->FindObject("hist");
    h[icomb]->SetTitle(histo_title[icomb]);
    h[icomb]->GetXaxis()->SetTitle("Analog Upstream BCM ");
    h[icomb]->GetYaxis()->SetTitleOffset(1.5);
    h[icomb]->GetYaxis()->SetTitle("SAMs");
    h[icomb]->Draw("COLZ");

    c6->cd(icomb+4);
    R->Draw(Form("0.25*(Ab3_%d+Ab7_%d+Ab4_%d+Ab8_%d):Aq4_%d>>hist",icomb,icomb,icomb,icomb,icomb),cut_raw,"COLZ");
    h[icomb+3] = (TH2D*)gDirectory->FindObject("hist");
    h[icomb+3]->SetTitle(histo_title[icomb]);
    h[icomb+3]->GetXaxis()->SetTitle("Analog Downstream BCM ");
    h[icomb+3]->GetYaxis()->SetTitleOffset(1.5);
    h[icomb+3]->GetYaxis()->SetTitle("SAMs");
    h[icomb+3]->Draw("COLZ");
  }

  c6 ->SaveAs(Form("run%d_blumivsbcm.pdf",runnum));

  // blumi 37 or 48 vs bcm1,2,3,4
  for( int ibcm = 1 ;ibcm < 4; ibcm ++){  // ibcm =1  =>  to skip digital upstream bcm
    for(int iplot = 0; iplot <2; iplot++){
      for(int icomb =0; icomb<8; icomb++){
	//cut_raw = Form("ev_num>10 && abs(0.5*(Ab%d_%d+ Ab%d_%d))<%f && abs(Aq%d_%d)<%f && bcm%d>500",sam_a[iplot],icomb,sam_b[iplot],icomb,sam_cut,ibcm+1,icomb,bcm_cut,ibcm+1);
	c1->cd(1);
	TString tree_name = Form("0.5*(Ab%d_%d+Ab%d_%d):Aq%d_%d",sam_a[iplot],icomb,sam_b[iplot],icomb,ibcm+1,icomb,icomb);
	TString histo_name = Form("h1colz%d",icomb);
	tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb],"goff");

	hcolz = (TH2D*)gDirectory->FindObject(Form("%s",histo_name.Data()));
	hcolz->SetTitle(histo_title[icomb]);
	hcolz->GetXaxis()->SetTitle(x_title[ibcm]);
	hcolz->GetYaxis()->SetTitle(Form("Asym(SAM %d + SAM %d)/2",sam_a[iplot],sam_b[iplot]));
	hcolz->Draw("COLZ");
	
	
	c1->cd(2);
	tree_name = Form("0.5*(Ab%d_%d+Ab%d_%d)",sam_a[iplot],icomb,sam_b[iplot],icomb);
	histo_name = Form("h11%d",icomb);
	tree[icomb]->Draw(Form("%s >>%s",tree_name.Data(),histo_name.Data()),cut_text[icomb],"goff");

	h1 = (TH1D*)gDirectory->FindObject(Form("%s",histo_name.Data()));
	h1->SetTitle(Form("Asym(SAM %d + SAM %d)/2",sam_a[iplot],sam_b[iplot]));
	h1->GetXaxis()->SetTitle(Form("Asym(SAM %d + SAM %d)/2",sam_a[iplot],sam_b[iplot]));
	h1->Draw();
		
	c1->cd(3);
	tree_name = Form("Aq%d_%d",ibcm+1,icomb);
	histo_name = Form("h12%d",icomb);
	tree[icomb]->Draw(Form("%s>> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb],"goff");
	h2 = (TH1D*)gDirectory->FindObject(Form("%s",histo_name.Data()));
	h2->SetTitle(x_title[ibcm]);
	h2->GetXaxis()->SetTitle(x_title[ibcm]);
	h2->Draw();
	
	c1->cd(4);
	tree_name = Form("0.5*(Ab%d_%d+Ab%d_%d)-Aq%d_%d",sam_a[iplot],icomb,sam_b[iplot],icomb,ibcm+1,icomb);
	histo_name = Form("h13%d",icomb);
	tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb],"goff");
	h3 = (TH1D*)gDirectory->FindObject(Form("%s",histo_name.Data()));
	h3->SetTitle("Double Difference");
	h3->GetXaxis()->SetTitle("Double Difference");
	h3->Draw();

	FitGaussian(h1);
	FitGaussian(h2);
	FitGaussian(h3);
	
	c1->SaveAs(Form("h%d.pdf",icomb));
      }
      gSystem->Exec(Form("pdfunite h*.pdf run%d_blumi%d%d_bcm%d.pdf",runnum,sam_a[iplot],sam_b[iplot],ibcm+1));
      gSystem->Exec(Form("rm h*.pdf"));
    }
  }
  // blumi (3+7)-(4+8)
  
  for(int ibcm =1 ; ibcm < 4; ibcm++){
    for(int icomb =0; icomb<8; icomb++){
      //cut_raw = Form("ev_num>10 && abs(Ab3_%d-Ab4_%d+Ab7_%d-Ab8_%d)<%f && abs(Aq%d_%d)<%f && bcm%d > 500 ",icomb,icomb,icomb,icomb,sam_cut_asym*4,ibcm+1,icomb,bcm_cut,ibcm+1);
      c1->cd(1);
      tree_name = Form("0.25*(Ab3_%d+Ab7_%d-Ab4_%d-Ab8_%d):Aq%d_%d",icomb,icomb,icomb,icomb,ibcm+1,icomb);
      histo_name =Form("h2colz%d",icomb);
      tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb],"goff");
      
      hcolz = (TH2D*)gDirectory->FindObject(Form("%s",histo_name.Data()));
      hcolz->SetTitle(histo_title[icomb]);
      hcolz->GetXaxis()->SetTitle(x_title[ibcm]);
      hcolz->GetYaxis()->SetTitle("Asym_blumi(3+7-4-8)/4");
      hcolz->Draw("COLZ");

      c1->cd(2);
      tree_name = Form("0.25*(Ab3_%d+Ab7_%d-Ab4_%d-Ab8_%d)",icomb,icomb,icomb,icomb);
      histo_name = Form("h21%d",icomb);
      tree[icomb]->Draw(Form("%s>>%s",tree_name.Data(),histo_name.Data()),cut_text[icomb]);
      
      h1 = (TH1D*)gDirectory->FindObject(Form("%s",histo_name.Data()));
      h1->SetTitle("Asym_blumi(3+7-4-8)/4");
      h1->GetXaxis()->SetTitle("Asym_blumi(3+7-4-8)/4");
      h1->Draw();
      
      c1->cd(3);
      tree_name = Form("Aq%d_%d",ibcm+1,icomb);
      histo_name = Form("h22%d",icomb);
      tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb]);
      
      h2 = (TH1D*)gDirectory->FindObject(Form("%s",histo_name.Data()));
      h2->SetTitle(x_title[ibcm]);
      h2->GetXaxis()->SetTitle(x_title[ibcm]);
      h2->Draw();
			       
      c1->cd(4);
      tree_name = Form("0.25*(Ab3_%d+Ab7_%d-Ab4_%d-Ab8_%d)-Aq%d_%d",icomb,icomb,icomb,icomb,ibcm+1,icomb);
      histo_name = Form("h23%d",icomb);
      tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb]);
      
      h3 = (TH1D*)gDirectory->FindObject(Form("%s",histo_name.Data()));
      h3->SetTitle("Double Difference");
      h3->GetXaxis()->SetTitle("Double Difference");
      h3->Draw();

      FitGaussian(h1);
      FitGaussian(h2);
      FitGaussian(h3);
       
      c1->SaveAs(Form("h%d.pdf",icomb));
    }
    // Save and pdf unite
    gSystem->Exec(Form("pdfunite h*.pdf run%d_blumi37-48_bcm%d.pdf",runnum,ibcm+1));
    gSystem->Exec(Form("rm h*.pdf"));

    // TCanvas *c01= new TCanvas("c01","c01",500,500);
    // c01->cd();
    // P->Draw(Form("0.25*(Ab3_%d+Ab7_%d-Ab4_%d-Ab8_%d)>> h3",icomb,icomb,icomb,icomb),cut_parity);
    // h3->SetTitle("3+7-4-8");
    // h3->GetXaxis()->SetTitle("Double Difference");
    // FitGaussian(h3);
    // c01->SaveAs("3+7-4-8.pdf");
  }
  icomb =7;

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  c2->cd(1);  
  //cut_parity = Form("m_ev_num>10 && abs(Ab3_%d+Ab4_%d+Ab7_%d+Ab8_%d)<%f && abs(Aq%d_%d)<%f && avg_bcm%d > 500 ",icomb,icomb,icomb,icomb,sam_cut_asym*3,4,icomb,2500.0,4);
  P->Draw(Form("0.25*(Ab3_%d+Ab7_%d+Ab4_%d+Ab8_%d):Aq4_%d>>hcore1",icomb,icomb,icomb,icomb,icomb),cut_parity,"COLZ");
  hcore1->SetTitle("3+7+4+8");
  c2->cd(2);
  //cut_parity = Form("m_ev_num>10 && abs(Ab3_%d+Ab4_%d+Ab7_%d+Ab8_%d)<%f && abs(Aq%d_%d)<%f && avg_bcm%d > 500 ",icomb,icomb,icomb,icomb,sam_cut_asym*2,4,icomb,2500.0,4);
  P->Draw(Form("0.25*(Ab3_%d+Ab7_%d-Ab4_%d-Ab8_%d):Aq4_%d>>hcore2",icomb,icomb,icomb,icomb,icomb),cut_parity,"COLZ");
  hcore2->SetTitle("3+7-4-8");
  c2->SaveAs(Form("run%d_blumibcmcorrelation.pdf",runnum));
  
  // blumi (3+7)+(4+8)
  
  for(int ibcm =1 ; ibcm < 4; ibcm++){
    for(int icomb =0; icomb<8; icomb++){
      //cut_raw = Form("ev_num>10 && abs(Ab3_%d+Ab4_%d+Ab7_%d+Ab8_%d)<%f && abs(Aq%d_%d)<%f && bcm%d > 500 ",icomb,icomb,icomb,icomb,sam_cut*4,ibcm+1,icomb,bcm_cut,ibcm+1);
      c1->cd(1);
      tree_name = Form("0.25*(Ab3_%d+Ab7_%d+Ab4_%d+Ab8_%d):Aq%d_%d",icomb,icomb,icomb,icomb,ibcm+1,icomb);
      histo_name = Form("h3colz%d",icomb);
      tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb],"goff");

      hcolz = (TH2D*)gDirectory->FindObject(Form("%s",histo_name.Data()));
      hcolz->SetTitle(histo_title[icomb]);
      hcolz->GetXaxis()->SetTitle(x_title[ibcm]);
      hcolz->GetYaxis()->SetTitle("Asym_blumi(3+7+4+8)/4");
      hcolz->Draw("COLZ");

      c1->cd(2);
      tree_name = Form("0.25*(Ab3_%d+Ab7_%d+Ab4_%d+Ab8_%d)",icomb,icomb,icomb,icomb);
      histo_name = Form("h31%d",icomb);
      tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb]);

      h1 = (TH1D*)gDirectory->FindObject(Form("%s",histo_name.Data()));
      h1->SetTitle("Asym_blumi(3+7+4+8)/4");
      h1->GetXaxis()->SetTitle("Asym_blumi(3+7+4+8)/4");
      h1->Draw();
      
      c1->cd(3);
      tree_name = Form("Aq%d_%d",ibcm+1,icomb);
      histo_name = Form("h32%d",icomb);
      tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb]);
      
      h2 = (TH1D*)gDirectory->FindObject(Form("%s",histo_name.Data()));
      h2->SetTitle(x_title[ibcm]);
      h2->GetXaxis()->SetTitle(x_title[ibcm]);
      h2->Draw();
			       
      c1->cd(4);
      tree_name = Form("0.25*(Ab3_%d+Ab7_%d+Ab4_%d+Ab8_%d)-Aq%d_%d",icomb,icomb,icomb,icomb,ibcm+1,icomb);
      histo_name = Form("h33%d",icomb);
      tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb]);
      h3 = (TH1D*)gDirectory->FindObject(histo_name);
      h3->SetTitle("Double Difference");
      h3->GetXaxis()->SetTitle("Double Difference");

      FitGaussian(h1);
      FitGaussian(h2);
      FitGaussian(h3);
       
      c1->SaveAs(Form("h%d.pdf",icomb));
    }
    // Save and pdf unite
    gSystem->Exec(Form("pdfunite h*.pdf run%d_blumi37+48_bcm%d.pdf",runnum,ibcm+1));
    gSystem->Exec(Form("rm h*.pdf"));
  }

  
  // plot 37 vs 48
 
  for(int icomb =0; icomb<8; icomb ++){
    c1->cd(1);
    tree_name = Form("0.5*(Ab3_%d+Ab7_%d)-Aq4_%d:0.5*(Ab8_%d+Ab4_%d)-Aq4_%d",icomb,icomb,icomb,icomb,icomb,icomb);
    histo_name = Form("h4colz%d",icomb);
    tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb],"goff");
    
    hcolz = (TH2D*)gDirectory->FindObject(histo_name);
    hcolz->SetTitle(histo_title[icomb]);
    hcolz->GetXaxis()->SetTitle("Normalized blumi(4+8)/2");
    hcolz->GetYaxis()->SetTitle("Normalized blumi(3+7)/2");
    hcolz->Draw("COLZ");
    
    c1->cd(2);
    tree_name = Form("0.5*(Ab3_%d+Ab7_%d)-Aq4_%d",icomb,icomb,icomb);
    histo_name = Form("h41%d",icomb);
    tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb]);
    h1 = (TH1D*)gDirectory->FindObject(histo_name);
    h1->SetTitle("Normalized (3+7)/2");
    h1->Draw();
    
    c1->cd(3);
    tree_name = Form("0.5*(Ab4_%d+Ab8_%d)-Aq4_%d",icomb,icomb,icomb);
    histo_name = Form("h42%d",icomb);
    tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb]);

    h2 = (TH1D*)gDirectory->FindObject(histo_name);
    h2->SetTitle("Normalized (4+8)/2");
    h2->Draw();

    c1->cd(4);
    histo_name = Form("h43%d",icomb);
    tree_name = Form("0.5*(Ab3_%d+Ab7_%d)-0.5*(Ab4_%d+Ab8_%d)",icomb,icomb,icomb,icomb);
    tree[icomb]->Draw(Form("%s>>%s",tree_name.Data(),histo_name.Data()),cut_text[icomb]);    

    h3 = (TH1D*)gDirectory->FindObject(histo_name.Data());
    h3->SetTitle("Double Difference");
    h3->Draw();
    
    FitGaussian(h1);
    FitGaussian(h2);
    FitGaussian(h3);
    c1->SaveAs(Form("h%d.pdf",icomb));
  }

  gSystem->Exec(Form("pdfunite h*.pdf run%d_blumi37vs48.pdf",runnum));
  gSystem->Exec(Form("rm h*.pdf"));
  
  // plot double difference of Downstream - Upstream BCM
  for(int icomb =0; icomb<8; icomb++){
    //cut_raw = Form("ev_num>10 && abs(Aq3_%d)< %f && abs(Aq4_%d)<%f && bcm%d > 500",icomb,bcm_cut,icomb,bcm_cut,ibcm+1);
    c1->cd(1);
    tree_name = Form("Aq3_%d:Aq4_%d",icomb,icomb);
    histo_name = Form("h5colz%d",icomb);
    tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb],"goff");

    hcolz = (TH2D*)gDirectory->FindObject(histo_name);
    hcolz->SetTitle(histo_title[icomb]);
    hcolz->GetXaxis()->SetTitle(x_title[2]);
    hcolz->GetYaxis()->SetTitle(x_title[3]);
    hcolz->Draw("COLZ");

    c1->cd(2);
    tree_name = Form("Aq3_%d",icomb);
    histo_name = Form("h51%d",icomb);
    tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb]);
    h1 = (TH1D*)gDirectory->FindObject(histo_name);
    h1->SetTitle(x_title[2]);
    h1->Draw();
      
    c1->cd(3);
    tree_name = Form("Aq4_%d",icomb);
    histo_name = Form("h52%d",icomb);
    tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb]);
    h2 = (TH1D*)gDirectory->FindObject(histo_name);
    h2->SetTitle(x_title[3]);
    h2->Draw();

    c1->cd(4);
    tree_name = Form("Aq3_%d-Aq4_%d",icomb,icomb);
    histo_name = Form("h53%d",icomb);
    tree[icomb]->Draw(Form("%s >> %s",tree_name.Data(),histo_name.Data()),cut_text[icomb]);
    h3 =(TH1D*)gDirectory->FindObject(histo_name);
    h3->SetTitle("Double Difference");
    h3->Draw();

    FitGaussian(h1);
    FitGaussian(h2);
    FitGaussian(h3);
      
    c1->SaveAs(Form("h%d.pdf",icomb));
  }

  // Save and pdf unite
  gSystem->Exec(Form("pdfunite h*.pdf run%d_bcm_dd.pdf",runnum));
  gSystem->Exec(Form("rm h*.pdf"));
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
