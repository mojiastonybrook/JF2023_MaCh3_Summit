void plot_autocorrelation_comp(std::string filename1, std::string filename2){
  // open root file
  std::string filenames[2];
  filenames[0]=filename1;
  filenames[1]=filename2;

  TFile* files[2];
  for(int i =0; i<2; i++){
    files[i] = new TFile(filenames[i].c_str());
  }

  TCanvas *c1 = new TCanvas("c1"," ", 0, 0, 800,630);
  gStyle->SetOptStat(0);
  // prepare output name
  
  TString pdfname;
  pdfname = "corrfddet_autocorr_halfchain.pdf";

  c1->Print(pdfname+"[");

  // get to autocorrelation subdirectory in the root file
  files[0]->cd("Auto_corr");

  int nParam = gDirectory->GetListOfKeys()->GetSize();

  TLine *y_0_2 = NULL;
  for(int i =0; i < nParam; i++){
   
    files[0]->cd("Auto_corr");
    std::string plotname = gDirectory->GetListOfKeys()->At(i)->GetName();
    TH1D* temp1 = (TH1D*)gDirectory->Get(plotname.c_str())->Clone(); 
	
    y_0_2 = new TLine(temp1->GetXaxis()->GetBinLowEdge(1), 0.2, temp1->GetXaxis()->GetBinLowEdge(temp1->GetXaxis()->GetNbins()+1),0.2);
    y_0_2->SetLineColor(kRed);
    y_0_2->SetLineWidth(2);
    y_0_2->SetLineStyle(kDashed);
       
    files[1]->cd("Auto_corr");
    std::string plotname = gDirectory->GetListOfKeys()->At(i)->GetName();
    TH1D* temp2 = (TH1D*)gDirectory->Get(plotname.c_str())->Clone();
    temp2->SetLineColor(kBlack);
    double min = temp1->GetMinimum() < temp2->GetMinimum()? temp1->GetMinimum():temp2->GetMinimum();

    temp1->SetMinimum(min*1.3);
     
    temp1->Draw();
    y_0_2->Draw("same");
    temp2->Draw("same");
    
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(temp1,"corrfddet fStepSize = 0.05","l");
    leg->AddEntry(temp2,"corrfddet fStepSize = 0.02","l");
    leg->Draw("same");

    c1->Print(pdfname);
  }

  c1->Print(pdfname+"]");
  
}
