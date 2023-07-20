void read_autocorrelation(std::string filename){
  // open root file
  TFile* file = new TFile(filename.c_str());

  TCanvas *c1 = new TCanvas("c1"," ", 0, 0, 800,630);
  // prepare output name
  std::string sub_str = filename.substr(filename.find("output/"));
  TString pdfname = sub_str;
  pdfname.ReplaceAll("output/","");
  pdfname.ReplaceAll(".root","");
  pdfname += "_autocorr.pdf";

  c1->Print(pdfname+"[");

  // get to autocorrelation subdirectory in the root file
  file->cd("Auto_corr");

  int nParam = gDirectory->GetListOfKeys()->GetSize();
  TLine *y_0_2 = NULL;
  for(int i =0; i < nParam; i++){
    std::string plotname = gDirectory->GetListOfKeys()->At(i)->GetName();
    // select syst group of interest
    if(plotname.find("xsec_")==std::string::npos) continue;
    std::cout << plotname << std::endl;

    TH1D* temp = (TH1D*)gDirectory->Get(plotname.c_str());

    y_0_2 = new TLine(temp->GetXaxis()->GetBinLowEdge(1), 0.2, temp->GetXaxis()->GetBinLowEdge(temp->GetXaxis()->GetNbins()+1),0.2);
    y_0_2->SetLineColor(kRed);
    y_0_2->SetLineWidth(kDashed);

    temp->Draw();    
    y_0_2->Draw("same");

    c1->Print(pdfname);
  }

  c1->Print(pdfname+"]");
  
}
