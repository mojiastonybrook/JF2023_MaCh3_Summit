void read_diagMCMC(std::string filename,std::string systname){
  //root -l read_diagMCMC.cpp'("filename","systname")'
  //systanme: xsec_,atmfulx_,corrfddet_

  // open root file
  TFile* file = new TFile(filename.c_str(),"open");

  TCanvas *c1 = new TCanvas("c1"," ", 0, 0, 1620,630);
  TPad* pad1 = new TPad("pad1","Autocorrelation",0.,0.,0.5,1.0);
  TPad* pad2 = new TPad("pad2","Trace",0.5,0.,1.0,1.0);
  gStyle->SetOptStat(0);

  // prepare output name
  std::string sub_str = filename.substr(filename.find("output/"));
  TString pdfname = sub_str;
  pdfname.ReplaceAll("output/","");
  pdfname.ReplaceAll(".root","");
  pdfname = pdfname+"_"+systname.c_str()+".pdf";

  c1->Print(pdfname+"[");

  // get to autocorrelation subdirectory in the root file
  TDirectory* auto_dir = file->Get("Auto_corr");
  TDirectory* trace_dir = file->Get("Trace");

  int nParam = auto_dir->GetListOfKeys()->GetSize();
  TLine *y_0_2 = NULL;
  for(int i =0; i < nParam; i++){
    std::string auto_corr_name = auto_dir->GetListOfKeys()->At(i)->GetName();
    // select syst group of interest
    if(auto_corr_name.find(systname)==std::string::npos) continue;
    std::cout << auto_corr_name << std::endl;

    std::string trace_name = trace_dir->GetListOfKeys()->At(i)->GetName();
    std::cout << trace_name << std::endl;

    // pad1 for auto_corr
    c1->cd();
    pad1->Draw();
    pad1->cd();
    TH1D* temp = (TH1D*)auto_dir->Get(auto_corr_name.c_str());

    y_0_2 = new TLine(temp->GetXaxis()->GetBinLowEdge(1), 0.2, temp->GetXaxis()->GetBinLowEdge(temp->GetXaxis()->GetNbins()+1),0.2);
    y_0_2->SetLineColor(kRed);
    y_0_2->SetLineWidth(kDashed);

    temp->Draw();    
    y_0_2->Draw("same");
    // pad2 for trace
    c1->cd();
    pad2->Draw();
    pad2->cd();

    TH1D* temp2 = (TH1D*)trace_dir->Get(trace_name.c_str());
    temp2->Draw();

    //save
    c1->Print(pdfname);
  }

  c1->Print(pdfname+"]");
  
}
