void plot_corrfddet_diag(std::string filename){
    // Open the MCMC file
  TFile *f = new TFile(filename.c_str());
  TTree *t = (TTree*)f->Get("posteriors");

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetLeftMargin(canv->GetLeftMargin()*1.5);

  std::string sub_str = filename.substr(filename.find("output/"));
  TString pdfname = sub_str;
  pdfname.ReplaceAll("output/","");
  pdfname.ReplaceAll(".root", "");
  pdfname+="_corrfddet_hist.pdf";
  canv->Print(pdfname+"["); 
  // Loop over the branches
  int nbranches = t->GetListOfBranches()->GetSize();
  for (int i=0;i < nbranches; i++){
    std::string branchname = t->GetListOfBranches()->At(i)->GetName();
    if (branchname.find("corrfddet_") == std::string::npos) continue;
    if (branchname[0] != std::string("c")) continue;
    //if (branchname.find("_prop")==std::string::npos) continue;
    std::cout << branchname << std::endl;

    std::string drawcmd = branchname+Form(">>h%i", i);
    t->Draw(drawcmd.c_str());
    canv->Print(pdfname);
  } 
  canv->Print(pdfname+"]");
}
