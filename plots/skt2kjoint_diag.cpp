void skt2kjoint_diag(std::string filename) {
  // Open the MCMC file
  TFile *f = new TFile(filename.c_str());
  TTree *t = (TTree*)f->Get("posteriors");

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetLeftMargin(canv->GetLeftMargin()*1.5);

  std::string sub_str = filename.substr(filename.find("output/"));
  TString pdfname = sub_str;
  pdfname.ReplaceAll("output/","");
  pdfname.ReplaceAll(".root", "all.pdf");
  canv->Print(pdfname+"[");

  // Open up another TFile with the priors, lower bound and upper bound
  TFile *fprior = new TFile("xsec_covariance_2020a_IncAtmosphericModel_v15.root");
  TObjArray *xsec_param_names = (TObjArray*)fprior->Get("xsec_param_names");
  TVectorD *xsec_param_prior = (TVectorD*)fprior->Get("xsec_param_prior");
  TVectorD *xsec_param_lb = (TVectorD*)fprior->Get("xsec_param_lb");
  TVectorD *xsec_param_ub = (TVectorD*)fprior->Get("xsec_param_ub");

  // Loop over the branches
  int nbranches = t->GetListOfBranches()->GetSize();
  for (int i = 0; i < nbranches; ++i) {
    std::string branchname = t->GetListOfBranches()->At(i)->GetName();

    // Only look at cross-section parameters
    if (branchname.find("xsec_") == std::string::npos) continue;
    // Make sure they are all "xsec_" by looking at the first character
    if (branchname[0] != std::string("x")) continue;
    std::cout << branchname << std::endl;
    bool isxsec = false;
    // if xsec, get the prior term, upper bound, lower bound
    std::string realname;
    // Find which nth xsec parameter this is by looking at the branch name
    size_t index = branchname.find("xsec_");
    int nxsec = -1;
    if (index != std::string::npos) {
      isxsec = true;
      // extract which xsec parameter this is
      std::string paramno = branchname.substr(5, branchname.size());
      // convert to int
      nxsec = std::atoi(paramno.c_str());
      // Get the pretty name of the parameter from the xsec_covariance...root file
      realname = xsec_param_names->At(nxsec)->GetName();
      if (realname.find("b_") != std::string::npos) continue;
    }

    // The draw command (save histogram into object h1, h2, h3, etc
    std::cout << "branchname: " << branchname << std::endl;
    std::string drawcmd = branchname+Form(">>h%i", i);
    std::cout << "draw: " << drawcmd << std::endl;
    t->Draw(drawcmd.c_str());

    // Make some TLine's to draw over the posteriors for the upper, lower bounds and the prior central value
    TLine *ub = NULL;
    TLine *lb = NULL;
    TLine *prior = NULL;
    bool isbad = false;
    if (isxsec) {
      TH1D *temp = (TH1D*)gDirectory->Get(Form("h%i", i));
      // Check the integral at the upper and lower bound
      if (temp->Integral(temp->FindBin((*xsec_param_ub)[nxsec]), temp->GetXaxis()->GetNbins()) > 0 ||
          temp->Integral(0, temp->FindBin((*xsec_param_lb)[nxsec])) > 0) {
        isbad = true;
      }
      if (!isbad) continue;
      temp->GetXaxis()->SetTitle(realname.c_str());
      // Also draw some lines for the prior, upper bound, lower bound
      ub =    new TLine((*xsec_param_ub)[nxsec], 0, (*xsec_param_ub)[nxsec], temp->GetMaximum());
      lb =    new TLine((*xsec_param_lb)[nxsec], 0, (*xsec_param_lb)[nxsec], temp->GetMaximum());
      prior = new TLine((*xsec_param_prior)[nxsec], 0, (*xsec_param_prior)[nxsec], temp->GetMaximum());
      ub->SetLineColor(kRed);
      lb->SetLineColor(kRed);
      prior->SetLineColor(kBlack);
      prior->SetLineStyle(kDashed);
      ub->SetLineWidth(2);
      lb->SetLineWidth(2);
      prior->SetLineWidth(2);
      if (ub != NULL) {
        ub->Draw("same");
        lb->Draw("same");
        prior->Draw("same");
      }
      canv->SetLogy(true);
    }
    canv->Print(pdfname);

    // Draw the TH2D of the branch name against the step
    std::string drawcmd2 = branchname+Form(":step>>h2_%i", i);
    t->Draw(drawcmd2.c_str());
    if (isxsec) {
      TH2D *h2 = (TH2D*)gDirectory->Get(Form("h2_%i", i));
      h2->GetYaxis()->SetTitle(realname.c_str());
      h2->GetXaxis()->SetTitle("step");

      ub->SetY1(ub->GetX1());
      ub->SetY2(ub->GetX2());
      ub->SetX1(h2->GetXaxis()->GetBinLowEdge(1));
      ub->SetX2(h2->GetXaxis()->GetBinLowEdge(h2->GetXaxis()->GetNbins()+1));

      lb->SetY1(lb->GetX1());
      lb->SetY2(lb->GetX2());
      lb->SetX1(h2->GetXaxis()->GetBinLowEdge(1));
      lb->SetX2(h2->GetXaxis()->GetBinLowEdge(h2->GetXaxis()->GetNbins()+1));

      prior->SetY1(prior->GetX1());
      prior->SetY2(prior->GetX2());
      prior->SetX1(h2->GetXaxis()->GetBinLowEdge(1));
      prior->SetX2(h2->GetXaxis()->GetBinLowEdge(h2->GetXaxis()->GetNbins()+1));

      h2->Draw();
      ub->Draw("same");
      lb->Draw("same");
      prior->Draw("same");
      canv->SetLogy(false);
    }

    canv->Print(pdfname);

    if (isxsec) {
      delete ub;
      delete lb;
      delete prior;
    }
  }
  canv->Print(pdfname+"]");
}
