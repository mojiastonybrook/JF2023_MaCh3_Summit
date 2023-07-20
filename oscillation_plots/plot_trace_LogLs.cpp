void plot_trace_LogLs(int start, int number){
  const int chain_num = number;

  //std::string chains_dir = "/gpfs/alpine/proj-shared/phy171/Asimov_fit_chains";
  std::string chains_dir = "/gpfs/alpine/proj-shared/phy171/Data_fit_chains"; 
  std::string chains[chain_num];
  for(int i =0; i < chain_num; i++){
    int c_id = i + start;
    chains[i] = chains_dir+Form("/MaCh3_MCMC_chain_%i_Iter_0_7_reduced.root",c_id);
    std::cout << chains[i] << std::endl;
  }
  
  TH1D **traceHists;
  traceHists = new TH1D*[chain_num];
  // loop over chains 
  for(int i=0; i < chain_num; i++){
    
    TFile* file = new TFile(chains[i].c_str(),"OPEN");
    TChain* chain = new TChain("osc_posteriors","");
    chain->Add(chains[i].c_str());

    int nEntries = chain->GetEntries();
    std::cout << chains[i] << ": " << nEntries << " steps in total. " << std::endl;

    traceHists[i] = new TH1D(Form("h_%i",i),Form("LogL of Chain_%i",i+start), nEntries, 0, nEntries);
    traceHists[i]->GetXaxis()->SetTitle("step");
    traceHists[i]->GetYaxis()->SetTitle("-logL");

    double LogL_value=-999.99;
    // loop over entries
    for(int j=0; j < nEntries; j++){
      chain->SetBranchAddress("LogL",&LogL_value);
      chain->GetEntry(j);
      //fill histogram
      traceHists[i]->SetBinContent(j,LogL_value);
      //if(LogL_value > 10000){
        //std::cout << "step " << j << " " << LogL_value << std::endl;
      //}
    } 
    std::cout << "done with " << chains[i] << std::endl; 
  }

  // plotting
  TCanvas *c1 = new TCanvas("c1"," ",0,0, 800,630);
  gStyle->SetOptStat(0);
  std::string canvas_name = Form("datafit_trace_LogL_chain_%i.eps",start);
  //std::string canvas_name = "trace_LogL.eps";
//  c1->Print(canvas_name+"[");

  for(int i=0; i < chain_num; i++){
    if(i==0){
      traceHists[i]->SetMaximum(5500);
      traceHists[i]->SetTitle("LogL traces");
      traceHists[i]->Draw();

      continue;
    }
    traceHists[i]->Draw("same");

  }
  c1->Print(canvas_name.c_str());
//  c1->Print(canvas_name+"]");

//  delete[] traceHists;
  std::cout<<"DONE!"<<std::endl;
}
