#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <algorithm>
#include <math.h>

#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TBranch.h"
#include "TLine.h"
#include "TString.h"
#include "TStyle.h"

#include "TH1.h"
#include "TVectorD.h"

#include "TCanvas.h"
#include "TColor.h"

#include "TStopwatch.h"

void plot_trace_xsec_cov(std::string file_name){
  
  std::string sub_str = file_name.substr(file_name.find("output/"));
  TString base_name = sub_str;
  base_name.ReplaceAll(".root","");
  base_name.ReplaceAll("output/","");
  std::cout << base_name << std::endl;
   
  TFile* input_file = new TFile(file_name.c_str(),"OPEN");
  TChain* Chain = new TChain("posteriors","");
  Chain->Add(file_name.c_str());

  int nEntries = Chain->GetEntries();
  std::cout << "Total Steps Recorded: " << nEntries << std::endl;

  TObjArray* brlis = (TObjArray*)(Chain->GetListOfBranches());
  int nBranches = brlis->GetEntries();

  std::vector<TString> XsecName_v;
  int nXsec = 0;
  
  for(int i=0;i<nBranches;i++){
    TBranch* br=(TBranch*)brlis->At(i);
    TString brname=br->GetName();

    if(brname.BeginsWith("xsec_")){
      std::cout<< "xsec parameter: "<< brname << std::endl;
      XsecName_v.push_back(brname);
      nXsec++;
    }
  }
  std::cout << "found "<< nXsec << " xsec parameters. " << std::endl;
  std::cout << "initializing arrays... "<< std::endl; 
  double **XsecValues;
  XsecValues = new double*[nEntries]();
  for(int i =0; i<nEntries;i++){
    XsecValues[i] =new double[nXsec]();
    for(int j=0;j<nXsec;j++){
       XsecValues[i][j] = -999.99;
    }
  }
  // read input tree
  std::cout << "Reading input tree... " << std::endl;
  TStopwatch clock;
  clock.Start();
  // set all the branches to off except for the ones of interests
  Chain->SetBranchStatus("*",false);
  for(int i=0; i<nXsec; i++){
    Chain->SetBranchStatus(XsecName_v[i].Data(),true);
  }
  //loop over the entries
  for(int i=0; i<nEntries; i++){
    if (i%100 ==0){
       std::cout << i << "/" << nEntries << " (" << double(i)/double(nEntries)*100. << "%)" << std::endl;
    }
    // set the branch address for sample LLH
    for(int j=0; j<nXsec; j++){
      Chain->SetBranchAddress(XsecName_v[j].Data(), &XsecValues[i][j]);
    }
    // fill the arry
    Chain->GetEntry(i);
  }
  clock.Stop();
  std::cout << "Took " << clock.RealTime() << "s to finish for " << nEntries << " events" << std::endl;
  double *maxs = new double[nXsec]();
  double *mins = new double[nXsec]();
  for(int i=0;i<nXsec;i++){
    maxs[i]=Chain->GetMaximum(XsecName_v[i].Data());
    mins[i]=Chain->GetMinimum(XsecName_v[i].Data());
    std:: cout << "upper limit: " << maxs[i] << "; lower limit: " << mins[i] << std::endl;
  }
  delete Chain;
  //making plots
  TH1D **TraceXsecPlots;
  TH1D **HistXsecPlots;
  TraceXsecPlots = new TH1D*[nXsec];
  HistXsecPlots =new TH1D*[nXsec];
  //initialize histograms
  for(int i =0; i<nXsec;i++){
    TraceXsecPlots[i] = new TH1D(XsecName_v[i].Data(),XsecName_v[i].Data(),nEntries,0,nEntries);
    TraceXsecPlots[i]->GetXaxis()->SetTitle("step");
    TraceXsecPlots[i]->GetYaxis()->SetTitle("Parameter Variation");
    double up_fact;
    double low_fact;
    if(maxs[i]>0){
      up_fact = 1.02;
    }
    else{ up_fact = 0.98;}
    if(mins[i]>0){
      low_fact = 0.98;
    }
    else{ low_fact = 1.02;}
    TraceXsecPlots[i]->GetYaxis()->SetRangeUser(low_fact*mins[i], up_fact*maxs[i]);
    // distribution later

  }
  std::cout << "filling histograms... "<< std::endl;
  for(int i=0; i<nEntries; i++){
    for(int j=0; j<nXsec;j++){
      TraceXsecPlots[j]->SetBinContent(i,XsecValues[i][j]);
      // distribution later

    }
  }
  std::cout << "Done with histogram.Now plotting... " << std::endl;  
  // draw and save 
  TCanvas *c1 = new TCanvas("c1"," ",0,0, 800,630);
  gStyle->SetOptStat(0);
  TString canvas_name = base_name+"_DiagMCMC_xSec";
  c1->Print(canvas_name+".pdf[");
  for(int i=0; i<nXsec;i++){
   c1->cd();
   TraceXsecPlots[i]->Draw();
   c1->Print(canvas_name+".pdf");
   delete TraceXsecPlots[i];
   delete XsecValues[i];
  }

  c1->Print(canvas_name+".pdf]");

  delete[] TraceXsecPlots;
  delete[] XsecValues;
 
  std::cout << "DONE!" << std::endl;

}
