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

void plot_trace_sampleLLH(std::string file_name){

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

  std::vector<TString> SampleName_v;
  int nSamples =0;

  for(int i=0; i<nBranches; i++){
    //Get the TBranch and its name
    TBranch* br =(TBranch*)brlis->At(i);
    TString brname = br->GetName();

    //find the sample Log LH 
    if(brname.BeginsWith("LogL_sample_")){
      std::cout << "get logL from sample: " << brname <<std::endl;
      SampleName_v.push_back(brname); 
      nSamples++; 
    }
  }
  // initialize the array to hold sample LLH from each entry
  double **SampleValues;
  SampleValues = new double*[nEntries]();
  for (int i = 0; i < nEntries; i++){
    SampleValues[i] = new double[nSamples]();
    for (int j =0; j < nSamples; j++){
      SampleValues[i][j] = -999.99;
    }
  }
  // read input tree
  std::cout << "Reading input tree... " << std::endl;
  TStopwatch clock;
  clock.Start();
  // set all the branches to off except for the ones of interests
  Chain->SetBranchStatus("*",false);
  for(int i=0; i<nSamples; i++){
    Chain->SetBranchStatus(SampleName_v[i].Data(),true);
  }
  //loop over the entries
  for(int i=0; i<nEntries; i++){
    if (i%100 ==0){
       std::cout << i << "/" << nEntries << " (" << double(i)/double(nEntries)*100. << "%)" << std::endl;
    }
    // set the branch address for sample LLH
    for(int j=0; j<nSamples; j++){
      Chain->SetBranchAddress(SampleName_v[j].Data(), &SampleValues[i][j]);
    }
    // fill the arry
    Chain->GetEntry(i);
    if ( (30000 < i)&&(i < 40000)){
      std::cout << "step " << i << " sample_0 LLH: " << SampleValues[i][0] << std::endl;
    }
  }
  clock.Stop();
  std::cout << "Took " << clock.RealTime() << "s to finish for " << nEntries << " events" << std::endl;

  delete Chain;

  // make trace plots
  std::cout << "Making trace plots... " << std::endl;
  // histogram array
  TH1D **TraceSamplePlots;
  TraceSamplePlots = new TH1D*[nSamples];
  for(int i=0; i<nSamples; i++){
    TraceSamplePlots[i] = new TH1D(SampleName_v[i].Data(),SampleName_v[i].Data(),nEntries, 0,nEntries);
    TraceSamplePlots[i]->GetXaxis()->SetTitle("step");
    TraceSamplePlots[i]->GetYaxis()->SetTitle("Sample -logL");
  }
  // fill the histogram
  std::cout << "filling histograms... "<< std::endl;
  for(int i=0; i<nEntries; i++){
    for(int j=0; j<nSamples;j++){
      TraceSamplePlots[j]->SetBinContent(i,SampleValues[i][j]);
    }
  }
  std::cout << "Done with histogram.Now plotting... " << std::endl;  
  // draw and save 
  TCanvas *c1 = new TCanvas("c1"," ",0,0, 800,630);
  gStyle->SetOptStat(0);
  TString canvas_name = base_name+"_DiagMCMC_trace_LLH";
  c1->Print(canvas_name+".pdf[");
  for(int i=0; i<nSamples;i++){
   c1->cd();
   TraceSamplePlots[i]->Draw();
   c1->Print(canvas_name+".pdf");
   delete TraceSamplePlots[i];
   delete SampleValues[i];
  }

  c1->Print(canvas_name+".pdf]");

  delete[] TraceSamplePlots;
  delete[] SampleValues;
 
  std::cout << "DONE!" << std::endl;
}
