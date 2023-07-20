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

void plot_trace_corrfddet(std::string file_name){
  
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

  std::vector<TString> CorrfddetName_v;
  int nCorrfddet = 0;
  
  for(int i=0;i<nBranches;i++){
    TBranch* br=(TBranch*)brlis->At(i);
    std::string brname=br->GetName();
    
    if (brname.find("corrfddet_") == std::string::npos) continue;   
    if (brname.find("_prop")!= std::string::npos) continue;

    std::cout<< "xsec parameter: "<< brname << std::endl;
    CorrfddetName_v.push_back(TString(brname));
    nCorrfddet++;   
  }
  std::cout << "found "<< nCorrfddet << " xsec parameters. " << std::endl;
  std::cout << "initializing arrays... "<< std::endl; 
  double **CorrfddetValues;
  CorrfddetValues = new double*[nEntries]();
  for(int i =0; i<nEntries;i++){
    CorrfddetValues[i] =new double[nCorrfddet]();
    for(int j=0;j<nCorrfddet;j++){
       CorrfddetValues[i][j] = -999.99;
    }
  }
  // read input tree
  std::cout << "Reading input tree... " << std::endl;
  TStopwatch clock;
  clock.Start();
  // set all the branches to off except for the ones of interests
  Chain->SetBranchStatus("*",false);
  for(int i=0; i<nCorrfddet; i++){
    Chain->SetBranchStatus(CorrfddetName_v[i].Data(),true);
  }
  //loop over the entries
  for(int i=0; i<nEntries; i++){
    if (i%100 ==0){
       std::cout << i << "/" << nEntries << " (" << double(i)/double(nEntries)*100. << "%)" << std::endl;
    }
    // set the branch address for sample LLH
    for(int j=0; j<nCorrfddet; j++){
      Chain->SetBranchAddress(CorrfddetName_v[j].Data(), &CorrfddetValues[i][j]);
    }
    // fill the arry
    Chain->GetEntry(i);
  }
  clock.Stop();
  std::cout << "Took " << clock.RealTime() << "s to finish for " << nEntries << " events" << std::endl;
  double *maxs = new double[nCorrfddet]();
  double *mins = new double[nCorrfddet]();
  for(int i=0;i<nCorrfddet;i++){
    maxs[i]=Chain->GetMaximum(CorrfddetName_v[i].Data());
    mins[i]=Chain->GetMinimum(CorrfddetName_v[i].Data());
    if(i < 28){
      if (maxs[i] < 1.){
        maxs[i] = 1.;
      }
      if (mins[i] > -1.){
        mins[i] = -1.;
      }
    }
    else{
       if (maxs[i] < 2.){
        maxs[i] = 2.;
      }
      if (mins[i] > 0.){
        mins[i] = 0.;
      }   
    }
    std:: cout << "upper limit: " << maxs[i] << "; lower limit: " << mins[i] << std::endl;
  }
  delete Chain;
  //making plots
  TH1D **TraceCorrfddetPlots;
  TH1D **HistCorrfddetPlots;
  TraceCorrfddetPlots = new TH1D*[nCorrfddet];
  HistCorrfddetPlots =new TH1D*[nCorrfddet];
  //initialize histograms
  for(int i =0; i<nCorrfddet;i++){
    TraceCorrfddetPlots[i] = new TH1D(CorrfddetName_v[i].Data(),CorrfddetName_v[i].Data(),nEntries,0,nEntries);
    TraceCorrfddetPlots[i]->GetXaxis()->SetTitle("step");
    TraceCorrfddetPlots[i]->GetYaxis()->SetTitle("Parameter Variation");
    double range_fact = 1.04;
    TraceCorrfddetPlots[i]->GetYaxis()->SetRangeUser(range_fact*mins[i]-0.04, range_fact*maxs[i]);
    // distribution later

  }
  std::cout << "filling histograms... "<< std::endl;
  for(int i=0; i<nEntries; i++){
    for(int j=0; j<nCorrfddet;j++){
      TraceCorrfddetPlots[j]->SetBinContent(i,CorrfddetValues[i][j]);
      // distribution later

    }
  }
  std::cout << "Done with histogram.Now plotting... " << std::endl;  
  // draw and save 
  TCanvas *c1 = new TCanvas("c1"," ",0,0, 800,630);
  gStyle->SetOptStat(0);
  TString canvas_name = base_name+"_DiagMCMC_corrfddet";
  c1->Print(canvas_name+".pdf[");
  TLine* ub_1 = NULL;
  TLine* lb_1 = NULL;
  for(int i=0; i<nCorrfddet;i++){
   c1->cd();
   if ( i< 28){
     ub_1 = new TLine(TraceCorrfddetPlots[i]->GetXaxis()->GetBinLowEdge(1), 1., TraceCorrfddetPlots[i]->GetXaxis()->GetBinLowEdge(TraceCorrfddetPlots[i]->GetXaxis()->GetNbins()+1), 1.);
     lb_1 = new TLine(TraceCorrfddetPlots[i]->GetXaxis()->GetBinLowEdge(1), -1., TraceCorrfddetPlots[i]->GetXaxis()->GetBinLowEdge(TraceCorrfddetPlots[i]->GetXaxis()->GetNbins()+1), -1.);
   }
   else{
     ub_1 = new TLine(TraceCorrfddetPlots[i]->GetXaxis()->GetBinLowEdge(1), 2., TraceCorrfddetPlots[i]->GetXaxis()->GetBinLowEdge(TraceCorrfddetPlots[i]->GetXaxis()->GetNbins()+1), 2.);
     lb_1 = new TLine(TraceCorrfddetPlots[i]->GetXaxis()->GetBinLowEdge(1), 0., TraceCorrfddetPlots[i]->GetXaxis()->GetBinLowEdge(TraceCorrfddetPlots[i]->GetXaxis()->GetNbins()+1), 0.);  
   }
   ub_1->SetLineColor(kRed);
   ub_1->SetLineStyle(kDashed);
   ub_1->SetLineWidth(2);
   lb_1->SetLineColor(kRed);
   lb_1->SetLineStyle(kDashed);
   lb_1->SetLineWidth(2);

   TraceCorrfddetPlots[i]->Draw("COLZ");
   ub_1->Draw("same");
   lb_1->Draw("same");
   c1->Print(canvas_name+".pdf");
   delete TraceCorrfddetPlots[i];
   delete CorrfddetValues[i];
  }

  c1->Print(canvas_name+".pdf]");

  delete[] TraceCorrfddetPlots;
  delete[] CorrfddetValues;
 
  std::cout << "DONE!" << std::endl;

}
