#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TAxis.h>
#include <TSpline.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>

const double Pi = 3.14159265358979323846; // Pi

void MAReweight()
{
    size_t NumberOfBins = 20;
    
    int TrkID;
    float TrackTheta[5000];
    double HistogramWeight;
    
    
    TChain* Chain = new TChain("anatree");
    Chain -> Add("/home/christoph/anatrees/CCInclusiveNote/Hist_Track_pandoraNu_Vertex_pandoraNu_MA_v05_08_00_Mod.root");
    
    Chain -> SetBranchAddress("TrackCand", &TrkID);
    Chain -> SetBranchAddress("trktheta_pandoraNu", TrackTheta);
    Chain -> SetBranchAddress("eventWeight_MA", &HistogramWeight);
    
    TH1F* SelectionTheta = new TH1F("#theta-Angle","#theta-Angle",NumberOfBins,0,180);
    SelectionTheta -> SetStats(0);
    SelectionTheta -> GetXaxis() -> SetTitle("Muon #theta-Angle [#circ]");
    SelectionTheta -> GetYaxis() -> SetTitle("No. of events");
    
    for(unsigned int tree_index = 0; tree_index < Chain -> GetEntries(); tree_index++)
    {
        Chain -> GetEntry(tree_index);
        SelectionTheta -> Fill(TrackTheta[TrkID]/Pi*180, HistogramWeight);
//         std::cout << HistogramWeight << std:: endl;
    }
    
    SelectionTheta -> Draw();
}
