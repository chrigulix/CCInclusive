#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <deque>
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

//This defines our current settings for the fiducial volume
const double FVx = 256.35;
const double FVy = 233;
const double FVz = 1036.8;
const double borderx = 10.;
const double bordery = 20.;
const double borderz = 10.;
const double Avogadro = 6.022140857e23; //mol^-1
const double ArMass = 39.948; // u
const double NoNucleons = 40;
const double Density = 1.396; // g/cm^3
const double Pi = 3.14159265358979323846; // Pi

size_t NumberOfBins = 20;


// Main Function
void DrawCosmic()
{
    float NumberOfTargets = (FVx - 2*borderx) * (FVy - 2*bordery) * (FVz - 2*borderz) * Density * Avogadro/ArMass*NoNucleons;

    std::string InputFolder = ".";
    std::string OutputFolder = ".";

    // Cosmic Selection container
    std::vector<TH1F*> CosmicTrackRange;
    std::vector<TH1F*> CosmicCosTheta;
    std::vector<TH1F*> CosmicTheta;
    std::vector<TH1F*> CosmicPhi;
    std::vector<TH1F*> CosmicMomentum;
    std::vector<TH1F*> CosmicTrackLength;
    std::vector<TH1F*> CosmicXVtxPosition;
    std::vector<TH1F*> CosmicYVtxPosition;
    std::vector<TH1F*> CosmicZVtxPosition;
    
    size_t NumberOfBins = 20;

    // Read cosmic comparison histograms
    TFile* CosmicFile = new TFile((InputFolder+"/Cosmic_Distributions_Histograms_Mod.root").c_str(),"READ");

    // cd into cosmic file
    CosmicFile -> cd();

    // Cosmic histogram labels
    std::vector<std::string> CosmicHistLabels;
    CosmicHistLabels.push_back("Data Off-Beam BNBEXT All");
    CosmicHistLabels.push_back("In Time Corsika All");
    CosmicHistLabels.push_back("Cosmic Systematics All");
    
    // Read from cosmic file
    for(auto Label : CosmicHistLabels)
    {
        CosmicTrackRange.push_back( (TH1F*) CosmicFile->Get(("Track Range "+Label).c_str()) );
        CosmicCosTheta.push_back( (TH1F*) CosmicFile->Get(("cos#theta "+Label).c_str()) );
        CosmicTheta.push_back( (TH1F*) CosmicFile->Get(("#theta-Angle "+Label).c_str()) );
        CosmicPhi.push_back( (TH1F*) CosmicFile->Get(("#phi-Angle "+Label).c_str()) );
        CosmicMomentum.push_back( (TH1F*) CosmicFile->Get(("Momentum "+Label).c_str()) );
        CosmicTrackLength.push_back( (TH1F*) CosmicFile->Get(("Track Length "+Label).c_str()) );
        CosmicXVtxPosition.push_back( (TH1F*) CosmicFile->Get(("Vertex X position "+Label).c_str()) );
        CosmicYVtxPosition.push_back( (TH1F*) CosmicFile->Get(("Vertex Y position "+Label).c_str()) );
        CosmicZVtxPosition.push_back( (TH1F*) CosmicFile->Get(("Vertex Z position "+Label).c_str()) );
    }
    
    TLegend* Legend = new TLegend(0.5,0.70,0.85,0.85);
    Legend->SetLineStyle ( 0 );
    Legend->SetLineColorAlpha ( 0,0 );
    Legend->SetFillStyle ( 0 );
    Legend->SetMargin ( 0.2 );
//     Legend->SetTextFont ( 43 );
//     Legend->SetTextSize ( 35 );
//     Legend->SetHeader("Normalised per Event","C");

//     TLegendEntry *Header = (TLegendEntry*)Legend->GetListOfPrimitives()->First();
//     Header->SetTextSize(0.05);
    
    Legend->AddEntry( CosmicTrackRange.at(0), "Data: off-Beam BNBEXT","lep" );
    Legend->AddEntry( CosmicTrackRange.at(1), "MC: inTime CORSIKA","lep" );
    

    // BEGIN DRAW COSMIC DISTRIBUTIONS -------------------------------------------------------------------------------------------------------------------
    
    TCanvas *C0 = new TCanvas("C0", "C0", 1000, 1000);
    TPad *pad1TrackRange = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1TrackRange->SetBottomMargin(0);
    pad1TrackRange->Draw();
    pad1TrackRange->cd();
    CosmicTrackRange.at(0) -> SetTitle("PandoraCosmic Track Range");
    CosmicTrackRange.at(0) -> GetXaxis() -> SetLabelOffset(999);
    CosmicTrackRange.at(0) -> GetXaxis() -> SetLabelSize(0);
    CosmicTrackRange.at(0) -> GetXaxis() -> SetRangeUser(0,450);
    CosmicTrackRange.at(0) -> GetYaxis() -> SetRangeUser(-0.1,4);
    CosmicTrackRange.at(0) -> GetYaxis() -> SetTitle("PandoraCosmic tracks per event");
    CosmicTrackRange.at(0) -> SetLineColor(9);
    CosmicTrackRange.at(0) -> Draw("SAME");
    CosmicTrackRange.at(1) -> SetLineColor(46);
    CosmicTrackRange.at(1) -> Draw("SAME");
    Legend -> Draw();
    C0->cd();          // Go back to the main canvas before defining pad2
    
    TPad *pad2TrackRange = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2TrackRange->SetTopMargin(0);
    pad2TrackRange->SetBottomMargin(0.2);
    pad2TrackRange->Draw();
    pad2TrackRange->cd(); 
    CosmicTrackRange.at(2) -> SetLineColor(13);
    CosmicTrackRange.at(2) -> SetTitle("");
    CosmicTrackRange.at(2) -> GetXaxis() -> SetLabelSize(0.09);
    CosmicTrackRange.at(2) -> GetXaxis() -> SetTitleSize(0.1);
    CosmicTrackRange.at(2) -> GetXaxis() -> SetRangeUser(0,450);
    CosmicTrackRange.at(2) -> GetYaxis() -> SetTitle("Ratio (Data / MC)");
    CosmicTrackRange.at(2) -> GetYaxis() -> SetLabelSize(0.09);
    CosmicTrackRange.at(2) -> GetYaxis() -> SetTitleSize(0.1);
    CosmicTrackRange.at(2) -> GetYaxis() -> SetTitleOffset(0.5);
    CosmicTrackRange.at(2) -> Draw("SAME");
    TLine* Line1TrackRange = new TLine(0,1,450,1);
    Line1TrackRange -> SetLineStyle(7);
    Line1TrackRange -> Draw("SAME");

    C0->SaveAs("CosmicTrackRangePandoraCosmic.pdf");
    
    Legend -> SetX1NDC(0.35);
    Legend -> SetX2NDC(0.70);
    Legend -> SetY1NDC(0.10);
    Legend -> SetY2NDC(0.25);

    TCanvas *C1 = new TCanvas("C1", "C1", 1000, 1000);
    TPad *pad1CosTheta = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1CosTheta->SetBottomMargin(0);
    pad1CosTheta->Draw();
    pad1CosTheta->cd();
    CosmicCosTheta.at(0) -> SetTitle("PandoraCosmic cos(#theta)");
    CosmicCosTheta.at(0) -> GetXaxis() -> SetLabelOffset(999);
    CosmicCosTheta.at(0) -> GetXaxis() -> SetLabelSize(0);
    CosmicCosTheta.at(0) -> GetYaxis() -> SetRangeUser(-0.1,4);
    CosmicCosTheta.at(0) -> GetYaxis() -> SetTitle("PandoraCosmic tracks per event");
    CosmicCosTheta.at(0) -> SetLineColor(9);
    CosmicCosTheta.at(0) -> Draw("SAME");
    CosmicCosTheta.at(1) -> SetLineColor(46);
    CosmicCosTheta.at(1) -> Draw("SAME");
    Legend -> Draw();
    C1->cd();          // Go back to the main canvas before defining pad2

    TPad *pad2CosTheta = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2CosTheta->SetTopMargin(0);
    pad2CosTheta->SetBottomMargin(0.2);
    pad2CosTheta->Draw();
    pad2CosTheta->cd(); 
    CosmicCosTheta.at(2) -> SetLineColor(13);
    CosmicCosTheta.at(2) -> SetTitle("");
    CosmicCosTheta.at(2) -> GetYaxis() -> SetRangeUser(0.52,1.15);
    CosmicCosTheta.at(2) -> GetXaxis() -> SetLabelSize(0.09);
    CosmicCosTheta.at(2) -> GetXaxis() -> SetTitleSize(0.1);
    CosmicCosTheta.at(2) -> GetYaxis() -> SetTitle("Ratio (Data / MC)");
    CosmicCosTheta.at(2) -> GetYaxis() -> SetLabelSize(0.09);
    CosmicCosTheta.at(2) -> GetYaxis() -> SetTitleSize(0.1);
    CosmicCosTheta.at(2) -> GetYaxis() -> SetTitleOffset(0.5);
    CosmicCosTheta.at(2) -> Draw("SAME");
    TLine* Line1CosTheta = new TLine(-1,1,1,1);
    Line1CosTheta -> SetLineStyle(7);
    Line1CosTheta -> Draw("SAME");
    
    C1->SaveAs("CosmicCosThetaPandoraCosmic.pdf");
    
    Legend -> SetX1NDC(0.50);
    Legend -> SetX2NDC(0.85);
    Legend -> SetY1NDC(0.70);
    Legend -> SetY2NDC(0.85);
    
    //     (0.5,0.70,0.85,0.85);
    
    TCanvas *C2 = new TCanvas("C2", "C2", 1000, 1000);
    TPad *pad1Theta = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1Theta->SetBottomMargin(0);
    pad1Theta->Draw();
    pad1Theta->cd();
    CosmicTheta.at(0) -> SetTitle("PandoraCosmic #theta-Angle");
    CosmicTheta.at(0) -> GetXaxis() -> SetLabelOffset(999);
    CosmicTheta.at(0) -> GetXaxis() -> SetLabelSize(0);
    CosmicTheta.at(0) -> GetYaxis() -> SetRangeUser(-0.1,5);
    CosmicTheta.at(0) -> GetYaxis() -> SetTitle("PandoraCosmic tracks per event");
    CosmicTheta.at(0) -> SetLineColor(9);
    CosmicTheta.at(0) -> Draw("SAME");
    CosmicTheta.at(1) -> SetLineColor(46);
    CosmicTheta.at(1) -> Draw("SAME");
    Legend -> Draw();
    C2->cd();          // Go back to the main canvas before defining pad2

    TPad *pad2Theta = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2Theta->SetTopMargin(0);
    pad2Theta->SetBottomMargin(0.2);
    pad2Theta->Draw();
    pad2Theta->cd(); 
    CosmicTheta.at(2) -> SetLineColor(13);
    CosmicTheta.at(2) -> SetTitle("");
    CosmicTheta.at(2) -> GetYaxis() -> SetRangeUser(0.25,1.25);
    CosmicTheta.at(2) -> GetXaxis() -> SetLabelSize(0.09);
    CosmicTheta.at(2) -> GetXaxis() -> SetTitleSize(0.1);
    CosmicTheta.at(2) -> GetYaxis() -> SetTitle("Ratio (Data / MC)");
    CosmicTheta.at(2) -> GetYaxis() -> SetLabelSize(0.09);
    CosmicTheta.at(2) -> GetYaxis() -> SetTitleSize(0.1);
    CosmicTheta.at(2) -> GetYaxis() -> SetTitleOffset(0.5);
    CosmicTheta.at(2) -> Draw("SAME");
    TLine* Line1Theta = new TLine(0,1,180,1);
    Line1Theta -> SetLineStyle(7);
    Line1Theta -> Draw("SAME");
    
    C2->SaveAs("CosmicThetaPandoraCosmic.pdf");
    
    TCanvas *C3 = new TCanvas("C3", "C3", 1000, 1000);
    TPad *pad1Phi = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1Phi->SetBottomMargin(0);
    pad1Phi->Draw();
    pad1Phi->cd();
    CosmicPhi.at(0) -> SetTitle("PandoraCosmic #phi-Angle");
    CosmicPhi.at(0) -> GetXaxis() -> SetLabelOffset(999);
    CosmicPhi.at(0) -> GetXaxis() -> SetLabelSize(0);
    CosmicPhi.at(0) -> GetYaxis() -> SetRangeUser(-0.1,10);
    CosmicPhi.at(0) -> GetYaxis() -> SetTitle("PandoraCosmic tracks per event");
    CosmicPhi.at(0) -> SetLineColor(9);
    CosmicPhi.at(0) -> Draw("SAME");
    CosmicPhi.at(1) -> SetLineColor(46);
    CosmicPhi.at(1) -> Draw("SAME");
    Legend -> Draw();
    C3->cd();          // Go back to the main canvas before defining pad2

    TPad *pad2Phi = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2Phi->SetTopMargin(0);
    pad2Phi->SetBottomMargin(0.2);
    pad2Phi->Draw();
    pad2Phi->cd(); 
    CosmicPhi.at(2) -> SetLineColor(13);
    CosmicPhi.at(2) -> SetTitle("");
//     CosmicPhi.at(2) -> GetYaxis() -> SetRangeUser(0.25,1.25);
    CosmicPhi.at(2) -> GetXaxis() -> SetLabelSize(0.09);
    CosmicPhi.at(2) -> GetXaxis() -> SetTitleSize(0.1);
    CosmicPhi.at(2) -> GetYaxis() -> SetTitle("Ratio (Data / MC)");
    CosmicPhi.at(2) -> GetYaxis() -> SetLabelSize(0.09);
    CosmicPhi.at(2) -> GetYaxis() -> SetTitleSize(0.1);
    CosmicPhi.at(2) -> GetYaxis() -> SetTitleOffset(0.5);
    CosmicPhi.at(2) -> Draw("SAME");
    TLine* Line1Phi = new TLine(-180,1,180,1);
    Line1Phi -> SetLineStyle(7);
    Line1Phi -> Draw("SAME");
    
    C3->SaveAs("CosmicPhiPandoraCosmic.pdf");
    
    TCanvas *C4 = new TCanvas("C4", "C4", 1000, 1000);
    TPad *pad1Momentum = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1Momentum->SetBottomMargin(0);
    pad1Momentum->Draw();
    pad1Momentum->cd();
    CosmicMomentum.at(0) -> SetTitle("PandoraCosmic Momentum");
    CosmicMomentum.at(0) -> GetXaxis() -> SetLabelOffset(999);
    CosmicMomentum.at(0) -> GetXaxis() -> SetLabelSize(0);
    CosmicMomentum.at(0) -> GetXaxis() -> SetRangeUser(0,1.8);
    CosmicMomentum.at(0) -> GetYaxis() -> SetRangeUser(-1,30);
    CosmicMomentum.at(0) -> GetYaxis() -> SetTitle("PandoraCosmic tracks per event");
    CosmicMomentum.at(0) -> SetLineColor(9);
    CosmicMomentum.at(0) -> Draw("SAME");
    CosmicMomentum.at(1) -> SetLineColor(46);
    CosmicMomentum.at(1) -> Draw("SAME");
    Legend -> Draw();
    C4->cd();          // Go back to the main canvas before defining pad2

    TPad *pad2Momentum = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2Momentum->SetTopMargin(0);
    pad2Momentum->SetBottomMargin(0.2);
    pad2Momentum->Draw();
    pad2Momentum->cd(); 
    CosmicMomentum.at(2) -> SetLineColor(13);
    CosmicMomentum.at(2) -> SetTitle("");
//     CosmicMomentum.at(2) -> GetYaxis() -> SetRangeUser(0.25,1.25);
    CosmicMomentum.at(2) -> GetXaxis() -> SetLabelSize(0.09);
    CosmicMomentum.at(2) -> GetXaxis() -> SetTitleSize(0.1);
    CosmicMomentum.at(2) -> GetXaxis() -> SetRangeUser(0,1.8);
    CosmicMomentum.at(2) -> GetYaxis() -> SetTitle("Ratio (Data / MC)");
    CosmicMomentum.at(2) -> GetYaxis() -> SetLabelSize(0.09);
    CosmicMomentum.at(2) -> GetYaxis() -> SetTitleSize(0.1);
    CosmicMomentum.at(2) -> GetYaxis() -> SetTitleOffset(0.5);
    CosmicMomentum.at(2) -> Draw("SAME");
    TLine* Line1Momentum = new TLine(0,1,1.95,1);
    Line1Momentum -> SetLineStyle(7);
    Line1Momentum -> Draw("SAME");
    
    C4->SaveAs("CosmicMomentumPandoraCosmic.pdf");
    
    TCanvas *C5 = new TCanvas("C5", "C5", 1000, 1000);
    TPad *pad1TrackLength = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1TrackLength->SetBottomMargin(0);
    pad1TrackLength->Draw();
    pad1TrackLength->cd();
    CosmicTrackLength.at(0) -> SetTitle("PandoraCosmic Track Length");
    CosmicTrackLength.at(0) -> GetXaxis() -> SetLabelOffset(999);
    CosmicTrackLength.at(0) -> GetXaxis() -> SetLabelSize(0);
    CosmicTrackLength.at(0) -> GetXaxis() -> SetRangeUser(0,500);
    CosmicTrackLength.at(0) -> GetYaxis() -> SetRangeUser(-1,30);
    CosmicTrackLength.at(0) -> GetYaxis() -> SetTitle("PandoraCosmic tracks per event");
    CosmicTrackLength.at(0) -> SetLineColor(9);
    CosmicTrackLength.at(0) -> Draw("SAME");
    CosmicTrackLength.at(1) -> SetLineColor(46);
    CosmicTrackLength.at(1) -> Draw("SAME");
    Legend -> Draw();
    C5->cd();          // Go back to the main canvas before defining pad2

    TPad *pad2TrackLength = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2TrackLength->SetTopMargin(0);
    pad2TrackLength->SetBottomMargin(0.2);
    pad2TrackLength->Draw();
    pad2TrackLength->cd(); 
    CosmicTrackLength.at(2) -> SetLineColor(13);
    CosmicTrackLength.at(2) -> SetTitle("");
//     CosmicTrackLength.at(2) -> GetYaxis() -> SetRangeUser(0.25,1.25);
    CosmicTrackLength.at(2) -> GetXaxis() -> SetLabelSize(0.09);
    CosmicTrackLength.at(2) -> GetXaxis() -> SetTitleSize(0.1);
    CosmicTrackLength.at(2) -> GetXaxis() -> SetRangeUser(0,500);
    CosmicTrackLength.at(2) -> GetYaxis() -> SetTitle("Ratio (Data / MC)");
    CosmicTrackLength.at(2) -> GetYaxis() -> SetLabelSize(0.09);
    CosmicTrackLength.at(2) -> GetYaxis() -> SetTitleSize(0.1);
    CosmicTrackLength.at(2) -> GetYaxis() -> SetTitleOffset(0.5);
    CosmicTrackLength.at(2) -> Draw("SAME");
    TLine* Line1TrackLength = new TLine(0,1,520,1);
    Line1TrackLength -> SetLineStyle(7);
    Line1TrackLength -> Draw("SAME");
    
    C5->SaveAs("CosmicTrackLengthPandoraCosmic.pdf");
    
    Legend -> SetX1NDC(0.35);
    Legend -> SetX2NDC(0.70);
    Legend -> SetY1NDC(0.10);
    Legend -> SetY2NDC(0.25);
    
    TCanvas *C6 = new TCanvas("C6", "C6", 1000, 1000);
    TPad *pad1XVtxPosition = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1XVtxPosition->SetBottomMargin(0);
    pad1XVtxPosition->Draw();
    pad1XVtxPosition->cd();
    CosmicXVtxPosition.at(0) -> SetTitle("PandoraCosmic Vertex Position in X");
    CosmicXVtxPosition.at(0) -> GetXaxis() -> SetLabelOffset(999);
    CosmicXVtxPosition.at(0) -> GetXaxis() -> SetLabelSize(0);
//     CosmicXVtxPosition.at(0) -> GetXaxis() -> SetRangeUser(0,500);
    CosmicXVtxPosition.at(0) -> GetYaxis() -> SetRangeUser(-0.1,7);
    CosmicXVtxPosition.at(0) -> GetYaxis() -> SetTitle("PandoraNu tracks per event");
    CosmicXVtxPosition.at(0) -> SetLineColor(9);
    CosmicXVtxPosition.at(0) -> Draw("SAME");
    CosmicXVtxPosition.at(1) -> SetLineColor(46);
    CosmicXVtxPosition.at(1) -> Draw("SAME");
    Legend -> Draw();
    C6->cd();          // Go back to the main canvas before defining pad2

    TPad *pad2XVtxPosition = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2XVtxPosition->SetTopMargin(0);
    pad2XVtxPosition->SetBottomMargin(0.2);
    pad2XVtxPosition->Draw();
    pad2XVtxPosition->cd(); 
    CosmicXVtxPosition.at(2) -> SetLineColor(13);
    CosmicXVtxPosition.at(2) -> SetTitle("");
//     CosmicXVtxPosition.at(2) -> GetYaxis() -> SetRangeUser(0.25,1.25);
    CosmicXVtxPosition.at(2) -> GetXaxis() -> SetLabelSize(0.09);
    CosmicXVtxPosition.at(2) -> GetXaxis() -> SetTitleSize(0.1);
//     CosmicXVtxPosition.at(2) -> GetXaxis() -> SetRangeUser(0,500);
    CosmicXVtxPosition.at(2) -> GetYaxis() -> SetTitle("Ratio (Data / MC)");
    CosmicXVtxPosition.at(2) -> GetYaxis() -> SetLabelSize(0.09);
    CosmicXVtxPosition.at(2) -> GetYaxis() -> SetTitleSize(0.1);
    CosmicXVtxPosition.at(2) -> GetYaxis() -> SetTitleOffset(0.5);
    CosmicXVtxPosition.at(2) -> Draw("SAME");
    TLine* Line1XVtxPosition = new TLine(0,1,256,1);
    Line1XVtxPosition -> SetLineStyle(7);
    Line1XVtxPosition -> Draw("SAME");
    
    C6->SaveAs("CosmicXVtxPositionPandoraCosmic.pdf");
    
    TCanvas *C7 = new TCanvas("C7", "C7", 1000, 1000);
    TPad *pad1YVtxPosition = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1YVtxPosition->SetBottomMargin(0);
    pad1YVtxPosition->Draw();
    pad1YVtxPosition->cd();
    CosmicYVtxPosition.at(0) -> SetTitle("PandoraCosmic Vertex Position in Y");
    CosmicYVtxPosition.at(0) -> GetXaxis() -> SetLabelOffset(999);
    CosmicYVtxPosition.at(0) -> GetXaxis() -> SetLabelSize(0);
//     CosmicYVtxPosition.at(0) -> GetXaxis() -> SetRangeUser(0,500);
    CosmicYVtxPosition.at(0) -> GetYaxis() -> SetRangeUser(-0.5,12);
    CosmicYVtxPosition.at(0) -> GetYaxis() -> SetTitle("PandoraNu tracks per event");
    CosmicYVtxPosition.at(0) -> SetLineColor(9);
    CosmicYVtxPosition.at(0) -> Draw("SAME");
    CosmicYVtxPosition.at(1) -> SetLineColor(46);
    CosmicYVtxPosition.at(1) -> Draw("SAME");
    Legend -> Draw();
    C7->cd();          // Go back to the main canvas before defining pad2

    TPad *pad2YVtxPosition = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2YVtxPosition->SetTopMargin(0);
    pad2YVtxPosition->SetBottomMargin(0.2);
    pad2YVtxPosition->Draw();
    pad2YVtxPosition->cd(); 
    CosmicYVtxPosition.at(2) -> SetLineColor(13);
    CosmicYVtxPosition.at(2) -> SetTitle("");
//     CosmicYVtxPosition.at(2) -> GetYaxis() -> SetRangeUser(0.25,1.25);
    CosmicYVtxPosition.at(2) -> GetXaxis() -> SetLabelSize(0.09);
    CosmicYVtxPosition.at(2) -> GetXaxis() -> SetTitleSize(0.1);
//     CosmicYVtxPosition.at(2) -> GetXaxis() -> SetRangeUser(0,500);
    CosmicYVtxPosition.at(2) -> GetYaxis() -> SetTitle("Ratio (Data / MC)");
    CosmicYVtxPosition.at(2) -> GetYaxis() -> SetLabelSize(0.09);
    CosmicYVtxPosition.at(2) -> GetYaxis() -> SetTitleSize(0.1);
    CosmicYVtxPosition.at(2) -> GetYaxis() -> SetTitleOffset(0.5);
    CosmicYVtxPosition.at(2) -> Draw("SAME");
    TLine* Line1YVtxPosition = new TLine(-233/2,1,233/2,1);
    Line1YVtxPosition -> SetLineStyle(7);
    Line1YVtxPosition -> Draw("SAME");
    
    C7->SaveAs("CosmicYVtxPositionPandoraCosmic.pdf");
    
    TCanvas *C8 = new TCanvas("C8", "C8", 1000, 1000);
    TPad *pad1ZVtxPosition = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1ZVtxPosition->SetBottomMargin(0);
    pad1ZVtxPosition->Draw();
    pad1ZVtxPosition->cd();
    CosmicZVtxPosition.at(0) -> SetTitle("PandoraCosmic Vertex Position in Z");
    CosmicZVtxPosition.at(0) -> GetXaxis() -> SetLabelOffset(999);
    CosmicZVtxPosition.at(0) -> GetXaxis() -> SetLabelSize(0);
//     CosmicZVtxPosition.at(0) -> GetXaxis() -> SetRangeUser(0,500);
    CosmicZVtxPosition.at(0) -> GetYaxis() -> SetRangeUser(-0.5,10);
    CosmicZVtxPosition.at(0) -> GetYaxis() -> SetTitle("PandoraNu tracks per event");
    CosmicZVtxPosition.at(0) -> SetLineColor(9);
    CosmicZVtxPosition.at(0) -> Draw("SAME");
    CosmicZVtxPosition.at(1) -> SetLineColor(46);
    CosmicZVtxPosition.at(1) -> Draw("SAME");
    Legend -> Draw();
    C8->cd();          // Go back to the main canvas before defining pad2

    TPad *pad2ZVtxPosition = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2ZVtxPosition->SetTopMargin(0);
    pad2ZVtxPosition->SetBottomMargin(0.2);
    pad2ZVtxPosition->Draw();
    pad2ZVtxPosition->cd(); 
    CosmicZVtxPosition.at(2) -> SetLineColor(13);
    CosmicZVtxPosition.at(2) -> SetTitle("");
    CosmicZVtxPosition.at(2) -> GetYaxis() -> SetRangeUser(0.75,1.45);
    CosmicZVtxPosition.at(2) -> GetXaxis() -> SetLabelSize(0.09);
    CosmicZVtxPosition.at(2) -> GetXaxis() -> SetTitleSize(0.1);
//     CosmicZVtxPosition.at(2) -> GetXaxis() -> SetRangeUser(0,500);
    CosmicZVtxPosition.at(2) -> GetYaxis() -> SetTitle("Ratio (Data / MC)");
    CosmicZVtxPosition.at(2) -> GetYaxis() -> SetLabelSize(0.09);
    CosmicZVtxPosition.at(2) -> GetYaxis() -> SetTitleSize(0.1);
    CosmicZVtxPosition.at(2) -> GetYaxis() -> SetTitleOffset(0.5);
    CosmicZVtxPosition.at(2) -> Draw("SAME");
    TLine* Line1ZVtxPosition = new TLine(0,1,1036,1);
    Line1ZVtxPosition -> SetLineStyle(7);
    Line1ZVtxPosition -> Draw("SAME");
    
    C8->SaveAs("CosmicZVtxPositionPandoraCosmic.pdf");

    // END DRAW COSMIC DISTRIBUTIONS ---------------------------------------------------------------------------------------------------------------------
}
