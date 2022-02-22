#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <utility>
#include <vector>
#include <map>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TLine.h>
#include <TTree.h>

void ReverseScan(TH1F* RawHisto);
void Scan(TH1F* RawHisto);
void CalculateSign(TH1F* SignalHist, TH1F* BgrHist);
void MovingAverage(TH1F* SignalHist, int AverageWidth);

void DrawCutOptimisation()
{
//     std::string Folder = "/home/christoph/anatrees/CCInclusiveNote";
    std::string Folder = "./CCInclFiles";
    
    TFile* InputFile = new TFile((Folder+"/Cut_Optimizer_prodgenie_bnb_nu_cosmic_uboone_v05_08_00_Mod.root").c_str());
//     TFile* InputFile = new TFile("rootfiles/Cut_Optimizer_prodgenie_bnb_nu_cosmic_uboone_v05_08_00_Mod.root");
    
    TH1F* FlashSignal = (TH1F*)InputFile->Get("FlashSignal");
    TH1F* FlashBGR = (TH1F*)InputFile->Get("FlashBGR");
//     FlashSignal->Rebin(4);
//     FlashBGR->Rebin(4);
    ReverseScan(FlashSignal);
    ReverseScan(FlashBGR);
           
    TH1F* VtxDistanceSignal = (TH1F*)InputFile->Get("VtxDistanceSignal");
    TH1F* VtxDistanceBGR = (TH1F*)InputFile->Get("VtxDistanceBGR");
    Scan(VtxDistanceSignal);
    Scan(VtxDistanceBGR);
    
    TH1F* TrueVtxDistanceSignal = (TH1F*)InputFile->Get("TrueVtxDistanceSignal");
    TH1F* TrueVtxDistanceBGR = (TH1F*)InputFile->Get("TrueVtxDistanceBGR");
    Scan(TrueVtxDistanceSignal);
    Scan(TrueVtxDistanceBGR);
          
    TH1F* XVtxPosSignal = (TH1F*)InputFile->Get("XVtxPosSignal");
    TH1F* XVtxPosBGR = (TH1F*)InputFile->Get("XVtxPosBGR");
    ReverseScan(XVtxPosSignal);
    ReverseScan(XVtxPosBGR);    

    TH1F* YVtxPosSignal = (TH1F*)InputFile->Get("YVtxPosSignal");
    TH1F* YVtxPosBGR = (TH1F*)InputFile->Get("YVtxPosBGR");
    ReverseScan(YVtxPosSignal); 
    ReverseScan(YVtxPosBGR); 
    
    TH1F* ZVtxPosSignal = (TH1F*)InputFile->Get("ZVtxPosSignal");
    TH1F* ZVtxPosBGR = (TH1F*)InputFile->Get("ZVtxPosBGR");
    ReverseScan(ZVtxPosSignal); 
    ReverseScan(ZVtxPosBGR); 
           
    TH1F* FlashDistSignal = (TH1F*)InputFile->Get("FlashDistSignal");
    TH1F* FlashDistBGR = (TH1F*)InputFile->Get("FlashDistBGR");
    Scan(FlashDistSignal);
    Scan(FlashDistBGR);

    TH1F* TrackRangeSignal = (TH1F*)InputFile->Get("TrackRangeSignal");
    TH1F* TrackRangeBGR = (TH1F*)InputFile->Get("TrackRangeBGR");
    ReverseScan(TrackRangeSignal); 
    ReverseScan(TrackRangeBGR); 
    
    CalculateSign(FlashSignal,FlashBGR);
    CalculateSign(VtxDistanceSignal,VtxDistanceBGR);
    CalculateSign(TrueVtxDistanceSignal,TrueVtxDistanceBGR);
    CalculateSign(XVtxPosSignal,XVtxPosBGR);
    CalculateSign(YVtxPosSignal,YVtxPosBGR);
    CalculateSign(ZVtxPosSignal,ZVtxPosBGR);
    CalculateSign(FlashDistSignal,FlashDistBGR);
    CalculateSign(TrackRangeSignal,TrackRangeBGR);
    
//     MovingAverage(TrackRangeSignal,5);
    
    TCanvas* C0 = new TCanvas("Flash", "Flash", 1400, 1000);
    C0->cd();
    FlashSignal->Draw();
    
    TCanvas* C1 = new TCanvas("Vertex Dist", "Vertex Dist", 1400, 1000);
    C1->cd();
    VtxDistanceSignal->Draw();
    
    TCanvas* C2 = new TCanvas("True Vertex Dist", "True Vertex Dist", 1400, 1000);
    C2->cd();
    TrueVtxDistanceSignal->Draw();
    
    TCanvas* C3 = new TCanvas("Vertex Pos X", "Vertex Pos X", 1400, 1000);
    C3->cd();
    XVtxPosSignal->Draw();
    
    TLine* YPosLine = new TLine(20.0,0.0,20.0,84.0);
    YPosLine -> SetLineWidth(4);
    YPosLine -> SetLineColor(12);
    TLatex *YPosText = new TLatex(22.0,40.0,"Max at 20 cm");
    
    TCanvas* C4 = new TCanvas("Vertex Pos Y", "Vertex Pos Y", 1400, 1000);
    C4->cd();
    YVtxPosSignal -> SetLineColor(38);
    YVtxPosSignal -> SetLineWidth(4);
    YVtxPosSignal -> SetTitle("Fiducial Cut Optimisation");
    YVtxPosSignal -> GetXaxis() -> SetRangeUser(0,100);
    YVtxPosSignal -> GetYaxis() -> SetRangeUser(0,100);
    YVtxPosSignal -> GetXaxis() -> SetTitle("Y-Axis Fiducial Cut Value [cm]");
    YVtxPosSignal -> GetYaxis() -> SetTitle("S / #sqrt{S + B} [ ]");
    YVtxPosSignal -> GetYaxis() -> SetTitleOffset(0.95);
    YVtxPosSignal -> Draw();
    YPosLine -> Draw("SAME");
    YPosText -> Draw("SAME");
    gPad->RedrawAxis();
    C4 -> SaveAs("../images/FirstCCInclusive/CutOptFiducialY.pdf");
    
    TCanvas* C5 = new TCanvas("Vertex Pos Z", "Vertex Pos Z", 1400, 1000);
    C5->cd();
    ZVtxPosSignal->Draw();
    
    TCanvas* C6 = new TCanvas("Flash Dist", "Flash Dist", 1400, 1000);
    C6->cd();
    FlashDistSignal->Draw();
    
    TLine* LengthLine = new TLine(75.0,0.0,75.0,60.0);
    LengthLine -> SetLineWidth(4);
    LengthLine -> SetLineColor(12);
    TLatex *LengthText = new TLatex(80.0,30.0,"Max at 75 cm");
    
    TCanvas* C7 = new TCanvas("Track Range", "Track Range", 1400, 1000);
    C7->cd();
    TrackRangeSignal -> Rebin(5);
    TrackRangeSignal -> Scale(0.2);
    TrackRangeSignal -> SetLineColor(38);
    TrackRangeSignal -> SetLineWidth(4);
    TrackRangeSignal -> SetTitle("Track Range Cut Optimisation");
    TrackRangeSignal -> GetXaxis() -> SetRangeUser(0,500);
    TrackRangeSignal -> GetYaxis() -> SetRangeUser(0,70);
    TrackRangeSignal -> GetXaxis() -> SetTitle("Track Range Cut Value [cm]");
    TrackRangeSignal -> GetYaxis() -> SetTitle("S / #sqrt{S + B} [ ]");
    TrackRangeSignal -> GetYaxis() -> SetTitleOffset(0.95);
    TrackRangeSignal->Draw("HIST");
    LengthLine -> Draw("SAME");
    LengthText -> Draw("SAME");
    gPad->RedrawAxis();
    C7 -> SaveAs("../images/FirstCCInclusive/CutOptTrackRange.pdf");
}

void ReverseScan ( TH1F* RawHisto )
{
    float Sum = 0.0;
    
    for(unsigned int bin = RawHisto->GetNbinsX(); bin > 0 ; --bin)
    {
        Sum += RawHisto->GetBinContent(bin);
        
        RawHisto->SetBinContent(bin,Sum);
    }
}

void Scan ( TH1F* RawHisto )
{
    float Sum = 0.0;
    
    for(unsigned int bin = 1; bin <= RawHisto->GetNbinsX(); ++bin)
    {
        Sum += RawHisto->GetBinContent(bin);
        
        RawHisto->SetBinContent(bin,Sum);
    }
}

void CalculateSign ( TH1F* SignalHist, TH1F* BgrHist )
{
    BgrHist->Add(SignalHist);
    
    for(unsigned bin = 1; bin <= BgrHist->GetNbinsX(); ++bin)
    {
        float TempBinContent = BgrHist->GetBinContent(bin);
        BgrHist->SetBinContent(bin,std::sqrt(TempBinContent));
    }
    
    SignalHist->Divide(BgrHist);
    
    std::cout << "Cut maximum bin: " << SignalHist->GetBinCenter(SignalHist->GetMaximumBin()) << " Maximum Value: " << SignalHist->GetBinContent(SignalHist->GetMaximumBin()) << std::endl;
}

void MovingAverage(TH1F* SignalHist, int AverageWidth)
{
    unsigned int NumberOfBins = SignalHist -> GetNbinsX();
    
    std::cout << NumberOfBins << std::endl;
    
    TH1F* HistCopy = (TH1F*) SignalHist -> Clone();
    
    for(unsigned int bin_no = 1; bin_no <= NumberOfBins; bin_no++)
    {
        if((bin_no - AverageWidth) > 0 && (bin_no + AverageWidth) <= NumberOfBins)
        {
            float Filtered = 0.0;
            unsigned int BinCount = 0;
            for(int filter_bin = -AverageWidth; filter_bin <= AverageWidth; filter_bin++)
            {
                Filtered += HistCopy -> GetBinContent(filter_bin+bin_no);
                BinCount++;
            }
            Filtered /= (float)BinCount;
            
            SignalHist -> SetBinContent(bin_no, Filtered);
        }
    }
}
