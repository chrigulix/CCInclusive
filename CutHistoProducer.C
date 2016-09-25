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

void CutHistoProducer()
{
//     TFile* InputFile = new TFile("/lheppc46/data/uBData/anatrees/Cut_Optimizer_prodgenie_bnb_nu_cosmic_uboone_v05_08_00_Mod.root");
    TFile* InputFile = new TFile("rootfiles/Cut_Optimizer_prodgenie_bnb_nu_cosmic_uboone_v05_08_00_Mod.root");
    
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
    
    TCanvas* C0 = new TCanvas("Flash", "Flash", 1400, 1000);
    C0->cd();
    FlashSignal->Draw();
    
    TCanvas* C1 = new TCanvas("Vertex Dist", "Vertex Dist", 1400, 1000);
    C1->cd();
    VtxDistanceSignal->Draw();
    
    TCanvas* C1a = new TCanvas("True Vertex Dist", "True Vertex Dist", 1400, 1000);
    C1a->cd();
    TrueVtxDistanceSignal->Draw();
    
    TCanvas* C2 = new TCanvas("Vertex Pos X", "Vertex Pos X", 1400, 1000);
    C2->cd();
    XVtxPosSignal->Draw();
    
    TCanvas* C3 = new TCanvas("Vertex Pos Y", "Vertex Pos Y", 1400, 1000);
    C3->cd();
    YVtxPosSignal->Draw();
    
    TCanvas* C4 = new TCanvas("Vertex Pos Z", "Vertex Pos Z", 1400, 1000);
    C4->cd();
    ZVtxPosSignal->Draw();
    
    TCanvas* C5 = new TCanvas("Flash Dist", "Flash Dist", 1400, 1000);
    C5->cd();
    FlashDistSignal->Draw();
    
    TCanvas* C6 = new TCanvas("Track Range", "Track Range", 1400, 1000);
    C6->cd();
    TrackRangeSignal->Draw();
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
    
    std::cout << "Cut maximum bin: " << SignalHist->GetBinCenter(SignalHist->GetMaximumBin()) << std::endl;
}
