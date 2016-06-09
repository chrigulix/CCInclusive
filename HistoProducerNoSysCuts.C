#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <iostream>

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

//This defines our current settings for the fiducial volume
double FVx = 256.35;
double FVy = 233;
double FVz = 1036.8;
double borderx = 10.;
double bordery = 20.;
double borderz = 10.;

float GetMaximum(const std::vector<TH1F*>& HistVector);
void AddFirstTwoHistograms(std::vector<TH1F*>& HistVector, float Weight);
void AddFirstTwoHistograms2D(std::vector<TH2F*>& HistVector, float Weight);
float CalcLength(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2);
double FlashTrackDist(double flash, double start, double end);
bool inDeadRegion(double y, double z);
std::vector<TSpline5> Systematics();
void AdjustSysError(std::vector<TH1F*>& HistVector);
bool inFV(double x, double y, double z);

void HistoProducerNoSysCuts()
{
    TGaxis::SetMaxDigits(4);
    
    //     std::string TrackProdName="pandoraNuKHit";
//     std::string TrackProdName = "pandoraCosmic";
    std::string TrackProdName="pandoraNu";
//     std::string TrackProdName="pmtrack";
//     std::string TrackProdName="pandoraNuPMA";
//     std::string TrackProdName="trackkalmanhit";

//     std::string  VertexProdName="nuvtx";
//     std::string VertexProdName="pandoraCosmic";
    std::string VertexProdName = "pandoraNu";
//     std::string VertexProdName = "pmtrack";

//     std::string SelectionLabel = "_Old";
    std::string SelectionLabel = "_Mod";
//     std::string SelectionLabel = "_New";
    
    std::string FileType = "png";
//     std::string FileType = "pdf";
    
    std::vector<TChain*> ChainVec;

    std::vector<std::string> DataLabel;
    std::vector<std::string> MCLabel;
    std::vector<std::string> GenLabel;
    std::vector<std::string> BgrLabel;

    std::vector<float> ScalingFactors;
    ScalingFactors.push_back(1);
    ScalingFactors.push_back(1.28711);
    ScalingFactors.push_back(1);

    // Binning
//     unsigned int NumberOfBins = 2000000;
    unsigned int NumberOfBins = 20;
    unsigned int NumberOf2DBins = 10;

    std::vector<TH1F*> SelectionTrackRange;
    std::vector<TH1F*> SelectionEnergy;
    std::vector<TH1F*> SelectionMomentum;
    std::vector<TH1F*> SelectionTheta;
    std::vector<TH1F*> SelectionCosTheta;
    std::vector<TH1F*> SelectionPhi;

    std::vector<TH1F*> SelXTrackStartEnd;
    std::vector<TH1F*> SelYTrackStartEnd;
    std::vector<TH1F*> SelZTrackStartEnd;

    std::vector<TH1F*> SelXVtxPosition;
    std::vector<TH1F*> SelYVtxPosition;
    std::vector<TH1F*> SelZVtxPosition;

    std::vector<TH1F*> BgrTrackRange;
    std::vector<TH1F*> BgrEnergy;
    std::vector<TH1F*> BgrMomentum;
    std::vector<TH1F*> BgrTheta;
    std::vector<TH1F*> BgrCosTheta;
    std::vector<TH1F*> BgrPhi;

    std::vector<TH1F*> BgrXTrackStartEnd;
    std::vector<TH1F*> BgrYTrackStartEnd;
    std::vector<TH1F*> BgrZTrackStartEnd;

    std::vector<TH1F*> BgrXVtxPosition;
    std::vector<TH1F*> BgrYVtxPosition;
    std::vector<TH1F*> BgrZVtxPosition;

    std::vector<TH2F*> PhiVsTheta;
    std::vector<TH2F*> PhiVsXPos;
    std::vector<TH2F*> PhiVsYPos;
    std::vector<TH2F*> PhiVsZPos;
    std::vector<TH2F*> XPosVsYPos;
    std::vector<TH2F*> ZPosVsYPos;
    std::vector<TH2F*> RangeVsPE;
    std::vector<TH2F*> RangeVsYPos;
    std::vector<TH2F*> PhiVsFlashTrackDist;

    TF1* SinTheta = new TF1("const","sin(x)",0,3.142);

    THStack* StackBgrTrackRange = new THStack("Bgr Track Range","Bgr Track Range");
    THStack* StackBgrEnergy = new THStack("Bgr Energy","Bgr Energy");
    THStack* StackBgrMomentum = new THStack("Bgr Momentum","Bgr Momentum");
    THStack* StackBgrTheta = new THStack("Bgr Theta","Bgr Theta");
    THStack* StackBgrCosTheta = new THStack("Bgr Theta","Bgr Theta");
    THStack* StackBgrPhi = new THStack("Bgr Phi","Bgr Phi");
    THStack* StackBgrXTrackStartEnd = new THStack("Bgr X Start","Bgr X Start");
    THStack* StackBgrYTrackStartEnd = new THStack("Bgr Y Start","Bgr Y Start");
    THStack* StackBgrZTrackStartEnd = new THStack("Bgr Z Start","Bgr Z Start");
    THStack* StackBgrXVtxPosition = new THStack("Bgr X Vertex","Bgr X Vertex");
    THStack* StackBgrYVtxPosition = new THStack("Bgr X Vertex","Bgr X Vertex");
    THStack* StackBgrZVtxPosition = new THStack("Bgr X Vertex","Bgr X Vertex");

    TLegend* LegendData = new TLegend(0.7,0.8,0.9,0.9);
    LegendData->SetHeader("Data Sample");

    DataLabel.push_back("On-Beam BNB");
    DataLabel.push_back("Off-Beam BNB ext");

    TLegend* LegendMC = new TLegend(0.6,0.72,0.9,0.9);
//     LegendMC->SetHeader("Data Type");

    MCLabel.push_back("On-Beam Minus Off-Beam Sample");
    MCLabel.push_back("Selection on MC BNB+Cosmic with Stat. Error");
//     MCLabel.push_back("MC BNB+Cosmic Sys. Error");

    TLegend* FlashLabel = new TLegend(0.7,0.7,0.9,0.9);
//     FlashLabel->SetHeader("Generator Type");

    GenLabel.push_back("Data On-Beam BNB");
    GenLabel.push_back("Data Off-Beam BNBEXT");
    GenLabel.push_back("MC Prodgenie BNB Nu Cosmic");
//     GenLabel.push_back("MC Systematic Errors");
//     GenLabel.push_back("MC Prodgenie BNB Nu");
//     GenLabel.push_back("MC Prodcosmic Corsika in-Time");

    BgrLabel.push_back("Bgr #bar{#nu}_{#mu} Events MC BNB+Cosmic");
    BgrLabel.push_back("Bgr #nu_{e} Events MC BNB+Cosmic");
    BgrLabel.push_back("Bgr #nu_{#mu} CC Out of FV MC BNB+Cosmic");
    BgrLabel.push_back("Bgr NC Events MC BNB+Cosmic");
    BgrLabel.push_back("Bgr Cosmic Events MC BNB+Cosmic");

    std::vector<TSpline5> SystematicErrors = Systematics();

    std::vector<unsigned int> ColorMap = {13,28,42,30,38};

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add(("/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_onbeam_bnb_v05_08_00_1"+ SelectionLabel +".root").c_str());
    ChainVec.back() -> Add(("/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_onbeam_bnb_v05_08_00_2"+ SelectionLabel +".root").c_str());
    ChainVec.back() -> Add(("/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_onbeam_bnb_v05_08_00_3"+ SelectionLabel +".root").c_str());

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add(("/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_offbeam_bnbext_v05_08_00_1"+ SelectionLabel +".root").c_str());
    ChainVec.back() -> Add(("/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_offbeam_bnbext_v05_08_00_2"+ SelectionLabel +".root").c_str());

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add(("/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_prodgenie_bnb_nu_cosmic_uboone_v05_08_00"+ SelectionLabel +".root").c_str());

    for(const auto& Label : GenLabel)
    {
        SelectionTrackRange.push_back(new TH1F(("Track Range"+Label).c_str(),"Track Range of Selected Track",NumberOfBins,0,1036.8));
        SelectionTrackRange.back()->SetStats(0);
        SelectionTrackRange.back()->GetXaxis()->SetTitle("Track Range [cm]");
        SelectionTrackRange.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");

        SelectionTheta.push_back(new TH1F(("#theta-Angle"+Label).c_str(),"#theta-Angle of Selected Track",NumberOfBins,0,3.142));
        SelectionTheta.back()->SetStats(0);
        SelectionTheta.back()->GetXaxis()->SetTitle("#theta [rad]");
        SelectionTheta.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{d#theta}");
        SelectionTheta.back()->GetYaxis()->SetTitleOffset(1.3);

        SelectionCosTheta.push_back(new TH1F(("cos#theta-Angle"+Label).c_str(),"cos#theta of Selected Track",NumberOfBins,-1,1));
        SelectionCosTheta.back()->SetStats(0);
        SelectionCosTheta.back()->GetXaxis()->SetTitle("cos#theta [ ]");
        SelectionCosTheta.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{d(cos#theta)}");
        SelectionCosTheta.back()->GetYaxis()->SetTitleOffset(1.3);

        SelectionPhi.push_back(new TH1F(("#phi-Angle"+Label).c_str(),"#phi-Angle of Selected Track",NumberOfBins,-3.142,3.142));
        SelectionPhi.back()->SetStats(0);
        SelectionPhi.back()->GetXaxis()->SetTitle("#phi angle [rad]");
        SelectionPhi.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{d#phi}");
        SelectionPhi.back()->GetYaxis()->SetTitleOffset(1.3);

        SelectionEnergy.push_back(new TH1F(("Energy"+Label).c_str(),"Energy of Selected Track",NumberOfBins,0,3));
        SelectionEnergy.back()->SetStats(0);
        SelectionEnergy.back()->GetXaxis()->SetTitle("Muon Kinetic Energy [MeV]");
        SelectionEnergy.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dE}");
        SelectionEnergy.back()->GetYaxis()->SetTitleOffset(1.3);
        
        SelectionMomentum.push_back(new TH1F(("Momentum"+Label).c_str(),"Momentum of Selected Track",NumberOfBins,0,3));
        SelectionMomentum.back()->SetStats(0);
        SelectionMomentum.back()->GetXaxis()->SetTitle("Muon Momentum [GeV/c]");
        SelectionMomentum.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dp}");
        SelectionMomentum.back()->GetYaxis()->SetTitleOffset(1.3);

        SelXTrackStartEnd.push_back(new TH1F(("XTrack"+Label).c_str(),"X Track Start & End Positions",NumberOfBins,0,256));
        SelXTrackStartEnd.back()->SetStats(0);
        SelXTrackStartEnd.back()->GetXaxis()->SetTitle("x [cm]");
        SelXTrackStartEnd.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
        SelXTrackStartEnd.back()->GetYaxis()->SetTitleOffset(1.3);

        SelYTrackStartEnd.push_back(new TH1F(("YTrack"+Label).c_str(),"Y Track Start & End Positions",NumberOfBins,-233/2,233/2));
        SelYTrackStartEnd.back()->SetStats(0);
        SelYTrackStartEnd.back()->GetXaxis()->SetTitle("y [cm]");
        SelYTrackStartEnd.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dy}");
        SelYTrackStartEnd.back()->GetYaxis()->SetTitleOffset(1.3);

        SelZTrackStartEnd.push_back(new TH1F(("ZTrack"+Label).c_str(),"Z Track Start & End Positions",NumberOfBins,0,1036.8));
        SelZTrackStartEnd.back()->SetStats(0);
        SelZTrackStartEnd.back()->GetXaxis()->SetTitle("z [cm]");
        SelZTrackStartEnd.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dz}");
        SelZTrackStartEnd.back()->GetYaxis()->SetTitleOffset(1.3);

        SelXVtxPosition.push_back(new TH1F(("XVertex"+Label).c_str(),"X Vertex Position",NumberOfBins,0,256));
        SelXVtxPosition.back()->SetStats(0);
        SelXVtxPosition.back()->GetXaxis()->SetTitle("x [cm]");
        SelXVtxPosition.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
        SelXVtxPosition.back()->GetYaxis()->SetTitleOffset(1.3);

        SelYVtxPosition.push_back(new TH1F(("YVertex"+Label).c_str(),"Y Vertex Position",NumberOfBins,-233/2,233/2));
        SelYVtxPosition.back()->SetStats(0);
        SelYVtxPosition.back()->GetXaxis()->SetTitle("y [cm]");
        SelYVtxPosition.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dy}");
        SelYVtxPosition.back()->GetYaxis()->SetTitleOffset(1.3);

        SelZVtxPosition.push_back(new TH1F(("ZVertex"+Label).c_str(),"Z Vertex Position",NumberOfBins,0,1036.8));
        SelZVtxPosition.back()->SetStats(0);
        SelZVtxPosition.back()->GetXaxis()->SetTitle("z [cm]");
        SelZVtxPosition.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dz}");
        SelZVtxPosition.back()->GetYaxis()->SetTitleOffset(1.3);

        PhiVsTheta.push_back(new TH2F(("PhiVsTheta"+Label).c_str(),"Phi Vs. Theta",NumberOf2DBins,-3.142,3.142,NumberOf2DBins,0,3.142));
        PhiVsTheta.back()->SetStats(0);
        PhiVsTheta.back()->GetXaxis()->SetTitle("#phi [rad]");
        PhiVsTheta.back()->GetYaxis()->SetTitle("#theta [rad]");

        PhiVsXPos.push_back(new TH2F(("PhiVsXPos"+Label).c_str(),"Phi Vs. X-Position",NumberOf2DBins,-3.142,3.142,NumberOf2DBins,0,256));
        PhiVsXPos.back()->SetStats(0);
        PhiVsXPos.back()->GetXaxis()->SetTitle("#phi [rad]");
        PhiVsXPos.back()->GetYaxis()->SetTitle("x [cm]");

        PhiVsYPos.push_back(new TH2F(("PhiVsYPos"+Label).c_str(),"Phi Vs. Y-Position",NumberOf2DBins,-3.142,3.142,NumberOf2DBins,-233/2,233/2));
        PhiVsYPos.back()->SetStats(0);
        PhiVsYPos.back()->GetXaxis()->SetTitle("#phi [rad]");
        PhiVsYPos.back()->GetYaxis()->SetTitle("y [cm]");

        PhiVsZPos.push_back(new TH2F(("PhiVsZPos"+Label).c_str(),"Phi Vs. Z-Position",NumberOf2DBins,-3.142,3.142,NumberOf2DBins,0,1036.8));
        PhiVsZPos.back()->SetStats(0);
        PhiVsZPos.back()->GetXaxis()->SetTitle("#phi [rad]");
        PhiVsZPos.back()->GetYaxis()->SetTitle("z [cm]");

        XPosVsYPos.push_back(new TH2F(("XPosVsYPos"+Label).c_str(),"X-Position Vs. Y-Position",NumberOf2DBins,0,256,NumberOf2DBins,-233/2,233/2));
        XPosVsYPos.back()->SetStats(0);
        XPosVsYPos.back()->GetXaxis()->SetTitle("x [cm]");
        XPosVsYPos.back()->GetYaxis()->SetTitle("y [cm]");

        ZPosVsYPos.push_back(new TH2F(("ZPosVsYPos"+Label).c_str(),"Z-Position Vs. Y-Position",NumberOf2DBins,0,1036.8,NumberOf2DBins,-233/2,233/2));
        ZPosVsYPos.back()->SetStats(0);
        ZPosVsYPos.back()->GetXaxis()->SetTitle("z [cm]");
        ZPosVsYPos.back()->GetYaxis()->SetTitle("y [cm]");

        RangeVsPE.push_back(new TH2F(("RangeVsPE"+Label).c_str(),"Track Range Vs. PE",NumberOf2DBins,0,700,NumberOf2DBins,50,500));
        RangeVsPE.back()->SetStats(0);
        RangeVsPE.back()->GetXaxis()->SetTitle("Track Range [cm]");
        RangeVsPE.back()->GetYaxis()->SetTitle("# PE [ ]");

        RangeVsYPos.push_back(new TH2F(("RangeVsYPos"+Label).c_str(),"Y-Position Vs. Track Range",NumberOf2DBins,0,700,NumberOf2DBins,-233/2,233/2));
        RangeVsYPos.back()->SetStats(0);
        RangeVsYPos.back()->GetXaxis()->SetTitle("Track Range [cm]");
        RangeVsYPos.back()->GetYaxis()->SetTitle("y [cm]");

        PhiVsFlashTrackDist.push_back(new TH2F(("PhiVsFlashTrackDist"+Label).c_str(),"Phi angle Vs. Track to Flash Distance",NumberOf2DBins,-3.142,3.142,NumberOf2DBins,0,14));
        PhiVsFlashTrackDist.back()->SetStats(0);
        PhiVsFlashTrackDist.back()->GetXaxis()->SetTitle("#phi angle [rad]");
        PhiVsFlashTrackDist.back()->GetYaxis()->SetTitle("Track to Flash Distance [cm]");
    }

    unsigned int BgrCount = 0;
    for(const auto& Label : BgrLabel)
    {
        BgrTrackRange.push_back(new TH1F(("Track Range"+Label).c_str(),"Track Range of Selected Track",NumberOfBins,0,1000));
        BgrTrackRange.back()->SetStats(0);
        BgrTrackRange.back()->SetFillColor(ColorMap.at(BgrCount));
        BgrTrackRange.back()->GetXaxis()->SetTitle("Track Range [cm]");
        BgrTrackRange.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");

        BgrTheta.push_back(new TH1F(("#theta-Angle"+Label).c_str(),"#theta of Selected Track",NumberOfBins,0,3.142));
        BgrTheta.back()->SetStats(0);
        BgrTheta.back()->SetFillColor(ColorMap.at(BgrCount));
        BgrTheta.back()->GetXaxis()->SetTitle("#theta [rad]");
        BgrTheta.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{d#theta}");
        BgrTheta.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrCosTheta.push_back(new TH1F(("cos#theta-Angle"+Label).c_str(),"cos#theta of Selected Track",NumberOfBins,-1,1));
        BgrCosTheta.back()->SetStats(0);
        BgrCosTheta.back()->SetFillColor(ColorMap.at(BgrCount));
        BgrCosTheta.back()->GetXaxis()->SetTitle("cos#theta [ ]");
        BgrCosTheta.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{d(cos#theta)}");
        BgrCosTheta.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrPhi.push_back(new TH1F(("#phi-Angle"+Label).c_str(),"#phi-Angle of Selected Track",NumberOfBins,-3.142,3.142));
        BgrPhi.back()->SetStats(0);
        BgrPhi.back()->SetFillColor(ColorMap.at(BgrCount));
        BgrPhi.back()->GetXaxis()->SetTitle("#phi angle [rad]");
        BgrPhi.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{d#phi}");
        BgrPhi.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrEnergy.push_back(new TH1F(("Energy"+Label).c_str(),"Energy of Selected Track",NumberOfBins,0,3));
        BgrEnergy.back()->SetStats(0);
        BgrEnergy.back()->SetFillColor(ColorMap.at(BgrCount));
        BgrEnergy.back()->GetXaxis()->SetTitle("Muon Kinetic Energy [MeV]");
        BgrEnergy.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dE}");
        BgrEnergy.back()->GetYaxis()->SetTitleOffset(1.3);
        
        BgrMomentum.push_back(new TH1F(("Momentum"+Label).c_str(),"Momentum of Selected Track",NumberOfBins,0,3));
        BgrMomentum.back()->SetStats(0);
        BgrMomentum.back()->SetFillColor(ColorMap.at(BgrCount));
        BgrMomentum.back()->GetXaxis()->SetTitle("Muon Momentum [GeV/c]");
        BgrMomentum.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dp}");
        BgrMomentum.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrXTrackStartEnd.push_back(new TH1F(("XTrack"+Label).c_str(),"X Track Start & End Positions",NumberOfBins,0,256));
        BgrXTrackStartEnd.back()->SetStats(0);
        BgrXTrackStartEnd.back()->SetFillColor(ColorMap.at(BgrCount));
        BgrXTrackStartEnd.back()->GetXaxis()->SetTitle("x [cm]");
        BgrXTrackStartEnd.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
        BgrXTrackStartEnd.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrYTrackStartEnd.push_back(new TH1F(("YTrack"+Label).c_str(),"Y Track Start & End Positions",NumberOfBins,-233/2,233/2));
        BgrYTrackStartEnd.back()->SetStats(0);
        BgrYTrackStartEnd.back()->SetFillColor(ColorMap.at(BgrCount));
        BgrYTrackStartEnd.back()->GetXaxis()->SetTitle("y [cm]");
        BgrYTrackStartEnd.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dy}");
        BgrYTrackStartEnd.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrZTrackStartEnd.push_back(new TH1F(("ZTrack"+Label).c_str(),"Z Track Start & End Positions",NumberOfBins,0,1036.8));
        BgrZTrackStartEnd.back()->SetStats(0);
        BgrZTrackStartEnd.back()->SetFillColor(ColorMap.at(BgrCount));
        BgrZTrackStartEnd.back()->GetXaxis()->SetTitle("z [cm]");
        BgrZTrackStartEnd.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dz}");
        BgrZTrackStartEnd.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrXVtxPosition.push_back(new TH1F(("XVertex"+Label).c_str(),"X Vertex Position",NumberOfBins,0,256));
        BgrXVtxPosition.back()->SetStats(0);
        BgrXVtxPosition.back()->SetFillColor(ColorMap.at(BgrCount));
        BgrXVtxPosition.back()->GetXaxis()->SetTitle("x [cm]");
        BgrXVtxPosition.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
        BgrXVtxPosition.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrYVtxPosition.push_back(new TH1F(("YVertex"+Label).c_str(),"Y Vertex Position",NumberOfBins,-233/2,233/2));
        BgrYVtxPosition.back()->SetStats(0);
        BgrYVtxPosition.back()->SetFillColor(ColorMap.at(BgrCount));
        BgrYVtxPosition.back()->GetXaxis()->SetTitle("y [cm]");
        BgrYVtxPosition.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dy}");
        BgrYVtxPosition.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrZVtxPosition.push_back(new TH1F(("ZVertex"+Label).c_str(),"Z Vertex Position",NumberOfBins,0,1036.8));
        BgrZVtxPosition.back()->SetStats(0);
        BgrZVtxPosition.back()->SetFillColor(ColorMap.at(BgrCount));
        BgrZVtxPosition.back()->GetXaxis()->SetTitle("z [cm]");
        BgrZVtxPosition.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dz}");
        BgrZVtxPosition.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrCount++;
    }

    int Run;
    int Subrun;
    int Event;

    int TrkID;
    int VtxID;

    float FlashPE[5000];
    int NumberOfFlashes;
    float FlashTime[5000];
    float ZFlashCenter[5000];

    int MCTrkID;
    int CCNCFlag[10];
    int TruthMode[10];
    int PDGTruth[5000];
    float NuEnergyTruth[10];
    float nuvtxx_truth[10]; //true vertex x (in cm)
    float nuvtxy_truth[10];
    float nuvtxz_truth[10];

    short TrkBestPlane[5000];
    short TrkOrigin[5000][3];

    float TrackTheta[5000];
    float TrackPhi[5000];
    float TrackMomentum[5000];

    float XTrackStart[5000];
    float YTrackStart[5000];
    float ZTrackStart[5000];

    float XTrackEnd[5000];
    float YTrackEnd[5000];
    float ZTrackEnd[5000];

    float XVertexPosition[500];
    float YVertexPosition[500];
    float ZVertexPosition[500];

    float KineticEnergy[5000][3];

    double beammin;
    double beammax;

    std::ofstream DataToLookAt("SelectedData.txt",std::ios::trunc);
    
    unsigned int OnBeamRemain = 0;
    unsigned int OffBeamRemain = 0;
    
    float OnBeamFactor;
    float OffBeamFactor;

    for(unsigned int file_no = 0; file_no < ChainVec.size(); file_no++)
    {
        if(!file_no)
        {
            beammin = 3.55 - 0.36;
            beammax = 5.15 - 0.36;
        }
        else
        {
            beammin = 3.55;
            beammax = 5.15;
        }

        ChainVec.at(file_no) -> SetBranchAddress("run", &Run);
        ChainVec.at(file_no) -> SetBranchAddress("subrun", &Subrun);
        ChainVec.at(file_no) -> SetBranchAddress("event", &Event);

        ChainVec.at(file_no) -> SetBranchAddress("TrackCand", &TrkID);
        ChainVec.at(file_no) -> SetBranchAddress("VertexCand", &VtxID);

        ChainVec.at(file_no) -> SetBranchAddress("flash_pe", FlashPE);
        ChainVec.at(file_no) -> SetBranchAddress("no_flashes", &NumberOfFlashes);
        ChainVec.at(file_no) -> SetBranchAddress("flash_time", FlashTime);
        ChainVec.at(file_no) -> SetBranchAddress("flash_zcenter", ZFlashCenter);

        ChainVec.at(file_no) -> SetBranchAddress("MCTrackCand", &MCTrkID);
        ChainVec.at(file_no) -> SetBranchAddress("ccnc_truth", CCNCFlag);
        ChainVec.at(file_no) -> SetBranchAddress("mode_truth", TruthMode);
        ChainVec.at(file_no) -> SetBranchAddress("pdg", PDGTruth);
        ChainVec.at(file_no) -> SetBranchAddress("enu_truth", NuEnergyTruth);
        ChainVec.at(file_no) -> SetBranchAddress("nuvtxx_truth", nuvtxx_truth);
        ChainVec.at(file_no) -> SetBranchAddress("nuvtxy_truth", nuvtxy_truth);
        ChainVec.at(file_no) -> SetBranchAddress("nuvtxz_truth", nuvtxz_truth);
        ChainVec.at(file_no) -> SetBranchAddress(("trkorigin_"+TrackProdName).c_str(), TrkOrigin);
        ChainVec.at(file_no) -> SetBranchAddress(("trkpidbestplane_"+TrackProdName).c_str(), TrkBestPlane);

        ChainVec.at(file_no) -> SetBranchAddress(("trkke_"+TrackProdName).c_str(), KineticEnergy);
        ChainVec.at(file_no) -> SetBranchAddress(("trkmomrange_"+TrackProdName).c_str(), TrackMomentum);
        ChainVec.at(file_no) -> SetBranchAddress(("trktheta_"+TrackProdName).c_str(), TrackTheta);
        ChainVec.at(file_no) -> SetBranchAddress(("trkphi_"+TrackProdName).c_str(),TrackPhi);

        ChainVec.at(file_no) -> SetBranchAddress(("trkstartx_"+TrackProdName).c_str(),XTrackStart);
        ChainVec.at(file_no) -> SetBranchAddress(("trkstarty_"+TrackProdName).c_str(),YTrackStart);
        ChainVec.at(file_no) -> SetBranchAddress(("trkstartz_"+TrackProdName).c_str(),ZTrackStart);

        ChainVec.at(file_no) -> SetBranchAddress(("trkendx_"+TrackProdName).c_str(),XTrackEnd);
        ChainVec.at(file_no) -> SetBranchAddress(("trkendy_"+TrackProdName).c_str(),YTrackEnd);
        ChainVec.at(file_no) -> SetBranchAddress(("trkendz_"+TrackProdName).c_str(),ZTrackEnd);

        if(VertexProdName != "nuvtx")
        {
            ChainVec.at(file_no) -> SetBranchAddress(("vtxx_"+VertexProdName).c_str(), XVertexPosition);
            ChainVec.at(file_no) -> SetBranchAddress(("vtxy_"+VertexProdName).c_str(), YVertexPosition);
            ChainVec.at(file_no) -> SetBranchAddress(("vtxz_"+VertexProdName).c_str(), ZVertexPosition);
        }
        else
        {
            ChainVec.at(file_no) -> SetBranchAddress("nuvtxx", XVertexPosition);
            ChainVec.at(file_no) -> SetBranchAddress("nuvtxy", YVertexPosition);
            ChainVec.at(file_no) -> SetBranchAddress("nuvtxz", ZVertexPosition);
        }

        unsigned int nubar = 0;
        unsigned int nue = 0;
        unsigned int outFV = 0;
        unsigned int NCnu = 0;
        unsigned int Cosmic = 0;
        unsigned int UnknownOrigin = 0;
        unsigned int Signal = 0;

        unsigned int nuQE = 0;
        unsigned int nuRES = 0;
        unsigned int nuDIS = 0;
        unsigned int nuCOH = 0;

        unsigned int negPhi = 0;
        unsigned int posPhi = 0;

        float XFVCutValue = 10; //10
        float YFVCutValue = 20; //20
        float ZFVCutValue = 10; //10
        float FlashTrackCut = 80; //80

        for(unsigned int tree_index = 0; tree_index < ChainVec.at(file_no)->GetEntries(); tree_index++)
        {
            if(!(tree_index % 100)) std::cout << "Event\t" << tree_index << "\t of \t" << ChainVec.at(file_no)->GetEntries() << std::endl;

            ChainVec.at(file_no)->GetEntry(tree_index);

            float FlashMax = 0.0;
            float ZFlashCenterMax = 0.0;

            for(int flash_no = 0; flash_no < NumberOfFlashes; flash_no++)
            {
                if( FlashTime[flash_no] > beammin && FlashTime[flash_no] < beammax && FlashPE[flash_no] > FlashMax)
                {
                    FlashMax = FlashPE[flash_no];
                    ZFlashCenterMax = ZFlashCenter[flash_no];
                }
            }
            
            Signal++;

//             if(FlashTrackDist(ZFlashCenterMax,ZTrackStart[TrkID],ZTrackEnd[TrkID]) < FlashTrackCut  && YTrackStart[TrkID] < (233./2.-30) && YTrackStart[TrkID] > (-233./2.+30) && YTrackEnd[TrkID] < (233./2.-30) && YTrackEnd[TrkID] > (-233./2.+30))
//             if( YTrackStart[TrkID] < (233./2.-YFVCutValue) && YTrackStart[TrkID] > (-233./2.+YFVCutValue) && YTrackEnd[TrkID] < (233./2.-YFVCutValue) && YTrackEnd[TrkID] > (-233./2.+YFVCutValue) &&
//                     ZTrackStart[TrkID] < (1036.8-ZFVCutValue) && ZTrackStart[TrkID] > ZFVCutValue && ZTrackEnd[TrkID] < (1036.8-ZFVCutValue) && ZTrackEnd[TrkID] > ZFVCutValue &&
//                     FlashTrackDist(ZFlashCenterMax,ZTrackStart[TrkID],ZTrackEnd[TrkID]) < FlashTrackCut )
//             if(!inDeadRegion(YTrackStart[TrkID],ZTrackStart[TrkID]) && !inDeadRegion(YTrackEnd[TrkID],ZTrackEnd[TrkID]))
            if(!(Run < 5750 && Run > 5650) || file_no > 1)
//             if(true)
            {
                if(file_no == 0)
                {
                    DataToLookAt << Run << " " << Subrun << " " << Event << "\n";
                    OnBeamRemain++;
                }
                else if(file_no == 1)
                {
                    OffBeamRemain++;
                }

                RangeVsPE.at(file_no)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),FlashMax);

                if(FlashTrackDist(ZFlashCenterMax,ZTrackStart[TrkID],ZTrackEnd[TrkID]))
                {
                    PhiVsFlashTrackDist.at(file_no)->Fill(TrackPhi[TrkID],FlashTrackDist(ZFlashCenterMax,ZTrackStart[TrkID],ZTrackEnd[TrkID]));
                }

                SelectionTrackRange.at(file_no)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                SelectionTheta.at(file_no)->Fill(TrackTheta[TrkID]);
                SelectionCosTheta.at(file_no)->Fill(cos(TrackTheta[TrkID]));
                SelectionPhi.at(file_no)->Fill(TrackPhi[TrkID]);
                SelectionEnergy.at(file_no)->Fill(KineticEnergy[TrkID][2]/1000);
                SelectionMomentum.at(file_no)->Fill(TrackMomentum[TrkID]);

                SelXTrackStartEnd.at(file_no)->Fill(XTrackStart[TrkID]);
                SelXTrackStartEnd.at(file_no)->Fill(XTrackEnd[TrkID]);
                SelYTrackStartEnd.at(file_no)->Fill(YTrackStart[TrkID]);
                SelYTrackStartEnd.at(file_no)->Fill(YTrackEnd[TrkID]);
                SelZTrackStartEnd.at(file_no)->Fill(ZTrackStart[TrkID]);
                SelZTrackStartEnd.at(file_no)->Fill(ZTrackEnd[TrkID]);
                SelXVtxPosition.at(file_no)->Fill(XVertexPosition[VtxID]);
                SelYVtxPosition.at(file_no)->Fill(YVertexPosition[VtxID]);
                SelZVtxPosition.at(file_no)->Fill(ZVertexPosition[VtxID]);

                PhiVsTheta.at(file_no)->Fill(TrackPhi[TrkID],TrackTheta[TrkID]);
                PhiVsXPos.at(file_no)->Fill(TrackPhi[TrkID],XTrackStart[TrkID]);
                PhiVsXPos.at(file_no)->Fill(TrackPhi[TrkID],XTrackEnd[TrkID]);
                PhiVsYPos.at(file_no)->Fill(TrackPhi[TrkID],YTrackStart[TrkID]);
                PhiVsYPos.at(file_no)->Fill(TrackPhi[TrkID],YTrackEnd[TrkID]);
                PhiVsZPos.at(file_no)->Fill(TrackPhi[TrkID],ZTrackStart[TrkID]);
                PhiVsZPos.at(file_no)->Fill(TrackPhi[TrkID],ZTrackEnd[TrkID]);

                if(TrackPhi[TrkID] < 0) negPhi++;
                else posPhi++;

                if(TrackTheta[TrkID] > 0.8 && TrackTheta[TrkID] < 1.5 && TrackPhi[TrkID] > -2.0 && TrackPhi[TrkID] < -0.7)
                {
                    XPosVsYPos.at(file_no)->Fill(XTrackStart[TrkID],YTrackStart[TrkID]);
                    XPosVsYPos.at(file_no)->Fill(XTrackEnd[TrkID],YTrackEnd[TrkID]);
                    ZPosVsYPos.at(file_no)->Fill(ZTrackStart[TrkID],YTrackStart[TrkID]);
                    ZPosVsYPos.at(file_no)->Fill(ZTrackEnd[TrkID],YTrackEnd[TrkID]);
                    RangeVsYPos.at(file_no)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),YTrackStart[TrkID]);
                    RangeVsYPos.at(file_no)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),YTrackEnd[TrkID]);
                }

                // Fill Bgr
                if(file_no == 2 && CCNCFlag[0] == 0 && TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1)
                {
                    if(MCTrkID > -1 && PDGTruth[MCTrkID] == -13)
                    {
                        nubar++;
                        BgrTrackRange.at(0)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                        BgrTheta.at(0)->Fill(TrackTheta[TrkID]);
                        BgrCosTheta.at(0)->Fill(cos(TrackTheta[TrkID]));
                        BgrPhi.at(0)->Fill(TrackPhi[TrkID]);
                        BgrEnergy.at(0)->Fill(KineticEnergy[TrkID][2]/1000);
                        BgrMomentum.at(0)->Fill(TrackMomentum[TrkID]);
                        BgrXTrackStartEnd.at(0)->Fill(XTrackStart[TrkID]);
                        BgrXTrackStartEnd.at(0)->Fill(XTrackEnd[TrkID]);
                        BgrYTrackStartEnd.at(0)->Fill(YTrackStart[TrkID]);
                        BgrYTrackStartEnd.at(0)->Fill(YTrackEnd[TrkID]);
                        BgrZTrackStartEnd.at(0)->Fill(ZTrackStart[TrkID]);
                        BgrZTrackStartEnd.at(0)->Fill(ZTrackEnd[TrkID]);
                        BgrXVtxPosition.at(0)->Fill(XVertexPosition[VtxID]);
                        BgrYVtxPosition.at(0)->Fill(YVertexPosition[VtxID]);
                        BgrZVtxPosition.at(0)->Fill(ZVertexPosition[VtxID]);
                    }
                    else if(MCTrkID > -1 && abs(PDGTruth[MCTrkID]) == 11)
                    {
                        nue++;
                        BgrTrackRange.at(1)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                        BgrTheta.at(1)->Fill(TrackTheta[TrkID]);
                        BgrCosTheta.at(1)->Fill(cos(TrackTheta[TrkID]));
                        BgrPhi.at(1)->Fill(TrackPhi[TrkID]);
                        BgrEnergy.at(1)->Fill(KineticEnergy[TrkID][2]/1000);
                        BgrMomentum.at(1)->Fill(TrackMomentum[TrkID]);
                        BgrXTrackStartEnd.at(1)->Fill(XTrackStart[TrkID]);
                        BgrXTrackStartEnd.at(1)->Fill(XTrackEnd[TrkID]);
                        BgrYTrackStartEnd.at(1)->Fill(YTrackStart[TrkID]);
                        BgrYTrackStartEnd.at(1)->Fill(YTrackEnd[TrkID]);
                        BgrZTrackStartEnd.at(1)->Fill(ZTrackStart[TrkID]);
                        BgrZTrackStartEnd.at(1)->Fill(ZTrackEnd[TrkID]);
                        BgrXVtxPosition.at(1)->Fill(XVertexPosition[VtxID]);
                        BgrYVtxPosition.at(1)->Fill(YVertexPosition[VtxID]);
                        BgrZVtxPosition.at(1)->Fill(ZVertexPosition[VtxID]);
                    }
                    else if(!inFV(nuvtxx_truth[0],nuvtxy_truth[0],nuvtxz_truth[0]))
                    {
                        outFV++;
                        BgrTrackRange.at(2)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                        BgrTheta.at(2)->Fill(TrackTheta[TrkID]);
                        BgrCosTheta.at(2)->Fill(cos(TrackTheta[TrkID]));
                        BgrPhi.at(2)->Fill(TrackPhi[TrkID]);
                        BgrEnergy.at(2)->Fill(KineticEnergy[TrkID][2]/1000);
                        BgrMomentum.at(2)->Fill(TrackMomentum[TrkID]);
                        BgrXTrackStartEnd.at(2)->Fill(XTrackStart[TrkID]);
                        BgrXTrackStartEnd.at(2)->Fill(XTrackEnd[TrkID]);
                        BgrYTrackStartEnd.at(2)->Fill(YTrackStart[TrkID]);
                        BgrYTrackStartEnd.at(2)->Fill(YTrackEnd[TrkID]);
                        BgrZTrackStartEnd.at(2)->Fill(ZTrackStart[TrkID]);
                        BgrZTrackStartEnd.at(2)->Fill(ZTrackEnd[TrkID]);
                        BgrXVtxPosition.at(2)->Fill(XVertexPosition[VtxID]);
                        BgrYVtxPosition.at(2)->Fill(YVertexPosition[VtxID]);
                        BgrZVtxPosition.at(2)->Fill(ZVertexPosition[VtxID]);
                    }
                    else if(MCTrkID > -1 && file_no == 2 && TruthMode[0] == 0)
                    {
                        nuQE++;
                    }
                    else if(MCTrkID > -1 && file_no == 2 && TruthMode[0] == 1)
                    {
                        nuRES++;
                    }
                    else if(MCTrkID > -1 && file_no == 2 && TruthMode[0] == 2)
                    {
                        nuDIS++;
                    }
                    else if(MCTrkID > -1 && file_no == 2 && TruthMode[0] == 3)
                    {
                        nuCOH++;
                    }
                }
                else if(file_no == 2 && CCNCFlag[0] == 1 && TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1)
                {
                    NCnu++;
                    BgrTrackRange.at(3)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                    BgrTheta.at(3)->Fill(TrackTheta[TrkID]);
                    BgrCosTheta.at(3)->Fill(cos(TrackTheta[TrkID]));
                    BgrPhi.at(3)->Fill(TrackPhi[TrkID]);
                    BgrEnergy.at(3)->Fill(KineticEnergy[TrkID][2]/1000);
                    BgrMomentum.at(3)->Fill(TrackMomentum[TrkID]);
                    BgrXTrackStartEnd.at(3)->Fill(XTrackStart[TrkID]);
                    BgrXTrackStartEnd.at(3)->Fill(XTrackEnd[TrkID]);
                    BgrYTrackStartEnd.at(3)->Fill(YTrackStart[TrkID]);
                    BgrYTrackStartEnd.at(3)->Fill(YTrackEnd[TrkID]);
                    BgrZTrackStartEnd.at(3)->Fill(ZTrackStart[TrkID]);
                    BgrZTrackStartEnd.at(3)->Fill(ZTrackEnd[TrkID]);
                    BgrXVtxPosition.at(3)->Fill(XVertexPosition[VtxID]);
                    BgrYVtxPosition.at(3)->Fill(YVertexPosition[VtxID]);
                    BgrZVtxPosition.at(3)->Fill(ZVertexPosition[VtxID]);
                }
                else if(file_no == 2 && TrkOrigin[TrkID][TrkBestPlane[TrkID]] != 1)
                {
                    Cosmic++;
                    BgrTrackRange.at(4)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                    BgrTheta.at(4)->Fill(TrackTheta[TrkID]);
                    BgrCosTheta.at(4)->Fill(cos(TrackTheta[TrkID]));
                    BgrPhi.at(4)->Fill(TrackPhi[TrkID]);
                    BgrEnergy.at(4)->Fill(KineticEnergy[TrkID][2]/1000);
                    BgrMomentum.at(4)->Fill(TrackMomentum[TrkID]);
                    BgrXTrackStartEnd.at(4)->Fill(XTrackStart[TrkID]);
                    BgrXTrackStartEnd.at(4)->Fill(XTrackEnd[TrkID]);
                    BgrYTrackStartEnd.at(4)->Fill(YTrackStart[TrkID]);
                    BgrYTrackStartEnd.at(4)->Fill(YTrackEnd[TrkID]);
                    BgrZTrackStartEnd.at(4)->Fill(ZTrackStart[TrkID]);
                    BgrZTrackStartEnd.at(4)->Fill(ZTrackEnd[TrkID]);
                    BgrXVtxPosition.at(4)->Fill(XVertexPosition[VtxID]);
                    BgrYVtxPosition.at(4)->Fill(YVertexPosition[VtxID]);
                    BgrZVtxPosition.at(4)->Fill(ZVertexPosition[VtxID]);

                    if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] == -1)
                    {
                        UnknownOrigin++;
                    }
                }
            }
        }
        std::cout << Signal << " " << nubar << " " << nue << " " << outFV << " " << NCnu << " " << Cosmic << " " << UnknownOrigin << std::endl;
        std::cout << nuQE << " " << nuRES << " " << nuDIS << " " << nuCOH << std::endl;

        std::cout << "Number of negative phi in " << GenLabel.at(file_no) << " : " << negPhi << std::endl;
        std::cout << "Number of positive phi in " << GenLabel.at(file_no) << " : " << posPhi << std::endl;
        
        if(file_no == 0)
        {
            OnBeamFactor = (float)OnBeamRemain/(float)Signal;
        }
        else if(file_no == 1)
        {
            OffBeamFactor = (float)OffBeamRemain/(float)Signal;
        }
        
        std::cout << "OnBeam  : " << OnBeamFactor << " " << OnBeamRemain << std::endl;
        std::cout << "OffBeam  : " << OffBeamFactor << " " << OffBeamRemain << std::endl;

        ChainVec.at(file_no)->ResetBranchAddresses();
    }

    DataToLookAt.close();
    
    ScalingFactors.at(1) *= OnBeamFactor/OffBeamFactor;
    
    std::cout << SelectionTrackRange.at(2)->Integral() << " " << BgrTrackRange.at(4)->Integral() << std::endl;

    for(unsigned int bgrhist_no = 0; bgrhist_no < BgrLabel.size(); bgrhist_no++)
    {
        BgrTrackRange.at(bgrhist_no)->Scale(1/SelectionTrackRange.at(2)->Integral());
        BgrTheta.at(bgrhist_no)->Scale(1/SelectionTheta.at(2)->Integral());
        BgrCosTheta.at(bgrhist_no)->Scale(1/SelectionCosTheta.at(2)->Integral());
        BgrPhi.at(bgrhist_no)->Scale(1/SelectionPhi.at(2)->Integral());
        BgrEnergy.at(bgrhist_no)->Scale(1/SelectionEnergy.at(2)->Integral());
        BgrMomentum.at(bgrhist_no)->Scale(1/SelectionMomentum.at(2)->Integral());
        BgrXTrackStartEnd.at(bgrhist_no)->Scale(1/SelXTrackStartEnd.at(2)->Integral());
        BgrYTrackStartEnd.at(bgrhist_no)->Scale(1/SelYTrackStartEnd.at(2)->Integral());
        BgrZTrackStartEnd.at(bgrhist_no)->Scale(1/SelZTrackStartEnd.at(2)->Integral());
        BgrXVtxPosition.at(bgrhist_no)->Scale(1/SelXVtxPosition.at(2)->Integral());
        BgrYVtxPosition.at(bgrhist_no)->Scale(1/SelYVtxPosition.at(2)->Integral());
        BgrZVtxPosition.at(bgrhist_no)->Scale(1/SelZVtxPosition.at(2)->Integral());

        StackBgrTrackRange->Add(BgrTrackRange.at(bgrhist_no));
        StackBgrTheta->Add(BgrTheta.at(bgrhist_no));
        StackBgrCosTheta->Add(BgrCosTheta.at(bgrhist_no));
        StackBgrPhi->Add(BgrPhi.at(bgrhist_no));
        StackBgrEnergy->Add(BgrEnergy.at(bgrhist_no));
        StackBgrMomentum->Add(BgrMomentum.at(bgrhist_no));
        StackBgrXTrackStartEnd->Add(BgrXTrackStartEnd.at(bgrhist_no));
        StackBgrYTrackStartEnd->Add(BgrYTrackStartEnd.at(bgrhist_no));
        StackBgrZTrackStartEnd->Add(BgrZTrackStartEnd.at(bgrhist_no));
        StackBgrXVtxPosition->Add(BgrXVtxPosition.at(bgrhist_no));
        StackBgrYVtxPosition->Add(BgrYVtxPosition.at(bgrhist_no));
        StackBgrZVtxPosition->Add(BgrZVtxPosition.at(bgrhist_no));
    }

    for(unsigned int file_no = 0; file_no < ScalingFactors.size(); file_no++)
    {
        SelectionTrackRange.at(file_no)->Sumw2();
        SelectionTheta.at(file_no)->Sumw2();
        SelectionCosTheta.at(file_no)->Sumw2();
        SelectionPhi.at(file_no)->Sumw2();
        SelectionEnergy.at(file_no)->Sumw2();
        SelectionMomentum.at(file_no)->Sumw2();
        SelXTrackStartEnd.at(file_no)->Sumw2();
        SelYTrackStartEnd.at(file_no)->Sumw2();
        SelZTrackStartEnd.at(file_no)->Sumw2();
        SelXVtxPosition.at(file_no)->Sumw2();
        SelYVtxPosition.at(file_no)->Sumw2();
        SelZVtxPosition.at(file_no)->Sumw2();

        if(file_no < 2)
        {
            SelectionTrackRange.at(file_no)->Scale(ScalingFactors.at(file_no));
            SelectionTheta.at(file_no)->Scale(ScalingFactors.at(file_no));
            SelectionCosTheta.at(file_no)->Scale(ScalingFactors.at(file_no));
            SelectionPhi.at(file_no)->Scale(ScalingFactors.at(file_no));
            SelectionEnergy.at(file_no)->Scale(ScalingFactors.at(file_no));
            SelectionMomentum.at(file_no)->Scale(ScalingFactors.at(file_no));
            SelXTrackStartEnd.at(file_no)->Scale(ScalingFactors.at(file_no));
            SelYTrackStartEnd.at(file_no)->Scale(ScalingFactors.at(file_no));
            SelZTrackStartEnd.at(file_no)->Scale(ScalingFactors.at(file_no));
            SelXVtxPosition.at(file_no)->Scale(ScalingFactors.at(file_no));
            SelYVtxPosition.at(file_no)->Scale(ScalingFactors.at(file_no));
            SelZVtxPosition.at(file_no)->Scale(ScalingFactors.at(file_no));
        }

        PhiVsTheta.at(file_no)->Scale(ScalingFactors.at(file_no));
        PhiVsXPos.at(file_no)->Scale(ScalingFactors.at(file_no));
        PhiVsYPos.at(file_no)->Scale(ScalingFactors.at(file_no));
        PhiVsZPos.at(file_no)->Scale(ScalingFactors.at(file_no));
        RangeVsPE.at(file_no)->Scale(ScalingFactors.at(file_no));
        XPosVsYPos.at(file_no)->Scale(ScalingFactors.at(file_no));
        ZPosVsYPos.at(file_no)->Scale(ScalingFactors.at(file_no));
        RangeVsPE.at(file_no)->Scale(ScalingFactors.at(file_no));
        RangeVsYPos.at(file_no)->Scale(ScalingFactors.at(file_no));
        PhiVsFlashTrackDist.at(file_no)->Scale(ScalingFactors.at(file_no));
    }

    for(unsigned int hist_no = 0; hist_no < DataLabel.size(); hist_no++)
    {
        LegendData->AddEntry( SelectionTrackRange.at(hist_no), (DataLabel.at(hist_no)).c_str(),"l" );
    }

    TCanvas *Canvas1 = new TCanvas("Track Range of Selected Track", "Track Range of Selected Track", 1400, 1000);
    Canvas1->cd();
    SelectionTrackRange.at(0)->SetMaximum(1.1*GetMaximum(SelectionTrackRange));
    SelectionTrackRange.at(0)->Draw();
    SelectionTrackRange.at(1)->SetLineColor(2);
    SelectionTrackRange.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas1->SaveAs(("DataSelRange"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas2 = new TCanvas("Theta-Angle of Selected Track", "Theta-Angle of Selected Track", 1400, 1000);
    Canvas2->cd();
    SelectionTheta.at(0)->SetMaximum(1.1*GetMaximum(SelectionTheta));
    SelectionTheta.at(0)->Draw();
    SelectionTheta.at(1)->SetLineColor(2);
    SelectionTheta.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas2->SaveAs(("DataSelTheta"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas2a = new TCanvas("Cos Theta-Angle of Selected Track", "Cos Theta-Angle of Selected Track", 1400, 1000);
    Canvas2a->cd();
    SelectionCosTheta.at(0)->SetMaximum(1.3*GetMaximum(SelectionCosTheta));
    SelectionCosTheta.at(0)->Draw();
    SelectionCosTheta.at(1)->SetLineColor(2);
    SelectionCosTheta.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas2a->SaveAs(("DataSelCosTheta"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas3 = new TCanvas("Phi-Angle of Selected Track", "Phi-Angle of Selected Track", 1400, 1000);
    Canvas3->cd();
    SelectionPhi.at(0)->SetMaximum(1.1*GetMaximum(SelectionPhi));
    SelectionPhi.at(0)->SetMinimum(0);
    SelectionPhi.at(0)->Draw();
    SelectionPhi.at(1)->SetLineColor(2);
    SelectionPhi.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas3->SaveAs(("DataSelPhi"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas4 = new TCanvas("Energy of Selected Track", "Energy of Selected Track", 1400, 1000);
    Canvas4->cd();
    SelectionEnergy.at(0)->SetMaximum(1.1*GetMaximum(SelectionEnergy));
    SelectionEnergy.at(0)->Draw();
    SelectionEnergy.at(1)->SetLineColor(2);
    SelectionEnergy.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas4->SaveAs(("DataSelEnergy"+SelectionLabel+"."+FileType).c_str());
    
    TCanvas *Canvas4a = new TCanvas("Momentum of Selected Track", "Momentum of Selected Track", 1400, 1000);
    Canvas4a->cd();
    SelectionMomentum.at(0)->SetMaximum(1.1*GetMaximum(SelectionMomentum));
    SelectionMomentum.at(0)->Draw();
    SelectionMomentum.at(1)->SetLineColor(2);
    SelectionMomentum.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas4a->SaveAs(("DataSelMomentum"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas5 = new TCanvas("X Start & End Point Selected Track", "X Start & End Point Selected Track", 1400, 1000);
    Canvas5->cd();
    SelXTrackStartEnd.at(0)->SetMaximum(1.3*GetMaximum(SelXTrackStartEnd));
    SelXTrackStartEnd.at(0)->SetMinimum(0);
    SelXTrackStartEnd.at(0)->Draw();
    SelXTrackStartEnd.at(1)->SetLineColor(2);
    SelXTrackStartEnd.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas5->SaveAs(("DataSelXTrack"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas6 = new TCanvas("Y Start & End Point Selected Track", "Y Start & End Point Selected Track", 1400, 1000);
    Canvas6->cd();
    SelYTrackStartEnd.at(0)->SetMaximum(1.3*GetMaximum(SelYTrackStartEnd));
    SelYTrackStartEnd.at(0)->SetMinimum(0);
    SelYTrackStartEnd.at(0)->Draw();
    SelYTrackStartEnd.at(1)->SetLineColor(2);
    SelYTrackStartEnd.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas6->SaveAs(("DataSelYTrack"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas7 = new TCanvas("Z Start & End Point Selected Track", "Z Start & End Point Selected Track", 1400, 1000);
    Canvas7->cd();
    SelZTrackStartEnd.at(0)->SetMaximum(1.2*GetMaximum(SelZTrackStartEnd));
    SelZTrackStartEnd.at(0)->SetMinimum(0);
    SelZTrackStartEnd.at(0)->Draw();
    SelZTrackStartEnd.at(1)->SetLineColor(2);
    SelZTrackStartEnd.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas7->SaveAs(("DataSelZTrack"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas8 = new TCanvas("X Vertex Postion", "X Vertex Postion", 1400, 1000);
    Canvas8->cd();
    SelXVtxPosition.at(0)->SetMaximum(1.3*GetMaximum(SelXVtxPosition));
    SelXVtxPosition.at(0)->SetMinimum(0);
    SelXVtxPosition.at(0)->Draw();
    SelXVtxPosition.at(1)->SetLineColor(2);
    SelXVtxPosition.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas8->SaveAs(("DataSelXVertex"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas9 = new TCanvas("Y Vertex Postion", "Y Vertex Postion", 1400, 1000);
    Canvas9->cd();
    SelYVtxPosition.at(0)->SetMaximum(1.3*GetMaximum(SelYVtxPosition));
    SelYVtxPosition.at(0)->SetMinimum(0);
    SelYVtxPosition.at(0)->Draw();
    SelYVtxPosition.at(1)->SetLineColor(2);
    SelYVtxPosition.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas9->SaveAs(("DataSelYVertex"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas10 = new TCanvas("Z Vertex Postion", "Z Vertex Postion", 1400, 1000);
    Canvas10->cd();
    SelZVtxPosition.at(0)->SetMaximum(1.3*GetMaximum(SelZVtxPosition));
    SelZVtxPosition.at(0)->SetMinimum(0);
    SelZVtxPosition.at(0)->Draw();
    SelZVtxPosition.at(1)->SetLineColor(2);
    SelZVtxPosition.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas10->SaveAs(("DataSelZVertex"+SelectionLabel+"."+FileType).c_str());

    AddFirstTwoHistograms(SelectionTrackRange,-1.);
    AddFirstTwoHistograms(SelectionTheta,-1.);
    AddFirstTwoHistograms(SelectionCosTheta,-1.);
    AddFirstTwoHistograms(SelectionPhi,-1.);
    AddFirstTwoHistograms(SelectionEnergy,-1);
    AddFirstTwoHistograms(SelectionMomentum,-1);
    AddFirstTwoHistograms(SelXTrackStartEnd,-1);
    AddFirstTwoHistograms(SelYTrackStartEnd,-1);
    AddFirstTwoHistograms(SelZTrackStartEnd,-1);
    AddFirstTwoHistograms(SelXVtxPosition,-1);
    AddFirstTwoHistograms(SelYVtxPosition,-1);
    AddFirstTwoHistograms(SelZVtxPosition,-1);

    AddFirstTwoHistograms2D(PhiVsTheta,-1);
    AddFirstTwoHistograms2D(PhiVsXPos,-1);
    AddFirstTwoHistograms2D(PhiVsYPos,-1);
    AddFirstTwoHistograms2D(PhiVsZPos,-1);
    AddFirstTwoHistograms2D(RangeVsPE,-1);
    AddFirstTwoHistograms2D(XPosVsYPos,-1);
    AddFirstTwoHistograms2D(ZPosVsYPos,-1);
    AddFirstTwoHistograms2D(RangeVsYPos,-1);
    AddFirstTwoHistograms2D(PhiVsFlashTrackDist,-1);

    std::cout << "Theta KS significance: " << SelectionTheta.at(1)->KolmogorovTest(SelectionTheta.at(0)) << std::endl;
    std::cout << "Phi KS significance: " << SelectionPhi.at(1)->KolmogorovTest(SelectionPhi.at(0)) << std::endl;

    LegendMC->AddEntry( SelectionTrackRange.at(0), (MCLabel.at(0)).c_str(),"lep" );
    LegendMC->AddEntry( SelectionTrackRange.at(1), (MCLabel.at(1)).c_str(),"f" );
    for(int bgrhist_no = BgrLabel.size()-1; bgrhist_no > -1 ; bgrhist_no--)
    {
        LegendMC->AddEntry( BgrTrackRange.at(bgrhist_no), (BgrLabel.at(bgrhist_no)).c_str(),"f" );
    }

    TCanvas *Canvas11 = new TCanvas("OnBeam Minus OffBeam Track Range", "OnBeam Minus OffBeam Track Range", 1400, 1000);
    Canvas11->cd();
    SelectionTrackRange.at(1)->SetMaximum(1.4*SelectionTrackRange.at(1)->GetBinContent(SelectionTrackRange.at(1)->GetMaximumBin()));
    SelectionTrackRange.at(1)->SetMinimum(0.0);
    SelectionTrackRange.at(1)->SetFillColor(46);
    SelectionTrackRange.at(1)->DrawNormalized("E2");
    StackBgrTrackRange->Draw("SAME");
    SelectionTrackRange.at(0)->SetLineWidth(2);
    SelectionTrackRange.at(0)->SetLineColor(1);
    SelectionTrackRange.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas11->SaveAs(("On-OffBeamSelRange"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas12 = new TCanvas("OnBeam Minus OffBeam Theta-Angle", "OnBeam Minus OffBeam Theta-Angle", 1400, 1000);
    Canvas12->cd();
    SelectionTheta.at(1)->SetMaximum(1.5*SelectionTheta.at(1)->GetBinContent(SelectionTheta.at(1)->GetMaximumBin()));
    SelectionTheta.at(1)->SetMinimum(0.0);
    SelectionTheta.at(1)->SetFillColor(46);
    SelectionTheta.at(1)->DrawNormalized("E2");
    StackBgrTheta->Draw("SAME");
    SelectionTheta.at(0)->SetLineWidth(2);
    SelectionTheta.at(0)->SetLineColor(1);
    SelectionTheta.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas12->SaveAs(("On-OffBeamSelTheta"+SelectionLabel+"."+FileType).c_str());


    for(auto& BgrThetaHist :BgrTheta)
    {
        BgrThetaHist->Scale(SelectionTheta.at(1)->Integral());
    }

    for(auto& ThetaHistogram : SelectionTheta)
    {
        ThetaHistogram->Divide(SinTheta,1.);
    }

    for(auto& BgrThetaHist :BgrTheta)
    {
        BgrThetaHist->Divide(SinTheta,1.);
        BgrThetaHist->Scale(1/SelectionTheta.at(1)->Integral());
    }
    StackBgrTheta->Modified();

    TCanvas *Canvas12a = new TCanvas("OnBeam Minus OffBeam Theta-Angle Omega", "OnBeam Minus OffBeam Theta-Angle Omega", 1400, 1000);
    Canvas12a->cd();
    SelectionTheta.at(1)->SetMaximum(1.5*SelectionTheta.at(1)->GetBinContent(SelectionTheta.at(1)->GetMaximumBin()));
    SelectionTheta.at(1)->SetMinimum(0.0);
    SelectionTheta.at(1)->GetYaxis()->SetTitle("Weighted #frac{dn}{d#Omega}");
    SelectionTheta.at(1)->SetFillColor(46);
    SelectionTheta.at(1)->DrawNormalized("E2");
    StackBgrTheta->Draw("SAME");
    SelectionTheta.at(0)->SetLineWidth(2);
    SelectionTheta.at(0)->SetLineColor(1);
    SelectionTheta.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas12a->SaveAs(("On-OffBeamSelThetaOmega_Cuts"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas12b = new TCanvas("OnBeam Minus OffBeam Cos Theta-Angle", "OnBeam Minus OffBeam Cos Theta-Angle", 1400, 1000);
    Canvas12b->cd();
    SelectionCosTheta.at(1)->SetMaximum(1.5*SelectionCosTheta.at(1)->GetBinContent(SelectionCosTheta.at(1)->GetMaximumBin()));
    SelectionCosTheta.at(1)->SetMinimum(0.0);
    SelectionCosTheta.at(1)->SetFillColor(46);
    SelectionCosTheta.at(1)->DrawNormalized("E2");
    StackBgrCosTheta->Draw("SAME");
    SelectionCosTheta.at(0)->SetLineWidth(2);
    SelectionCosTheta.at(0)->SetLineColor(1);
    SelectionCosTheta.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas12b->SaveAs(("On-OffBeamSelCosTheta_Cuts"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas13 = new TCanvas("OnBeam Minus OffBeam Phi-Angle", "OnBeam Minus OffBeam Phi-Angle", 1400, 1000);
    Canvas13->cd();
    SelectionPhi.at(1)->SetMaximum(1.9*SelectionPhi.at(1)->GetBinContent(SelectionPhi.at(1)->GetMaximumBin()));
    SelectionPhi.at(1)->SetMinimum(0.0);
    SelectionPhi.at(1)->SetFillColor(46);
    SelectionPhi.at(1)->DrawNormalized("E2");
    StackBgrPhi->Draw("SAME");
    SelectionPhi.at(0)->SetLineWidth(2);
    SelectionPhi.at(0)->SetLineColor(1);
    SelectionPhi.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas13->SaveAs(("On-OffBeamSelPhi_Cuts"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas14 = new TCanvas("Energy", "Energy", 1400, 1000);
    Canvas14->cd();
    SelectionEnergy.at(1)->SetMaximum(1.2*SelectionEnergy.at(1)->GetBinContent(SelectionEnergy.at(1)->GetMaximumBin()));
    SelectionEnergy.at(1)->SetMinimum(0.0);
    SelectionEnergy.at(1)->SetFillColor(46);
    SelectionEnergy.at(1)->DrawNormalized("E2");
    StackBgrEnergy->Draw("SAME");
    SelectionEnergy.at(0)->SetLineWidth(2);
    SelectionEnergy.at(0)->SetLineColor(1);
    SelectionEnergy.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas14->SaveAs(("On-OffBeamSelEnergy_Cuts"+SelectionLabel+"."+FileType).c_str());
    
    TCanvas *Canvas14a = new TCanvas("Momentum", "Momentum", 1400, 1000);
    Canvas14a->cd();
    SelectionMomentum.at(1)->SetMaximum(1.2*SelectionMomentum.at(1)->GetBinContent(SelectionMomentum.at(1)->GetMaximumBin()));
    SelectionMomentum.at(1)->SetMinimum(0.0);
    SelectionMomentum.at(1)->SetFillColor(46);
    SelectionMomentum.at(1)->DrawNormalized("E2");
    StackBgrMomentum->Draw("SAME");
    SelectionMomentum.at(0)->SetLineWidth(2);
    SelectionMomentum.at(0)->SetLineColor(1);
    SelectionMomentum.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas14a->SaveAs(("On-OffBeamSelMomentum_Cuts"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas15 = new TCanvas("OnBeam Minus OffBeam X Start & End Point ", "OnBeam Minus OffBeam X Start & End Point ", 1400, 1000);
    Canvas15->cd();
    SelXTrackStartEnd.at(1)->SetMaximum(1.7*SelXTrackStartEnd.at(1)->GetBinContent(SelXTrackStartEnd.at(1)->GetMaximumBin()));
    SelXTrackStartEnd.at(1)->SetMinimum(0.0);
    SelXTrackStartEnd.at(1)->SetFillColor(46);
    SelXTrackStartEnd.at(1)->DrawNormalized("E2");
    StackBgrXTrackStartEnd->Draw("SAME");
    SelXTrackStartEnd.at(0)->SetLineWidth(2);
    SelXTrackStartEnd.at(0)->SetLineColor(1);
    SelXTrackStartEnd.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas15->SaveAs(("On-OffBeamSelXTrack_Cuts"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas16 = new TCanvas("OnBeam Minus OffBeam Y Start & End Point ", "OnBeam Minus OffBeam Y Start & End Point ", 1400, 1000);
    Canvas16->cd();
    SelYTrackStartEnd.at(1)->SetMaximum(1.7*SelYTrackStartEnd.at(1)->GetBinContent(SelYTrackStartEnd.at(1)->GetMaximumBin()));
    SelYTrackStartEnd.at(1)->SetMinimum(0.0);
    SelYTrackStartEnd.at(1)->SetFillColor(46);
    SelYTrackStartEnd.at(1)->DrawNormalized("E2");
    StackBgrYTrackStartEnd->Draw("SAME");
    SelYTrackStartEnd.at(0)->SetLineWidth(2);
    SelYTrackStartEnd.at(0)->SetLineColor(1);
    SelYTrackStartEnd.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas16->SaveAs(("On-OffBeamSelYTrack_Cuts"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas17 = new TCanvas("OnBeam Minus OffBeam Z Start & End Point ", "OnBeam Minus OffBeam Z Start & End Point ", 1400, 1000);
    Canvas17->cd();
    SelZTrackStartEnd.at(1)->SetMaximum(1.5*SelZTrackStartEnd.at(1)->GetBinContent(SelZTrackStartEnd.at(1)->GetMaximumBin()));
    SelZTrackStartEnd.at(1)->SetMinimum(0.0);
    SelZTrackStartEnd.at(1)->SetFillColor(46);
    SelZTrackStartEnd.at(1)->DrawNormalized("E2");
    StackBgrZTrackStartEnd->Draw("SAME");
    SelZTrackStartEnd.at(0)->SetLineWidth(2);
    SelZTrackStartEnd.at(0)->SetLineColor(1);
    SelZTrackStartEnd.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas17->SaveAs(("On-OffBeamSelZTrack_Cuts"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas18 = new TCanvas("OnBeam Minus OffBeam X Vertex Postion", "OnBeam Minus OffBeam X Vertex Postion", 1400, 1000);
    Canvas18->cd();
    SelXVtxPosition.at(1)->SetMaximum(1.7*SelXVtxPosition.at(1)->GetBinContent(SelXVtxPosition.at(1)->GetMaximumBin()));
    SelXVtxPosition.at(1)->SetMinimum(0.0);
    SelXVtxPosition.at(1)->SetFillColor(46);
    SelXVtxPosition.at(1)->DrawNormalized("E2");
    StackBgrXVtxPosition->Draw("SAME");
    SelXVtxPosition.at(0)->SetLineWidth(2);
    SelXVtxPosition.at(0)->SetLineColor(1);
    SelXVtxPosition.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas18->SaveAs(("On-OffBeamSelXVertex_Cuts"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas19 = new TCanvas("OnBeam Minus OffBeam Y Vertex Postion", "OnBeam Minus OffBeam Y Vertex Postion", 1400, 1000);
    Canvas19->cd();
    SelYVtxPosition.at(1)->SetMaximum(1.7*SelYVtxPosition.at(1)->GetBinContent(SelYVtxPosition.at(1)->GetMaximumBin()));
    SelYVtxPosition.at(1)->SetMinimum(0.0);
    SelYVtxPosition.at(1)->SetFillColor(46);
    SelYVtxPosition.at(1)->DrawNormalized("E2");
    StackBgrYVtxPosition->Draw("SAME");
    SelYVtxPosition.at(0)->SetLineWidth(2);
    SelYVtxPosition.at(0)->SetLineColor(1);
    SelYVtxPosition.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas19->SaveAs(("On-OffBeamSelYVertex_Cuts"+SelectionLabel+"."+FileType).c_str());

    TCanvas *Canvas20 = new TCanvas("OnBeam Minus OffBeam Z Vertex Postion", "OnBeam Minus OffBeam Z Vertex Postion", 1400, 1000);
    Canvas20->cd();
    SelZVtxPosition.at(1)->SetMaximum(1.7*SelZVtxPosition.at(1)->GetBinContent(SelZVtxPosition.at(1)->GetMaximumBin()));
    SelZVtxPosition.at(1)->SetMinimum(0.0);
    SelZVtxPosition.at(1)->SetFillColor(46);
    SelZVtxPosition.at(1)->DrawNormalized("E2");
    StackBgrZVtxPosition->Draw("SAME");
    SelZVtxPosition.at(0)->SetLineWidth(2);
    SelZVtxPosition.at(0)->SetLineColor(1);
    SelZVtxPosition.at(0)->DrawNormalized("SAME");
    LegendMC->Draw();
    Canvas20->SaveAs(("On-OffBeamSelZVertex_Cuts"+SelectionLabel+"."+FileType).c_str());

}

float GetMaximum(const std::vector<TH1F*>& HistVector)
{
    float Maximum = 0;

    for(unsigned int hist_no = 0; hist_no < 2; hist_no++)
    {
        float TempMax = HistVector.at(hist_no)->GetBinContent(HistVector.at(hist_no)->GetMaximumBin());

        if(TempMax > Maximum)
        {
            Maximum = TempMax;
        }
    }

    return Maximum;
}

void AddFirstTwoHistograms(std::vector<TH1F*>& HistVector, float Weight)
{
    if(HistVector.size() > 1)
    {
        HistVector.at(0)->Add(HistVector.at(1), Weight);
        delete HistVector.at(1);
        HistVector.erase(HistVector.begin()+1);
    }
    else
    {
        std::cout << "Histograms not added!" << std::endl;
    }
}

void AddFirstTwoHistograms2D(std::vector<TH2F*>& HistVector, float Weight)
{
    if(HistVector.size() > 1)
    {
        HistVector.at(0)->Add(HistVector.at(1), Weight);
        delete HistVector.at(1);
        HistVector.erase(HistVector.begin()+1);
    }
    else
    {
        std::cout << "Histograms not added!" << std::endl;
    }
}

float CalcLength(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2)
{
    return sqrt(pow(x_1-x_2, 2) + pow(y_1-y_2, 2) + pow(z_1-z_2, 2));
}

double FlashTrackDist(double flash, double start, double end)
{
    if(end >= start)
    {
        if(flash < end && flash > start) return 0;
        else return TMath::Min(fabs(flash-start), fabs(flash-end));
    }
    else
    {
        if(flash > end && flash < start) return 0;
        else return TMath::Min(fabs(flash-start), fabs(flash-end));
    }
}

bool inDeadRegion(double y, double z)
{
    if((y < (0.63*z+20)) && (y > (0.63*z-130))) return true;
    else if((y < (0.63*z-185)) && (y > (0.63*z-232.3))) return true;
    else if((y > (-0.63*z+429.3)) && (y < (-0.63*z+476.5))) return true;
    else if(z > 700 && z < 750) return true;
    else return false;
}

std::vector<TSpline5> Systematics()
{
    // Number of columns in file
    unsigned short NumberOfColumns = 5;

    // Initialize data structure
    std::vector <std::vector<float>> BeamSystematics;
    BeamSystematics.resize(NumberOfColumns);

    // Line and cell string for ifstream
    std::string FileLine;
    std::string Cell;

    // Open field systematic error file
    std::ifstream SysFile ("bnb_sys_error_uboone.txt");

    // check file
    if(SysFile.bad())
    {
        std::cout << "No such file or directory: " << "bnb_sys_error_uboone.txt" << std::endl;
        exit(-1);
    }

    // Loop over lines until files end
    while(std::getline(SysFile,FileLine))
    {
        // If not a header line
        if(FileLine[0] != 'E')
        {
            // Loop over all columns
            for(unsigned column_no = 0; column_no < NumberOfColumns; column_no++)
            {
                // Only read data if data stream works
                if(SysFile >> Cell)
                {
                    // Fill systematic error data into data structure
                    BeamSystematics.at(column_no).push_back(std::stof(Cell));
                }
            } // End of column loop
        } // if not header
    } // line loop

    // Initialize Graph vector
    std::vector<TGraph*> GraphVector;

    // Fill graphs with beam systematic data
    for(unsigned int entry_no = 1; entry_no < NumberOfColumns; entry_no++)
    {
        GraphVector.push_back( new TGraph(BeamSystematics.at(0).size(),BeamSystematics.at(0).data(),BeamSystematics.at(entry_no).data()) );
    }

    // Initialize spline vector
    std::vector<TSpline5> SplineVector;

    // Produce spline vector by fitting all graphs
    for(const auto& Graph : GraphVector)
    {
        SplineVector.push_back(TSpline5("",Graph));
        delete Graph;
    }

    return SplineVector;
}

void AdjustSysError(std::vector<TH1F*>& HistVector)
{
    for(unsigned int bin_no = 1; bin_no < HistVector.back()->GetNbinsX()+1; bin_no++)
    {
        HistVector.back()->SetBinError( bin_no, HistVector.back()->GetBinContent(bin_no) - HistVector.at(2)->GetBinContent(bin_no) + HistVector.at(2)->GetBinError(bin_no) );
        HistVector.back()->SetBinContent( bin_no, HistVector.at(2)->GetBinContent(bin_no) );
    }
}

bool inFV(double x, double y, double z) 
{
    if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)) return true;
    else return false;
}
