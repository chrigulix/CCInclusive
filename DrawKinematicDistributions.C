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

TSpline3* KEvsRSpline; // Global spline for momentum calculation

// Function which calculates the distance between two points
float CalcRange(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2);

// Function which checks if a point is in the FV
bool inFV(double x, double y, double z);

// Function which checks if a point is in the TPC
bool inTPC(double x, double y, double z);

// Add two histogramms with indices First and Last and weight
void AddHistograms(std::vector<TH1F*>& HistVector, unsigned int First, unsigned int Last, float Weight, bool EraseLast = false);

TH1F* AddToNewHist(std::vector<TH1F*>& HistVector, unsigned int First, unsigned int Last, float Weight);

// Subtract background histogram from selection histogram
void SubtractBgr(std::vector<TH1F*>& HistVector, std::vector<std::vector<TH1F*>>& BgrVector, unsigned int First, unsigned int Last, float Weight);

// Normalize Matrix by row
void NormMatrixByColumn(TH2F* UMatrix);

// Unsmearing of selected events
void SelectionUnsmearing(TH2F*& UMatrix, TH1F*& SVector);

// Momentum calculation
void MomentumSplinePreparation();

// Get Momentum
float GetMomentum(float TrackLength);

void CalcSigEfficiency(std::vector<TH1F*>& HistVector);

// Main Function
void DrawKinematicDistributions()
{
    float NumberOfTargets = (FVx - 2*borderx) * (FVy - 2*bordery) * (FVz - 2*borderz) * Density * Avogadro/ArMass*NoNucleons;

    // Fill momentum calculation spline
    MomentumSplinePreparation();

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

    // Selection Histogram Vectors
    std::vector<TH1F*> SelectionTrackRange;
    std::vector<TH1F*> SelectionCosTheta;
    std::vector<TH1F*> SelectionTheta;
    std::vector<TH1F*> SelectionPhi;
    std::vector<TH1F*> SelectionMomentum;
    std::vector<TH1F*> SelectionTrackLength;
    std::vector<TH1F*> SelXVtxPosition;
    std::vector<TH1F*> SelYVtxPosition;
    std::vector<TH1F*> SelZVtxPosition;

    // Background Histogram Vectors
    std::vector<std::vector<TH1F*>> BgrTrackRange;
    std::vector<std::vector<TH1F*>> BgrCosTheta;
    std::vector<std::vector<TH1F*>> BgrTheta;
    std::vector<std::vector<TH1F*>> BgrPhi;
    std::vector<std::vector<TH1F*>> BgrMomentum;
    std::vector<std::vector<TH1F*>> BgrTrackLength;
    std::vector<std::vector<TH1F*>> BgrXVtxPosition;
    std::vector<std::vector<TH1F*>> BgrYVtxPosition;
    std::vector<std::vector<TH1F*>> BgrZVtxPosition;

    // Produce beam systematics Histogram
    std::vector<std::deque<TH1F*>> TrackRangeBeamSys;
    std::vector<std::deque<TH1F*>> CosThetaBeamSys;
    std::vector<std::deque<TH1F*>> ThetaBeamSys;
    std::vector<std::deque<TH1F*>> PhiBeamSys;
    std::vector<std::deque<TH1F*>> MomentumBeamSys;
    std::vector<std::deque<TH1F*>> TrackLengthBeamSys;
    std::vector<std::deque<TH1F*>> XVtxPositionBeamSys;
    std::vector<std::deque<TH1F*>> YVtxPositionBeamSys;
    std::vector<std::deque<TH1F*>> ZVtxPositionBeamSys;

    // Efficiencies
    std::vector<TEfficiency*> EffTrackRange;
    std::vector<TEfficiency*> EffCosTheta;
    std::vector<TEfficiency*> EffTheta;
    std::vector<TEfficiency*> EffPhi;
    std::vector<TEfficiency*> EffMomentum;
    std::vector<TEfficiency*> EffXVtxPosition;
    std::vector<TEfficiency*> EffYVtxPosition;
    std::vector<TEfficiency*> EffZVtxPosition;

    // Purities
    std::vector<TEfficiency*> PurTrackRange;
    std::vector<TEfficiency*> PurCosTheta;
    std::vector<TEfficiency*> PurTheta;
    std::vector<TEfficiency*> PurPhi;
    std::vector<TEfficiency*> PurMomentum;
    std::vector<TEfficiency*> PurTrackLength;
    std::vector<TEfficiency*> PurXVtxPosition;
    std::vector<TEfficiency*> PurYVtxPosition;
    std::vector<TEfficiency*> PurZVtxPosition;

    // Unsemaring Matrix
    TH2F* UMatrixTrackRange;
    TH2F* UMatrixCosTheta;
    TH2F* UMatrixTheta;
    TH2F* UMatrixPhi;
    TH2F* UMatrixMomentum;
    TH2F* UMatrixXVtxPosition;
    TH2F* UMatrixYVtxPosition;
    TH2F* UMatrixZVtxPosition;

    size_t NumberOfBins = 20;

//     double MCPOT = 2.3e20/191362*92498;
//     double MCPOT = 2.304e20/141*62;
//     double TruthPOT = 5.451e19;
    double DataPOT = 4.950e19;
    double MCPOT = 2.300468e20; // 2.300468e20
    double MCMaPOT  = 2.000320e20; // 2.000320e20
    double MCMECPOT = 2.075599e20; // 2.075599e20
    double MCTEMPOT = 2.317981e20; // 2.317981e20


    // Read cosmic comparison histograms
    TFile* CosmicFile = new TFile((InputFolder+"/Cosmic_Distributions_Histograms_Mod.root").c_str(),"READ");

    // cd into cosmic file
    CosmicFile -> cd();

    // Cosmic histogram labels
    std::vector<std::string> CosmicHistLabels;
    CosmicHistLabels.push_back("Data Off-Beam BNBEXT All");
    CosmicHistLabels.push_back("In Time Corsika All");
    CosmicHistLabels.push_back("Cosmic Systematics All");

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

    // Read cosmic comparison histograms
    TFile* SelectionFile = new TFile((InputFolder+"/Selection_Histograms_Mod.root").c_str(),"READ");

    // cd into cosmic file
    SelectionFile -> cd();
//     SelectionFile -> ls();

    // Selection generator labels
    std::vector<std::pair<std::string,unsigned int>> GenLabel;

    // Scaling vector
    std::vector<float> ScalingFactors;

    GenLabel.push_back(std::make_pair("Data On-Beam BNB",1));
    ScalingFactors.push_back(1);

    GenLabel.push_back(std::make_pair("Data Off-Beam BNBEXT",9));
    ScalingFactors.push_back(1.2300);

    GenLabel.push_back(std::make_pair("MC Selection",46));
    ScalingFactors.push_back(DataPOT/MCPOT);

    GenLabel.push_back(std::make_pair("MA Adjusted Selection",28));
    ScalingFactors.push_back(DataPOT/MCMaPOT);

    GenLabel.push_back(std::make_pair("TEM Selection",30));
    ScalingFactors.push_back(DataPOT/MCTEMPOT);

    GenLabel.push_back(std::make_pair("MEC Selection",38));
    ScalingFactors.push_back(DataPOT/MCMECPOT);

    GenLabel.push_back(std::make_pair("MC Truth",1));
    ScalingFactors.push_back(DataPOT/MCPOT);

    // BEGIN READ --------------------------------------------------------------------------------------------------------------------------------------------

    // Fill selection histograms
    for(auto Label : GenLabel)
    {
        SelectionTrackRange.push_back( (TH1F*) SelectionFile->Get(("Track Range "+Label.first).c_str()) );
        SelectionCosTheta.push_back( (TH1F*) SelectionFile->Get(("cos#theta "+Label.first).c_str()) );
        SelectionTheta.push_back( (TH1F*) SelectionFile->Get(("#theta-Angle "+Label.first).c_str()) );
        SelectionPhi.push_back( (TH1F*) SelectionFile->Get(("#phi-Angle "+Label.first).c_str()) );
        SelectionMomentum.push_back( (TH1F*) SelectionFile->Get(("Momentum "+Label.first).c_str()) );
        SelectionTrackLength.push_back( (TH1F*) SelectionFile->Get(("Track Length "+Label.first).c_str()) );
        SelXVtxPosition.push_back( (TH1F*) SelectionFile->Get(("Vertex X position "+Label.first).c_str()) );
        SelYVtxPosition.push_back( (TH1F*) SelectionFile->Get(("Vertex Y position "+Label.first).c_str()) );
        SelZVtxPosition.push_back( (TH1F*) SelectionFile->Get(("Vertex Z position "+Label.first).c_str()) );

        // Set colour
        SelectionTrackRange.back()->SetFillColor(Label.second);
        SelectionTrackRange.back()->SetLineColor(Label.second);
        SelectionTrackRange.back()->SetMarkerColor(Label.second);

        SelectionCosTheta.back()->SetFillColor(Label.second);
        SelectionCosTheta.back()->SetLineColor(Label.second);
        SelectionCosTheta.back()->SetMarkerColor(Label.second);

        SelectionTheta.back()->SetFillColor(Label.second);
        SelectionTheta.back()->SetLineColor(Label.second);
        SelectionTheta.back()->SetMarkerColor(Label.second);

        SelectionPhi.back()->SetFillColor(Label.second);
        SelectionPhi.back()->SetLineColor(Label.second);
        SelectionPhi.back()->SetMarkerColor(Label.second);

        SelectionMomentum.back()->SetFillColor(Label.second);
        SelectionMomentum.back()->SetLineColor(Label.second);
        SelectionMomentum.back()->SetMarkerColor(Label.second);

        SelectionTrackLength.back()->SetFillColor(Label.second);
        SelectionTrackLength.back()->SetLineColor(Label.second);
        SelectionTrackLength.back()->SetMarkerColor(Label.second);

        SelXVtxPosition.back()->SetFillColor(Label.second);
        SelXVtxPosition.back()->SetLineColor(Label.second);
        SelXVtxPosition.back()->SetMarkerColor(Label.second);

        SelYVtxPosition.back()->SetFillColor(Label.second);
        SelYVtxPosition.back()->SetLineColor(Label.second);
        SelYVtxPosition.back()->SetMarkerColor(Label.second);

        SelZVtxPosition.back()->SetFillColor(Label.second);
        SelZVtxPosition.back()->SetLineColor(Label.second);
        SelZVtxPosition.back()->SetMarkerColor(Label.second);
    }

    SelectionTrackRange.at(0)->SetMarkerStyle(9);
    SelectionCosTheta.at(0)->SetMarkerStyle(9);
    SelectionTheta.at(0)->SetMarkerStyle(9);
    SelectionPhi.at(0)->SetMarkerStyle(9);
    SelectionMomentum.at(0)->SetMarkerStyle(9);
    SelectionTrackLength.at(0)->SetMarkerStyle(9);
    SelXVtxPosition.at(0)->SetMarkerStyle(9);
    SelYVtxPosition.at(0)->SetMarkerStyle(9);
    SelZVtxPosition.at(0)->SetMarkerStyle(9);

    // MC Background
    std::vector<std::pair<std::string,unsigned int>> BgrLabel;
    BgrLabel.push_back(std::make_pair("All",1));
    BgrLabel.push_back(std::make_pair("cosmic",38));
    BgrLabel.push_back(std::make_pair("dirt",28));
    BgrLabel.push_back(std::make_pair("outFV",42));
    BgrLabel.push_back(std::make_pair("anti nu_mu",kOrange-3));
    BgrLabel.push_back(std::make_pair("n_e-like",13));
    BgrLabel.push_back(std::make_pair("nu_NC",30));
    BgrLabel.push_back(std::make_pair("PureSelected",46));

//     std::vector<unsigned int> ColorMap = {13,28,42,30,38};

    BgrTrackRange.resize(4);
    BgrCosTheta.resize(4);
    BgrTheta.resize(4);
    BgrPhi.resize(4);
    BgrMomentum.resize(4);
    BgrTrackLength.resize(4);
    BgrXVtxPosition.resize(4);
    BgrYVtxPosition.resize(4);
    BgrZVtxPosition.resize(4);

    // Fill background histograms
    for(unsigned int file_no = 0; file_no < 4; file_no++)
    {
        for(auto Label : BgrLabel)
        {
            BgrTrackRange.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label.first+"Background Range "+std::to_string(file_no)).c_str()) );
            BgrCosTheta.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label.first+"Background cos#theta "+std::to_string(file_no)).c_str()) );
            BgrTheta.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label.first+"Background #theta "+std::to_string(file_no)).c_str()) );
            BgrPhi.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label.first+"Background #phi "+std::to_string(file_no)).c_str()) );
            BgrMomentum.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label.first+"Background Momentum "+std::to_string(file_no)).c_str()) );
            BgrTrackLength.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label.first+"Background Length "+std::to_string(file_no)).c_str()) );
            BgrXVtxPosition.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label.first+"Background XVtx "+std::to_string(file_no)).c_str()) );
            BgrYVtxPosition.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label.first+"Background YVtx "+std::to_string(file_no)).c_str()) );
            BgrZVtxPosition.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label.first+"Background ZVtx "+std::to_string(file_no)).c_str()) );

            // Set colours
            if(Label.second != BgrLabel.back().second)BgrTrackRange.at(file_no).back()->SetFillColor(Label.second);
            BgrTrackRange.at(file_no).back()->SetLineColor(Label.second);
            BgrTrackRange.at(file_no).back()->SetMarkerColor(Label.second);

            if(Label.second != BgrLabel.back().second)BgrCosTheta.at(file_no).back()->SetFillColor(Label.second);
            BgrCosTheta.at(file_no).back()->SetLineColor(Label.second);
            BgrCosTheta.at(file_no).back()->SetMarkerColor(Label.second);

            if(Label.second != BgrLabel.back().second)BgrTheta.at(file_no).back()->SetFillColor(Label.second);
            BgrTheta.at(file_no).back()->SetLineColor(Label.second);
            BgrTheta.at(file_no).back()->SetMarkerColor(Label.second);

            if(Label.second != BgrLabel.back().second)BgrPhi.at(file_no).back()->SetFillColor(Label.second);
            BgrPhi.at(file_no).back()->SetLineColor(Label.second);
            BgrPhi.at(file_no).back()->SetMarkerColor(Label.second);

            if(Label.second != BgrLabel.back().second)BgrMomentum.at(file_no).back()->SetFillColor(Label.second);
            BgrMomentum.at(file_no).back()->SetLineColor(Label.second);
            BgrMomentum.at(file_no).back()->SetMarkerColor(Label.second);

            if(Label.second != BgrLabel.back().second)BgrTrackLength.at(file_no).back()->SetFillColor(Label.second);
            BgrTrackLength.at(file_no).back()->SetLineColor(Label.second);
            BgrTrackLength.at(file_no).back()->SetMarkerColor(Label.second);

            if(Label.second != BgrLabel.back().second)BgrXVtxPosition.at(file_no).back()->SetFillColor(Label.second);
            if(Label.second != BgrLabel.back().second)BgrXVtxPosition.at(file_no).back()->SetLineColor(Label.second);
            if(Label.second != BgrLabel.back().second)BgrXVtxPosition.at(file_no).back()->SetMarkerColor(Label.second);

            if(Label.second != BgrLabel.back().second)BgrYVtxPosition.at(file_no).back()->SetFillColor(Label.second);
            if(Label.second != BgrLabel.back().second)BgrYVtxPosition.at(file_no).back()->SetLineColor(Label.second);
            if(Label.second != BgrLabel.back().second)BgrYVtxPosition.at(file_no).back()->SetMarkerColor(Label.second);

            if(Label.second != BgrLabel.back().second)BgrZVtxPosition.at(file_no).back()->SetFillColor(Label.second);
            if(Label.second != BgrLabel.back().second)BgrZVtxPosition.at(file_no).back()->SetLineColor(Label.second);
            if(Label.second != BgrLabel.back().second)BgrZVtxPosition.at(file_no).back()->SetMarkerColor(Label.second);
        }
    }

    // Systematic labels
    std::vector<std::string> SystLabel;
    SystLabel.push_back("dirt");
    SystLabel.push_back("outFV");
    SystLabel.push_back("anti nu_mu");
    SystLabel.push_back("n_e-like");
    SystLabel.push_back("nu_NC");
    SystLabel.push_back("PureSelected");

    // Initialize systematics vector
    TrackRangeBeamSys.resize(4);
    CosThetaBeamSys.resize(4);
    ThetaBeamSys.resize(4);
    PhiBeamSys.resize(4);
    MomentumBeamSys.resize(4);
    TrackLengthBeamSys.resize(4);
    XVtxPositionBeamSys.resize(4);
    YVtxPositionBeamSys.resize(4);
    ZVtxPositionBeamSys.resize(4);

    // Fill beam systematic histograms
    for(unsigned int file_no = 0; file_no < 4; file_no++)
    {
        for(auto Label : SystLabel)
        {
            TrackRangeBeamSys.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label+"Systematics Range "+std::to_string(file_no)).c_str()) );
            CosThetaBeamSys.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label+"Systematics cos#theta "+std::to_string(file_no)).c_str()) );
            ThetaBeamSys.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label+"Systematics #theta "+std::to_string(file_no)).c_str()) );
            PhiBeamSys.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label+"Systematics #phi "+std::to_string(file_no)).c_str()) );
            MomentumBeamSys.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label+"Systematics Momentum "+std::to_string(file_no)).c_str()) );
            TrackLengthBeamSys.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label+"Systematics Length "+std::to_string(file_no)).c_str()) );
            XVtxPositionBeamSys.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label+"Systematics XVtx "+std::to_string(file_no)).c_str()) );
            YVtxPositionBeamSys.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label+"Systematics YVtx "+std::to_string(file_no)).c_str()) );
            ZVtxPositionBeamSys.at(file_no).push_back( (TH1F*) SelectionFile->Get((Label+"Systematics ZVtx "+std::to_string(file_no)).c_str()) );
        }
    }

    // Fill smearing matrices
    UMatrixTrackRange = (TH2F*) SelectionFile->Get("UMatrixTrackRange");
    UMatrixCosTheta = (TH2F*) SelectionFile->Get("UMatrixCosTheta");
    UMatrixTheta = (TH2F*) SelectionFile->Get("UMatrixTheta");
    UMatrixPhi = (TH2F*) SelectionFile->Get("UMatrixPhi");
    UMatrixMomentum = (TH2F*) SelectionFile->Get("UMatrixMomentum");
    UMatrixXVtxPosition = (TH2F*) SelectionFile->Get("UMatrixXVtxPosition");
    UMatrixYVtxPosition = (TH2F*) SelectionFile->Get("UMatrixYVtxPosition");
    UMatrixZVtxPosition = (TH2F*) SelectionFile->Get("UMatrixZVtxPosition");

    // END READ ----------------------------------------------------------------------------------------------------------------------------------------------

    // BEGIN HISTOGRAM CALCULATIONS --------------------------------------------------------------------------------------------------------------------------

    // Scale histograms
    for(unsigned int scale_no = 0; scale_no < ScalingFactors.size(); scale_no++)
    {
        SelectionTrackRange.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionCosTheta.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionTheta.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionPhi.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionMomentum.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionTrackLength.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelXVtxPosition.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelYVtxPosition.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelZVtxPosition.at(scale_no)->Scale(ScalingFactors.at(scale_no));

        if(scale_no > 1 && scale_no < 6)
        {
            for(unsigned int bgr_no = 0; bgr_no < BgrLabel.size(); bgr_no++)
            {
                BgrTrackRange.at(scale_no-2).at(bgr_no) -> Scale(ScalingFactors.at(scale_no));
                BgrCosTheta.at(scale_no-2).at(bgr_no) -> Scale(ScalingFactors.at(scale_no));
                BgrTheta.at(scale_no-2).at(bgr_no) -> Scale(ScalingFactors.at(scale_no));
                BgrPhi.at(scale_no-2).at(bgr_no) -> Scale(ScalingFactors.at(scale_no));
                BgrMomentum.at(scale_no-2).at(bgr_no) -> Scale(ScalingFactors.at(scale_no));
                BgrTrackLength.at(scale_no-2).at(bgr_no) -> Scale(ScalingFactors.at(scale_no));
                BgrXVtxPosition.at(scale_no-2).at(bgr_no) -> Scale(ScalingFactors.at(scale_no));
                BgrYVtxPosition.at(scale_no-2).at(bgr_no) -> Scale(ScalingFactors.at(scale_no));
                BgrZVtxPosition.at(scale_no-2).at(bgr_no) -> Scale(ScalingFactors.at(scale_no));
            }

            BgrTrackRange.at(scale_no-2).back() -> Add(BgrTrackRange.at(scale_no-2).front());
            BgrCosTheta.at(scale_no-2).back() -> Add(BgrCosTheta.at(scale_no-2).front());
            BgrTheta.at(scale_no-2).back() -> Add(BgrTheta.at(scale_no-2).front());
            BgrPhi.at(scale_no-2).back() -> Add(BgrPhi.at(scale_no-2).front());
            BgrMomentum.at(scale_no-2).back() -> Add(BgrMomentum.at(scale_no-2).front());
            BgrTrackLength.at(scale_no-2).back() -> Add(BgrTrackLength.at(scale_no-2).front());
            BgrXVtxPosition.at(scale_no-2).back() -> Add(BgrXVtxPosition.at(scale_no-2).front());
            BgrYVtxPosition.at(scale_no-2).back() -> Add(BgrYVtxPosition.at(scale_no-2).front());
            BgrZVtxPosition.at(scale_no-2).back() -> Add(BgrZVtxPosition.at(scale_no-2).front());

            BgrTrackRange.at(scale_no-2).back() -> Add(SelectionTrackRange.at(1));
            BgrCosTheta.at(scale_no-2).back() -> Add(SelectionCosTheta.at(1));
            BgrTheta.at(scale_no-2).back() -> Add(SelectionTheta.at(1));
            BgrPhi.at(scale_no-2).back() -> Add(SelectionPhi.at(1));
            BgrMomentum.at(scale_no-2).back() -> Add(SelectionMomentum.at(1));
            BgrTrackLength.at(scale_no-2).back() -> Add(SelectionTrackLength.at(1));
            BgrXVtxPosition.at(scale_no-2).back() -> Add(SelXVtxPosition.at(1));
            BgrYVtxPosition.at(scale_no-2).back() -> Add(SelYVtxPosition.at(1));
            BgrZVtxPosition.at(scale_no-2).back() -> Add(SelZVtxPosition.at(1));

            for(unsigned int syst_no = 0; syst_no < SystLabel.size(); syst_no++)
            {
                TrackRangeBeamSys.at(scale_no-2).at(syst_no) -> Scale(std::pow(ScalingFactors.at(scale_no),2));
                CosThetaBeamSys.at(scale_no-2).at(syst_no) -> Scale(std::pow(ScalingFactors.at(scale_no),2));
                ThetaBeamSys.at(scale_no-2).at(syst_no) -> Scale(std::pow(ScalingFactors.at(scale_no),2));
                PhiBeamSys.at(scale_no-2).at(syst_no) -> Scale(std::pow(ScalingFactors.at(scale_no),2));
                MomentumBeamSys.at(scale_no-2).at(syst_no) -> Scale(std::pow(ScalingFactors.at(scale_no),2));
                TrackLengthBeamSys.at(scale_no-2).at(syst_no) -> Scale(std::pow(ScalingFactors.at(scale_no),2));
                XVtxPositionBeamSys.at(scale_no-2).at(syst_no) -> Scale(std::pow(ScalingFactors.at(scale_no),2));
                YVtxPositionBeamSys.at(scale_no-2).at(syst_no) -> Scale(std::pow(ScalingFactors.at(scale_no),2));
                ZVtxPositionBeamSys.at(scale_no-2).at(syst_no) -> Scale(std::pow(ScalingFactors.at(scale_no),2));
            }
        }

    }

    // Loop over all MC files
    for(unsigned int file_no = 0; file_no < 4; file_no++)
    {
        // Change dirt relative systematic uncertainties to 100%
        TrackRangeBeamSys.at(file_no).at(1) = (TH1F*)BgrTrackRange.at(file_no).at(2)->Clone();
        CosThetaBeamSys.at(file_no).at(1) = (TH1F*)BgrCosTheta.at(file_no).at(2)->Clone();
        ThetaBeamSys.at(file_no).at(1) = (TH1F*)BgrTheta.at(file_no).at(2)->Clone();
        PhiBeamSys.at(file_no).at(1) = (TH1F*)BgrPhi.at(file_no).at(2)->Clone();
        MomentumBeamSys.at(file_no).at(1) = (TH1F*)BgrMomentum.at(file_no).at(2)->Clone();
        TrackLengthBeamSys.at(file_no).at(1) = (TH1F*)BgrTrackLength.at(file_no).at(2)->Clone();
        XVtxPositionBeamSys.at(file_no).at(1) = (TH1F*)BgrXVtxPosition.at(file_no).at(2)->Clone();
        YVtxPositionBeamSys.at(file_no).at(1) = (TH1F*)BgrYVtxPosition.at(file_no).at(2)->Clone();
        ZVtxPositionBeamSys.at(file_no).at(1) = (TH1F*)BgrZVtxPosition.at(file_no).at(2)->Clone();
        // Square entries
        TrackRangeBeamSys.at(file_no).at(1) -> Multiply(BgrTrackRange.at(file_no).at(2));
        CosThetaBeamSys.at(file_no).at(1) -> Multiply(BgrCosTheta.at(file_no).at(2));
        ThetaBeamSys.at(file_no).at(1) -> Multiply(BgrTheta.at(file_no).at(2));
        PhiBeamSys.at(file_no).at(1) -> Multiply(BgrPhi.at(file_no).at(2));
        MomentumBeamSys.at(file_no).at(1) -> Multiply(BgrMomentum.at(file_no).at(2));
        TrackLengthBeamSys.at(file_no).at(1) -> Multiply(BgrTrackLength.at(file_no).at(2));
        XVtxPositionBeamSys.at(file_no).at(1) -> Multiply(BgrXVtxPosition.at(file_no).at(2));
        YVtxPositionBeamSys.at(file_no).at(1) -> Multiply(BgrYVtxPosition.at(file_no).at(2));
        ZVtxPositionBeamSys.at(file_no).at(1) -> Multiply(BgrZVtxPosition.at(file_no).at(2));

        // First clone cosmic relative variance
//         TrackRangeBeamSys.at(file_no).push_front( (TH1F*)CosmicTrackRange.back()->Clone() );
//         CosThetaBeamSys.at(file_no).push_front( (TH1F*)CosmicCosTheta.back()->Clone() );
//         ThetaBeamSys.at(file_no).push_front( (TH1F*)CosmicTheta.back()->Clone() );
//         PhiBeamSys.at(file_no).push_front( (TH1F*)CosmicPhi.back()->Clone() );
//         MomentumBeamSys.at(file_no).push_front( (TH1F*)CosmicMomentum.back()->Clone() );
//         TrackLengthBeamSys.at(file_no).push_front( (TH1F*)CosmicTrackLength.back()->Clone() );
//         XVtxPositionBeamSys.at(file_no).push_front( (TH1F*)CosmicXVtxPosition.back()->Clone() );
//         YVtxPositionBeamSys.at(file_no).push_front( (TH1F*)CosmicYVtxPosition.back()->Clone() );
//         ZVtxPositionBeamSys.at(file_no).push_front( (TH1F*)CosmicZVtxPosition.back()->Clone() );

        // First clone cosmic relative variance
        TrackRangeBeamSys.at(file_no).push_front( (TH1F*)BgrTrackRange.at(file_no).at(1)->Clone() );
        CosThetaBeamSys.at(file_no).push_front( (TH1F*)BgrCosTheta.at(file_no).at(1)->Clone() );
        ThetaBeamSys.at(file_no).push_front( (TH1F*)BgrTheta.at(file_no).at(1)->Clone() );
        PhiBeamSys.at(file_no).push_front( (TH1F*)BgrPhi.at(file_no).at(1)->Clone() );
        MomentumBeamSys.at(file_no).push_front( (TH1F*)BgrMomentum.at(file_no).at(1)->Clone() );
        TrackLengthBeamSys.at(file_no).push_front( (TH1F*)BgrTrackLength.at(file_no).at(1)->Clone() );
        XVtxPositionBeamSys.at(file_no).push_front( (TH1F*)BgrXVtxPosition.at(file_no).at(1)->Clone() );
        YVtxPositionBeamSys.at(file_no).push_front( (TH1F*)BgrYVtxPosition.at(file_no).at(1)->Clone() );
        ZVtxPositionBeamSys.at(file_no).push_front( (TH1F*)BgrZVtxPosition.at(file_no).at(1)->Clone() );

        // Multiply squared background to get variance
        TrackRangeBeamSys.at(file_no).front() -> Multiply(BgrTrackRange.at(file_no).at(1));
//         TrackRangeBeamSys.at(file_no).front() -> Multiply(BgrTrackRange.at(file_no).at(1));
        CosThetaBeamSys.at(file_no).front() -> Multiply(BgrCosTheta.at(file_no).at(1));
//         CosThetaBeamSys.at(file_no).front() -> Multiply(BgrCosTheta.at(file_no).at(1));
        ThetaBeamSys.at(file_no).front() -> Multiply(BgrTheta.at(file_no).at(1));
//         ThetaBeamSys.at(file_no).front() -> Multiply(BgrTheta.at(file_no).at(1));
        PhiBeamSys.at(file_no).front() -> Multiply(BgrPhi.at(file_no).at(1));
//         PhiBeamSys.at(file_no).front() -> Multiply(BgrPhi.at(file_no).at(1));
        MomentumBeamSys.at(file_no).front() -> Multiply(BgrMomentum.at(file_no).at(1));
//         MomentumBeamSys.at(file_no).front() -> Multiply(BgrMomentum.at(file_no).at(1));
        TrackLengthBeamSys.at(file_no).front() -> Multiply(BgrTrackLength.at(file_no).at(1));
//         TrackLengthBeamSys.at(file_no).front() -> Multiply(BgrTrackLength.at(file_no).at(1));
        XVtxPositionBeamSys.at(file_no).front() -> Multiply(BgrXVtxPosition.at(file_no).at(1));
//         XVtxPositionBeamSys.at(file_no).front() -> Multiply(BgrXVtxPosition.at(file_no).at(1));
        YVtxPositionBeamSys.at(file_no).front() -> Multiply(BgrYVtxPosition.at(file_no).at(1));
//         YVtxPositionBeamSys.at(file_no).front() -> Multiply(BgrYVtxPosition.at(file_no).at(1));
        ZVtxPositionBeamSys.at(file_no).front() -> Multiply(BgrZVtxPosition.at(file_no).at(1));
//         ZVtxPositionBeamSys.at(file_no).front() -> Multiply(BgrZVtxPosition.at(file_no).at(1));

        // Add all variances up
        TrackRangeBeamSys.at(file_no).push_front( (TH1F*)TrackRangeBeamSys.at(file_no).front()->Clone() );
        CosThetaBeamSys.at(file_no).push_front( (TH1F*)CosThetaBeamSys.at(file_no).front()->Clone() );
        ThetaBeamSys.at(file_no).push_front( (TH1F*)ThetaBeamSys.at(file_no).front()->Clone() );
        PhiBeamSys.at(file_no).push_front( (TH1F*)PhiBeamSys.at(file_no).front()->Clone() );
        MomentumBeamSys.at(file_no).push_front( (TH1F*)MomentumBeamSys.at(file_no).front()->Clone() );
        TrackLengthBeamSys.at(file_no).push_front( (TH1F*)TrackLengthBeamSys.at(file_no).front()->Clone() );
        XVtxPositionBeamSys.at(file_no).push_front( (TH1F*)XVtxPositionBeamSys.at(file_no).front()->Clone() );
        YVtxPositionBeamSys.at(file_no).push_front( (TH1F*)YVtxPositionBeamSys.at(file_no).front()->Clone() );
        ZVtxPositionBeamSys.at(file_no).push_front( (TH1F*)ZVtxPositionBeamSys.at(file_no).front()->Clone() );

        // Create a total systematic entry at the front of the vector
//         TrackRangeBeamSys.at(file_no).push_front(new TH1F(("Total Systematics Range "+std::to_string(file_no)).c_str(),"Range",NumberOfBins,0,700));
//         CosThetaBeamSys.at(file_no).push_front(new TH1F(("Total Systematics cos#theta "+std::to_string(file_no)).c_str(),"cos#theta",NumberOfBins,-1,1));
//         ThetaBeamSys.at(file_no).push_front(new TH1F(("Total Systematics #theta "+std::to_string(file_no)).c_str(),"#theta",NumberOfBins,0,180));
//         PhiBeamSys.at(file_no).push_front(new TH1F(("Total Systematics #phi "+std::to_string(file_no)).c_str(),"#varphi",NumberOfBins,-180,180));
//         MomentumBeamSys.at(file_no).push_front(new TH1F(("Total Systematics Momentum "+std::to_string(file_no)).c_str(),"Momentum",NumberOfBins,0,3));
//         TrackLengthBeamSys.at(file_no).push_front(new TH1F(("Total Systematics Length "+std::to_string(file_no)).c_str(),"Lenght",NumberOfBins,0,800));
//         XVtxPositionBeamSys.at(file_no).push_front(new TH1F(("Total Systematics XVtx "+std::to_string(file_no)).c_str(),"XVtx",NumberOfBins,0,256.35));
//         YVtxPositionBeamSys.at(file_no).push_front(new TH1F(("Total Systematics YVtx "+std::to_string(file_no)).c_str(),"YVtx",NumberOfBins,-233*0.5,233*0.5));
//         ZVtxPositionBeamSys.at(file_no).push_front(new TH1F(("Total Systematics ZVtx "+std::to_string(file_no)).c_str(),"ZVtx",NumberOfBins,0,1036.8));

        // Add all the variances
        for(unsigned int sys_no = 1; sys_no < TrackRangeBeamSys.size(); sys_no++)
        {
            TrackRangeBeamSys.at(file_no).front()->Add(TrackRangeBeamSys.at(file_no).at(sys_no));
            CosThetaBeamSys.at(file_no).front()->Add(CosThetaBeamSys.at(file_no).at(sys_no));
            ThetaBeamSys.at(file_no).front()->Add(ThetaBeamSys.at(file_no).at(sys_no));
            PhiBeamSys.at(file_no).front()->Add(PhiBeamSys.at(file_no).at(sys_no));
            MomentumBeamSys.at(file_no).front()->Add(MomentumBeamSys.at(file_no).at(sys_no));
            TrackLengthBeamSys.at(file_no).front()->Add(TrackLengthBeamSys.at(file_no).at(sys_no));
            XVtxPositionBeamSys.at(file_no).front()->Add(XVtxPositionBeamSys.at(file_no).at(sys_no));
            YVtxPositionBeamSys.at(file_no).front()->Add(YVtxPositionBeamSys.at(file_no).at(sys_no));
            ZVtxPositionBeamSys.at(file_no).front()->Add(ZVtxPositionBeamSys.at(file_no).at(sys_no));
        }

        // Add Off-Beam sample to MC prediction standard (for all MC predictions, put this into loop with file_no+2 index)
        SelectionTrackRange.at(file_no+2)-> Add(SelectionTrackRange.at(1));
        SelectionCosTheta.at(file_no+2)-> Add(SelectionCosTheta.at(1));
        SelectionTheta.at(file_no+2)-> Add(SelectionTheta.at(1));
        SelectionPhi.at(file_no+2)-> Add(SelectionPhi.at(1));
        SelectionMomentum.at(file_no+2)-> Add(SelectionMomentum.at(1));
        SelectionTrackLength.at(file_no+2)-> Add(SelectionTrackLength.at(1));
        SelXVtxPosition.at(file_no+2)-> Add(SelXVtxPosition.at(1));
        SelYVtxPosition.at(file_no+2)-> Add(SelYVtxPosition.at(1));
        SelZVtxPosition.at(file_no+2)-> Add(SelZVtxPosition.at(1));

        // Add systematics to error bars: sqrt(systematics^2 + statistics^2)
        for(unsigned int bin_no = 1; bin_no <= NumberOfBins; bin_no++)
        {
            SelectionTrackRange.at(file_no+2) -> SetBinError( bin_no, std::sqrt( std::pow(SelectionTrackRange.at(file_no+2)->GetBinError(bin_no),2) + TrackRangeBeamSys.at(file_no).front()->GetBinContent(bin_no) ) );

            SelectionCosTheta.at(file_no+2) -> SetBinError( bin_no, std::sqrt( std::pow(SelectionCosTheta.at(file_no+2)->GetBinError(bin_no),2) + CosThetaBeamSys.at(file_no).front()->GetBinContent(bin_no) ) );

            SelectionTheta.at(file_no+2) -> SetBinError( bin_no, std::sqrt( std::pow(SelectionTheta.at(file_no+2)->GetBinError(bin_no),2) + ThetaBeamSys.at(file_no).front()->GetBinContent(bin_no) ) );

            SelectionPhi.at(file_no+2) -> SetBinError( bin_no, std::sqrt( std::pow(SelectionPhi.at(file_no+2)->GetBinError(bin_no),2) + PhiBeamSys.at(file_no).front()->GetBinContent(bin_no) ) );

            SelectionMomentum.at(file_no+2) -> SetBinError( bin_no, std::sqrt( std::pow(SelectionMomentum.at(file_no+2)->GetBinError(bin_no),2) + MomentumBeamSys.at(file_no).front()->GetBinContent(bin_no) ) );

            SelectionTrackLength.at(file_no+2) -> SetBinError( bin_no, std::sqrt( std::pow(SelectionTrackLength.at(file_no+2)->GetBinError(bin_no),2) + TrackLengthBeamSys.at(file_no).front()->GetBinContent(bin_no) ) );

            SelXVtxPosition.at(file_no+2) -> SetBinError( bin_no, std::sqrt( std::pow(SelXVtxPosition.at(file_no+2)->GetBinError(bin_no),2) + XVtxPositionBeamSys.at(file_no).front()->GetBinContent(bin_no) ) );

            SelYVtxPosition.at(file_no+2) -> SetBinError( bin_no, std::sqrt( std::pow(SelYVtxPosition.at(file_no+2)->GetBinError(bin_no),2) + YVtxPositionBeamSys.at(file_no).front()->GetBinContent(bin_no) ) );

            SelZVtxPosition.at(file_no+2) -> SetBinError (bin_no, std::sqrt( std::pow(SelZVtxPosition.at(file_no+2)->GetBinError(bin_no),2) + ZVtxPositionBeamSys.at(file_no).front()->GetBinContent(bin_no) ) );
        }
    } // end MC-File loop
    
    // Test for chi2
    double MomentumChi2 = 0.0;
    double NDF = SelectionTrackRange.at(0)->Chi2Test(SelectionTrackRange.at(2),"CHI2 UW") / SelectionTrackRange.at(0)->Chi2Test(SelectionTrackRange.at(2),"CHI2/NDF UW");
    
    for(int bin = 1; bin <= SelectionCosTheta.at(0)->GetNbinsX(); bin++)
    {
        MomentumChi2 += std::pow( SelectionCosTheta.at(0)->GetBinContent(bin) - SelectionCosTheta.at(2)->GetBinContent(bin) ,2 ) /  ( std::pow(SelectionCosTheta.at(0)->GetBinError(bin),2) + std::pow(SelectionCosTheta.at(2)->GetBinError(bin),2));
    }
    
    std::cout << "Chi2 " << MomentumChi2 << std::endl;
    std::cout << "NDF " <<  NDF << std::endl;
    std::cout << "Chi2/NDF " << MomentumChi2/NDF << std::endl;
   

    // END HISTOGRAM CALCULATIONS ------------------------------------------------------------------------------------------------------------------------

    // BEGIN CHI-SQUARED ANALYSIS ------------------------------------------------------------------------------------------------------------------------

    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << " Distribution | \t MA = 0.99 \t MA = 1.35 \t TEM \t \t MEC" << std::endl;
    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << "Track Range |\t \t" << SelectionTrackRange.at(0)->Chi2Test(SelectionTrackRange.at(2),"CHI2/NDF") << "\t" << SelectionTrackRange.at(0)->Chi2Test(SelectionTrackRange.at(3),"CHI2/NDF") << "\t\t" << SelectionTrackRange.at(0)->Chi2Test(SelectionTrackRange.at(4),"CHI2/NDF") << "\t" << SelectionTrackRange.at(0)->Chi2Test(SelectionTrackRange.at(5),"CHI2/NDF") << std::endl;
    std::cout << "Cos Theta |\t \t" << SelectionCosTheta.at(0)->Chi2Test(SelectionCosTheta.at(2),"UW CHI2/NDF") << "\t" << SelectionCosTheta.at(0)->Chi2Test(SelectionCosTheta.at(3),"CHI2/NDF") << "\t" << SelectionCosTheta.at(0)->Chi2Test(SelectionCosTheta.at(4),"CHI2/NDF") << "\t" << SelectionCosTheta.at(0)->Chi2Test(SelectionCosTheta.at(5),"CHI2/NDF") << std::endl;
    std::cout << "Theta |\t \t \t" << SelectionTheta.at(0)->Chi2Test(SelectionTheta.at(2),"CHI2/NDF") << "\t" << SelectionTheta.at(0)->Chi2Test(SelectionTheta.at(3),"CHI2/NDF") << "\t" << SelectionTheta.at(0)->Chi2Test(SelectionTheta.at(4),"CHI2/NDF") << "\t\t" << SelectionTheta.at(0)->Chi2Test(SelectionTheta.at(5),"CHI2/NDF") << std::endl;
    std::cout << "Phi |\t \t \t" << SelectionPhi.at(0)->Chi2Test(SelectionPhi.at(2),"CHI2/NDF") << "\t\t" << SelectionPhi.at(0)->Chi2Test(SelectionPhi.at(3),"CHI2/NDF") << "\t" << SelectionPhi.at(0)->Chi2Test(SelectionPhi.at(4),"CHI2/NDF") << "\t" << SelectionPhi.at(0)->Chi2Test(SelectionPhi.at(5),"CHI2/NDF") << std::endl;
    std::cout << "Momentum |\t \t" << SelectionMomentum.at(0)->Chi2Test(SelectionMomentum.at(2),"CHI2/NDF") << "\t\t" << SelectionMomentum.at(0)->Chi2Test(SelectionMomentum.at(3),"CHI2/NDF") << "\t\t" << SelectionMomentum.at(0)->Chi2Test(SelectionMomentum.at(4),"CHI2/NDF") << "\t\t" << SelectionMomentum.at(0)->Chi2Test(SelectionMomentum.at(5),"CHI2/NDF") << std::endl;
    std::cout << "Track Length |\t \t" << SelectionTrackLength.at(0)->Chi2Test(SelectionTrackLength.at(2),"CHI2/NDF") << "\t" << SelectionTrackLength.at(0)->Chi2Test(SelectionTrackLength.at(3),"CHI2/NDF") << "\t" << SelectionTrackLength.at(0)->Chi2Test(SelectionTrackLength.at(4),"CHI2/NDF") << "\t" << SelectionTrackLength.at(0)->Chi2Test(SelectionTrackLength.at(5),"CHI2/NDF") << std::endl;
    std::cout << "X Vtx |\t \t \t" << SelXVtxPosition.at(0)->Chi2Test(SelXVtxPosition.at(2),"CHI2/NDF") << "\t\t" << SelXVtxPosition.at(0)->Chi2Test(SelXVtxPosition.at(3),"CHI2/NDF") << "\t\t" << SelXVtxPosition.at(0)->Chi2Test(SelXVtxPosition.at(4),"CHI2/NDF") << "\t\t" << SelXVtxPosition.at(0)->Chi2Test(SelXVtxPosition.at(5),"CHI2/NDF") << std::endl;
    std::cout << "Y Vtx |\t \t \t" << SelYVtxPosition.at(0)->Chi2Test(SelYVtxPosition.at(2),"CHI2/NDF") << "\t" << SelYVtxPosition.at(0)->Chi2Test(SelYVtxPosition.at(3),"CHI2/NDF") << "\t" << SelYVtxPosition.at(0)->Chi2Test(SelYVtxPosition.at(4),"CHI2/NDF") << "\t" << SelYVtxPosition.at(0)->Chi2Test(SelYVtxPosition.at(5),"CHI2/NDF") << std::endl;
    std::cout << "Z Vtx |\t \t \t" << SelZVtxPosition.at(0)->Chi2Test(SelZVtxPosition.at(2),"CHI2/NDF") << "\t\t" << SelZVtxPosition.at(0)->Chi2Test(SelZVtxPosition.at(3),"CHI2/NDF") << "\t\t" << SelZVtxPosition.at(0)->Chi2Test(SelZVtxPosition.at(4),"CHI2/NDF") << "\t\t" << SelZVtxPosition.at(0)->Chi2Test(SelZVtxPosition.at(5),"CHI2/NDF") << std::endl;
    std::cout << "----------------------------------------------------------------------------------" << std::endl;

    // END CHI-SQUARED ANALYSIS --------------------------------------------------------------------------------------------------------------------------

    std::cout << "Cosmic Syst Variance Range: " << TrackRangeBeamSys.at(0).at(1) -> Integral() << std::endl;
    std::cout << "Cosmic Syst Variance CosTheta: " << CosThetaBeamSys.at(0).at(1) -> Integral() << std::endl;
    std::cout << "Cosmic Syst Variance Theta: " << ThetaBeamSys.at(0).at(1) -> Integral() << std::endl;
    std::cout << "Cosmic Syst Variance phi: " << PhiBeamSys.at(0).at(1) -> Integral() << std::endl;
    std::cout << "Cosmic Syst Variance Mom: " << MomentumBeamSys.at(0).at(1) -> Integral() << std::endl;
    std::cout << "Cosmic Syst Variance Length: " << TrackLengthBeamSys.at(0).at(1) -> Integral() << std::endl;
    std::cout << "Cosmic Syst Variance X: " << XVtxPositionBeamSys.at(0).at(1) -> Integral() << std::endl;
    std::cout << "Cosmic Syst Variance Y: " << YVtxPositionBeamSys.at(0).at(1) -> Integral() << std::endl;
    std::cout << "Cosmic Syst Variance Z: " << ZVtxPositionBeamSys.at(0).at(1) -> Integral() << std::endl;

    // BEGIN BGR STACKING --------------------------------------------------------------------------------------------------------------------------------
    THStack* StackBgrTrackRange = new THStack("Bgr Track Range","Bgr Track Range");
    THStack* StackBgrCosTheta = new THStack("Bgr Cos Theta","Bgr Cos Theta");
    THStack* StackBgrTheta = new THStack("Bgr Theta","Bgr Theta");
    THStack* StackBgrPhi = new THStack("Bgr Phi","Bgr Phi");
    THStack* StackBgrMomentum = new THStack("Bgr Momentum","Bgr Momentum");
    THStack* StackBgrTrackLength = new THStack("Bgr Track Length","Bgr Track Length");
    THStack* StackBgrXVtxPosition = new THStack("Bgr X Vertex","Bgr X Vertex");
    THStack* StackBgrYVtxPosition = new THStack("Bgr X Vertex","Bgr X Vertex");
    THStack* StackBgrZVtxPosition = new THStack("Bgr X Vertex","Bgr X Vertex");

    // First entry is the BNB EXT Cosmic
    StackBgrTrackRange->Add(SelectionTrackRange.at(1));
    StackBgrCosTheta->Add(SelectionCosTheta.at(1));
    StackBgrTheta->Add(SelectionTheta.at(1));
    StackBgrPhi->Add(SelectionPhi.at(1));
    StackBgrMomentum->Add(SelectionMomentum.at(1));
    StackBgrTrackLength->Add(SelectionTrackLength.at(1));
    StackBgrXVtxPosition->Add(SelXVtxPosition.at(1));
    StackBgrYVtxPosition->Add(SelYVtxPosition.at(1));
    StackBgrZVtxPosition->Add(SelZVtxPosition.at(1));

    for(unsigned int bgr_no = 1; bgr_no < BgrLabel.size()-1; bgr_no++)
    {
        StackBgrTrackRange->Add(BgrTrackRange.at(0).at(bgr_no));
        StackBgrCosTheta->Add(BgrCosTheta.at(0).at(bgr_no));
        StackBgrTheta->Add(BgrTheta.at(0).at(bgr_no));
        StackBgrPhi->Add(BgrPhi.at(0).at(bgr_no));
        StackBgrMomentum->Add(BgrMomentum.at(0).at(bgr_no));
        StackBgrTrackLength->Add(BgrTrackLength.at(0).at(bgr_no));
        StackBgrXVtxPosition->Add(BgrXVtxPosition.at(0).at(bgr_no));
        StackBgrYVtxPosition->Add(BgrYVtxPosition.at(0).at(bgr_no));
        StackBgrZVtxPosition->Add(BgrZVtxPosition.at(0).at(bgr_no));
    }

    // END BGR STACKING ----------------------------------------------------------------------------------------------------------------------------------

    // BEGIN LEGENDS -------------------------------------------------------------------------------------------------------------------------------------

    // Forward folding with backgrounds legend
    TLegend* Legend = new TLegend(0.48,0.41,0.85,0.85);
    Legend->SetLineStyle ( 0 );
    Legend->SetLineColorAlpha ( 0,0 );
    Legend->SetFillStyle ( 0 );
    Legend->SetMargin ( 0.2 );
//     Legend->SetTextFont ( 43 );
//     Legend->SetTextSize ( 35 );
    Legend->SetHeader("Normalised to 4.95 #times 10^{19} POT","C");

    TLegendEntry *Header = (TLegendEntry*)Legend->GetListOfPrimitives()->First();
    Header->SetTextSize(0.05);

    std::vector<std::string> LegendLabel;
    LegendLabel.push_back("Off-Beam Cosmic (Data)");
    LegendLabel.push_back("Cosmic BGR (MC)");
    LegendLabel.push_back("Dirt BGR (MC)");
    LegendLabel.push_back("Out of FV BGR (MC)");
    LegendLabel.push_back("#bar{#nu}_{#mu} CC BGR (MC)");
    LegendLabel.push_back("#nu_{e} & #bar{#nu}_{e} CC BGR (MC)");
    LegendLabel.push_back("#nu_{x} NC BGR events (MC)");

    Legend->AddEntry( SelectionTrackRange.at(0), "On-Beam BNB (Data)","lep" );
    Legend->AddEntry( SelectionTrackRange.at(2), "Selected #nu_{#mu} CC + BGR (MC)","lf" );
    for(int bgr_no = LegendLabel.size()-1; bgr_no > 0 ; bgr_no--)
    {
        Legend->AddEntry( BgrTrackRange.at(0).at(bgr_no), (LegendLabel.at(bgr_no)).c_str(),"f" );
    }
    Legend->AddEntry( SelectionTrackRange.at(1), (LegendLabel.front()).c_str(),"f" );

    // ------------------------------------------------------------------------------------------------

    // Legend for different models
    TLegend* ModelLegend = new TLegend(0.48,0.60,0.85,0.85);
    ModelLegend->SetLineStyle ( 0 );
    ModelLegend->SetLineColorAlpha ( 0,0 );
    ModelLegend->SetFillStyle ( 0 );
    ModelLegend->SetMargin ( 0.2 );
    ModelLegend->SetHeader("Normalised to 4.95 #times 10^{19} POT","C");

    TLegendEntry *ModelHeader = (TLegendEntry*)ModelLegend->GetListOfPrimitives()->First();
    ModelHeader->SetTextSize(0.05);

    std::vector<std::string> ModelLegendLabel;
    ModelLegendLabel.push_back("MC M_{A} = 0.99 GeV");
    ModelLegendLabel.push_back("MC M_{A} = 1.35 GeV");
    ModelLegendLabel.push_back("MC ESF & TEM");
    ModelLegendLabel.push_back("MC MEC");

    ModelLegend->AddEntry( SelectionTrackRange.at(0), "Data On-Beam BNB","lep" );
    for(unsigned int bgr_no = 0; bgr_no < ModelLegendLabel.size() ; bgr_no++)
    {
        ModelLegend->AddEntry( SelectionTrackRange.at(bgr_no+2), (ModelLegendLabel.at(bgr_no)).c_str(),"lf" );
    }

    // END LEGENDS ---------------------------------------------------------------------------------------------------------------------------------------

    // BEGIN Draw SMEARING MATRICES ----------------------------------------------------------------------------------------------------------------------

    TCanvas *Smearing0 = new TCanvas("Smearing0", "Smearing0", 1100, 1000);
    Smearing0 -> cd();
    TBox* TrackRangeBack = new TBox(0, 0, 700, 700);
    TrackRangeBack -> SetFillColor(kBlue+3);
    UMatrixTrackRange -> SetTitle("Track Range Smearing Matrix");
    UMatrixTrackRange -> Draw("COLZ");
    TrackRangeBack -> Draw("SAME");
    UMatrixTrackRange -> Draw("COLZ SAME");
    gPad->RedrawAxis();
    Smearing0->SaveAs("SmearingMatrixTrackRange.pdf");

    TCanvas *Smearing1 = new TCanvas("Smearing1", "Smearing1", 1100, 1000);
    Smearing1 -> cd();
    TBox* CosThetaBack = new TBox(-1, -1, 1, 1);
    CosThetaBack -> SetFillColor(kBlue+3);
    UMatrixCosTheta -> SetTitle("cos(#theta) Smearing Matrix");
    UMatrixCosTheta -> Draw("COLZ");
    CosThetaBack -> Draw("SAME");
    UMatrixCosTheta -> Draw("COLZ SAME");
    gPad->RedrawAxis();
    Smearing1->SaveAs("SmearingMatrixCosTheta.pdf");

    TCanvas *Smearing2 = new TCanvas("Smearing2", "Smearing2", 1100, 1000);
    Smearing2 -> cd();
    TBox* ThetaBack = new TBox(0, 0, 180, 180);
    ThetaBack -> SetFillColor(kBlue+3);
    UMatrixTheta -> SetTitle("#theta-Angle Smearing Matrix");
    UMatrixTheta -> Draw("COLZ");
    ThetaBack -> Draw("SAME");
    UMatrixTheta -> Draw("COLZ SAME");
    gPad->RedrawAxis();
    Smearing2->SaveAs("SmearingMatrixTheta.pdf");

    TCanvas *Smearing3 = new TCanvas("Smearing3", "Smearing3", 1100, 1000);
    Smearing3 -> cd();
    TBox* PhiBack = new TBox(-180, -180, 180, 180);
    PhiBack -> SetFillColor(kBlue+3);
    UMatrixPhi -> SetTitle("#varphi-Angle Smearing Matrix");
    UMatrixPhi -> Draw("COLZ");
    PhiBack -> Draw("SAME");
    UMatrixPhi -> Draw("COLZ SAME");
    gPad->RedrawAxis();
    Smearing3->SaveAs("SmearingMatrixPhi.pdf");

    TCanvas *Smearing4 = new TCanvas("Smearing4", "Smearing4", 1100, 1000);
    Smearing4 -> cd();
    TBox* MomentumBack = new TBox(0, 0, 3, 3);
    MomentumBack -> SetFillColor(kBlue+3);
    UMatrixMomentum -> SetTitle("Momentum Smearing Matrix");
    UMatrixMomentum -> Draw("COLZ");
    MomentumBack -> Draw("SAME");
    UMatrixMomentum -> Draw("COLZ SAME");
    gPad->RedrawAxis();
    Smearing4->SaveAs("SmearingMatrixMomentum.pdf");

    TCanvas *Smearing5 = new TCanvas("Smearing5", "Smearing5", 1100, 1000);
    Smearing5 -> cd();
    TBox* XVtxPositionBack = new TBox(0, 0, 256.35, 256.35);
    XVtxPositionBack -> SetFillColor(kBlue+3);
    UMatrixXVtxPosition -> SetTitle("X-Vtx Smearing Matrix");
    UMatrixXVtxPosition -> Draw("COLZ");
    XVtxPositionBack -> Draw("SAME");
    UMatrixXVtxPosition -> Draw("COLZ SAME");
    gPad->RedrawAxis();
    Smearing5->SaveAs("SmearingMatrixXVtxPosition.pdf");

    TCanvas *Smearing6 = new TCanvas("Smearing6", "Smearing6", 1100, 1000);
    Smearing6 -> cd();
    TBox* YVtxPositionBack = new TBox(-233*0.5, -233*0.5, 233*0.5, 233*0.5);
    YVtxPositionBack -> SetFillColor(kBlue+3);
    UMatrixYVtxPosition -> SetTitle("Y-Vtx Smearing Matrix");
    UMatrixYVtxPosition -> Draw("COLZ");
    YVtxPositionBack -> Draw("SAME");
    UMatrixYVtxPosition -> Draw("COLZ SAME");
    gPad->RedrawAxis();
    Smearing6->SaveAs("SmearingMatrixYVtxPosition.pdf");

    TCanvas *Smearing7 = new TCanvas("Smearing7", "Smearing7", 1100, 1000);
    Smearing7 -> cd();
    TBox* ZVtxPositionBack = new TBox(0, 0, 1036.8, 1036.8);
    ZVtxPositionBack -> SetFillColor(kBlue+3);
    UMatrixZVtxPosition -> SetTitle("Z-Vtx Smearing Matrix");
    UMatrixZVtxPosition -> Draw("COLZ");
    ZVtxPositionBack -> Draw("SAME");
    UMatrixZVtxPosition -> Draw("COLZ SAME");
    gPad->RedrawAxis();
    Smearing7->SaveAs("SmearingMatrixZVtxPosition.pdf");

    // END Draw SMEARING MATRICES ̣------------------------------------------------------------------------------------------------------------------------

    // BEGIN DRAW FORWARD FOLDING ------------------------------------------------------------------------------------------------------------------------

    // Ratios
    std::vector<TH1F*> RatioTrackRange;
    std::vector<TH1F*> RatioCosTheta;
    std::vector<TH1F*> RatioTheta;
    std::vector<TH1F*> RatioPhi;
    std::vector<TH1F*> RatioMomentum;
    std::vector<TH1F*> RatioTrackLength;
    std::vector<TH1F*> RatioXVtxPosition;
    std::vector<TH1F*> RatioYVtxPosition;
    std::vector<TH1F*> RatioZVtxPosition;

    for(unsigned int model_no = 0; model_no < 4; model_no++)
    {
        RatioTrackRange.push_back( (TH1F*) SelectionTrackRange.at(0) -> Clone() );
        RatioTrackRange.back() -> Divide(SelectionTrackRange.at(2+model_no));
        RatioTrackRange.back() -> SetMarkerStyle(9);
        RatioTrackRange.back() -> SetMarkerColor(13);
        RatioTrackRange.back() -> SetLineColor(13);

        RatioCosTheta.push_back( (TH1F*) SelectionCosTheta.at(0) -> Clone() );
        RatioCosTheta.back() -> Divide(SelectionCosTheta.at(2+model_no));
        RatioCosTheta.back() -> SetMarkerStyle(9);
        RatioCosTheta.back() -> SetMarkerColor(13);
        RatioCosTheta.back() -> SetLineColor(13);

        RatioTheta.push_back( (TH1F*) SelectionTheta.at(0) -> Clone() );
        RatioTheta.back() -> Divide(SelectionTheta.at(2+model_no));
        RatioTheta.back() -> SetMarkerStyle(9);
        RatioTheta.back() -> SetMarkerColor(13);
        RatioTheta.back() -> SetLineColor(13);

        RatioPhi.push_back( (TH1F*) SelectionPhi.at(0) -> Clone() );
        RatioPhi.back() -> Divide(SelectionPhi.at(2+model_no));
        RatioPhi.back() -> SetMarkerStyle(9);
        RatioPhi.back() -> SetMarkerColor(13);
        RatioPhi.back() -> SetLineColor(13);

        RatioMomentum.push_back( (TH1F*) SelectionMomentum.at(0) -> Clone() );
        RatioMomentum.back() -> Divide(SelectionMomentum.at(2+model_no));
        RatioMomentum.back() -> SetMarkerStyle(9);
        RatioMomentum.back() -> SetMarkerColor(13);
        RatioMomentum.back() -> SetLineColor(13);

        RatioTrackLength.push_back( (TH1F*) SelectionTrackLength.at(0) -> Clone() );
        RatioTrackLength.back() -> Divide(SelectionTrackLength.at(2+model_no));
        RatioTrackLength.back() -> SetMarkerStyle(9);
        RatioTrackLength.back() -> SetMarkerColor(13);
        RatioTrackLength.back() -> SetLineColor(13);

        RatioXVtxPosition.push_back( (TH1F*) SelXVtxPosition.at(0) -> Clone() );
        RatioXVtxPosition.back() -> Divide(SelXVtxPosition.at(2+model_no));
        RatioXVtxPosition.back() -> SetMarkerStyle(9);
        RatioXVtxPosition.back() -> SetMarkerColor(13);
        RatioXVtxPosition.back() -> SetLineColor(13);

        RatioYVtxPosition.push_back( (TH1F*) SelYVtxPosition.at(0) -> Clone() );
        RatioYVtxPosition.back() -> Divide(SelYVtxPosition.at(2+model_no));
        RatioYVtxPosition.back() -> SetMarkerStyle(9);
        RatioYVtxPosition.back() -> SetMarkerColor(13);
        RatioYVtxPosition.back() -> SetLineColor(13);

        RatioZVtxPosition.push_back( (TH1F*) SelZVtxPosition.at(0) -> Clone() );
        RatioZVtxPosition.back() -> Divide(SelZVtxPosition.at(2+model_no));
        RatioZVtxPosition.back() -> SetMarkerStyle(9);
        RatioZVtxPosition.back() -> SetMarkerColor(13);
        RatioZVtxPosition.back() -> SetLineColor(13);
    }

    TCanvas *C0 = new TCanvas("C0", "C0", 1000, 1000);
    TPad *pad1TrackRange = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1TrackRange->SetBottomMargin(0);
    pad1TrackRange->Draw();
    pad1TrackRange->cd();
    SelectionTrackRange.at(2) -> GetYaxis() -> SetRangeUser(-20,900);
    SelectionTrackRange.at(2) -> SetFillColorAlpha(46,0.3);
    SelectionTrackRange.at(2) -> Draw("E2 SAME");
    BgrTrackRange.at(0).back() -> SetLineColor(46);
    BgrTrackRange.at(0).back() -> Draw("HIST ][ SAME");
    StackBgrTrackRange -> Draw("HIST SAME");
    SelectionTrackRange.at(0) -> SetLineColor(1);
    SelectionTrackRange.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C0->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2TrackRange = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2TrackRange->SetTopMargin(0);
    pad2TrackRange->SetBottomMargin(0.2);
    pad2TrackRange->Draw();
    pad2TrackRange->cd();
    RatioTrackRange.at(0) -> SetLineColor(13);
    RatioTrackRange.at(0) -> SetTitle("");
    RatioTrackRange.at(0) -> GetXaxis() -> SetLabelSize(0.09);
    RatioTrackRange.at(0) -> GetXaxis() -> SetTitleSize(0.1);
    RatioTrackRange.at(0) -> GetYaxis() -> SetRangeUser(-0.1,1.7);
    RatioTrackRange.at(0) -> GetYaxis() -> SetTitle("OnBeam/(MC + BGR)");
    RatioTrackRange.at(0) -> GetYaxis() -> SetLabelSize(0.09);
    RatioTrackRange.at(0) -> GetYaxis() -> SetTitleSize(0.1);
    RatioTrackRange.at(0) -> GetYaxis() -> SetTitleOffset(0.5);
    RatioTrackRange.at(0) -> Draw("SAME");
    TLine* Line1TrackRange = new TLine(0,1,700,1);
    Line1TrackRange -> SetLineStyle(7);
    Line1TrackRange -> Draw("SAME");
    C0->SaveAs("ForwardFoldedTrackRange.pdf");

    Legend -> SetX1NDC(0.15);
    Legend -> SetX2NDC(0.52);

    TCanvas *C1 = new TCanvas("C1", "C1", 1000, 1000);
    TPad *pad1CosTheta = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1CosTheta->SetBottomMargin(0);
    pad1CosTheta->Draw();
    pad1CosTheta->cd();
    SelectionCosTheta.at(2) -> GetYaxis() -> SetRangeUser(-50,1400);
    SelectionCosTheta.at(2) -> SetFillColorAlpha(46,0.3);
    SelectionCosTheta.at(2) -> Draw("E2 SAME");
    BgrCosTheta.at(0).back() -> SetLineColor(46);
    BgrCosTheta.at(0).back() -> Draw("HIST ][ SAME");
    StackBgrCosTheta -> Draw("HIST SAME");
    SelectionCosTheta.at(0) -> SetLineColor(1);
    SelectionCosTheta.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C1->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2CosTheta = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2CosTheta->SetTopMargin(0);
    pad2CosTheta->SetBottomMargin(0.2);
    pad2CosTheta->Draw();
    pad2CosTheta->cd();
    RatioCosTheta.at(0) -> SetLineColor(13);
    RatioCosTheta.at(0) -> SetTitle("");
    RatioCosTheta.at(0) -> GetXaxis() -> SetLabelSize(0.09);
    RatioCosTheta.at(0) -> GetXaxis() -> SetTitleSize(0.1);
//     RatioCosTheta.at(0) -> GetYaxis() -> SetRangeUser(-0.1,1.7);
    RatioCosTheta.at(0) -> GetYaxis() -> SetTitle("OnBeam/(MC + BGR)");
    RatioCosTheta.at(0) -> GetYaxis() -> SetLabelSize(0.09);
    RatioCosTheta.at(0) -> GetYaxis() -> SetTitleSize(0.1);
    RatioCosTheta.at(0) -> GetYaxis() -> SetTitleOffset(0.5);
    RatioCosTheta.at(0) -> Draw("SAME");
    TLine* Line1CosTheta = new TLine(-1,1,1,1);
    Line1CosTheta -> SetLineStyle(7);
    Line1CosTheta -> Draw("SAME");
    C1->SaveAs("ForwardFoldedCosTheta.pdf");

    Legend -> SetX1NDC(0.48);
    Legend -> SetX2NDC(0.85);

    TCanvas *C2 = new TCanvas("C2", "C2", 1000, 1000);
    TPad *pad1Theta = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1Theta->SetBottomMargin(0);
    pad1Theta->Draw();
    pad1Theta->cd();
    SelectionTheta.at(2) -> GetYaxis() -> SetRangeUser(-20,700);
    SelectionTheta.at(2) -> SetFillColorAlpha(46,0.3);
    SelectionTheta.at(2) -> Draw("E2 SAME");
    BgrTheta.at(0).back() -> SetLineColor(46);
    BgrTheta.at(0).back() -> Draw("HIST ][ SAME");
    StackBgrTheta -> Draw("HIST SAME");
    SelectionTheta.at(0) -> SetLineColor(1);
    SelectionTheta.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C2->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2Theta = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2Theta->SetTopMargin(0);
    pad2Theta->SetBottomMargin(0.2);
    pad2Theta->Draw();
    pad2Theta->cd();
    RatioTheta.at(0) -> SetLineColor(13);
    RatioTheta.at(0) -> SetTitle("");
    RatioTheta.at(0) -> GetXaxis() -> SetLabelSize(0.09);
    RatioTheta.at(0) -> GetXaxis() -> SetTitleSize(0.1);
//     RatioTheta.at(0) -> GetYaxis() -> SetRangeUser(-0.1,1.7);
    RatioTheta.at(0) -> GetYaxis() -> SetTitle("OnBeam/(MC + BGR)");
    RatioTheta.at(0) -> GetYaxis() -> SetLabelSize(0.09);
    RatioTheta.at(0) -> GetYaxis() -> SetTitleSize(0.1);
    RatioTheta.at(0) -> GetYaxis() -> SetTitleOffset(0.5);
    RatioTheta.at(0) -> Draw("SAME");
    TLine* Line1Theta = new TLine(0,1,180,1);
    Line1Theta -> SetLineStyle(7);
    Line1Theta -> Draw("SAME");
    C2->SaveAs("ForwardFoldedTheta.pdf");

//     Legend -> SetX1NDC(0.40);
//     Legend -> SetX2NDC(0.77);

    TCanvas *C3 = new TCanvas("C3", "C3", 1000, 1000);
    TPad *pad1Phi = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1Phi->SetBottomMargin(0);
    pad1Phi->Draw();
    pad1Phi->cd();
    SelectionPhi.at(2) -> GetYaxis() -> SetRangeUser(-10,400);
    SelectionPhi.at(2) -> SetFillColorAlpha(46,0.3);
    SelectionPhi.at(2) -> Draw("E2 SAME");
    BgrPhi.at(0).back() -> SetLineColor(46);
    BgrPhi.at(0).back() -> Draw("HIST ][ SAME");
    StackBgrPhi -> Draw("HIST SAME");
    SelectionPhi.at(0) -> SetLineColor(1);
    SelectionPhi.at(0) -> Draw("SAME");
//     Legend -> Draw();
    gPad->RedrawAxis();
    C3->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2Phi = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2Phi->SetTopMargin(0);
    pad2Phi->SetBottomMargin(0.2);
    pad2Phi->Draw();
    pad2Phi->cd();
    RatioPhi.at(0) -> SetLineColor(13);
    RatioPhi.at(0) -> SetTitle("");
    RatioPhi.at(0) -> GetXaxis() -> SetLabelSize(0.09);
    RatioPhi.at(0) -> GetXaxis() -> SetTitleSize(0.1);
//     RatioPhi.at(0) -> GetYaxis() -> SetRangeUser(-0.1,1.7);
    RatioPhi.at(0) -> GetYaxis() -> SetTitle("OnBeam/(MC + BGR)");
    RatioPhi.at(0) -> GetYaxis() -> SetLabelSize(0.09);
    RatioPhi.at(0) -> GetYaxis() -> SetTitleSize(0.1);
    RatioPhi.at(0) -> GetYaxis() -> SetTitleOffset(0.5);
    RatioPhi.at(0) -> Draw("SAME");
    TLine* Line1Phi = new TLine(-180,1,180,1);
    Line1Phi -> SetLineStyle(7);
    Line1Phi -> Draw("SAME");
    gPad->RedrawAxis();
    C3->SaveAs("ForwardFoldedPhi.pdf");

    TCanvas *C4 = new TCanvas("C4", "C4", 1000, 1000);
    TPad *pad1Momentum = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1Momentum->SetBottomMargin(0);
    pad1Momentum->Draw();
    pad1Momentum->cd();
    SelectionMomentum.at(2) -> GetYaxis() -> SetRangeUser(-50,1400);
    SelectionMomentum.at(2) -> SetFillColorAlpha(46,0.3);
    SelectionMomentum.at(2) -> Draw("E2 SAME");
    BgrMomentum.at(0).back() -> SetLineColor(46);
    BgrMomentum.at(0).back() -> Draw("HIST ][ SAME");
    StackBgrMomentum -> Draw("HIST SAME");
    SelectionMomentum.at(0) -> SetLineColor(1);
    SelectionMomentum.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C4->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2Momentum = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2Momentum->SetTopMargin(0);
    pad2Momentum->SetBottomMargin(0.2);
    pad2Momentum->Draw();
    pad2Momentum->cd();
    RatioMomentum.at(0) -> SetLineColor(13);
    RatioMomentum.at(0) -> SetTitle("");
    RatioMomentum.at(0) -> GetXaxis() -> SetLabelSize(0.09);
    RatioMomentum.at(0) -> GetXaxis() -> SetTitleSize(0.1);
//     RatioMomentum.at(0) -> GetYaxis() -> SetRangeUser(-0.1,1.7);
    RatioMomentum.at(0) -> GetYaxis() -> SetTitle("OnBeam/(MC + BGR)");
    RatioMomentum.at(0) -> GetYaxis() -> SetLabelSize(0.09);
    RatioMomentum.at(0) -> GetYaxis() -> SetTitleSize(0.1);
    RatioMomentum.at(0) -> GetYaxis() -> SetTitleOffset(0.5);
    RatioMomentum.at(0) -> Draw("SAME");
    TLine* Line1Momentum = new TLine(0,1,3,1);
    Line1Momentum -> SetLineStyle(7);
    Line1Momentum -> Draw("SAME");
    gPad->RedrawAxis();
    C4->SaveAs("ForwardFoldedMomentum.pdf");

    TCanvas *C5 = new TCanvas("C5", "C5", 1000, 1000);
    TPad *pad1TrackLength = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1TrackLength->SetBottomMargin(0);
    pad1TrackLength->Draw();
    pad1TrackLength->cd();
    SelectionTrackLength.at(2) -> GetYaxis() -> SetRangeUser(-50,1000);
    SelectionTrackLength.at(2) -> SetFillColorAlpha(46,0.3);
    SelectionTrackLength.at(2) -> Draw("E2 SAME");
    BgrTrackLength.at(0).back() -> SetLineColor(46);
    BgrTrackLength.at(0).back() -> Draw("HIST ][ SAME");
    StackBgrTrackLength -> Draw("HIST SAME");
    SelectionTrackLength.at(0) -> SetLineColor(1);
    SelectionTrackLength.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C5->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2TrackLength = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2TrackLength->SetTopMargin(0);
    pad2TrackLength->SetBottomMargin(0.2);
    pad2TrackLength->Draw();
    pad2TrackLength->cd();
    RatioTrackLength.at(0) -> SetLineColor(13);
    RatioTrackLength.at(0) -> SetTitle("");
    RatioTrackLength.at(0) -> GetXaxis() -> SetLabelSize(0.09);
    RatioTrackLength.at(0) -> GetXaxis() -> SetTitleSize(0.1);
//     RatioTrackLength.at(0) -> GetYaxis() -> SetRangeUser(-0.1,1.7);
    RatioTrackLength.at(0) -> GetYaxis() -> SetTitle("OnBeam/(MC + BGR)");
    RatioTrackLength.at(0) -> GetYaxis() -> SetLabelSize(0.09);
    RatioTrackLength.at(0) -> GetYaxis() -> SetTitleSize(0.1);
    RatioTrackLength.at(0) -> GetYaxis() -> SetTitleOffset(0.5);
    RatioTrackLength.at(0) -> Draw("SAME");
    TLine* Line1TrackLength = new TLine(0,1,800,1);
    Line1TrackLength -> SetLineStyle(7);
    Line1TrackLength -> Draw("SAME");
    gPad->RedrawAxis();
    C5->SaveAs("ForwardFoldedTrackLength.pdf");

    Legend -> SetNColumns(2);
    Legend -> SetX1NDC(0.13);
    Legend -> SetX2NDC(0.87);
    Legend -> SetY1NDC(0.62);
    Legend -> SetY2NDC(0.87);

    TCanvas *C6 = new TCanvas("C6", "C6", 1000, 1000);
    TPad *pad1XVtxPosition = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1XVtxPosition->SetBottomMargin(0);
    pad1XVtxPosition->Draw();
    pad1XVtxPosition->cd();
    SelXVtxPosition.at(2) -> GetYaxis() -> SetRangeUser(-10,400);
    SelXVtxPosition.at(2) -> SetFillColorAlpha(46,0.3);
    SelXVtxPosition.at(2) -> Draw("E2 SAME");
    BgrXVtxPosition.at(0).back() -> SetLineColor(46);
    BgrXVtxPosition.at(0).back() -> Draw("HIST ][ SAME");
    StackBgrXVtxPosition -> Draw("HIST SAME");
    SelXVtxPosition.at(0) -> SetLineColor(1);
    SelXVtxPosition.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C6->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2XVtxPosition = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2XVtxPosition->SetTopMargin(0);
    pad2XVtxPosition->SetBottomMargin(0.2);
    pad2XVtxPosition->Draw();
    pad2XVtxPosition->cd();
    RatioXVtxPosition.at(0) -> SetLineColor(13);
    RatioXVtxPosition.at(0) -> SetTitle("");
    RatioXVtxPosition.at(0) -> GetXaxis() -> SetLabelSize(0.09);
    RatioXVtxPosition.at(0) -> GetXaxis() -> SetTitleSize(0.1);
//     RatioXVtxPosition.at(0) -> GetYaxis() -> SetRangeUser(-0.1,1.7);
    RatioXVtxPosition.at(0) -> GetYaxis() -> SetTitle("OnBeam/(MC + BGR)");
    RatioXVtxPosition.at(0) -> GetYaxis() -> SetLabelSize(0.09);
    RatioXVtxPosition.at(0) -> GetYaxis() -> SetTitleSize(0.1);
    RatioXVtxPosition.at(0) -> GetYaxis() -> SetTitleOffset(0.5);
    RatioXVtxPosition.at(0) -> Draw("SAME");
    TLine* Line1XVtxPosition = new TLine(0,1,256.35,1);
    Line1XVtxPosition -> SetLineStyle(7);
    Line1XVtxPosition -> Draw("SAME");
    gPad->RedrawAxis();
    C6->SaveAs("ForwardFoldedXVtxPosition.pdf");

    TCanvas *C7 = new TCanvas("C7", "C7", 1000, 1000);
    TPad *pad1YVtxPosition = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1YVtxPosition->SetBottomMargin(0);
    pad1YVtxPosition->Draw();
    pad1YVtxPosition->cd();
    SelYVtxPosition.at(2) -> GetYaxis() -> SetRangeUser(-10,400);
    SelYVtxPosition.at(2) -> SetFillColorAlpha(46,0.3);
    SelYVtxPosition.at(2) -> Draw("E2 SAME");
    BgrYVtxPosition.at(0).back() -> SetLineColor(46);
    BgrYVtxPosition.at(0).back() -> Draw("HIST ][ SAME");
    StackBgrYVtxPosition -> Draw("HIST SAME");
    SelYVtxPosition.at(0) -> SetLineColor(1);
    SelYVtxPosition.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C7->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2YVtxPosition = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2YVtxPosition->SetTopMargin(0);
    pad2YVtxPosition->SetBottomMargin(0.2);
    pad2YVtxPosition->Draw();
    pad2YVtxPosition->cd();
    RatioYVtxPosition.at(0) -> SetLineColor(13);
    RatioYVtxPosition.at(0) -> SetTitle("");
    RatioYVtxPosition.at(0) -> GetXaxis() -> SetLabelSize(0.09);
    RatioYVtxPosition.at(0) -> GetXaxis() -> SetTitleSize(0.1);
    RatioYVtxPosition.at(0) -> GetYaxis() -> SetRangeUser(0.7,1.45);
    RatioYVtxPosition.at(0) -> GetYaxis() -> SetTitle("OnBeam/(MC + BGR)");
    RatioYVtxPosition.at(0) -> GetYaxis() -> SetLabelSize(0.09);
    RatioYVtxPosition.at(0) -> GetYaxis() -> SetTitleSize(0.1);
    RatioYVtxPosition.at(0) -> GetYaxis() -> SetTitleOffset(0.5);
    RatioYVtxPosition.at(0) -> Draw("SAME");
    TLine* Line1YVtxPosition = new TLine(-233/2,1,233/2,1);
    Line1YVtxPosition -> SetLineStyle(7);
    Line1YVtxPosition -> Draw("SAME");
    gPad->RedrawAxis();
    C7->SaveAs("ForwardFoldedYVtxPosition.pdf");

    TCanvas *C8 = new TCanvas("C8", "C8", 1000, 1000);
    TPad *pad1ZVtxPosition = new TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1ZVtxPosition->SetBottomMargin(0);
    pad1ZVtxPosition->Draw();
    pad1ZVtxPosition->cd();
    SelZVtxPosition.at(2) -> GetYaxis() -> SetRangeUser(-10,420);
    SelZVtxPosition.at(2) -> SetFillColorAlpha(46,0.3);
    SelZVtxPosition.at(2) -> Draw("E2 SAME");
    BgrZVtxPosition.at(0).back() -> SetLineColor(46);
    BgrZVtxPosition.at(0).back() -> Draw("HIST ][ SAME");
    StackBgrZVtxPosition -> Draw("HIST SAME");
    SelZVtxPosition.at(0) -> SetLineColor(1);
    SelZVtxPosition.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C8->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2ZVtxPosition = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2ZVtxPosition->SetTopMargin(0);
    pad2ZVtxPosition->SetBottomMargin(0.2);
    pad2ZVtxPosition->Draw();
    pad2ZVtxPosition->cd();
    RatioZVtxPosition.at(0) -> SetLineColor(13);
    RatioZVtxPosition.at(0) -> SetTitle("");
    RatioZVtxPosition.at(0) -> GetXaxis() -> SetLabelSize(0.09);
    RatioZVtxPosition.at(0) -> GetXaxis() -> SetTitleSize(0.1);
//     RatioZVtxPosition.at(0) -> GetYaxis() -> SetRangeUser(-0.1,1.7);
    RatioZVtxPosition.at(0) -> GetYaxis() -> SetTitle("OnBeam/(MC + BGR)");
    RatioZVtxPosition.at(0) -> GetYaxis() -> SetLabelSize(0.09);
    RatioZVtxPosition.at(0) -> GetYaxis() -> SetTitleSize(0.1);
    RatioZVtxPosition.at(0) -> GetYaxis() -> SetTitleOffset(0.5);
    RatioZVtxPosition.at(0) -> Draw("SAME");
    TLine* Line1ZVtxPosition = new TLine(0,1,1036.8,1);
    Line1ZVtxPosition -> SetLineStyle(7);
    Line1ZVtxPosition -> Draw("SAME");
    gPad->RedrawAxis();
    C8->SaveAs("ForwardFoldedZVtxPosition.pdf");

    // END DRAW FORWARD FOLDING --------------------------------------------------------------------------------------------------------------------------

    // BEGIN DRAW Models--------- ------------------------------------------------------------------------------------------------------------------------

    TCanvas *Canvas0 = new TCanvas("Range", "Range", 1400, 1000);
    Canvas0->cd();
    SelectionTrackRange.at(2) ->  GetYaxis() -> SetRangeUser(0,900);
    SelectionTrackRange.at(2) -> SetFillColorAlpha(46,0.3);
    SelectionTrackRange.at(2) -> Draw("E2 SAME");
    BgrTrackRange.at(0).back() -> SetLineColor(46);
    BgrTrackRange.at(0).back() -> Draw("HIST SAME");
    SelectionTrackRange.at(3) -> SetFillColorAlpha(28,0.3);
    SelectionTrackRange.at(3) -> SetMarkerColor(28);
    SelectionTrackRange.at(3) -> Draw("E2 SAME");
    BgrTrackRange.at(1).back() -> SetLineColor(28);
    BgrTrackRange.at(1).back() -> Draw("HIST SAME");
    SelectionTrackRange.at(4) -> SetFillColorAlpha(30,0.3);
    SelectionTrackRange.at(4) -> Draw("E2 SAME");
    BgrTrackRange.at(2).back() -> SetLineColor(30);
    BgrTrackRange.at(2).back() -> Draw("HIST SAME");
    SelectionTrackRange.at(5) -> SetFillColorAlpha(38,0.3);
    SelectionTrackRange.at(5) -> Draw("E2 SAME");
    BgrTrackRange.at(3).back() -> SetLineColor(38);
    BgrTrackRange.at(3).back() -> Draw("HIST SAME");
    SelectionTrackRange.at(0)->Draw("SAME");
    ModelLegend -> Draw();
    gPad->RedrawAxis();
    Canvas0->SaveAs("ModelComparisonTrackRange.pdf");

    ModelLegend -> SetX1NDC(0.15);
    ModelLegend -> SetX2NDC(0.52);

    TCanvas *Canvas1 = new TCanvas("CosTheta", "CosTheta", 1400, 1000);
    Canvas1->cd();
    SelectionCosTheta.at(2) ->  GetYaxis() -> SetRangeUser(0,1400);
    SelectionCosTheta.at(2) -> SetFillColorAlpha(46,0.3);
    SelectionCosTheta.at(2) -> Draw("E2 SAME");
    BgrCosTheta.at(0).back() -> SetLineColor(46);
    BgrCosTheta.at(0).back() -> Draw("HIST SAME");
    SelectionCosTheta.at(3) -> SetFillColorAlpha(28,0.3);
    SelectionCosTheta.at(3) -> SetMarkerColor(28);
    SelectionCosTheta.at(3) -> Draw("E2 SAME");
    BgrCosTheta.at(1).back() -> SetLineColor(28);
    BgrCosTheta.at(1).back() -> Draw("HIST SAME");
    SelectionCosTheta.at(4) -> SetFillColorAlpha(30,0.3);
    SelectionCosTheta.at(4) -> Draw("E2 SAME");
    BgrCosTheta.at(2).back() -> SetLineColor(30);
    BgrCosTheta.at(2).back() -> Draw("HIST SAME");
    SelectionCosTheta.at(5) -> SetFillColorAlpha(38,0.3);
    SelectionCosTheta.at(5) -> Draw("E2 SAME");
    BgrCosTheta.at(3).back() -> SetLineColor(38);
    BgrCosTheta.at(3).back() -> Draw("HIST SAME");
    SelectionCosTheta.at(0)->Draw("SAME");
    ModelLegend -> Draw();
    gPad->RedrawAxis();
    Canvas1->SaveAs("ModelComparisonCosTheta.pdf");

    ModelLegend -> SetX1NDC(0.48);
    ModelLegend -> SetX2NDC(0.85);

    TCanvas *Canvas2 = new TCanvas("Theta", "Theta", 1400, 1000);
    Canvas2->cd();
    SelectionTheta.at(2) ->  GetYaxis() -> SetRangeUser(0,700);
    SelectionTheta.at(2) -> SetFillColorAlpha(46,0.3);
    SelectionTheta.at(2) -> Draw("E2 SAME");
    BgrTheta.at(0).back() -> SetLineColor(46);
    BgrTheta.at(0).back() -> Draw("HIST SAME");
    SelectionTheta.at(3) -> SetFillColorAlpha(28,0.3);
    SelectionTheta.at(3) -> SetMarkerColor(28);
    SelectionTheta.at(3) -> Draw("E2 SAME");
    BgrTheta.at(1).back() -> SetLineColor(28);
    BgrTheta.at(1).back() -> Draw("HIST SAME");
    SelectionTheta.at(4) -> SetFillColorAlpha(30,0.3);
    SelectionTheta.at(4) -> Draw("E2 SAME");
    BgrTheta.at(2).back() -> SetLineColor(30);
    BgrTheta.at(2).back() -> Draw("HIST SAME");
    SelectionTheta.at(5) -> SetFillColorAlpha(38,0.3);
    SelectionTheta.at(5) -> Draw("E2 SAME");
    BgrTheta.at(3).back() -> SetLineColor(38);
    BgrTheta.at(3).back() -> Draw("HIST SAME");
    SelectionTheta.at(0)->Draw("SAME");
    ModelLegend -> Draw();
    gPad->RedrawAxis();
    Canvas2->SaveAs("ModelComparisonTheta.pdf");

    ModelLegend -> SetX1NDC(0.38);
    ModelLegend -> SetX2NDC(0.75);

    TCanvas *Canvas3 = new TCanvas("Phi", "Phi", 1400, 1000);
    Canvas3->cd();
    SelectionPhi.at(2) ->  GetYaxis() -> SetRangeUser(0,400);
    SelectionPhi.at(2) -> SetFillColorAlpha(46,0.3);
    SelectionPhi.at(2) -> Draw("E2 SAME");
    BgrPhi.at(0).back() -> SetLineColor(46);
    BgrPhi.at(0).back() -> Draw("HIST SAME");
    SelectionPhi.at(3) -> SetFillColorAlpha(28,0.3);
    SelectionPhi.at(3) -> SetMarkerColor(28);
    SelectionPhi.at(3) -> Draw("E2 SAME");
    BgrPhi.at(1).back() -> SetLineColor(28);
    BgrPhi.at(1).back() -> Draw("HIST SAME");
    SelectionPhi.at(4) -> SetFillColorAlpha(30,0.3);
    SelectionPhi.at(4) -> Draw("E2 SAME");
    BgrPhi.at(2).back() -> SetLineColor(30);
    BgrPhi.at(2).back() -> Draw("HIST SAME");
    SelectionPhi.at(5) -> SetFillColorAlpha(38,0.3);
    SelectionPhi.at(5) -> Draw("E2 SAME");
    BgrPhi.at(3).back() -> SetLineColor(38);
    BgrPhi.at(3).back() -> Draw("HIST SAME");
    SelectionPhi.at(0)->Draw("SAME");
    ModelLegend -> Draw();
    gPad->RedrawAxis();
    Canvas3->SaveAs("ModelComparisonPhi.pdf");

    ModelLegend -> SetX1NDC(0.48);
    ModelLegend -> SetX2NDC(0.85);

    TCanvas *Canvas4 = new TCanvas("Momentum", "Momentum", 1400, 1000);
    Canvas4->cd();
    SelectionMomentum.at(2) ->  GetYaxis() -> SetRangeUser(0,1400);
    SelectionMomentum.at(2) -> SetFillColorAlpha(46,0.3);
    SelectionMomentum.at(2) -> Draw("E2 SAME");
    BgrMomentum.at(0).back() -> SetLineColor(46);
    BgrMomentum.at(0).back() -> Draw("HIST SAME");
    SelectionMomentum.at(3) -> SetFillColorAlpha(28,0.3);
    SelectionMomentum.at(3) -> SetMarkerColor(28);
    SelectionMomentum.at(3) -> Draw("E2 SAME");
    BgrMomentum.at(1).back() -> SetLineColor(28);
    BgrMomentum.at(1).back() -> Draw("HIST SAME");
    SelectionMomentum.at(4) -> SetFillColorAlpha(30,0.3);
    SelectionMomentum.at(4) -> Draw("E2 SAME");
    BgrMomentum.at(2).back() -> SetLineColor(30);
    BgrMomentum.at(2).back() -> Draw("HIST SAME");
    SelectionMomentum.at(5) -> SetFillColorAlpha(38,0.3);
    SelectionMomentum.at(5) -> Draw("E2 SAME");
    BgrMomentum.at(3).back() -> SetLineColor(38);
    BgrMomentum.at(3).back() -> Draw("HIST SAME");
    SelectionMomentum.at(0)->Draw("SAME");
    ModelLegend -> Draw();
    gPad->RedrawAxis();
    Canvas4->SaveAs("ModelComparisonMomentum.pdf");

    TCanvas *Canvas5 = new TCanvas("Length", "Length", 1400, 1000);
    Canvas5->cd();
    SelectionTrackLength.at(2) ->  GetYaxis() -> SetRangeUser(0,1000);
    SelectionTrackLength.at(2) -> SetFillColorAlpha(46,0.3);
    SelectionTrackLength.at(2) -> Draw("E2 SAME");
    BgrTrackLength.at(0).back() -> SetLineColor(46);
    BgrTrackLength.at(0).back() -> Draw("HIST SAME");
    SelectionTrackLength.at(3) -> SetFillColorAlpha(28,0.3);
    SelectionTrackLength.at(3) -> SetMarkerColor(28);
    SelectionTrackLength.at(3) -> Draw("E2 SAME");
    BgrTrackLength.at(1).back() -> SetLineColor(28);
    BgrTrackLength.at(1).back() -> Draw("HIST SAME");
    SelectionTrackLength.at(4) -> SetFillColorAlpha(30,0.3);
    SelectionTrackLength.at(4) -> Draw("E2 SAME");
    BgrTrackLength.at(2).back() -> SetLineColor(30);
    BgrTrackLength.at(2).back() -> Draw("HIST SAME");
    SelectionTrackLength.at(5) -> SetFillColorAlpha(38,0.3);
    SelectionTrackLength.at(5) -> Draw("E2 SAME");
    BgrTrackLength.at(3).back() -> SetLineColor(38);
    BgrTrackLength.at(3).back() -> Draw("HIST SAME");
    SelectionTrackLength.at(0)->Draw("SAME");
    ModelLegend -> Draw();
    gPad->RedrawAxis();
    Canvas5->SaveAs("ModelComparisonTrackLength.pdf");

    ModelLegend -> SetX1NDC(0.15);
    ModelLegend -> SetX2NDC(0.52);

    TCanvas *Canvas6 = new TCanvas("XVtx", "XVtx", 1400, 1000);
    Canvas6->cd();
    SelXVtxPosition.at(2) ->  GetYaxis() -> SetRangeUser(0,400);
    SelXVtxPosition.at(2) -> SetFillColorAlpha(46,0.3);
    SelXVtxPosition.at(2) -> Draw("E2 SAME");
    BgrXVtxPosition.at(0).back() -> SetLineColor(46);
    BgrXVtxPosition.at(0).back() -> Draw("HIST SAME");
    SelXVtxPosition.at(3) -> SetFillColorAlpha(28,0.3);
    SelXVtxPosition.at(3) -> SetMarkerColor(28);
    SelXVtxPosition.at(3) -> Draw("E2 SAME");
    BgrXVtxPosition.at(1).back() -> SetLineColor(28);
    BgrXVtxPosition.at(1).back() -> Draw("HIST SAME");
    SelXVtxPosition.at(4) -> SetFillColorAlpha(30,0.3);
    SelXVtxPosition.at(4) -> Draw("E2 SAME");
    BgrXVtxPosition.at(2).back() -> SetLineColor(30);
    BgrXVtxPosition.at(2).back() -> Draw("HIST SAME");
    SelXVtxPosition.at(5) -> SetFillColorAlpha(38,0.3);
    SelXVtxPosition.at(5) -> Draw("E2 SAME");
    BgrXVtxPosition.at(3).back() -> SetLineColor(38);
    BgrXVtxPosition.at(3).back() -> Draw("HIST SAME");
    SelXVtxPosition.at(0)->Draw("SAME");
    ModelLegend -> Draw();
    gPad->RedrawAxis();
    Canvas6->SaveAs("ModelComparisonXVtx.pdf");

    TCanvas *Canvas7 = new TCanvas("YVtx", "YVtx", 1400, 1000);
    Canvas7->cd();
    SelYVtxPosition.at(2) ->  GetYaxis() -> SetRangeUser(0,400);
    SelYVtxPosition.at(2) -> SetFillColorAlpha(46,0.3);
    SelYVtxPosition.at(2) -> Draw("E2 SAME");
    BgrYVtxPosition.at(0).back() -> SetLineColor(46);
    BgrYVtxPosition.at(0).back() -> Draw("HIST SAME");
    SelYVtxPosition.at(3) -> SetFillColorAlpha(28,0.3);
    SelYVtxPosition.at(3) -> SetMarkerColor(28);
    SelYVtxPosition.at(3) -> Draw("E2 SAME");
    BgrYVtxPosition.at(1).back() -> SetLineColor(28);
    BgrYVtxPosition.at(1).back() -> Draw("HIST SAME");
    SelYVtxPosition.at(4) -> SetFillColorAlpha(30,0.3);
    SelYVtxPosition.at(4) -> Draw("E2 SAME");
    BgrYVtxPosition.at(2).back() -> SetLineColor(30);
    BgrYVtxPosition.at(2).back() -> Draw("HIST SAME");
    SelYVtxPosition.at(5) -> SetFillColorAlpha(38,0.3);
    SelYVtxPosition.at(5) -> Draw("E2 SAME");
    BgrYVtxPosition.at(3).back() -> SetLineColor(38);
    BgrYVtxPosition.at(3).back() -> Draw("HIST SAME");
    SelYVtxPosition.at(0)->Draw("SAME");
    ModelLegend -> Draw();
    gPad->RedrawAxis();
    Canvas7->SaveAs("ModelComparisonYVtx.pdf");

    TCanvas *Canvas8 = new TCanvas("ZVtx", "ZVtx", 1400, 1000);
    Canvas8->cd();
    SelZVtxPosition.at(2) ->  GetYaxis() -> SetRangeUser(0,400);
    SelZVtxPosition.at(2) -> SetFillColorAlpha(46,0.3);
    SelZVtxPosition.at(2) -> Draw("E2 SAME");
    BgrZVtxPosition.at(0).back() -> SetLineColor(46);
    BgrZVtxPosition.at(0).back() -> Draw("HIST SAME");
    SelZVtxPosition.at(3) -> SetFillColorAlpha(28,0.3);
    SelZVtxPosition.at(3) -> SetMarkerColor(28);
    SelZVtxPosition.at(3) -> Draw("E2 SAME");
    BgrZVtxPosition.at(1).back() -> SetLineColor(28);
    BgrZVtxPosition.at(1).back() -> Draw("HIST SAME");
    SelZVtxPosition.at(4) -> SetFillColorAlpha(30,0.3);
    SelZVtxPosition.at(4) -> Draw("E2 SAME");
    BgrZVtxPosition.at(2).back() -> SetLineColor(30);
    BgrZVtxPosition.at(2).back() -> Draw("HIST SAME");
    SelZVtxPosition.at(5) -> SetFillColorAlpha(38,0.3);
    SelZVtxPosition.at(5) -> Draw("E2 SAME");
    BgrZVtxPosition.at(3).back() -> SetLineColor(38);
    BgrZVtxPosition.at(3).back() -> Draw("HIST SAME");
    SelZVtxPosition.at(0)->Draw("SAME");
    ModelLegend -> Draw();
    gPad->RedrawAxis();
    Canvas8->SaveAs("ModelComparisonZVtx.pdf");

//     TCanvas *C1 = new TCanvas("C1", "C1", 1400, 1000);
//     MomentumBeamSys.at(1).at(0) -> Draw();
}


float CalcRange(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2)
{
    return sqrt(pow(x_1-x_2, 2) + pow(y_1-y_2, 2) + pow(z_1-z_2, 2));
}

bool inFV(double x, double y, double z)
{
    if(x < (FVx - borderx) && x > borderx && y < (FVy/2. - bordery) && y > (-FVy/2. + bordery) && z < (FVz - borderz) && z > borderz) return true;
    else return false;
}

bool inTPC(double x, double y, double z)
{
    if(x < FVx && x > 0 && y < FVy/2. && y > -FVy/2. && z < FVz && z > 0) return true;
    else return false;
}

void AddHistograms(std::vector<TH1F*>& HistVector, unsigned int First, unsigned int Last, float Weight, bool EraseLast)
{
    // Check if there is something to be added
    if (HistVector.size() > Last)
    {
        // Add histograms
        HistVector.at(First) -> Add(HistVector.at(Last), Weight);

        // Erase last histogram if flag is set
        if(EraseLast)
        {
            HistVector.at(Last)->Delete();
            HistVector.erase(HistVector.begin() + Last);
        }
    }
    else // if nothing can be added
    {
        std::cout << "Histograms not added!" << std::endl;
    }
}

TH1F* AddToNewHist(std::vector<TH1F*>& HistVector, unsigned int First, unsigned int Last, float Weight)
{
    TH1F *OutputHist = (TH1F*) HistVector.at(First)->Clone();

    // Check if there is something to be added
    if (HistVector.size() > Last)
    {
        // Add histograms
        OutputHist -> Add(HistVector.at(Last), Weight);
    }
    else // if nothing can be added
    {
        std::cout << "Histograms not added!" << std::endl;
    }
    return OutputHist;
}

void SubtractBgr(std::vector<TH1F*>& HistVector, std::vector<std::vector<TH1F*>>& BgrVector, unsigned int First, unsigned int Last, float Weight)
{
    // Check if there is something to be added
    if (HistVector.size() > First && BgrVector.size() > Last)
    {
        // Add histograms
        HistVector.at(First) -> Add(BgrVector.at(Last).at(0), -Weight);
    }
    else // if nothing can be added
    {
        std::cout << "Histograms not added!" << std::endl;
    }
}

void NormMatrixByColumn(TH2F* UMatrix)
{
    // loop over xbins of the smearing matrices
    for(unsigned int xbin = 1; xbin <= UMatrix->GetNbinsX(); xbin++)
    {
        float NormFact = 0;

        // loop over ybins (column)
        for(unsigned int ybin = 1; ybin <= UMatrix->GetNbinsY(); ybin++)
        {
            // Add column entry to normalization factor
            NormFact += UMatrix->GetBinContent(xbin,ybin);
        } // ybin loop

        // loop over xbins (column)
        for(unsigned int ybin = 1; ybin <= UMatrix->GetNbinsY(); ybin++)
        {
            // Normalize entire row of the matrix
            if(NormFact) UMatrix->SetBinContent(xbin,ybin,UMatrix->GetBinContent(xbin,ybin)/NormFact) ;
        }// ybin loop
    }// xbin loop
}

void SelectionUnsmearing(TH2F*& UMatrix, TH1F*& SVector)
{
    TH1F* CloneVector = (TH1F*)SVector->Clone();

    // loop over xbins (row)
    for(unsigned int xbin = 1; xbin <= UMatrix->GetNbinsX(); xbin++)
    {
        // Unsmeared bin content
        float UnsmearedContent = 0;

        // loop over ybins of the unsmearing matrices
        for(unsigned int ybin = 1; ybin <= UMatrix->GetNbinsY(); ybin++)
        {
            // Sum up the vertical contribution of the unsmeared vector
            UnsmearedContent += CloneVector->GetBinContent(ybin)*UMatrix->GetBinContent(xbin,ybin);
        }
        // Fill unsmeared content into vector again
        SVector->SetBinContent(xbin,UnsmearedContent);
    }
}

void MomentumSplinePreparation()
{
    float RangeGramPerCM[29] = {9.833E-1, 1.786E0, 3.321E0, 6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1, 1.063E2, 1.725E2,
                                2.385E2, 4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3, 2.297E3, 4.359E3, 5.354E3, 7.298E3,
                                1.013E4, 1.469E4, 1.910E4, 3.558E4, 4.326E4, 5.768E4, 7.734E4, 1.060E5, 1.307E5
                               };

    float KEMeV[29] = {10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000, 1400, 2000, 3000, 4000,
                       8000, 10000, 14000, 20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000
                      };

    // convert to cm
    for(auto & RangePoint : RangeGramPerCM)
    {
        RangePoint /= Density;
    }

    TGraph* KEvsR = new TGraph(29, RangeGramPerCM, KEMeV);
//     KEvsR -> Draw();

    KEvsRSpline = new TSpline3("KEvsRS",KEvsR);

    delete KEvsR;
}

float GetMomentum(float TrackLength)
{
    float MuonMass = 105.7; //MeV

    // Change Track length to kinetic energy
    TrackLength = KEvsRSpline->Eval(TrackLength);

    // Convert kinetic energy to momentum
    TrackLength = std::sqrt( std::pow(TrackLength,2) + 2*TrackLength*MuonMass );

    // Convert MeV to GeV
    TrackLength /= 1000;

    return TrackLength;
}

void CalcSigEfficiency (std::vector<TH1F*>& HistVector)
{
    for(unsigned bin_no = 1; bin_no <= HistVector.at(3)->GetNbinsX(); bin_no++)
    {
        float SignalBinContent = HistVector.at(2)->GetBinContent(bin_no);

        // This calculates the per bin efficiency of signal compared to cosmic contamination
        float Efficiency = SignalBinContent / ( SignalBinContent + HistVector.at(3)->GetBinContent(bin_no) );
        HistVector.at(3)->SetBinContent(bin_no,Efficiency);
    }
}

