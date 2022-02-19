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
void DrawCCInclusive()
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

    GenLabel.push_back(std::make_pair("MA Adjusted Selection",38));
    ScalingFactors.push_back(DataPOT/MCMaPOT);

    GenLabel.push_back(std::make_pair("TEM Selection",30));
    ScalingFactors.push_back(DataPOT/MCTEMPOT);

    GenLabel.push_back(std::make_pair("MEC Selection",42));
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
            BgrTrackRange.at(file_no).back()->SetFillColor(Label.second);
            BgrTrackRange.at(file_no).back()->SetLineColor(Label.second);
            BgrTrackRange.at(file_no).back()->SetMarkerColor(Label.second);

            BgrCosTheta.at(file_no).back()->SetFillColor(Label.second);
            BgrCosTheta.at(file_no).back()->SetLineColor(Label.second);
            BgrCosTheta.at(file_no).back()->SetMarkerColor(Label.second);

            BgrTheta.at(file_no).back()->SetFillColor(Label.second);
            BgrTheta.at(file_no).back()->SetLineColor(Label.second);
            BgrTheta.at(file_no).back()->SetMarkerColor(Label.second);

            BgrPhi.at(file_no).back()->SetFillColor(Label.second);
            BgrPhi.at(file_no).back()->SetLineColor(Label.second);
            BgrPhi.at(file_no).back()->SetMarkerColor(Label.second);

            BgrMomentum.at(file_no).back()->SetFillColor(Label.second);
            BgrMomentum.at(file_no).back()->SetLineColor(Label.second);
            BgrMomentum.at(file_no).back()->SetMarkerColor(Label.second);

            BgrTrackLength.at(file_no).back()->SetFillColor(Label.second);
            BgrTrackLength.at(file_no).back()->SetLineColor(Label.second);
            BgrTrackLength.at(file_no).back()->SetMarkerColor(Label.second);

            BgrXVtxPosition.at(file_no).back()->SetFillColor(Label.second);
            BgrXVtxPosition.at(file_no).back()->SetLineColor(Label.second);
            BgrXVtxPosition.at(file_no).back()->SetMarkerColor(Label.second);

            BgrYVtxPosition.at(file_no).back()->SetFillColor(Label.second);
            BgrYVtxPosition.at(file_no).back()->SetLineColor(Label.second);
            BgrYVtxPosition.at(file_no).back()->SetMarkerColor(Label.second);

            BgrZVtxPosition.at(file_no).back()->SetFillColor(Label.second);
            BgrZVtxPosition.at(file_no).back()->SetLineColor(Label.second);
            BgrZVtxPosition.at(file_no).back()->SetMarkerColor(Label.second);
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

    // Generate New set of histograms containing the distributions with systematics + statistics
    std::vector<TH1F*> SystTrackRange;
    std::vector<TH1F*> SystCosTheta;
    std::vector<TH1F*> SystTheta;
    std::vector<TH1F*> SystPhi;
    std::vector<TH1F*> SystMomentum;
    std::vector<TH1F*> SystTrackLength;
    std::vector<TH1F*> SystXVtxPosition;
    std::vector<TH1F*> SystYVtxPosition;
    std::vector<TH1F*> SystZVtxPosition;

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
        TrackRangeBeamSys.at(file_no).push_front( (TH1F*)CosmicTrackRange.back()->Clone() );
        CosThetaBeamSys.at(file_no).push_front( (TH1F*)CosmicCosTheta.back()->Clone() );
        ThetaBeamSys.at(file_no).push_front( (TH1F*)CosmicTheta.back()->Clone() );
        PhiBeamSys.at(file_no).push_front( (TH1F*)CosmicPhi.back()->Clone() );
        MomentumBeamSys.at(file_no).push_front( (TH1F*)CosmicMomentum.back()->Clone() );
        TrackLengthBeamSys.at(file_no).push_front( (TH1F*)CosmicTrackLength.back()->Clone() );
        XVtxPositionBeamSys.at(file_no).push_front( (TH1F*)CosmicXVtxPosition.back()->Clone() );
        YVtxPositionBeamSys.at(file_no).push_front( (TH1F*)CosmicYVtxPosition.back()->Clone() );
        ZVtxPositionBeamSys.at(file_no).push_front( (TH1F*)CosmicZVtxPosition.back()->Clone() );

        // Multiply squared background to get variance
        TrackRangeBeamSys.at(file_no).front() -> Multiply(BgrTrackRange.at(file_no).at(1));
        TrackRangeBeamSys.at(file_no).front() -> Multiply(BgrTrackRange.at(file_no).at(1));
        CosThetaBeamSys.at(file_no).front() -> Multiply(BgrCosTheta.at(file_no).at(1));
        CosThetaBeamSys.at(file_no).front() -> Multiply(BgrCosTheta.at(file_no).at(1));
        ThetaBeamSys.at(file_no).front() -> Multiply(BgrTheta.at(file_no).at(1));
        ThetaBeamSys.at(file_no).front() -> Multiply(BgrTheta.at(file_no).at(1));
        PhiBeamSys.at(file_no).front() -> Multiply(BgrPhi.at(file_no).at(1));
        PhiBeamSys.at(file_no).front() -> Multiply(BgrPhi.at(file_no).at(1));
        MomentumBeamSys.at(file_no).front() -> Multiply(BgrMomentum.at(file_no).at(1));
        MomentumBeamSys.at(file_no).front() -> Multiply(BgrMomentum.at(file_no).at(1));
        TrackLengthBeamSys.at(file_no).front() -> Multiply(BgrTrackLength.at(file_no).at(1));
        TrackLengthBeamSys.at(file_no).front() -> Multiply(BgrTrackLength.at(file_no).at(1));
        XVtxPositionBeamSys.at(file_no).front() -> Multiply(BgrXVtxPosition.at(file_no).at(1));
        XVtxPositionBeamSys.at(file_no).front() -> Multiply(BgrXVtxPosition.at(file_no).at(1));
        YVtxPositionBeamSys.at(file_no).front() -> Multiply(BgrYVtxPosition.at(file_no).at(1));
        YVtxPositionBeamSys.at(file_no).front() -> Multiply(BgrYVtxPosition.at(file_no).at(1));
        ZVtxPositionBeamSys.at(file_no).front() -> Multiply(BgrZVtxPosition.at(file_no).at(1));
        ZVtxPositionBeamSys.at(file_no).front() -> Multiply(BgrZVtxPosition.at(file_no).at(1));

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
        TrackRangeBeamSys.at(file_no).push_front(new TH1F(("Total Systematics Range "+std::to_string(file_no)).c_str(),"Range",NumberOfBins,0,700));
        CosThetaBeamSys.at(file_no).push_front(new TH1F(("Total Systematics cos#theta "+std::to_string(file_no)).c_str(),"cos#theta",NumberOfBins,-1,1));
        ThetaBeamSys.at(file_no).push_front(new TH1F(("Total Systematics #theta "+std::to_string(file_no)).c_str(),"#theta",NumberOfBins,0,180));
        PhiBeamSys.at(file_no).push_front(new TH1F(("Total Systematics #phi "+std::to_string(file_no)).c_str(),"#phi",NumberOfBins,-180,180));
        MomentumBeamSys.at(file_no).push_front(new TH1F(("Total Systematics Momentum "+std::to_string(file_no)).c_str(),"Momentum",NumberOfBins,0,3));
        TrackLengthBeamSys.at(file_no).push_front(new TH1F(("Total Systematics Length "+std::to_string(file_no)).c_str(),"Lenght",NumberOfBins,0,800));
        XVtxPositionBeamSys.at(file_no).push_front(new TH1F(("Total Systematics XVtx "+std::to_string(file_no)).c_str(),"XVtx",NumberOfBins,0,256.35));
        YVtxPositionBeamSys.at(file_no).push_front(new TH1F(("Total Systematics YVtx "+std::to_string(file_no)).c_str(),"YVtx",NumberOfBins,-233*0.5,233*0.5));
        ZVtxPositionBeamSys.at(file_no).push_front(new TH1F(("Total Systematics ZVtx "+std::to_string(file_no)).c_str(),"ZVtx",NumberOfBins,0,1036.8));

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

        // Calculate standard deviation
        for(unsigned int bin_no = 1; bin_no <= NumberOfBins; bin_no++)
        {
            TrackRangeBeamSys.at(file_no).front()->SetBinContent(bin_no,std::sqrt(TrackRangeBeamSys.at(file_no).front()->GetBinContent(bin_no)));
            CosThetaBeamSys.at(file_no).front()->SetBinContent(bin_no,std::sqrt(CosThetaBeamSys.at(file_no).front()->GetBinContent(bin_no)));
            ThetaBeamSys.at(file_no).front()->SetBinContent(bin_no,std::sqrt(ThetaBeamSys.at(file_no).front()->GetBinContent(bin_no)));
            PhiBeamSys.at(file_no).front()->SetBinContent(bin_no,std::sqrt(PhiBeamSys.at(file_no).front()->GetBinContent(bin_no)));
            MomentumBeamSys.at(file_no).front()->SetBinContent(bin_no,std::sqrt(MomentumBeamSys.at(file_no).front()->GetBinContent(bin_no)));
            TrackLengthBeamSys.at(file_no).front()->SetBinContent(bin_no,std::sqrt(TrackLengthBeamSys.at(file_no).front()->GetBinContent(bin_no)));
            XVtxPositionBeamSys.at(file_no).front()->SetBinContent(bin_no,std::sqrt(XVtxPositionBeamSys.at(file_no).front()->GetBinContent(bin_no)));
            YVtxPositionBeamSys.at(file_no).front()->SetBinContent(bin_no,std::sqrt(YVtxPositionBeamSys.at(file_no).front()->GetBinContent(bin_no)));
            ZVtxPositionBeamSys.at(file_no).front()->SetBinContent(bin_no,std::sqrt(ZVtxPositionBeamSys.at(file_no).front()->GetBinContent(bin_no)));
        }

        // Fill data point with statistic uncertainties first
        SystTrackRange.push_back( (TH1F*)SelectionTrackRange.at(file_no+2)->Clone() );
        SystCosTheta.push_back( (TH1F*)SelectionCosTheta.at(file_no+2)->Clone() );
        SystTheta.push_back( (TH1F*)SelectionTheta.at(file_no+2)->Clone() );
        SystPhi.push_back( (TH1F*)SelectionPhi.at(file_no+2)->Clone() );
        SystMomentum.push_back( (TH1F*)SelectionMomentum.at(file_no+2)->Clone() );
        SystTrackLength.push_back( (TH1F*)SelectionTrackLength.at(file_no+2)->Clone() );
        SystXVtxPosition.push_back( (TH1F*)SelXVtxPosition.at(file_no+2)->Clone() );
        SystYVtxPosition.push_back( (TH1F*)SelYVtxPosition.at(file_no+2)->Clone() );
        SystZVtxPosition.push_back( (TH1F*)SelZVtxPosition.at(file_no+2)->Clone() );

        SystTrackRange.back() -> Add(SelectionTrackRange.at(1));
        SystCosTheta.back() -> Add(SelectionCosTheta.at(1));
        SystTheta.back() -> Add(SelectionTheta.at(1));
        SystPhi.back() -> Add(SelectionPhi.at(1));
        SystMomentum.back() -> Add(SelectionMomentum.at(1));
        SystTrackLength.back() -> Add(SelectionTrackLength.at(1));
        SystXVtxPosition.back() -> Add(SelXVtxPosition.at(1));
        SystYVtxPosition.back() -> Add(SelYVtxPosition.at(1));
        SystZVtxPosition.back() -> Add(SelZVtxPosition.at(1));

        // Generate New set of histograms containing the distributions with systematics + statistics
        for(unsigned int bin_no = 1; bin_no <= NumberOfBins; bin_no++)
        {
            SystTrackRange.back()->SetBinError(bin_no,SystTrackRange.back()->GetBinError(bin_no)+TrackRangeBeamSys.at(file_no).front()->GetBinContent(bin_no));
            SystCosTheta.back()->SetBinError(bin_no,SystCosTheta.back()->GetBinError(bin_no)+CosThetaBeamSys.at(file_no).front()->GetBinContent(bin_no));
            SystTheta.back()->SetBinError(bin_no,SystTheta.back()->GetBinError(bin_no)+ThetaBeamSys.at(file_no).front()->GetBinContent(bin_no));
            SystPhi.back()->SetBinError(bin_no,SystPhi.back()->GetBinError(bin_no)+PhiBeamSys.at(file_no).front()->GetBinContent(bin_no));
            SystMomentum.back()->SetBinError(bin_no,SystMomentum.back()->GetBinError(bin_no)+MomentumBeamSys.at(file_no).front()->GetBinContent(bin_no));
            SystTrackLength.back()->SetBinError(bin_no,SystTrackLength.back()->GetBinError(bin_no)+TrackLengthBeamSys.at(file_no).front()->GetBinContent(bin_no));
            SystXVtxPosition.back()->SetBinError(bin_no,SystXVtxPosition.back()->GetBinError(bin_no)+XVtxPositionBeamSys.at(file_no).front()->GetBinContent(bin_no));
            SystYVtxPosition.back()->SetBinError(bin_no,SystYVtxPosition.back()->GetBinError(bin_no)+YVtxPositionBeamSys.at(file_no).front()->GetBinContent(bin_no));
            SystZVtxPosition.back()->SetBinError(bin_no,SystZVtxPosition.back()->GetBinError(bin_no)+ZVtxPositionBeamSys.at(file_no).front()->GetBinContent(bin_no));
        }
    }

    // Add Off-Beam sample to MC prediction standard (for all MC predictions, put this into loop with file_no+2 index)
    for(unsigned int file_no = 0; file_no < 4; file_no++)
    {
        SelectionTrackRange.at(file_no+2)-> Add(SelectionTrackRange.at(1));
        SelectionCosTheta.at(file_no+2)-> Add(SelectionCosTheta.at(1));
        SelectionTheta.at(file_no+2)-> Add(SelectionTheta.at(1));
        SelectionPhi.at(file_no+2)-> Add(SelectionPhi.at(1));
        SelectionMomentum.at(file_no+2)-> Add(SelectionMomentum.at(1));
        SelectionTrackLength.at(file_no+2)-> Add(SelectionTrackLength.at(1));
        SelXVtxPosition.at(file_no+2)-> Add(SelXVtxPosition.at(1));
        SelYVtxPosition.at(file_no+2)-> Add(SelYVtxPosition.at(1));
        SelZVtxPosition.at(file_no+2)-> Add(SelZVtxPosition.at(1));
    }
    // END HISTOGRAM CALCULATIONS ------------------------------------------------------------------------------------------------------------------------

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
    
    // BEGIN LEGEND --------------------------------------------------------------------------------------------------------------------------------------
    TLegend* Legend = new TLegend(0.48,0.41,0.85,0.85);
    Legend->SetLineStyle ( 0 );
    Legend->SetLineColorAlpha ( 0,0 );
    Legend->SetFillStyle ( 0 );
    Legend->SetMargin ( 0.2 );
//     Legend->SetTextFont ( 43 );
//     Legend->SetTextSize ( 35 );
    Legend->SetHeader("Normalised to 4.95 #times 10^{19} POT","C");
    
    TLegendEntry *header = (TLegendEntry*)Legend->GetListOfPrimitives()->First();
//     header->SetTextAlign(22);
//     header->SetTextColor(2);
    header->SetTextSize(0.05);
    
    
    std::vector<std::string> LegendLabel;
    LegendLabel.push_back("Off-Beam Cosmic (Data)");
    LegendLabel.push_back("Cosmic BGR (MC)");
    LegendLabel.push_back("Dirt BGR (MC)");
    LegendLabel.push_back("Out of FV BGR (MC)");
    LegendLabel.push_back("#bar{#nu}_{#mu} BGR (MC)");
    LegendLabel.push_back("#nu_{e} & #bar{#nu}_{e} BGR (MC)");
    LegendLabel.push_back("#nu_{x} NC BGR events (MC)");
    LegendLabel.push_back("Selected #nu_{#mu} CC + BGR (MC)");
    
    Legend->AddEntry( SelectionTrackRange.at(0), "On-Beam BNB (Data)","lep" );
    
//     Legend->AddEntry( SelectionTrackRange.at(2), "Selected #nu_{#mu} CC signal & bgr","f" );
    for(int bgr_no = LegendLabel.size()-1; bgr_no > 0 ; bgr_no--)
    {
        Legend->AddEntry( BgrTrackRange.at(0).at(bgr_no), (LegendLabel.at(bgr_no)).c_str(),"f" );
    }
    Legend->AddEntry( SelectionTrackRange.at(1), (LegendLabel.front()).c_str(),"f" );
    
    // END LEGEND ----------------------------------------------------------------------------------------------------------------------------------------
//     std::cout << CosmicMomentum.back() << std::endl;

//     TrackRangeBeamSys
//     CosThetaBeamSys
//     ThetaBeamSys
//     PhiBeamSys
//     MomentumBeamSys
//     TrackLengthBeamSys
//     XVtxPositionBeamSys
//     YVtxPositionBeamSys
//     ZVtxPositionBeamSys

//     0 SystLabel.push_back("all");
//     1 SystLabel.push_back("cosmic");
//     2 SystLabel.push_back("dirt");
//     3 SystLabel.push_back("outFV");
//     4 SystLabel.push_back("anti nu_mu");
//     5 SystLabel.push_back("n_e-like");
//     6 SystLabel.push_back("nu_NC");
//     7 SystLabel.push_back("PureSelected");

    TCanvas *C0 = new TCanvas("C0", "C0", 700, 500);
    SystTrackRange.at(0) -> SetFillColorAlpha(46,0.35);
    SystTrackRange.at(0) -> GetYaxis() -> SetRangeUser(0,900);
    SystTrackRange.at(0) -> Draw("E2 SAME");
    StackBgrTrackRange -> Draw("HIST SAME");
    SelectionTrackRange.at(2) -> Draw("E2 SAME");
    SelectionTrackRange.at(0) -> SetLineColor(1);
    SelectionTrackRange.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C0->SaveAs("ForwardFoldedTrackRange.pdf");
    
    Legend -> SetX1NDC(0.15);
    Legend -> SetX2NDC(0.52);

    TCanvas *C1 = new TCanvas("C1", "C1", 700, 500);
    SystCosTheta.at(0) -> SetFillColorAlpha(46,0.35);
    SystCosTheta.at(0) -> GetYaxis() -> SetRangeUser(0,1400);
    SystCosTheta.at(0) -> Draw("E2 SAME");
    StackBgrCosTheta -> Draw("HIST SAME");
    SelectionCosTheta.at(2) -> Draw("E2 SAME");
    SelectionCosTheta.at(0) -> SetLineColor(1);
    SelectionCosTheta.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C1->SaveAs("ForwardFoldedCosTheta.pdf");
    
    Legend -> SetX1NDC(0.48);
    Legend -> SetX2NDC(0.85);

    TCanvas *C2 = new TCanvas("C2", "C2", 700, 500);
    SystTheta.at(0) -> SetFillColorAlpha(46,0.35);
    SystTheta.at(0) -> GetYaxis() -> SetRangeUser(0,700);
    SystTheta.at(0) -> Draw("E2 SAME");
    StackBgrTheta -> Draw("HIST SAME");
    SelectionTheta.at(2) -> Draw("E2 SAME");
    SelectionTheta.at(0) -> SetLineColor(1);
    SelectionTheta.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C2->SaveAs("ForwardFoldedTheta.pdf");

//     Legend -> SetX1NDC(0.40);
//     Legend -> SetX2NDC(0.77);
    
    TCanvas *C3 = new TCanvas("C3", "C3", 700, 500);
    SystPhi.at(0) -> SetFillColorAlpha(46,0.35);
    SystPhi.at(0) -> SetLineColorAlpha(46,0.35);
    SystPhi.at(0) -> GetYaxis() -> SetRangeUser(0,400);
    SystPhi.at(0) -> Draw("E2 SAME");
    StackBgrPhi-> Draw("HIST SAME");
    SelectionPhi.at(2) -> Draw("E2 SAME");
    SelectionPhi.at(0) -> SetLineColor(1);
    SelectionPhi.at(0) -> Draw("SAME");
//     Legend -> Draw();
    gPad->RedrawAxis();
    C3->SaveAs("ForwardFoldedPhi.pdf");

    TCanvas *C4 = new TCanvas("C4", "C4", 700, 500);
    SystMomentum.at(0) -> SetFillColorAlpha(46,0.35);
    SystMomentum.at(0) -> GetYaxis() -> SetRangeUser(0,1400);
    SystMomentum.at(0) -> Draw("E2 SAME");
    StackBgrMomentum -> Draw("HIST SAME");
    SelectionMomentum.at(2) -> Draw("E2 SAME");
    SelectionMomentum.at(0) -> SetLineColor(1);
    SelectionMomentum.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C4->SaveAs("ForwardFoldedMomentum.pdf");

    TCanvas *C5 = new TCanvas("C5", "C5", 700, 500);
    SystTrackLength.at(0) -> SetFillColorAlpha(46,0.35);
    SystTrackLength.at(0) -> GetYaxis() -> SetRangeUser(0,1000);
    SystTrackLength.at(0) -> Draw("E2 SAME");
    StackBgrTrackLength -> Draw("HIST SAME");
    SelectionTrackLength.at(2) -> Draw("E2 SAME");
    SelectionTrackLength.at(0) -> SetLineColor(1);
    SelectionTrackLength.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C5->SaveAs("ForwardFoldedTrackLength.pdf");
    
    Legend -> SetNColumns(2);
    Legend -> SetX1NDC(0.13);
    Legend -> SetX2NDC(0.87);
    Legend -> SetY1NDC(0.60);
    Legend -> SetY2NDC(0.85);
    
    TCanvas *C6 = new TCanvas("C6", "C6", 700, 500);
    SystXVtxPosition.at(0) -> SetFillColorAlpha(46,0.35);
    SystXVtxPosition.at(0) -> GetYaxis() -> SetRangeUser(0,400);
    SystXVtxPosition.at(0) -> Draw("E2 SAME");
    StackBgrXVtxPosition -> Draw("HIST SAME");
    SelXVtxPosition.at(2) -> Draw("E2 SAME");
    SelXVtxPosition.at(0) -> SetLineColor(1);
    SelXVtxPosition.at(0) -> Draw("SAME");
    Legend -> Draw();
    gPad->RedrawAxis();
    C6->SaveAs("ForwardFoldedXVtxPosition.pdf");

    TCanvas *C7 = new TCanvas("C7", "C7", 700, 500);
    SystYVtxPosition.at(0) -> SetFillColorAlpha(46,0.35);
    SystYVtxPosition.at(0) -> GetYaxis() -> SetRangeUser(0,400);
    SystYVtxPosition.at(0) -> Draw("E2 SAME");
    StackBgrYVtxPosition -> Draw("HIST SAME");
    SelYVtxPosition.at(2) -> Draw("E2 SAME");
    SelYVtxPosition.at(0) -> SetLineColor(1);
    SelYVtxPosition.at(0) -> Draw("SAME");
    gPad->RedrawAxis();
    C7->SaveAs("ForwardFoldedYVtxPosition.pdf");

    TCanvas *C8 = new TCanvas("C8", "C8", 700, 500);
    SystZVtxPosition.at(0) -> SetFillColorAlpha(46,0.35);
    SystZVtxPosition.at(0) -> GetYaxis() -> SetRangeUser(0,400);
    SystZVtxPosition.at(0) -> Draw("E2 SAME");
    StackBgrZVtxPosition -> Draw("HIST SAME");
    SelZVtxPosition.at(2) -> Draw("E2 SAME");
    SelZVtxPosition.at(0) -> SetLineColor(1);
    SelZVtxPosition.at(0) -> Draw("SAME");
    gPad->RedrawAxis();
    C8->SaveAs("ForwardFoldedZVtxPosition.pdf");

    //     SelectionTrackRange.at(0) -> Draw(); // 1200
//     SelectionCosTheta.at(0) -> Draw(); // 1400
//     SelectionTheta.at(0) -> Draw(); // 700
//     SelectionPhi.at(0) -> Draw(); // 400
//     SelectionMomentum.at(0) -> Draw(); // 1400
//     SelectionTrackLength.at(0) -> Draw(); // 1000
//     SelXVtxPosition.at(0) -> Draw(); // 350
//     SelYVtxPosition.at(0) -> Draw(); // 350
//     SelZVtxPosition.at(0) -> Draw(); // 350


//     TCanvas *C1 = new TCanvas("C1", "C1", 700, 500);
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
