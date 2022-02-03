#include <algorithm>
#include <array>
#include <fstream>
#include <sstream>
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

TSpline3* KEvsRSpline; // Global spline for momentum calculation

// Read Beam systematics for nu_mu, anti-nu_mu, nu_e, and anti-nu_e into a TGraph
std::vector<TGraph*> ReadFluxSystematics(const std::string &PathToFile);

// Function which calculates the distance between two points
float CalcRange(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2);

// Function which checks if a point is in the FV
bool inFV(double x, double y, double z);

// Function which checks if a point is in the TPC
bool inTPC(double x, double y, double z);

// Normalize Matrix by row
void NormMatrixByColumn(TH2F* UMatrix);

// Prepare Tracklength -> Momentum conversion Spline
void MomentumSplinePreparation();

// Get Momentum
float GetMomentum(float TrackLength);

void CalcSigEfficiency(std::vector<TH1F*>& HistVector);

// CC inlcusive cross section function (main)
void StoreHistograms()
{
    float NumberOfTargets = (FVx - 2*borderx) * (FVy - 2*bordery) * (FVz - 2*borderz) * Density * Avogadro/ArMass*NoNucleons;

    std::string Folder = "/home/christoph/anatrees/CCInclusiveNote";
//     std::string Folder = "/home/christoph/anatrees/ThesisSelection";

    // Output file file type
    std::string FileType = "pdf";
//     std::string FileType = "png";

    // Data input file vector
    std::vector<TChain*> ChainVec;

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

    // Efficiency graphs
    //TH1F* EffTrackRange;
    //TH1F* EffCosTheta;
    //TH1F* EffPhi;
    //TH1F* EffMomentum;

    size_t NumberOfBins = 20;

//     double MCPOT = 2.3e20/191362*92498;
    double MCPOT = 2.304e20; ///141*62;
    double TruthPOT = 5.451e19;
    double DataPOT = 4.950e19;

    double IntegratedFlux;

    int TrkID;
    int VtxID;
    int MCTrkID;
    int MCVtxID;

    int CCNCFlag[10];
    int TruthMode[10];
    float TrueLeptonMomentum[10];
    float NuEnergyTruth[10];
    float nuvtxx_truth[10]; //true vertex x (in cm)
    float nuvtxy_truth[10];
    float nuvtxz_truth[10];

    int PDGTruth[5000];
    short TrkBestPlane[5000];
    short TrkOrigin[5000][3];
    int TrackIDTruth[5000][3];

    float TrackTheta[5000];
    float TrackPhi[5000];
//     float TrackMomentum[5000];
    float TrackLength[5000];

    float XTrackStart[5000];
    float YTrackStart[5000];
    float ZTrackStart[5000];

    float XTrackEnd[5000];
    float YTrackEnd[5000];
    float ZTrackEnd[5000];

    float XVertexPosition[500];
    float YVertexPosition[500];
    float ZVertexPosition[500];

    //MC truth
    int mcevts_truth; //neutrino interactions per event
    float XnuVtxTruth[10]; //true vertex x (in cm)
    float YnuVtxTruth[10];
    float ZnuVtxTruth[10];
    int nuPDGTruth[10]; //true neutrino pdg code. numu = 14

    int NumberOfMCTracks;
    int MCTrackID[5000];
    int MCTrueIndex[5000];

    float XMCTrackStart[5000];
    float YMCTrackStart[5000];
    float ZMCTrackStart[5000];

    float XMCTrackEnd[5000];
    float YMCTrackEnd[5000];
    float ZMCTrackEnd[5000];

    float MCTheta[5000];
    float MCPhi[5000];
    float MCEnergy[5000];

    double HistogramWeight;

    // Calculating total efficiency
    unsigned int ExpectedEvents = 0;
    unsigned int SelectedEvents = 0;
    unsigned int CheckIfSane = 0;

    // Check MA normalization
    unsigned int MASamples = 0;
    double MASamplesCorr = 0;

    TH1D* NuMuFlux;
//     TH1F* NuMuFlux = new TH1F("NuMuFlux","NuMuFlux",NumberOfBins,0,3);

    TFile* BNBFlux = new TFile("/home/christoph/anatrees/BNBFlux/numode_bnb_470m_r200.root");

    NuMuFlux = (TH1D*) BNBFlux->Get("numu");

    IntegratedFlux = NuMuFlux->Integral()*4.95e19/1e20;

//     std::cout << TempNuMuFlux->Integral(5,NuMuFlux->GetNbinsX())*4.95e19/1e20 << std::endl;
    std::cout << "Integrated flux corresponding to 4.95e19 POT: " << IntegratedFlux << " cm^-2 s^-1"<< std::endl;
    
    std::string FluxSysFile = "bnb_sys_error_uboone.txt";
    // Read beam systematics
    std::vector<TGraph*> Test = ReadFluxSystematics(FluxSysFile);
    
//     TempNuMuFlux->Rebin(3);

//     for(unsigned int bin_no = 1; bin_no <= NuMuFlux->GetNbinsX(); bin_no++)
//     {
//         NuMuFlux->SetBinContent(bin_no, TempNuMuFlux->GetBinContent(bin_no));
//         NuMuFlux->SetBinError(bin_no, TempNuMuFlux->GetBinError(bin_no));
//     }

//     NuMuFlux->Scale(DataPOT/1e20);

    // Number of Events
    std::vector<unsigned int> NumberOfEvents;

    // Name vector
    std::vector<std::string> GenLabel;

    // Scaling vector
    std::vector<float> ScalingFactors;

    // Fill momentum calculation spline
    MomentumSplinePreparation();

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_1_Mod.root").c_str());
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_2_Mod.root").c_str());
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_3_Mod.root").c_str());
    NumberOfEvents.push_back(ChainVec.back() -> GetEntries());
    GenLabel.push_back("Data On-Beam BNB");
    ScalingFactors.push_back(1);

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_1_Mod.root").c_str());
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_2_Mod.root").c_str());
    NumberOfEvents.push_back(ChainVec.back() -> GetEntries());
    GenLabel.push_back("Data Off-Beam BNBEXT");
    ScalingFactors.push_back(1.2300);

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_prodgenie_bnb_nu_cosmic_uboone_v05_08_00_Mod.root").c_str());
//     ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_prodgenie_bnb_nu_cosmic_uboone_field_v05_08_00_Mod.root").c_str());
    NumberOfEvents.push_back(ChainVec.back() -> GetEntries());
    GenLabel.push_back("MC Selection");
    ScalingFactors.push_back(DataPOT/MCPOT);

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_MA_v05_08_00_Mod.root").c_str());
    NumberOfEvents.push_back(ChainVec.back() -> GetEntries());
    GenLabel.push_back("MA Adjusted Selection");
    ScalingFactors.push_back(DataPOT/MCPOT*NumberOfEvents.at(2)/NumberOfEvents.back());

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_TEM_v05_08_00_Mod.root").c_str());
    NumberOfEvents.push_back(ChainVec.back() -> GetEntries());
    GenLabel.push_back("TEM Selection");
    ScalingFactors.push_back(DataPOT/MCPOT);

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_MEC_v05_08_00_Mod.root").c_str());
    NumberOfEvents.push_back(ChainVec.back() -> GetEntries());
    GenLabel.push_back("MEC Selection");
    ScalingFactors.push_back(DataPOT/MCPOT);

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_MC_Truth_prodgenie_bnb_nu_cosmic_uboone_v05_08_00.root").c_str());
    NumberOfEvents.push_back(ChainVec.back() -> GetEntries());
    GenLabel.push_back("MC Truth");
    ScalingFactors.push_back(DataPOT/MCPOT);
//     ScalingFactors.push_back(DataPOT/TruthPOT);

//     GenLabel.push_back("Efficiency");
//     ScalingFactors.push_back(1);

    // Loop over all generation labels
    for(const auto& Label : GenLabel)
    {
        SelectionTrackRange.push_back(new TH1F(("Track Range "+Label).c_str(),"Muon Track Range",NumberOfBins,0,1036.8));
        SelectionTrackRange.back() -> SetStats(0);
        SelectionTrackRange.back() -> GetXaxis() -> SetTitle("Muon track range [cm]");
        SelectionTrackRange.back() -> GetYaxis() -> SetTitle("No. of events");
//         SelectionTrackRange.back() -> GetYaxis() -> SetTitle("d#sigma/dl [cm^{2}/cm]");

        SelectionCosTheta.push_back(new TH1F(("cos#theta "+Label).c_str(),"Cosine of #theta-Angle",NumberOfBins,-1,1));
        SelectionCosTheta.back() -> SetStats(0);
        SelectionCosTheta.back() -> GetXaxis() -> SetTitle("Muon cos(#theta)");
        SelectionCosTheta.back() -> GetYaxis() -> SetTitle("No. of events");
//         SelectionCosTheta.back() -> GetYaxis() -> SetTitle("d#sigma/d(cos#theta) [cm^{2}/cos(#theta)]");

        SelectionTheta.push_back(new TH1F(("#theta-Angle "+Label).c_str(),"#theta-Angle",NumberOfBins,0,180));
        SelectionTheta.back() -> SetStats(0);
        SelectionTheta.back() -> GetXaxis() -> SetTitle("Muon #theta-Angle [#circ]");
        SelectionTheta.back() -> GetYaxis() -> SetTitle("No. of events");
//         SelectionTheta.back() -> GetYaxis() -> SetTitle("d#sigma/d#theta [cm^{2}/rad]");

        SelectionPhi.push_back(new TH1F(("#phi-Angle "+Label).c_str(),"#varphi-Angle",NumberOfBins,-180,180));
        SelectionPhi.back() -> SetStats(0);
        SelectionPhi.back() -> GetXaxis() -> SetTitle("Muon #varphi-Angle [#circ]");
        SelectionPhi.back() -> GetYaxis() -> SetTitle("No. of events");
//         SelectionPhi.back() -> GetYaxis() -> SetTitle("d#sigma/d#phi [cm^{2}/#circ]");

        SelectionMomentum.push_back(new TH1F(("Momentum "+Label).c_str(),"Muon Momentum",NumberOfBins,0,3));
        SelectionMomentum.back() -> SetStats(0);
        SelectionMomentum.back() -> GetXaxis() -> SetTitle("Muon Momentum p_{#mu} [GeV/c]");
        SelectionMomentum.back() -> GetYaxis() -> SetTitle("No. of events");
//         SelectionMomentum.back() -> GetYaxis() -> SetTitle("d#sigma/dp [cm^{2}/(GeV/c)]");

        SelectionTrackLength.push_back(new TH1F(("Track Length "+Label).c_str(),"Candidate Track Length",NumberOfBins,0,1036.8));
        SelectionTrackLength.back() -> SetStats(0);
        SelectionTrackLength.back() -> GetXaxis() -> SetTitle("Muon Track Length l_ [cm]");
        SelectionTrackLength.back() -> GetYaxis() -> SetTitle("No. of events");

        SelXVtxPosition.push_back(new TH1F(("Vertex X position "+Label).c_str(),"Vertex Position in X",NumberOfBins,0,256.35));
        SelXVtxPosition.back() -> SetStats(0);
        SelXVtxPosition.back() -> GetXaxis() -> SetTitle("x-coordinate [cm]");
        SelXVtxPosition.back() -> GetYaxis() -> SetTitle("No. of events");

        SelYVtxPosition.push_back(new TH1F(("Vertex Y position "+Label).c_str(),"Vertex Position in Y",NumberOfBins,-233*0.5,233*0.5));
        SelYVtxPosition.back() -> SetStats(0);
        SelYVtxPosition.back() -> GetXaxis() -> SetTitle("y-coordinate [cm]");
        SelYVtxPosition.back() -> GetYaxis() -> SetTitle("No. of events");

        SelZVtxPosition.push_back(new TH1F(("Vertex Z position "+Label).c_str(),"Vertex Position in Z",NumberOfBins,0,1036.8));
        SelZVtxPosition.back() -> SetStats(0);
        SelZVtxPosition.back() -> GetXaxis() -> SetTitle("z-coordinate [cm]");
        SelZVtxPosition.back() -> GetYaxis() -> SetTitle("No. of events");

        // Calculate standard deviation for all histograms
        SelectionTrackRange.back()->Sumw2();
        SelectionCosTheta.back()->Sumw2();
        SelectionTheta.back()->Sumw2();
        SelectionPhi.back()->Sumw2();
        SelectionMomentum.back()->Sumw2();
        SelectionTrackLength.back()->Sumw2();
        SelXVtxPosition.back()->Sumw2();
        SelYVtxPosition.back()->Sumw2();
        SelZVtxPosition.back()->Sumw2();


    } // loop over generation label

    // Initialize smearing matrices
    UMatrixTrackRange = new TH2F("Unsmering Matrix Track Range","Smearing Matrix Track Range",NumberOfBins,0,1036.8,NumberOfBins,0,1036.8);
    UMatrixTrackRange -> GetXaxis() -> SetTitle("Muon track length (true) [cm]");
    UMatrixTrackRange -> GetYaxis() -> SetTitle("Muon track length (reco) [cm]");
    UMatrixCosTheta = new TH2F("Smearing Matrix CosTheta","Smearing Matrix cos(#theta)",NumberOfBins,-1,1,NumberOfBins,-1,1);
    UMatrixCosTheta -> GetXaxis() -> SetTitle("Muon cos(#theta) (true)");
    UMatrixCosTheta -> GetYaxis() -> SetTitle("Muon cos(#theta) (reco)");
    UMatrixTheta = new TH2F("Smearing Matrix Theta","Smearing Matrix Theta",NumberOfBins,0,180,NumberOfBins,0,180);
    UMatrixTheta -> GetXaxis() -> SetTitle("Muon #theta-Angle (true) [#circ]");
    UMatrixTheta -> GetYaxis() -> SetTitle("Muon #theta-Angle (reco) [#circ]");
    UMatrixPhi = new TH2F("Smearing Matrix Phi","Smearing Matrix Phi",NumberOfBins,-180,180,NumberOfBins,-180,180);
    UMatrixPhi -> GetXaxis() -> SetTitle("Muon #varphi-Angle (true) [#circ]");
    UMatrixPhi -> GetYaxis() -> SetTitle("Muon #varphi-Angle (reco) [#circ]");
    UMatrixMomentum = new TH2F("Smearing Matrix Momentum","Smearing Matrix Momentum",NumberOfBins,0,3,NumberOfBins,0,3);
    UMatrixMomentum -> GetXaxis() -> SetTitle("Muon Momentum (true) [GeV/c]");
    UMatrixMomentum -> GetYaxis() -> SetTitle("Muon Momentum (reco) [GeV/c]");
    UMatrixXVtxPosition = new TH2F("Unsmering Matrix XVtx","Smearing Matrix XVtx",NumberOfBins,0,256.35,NumberOfBins,0,256.35);
    UMatrixXVtxPosition -> GetXaxis() -> SetTitle("Vertex x-coordinate (true) [cm]");
    UMatrixXVtxPosition -> GetYaxis() -> SetTitle("Vertex x-coordinate (reco) [cm]");
    UMatrixYVtxPosition = new TH2F("Unsmering Matrix YVtx","Smearing Matrix YVtx",NumberOfBins,-233*0.5,233*0.5,NumberOfBins,-233*0.5,233*0.5);
    UMatrixYVtxPosition -> GetXaxis() -> SetTitle("Vertex y-coordinate (true) [cm]");
    UMatrixYVtxPosition -> GetYaxis() -> SetTitle("Vertex y-coordinate (reco) [cm]");
    UMatrixZVtxPosition = new TH2F("Unsmering Matrix ZVtx","Smearing Matrix ZVtx",NumberOfBins,0,1036.8,NumberOfBins,0,1036.8);
    UMatrixZVtxPosition -> GetXaxis() -> SetTitle("Vertex z-coordinate (true) [cm]");
    UMatrixZVtxPosition -> GetYaxis() -> SetTitle("Vertex z-coordinate (reco) [cm]");

    // MC Background
    std::vector<std::string> BgrLabel;
    BgrLabel.push_back("All");
    BgrLabel.push_back("cosmic");
    BgrLabel.push_back("dirt");
    BgrLabel.push_back("outFV");
    BgrLabel.push_back("anti nu_mu");
    BgrLabel.push_back("n_e-like");
    BgrLabel.push_back("nu_NC");
    BgrLabel.push_back("PureSelected");

    // Initialize vector
    BgrTrackRange.resize(4);
    BgrCosTheta.resize(4);
    BgrTheta.resize(4);
    BgrPhi.resize(4);
    BgrMomentum.resize(4);
    BgrTrackLength.resize(4);
    BgrXVtxPosition.resize(4);
    BgrYVtxPosition.resize(4);
    BgrZVtxPosition.resize(4);

    for(unsigned int file_no = 0; file_no < BgrTrackRange.size(); file_no++)
    {
        for(auto Label : BgrLabel)
        {
            BgrTrackRange.at(file_no).push_back(new TH1F((Label+"Background Range "+std::to_string(file_no)).c_str(),"Range",NumberOfBins,0,1036.8));
            BgrCosTheta.at(file_no).push_back(new TH1F((Label+"Background cos#theta "+std::to_string(file_no)).c_str(),"cos#theta",NumberOfBins,-1,1));
            BgrTheta.at(file_no).push_back(new TH1F((Label+"Background #theta "+std::to_string(file_no)).c_str(),"#theta",NumberOfBins,0,180));
            BgrPhi.at(file_no).push_back(new TH1F((Label+"Background #phi "+std::to_string(file_no)).c_str(),"#phi",NumberOfBins,-180,180));
            BgrMomentum.at(file_no).push_back(new TH1F((Label+"Background Momentum "+std::to_string(file_no)).c_str(),"Momentum",NumberOfBins,0,3));
            BgrTrackLength.at(file_no).push_back(new TH1F((Label+"Background Length "+std::to_string(file_no)).c_str(),"Lenght",NumberOfBins,0,1036.8));
            BgrXVtxPosition.at(file_no).push_back(new TH1F((Label+"Background XVtx "+std::to_string(file_no)).c_str(),"XVtx",NumberOfBins,0,256.35));
            BgrYVtxPosition.at(file_no).push_back(new TH1F((Label+"Background YVtx "+std::to_string(file_no)).c_str(),"YVtx",NumberOfBins,-233*0.5,233*0.5));
            BgrZVtxPosition.at(file_no).push_back(new TH1F((Label+"Background ZVtx "+std::to_string(file_no)).c_str(),"ZVtx",NumberOfBins,0,1036.8));

            BgrTrackRange.at(file_no).back() -> Sumw2();
            BgrCosTheta.at(file_no).back() -> Sumw2();
            BgrTheta.at(file_no).back() -> Sumw2();
            BgrPhi.at(file_no).back() -> Sumw2();
            BgrMomentum.at(file_no).back() -> Sumw2();
            BgrTrackLength.at(file_no).back() -> Sumw2();
            BgrXVtxPosition.at(file_no).back() -> Sumw2();
            BgrYVtxPosition.at(file_no).back() -> Sumw2();
            BgrZVtxPosition.at(file_no).back() -> Sumw2();
        }
    }

    // Loop over all files
    for(unsigned int file_no = 0; file_no < ChainVec.size(); file_no++)
    {
        std::cout << "-------------File Progress--------------" << std::endl;
        std::cout << "File \t \t" << file_no+1 << " of " << ChainVec.size() << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        // Reco entities for all files except truth
        if(file_no < 6)
        {
            ChainVec.at(file_no) -> SetBranchAddress("TrackCand", &TrkID);
            ChainVec.at(file_no) -> SetBranchAddress("VertexCand", &VtxID);

//             ChainVec.at(file_no) -> SetBranchAddress("trkmomrange_pandoraNu", TrackMomentum);
            ChainVec.at(file_no) -> SetBranchAddress("trklen_pandoraNu", TrackLength);
            ChainVec.at(file_no) -> SetBranchAddress("trktheta_pandoraNu", TrackTheta);
            ChainVec.at(file_no) -> SetBranchAddress("trkphi_pandoraNu",TrackPhi);

            ChainVec.at(file_no) -> SetBranchAddress("trkstartx_pandoraNu",XTrackStart);
            ChainVec.at(file_no) -> SetBranchAddress("trkstarty_pandoraNu",YTrackStart);
            ChainVec.at(file_no) -> SetBranchAddress("trkstartz_pandoraNu",ZTrackStart);

            ChainVec.at(file_no) -> SetBranchAddress("trkendx_pandoraNu",XTrackEnd);
            ChainVec.at(file_no) -> SetBranchAddress("trkendy_pandoraNu",YTrackEnd);
            ChainVec.at(file_no) -> SetBranchAddress("trkendz_pandoraNu",ZTrackEnd);

            ChainVec.at(file_no) -> SetBranchAddress("vtxx_pandoraNu", XVertexPosition);
            ChainVec.at(file_no) -> SetBranchAddress("vtxy_pandoraNu", YVertexPosition);
            ChainVec.at(file_no) -> SetBranchAddress("vtxz_pandoraNu", ZVertexPosition);
        }

        // MC entities just for non-data files
        if(file_no > 1)
        {
            ChainVec.at(file_no) -> SetBranchAddress("ccnc_truth", CCNCFlag);
            ChainVec.at(file_no) -> SetBranchAddress("mode_truth", TruthMode);
            ChainVec.at(file_no) -> SetBranchAddress("pdg", PDGTruth);
            ChainVec.at(file_no) -> SetBranchAddress("enu_truth", NuEnergyTruth);
            ChainVec.at(file_no) -> SetBranchAddress("lep_mom_truth", TrueLeptonMomentum);
            ChainVec.at(file_no) -> SetBranchAddress("mcevts_truth", &mcevts_truth);
            ChainVec.at(file_no) -> SetBranchAddress("nuvtxx_truth", XnuVtxTruth);
            ChainVec.at(file_no) -> SetBranchAddress("nuvtxy_truth", YnuVtxTruth);
            ChainVec.at(file_no) -> SetBranchAddress("nuvtxz_truth", ZnuVtxTruth);
            ChainVec.at(file_no) -> SetBranchAddress("nuPDG_truth", nuPDGTruth);
            ChainVec.at(file_no) -> SetBranchAddress("geant_list_size", &NumberOfMCTracks);
            ChainVec.at(file_no) -> SetBranchAddress("MCTruthIndex", MCTrueIndex);
            ChainVec.at(file_no) -> SetBranchAddress("StartPointx", XMCTrackStart);
            ChainVec.at(file_no) -> SetBranchAddress("StartPointy", YMCTrackStart);
            ChainVec.at(file_no) -> SetBranchAddress("StartPointz", ZMCTrackStart);
            ChainVec.at(file_no) -> SetBranchAddress("EndPointx", XMCTrackEnd);
            ChainVec.at(file_no) -> SetBranchAddress("EndPointy", YMCTrackEnd);
            ChainVec.at(file_no) -> SetBranchAddress("EndPointz", ZMCTrackEnd);
            ChainVec.at(file_no) -> SetBranchAddress("enu_truth", NuEnergyTruth);
            ChainVec.at(file_no) -> SetBranchAddress("theta", MCTheta);
            ChainVec.at(file_no) -> SetBranchAddress("phi", MCPhi);
            ChainVec.at(file_no) -> SetBranchAddress("Eng", MCEnergy);

            if(file_no < 6)
            {
                ChainVec.at(file_no) -> SetBranchAddress("TrackId", MCTrackID);
                ChainVec.at(file_no) -> SetBranchAddress("trkorigin_pandoraNu", TrkOrigin);
                ChainVec.at(file_no) -> SetBranchAddress("trkidtruth_pandoraNu",TrackIDTruth);
                ChainVec.at(file_no) -> SetBranchAddress("trkpidbestplane_pandoraNu", TrkBestPlane);
            }
            else // if(file_no == 6)
            {
                ChainVec.at(file_no) -> SetBranchAddress("MCTrackCand", &MCTrkID);
                ChainVec.at(file_no) -> SetBranchAddress("MCVertexCand", &MCVtxID);
            }
        }

        // MA re-weight factor
        if(file_no == 3)
        {
            ChainVec.at(file_no) -> SetBranchAddress("eventWeight_MA", &HistogramWeight);
        }
        else
        {
            HistogramWeight = 1.0;
        }

        // Background counters and their errors
        unsigned int allbgr = 0;
        unsigned int cosmic = 0;
        unsigned int dirt = 0;
        unsigned int nuOutFV = 0;
        unsigned int antinu_mu = 0;
        unsigned int nu_e = 0;
        unsigned int nuNC = 0;

        // Calculating total efficiency
//         unsigned int ExpectedEvents = 0;
//         unsigned int SelectedEvents = 0;

        // Loop over all events
        for(unsigned int tree_index = 0; tree_index < ChainVec.at(file_no) -> GetEntries(); tree_index++)
        {
            // Progress indicator
            if(!(tree_index % 1000)) std::cout << "Event\t" << tree_index << "\t of \t" << ChainVec.at(file_no) -> GetEntries() << std::endl;

            // Skip corrupted events in files file
            if(file_no == 5 && tree_index == 5455) continue;
            if(file_no == 6 && (tree_index == 11602 || tree_index == 11675 || tree_index == 13510 || tree_index == 33027 || tree_index == 33070 || tree_index == 36239 || tree_index == 44078)) continue;

            // Get tree entry for this event
            ChainVec.at(file_no) -> GetEntry(tree_index);

            // if there are reco products
            if(file_no < 6)
            {
                // Fill histograms as usual for all files
                SelectionTrackRange.at(file_no) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),HistogramWeight);
                SelectionCosTheta.at(file_no) -> Fill(std::cos(TrackTheta[TrkID]),HistogramWeight);
                SelectionTheta.at(file_no) -> Fill(TrackTheta[TrkID]/Pi*180,HistogramWeight);
                SelectionPhi.at(file_no) -> Fill(TrackPhi[TrkID]/Pi*180,HistogramWeight);
                SelectionMomentum.at(file_no) -> Fill(GetMomentum(TrackLength[TrkID]),HistogramWeight);
                SelectionTrackLength.at(file_no) -> Fill(TrackLength[TrkID],HistogramWeight);
                SelXVtxPosition.at(file_no) -> Fill(XVertexPosition[VtxID],HistogramWeight);
                SelYVtxPosition.at(file_no) -> Fill(YVertexPosition[VtxID],HistogramWeight);
                SelZVtxPosition.at(file_no) -> Fill(ZVertexPosition[VtxID],HistogramWeight);
            }

            // if file is a MC Selection
            if(file_no > 1 && file_no < 6)
            {
                // Match track and vertex to true information
                int MCTrkID = -1;
                int MCVtxID = -1;

                // If the candidate is a neutrino track
                if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1)
                {
                    // Loop over all MCTracks
                    for(unsigned track_no = 0; track_no < NumberOfMCTracks; track_no++)
                    {
                        // If the Track ID of the neutrino track is the same as the MC truth ID
                        if(MCTrackID[track_no] == TrackIDTruth[TrkID][TrkBestPlane[TrkID]])
                        {
                            // Store MC truth vertex and track ID
                            MCTrkID = track_no;
                            MCVtxID = MCTrueIndex[track_no];
                        }
                    }
                }

                // Count Histogram weight
                if(file_no == 3)
                {
                    MASamples++;
                    MASamplesCorr += HistogramWeight;
                }

                // If track origin isn't a neutrino or not nu_mu or an NC event or not in FV it is background
                if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] != 1 || nuPDGTruth[MCVtxID] != 14 || CCNCFlag[MCVtxID] == 1 || !inFV(XnuVtxTruth[MCVtxID],YnuVtxTruth[MCVtxID],ZnuVtxTruth[MCVtxID]))
                {
                    BgrTrackRange.at(file_no-2).at(0) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),HistogramWeight);
                    BgrCosTheta.at(file_no-2).at(0) -> Fill(std::cos(TrackTheta[TrkID]),HistogramWeight);
                    BgrTheta.at(file_no-2).at(0) -> Fill(TrackTheta[TrkID]/Pi*180,HistogramWeight);
                    BgrPhi.at(file_no-2).at(0) -> Fill(TrackPhi[TrkID]/Pi*180);
                    BgrMomentum.at(file_no-2).at(0) -> Fill(GetMomentum(TrackLength[TrkID]),HistogramWeight);
                    BgrTrackLength.at(file_no-2).at(0)-> Fill(TrackLength[TrkID],HistogramWeight);;
                    BgrXVtxPosition.at(file_no-2).at(0) -> Fill(XVertexPosition[VtxID],HistogramWeight);
                    BgrYVtxPosition.at(file_no-2).at(0) -> Fill(YVertexPosition[VtxID],HistogramWeight);
                    BgrZVtxPosition.at(file_no-2).at(0) -> Fill(ZVertexPosition[VtxID],HistogramWeight);

                    allbgr++;
                }

                // If track origin is not neutrino
                if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] != 1)
                {
                    BgrTrackRange.at(file_no-2).at(1) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),HistogramWeight);
                    BgrCosTheta.at(file_no-2).at(1) -> Fill(std::cos(TrackTheta[TrkID]),HistogramWeight);
                    BgrTheta.at(file_no-2).at(1) -> Fill(TrackTheta[TrkID]/Pi*180,HistogramWeight);
                    BgrPhi.at(file_no-2).at(1) -> Fill(TrackPhi[TrkID]/Pi*180,HistogramWeight);
                    BgrMomentum.at(file_no-2).at(1) -> Fill(GetMomentum(TrackLength[TrkID]),HistogramWeight);
                    BgrTrackLength.at(file_no-2).at(1)-> Fill(TrackLength[TrkID],HistogramWeight);;
                    BgrXVtxPosition.at(file_no-2).at(1) -> Fill(XVertexPosition[VtxID],HistogramWeight);
                    BgrYVtxPosition.at(file_no-2).at(1) -> Fill(YVertexPosition[VtxID],HistogramWeight);
                    BgrZVtxPosition.at(file_no-2).at(1) -> Fill(ZVertexPosition[VtxID],HistogramWeight);

                    cosmic++;
                }
                // else if neutrino but not in TPC
                else if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1 && !inTPC(XnuVtxTruth[MCVtxID],YnuVtxTruth[MCVtxID],ZnuVtxTruth[MCVtxID]))
                {
                    BgrTrackRange.at(file_no-2).at(2) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),HistogramWeight);
                    BgrCosTheta.at(file_no-2).at(2) -> Fill(std::cos(TrackTheta[TrkID]),HistogramWeight);
                    BgrTheta.at(file_no-2).at(2) -> Fill(TrackTheta[TrkID]/Pi*180,HistogramWeight);
                    BgrPhi.at(file_no-2).at(2) -> Fill(TrackPhi[TrkID]/Pi*180,HistogramWeight);
                    BgrMomentum.at(file_no-2).at(2) -> Fill(GetMomentum(TrackLength[TrkID]),HistogramWeight);
                    BgrTrackLength.at(file_no-2).at(2)-> Fill(TrackLength[TrkID],HistogramWeight);;
                    BgrXVtxPosition.at(file_no-2).at(2) -> Fill(XVertexPosition[VtxID],HistogramWeight);
                    BgrYVtxPosition.at(file_no-2).at(2) -> Fill(YVertexPosition[VtxID],HistogramWeight);
                    BgrZVtxPosition.at(file_no-2).at(2) -> Fill(ZVertexPosition[VtxID],HistogramWeight);

                    dirt++;
                }
                // else if not in FV (excluding out of TPC because of else if)
                else if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1 && !inFV(XnuVtxTruth[MCVtxID],YnuVtxTruth[MCVtxID],ZnuVtxTruth[MCVtxID]))
                {
                    BgrTrackRange.at(file_no-2).at(3) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),HistogramWeight);
                    BgrCosTheta.at(file_no-2).at(3) -> Fill(std::cos(TrackTheta[TrkID]),HistogramWeight);
                    BgrTheta.at(file_no-2).at(3) -> Fill(TrackTheta[TrkID]/Pi*180,HistogramWeight);
                    BgrPhi.at(file_no-2).at(3) -> Fill(TrackPhi[TrkID]/Pi*180,HistogramWeight);
                    BgrMomentum.at(file_no-2).at(3) -> Fill(GetMomentum(TrackLength[TrkID]),HistogramWeight);
                    BgrTrackLength.at(file_no-2).at(3)-> Fill(TrackLength[TrkID],HistogramWeight);;
                    BgrXVtxPosition.at(file_no-2).at(3) -> Fill(XVertexPosition[VtxID],HistogramWeight);
                    BgrYVtxPosition.at(file_no-2).at(3) -> Fill(YVertexPosition[VtxID],HistogramWeight);
                    BgrZVtxPosition.at(file_no-2).at(3) -> Fill(ZVertexPosition[VtxID],HistogramWeight);

                    nuOutFV++;
                }
                // If Origin is neutrino & CC event & interaction product is anti-nu_mu
                else if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1 && CCNCFlag[MCVtxID] == 0 && nuPDGTruth[MCVtxID] == -14)
                {
                    BgrTrackRange.at(file_no-2).at(4) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),HistogramWeight);
                    BgrCosTheta.at(file_no-2).at(4) -> Fill(std::cos(TrackTheta[TrkID]),HistogramWeight);
                    BgrTheta.at(file_no-2).at(4) -> Fill(TrackTheta[TrkID]/Pi*180,HistogramWeight);
                    BgrPhi.at(file_no-2).at(4) -> Fill(TrackPhi[TrkID]/Pi*180,HistogramWeight);
                    BgrMomentum.at(file_no-2).at(4) -> Fill(GetMomentum(TrackLength[TrkID]),HistogramWeight);
                    BgrTrackLength.at(file_no-2).at(4)-> Fill(TrackLength[TrkID],HistogramWeight);;
                    BgrXVtxPosition.at(file_no-2).at(4) -> Fill(XVertexPosition[VtxID],HistogramWeight);
                    BgrYVtxPosition.at(file_no-2).at(4) -> Fill(YVertexPosition[VtxID],HistogramWeight);
                    BgrZVtxPosition.at(file_no-2).at(4) -> Fill(ZVertexPosition[VtxID],HistogramWeight);

                    antinu_mu++;
                }
                // else if nu_e like event
                else if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1 && CCNCFlag[MCVtxID] == 0 && std::abs(nuPDGTruth[MCVtxID]) == 12)
                {
                    BgrTrackRange.at(file_no-2).at(5) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),HistogramWeight);
                    BgrCosTheta.at(file_no-2).at(5) -> Fill(std::cos(TrackTheta[TrkID]),HistogramWeight);
                    BgrTheta.at(file_no-2).at(5) -> Fill(TrackTheta[TrkID]/Pi*180,HistogramWeight);
                    BgrPhi.at(file_no-2).at(5) -> Fill(TrackPhi[TrkID]/Pi*180,HistogramWeight);
                    BgrMomentum.at(file_no-2).at(5) -> Fill(GetMomentum(TrackLength[TrkID]),HistogramWeight);
                    BgrTrackLength.at(file_no-2).at(5)-> Fill(TrackLength[TrkID],HistogramWeight);;
                    BgrXVtxPosition.at(file_no-2).at(5) -> Fill(XVertexPosition[VtxID],HistogramWeight);
                    BgrYVtxPosition.at(file_no-2).at(5) -> Fill(YVertexPosition[VtxID],HistogramWeight);
                    BgrZVtxPosition.at(file_no-2).at(5) -> Fill(ZVertexPosition[VtxID],HistogramWeight);

                    nu_e++;
                }
                // else if neutral current event
                else if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1 && CCNCFlag[MCVtxID] == 1)
                {
                    BgrTrackRange.at(file_no-2).at(6) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),HistogramWeight);
                    BgrCosTheta.at(file_no-2).at(6) -> Fill(std::cos(TrackTheta[TrkID]),HistogramWeight);
                    BgrTheta.at(file_no-2).at(6) -> Fill(TrackTheta[TrkID]/Pi*180,HistogramWeight);
                    BgrPhi.at(file_no-2).at(6) -> Fill(TrackPhi[TrkID]/Pi*180,HistogramWeight);
                    BgrMomentum.at(file_no-2).at(6) -> Fill(GetMomentum(TrackLength[TrkID]),HistogramWeight);
                    BgrTrackLength.at(file_no-2).at(6)-> Fill(TrackLength[TrkID],HistogramWeight);;
                    BgrXVtxPosition.at(file_no-2).at(6) -> Fill(XVertexPosition[VtxID],HistogramWeight);
                    BgrYVtxPosition.at(file_no-2).at(6) -> Fill(YVertexPosition[VtxID],HistogramWeight);
                    BgrZVtxPosition.at(file_no-2).at(6) -> Fill(ZVertexPosition[VtxID],HistogramWeight);

                    nuNC++;
                }
                // everything that is not background
                else
                {
                    // This background vector entry contains the pure signal, without background
                    BgrTrackRange.at(file_no-2).at(7) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),HistogramWeight);
                    BgrCosTheta.at(file_no-2).at(7) -> Fill(std::cos(TrackTheta[TrkID]),HistogramWeight);
                    BgrTheta.at(file_no-2).at(7) -> Fill(TrackTheta[TrkID]/Pi*180,HistogramWeight);
                    BgrPhi.at(file_no-2).at(7) -> Fill(TrackPhi[TrkID]/Pi*180,HistogramWeight);
                    BgrMomentum.at(file_no-2).at(7) -> Fill(GetMomentum(TrackLength[TrkID]),HistogramWeight);
                    BgrTrackLength.at(file_no-2).at(7)-> Fill(TrackLength[TrkID],HistogramWeight);;
                    BgrXVtxPosition.at(file_no-2).at(7) -> Fill(XVertexPosition[VtxID],HistogramWeight);
                    BgrYVtxPosition.at(file_no-2).at(7) -> Fill(YVertexPosition[VtxID],HistogramWeight);
                    BgrZVtxPosition.at(file_no-2).at(7) -> Fill(ZVertexPosition[VtxID],HistogramWeight);

                    // only if file 2 and nu_mu CC event
                    if(file_no == 2 && TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1 && nuPDGTruth[MCVtxID] == 14 && inFV(XnuVtxTruth[MCVtxID],YnuVtxTruth[MCVtxID],ZnuVtxTruth[MCVtxID]))
                    {
                        // Fill searing matrices
                        UMatrixTrackRange -> Fill( CalcRange(XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID]), CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]) );
                        UMatrixCosTheta -> Fill( std::cos(MCTheta[MCTrkID]), std::cos(TrackTheta[TrkID]) );
                        UMatrixTheta -> Fill( MCTheta[MCTrkID]/Pi*180, TrackTheta[TrkID]/Pi*180 );
                        UMatrixPhi -> Fill( MCPhi[MCTrkID]/Pi*180, TrackPhi[TrkID]/Pi*180 );
                        UMatrixMomentum -> Fill( TrueLeptonMomentum[MCVtxID], GetMomentum(TrackLength[TrkID]) );
                        UMatrixXVtxPosition -> Fill( XnuVtxTruth[MCVtxID], XVertexPosition[VtxID]);
                        UMatrixYVtxPosition -> Fill( YnuVtxTruth[MCVtxID], YVertexPosition[VtxID]);
                        UMatrixZVtxPosition -> Fill( ZnuVtxTruth[MCVtxID], ZVertexPosition[VtxID]);

                        SelectedEvents++;
                    }
                    if(file_no == 2) CheckIfSane++;
                }
            }
            // if truth selection file
            else if(file_no == 6 && MCTrkID >= 0 && nuPDGTruth[MCVtxID] == 14 && inFV(XnuVtxTruth[MCVtxID],YnuVtxTruth[MCVtxID],ZnuVtxTruth[MCVtxID]))
            {
                // Fill background histograms
                SelectionTrackRange.at(file_no) -> Fill(CalcRange(XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID]));
                SelectionCosTheta.at(file_no) -> Fill(std::cos(MCTheta[MCTrkID]));
                SelectionTheta.at(file_no) -> Fill(MCTheta[MCTrkID]/Pi*180);
                SelectionPhi.at(file_no) -> Fill(MCPhi[MCTrkID]/Pi*180);
                SelectionMomentum.at(file_no) -> Fill(TrueLeptonMomentum[MCVtxID]);
                SelXVtxPosition.at(file_no) -> Fill(XnuVtxTruth[MCVtxID]);
                SelYVtxPosition.at(file_no) -> Fill(YnuVtxTruth[MCVtxID]);
                SelZVtxPosition.at(file_no) -> Fill(ZnuVtxTruth[MCVtxID]);

                ExpectedEvents++;
            }

            // Selection sanity check
            if(CheckIfSane != SelectedEvents)
            {
                std::cout << "----------------------------------------" << std::endl;
                std::cout << "!----Background Sanity Check Failed----!" << std::endl;
                std::cout << "----------------------------------------" << std::endl;
                std::cout << "File: " << file_no << ", Event: " << tree_index << ", Selected: " << SelectedEvents << ", Selection Check: " << CheckIfSane << std::endl;
            }
        } // Event loop

        // Fill purities
        if(file_no > 1 && file_no < 6)
        {
            PurTrackRange.push_back( new TEfficiency(*BgrTrackRange.at(file_no-2).back(), *SelectionTrackRange.at(file_no)) );
            PurCosTheta.push_back( new TEfficiency(*BgrCosTheta.at(file_no-2).back(), *SelectionCosTheta.at(file_no)) );
            PurTheta.push_back( new TEfficiency(*BgrTheta.at(file_no-2).back(), *SelectionTheta.at(file_no)) );
            PurPhi.push_back( new TEfficiency(*BgrPhi.at(file_no-2).back(), *SelectionPhi.at(file_no)) );
            PurMomentum.push_back( new TEfficiency(*BgrMomentum.at(file_no-2).back(), *SelectionMomentum.at(file_no)) );
            PurTrackLength.push_back( new TEfficiency(*BgrTrackLength.at(file_no-2).back(), *SelectionTrackRange.at(file_no)) );
            PurXVtxPosition.push_back( new TEfficiency(*BgrXVtxPosition.at(file_no-2).back(), *SelXVtxPosition.at(file_no)) );
            PurYVtxPosition.push_back( new TEfficiency(*BgrYVtxPosition.at(file_no-2).back(), *SelYVtxPosition.at(file_no)) );
            PurZVtxPosition.push_back( new TEfficiency(*BgrZVtxPosition.at(file_no-2).back(), *SelZVtxPosition.at(file_no)) );
        }

        // Calculate Background errors
        double allbgrErr = std::sqrt(allbgr);
        double cosmicErr = std::sqrt(cosmic);
        double dirtErr = std::sqrt(dirt);
        double nuOutFVErr = std::sqrt(nuOutFV);
        double antinu_muErr = std::sqrt(antinu_mu);
        double nu_eErr = std::sqrt(nu_e);
        double nuNCErr = std::sqrt(nuNC);

        // Print background summary
        std::cout << "-----------Background Summary-----------" << std::endl;
        std::cout << "Total Bgr: \t" << allbgr << " ± " << allbgrErr << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "Cosmic bgr: \t" << cosmic << " ± " << cosmicErr << std::endl;
        std::cout << "Dirt event: \t" << dirt << " ± " << dirtErr << std::endl;
        std::cout << "Out of FV: \t" << nuOutFV << " ± " << nuOutFVErr << std::endl;
        std::cout << "anti nu_mu: \t" << antinu_mu << " ± " << antinu_muErr << std::endl;
        std::cout << "nu_e like: \t" << nu_e << " ± " << nu_eErr << std::endl;
        std::cout << "NC event: \t" << nuNC << " ± " << nuNCErr << std::endl;
        std::cout << "----------------------------------------" << std::endl;

        // Check if backgrounds add up correctly
        if(allbgr != cosmic + dirt + nuOutFV + antinu_mu + nu_e + nuNC)
        {
            std::cout << "WARNING: There is a difference in background sums by " << allbgr - cosmic - dirt - nuOutFV - antinu_mu - nu_e - nuNC << " events!" << std::endl;
        }

        // Reset branch addresses to avoid problems
        ChainVec.at(file_no) -> ResetBranchAddresses();

        // Delete chain from memory (preventing overflow)
        ChainVec.at(file_no) -> Delete();
    } // file loop

    std::cout << "---------MA Normalization Check---------" << std::endl;
    std::cout << "MA count normal \t" << MASamples  <<  std::endl;
    std::cout << "MA weighted count \t" << MASamplesCorr << std::endl;
    std::cout << "MA weight check \t" << SelectionTrackRange.at(6)->Integral() << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    double TotalEfficiency = (double)SelectedEvents/(double)ExpectedEvents;
    double EfficiencyError = TotalEfficiency*(std::sqrt((double)SelectedEvents)/(double)SelectedEvents +  std::sqrt((double)ExpectedEvents)/(double)ExpectedEvents);

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Overall selected events \t" << SelectedEvents << " ± " << std::sqrt((double)SelectedEvents)  <<  std::endl;
    std::cout << "Overall expected events \t" << ExpectedEvents << " ± " << std::sqrt((double)ExpectedEvents) << std::endl;
    std::cout << "Overall selection efficiency \t" << TotalEfficiency << " ± " << EfficiencyError << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    // Fill efficiency (could be a loop if there are more of them)
    {
        unsigned int file_no = 2;

        // Efficiencies
        EffTrackRange.push_back( new TEfficiency(*BgrTrackRange.at(file_no-2).back(), *SelectionTrackRange.back()) );
        EffCosTheta.push_back( new TEfficiency(*BgrCosTheta.at(file_no-2).back(), *SelectionCosTheta.back()) );
        EffTheta.push_back( new TEfficiency(*BgrTheta.at(file_no-2).back(), *SelectionTheta.back()) );
        EffPhi.push_back( new TEfficiency(*BgrPhi.at(file_no-2).back(), *SelectionPhi.back()) );
        EffMomentum.push_back( new TEfficiency(*BgrMomentum.at(file_no-2).back(), *SelectionMomentum.back()) );
        EffXVtxPosition.push_back( new TEfficiency(*BgrXVtxPosition.at(file_no-2).back(), *SelXVtxPosition.back()) );
        EffYVtxPosition.push_back( new TEfficiency(*BgrYVtxPosition.at(file_no-2).back(), *SelYVtxPosition.back()) );
        EffZVtxPosition.push_back( new TEfficiency(*BgrZVtxPosition.at(file_no-2).back(), *SelZVtxPosition.back()) );

        // Maybe more?
    }

    // loop over all histograms

    // Fill efficiency
//     SelectionTrackRange.back()->Divide(SelectionTrackRange.at(2),SelectionTrackRange.at(6));
//     SelectionCosTheta.back()->Divide(SelectionCosTheta.at(2),SelectionCosTheta.at(6));
//     SelectionTheta.back()->Divide(SelectionTheta.at(2),SelectionTheta.at(6));
//     SelectionPhi.back()->Divide(SelectionPhi.at(2),SelectionPhi.at(6));
//     SelectionMomentum.back()->Divide(SelectionMomentum.at(2),SelectionMomentum.at(6));

    // Normalize matrices by row
    NormMatrixByColumn(UMatrixTrackRange);
    NormMatrixByColumn(UMatrixCosTheta);
    NormMatrixByColumn(UMatrixTheta);
    NormMatrixByColumn(UMatrixPhi);
    NormMatrixByColumn(UMatrixMomentum);
    NormMatrixByColumn(UMatrixXVtxPosition);
    NormMatrixByColumn(UMatrixYVtxPosition);
    NormMatrixByColumn(UMatrixZVtxPosition);

    // Open output file
    TFile* OutputFile = new TFile("Selection_Histograms_Mod.root","RECREATE");

    // switch to output file
    OutputFile->cd();

    // TODO Find an other way to put the histograms in file!

    for(unsigned int histo_no = 0; histo_no < GenLabel.size(); histo_no++)
    {
        OutputFile->WriteObject(SelectionTrackRange.at(histo_no), SelectionTrackRange.at(histo_no)->GetName());
        OutputFile->WriteObject(SelectionCosTheta.at(histo_no), SelectionCosTheta.at(histo_no)->GetName());
        OutputFile->WriteObject(SelectionTheta.at(histo_no), SelectionTheta.at(histo_no)->GetName());
        OutputFile->WriteObject(SelectionPhi.at(histo_no), SelectionPhi.at(histo_no)->GetName());
        OutputFile->WriteObject(SelectionMomentum.at(histo_no), SelectionMomentum.at(histo_no)->GetName());
        OutputFile->WriteObject(SelectionTrackLength.at(histo_no), SelectionTrackLength.at(histo_no)->GetName());
        OutputFile->WriteObject(SelXVtxPosition.at(histo_no), SelXVtxPosition.at(histo_no)->GetName());
        OutputFile->WriteObject(SelYVtxPosition.at(histo_no), SelYVtxPosition.at(histo_no)->GetName());
        OutputFile->WriteObject(SelZVtxPosition.at(histo_no), SelZVtxPosition.at(histo_no)->GetName());
    }

    for(unsigned int file_no = 1; file_no < BgrTrackRange.size(); file_no++)
    {
        for(unsigned int histo_no = 1; histo_no < BgrLabel.size(); histo_no++)
        {
            OutputFile->WriteObject(BgrTrackRange.at(file_no).at(histo_no), BgrTrackRange.at(file_no).at(histo_no)->GetName());
            OutputFile->WriteObject(BgrCosTheta.at(file_no).at(histo_no), BgrCosTheta.at(file_no).at(histo_no)->GetName());
            OutputFile->WriteObject(BgrTheta.at(file_no).at(histo_no), BgrTheta.at(file_no).at(histo_no)->GetName());
            OutputFile->WriteObject(BgrPhi.at(file_no).at(histo_no), BgrPhi.at(file_no).at(histo_no)->GetName());
            OutputFile->WriteObject(BgrMomentum.at(file_no).at(histo_no), BgrMomentum.at(file_no).at(histo_no)->GetName());
            OutputFile->WriteObject(BgrTrackLength.at(file_no).at(histo_no), BgrTrackLength.at(file_no).at(histo_no)->GetName());
            OutputFile->WriteObject(BgrXVtxPosition.at(file_no).at(histo_no), BgrXVtxPosition.at(file_no).at(histo_no)->GetName());
            OutputFile->WriteObject(BgrYVtxPosition.at(file_no).at(histo_no), BgrYVtxPosition.at(file_no).at(histo_no)->GetName());
            OutputFile->WriteObject(BgrZVtxPosition.at(file_no).at(histo_no), BgrZVtxPosition.at(file_no).at(histo_no)->GetName());
        }

    }

    // Write efficiency (-ies)
    for(unsigned int eff_no = 0; eff_no < EffTrackRange.size(); eff_no++)
    {
        OutputFile->WriteObject(EffTrackRange.at(eff_no), ("EffTrackRange "+ std::to_string(eff_no)).c_str());
        OutputFile->WriteObject(EffCosTheta.at(eff_no), ("EffCosTheta "+ std::to_string(eff_no)).c_str());
        OutputFile->WriteObject(EffTheta.at(eff_no), ("EffTheta "+ std::to_string(eff_no)).c_str());
        OutputFile->WriteObject(EffPhi.at(eff_no), ("EffPhi "+ std::to_string(eff_no)).c_str());
        OutputFile->WriteObject(EffMomentum.at(eff_no), ("EffMomentum "+ std::to_string(eff_no)).c_str());
        OutputFile->WriteObject(EffXVtxPosition.at(eff_no), ("EffXVtxPosition "+ std::to_string(eff_no)).c_str());
        OutputFile->WriteObject(EffYVtxPosition.at(eff_no), ("EffYVtxPosition "+ std::to_string(eff_no)).c_str());
        OutputFile->WriteObject(EffZVtxPosition.at(eff_no), ("EffZVtxPosition "+ std::to_string(eff_no)).c_str());
    }

    // Write purities
    for(unsigned int pur_no = 0; pur_no < PurTrackRange.size(); pur_no++)
    {
        OutputFile->WriteObject(PurTrackRange.at(pur_no),("PurTrackRange "+ std::to_string(pur_no)).c_str());
        OutputFile->WriteObject(PurCosTheta.at(pur_no),("PurCosTheta "+ std::to_string(pur_no)).c_str());
        OutputFile->WriteObject(PurTheta.at(pur_no),("PurTheta "+ std::to_string(pur_no)).c_str());
        OutputFile->WriteObject(PurPhi.at(pur_no),("PurPhi "+ std::to_string(pur_no)).c_str());
        OutputFile->WriteObject(PurMomentum.at(pur_no),("PurMomentum "+ std::to_string(pur_no)).c_str());
        OutputFile->WriteObject(PurTrackLength.at(pur_no),("PurTrackLength "+ std::to_string(pur_no)).c_str());
        OutputFile->WriteObject(PurXVtxPosition.at(pur_no),("PurXVtxPosition "+ std::to_string(pur_no)).c_str());
        OutputFile->WriteObject(PurYVtxPosition.at(pur_no),("PurYVtxPosition "+ std::to_string(pur_no)).c_str());
        OutputFile->WriteObject(PurZVtxPosition.at(pur_no),("PurZVtxPosition "+ std::to_string(pur_no)).c_str());
    }

    OutputFile->WriteObject(UMatrixTrackRange, "UMatrixTrackRange");
    OutputFile->WriteObject(UMatrixCosTheta, "UMatrixCosTheta");
    OutputFile->WriteObject(UMatrixTheta, "UMatrixTheta");
    OutputFile->WriteObject(UMatrixPhi, "UMatrixPhi");
    OutputFile->WriteObject(UMatrixMomentum, "UMatrixMomentum");
    OutputFile->WriteObject(UMatrixXVtxPosition, "UMatrixXVtxPosition");
    OutputFile->WriteObject(UMatrixYVtxPosition, "UMatrixYVtxPosition");
    OutputFile->WriteObject(UMatrixZVtxPosition, "UMatrixZVtxPosition");

    OutputFile->Close();
}

std::vector<TGraph*> ReadFluxSystematics(const std::string &PathToFile)
{
    std::ifstream InputFile(PathToFile);
    
    // Readout variables
    double Energy;
    double nu_mu;
    double nu_mubar;
    double nu_e;
    double nu_ebar;
    std::string Line;
    
    // Initialize data vectors
    std::vector<double> EnergyData;
    std::vector<double> nu_muData;
    std::vector<double> nu_mubarData;
    std::vector<double> nu_eData;
    std::vector<double> nu_ebarData;
    
    // Skip header
    std::getline(InputFile,Line);
    
    // Read variables
    while(InputFile >> Energy >> nu_mu >> nu_mubar >> nu_e >> nu_ebar)
    {
        EnergyData.push_back(Energy);
        nu_muData.push_back(nu_mu);
        nu_mubarData.push_back(nu_mubar);
        nu_eData.push_back(nu_e);
        nu_ebarData.push_back(nu_ebar);
    }
    
    // Convert, because root sucks
    double* EnergyEntry = EnergyData.data();
    double* nu_muEntry =nu_muData.data();
    double* nu_mubarEntry = nu_mubarData.data();
    double* nu_eEntry = nu_eData.data();
    double* nu_ebarEntry =nu_ebarData.data();
    
    std::vector<TGraph*> Graphs;
    Graphs.push_back(new TGraph(EnergyData.size(),EnergyEntry,nu_muEntry));
    Graphs.push_back(new TGraph(EnergyData.size(),EnergyEntry,nu_mubarEntry));
    Graphs.push_back(new TGraph(EnergyData.size(),EnergyEntry,nu_eEntry));
    Graphs.push_back(new TGraph(EnergyData.size(),EnergyEntry,nu_ebarEntry));
    
    return Graphs;
}

float CalcRange(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2)
{
    return std::sqrt(pow(x_1-x_2, 2) + pow(y_1-y_2, 2) + pow(z_1-z_2, 2));
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
