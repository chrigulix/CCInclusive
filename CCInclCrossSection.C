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

// Function which calculates the distance between two points
float CalcRange(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2);

// Function which checks if a point is in the FV
bool inFV(double x, double y, double z);

// Add two histogramms with indices First and Last and weight
void AddHistograms(std::vector<TH1F*>& HistVector, unsigned int First, unsigned int Last, float Weight, bool EraseLast = false);

// Normalize Matrix by row
void NormMatrixByColumn(TH2F* UMatrix);

// Unsmearing of selected events
void SelectionUnsmearing(TH2F*& UMatrix, TH1F*& SVector);

// Momentum calculation
void MomentumSplinePreparation();

// Get Momentum
float GetMomentum(float TrackLength);

void CalcSigEfficiency(std::vector<TH1F*>& HistVector);

// CC inlcusive cross section function (main) 
void CCInclCrossSection()
{
    float NumberOfTargets = (FVx - 2*borderx) * (FVy - 2*bordery) * (FVz - 2*borderz) * Density * Avogadro/ArMass*NoNucleons;
    
    std::string Folder = "/home/christoph/anatrees/CCInclusiveNote";
//     std::string Folder = "/home/christoph/anatrees/ThesisSelection";
    
    // Output file file type
    std::string FileType = "pdf";
//     std::string FileType = "png";

    // Data input file vector
    std::vector<TChain*> ChainVec;

    // Histogram Vectors
    std::vector<TH1F*> SelectionTrackRange;
    std::vector<TH1F*> SelectionCosTheta;
    std::vector<TH1F*> SelectionTheta;
    std::vector<TH1F*> SelectionPhi;
    std::vector<TH1F*> SelectionMomentum;

    // Unsemaring Matrix
    TH2F* UMatrixTrackRange;
    TH2F* UMatrixCosTheta;
    TH2F* UMatrixTheta;
    TH2F* UMatrixPhi;
    TH2F* UMatrixMomentum;
    
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
    int NuPDGTruth[10];
    float TrueLeptonMomentum[10];
    float NuEnergyTruth[10];
    float nuvtxx_truth[10]; //true vertex x (in cm)
    float nuvtxy_truth[10];
    float nuvtxz_truth[10];

    int PDGTruth[5000];
    short TrkBestPlane[5000];
    short TrkOrigin[5000][3];

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
    
    // Check MA normalization
    unsigned int MASamples = 0;
    double MASamplesCorr = 0;
    
    TH1D* NuMuFlux;
//     TH1F* NuMuFlux = new TH1F("NuMuFlux","NuMuFlux",NumberOfBins,0,3);
    
    TFile* BNBFlux = new TFile("/home/christoph/anatrees/BNBFlux/numode_bnb_470m_r200.root");
    
    NuMuFlux = (TH1D*) BNBFlux->Get("numu");
    
    IntegratedFlux = NuMuFlux->Integral()*4.95e19/1e20;
    
//     std::cout << TempNuMuFlux->Integral(5,NuMuFlux->GetNbinsX())*4.95e19/1e20 << std::endl;
    std::cout << IntegratedFlux << std::endl;
    
//     TempNuMuFlux->Rebin(3);
    
//     for(unsigned int bin_no = 1; bin_no <= NuMuFlux->GetNbinsX(); bin_no++)
//     {
//         NuMuFlux->SetBinContent(bin_no, TempNuMuFlux->GetBinContent(bin_no));
//         NuMuFlux->SetBinError(bin_no, TempNuMuFlux->GetBinError(bin_no));
//     }
    
//     NuMuFlux->Scale(DataPOT/1e20);
    
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
    GenLabel.push_back("Data On-Beam BNB");
    ScalingFactors.push_back(1);

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_1_Mod.root").c_str());
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_2_Mod.root").c_str());
    GenLabel.push_back("Data Off-Beam BNBEXT");
    ScalingFactors.push_back(1.2300);

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_prodgenie_bnb_nu_cosmic_uboone_v05_08_00_Mod.root").c_str());
//     ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_prodgenie_bnb_nu_cosmic_uboone_field_v05_08_00_Mod.root").c_str());
    GenLabel.push_back("MC Selection");
    ScalingFactors.push_back(DataPOT/MCPOT);
    
    // MC selection categories
    GenLabel.push_back("MC Cosmics");
    ScalingFactors.push_back(DataPOT/MCPOT);
    GenLabel.push_back("MC Beam Backgrounds");
    ScalingFactors.push_back(DataPOT/MCPOT);
    GenLabel.push_back("MC True Selection");
    ScalingFactors.push_back(DataPOT/MCPOT);
    
    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_MA_v05_08_00_Mod.root").c_str());
    GenLabel.push_back("MA Adjusted Selection");
    ScalingFactors.push_back(DataPOT/MCPOT*6920/6194);
    
    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_TEM_v05_08_00_Mod.root").c_str());
    GenLabel.push_back("TEM Selection");
    ScalingFactors.push_back(DataPOT/MCPOT);
    
    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_MEC_v05_08_00_Mod.root").c_str());
    GenLabel.push_back("MEC Selection");
    ScalingFactors.push_back(DataPOT/MCPOT);
    
    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((Folder+"/Hist_MC_Truth_prodgenie_bnb_nu_cosmic_uboone_v05_08_00.root").c_str());
    GenLabel.push_back("MC Truth");
    ScalingFactors.push_back(DataPOT/MCPOT);
//     ScalingFactors.push_back(DataPOT/TruthPOT);
    
    // Efficiency histograms
    GenLabel.push_back("Efficiency");
    ScalingFactors.push_back(1);

    // Loop over all generation labels
    for(const auto& Label : GenLabel)
    {
        SelectionTrackRange.push_back(new TH1F(("Track Range"+Label).c_str(),"Muon Track Range",NumberOfBins,0,1036.8));
        SelectionTrackRange.back() -> SetStats(0);
        SelectionTrackRange.back() -> GetXaxis() -> SetTitle("Muon track range [cm]");
        SelectionTrackRange.back() -> GetYaxis() -> SetTitle("No. of events");
//         SelectionTrackRange.back() -> GetYaxis() -> SetTitle("d#sigma/dl [cm^{2}/cm]");

        SelectionCosTheta.push_back(new TH1F(("cos#theta-Angle"+Label).c_str(),"Cosine of #theta-Angle",NumberOfBins,-1,1));
        SelectionCosTheta.back() -> SetStats(0);
        SelectionCosTheta.back() -> GetXaxis() -> SetTitle("Muon cos(#theta)");
        SelectionCosTheta.back() -> GetYaxis() -> SetTitle("No. of events");
//         SelectionCosTheta.back() -> GetYaxis() -> SetTitle("d#sigma/d(cos#theta) [cm^{2}/cos(#theta)]");
        
        SelectionTheta.push_back(new TH1F(("#theta-Angle"+Label).c_str(),"#theta-Angle",NumberOfBins,0,180));
        SelectionTheta.back() -> SetStats(0);
        SelectionTheta.back() -> GetXaxis() -> SetTitle("Muon #theta-Angle [#circ]");
        SelectionTheta.back() -> GetYaxis() -> SetTitle("No. of events");
//         SelectionTheta.back() -> GetYaxis() -> SetTitle("d#sigma/d#theta [cm^{2}/rad]");

        SelectionPhi.push_back(new TH1F(("#phi-Angle"+Label).c_str(),"#varphi-Angle",NumberOfBins,-180,180));
        SelectionPhi.back() -> SetStats(0);
        SelectionPhi.back() -> GetXaxis() -> SetTitle("Muon #varphi-Angle [#circ]");
        SelectionPhi.back() -> GetYaxis() -> SetTitle("No. of events");
//         SelectionPhi.back() -> GetYaxis() -> SetTitle("d#sigma/d#phi [cm^{2}/#circ]");

        SelectionMomentum.push_back(new TH1F(("Momentum"+Label).c_str(),"Muon Momentum",NumberOfBins,0,3));
        SelectionMomentum.back() -> SetStats(0);
        SelectionMomentum.back() -> GetXaxis() -> SetTitle("Muon Momentum p_{#mu} [GeV/c]");
        SelectionMomentum.back() -> GetYaxis() -> SetTitle("No. of events");
//         SelectionMomentum.back() -> GetYaxis() -> SetTitle("d#sigma/dp [cm^{2}/(GeV/c)]");
    } // loop over generation label
    
    // Initialize smearing matrices
    UMatrixTrackRange = new TH2F("Unsmering Matrix Track Range","Smearing Matrix Track Range",NumberOfBins,0,1036.8,NumberOfBins,0,1036.8);
    UMatrixTrackRange -> GetXaxis() -> SetTitle("Muon track length (truth) [cm]");
    UMatrixTrackRange -> GetYaxis() -> SetTitle("Muon track length (reco) [cm]");
    UMatrixCosTheta = new TH2F("Smearing Matrix CosTheta","Smearing Matrix cos(#theta)",NumberOfBins,0,-1,NumberOfBins,0,-1);
    UMatrixCosTheta -> GetXaxis() -> SetTitle("Muon cos(#theta) (truth)");
    UMatrixCosTheta -> GetYaxis() -> SetTitle("Muon cos(#theta) (reco)");
    UMatrixTheta = new TH2F("Smearing Matrix Theta","Smearing Matrix Theta",NumberOfBins,0,180,NumberOfBins,0,180);
    UMatrixTheta -> GetXaxis() -> SetTitle("Muon #theta-Angle (truth) [#circ]");
    UMatrixTheta -> GetYaxis() -> SetTitle("Muon #theta-Angle (reco) [#circ]");
    UMatrixPhi = new TH2F("Smearing Matrix Phi","Smearing Matrix Phi",NumberOfBins,0,-180,NumberOfBins,0,180);
    UMatrixPhi -> GetXaxis() -> SetTitle("Muon #varphi-Angle (truth) [#circ]");
    UMatrixPhi -> GetYaxis() -> SetTitle("Muon #varphi-Angle (reco) [#circ]");
    UMatrixMomentum = new TH2F("Smearing Matrix Momentum","Smearing Matrix Momentum",NumberOfBins,0,3,NumberOfBins,0,3);
    UMatrixMomentum -> GetXaxis() -> SetTitle("Muon Momentum (truth) [#circ]");
    UMatrixMomentum -> GetYaxis() -> SetTitle("Muon Momentum (reco) [#circ]");

    // Loop over all files
    for(unsigned int file_no = 0; file_no < ChainVec.size(); file_no++)
    {
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
        if(file_no == 2 || file_no == 6)
        {
            ChainVec.at(file_no) -> SetBranchAddress("MCTrackCand", &MCTrkID);
            ChainVec.at(file_no) -> SetBranchAddress("MCVertexCand", &MCVtxID);
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
            ChainVec.at(file_no) -> SetBranchAddress("StartPointx", XMCTrackStart);
            ChainVec.at(file_no) -> SetBranchAddress("StartPointy", YMCTrackStart);
            ChainVec.at(file_no) -> SetBranchAddress("StartPointz", ZMCTrackStart);
            ChainVec.at(file_no) -> SetBranchAddress("EndPointx", XMCTrackEnd);
            ChainVec.at(file_no) -> SetBranchAddress("EndPointy", YMCTrackEnd);
            ChainVec.at(file_no) -> SetBranchAddress("EndPointz", ZMCTrackEnd);
            ChainVec.at(file_no) -> SetBranchAddress("theta", MCTheta);
            ChainVec.at(file_no) -> SetBranchAddress("enu_truth", NuEnergyTruth);
            ChainVec.at(file_no) -> SetBranchAddress("phi", MCPhi);
            ChainVec.at(file_no) -> SetBranchAddress("Eng", MCEnergy);
            
//             if(file_no == 2 || file_no == 6)
//             {
//                 ChainVec.at(file_no) -> SetBranchAddress("MCVertexCand", &MCVtxID);
//             }
            
            // Backtracker information only if there are reco objects
            if(file_no < 6)
            {
                ChainVec.at(file_no) -> SetBranchAddress("trkorigin_pandoraNu", TrkOrigin);
                ChainVec.at(file_no) -> SetBranchAddress("trkpidbestplane_pandoraNu", TrkBestPlane);
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
        
        unsigned int cosmics = 0;
        unsigned int beambgr = 0;
        
        // Loop over all events
        for(unsigned int tree_index = 0; tree_index < ChainVec.at(file_no) -> GetEntries(); tree_index++)
        {
            // Progress indicator
            if(!(tree_index % 1000)) std::cout << "Event\t" << tree_index << "\t of \t" << ChainVec.at(file_no) -> GetEntries() << std::endl;

            // Skip corrupted events in truth file
            if(file_no == 6 && (tree_index == 11602 || tree_index == 11675 || tree_index == 13510 || tree_index == 33027 || tree_index == 33070 || tree_index == 36239 || tree_index == 44078)) continue;

            // Get tree entry for this event
            ChainVec.at(file_no) -> GetEntry(tree_index);

            // if there are reco products
            if(file_no < 3)
            {
                // Fill histograms as usual for all files
                SelectionTrackRange.at(file_no) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                SelectionCosTheta.at(file_no) -> Fill(std::cos(TrackTheta[TrkID]));
                SelectionTheta.at(file_no) -> Fill(TrackTheta[TrkID]/Pi*180);
                SelectionPhi.at(file_no) -> Fill(TrackPhi[TrkID]/Pi*180);
                SelectionMomentum.at(file_no) -> Fill(GetMomentum(TrackLength[TrkID]));
            }

            // if we are looking at the MC selection file
            if(file_no == 2)
            {                
                // If event is cosmic background
                if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] != 1)
                {
//                     if(file_no == 2)
//                     {
                        cosmics++;
                        // Fill cosmic background histograms
                        SelectionTrackRange.at(file_no+1) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                        SelectionCosTheta.at(file_no+1) -> Fill(std::cos(TrackTheta[TrkID]));
                        SelectionTheta.at(file_no+1) -> Fill(TrackTheta[TrkID]/Pi*180);
                        SelectionPhi.at(file_no+1) -> Fill(TrackPhi[TrkID]/Pi*180);
                        SelectionMomentum.at(file_no+1) -> Fill(GetMomentum(TrackLength[TrkID]));
//                     }
                }
                // else if event is other background
                else if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1 && ( nuPDGTruth[MCVtxID] != 14 || CCNCFlag[MCVtxID] == 1 || !inFV(XnuVtxTruth[MCVtxID],YnuVtxTruth[MCVtxID],ZnuVtxTruth[MCVtxID]) )/* || PDGTruth[MCTrkID] != 13*/)
                {
//                     if(file_no == 2)
//                     {
                        beambgr++;
                        // Fill beam related background histograms
                        SelectionTrackRange.at(file_no+2) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                        SelectionCosTheta.at(file_no+2) -> Fill(std::cos(TrackTheta[TrkID]));
                        SelectionTheta.at(file_no+2) -> Fill(TrackTheta[TrkID]/Pi*180);
                        SelectionPhi.at(file_no+2) -> Fill(TrackPhi[TrkID]/Pi*180);
                        SelectionMomentum.at(file_no+2) -> Fill(GetMomentum(TrackLength[TrkID]));
//                     }
                }
                else // if event is signal and truth
                {
                    // Fill background histograms
                    SelectionTrackRange.at(file_no+3) -> Fill(CalcRange(XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID]));
                    SelectionCosTheta.at(file_no+3) -> Fill(std::cos(MCTheta[MCTrkID]));
                    SelectionTheta.at(file_no+3) -> Fill(MCTheta[MCTrkID]/Pi*180);
                    SelectionPhi.at(file_no+3) -> Fill(MCPhi[MCTrkID]/Pi*180);
                    SelectionMomentum.at(file_no+3) -> Fill(TrueLeptonMomentum[MCVtxID]);
                    
//                     if(file_no == 2)
//                     {
                        // Fill smearing matrix
                        UMatrixTrackRange -> Fill( CalcRange(XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID]),CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]) );
                        UMatrixCosTheta -> Fill( std::cos(MCTheta[MCTrkID]),std::cos(TrackTheta[TrkID]) );
                        UMatrixTheta -> Fill( MCTheta[MCTrkID]/Pi*180,TrackTheta[TrkID]/Pi*180 );
                        UMatrixPhi -> Fill( MCPhi[MCTrkID]/Pi*180,TrackPhi[TrkID]/Pi*180 );
                        UMatrixMomentum -> Fill( TrueLeptonMomentum[MCVtxID],GetMomentum(TrackLength[TrkID]) );
//                     }
                    SelectedEvents++;
                }
            }
//             else if(file_no == 3)
//             {
//                 SelectionTrackRange.at(file_no+3) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),HistogramWeight);
//                 SelectionCosTheta.at(file_no+3) -> Fill(std::cos(TrackTheta[TrkID]),HistogramWeight);
//                 SelectionTheta.at(file_no+3) -> Fill(TrackTheta[TrkID]/Pi*180,HistogramWeight);
//                 SelectionPhi.at(file_no+3) -> Fill(TrackPhi[TrkID]/Pi*180,HistogramWeight);
//                 SelectionMomentum.at(file_no+3) -> Fill(GetMomentum(TrackLength[TrkID]),HistogramWeight);
//                 
//                 if(HistogramWeight!=1.0) std::cout << HistogramWeight << std:: endl;
//                 MASamples++;
//                 MASamplesCorr += HistogramWeight;
//             }
            else if(file_no > 2 && file_no < 6)
            {
                SelectionTrackRange.at(file_no+3) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]),HistogramWeight);
                SelectionCosTheta.at(file_no+3) -> Fill(std::cos(TrackTheta[TrkID]),HistogramWeight);
                SelectionTheta.at(file_no+3) -> Fill(TrackTheta[TrkID]/Pi*180,HistogramWeight);
                SelectionPhi.at(file_no+3) -> Fill(TrackPhi[TrkID]/Pi*180,HistogramWeight);
                SelectionMomentum.at(file_no+3) -> Fill(GetMomentum(TrackLength[TrkID]),HistogramWeight);
                
//                 if(HistogramWeight!=1.0) std::cout << HistogramWeight << std:: endl;
                
                if(file_no == 3)
                {
                    MASamples++;
                    MASamplesCorr += HistogramWeight;
                }
            }
            // if truth selection file
            else if(file_no == 6 && MCTrkID >= 0 && nuPDGTruth[MCVtxID] == 14 && inFV(XnuVtxTruth[MCVtxID],YnuVtxTruth[MCVtxID],ZnuVtxTruth[MCVtxID]))
            {
                // Fill background histograms
                SelectionTrackRange.at(file_no+3) -> Fill(CalcRange(XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID]));
                SelectionCosTheta.at(file_no+3) -> Fill(std::cos(MCTheta[MCTrkID]));
                SelectionTheta.at(file_no+3) -> Fill(MCTheta[MCTrkID]/Pi*180);
                SelectionPhi.at(file_no+3) -> Fill(MCPhi[MCTrkID]/Pi*180);
                SelectionMomentum.at(file_no+3) -> Fill(TrueLeptonMomentum[MCVtxID]);
                
                ExpectedEvents++;
            }
            
        } // Event loop

        // Reset branch addresses to avoid problems
        ChainVec.at(file_no) -> ResetBranchAddresses();
        
        // Delete chain from memory (preventing overflow)
        ChainVec.at(file_no) -> Delete();
        
        std::cout << "Cosmic count: " << cosmics << std::endl;
        std::cout << "Beam bgr count: " << beambgr << std::endl;

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
    
    // loop over all histograms
    for(unsigned int hist_no = 0; hist_no < GenLabel.size(); hist_no++)
    {
        // Calculate standard deviation for all histograms
        SelectionTrackRange.at(hist_no)->Sumw2();
        SelectionCosTheta.at(hist_no)->Sumw2();
        SelectionTheta.at(hist_no)->Sumw2();
        SelectionPhi.at(hist_no)->Sumw2();
        SelectionMomentum.at(hist_no)->Sumw2();
    } // histogram loop

    // loop over scaling factors 
    for(unsigned int scale_no = 0; scale_no < ScalingFactors.size(); scale_no++)
    {
        // Scale histograms
        SelectionTrackRange.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionCosTheta.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionTheta.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionPhi.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionMomentum.at(scale_no)->Scale(ScalingFactors.at(scale_no));
    } // scaling loop
    
    // Fill efficiency
    SelectionTrackRange.back()->Divide(SelectionTrackRange.at(5),SelectionTrackRange.at(9));
    SelectionCosTheta.back()->Divide(SelectionCosTheta.at(5),SelectionCosTheta.at(9));
    SelectionTheta.back()->Divide(SelectionTheta.at(5),SelectionTheta.at(9));
    SelectionPhi.back()->Divide(SelectionPhi.at(5),SelectionPhi.at(9));
    SelectionMomentum.back()->Divide(SelectionMomentum.at(5),SelectionMomentum.at(9));
    
    // Subtract offbeam from onbeam data, overwrite onbeam
    AddHistograms(SelectionTrackRange,0,1,-1);
    AddHistograms(SelectionCosTheta,0,1,-1);
    AddHistograms(SelectionTheta,0,1,-1);
    AddHistograms(SelectionPhi,0,1,-1);
    AddHistograms(SelectionMomentum,0,1,-1);
    
    // Normalize matrices by row
    NormMatrixByColumn(UMatrixTrackRange);
    NormMatrixByColumn(UMatrixCosTheta);
    NormMatrixByColumn(UMatrixTheta);
    NormMatrixByColumn(UMatrixPhi);
    NormMatrixByColumn(UMatrixMomentum);

    // Subtract cosmic background from data and overwrite data histogram
    AddHistograms(SelectionTrackRange,0,3,-1);
    AddHistograms(SelectionCosTheta,0,3,-1);
    AddHistograms(SelectionTheta,0,3,-1);
    AddHistograms(SelectionPhi,0,3,-1);
    AddHistograms(SelectionMomentum,0,3,-1);
    // Subtract beam background from data and overwrite data histogram
    AddHistograms(SelectionTrackRange,0,4,-1);
    AddHistograms(SelectionPhi,0,4,-1);
    AddHistograms(SelectionCosTheta,0,4,-1);
    AddHistograms(SelectionTheta,0,4,-1);
    AddHistograms(SelectionMomentum,0,4,-1);
    
    // Subtract cosmic background from mc selection
    AddHistograms(SelectionTrackRange,2,3,-1);
    AddHistograms(SelectionPhi,2,3,-1);
    AddHistograms(SelectionCosTheta,2,3,-1);
    AddHistograms(SelectionTheta,2,3,-1);
    AddHistograms(SelectionMomentum,2,3,-1);
    
    // Subtract beam backgrounds from mc selection
    AddHistograms(SelectionTrackRange,2,4,-1);
    AddHistograms(SelectionPhi,2,4,-1);
    AddHistograms(SelectionCosTheta,2,4,-1);
    AddHistograms(SelectionTheta,2,4,-1);
    AddHistograms(SelectionMomentum,2,4,-1);
    
    // Subtract cosmic background from mc selection
    AddHistograms(SelectionTrackRange,6,3,-1);
    AddHistograms(SelectionPhi,6,3,-1);
    AddHistograms(SelectionCosTheta,6,3,-1);
    AddHistograms(SelectionTheta,6,3,-1);
    AddHistograms(SelectionMomentum,6,3,-1);
    
    // Subtract beam backgrounds from mc selection
    AddHistograms(SelectionTrackRange,6,4,-1);
    AddHistograms(SelectionPhi,6,4,-1);
    AddHistograms(SelectionCosTheta,6,4,-1);
    AddHistograms(SelectionTheta,6,4,-1);
    AddHistograms(SelectionMomentum,6,4,-1);
    
    // Subtract cosmic background from mc selection
    AddHistograms(SelectionTrackRange,7,3,-1);
    AddHistograms(SelectionPhi,7,3,-1);
    AddHistograms(SelectionCosTheta,7,3,-1);
    AddHistograms(SelectionTheta,7,3,-1);
    AddHistograms(SelectionMomentum,7,3,-1);
    
    // Subtract beam backgrounds from mc selection
    AddHistograms(SelectionTrackRange,7,4,-1);
    AddHistograms(SelectionPhi,7,4,-1);
    AddHistograms(SelectionCosTheta,7,4,-1);
    AddHistograms(SelectionTheta,7,4,-1);
    AddHistograms(SelectionMomentum,7,4,-1);
    
    // Subtract cosmic background from mc selection
    AddHistograms(SelectionTrackRange,8,3,-1);
    AddHistograms(SelectionPhi,8,3,-1);
    AddHistograms(SelectionCosTheta,8,3,-1);
    AddHistograms(SelectionTheta,8,3,-1);
    AddHistograms(SelectionMomentum,8,3,-1);
    
    // Subtract beam backgrounds from mc selection
    AddHistograms(SelectionTrackRange,8,4,-1);
    AddHistograms(SelectionPhi,8,4,-1);
    AddHistograms(SelectionCosTheta,8,4,-1);
    AddHistograms(SelectionTheta,8,4,-1);
    AddHistograms(SelectionMomentum,8,4,-1);
    
//     CalcSigEfficiency(SelectionTrackRange);
//     CalcSigEfficiency(SelectionPhi);
//     CalcSigEfficiency(SelectionCosTheta);
//     CalcSigEfficiency(SelectionTheta);
//     CalcSigEfficiency(SelectionMomentum);
    
    // Draw histogram
    TCanvas *Canvas1a = new TCanvas("Range a", "Range", 1400, 1000);
    Canvas1a->cd();
    SelectionTrackRange.at(2)->SetMaximum(1.2*SelectionTrackRange.at(0)->GetBinContent(SelectionTrackRange.at(0)->GetMaximumBin()));
    SelectionTrackRange.at(2)->SetMinimum(0.0);
    SelectionTrackRange.at(2)->SetFillColor(46);
    SelectionTrackRange.at(2)->Draw("E2 SAME");
    SelectionTrackRange.at(6)->SetFillColor(28);
    SelectionTrackRange.at(6)->Draw("E2 SAME");
    SelectionTrackRange.at(7)->SetFillColor(30);
    SelectionTrackRange.at(7)->Draw("E2 SAME");
    SelectionTrackRange.at(8)->SetFillColor(38);
    SelectionTrackRange.at(8)->Draw("E2 SAME");
    SelectionTrackRange.at(0)->SetLineWidth(2);
    SelectionTrackRange.at(0)->SetLineColor(1);
    SelectionTrackRange.at(0)->SetMarkerColor(1);
    SelectionTrackRange.at(0)->Draw("SAME");
    Canvas1a->SaveAs(("ScaledOn-OffBeamSelRange."+FileType).c_str());
    
    TCanvas *Canvas2a = new TCanvas("CosTheta a", "CosTheta", 1400, 1000);
    Canvas2a->cd();
    SelectionCosTheta.at(2)->SetMaximum(1.2*SelectionCosTheta.at(0)->GetBinContent(SelectionCosTheta.at(0)->GetMaximumBin()));
    SelectionCosTheta.at(2)->SetMinimum(0.0);
    SelectionCosTheta.at(2)->SetFillColor(46);
    SelectionCosTheta.at(2)->Draw("E2");
    SelectionCosTheta.at(6)->SetFillColor(28);
    SelectionCosTheta.at(6)->Draw("E2 SAME");
    SelectionCosTheta.at(7)->SetFillColor(30);
    SelectionCosTheta.at(7)->Draw("E2 SAME");
    SelectionCosTheta.at(8)->SetFillColor(38);
    SelectionCosTheta.at(8)->Draw("E2 SAME");
    SelectionCosTheta.at(0)->SetLineWidth(2);
    SelectionCosTheta.at(0)->SetLineColor(1);
    SelectionCosTheta.at(0)->SetMarkerColor(1);
    SelectionCosTheta.at(0)->Draw("SAME");
    Canvas2a->SaveAs(("ScaledOn-OffBeamSelCosTheta."+FileType).c_str());
    
    TCanvas *Canvas3a = new TCanvas("Theta a", "Theta", 1400, 1000);
    Canvas3a->cd();
    SelectionTheta.at(2)->SetMaximum(1.2*SelectionTheta.at(0)->GetBinContent(SelectionTheta.at(0)->GetMaximumBin()));
    SelectionTheta.at(2)->SetMinimum(0.0);
    SelectionTheta.at(2)->SetFillColor(46);
    SelectionTheta.at(2)->Draw("E2");
    SelectionTheta.at(6)->SetFillColor(28);
    SelectionTheta.at(6)->Draw("E2 SAME");
    SelectionTheta.at(7)->SetFillColor(30);
    SelectionTheta.at(7)->Draw("E2 SAME");
    SelectionTheta.at(8)->SetFillColor(38);
    SelectionTheta.at(8)->Draw("E2 SAME");
    SelectionTheta.at(0)->SetLineWidth(2);
    SelectionTheta.at(0)->SetLineColor(1);
    SelectionTheta.at(0)->SetMarkerColor(1);
    SelectionTheta.at(0)->Draw("SAME");
    Canvas3a->SaveAs(("ScaledOn-OffBeamSelTheta."+FileType).c_str());
    
    TCanvas *Canvas4a = new TCanvas("Phi a", "Phi", 1400, 1000);
    Canvas4a->cd();
    SelectionPhi.at(2)->SetMaximum(1.2*SelectionPhi.at(0)->GetBinContent(SelectionPhi.at(0)->GetMaximumBin()));
    SelectionPhi.at(2)->SetMinimum(0.0);
    SelectionPhi.at(2)->SetFillColor(46);
    SelectionPhi.at(2)->Draw("E2");
    SelectionPhi.at(6)->SetFillColor(28);
    SelectionPhi.at(6)->Draw("E2 SAME");
    SelectionPhi.at(7)->SetFillColor(30);
    SelectionPhi.at(7)->Draw("E2 SAME");
    SelectionPhi.at(8)->SetFillColor(38);
    SelectionPhi.at(8)->Draw("E2 SAME");
    SelectionPhi.at(0)->SetLineWidth(2);
    SelectionPhi.at(0)->SetLineColor(1);
    SelectionPhi.at(0)->SetMarkerColor(1);
    SelectionPhi.at(0)->Draw("SAME");
    Canvas4a->SaveAs(("ScaledOn-OffBeamSelPhi."+FileType).c_str());
    
    TCanvas *Canvas5a = new TCanvas("Momentum a", "Momentum", 1400, 1000);
    Canvas5a->cd();
    SelectionMomentum.at(2)->SetMaximum(1.5*SelectionMomentum.at(0)->GetBinContent(SelectionMomentum.at(0)->GetMaximumBin()));
    SelectionMomentum.at(2)->SetMinimum(0.0);
    SelectionMomentum.at(2)->SetFillColor(46);
    SelectionMomentum.at(2)->Draw("E2");
    SelectionMomentum.at(6)->SetFillColor(28);
    SelectionMomentum.at(6)->Draw("E2 SAME");
    SelectionMomentum.at(7)->SetFillColor(30);
    SelectionMomentum.at(7)->Draw("E2 SAME");
    SelectionMomentum.at(8)->SetFillColor(38);
    SelectionMomentum.at(8)->Draw("E2 SAME");
    SelectionMomentum.at(0)->SetLineWidth(2);
    SelectionMomentum.at(0)->SetLineColor(1);
    SelectionMomentum.at(0)->SetMarkerColor(1);
    SelectionMomentum.at(0)->Draw("SAME");
    Canvas5a->SaveAs(("ScaledOn-OffBeamSelMomentum."+FileType).c_str());
    
    // Selection data/MC loop
    for(unsigned int hist_no = 0; hist_no < 3; hist_no++)
    {
        // Unsmearing of data
        SelectionUnsmearing(UMatrixTrackRange,SelectionTrackRange.at(hist_no));
        SelectionUnsmearing(UMatrixCosTheta,SelectionCosTheta.at(hist_no));
        SelectionUnsmearing(UMatrixTheta,SelectionTheta.at(hist_no));
        SelectionUnsmearing(UMatrixPhi,SelectionPhi.at(hist_no));
        SelectionUnsmearing(UMatrixMomentum,SelectionMomentum.at(hist_no));
        
        // Efficiency unfolding
        SelectionTrackRange.at(hist_no)->Divide(SelectionTrackRange.back());
        SelectionCosTheta.at(hist_no)->Divide(SelectionCosTheta.back());
        SelectionTheta.at(hist_no)->Divide(SelectionTheta.back());
        SelectionPhi.at(hist_no)->Divide(SelectionPhi.back());
        SelectionMomentum.at(hist_no)->Divide(SelectionMomentum.back());
    }
    
    // Efficiency unfolding
    SelectionTrackRange.at(0)->Multiply(SelectionTrackRange.at(3));
    SelectionCosTheta.at(0)->Multiply(SelectionCosTheta.at(3));
    SelectionTheta.at(0)->Multiply(SelectionTheta.at(3));
    SelectionPhi.at(0)->Multiply(SelectionPhi.at(3));
    SelectionMomentum.at(0)->Multiply(SelectionMomentum.at(3));
    
    // Draw histogram
    TCanvas *Canvas1b = new TCanvas("Range b", "Range", 1400, 1000);
    Canvas1b->cd();
    SelectionTrackRange.at(2)->SetMaximum(1.2*SelectionTrackRange.at(0)->GetBinContent(SelectionTrackRange.at(0)->GetMaximumBin()));
    SelectionTrackRange.at(2)->SetMinimum(0.0);
    SelectionTrackRange.at(2)->SetFillColor(46);
    SelectionTrackRange.at(2)->Draw("E2");
    SelectionTrackRange.at(0)->SetLineWidth(2);
    SelectionTrackRange.at(0)->SetLineColor(1);
    SelectionTrackRange.at(0)->SetMarkerColor(1);
    SelectionTrackRange.at(0)->Draw("SAME");
    Canvas1b->SaveAs(("UnsmearedNoBGRRange."+FileType).c_str());
    
    TCanvas *Canvas2b = new TCanvas("CosTheta b", "CosTheta", 1400, 1000);
    Canvas2b->cd();
    SelectionCosTheta.at(2)->SetMaximum(1.2*SelectionCosTheta.at(0)->GetBinContent(SelectionCosTheta.at(0)->GetMaximumBin()));
    SelectionCosTheta.at(2)->SetMinimum(0.0);
    SelectionCosTheta.at(2)->SetFillColor(46);
    SelectionCosTheta.at(2)->Draw("E2");
    SelectionCosTheta.at(0)->SetLineWidth(2);
    SelectionCosTheta.at(0)->SetLineColor(1);
    SelectionCosTheta.at(0)->SetMarkerColor(1);
    SelectionCosTheta.at(0)->Draw("SAME");
    Canvas2b->SaveAs(("UnsmearedNoBGRCosTheta."+FileType).c_str());
    
    TCanvas *Canvas3b = new TCanvas("Theta b", "Theta", 1400, 1000);
    Canvas3b->cd();
    SelectionTheta.at(2)->SetMaximum(1.2*SelectionTheta.at(0)->GetBinContent(SelectionTheta.at(0)->GetMaximumBin()));
    SelectionTheta.at(2)->SetMinimum(0.0);
    SelectionTheta.at(2)->SetFillColor(46);
    SelectionTheta.at(2)->Draw("E2");
    SelectionTheta.at(0)->SetLineWidth(2);
    SelectionTheta.at(0)->SetLineColor(1);
    SelectionTheta.at(0)->SetMarkerColor(1);
    SelectionTheta.at(0)->Draw("SAME");
    Canvas3b->SaveAs(("UnsmearedNoBGRTheta."+FileType).c_str());
    
    TCanvas *Canvas4b = new TCanvas("Phi b", "Phi", 1400, 1000);
    Canvas4b->cd();
    SelectionPhi.at(2)->SetMaximum(1.2*SelectionPhi.at(0)->GetBinContent(SelectionPhi.at(0)->GetMaximumBin()));
    SelectionPhi.at(2)->SetMinimum(0.0);
    SelectionPhi.at(2)->SetFillColor(46);
    SelectionPhi.at(2)->Draw("E2");
    SelectionPhi.at(0)->SetLineWidth(2);
    SelectionPhi.at(0)->SetLineColor(1);
    SelectionPhi.at(0)->SetMarkerColor(1);
    SelectionPhi.at(0)->Draw("SAME");
    Canvas4b->SaveAs(("UnsmearedNoBGRSelPhi."+FileType).c_str());
    
    TCanvas *Canvas5b = new TCanvas("Momentum b", "Momentum", 1400, 1000);
    Canvas5b->cd();
    SelectionMomentum.at(2)->SetMaximum(1.5*SelectionMomentum.at(0)->GetBinContent(SelectionMomentum.at(0)->GetMaximumBin()));
    SelectionMomentum.at(2)->SetMinimum(0.0);
    SelectionMomentum.at(2)->SetFillColor(46);
    SelectionMomentum.at(2)->Draw("E2");
    SelectionMomentum.at(0)->SetLineWidth(2);
    SelectionMomentum.at(0)->SetLineColor(1);
    SelectionMomentum.at(0)->SetMarkerColor(1);
    SelectionMomentum.at(0)->Draw("SAME");
    Canvas5b->SaveAs(("UnsmearedNoBGRMomentum."+FileType).c_str());
    
    for(unsigned int hist_no = 0; hist_no < 3; hist_no++)
    {
        // Scaling to flux number of target nucleons and bin width
        SelectionTrackRange.at(hist_no)->Scale(1/NumberOfTargets/SelectionTrackRange.at(hist_no)->GetBinWidth(1)/IntegratedFlux);
        SelectionCosTheta.at(hist_no)->Scale(1/NumberOfTargets/SelectionCosTheta.at(hist_no)->GetBinWidth(1)/IntegratedFlux);
        SelectionTheta.at(hist_no)->Scale(1/NumberOfTargets/SelectionTheta.at(hist_no)->GetBinWidth(1)/IntegratedFlux);
        SelectionPhi.at(hist_no)->Scale(1/NumberOfTargets/SelectionPhi.at(hist_no)->GetBinWidth(1)/IntegratedFlux);
        SelectionMomentum.at(hist_no)->Scale(1/NumberOfTargets/SelectionMomentum.at(hist_no)->GetBinWidth(1)/IntegratedFlux);
    }
    
    TCanvas *Canvas0 = new TCanvas("Test", "Test", 1400, 1000);
    Canvas0->cd();
    SelectionPhi.at(3)->Draw();
    
    TCanvas *Canvas0a = new TCanvas("Testa", "Testa", 1400, 1000);
    Canvas0a->cd();
    SelectionPhi.at(4)->Draw();

    // Draw histogram
    TCanvas *Canvas1 = new TCanvas("Range", "Range", 1400, 1000);
    Canvas1->cd();
    SelectionTrackRange.at(2)->SetMaximum(1.2*SelectionTrackRange.at(0)->GetBinContent(SelectionTrackRange.at(0)->GetMaximumBin()));
    SelectionTrackRange.at(2)->SetMinimum(0.0);
    SelectionTrackRange.at(2)->SetFillColor(46);
    SelectionTrackRange.at(2) -> GetYaxis() -> SetTitle("d#sigma/dl [cm^{2}/cm/Nucleon]");
    SelectionTrackRange.at(2)->Draw("E2");
    SelectionTrackRange.at(0)->SetLineWidth(2);
    SelectionTrackRange.at(0)->SetLineColor(1);
    SelectionTrackRange.at(0)->SetMarkerColor(1);
    SelectionTrackRange.at(0)->Draw("SAME");
    Canvas1->SaveAs(("DiffCrossSectionRange."+FileType).c_str());
    
    TCanvas *Canvas2 = new TCanvas("CosTheta", "CosTheta", 1400, 1000);
    Canvas2->cd();
    SelectionCosTheta.at(2)->SetMaximum(1.2*SelectionCosTheta.at(0)->GetBinContent(SelectionCosTheta.at(0)->GetMaximumBin()));
    SelectionCosTheta.at(2)->SetMinimum(0.0);
    SelectionCosTheta.at(2)->SetFillColor(46);
    SelectionCosTheta.at(2) -> GetYaxis() -> SetTitle("d#sigma/d(cos#theta) [cm^{2}/cos(#theta)/Nucleon]");
    SelectionCosTheta.at(2)->Draw("E2");
    SelectionCosTheta.at(0)->SetLineWidth(2);
    SelectionCosTheta.at(0)->SetLineColor(1);
    SelectionCosTheta.at(0)->SetMarkerColor(1);
    SelectionCosTheta.at(0)->Draw("SAME");
    Canvas2->SaveAs(("DiffCrossSectionCosTheta."+FileType).c_str());
    
    TCanvas *Canvas3 = new TCanvas("Theta", "Theta", 1400, 1000);
    Canvas3->cd();
    SelectionTheta.at(2)->SetMaximum(1.2*SelectionTheta.at(0)->GetBinContent(SelectionTheta.at(0)->GetMaximumBin()));
    SelectionTheta.at(2)->SetMinimum(0.0);
    SelectionTheta.at(2)->SetFillColor(46);
    SelectionTheta.at(2)->GetYaxis()->SetTitle("d#sigma/d#theta [cm^{2}/rad/Nucleon]");
    SelectionTheta.at(2)->Draw("E2");
    SelectionTheta.at(0)->SetLineWidth(2);
    SelectionTheta.at(0)->SetLineColor(1);
    SelectionTheta.at(0)->SetMarkerColor(1);
    SelectionTheta.at(0)->Draw("SAME");
    Canvas3->SaveAs(("DiffCrossSectionTheta."+FileType).c_str());
    
    TCanvas *Canvas4 = new TCanvas("Phi", "Phi", 1400, 1000);
    Canvas4->cd();
    SelectionPhi.at(2)->SetMaximum(1.2*SelectionPhi.at(0)->GetBinContent(SelectionPhi.at(0)->GetMaximumBin()));
    SelectionPhi.at(2)->SetMinimum(0.0);
    SelectionPhi.at(2)->SetFillColor(46);
    SelectionPhi.at(2) -> GetYaxis() -> SetTitle("d#sigma/d#phi [cm^{2}/rad/Nucleon]");
    SelectionPhi.at(2)->Draw("E2");
    SelectionPhi.at(0)->SetLineWidth(2);
    SelectionPhi.at(0)->SetLineColor(1);
    SelectionPhi.at(0)->SetMarkerColor(1);
    SelectionPhi.at(0)->Draw("SAME");
    Canvas4->SaveAs(("DiffCrossSectionPhi."+FileType).c_str());
    
    TCanvas *Canvas5 = new TCanvas("Momentum", "Momentum", 1400, 1000);
    Canvas5->cd();
    SelectionMomentum.at(2)->SetMaximum(1.5*SelectionMomentum.at(0)->GetBinContent(SelectionMomentum.at(0)->GetMaximumBin()));
    SelectionMomentum.at(2)->SetMinimum(0.0);
    SelectionMomentum.at(2)->SetFillColor(46);
    SelectionMomentum.at(2) -> GetYaxis() -> SetTitle("d#sigma/dp [cm^{2}/(GeV/c)/Nucleon]");
    SelectionMomentum.at(2)->Draw("E2");
    SelectionMomentum.at(0)->SetLineWidth(2);
    SelectionMomentum.at(0)->SetLineColor(1);
    SelectionMomentum.at(0)->SetMarkerColor(1);
    SelectionMomentum.at(0)->Draw("SAME");
    Canvas5->SaveAs(("DiffCrossSectionMomentum."+FileType).c_str());
    
    // Draw unsmearing matrix
    TCanvas *Canvas6 = new TCanvas("Unsmearing Range", "Unsmearing Range", 1400, 1000);
    Canvas6->cd();
    UMatrixTrackRange->Draw("colz");
    Canvas6->SaveAs(("UnsmearingRange."+FileType).c_str());
    
    TCanvas *Canvas7 = new TCanvas("Unsmearing CosTheta", "Unsmearing CosTheta", 1400, 1000);
    Canvas7->cd();
    UMatrixCosTheta->Draw("colz");
    Canvas7->SaveAs(("UnsmearingCosTheta."+FileType).c_str());
    
    TCanvas *Canvas8 = new TCanvas("Unsmearing Theta", "Unsmearing Theta", 1400, 1000);
    Canvas8->cd();
    UMatrixTheta->Draw("colz");
    Canvas8->SaveAs(("UnsmearingTheta."+FileType).c_str());
    
    TCanvas *Canvas9 = new TCanvas("Unsmearing Phi", "Unsmearing Phi", 1400, 1000);
    Canvas9->cd();
    UMatrixPhi->Draw("colz");
    Canvas9->SaveAs(("UnsmearingPhi."+FileType).c_str());
    
    TCanvas *Canvas10 = new TCanvas("Unsmearing Momentum", "Unsmearing Momentum", 1400, 1000);
    Canvas10->cd();
    UMatrixMomentum->Draw("colz");
    Canvas10->SaveAs(("UnsmearingMomentum."+FileType).c_str());

}

float CalcRange(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2)
{
    return sqrt(pow(x_1-x_2, 2) + pow(y_1-y_2, 2) + pow(z_1-z_2, 2));
}

bool inFV(double x, double y, double z)
{
    if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)) return true;
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
                               1.013E4, 1.469E4, 1.910E4, 3.558E4, 4.326E4, 5.768E4, 7.734E4, 1.060E5, 1.307E5};

    float KEMeV[29] = {10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000, 1400, 2000, 3000, 4000, 
                       8000, 10000, 14000, 20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000};
                       
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
