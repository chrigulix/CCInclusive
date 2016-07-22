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
const double Density = 1.4; // g/cm^3

// Function which calculates the distance between two points
float CalcRange(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2);

// Function which checks if a point is in the FV
bool inFV(double x, double y, double z);

// Add two histogramms with indices First and Last and weight
void AddHistograms(std::vector<TH1F*>& HistVector, unsigned int First, unsigned int Last, float Weight, bool EraseLast = false);

// Normalize Matrix by row
void NormMatrixByRow(TH2F* UMatrix);

// Unsmearing of selected events
void SelectionUnsmearing(TH2F*& UMatrix, TH1F*& SVector);


// CC inlcusive cross section function (main) 
void CCInclCrossSection()
{
    float NumberOfTargets = (FVx - 2*borderx) * (FVy - 2*bordery) * (FVz - 2*borderz) * Density * Avogadro/ArMass*NoNucleons;
    
    // Output file file type
    std::string FileType = "pdf";

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

    double MCPOT = 2.3e20;
    double DataPOT = 4.95e19;
    
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
    
    TH1D* NuMuFlux;
//     TH1F* NuMuFlux = new TH1F("NuMuFlux","NuMuFlux",NumberOfBins,0,3);
    
    TFile* BNBFlux = new TFile("/lheppc46/data/uBData/bnbflux/numode_bnb_470m_r200.root");
    
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

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add("/lheppc46/data/uBData/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_1_Mod.root");
    ChainVec.back() -> Add("/lheppc46/data/uBData/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_2_Mod.root");
    ChainVec.back() -> Add("/lheppc46/data/uBData/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_3_Mod.root");

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add("/lheppc46/data/uBData/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_1_Mod.root");
    ChainVec.back() -> Add("/lheppc46/data/uBData/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_2_Mod.root");

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add("/lheppc46/data/uBData/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_prodgenie_bnb_nu_cosmic_uboone_v05_08_00_Mod.root");

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add("/lheppc46/data/uBData/anatrees/Hist_MC_Truth_prodgenie_bnb_nu_cosmic_uboone_v05_08_00.root");

    GenLabel.push_back("Data On-Beam BNB");
    GenLabel.push_back("Data Off-Beam BNBEXT");
    GenLabel.push_back("MC Selection");
    GenLabel.push_back("MC Selection Backgrounds");
    GenLabel.push_back("MC Selection Truth");
    GenLabel.push_back("MC Truth");
    GenLabel.push_back("Efficiency");
    
    std::vector<float> ScalingFactors;
    ScalingFactors.push_back(1);
    ScalingFactors.push_back(1.2300);
    ScalingFactors.push_back(DataPOT/MCPOT);
    ScalingFactors.push_back(DataPOT/MCPOT);
    ScalingFactors.push_back(DataPOT/MCPOT);
    ScalingFactors.push_back(DataPOT/MCPOT);
    ScalingFactors.push_back(1);

    // Loop over all generation labels
    for(const auto& Label : GenLabel)
    {
        SelectionTrackRange.push_back(new TH1F(("Track Range"+Label).c_str(),"Track Range",NumberOfBins,0,1036.8));
        SelectionTrackRange.back() -> SetStats(0);
        SelectionTrackRange.back() -> GetXaxis() -> SetTitle("Track range [cm]");
        SelectionTrackRange.back() -> GetYaxis() -> SetTitle("No. of events");

        SelectionCosTheta.push_back(new TH1F(("cos#theta-Angle"+Label).c_str(),"cos#theta",NumberOfBins,-1,1));
        SelectionCosTheta.back() -> SetStats(0);
        SelectionCosTheta.back() -> GetXaxis() -> SetTitle("cos(#theta)");
        SelectionCosTheta.back() -> GetYaxis() -> SetTitle("No. of events");
        
        SelectionTheta.push_back(new TH1F(("#theta-Angle"+Label).c_str(),"#theta",NumberOfBins,0,3.142));
        SelectionTheta.back() -> SetStats(0);
        SelectionTheta.back() -> GetXaxis() -> SetTitle("#theta [rad]");
        SelectionTheta.back() -> GetYaxis() -> SetTitle("No. of events");

        SelectionPhi.push_back(new TH1F(("#phi-Angle"+Label).c_str(),"#phi-Angle",NumberOfBins,-3.142,3.142));
        SelectionPhi.back() -> SetStats(0);
        SelectionPhi.back() -> GetXaxis() -> SetTitle("#phi angle [rad]");
        SelectionPhi.back() -> GetYaxis() -> SetTitle("No. of events");

        SelectionMomentum.push_back(new TH1F(("Momentum"+Label).c_str(),"Momentum",NumberOfBins,0,3));
        SelectionMomentum.back() -> SetStats(0);
        SelectionMomentum.back() -> GetXaxis() -> SetTitle("Muon momentum [GeV/c]");
        SelectionMomentum.back() -> GetYaxis() -> SetTitle("No. of events");
    } // loop over generation label
    
    // Initialize unsmearing matrices
    UMatrixTrackRange = new TH2F("Unsmering Matrix Track Range","Unsmering Matrix Track Range",NumberOfBins,0,1036.8,NumberOfBins,0,1036.8);
    UMatrixCosTheta = new TH2F("Unsmering Matrix CosTheta","Unsmering Matrix CosTheta",NumberOfBins,0,-1,NumberOfBins,0,-1);
    UMatrixTheta = new TH2F("Unsmering Matrix Theta","Unsmering Matrix Theta",NumberOfBins,0,3.142,NumberOfBins,0,3.142);
    UMatrixPhi = new TH2F("Unsmering Matrix Phi","Unsmering Matrix Phi",NumberOfBins,0,-3.142,NumberOfBins,0,3.142);
    UMatrixMomentum = new TH2F("Unsmering Matrix Momentum","Unsmering Matrix Momentum",NumberOfBins,0,3,NumberOfBins,0,3);

    // Loop over all files
    for(unsigned int file_no = 0; file_no < ChainVec.size(); file_no++)
    {
        std::cout << "----------------------------------------" << std::endl;

        // Reco entities for all files except truth
        if(file_no < 3)
        {
            ChainVec.at(file_no) -> SetBranchAddress("TrackCand", &TrkID);
            ChainVec.at(file_no) -> SetBranchAddress("VertexCand", &VtxID);

            ChainVec.at(file_no) -> SetBranchAddress("trkmomrange_pandoraNu", TrackMomentum);
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
            ChainVec.at(file_no) -> SetBranchAddress("StartPointz", ZMCTrackStart);
            ChainVec.at(file_no) -> SetBranchAddress("EndPointx", XMCTrackEnd);
            ChainVec.at(file_no) -> SetBranchAddress("EndPointy", YMCTrackEnd);
            ChainVec.at(file_no) -> SetBranchAddress("EndPointz", ZMCTrackEnd);
            ChainVec.at(file_no) -> SetBranchAddress("theta", MCTheta);
            ChainVec.at(file_no) -> SetBranchAddress("enu_truth", NuEnergyTruth);
            ChainVec.at(file_no) -> SetBranchAddress("phi", MCPhi);
            ChainVec.at(file_no) -> SetBranchAddress("Eng", MCEnergy);
            
            // Backtracker information only if there are reco objects
            if(file_no < 3)
            {
                ChainVec.at(file_no) -> SetBranchAddress("trkorigin_pandoraNu", TrkOrigin);
                ChainVec.at(file_no) -> SetBranchAddress("trkpidbestplane_pandoraNu", TrkBestPlane);
            }
        }

        // Loop over all events
        for(unsigned int tree_index = 0; tree_index < ChainVec.at(file_no) -> GetEntries(); tree_index++)
        {
            // Skip corrupted events
            if(file_no == 3 && (tree_index == 11602 || tree_index == 11675 || tree_index == 13510 || tree_index == 33027 || tree_index == 33070 || tree_index == 36239 || tree_index == 44078)) continue;
            
            // Progress indicator
            if(!(tree_index % 1000)) std::cout << "Event\t" << tree_index << "\t of \t" << ChainVec.at(file_no) -> GetEntries() << std::endl;

            // Get tree entry for this event
            ChainVec.at(file_no) -> GetEntry(tree_index);

            // if there are reco products
            if(file_no < 3)
            {
                // Fill histograms as usual
                SelectionTrackRange.at(file_no) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                SelectionCosTheta.at(file_no) -> Fill(cos(TrackTheta[TrkID]));
                SelectionTheta.at(file_no) -> Fill(TrackTheta[TrkID]);
                SelectionPhi.at(file_no) -> Fill(TrackPhi[TrkID]);
                SelectionMomentum.at(file_no) -> Fill(TrackMomentum[TrkID]);
            }

            // if we are looking at the mc selection file
            if(file_no == 2)
            {
                // if event is background
                if(TrkID < 0 || MCVtxID < 0 || TrkOrigin[TrkID][TrkBestPlane[TrkID]] != 1 || nuPDGTruth[MCVtxID] != 14 || CCNCFlag[MCVtxID] == 1 || !inFV(XnuVtxTruth[MCVtxID],YnuVtxTruth[MCVtxID],ZnuVtxTruth[MCVtxID]) || PDGTruth[MCTrkID] != 13)
                {
                    // Fill background histograms
                    SelectionTrackRange.at(file_no+1) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                    SelectionCosTheta.at(file_no+1) -> Fill(cos(TrackTheta[TrkID]));
                    SelectionTheta.at(file_no+1) -> Fill(TrackTheta[TrkID]);
                    SelectionPhi.at(file_no+1) -> Fill(TrackPhi[TrkID]);
                    SelectionMomentum.at(file_no+1) -> Fill(TrackMomentum[TrkID]);
                }
                else // if event is signal and truth
                {
                    // Fill background histograms
                    SelectionTrackRange.at(file_no+2) -> Fill(CalcRange(XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID]));
                    SelectionCosTheta.at(file_no+2) -> Fill(cos(MCTheta[MCTrkID]));
                    SelectionTheta.at(file_no+2) -> Fill(MCTheta[MCTrkID]);
                    SelectionPhi.at(file_no+2) -> Fill(MCPhi[MCTrkID]);
                    SelectionMomentum.at(file_no+2) -> Fill(TrueLeptonMomentum[MCVtxID]);
                    
                    // Fill unsmearing matrix
                    UMatrixTrackRange -> Fill( CalcRange(XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID]),CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]) );
                    UMatrixCosTheta -> Fill( cos(MCTheta[MCTrkID]),cos(TrackTheta[TrkID]) );
                    UMatrixTheta -> Fill( MCTheta[MCTrkID],TrackTheta[TrkID] );
                    UMatrixPhi -> Fill( MCPhi[MCTrkID],TrackPhi[TrkID] );
                    UMatrixMomentum -> Fill( TrueLeptonMomentum[MCVtxID],TrackMomentum[TrkID] );
                }
            } // if mc selection file
            else if(file_no == 3 && MCTrkID > -1 && nuPDGTruth[MCVtxID] == 14)
            {
                // Fill background histograms
                SelectionTrackRange.at(file_no+2) -> Fill(CalcRange(XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID]));
                SelectionCosTheta.at(file_no+2) -> Fill(cos(MCTheta[MCTrkID]));
                SelectionTheta.at(file_no+2) -> Fill(MCTheta[MCTrkID]);
                SelectionPhi.at(file_no+2) -> Fill(MCPhi[MCTrkID]);
                SelectionMomentum.at(file_no+2) -> Fill(TrueLeptonMomentum[MCVtxID]);
            }
            
        } // Event loop

        // Reset branch addresses to avoid problems
        ChainVec.at(file_no) -> ResetBranchAddresses();

    } // file loop
    
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

    // loop over scaleing factors 
    for(unsigned int scale_no = 0; scale_no < ScalingFactors.size(); scale_no++)
    {
        // Scale histograms
        SelectionTrackRange.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionCosTheta.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionTheta.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionPhi.at(scale_no)->Scale(ScalingFactors.at(scale_no));
        SelectionMomentum.at(scale_no)->Scale(ScalingFactors.at(scale_no));
    } // scaling loop
    
    SelectionTrackRange.back()->Divide(SelectionTrackRange.at(4),SelectionTrackRange.at(5));
    SelectionCosTheta.back()->Divide(SelectionCosTheta.at(4),SelectionCosTheta.at(5));
    SelectionTheta.back()->Divide(SelectionTheta.at(4),SelectionTheta.at(5));
    SelectionPhi.back()->Divide(SelectionPhi.at(4),SelectionPhi.at(5));
    SelectionMomentum.back()->Divide(SelectionMomentum.at(4),SelectionMomentum.at(5));
    
//     TEfficiency* EffMomentum = new TEfficiency( *SelectionMomentum.at(4),*SelectionMomentum.at(5) );
    
//     EffTrackRange.BayesDivide(SelectionTrackRange.at(3),SelectionTrackRange.at(2));
//     EffCosTheta.BayesDivide(SelectionCosTheta.at(3),SelectionCosTheta.at(2));
//     EffPhi.BayesDivide(SelectionPhi.at(3),SelectionPhi.at(2));
//     EffMomentum.BayesDivide(SelectionMomentum.at(3),SelectionMomentum.at(2));
    
    // Normalize matrices by row
    NormMatrixByRow(UMatrixTrackRange);
    NormMatrixByRow(UMatrixCosTheta);
    NormMatrixByRow(UMatrixTheta);
    NormMatrixByRow(UMatrixPhi);
    NormMatrixByRow(UMatrixMomentum);

    // Subtract offbeam from onbeam data, overwrite onbeam- and erase offbeam histogram
    AddHistograms(SelectionTrackRange,0,1,-1,1);
    AddHistograms(SelectionCosTheta,0,1,-1,1);
    AddHistograms(SelectionTheta,0,1,-1,1);
    AddHistograms(SelectionPhi,0,1,-1,1);
    AddHistograms(SelectionMomentum,0,1,-1,1);
    
    // Subtract backgrounds from data
    AddHistograms(SelectionTrackRange,0,2,-1);
    AddHistograms(SelectionCosTheta,0,2,-1);
    AddHistograms(SelectionTheta,0,2,-1);
    AddHistograms(SelectionPhi,0,2,-1);
    AddHistograms(SelectionMomentum,0,2,-1);
    
    // Subtract backgrounds from mc selection and erase bgr
    AddHistograms(SelectionTrackRange,1,2,-1,1);
    AddHistograms(SelectionPhi,1,2,-1,1);
    AddHistograms(SelectionCosTheta,1,2,-1,1);
    AddHistograms(SelectionTheta,1,2,-1,1);
    AddHistograms(SelectionMomentum,1,2,-1,1);
    
    // Unsmearing loop
    for(unsigned int hist_no = 0; hist_no < 2; hist_no++)
    {
        SelectionUnsmearing(UMatrixTrackRange,SelectionTrackRange.at(hist_no));
        SelectionUnsmearing(UMatrixCosTheta,SelectionCosTheta.at(hist_no));
        SelectionUnsmearing(UMatrixTheta,SelectionTheta.at(hist_no));
        SelectionUnsmearing(UMatrixPhi,SelectionPhi.at(hist_no));
        SelectionUnsmearing(UMatrixMomentum,SelectionMomentum.at(hist_no));
        
        SelectionTrackRange.at(hist_no)->Divide(SelectionTrackRange.back());
        SelectionCosTheta.at(hist_no)->Divide(SelectionCosTheta.back());
        SelectionTheta.at(hist_no)->Divide(SelectionTheta.back());
        SelectionPhi.at(hist_no)->Divide(SelectionPhi.back());
        SelectionMomentum.at(hist_no)->Divide(SelectionMomentum.back());
        
        SelectionTrackRange.at(hist_no)->Scale(1/NumberOfTargets/SelectionTrackRange.at(hist_no)->GetBinWidth(1)/IntegratedFlux);
        SelectionCosTheta.at(hist_no)->Scale(1/NumberOfTargets/SelectionCosTheta.at(hist_no)->GetBinWidth(1)/IntegratedFlux);
        SelectionTheta.at(hist_no)->Scale(1/NumberOfTargets/SelectionTheta.at(hist_no)->GetBinWidth(1)/IntegratedFlux);
        SelectionPhi.at(hist_no)->Scale(1/NumberOfTargets/SelectionPhi.at(hist_no)->GetBinWidth(1)/IntegratedFlux);
        SelectionMomentum.at(hist_no)->Scale(1/NumberOfTargets/SelectionMomentum.at(hist_no)->GetBinWidth(1)/IntegratedFlux);
    }
    
    TCanvas *Canvas0 = new TCanvas("Test", "Test", 1400, 1000);
    Canvas0->cd();
    SelectionTheta.back()->Draw();

    // Draw histogram
    TCanvas *Canvas1 = new TCanvas("Range", "Range", 1400, 1000);
    Canvas1->cd();
    SelectionTrackRange.at(1)->SetMaximum(1.2*SelectionTrackRange.at(0)->GetBinContent(SelectionTrackRange.at(0)->GetMaximumBin()));
    SelectionTrackRange.at(1)->SetMinimum(0.0);
    SelectionTrackRange.at(1)->SetFillColor(46);
    SelectionTrackRange.at(1)->Draw("E2");
//     StackBgrTrackRange->Draw("SAME");
    SelectionTrackRange.at(0)->SetLineWidth(2);
    SelectionTrackRange.at(0)->SetLineColor(1);
    SelectionTrackRange.at(0)->SetMarkerColor(1);
    SelectionTrackRange.at(0)->Draw("SAME");
//     LegendMC->Draw();
//     LegendBGR->Draw();
//     TextPreliminary.Draw();
//     TextSelection.Draw();
    Canvas1->SaveAs(("ScaledOn-OffBeamSelRange."+FileType).c_str());
    
    TCanvas *Canvas2 = new TCanvas("CosTheta", "CosTheta", 1400, 1000);
    Canvas2->cd();
    SelectionCosTheta.at(1)->SetMaximum(1.2*SelectionCosTheta.at(0)->GetBinContent(SelectionCosTheta.at(0)->GetMaximumBin()));
    SelectionCosTheta.at(1)->SetMinimum(0.0);
    SelectionCosTheta.at(1)->SetFillColor(46);
    SelectionCosTheta.at(1)->Draw("E2");
//     StackBgrTrackRange->Draw("SAME");
    SelectionCosTheta.at(0)->SetLineWidth(2);
    SelectionCosTheta.at(0)->SetLineColor(1);
    SelectionCosTheta.at(0)->SetMarkerColor(1);
    SelectionCosTheta.at(0)->Draw("SAME");
//     LegendMC->Draw();
//     LegendBGR->Draw();
//     TextPreliminary.Draw();
//     TextSelection.Draw();
    Canvas2->SaveAs(("ScaledOn-OffBeamSelCosTheta."+FileType).c_str());
    
    TCanvas *Canvas3 = new TCanvas("Theta", "Theta", 1400, 1000);
    Canvas3->cd();
    SelectionTheta.at(1)->SetMaximum(1.2*SelectionTheta.at(0)->GetBinContent(SelectionTheta.at(0)->GetMaximumBin()));
    SelectionTheta.at(1)->SetMinimum(0.0);
    SelectionTheta.at(1)->SetFillColor(46);
    SelectionTheta.at(1)->Draw("E2");
//     StackBgrTrackRange->Draw("SAME");
    SelectionTheta.at(0)->SetLineWidth(2);
    SelectionTheta.at(0)->SetLineColor(1);
    SelectionTheta.at(0)->SetMarkerColor(1);
    SelectionTheta.at(0)->Draw("SAME");
//     LegendMC->Draw();
//     LegendBGR->Draw();
//     TextPreliminary.Draw();
//     TextSelection.Draw();
    Canvas3->SaveAs(("ScaledOn-OffBeamSelTheta."+FileType).c_str());
    
    TCanvas *Canvas4 = new TCanvas("Phi", "Phi", 1400, 1000);
    Canvas4->cd();
    SelectionPhi.at(1)->SetMaximum(1.2*SelectionPhi.at(0)->GetBinContent(SelectionPhi.at(0)->GetMaximumBin()));
    SelectionPhi.at(1)->SetMinimum(0.0);
    SelectionPhi.at(1)->SetFillColor(46);
    SelectionPhi.at(1)->Draw("E2");
//     StackBgrTrackRange->Draw("SAME");
    SelectionPhi.at(0)->SetLineWidth(2);
    SelectionPhi.at(0)->SetLineColor(1);
    SelectionPhi.at(0)->SetMarkerColor(1);
    SelectionPhi.at(0)->Draw("SAME");
//     LegendMC->Draw();
//     LegendBGR->Draw();
//     TextPreliminary.Draw();
//     TextSelection.Draw();
    Canvas4->SaveAs(("ScaledOn-OffBeamSelPhi."+FileType).c_str());
    
    TCanvas *Canvas5 = new TCanvas("Momentum", "Momentum", 1400, 1000);
    Canvas5->cd();
    SelectionMomentum.at(1)->SetMaximum(1.2*SelectionMomentum.at(0)->GetBinContent(SelectionMomentum.at(0)->GetMaximumBin()));
    SelectionMomentum.at(1)->SetMinimum(0.0);
    SelectionMomentum.at(1)->SetFillColor(46);
    SelectionMomentum.at(1)->Draw("E2");
//     StackBgrTrackRange->Draw("SAME");
    SelectionMomentum.at(0)->SetLineWidth(2);
    SelectionMomentum.at(0)->SetLineColor(1);
    SelectionMomentum.at(0)->SetMarkerColor(1);
    SelectionMomentum.at(0)->Draw("SAME");
//     LegendMC->Draw();
//     LegendBGR->Draw();
//     TextPreliminary.Draw();
//     TextSelection.Draw();
    Canvas5->SaveAs(("ScaledOn-OffBeamSelMomentum."+FileType).c_str());
    
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
    if (HistVector.size() >= Last)
    {
        // Add histograms
        HistVector.at(First) -> Add(HistVector.at(Last), Weight);

        // Erase last histogram if flag is set
        if(EraseLast)
        {
            delete HistVector.at(Last);
            HistVector.erase(HistVector.begin() + Last);
        }
    }
    else // if nothing can be added
    {
        std::cout << "Histograms not added!" << std::endl;
    }
}

void NormMatrixByRow(TH2F* UMatrix)
{
    // loop over ybins of the unsmearing matrices
    for(unsigned int ybin = 1; ybin <= UMatrix->GetNbinsY(); ybin++)
    {
        float NormFact = 0;
        
        // loop over xbins (row)
        for(unsigned int xbin = 1; xbin <= UMatrix->GetNbinsX(); xbin++)
        {
            // Add row entry to normalization factor
            NormFact += UMatrix->GetBinContent(xbin,ybin);
        } // xbin loop
        
        // loop over xbins (row)
        for(unsigned int xbin = 1; xbin <= UMatrix->GetNbinsX(); xbin++)
        {
            // Normalize entire row of the matrix
            if(NormFact) UMatrix->SetBinContent(xbin,ybin,UMatrix->GetBinContent(xbin,ybin)/NormFact) ;
        }// xbin loop
    }// ybin loop
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