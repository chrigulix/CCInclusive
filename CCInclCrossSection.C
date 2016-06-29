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

// Function which calculates the distance between two points
float CalcRange(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2);

// Function which checks if a point is in the FV
bool inFV(double x, double y, double z);

// Add two histogramms with indices First and Last and weight
void AddHistograms(std::vector<TH1F*>& HistVector, unsigned int First, unsigned int Last, float Weight, bool EraseLast = false);


void CCInclCrossSection()
{
    // Data input file vector
    std::vector<TChain*> ChainVec;

    // Histogram Vectors
    std::vector<TH1F*> SelectionTrackRange;
    std::vector<TH1F*> SelectionCosTheta;
    std::vector<TH1F*> SelectionPhi;
    std::vector<TH1F*> SelectionMomentum;

    // Unsemaring Matrix
    TH2F* UMatrixTrackRange;
    TH2F* UMatrixCosTheta;
    TH2F* UMatrixPhi;
    TH2F* UMatrixMomentum;

    size_t NumberOfBins = 20;

    double MCPOT = 2.3e23;
    double DataPOT = 4.95e19;

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

        SelectionPhi.push_back(new TH1F(("#phi-Angle"+Label).c_str(),"#phi-Angle",NumberOfBins,-3.142,3.142));
        SelectionPhi.back() -> SetStats(0);
        SelectionPhi.back() -> GetXaxis() -> SetTitle("#phi angle [rad]");
        SelectionPhi.back() -> GetYaxis() -> SetTitle("No. of events");

        SelectionMomentum.push_back(new TH1F(("Momentum"+Label).c_str(),"Momentum",NumberOfBins,0,3));
        SelectionMomentum.back() -> SetStats(0);
        SelectionMomentum.back() -> GetXaxis() -> SetTitle("Muon momentum [GeV/c]");
        SelectionMomentum.back() -> GetYaxis() -> SetTitle("No. of events");
    } // loop over generation label

    // Loop over all files
    for(unsigned int file_no = 0; file_no < ChainVec.size(); file_no++)
    {
        // Reco entities for all files except truth
        if(file_no != ChainVec.size()-1)
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
            ChainVec.at(file_no) -> SetBranchAddress("trkorigin_pandoraNu", TrkOrigin);
            ChainVec.at(file_no) -> SetBranchAddress("trkpidbestplane_pandoraNu", TrkBestPlane);
        }

        // Loop over all events
        for(unsigned int tree_index = 0; tree_index < ChainVec.at(file_no) -> GetEntries(); tree_index++)
        {
            // Progress indicator
            if(!(tree_index % 1000)) std::cout << "Event\t" << tree_index << "\t of \t" << ChainVec.at(file_no) -> GetEntries() << std::endl;

            // Get tree entry for this event
            ChainVec.at(file_no) -> GetEntry(tree_index);

            // if there are reco products
            if(file_no <= 2)
            {
                // Fill histograms as usual
                SelectionTrackRange.at(file_no) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                SelectionCosTheta.at(file_no) -> Fill(cos(TrackTheta[TrkID]));
                SelectionPhi.at(file_no) -> Fill(TrackPhi[TrkID]);
                SelectionMomentum.at(file_no) -> Fill(TrackMomentum[TrkID]);
            }

            // if we are looking at the mc selection file
            if(file_no == 2)
            {
                // if event is background
                if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] != 1 || nuPDGTruth[MCVtxID] != 14 || CCNCFlag[MCVtxID] == 1 || !inFV(XnuVtxTruth[MCVtxID],YnuVtxTruth[MCVtxID],ZnuVtxTruth[MCVtxID]) || PDGTruth[MCTrkID] != 13)
                {
                    // Fill background histograms
                    SelectionTrackRange.at(file_no+1) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                    SelectionCosTheta.at(file_no+1) -> Fill(cos(TrackTheta[TrkID]));
                    SelectionPhi.at(file_no+1) -> Fill(TrackPhi[TrkID]);
                    SelectionMomentum.at(file_no+1) -> Fill(TrackMomentum[TrkID]);
                }
                else // if event is signal and truth
                {
                    // Fill background histograms
                    SelectionTrackRange.at(file_no+2) -> Fill(CalcRange(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                    SelectionCosTheta.at(file_no+2) -> Fill(cos(TrackTheta[TrkID]));
                    SelectionPhi.at(file_no+2) -> Fill(TrackPhi[TrkID]);
                    SelectionMomentum.at(file_no+2) -> Fill(TrackMomentum[TrkID]);
                }
            } // if mc selection file
            else if(file_no == 3)
            {
                // Fill background histograms
                SelectionTrackRange.at(file_no+2) -> Fill(CalcRange(XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID]));
                SelectionCosTheta.at(file_no+2) -> Fill(cos(MCTheta[MCTrkID]));
                SelectionPhi.at(file_no+2) -> Fill(MCPhi[MCTrkID]);
                SelectionMomentum.at(file_no+2) -> Fill(TrueLeptonMomentum[MCVtxID]);
            }

        } // Event loop

        // Reset branch addresses to avoid problems
        ChainVec.at(file_no) -> ResetBranchAddresses();

    } // file loop


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
    if (HistVector.size() >= Last && First < Last)
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
