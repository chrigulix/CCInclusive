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

const double Pi = 3.14159265358979323846; // Pi

void MCTruthTest()
{
    size_t NumberOfBins = 20;
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

    TChain* Chain = new TChain("anatree");
    Chain -> Add("/home/christoph/anatrees/CCInclusiveNote/Hist_Track_pandoraNu_Vertex_pandoraNu_prodgenie_bnb_nu_cosmic_uboone_v05_08_00_Mod.root");
//     Chain -> Add("/home/christoph/anatrees/CCInclusiveNote/Hist_Track_pandoraNu_Vertex_pandoraNu_MA_v05_08_00_Mod.root");
//     Chain -> Add("/home/christoph/anatrees/CCInclusiveNote/Hist_Track_pandoraNu_Vertex_pandoraNu_TEM_v05_08_00_Mod.root");
//     Chain -> Add("/home/christoph/anatrees/CCInclusiveNote/Hist_Track_pandoraNu_Vertex_pandoraNu_MEC_v05_08_00_Mod.root");


    // Reco properties
    Chain -> SetBranchAddress("TrackCand", &TrkID);
    Chain -> SetBranchAddress("VertexCand", &VtxID);
    Chain -> SetBranchAddress("trklen_pandoraNu", TrackLength);
    Chain -> SetBranchAddress("trktheta_pandoraNu", TrackTheta);
    Chain -> SetBranchAddress("trkphi_pandoraNu",TrackPhi);
    Chain -> SetBranchAddress("trkstartx_pandoraNu",XTrackStart);
    Chain -> SetBranchAddress("trkstarty_pandoraNu",YTrackStart);
    Chain -> SetBranchAddress("trkstartz_pandoraNu",ZTrackStart);
    Chain -> SetBranchAddress("trkendx_pandoraNu",XTrackEnd);
    Chain -> SetBranchAddress("trkendy_pandoraNu",YTrackEnd);
    Chain -> SetBranchAddress("trkendz_pandoraNu",ZTrackEnd);
    Chain -> SetBranchAddress("vtxx_pandoraNu", XVertexPosition);
    Chain -> SetBranchAddress("vtxy_pandoraNu", YVertexPosition);
    Chain -> SetBranchAddress("vtxz_pandoraNu", ZVertexPosition);

    // Truth properties
    Chain -> SetBranchAddress("MCTrackCand", &MCTrkID);
    Chain -> SetBranchAddress("MCVertexCand", &MCVtxID);
    Chain -> SetBranchAddress("ccnc_truth", CCNCFlag);
    Chain -> SetBranchAddress("mode_truth", TruthMode);
    Chain -> SetBranchAddress("pdg", PDGTruth);
    Chain -> SetBranchAddress("enu_truth", NuEnergyTruth);
    Chain -> SetBranchAddress("lep_mom_truth", TrueLeptonMomentum);
    Chain -> SetBranchAddress("mcevts_truth", &mcevts_truth);
    Chain -> SetBranchAddress("nuvtxx_truth", XnuVtxTruth);
    Chain -> SetBranchAddress("nuvtxy_truth", YnuVtxTruth);
    Chain -> SetBranchAddress("nuvtxz_truth", ZnuVtxTruth);
    Chain -> SetBranchAddress("nuPDG_truth", nuPDGTruth);
    Chain -> SetBranchAddress("geant_list_size", &NumberOfMCTracks);
    Chain -> SetBranchAddress("TrackId", MCTrackID);
    Chain -> SetBranchAddress("MCTruthIndex", MCTrueIndex);
    Chain -> SetBranchAddress("StartPointx", XMCTrackStart);
    Chain -> SetBranchAddress("StartPointy", YMCTrackStart);
    Chain -> SetBranchAddress("StartPointz", ZMCTrackStart);
    Chain -> SetBranchAddress("EndPointx", XMCTrackEnd);
    Chain -> SetBranchAddress("EndPointy", YMCTrackEnd);
    Chain -> SetBranchAddress("EndPointz", ZMCTrackEnd);
    Chain -> SetBranchAddress("theta", MCTheta);
    Chain -> SetBranchAddress("enu_truth", NuEnergyTruth);
    Chain -> SetBranchAddress("phi", MCPhi);
    Chain -> SetBranchAddress("Eng", MCEnergy);

    Chain -> SetBranchAddress("trkorigin_pandoraNu", TrkOrigin);
    Chain -> SetBranchAddress("trkidtruth_pandoraNu",TrackIDTruth);
    Chain -> SetBranchAddress("trkpidbestplane_pandoraNu", TrkBestPlane);

    TH1F* SelectionTheta = new TH1F("#theta-Angle","#theta-Angle",NumberOfBins,0,180);
    SelectionTheta -> SetStats(0);
    SelectionTheta -> GetXaxis() -> SetTitle("Muon #theta-Angle [#circ]");
    SelectionTheta -> GetYaxis() -> SetTitle("No. of events");

    unsigned int cosmic = 0;

    for(unsigned int tree_index = 0; tree_index < Chain -> GetEntries(); tree_index++)
    {
//         if(tree_index == 5455) continue;

        Chain -> GetEntry(tree_index);

        // Match track and vertex to true information
        int MCTrackCandidate = -1;
        int MCVertexCandidate = -1;
        if(TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1)
        {
            for(unsigned track_no = 0; track_no < NumberOfMCTracks; track_no++)
            {
                if(MCTrackID[track_no] == TrackIDTruth[TrkID][TrkBestPlane[TrkID]])
                {
                    MCTrackCandidate = track_no;
                    MCVertexCandidate = MCTrueIndex[track_no];

                    if(MCTrackCandidate != MCTrkID || MCTrackCandidate < 0)
                    {
                        std::cout << "Event: " << tree_index << "\t| Track Candidate: \t (S) " << MCTrkID << " \t (C) " << MCTrackCandidate << "\tVert: " <<  MCVertexCandidate << std::endl;
                    }
                    if(MCVertexCandidate != 0 || MCVertexCandidate != MCVtxID)
                    {
                        std::cout << "Event: " << tree_index << "\t| Vertex Candidate: \t (S) " << MCVtxID << " \t (C) " << MCVertexCandidate << std::endl;
                    }
                }
            }
        }
        else cosmic++;
    }

    std::cout << cosmic << " with no neutrino origin!" << std::endl;
//     SelectionTheta -> Draw();
}
