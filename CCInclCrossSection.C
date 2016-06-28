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

void CCInclCrossSection()
{
    // Data input file vector
    std::vector<TChain*> ChainVec;
    
    // Histogram Vectors
    std::vector<TH1F*> SelectionTrackRange;
    std::vector<TH1F*> SelectionCosTheta;
    std::vector<TH1F*> SelectionPhi;
    std::vector<TH1F*> SelectionMomentum;

    size_t NumberOfBins = 20;
    
    int TrkID;
    int VtxID;
    int MCTrkID;
    int MCVtxID;
    
    int CCNCFlag[10];
    int TruthMode[10];
    int NuPDGTruth[10];
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
    GenLabel.push_back("MC Backgrounds");
    GenLabel.push_back("MC Truth");
    
    // Loop over all generated
    for(const auto& Label : GenLabel)
    {
        SelectionTrackRange.push_back(new TH1F(("Track Range"+Label).c_str(),"Track Range of Selected Track",NumberOfBins,0,1036.8));
        SelectionTrackRange.back()->SetStats(0);
        SelectionTrackRange.back()->GetXaxis()->SetTitle("Track range [cm]");
        SelectionTrackRange.back()->GetYaxis()->SetTitle("No. of events");
        
        SelectionCosTheta.push_back(new TH1F(("cos#theta-Angle"+Label).c_str(),"cos#theta of Selected Track",NumberOfBins,-1,1));
        SelectionCosTheta.back()->SetStats(0);
        SelectionCosTheta.back()->GetXaxis()->SetTitle("cos(#theta)");
        SelectionCosTheta.back()->GetYaxis()->SetTitle("No. of events");
        
        SelectionPhi.push_back(new TH1F(("#phi-Angle"+Label).c_str(),"#phi-Angle of Selected Track",NumberOfBins,-3.142,3.142));
        SelectionPhi.back()->SetStats(0);
        SelectionPhi.back()->GetXaxis()->SetTitle("#phi angle [rad]");
        SelectionPhi.back()->GetYaxis()->SetTitle("No. of events");
        
        SelectionMomentum.push_back(new TH1F(("Momentum"+Label).c_str(),"Momentum of Selected Track",NumberOfBins,0,3));
        SelectionMomentum.back()->SetStats(0);
        SelectionMomentum.back()->GetXaxis()->SetTitle("Muon momentum [GeV/c]");
        SelectionMomentum.back()->GetYaxis()->SetTitle("No. of events");
    }
    
    for(unsigned int file_no = 0; file_no < ChainVec.size(); file_no++)
    {
        ChainVec.at(file_no) -> SetBranchAddress("TrackCand", &TrkID);
        ChainVec.at(file_no) -> SetBranchAddress("VertexCand", &VtxID);

        ChainVec.at(file_no) -> SetBranchAddress("MCTrackCand", &MCTrkID);
        ChainVec.at(file_no) -> SetBranchAddress("MCVertexCand", &MCVtxID);
        ChainVec.at(file_no) -> SetBranchAddress("ccnc_truth", CCNCFlag);
        ChainVec.at(file_no) -> SetBranchAddress("mode_truth", TruthMode);
        ChainVec.at(file_no) -> SetBranchAddress("nuPDG_truth", NuPDGTruth);
        ChainVec.at(file_no) -> SetBranchAddress("pdg", PDGTruth);
        ChainVec.at(file_no) -> SetBranchAddress("enu_truth", NuEnergyTruth);
        ChainVec.at(file_no) -> SetBranchAddress("nuvtxx_truth", nuvtxx_truth);
        ChainVec.at(file_no) -> SetBranchAddress("nuvtxy_truth", nuvtxy_truth);
        ChainVec.at(file_no) -> SetBranchAddress("nuvtxz_truth", nuvtxz_truth);
        ChainVec.at(file_no) -> SetBranchAddress("trkorigin_pandoraNu", TrkOrigin);
        ChainVec.at(file_no) -> SetBranchAddress("trkpidbestplane_pandoraNu", TrkBestPlane);

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
    
    
}