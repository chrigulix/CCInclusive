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

// using namespace std;

//This defines our current settings for the fiducial volume
double FVx = 256.35;
double FVy = 233;
double FVz = 1036.8;
double borderx = 10.;
double bordery = 20.;
double borderz = 10.;

//This function returns if a 3D point is within the fiducial volume
bool inFV(double x, double y, double z) {
    if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)) return true;
    else return false;
}

// Main function
int MCSelection(std::string GeneratorName, unsigned int ThreadNumber, unsigned int NumberOfThreads)
{

    std::string Version = "v05_08_00";

//     std::string GeneratorName = "prodgenie_bnb_nu_cosmic";
//     std::string GeneratorName = "prodgenie_bnb_nu";
//     std::string GeneratorName = "prodcosmics_corsika_inTime";
//     std::string GeneratorName = "prodgenie_bnb_nu_cosmic_sc_uboone";
//     std::string GeneratorName = "data_onbeam_bnb";
//     std::string GeneratorName = "data_offbeam_bnbext";
//     std::string GeneratorName = "prodgenie_bnb_nu_cosmic_uboone";

    std::string FileNumberStr;

    if(NumberOfThreads > 1)
    {
        FileNumberStr = "_" + std::to_string(ThreadNumber);
    }
    else
    {
        FileNumberStr = "";
    }

    TChain *treenc = new TChain("analysistree/anatree");

    if(GeneratorName == "data_onbeam_bnb")
    {
        treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/onbeam_data_bnbSWtrigger/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
    }
    else if(GeneratorName == "data_offbeam_bnbext")
    {
        treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/offbeam_data_bnbSWtrigger/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
    }
    else if(GeneratorName == "prodgenie_bnb_nu_cosmic_uboone")
    {
        treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/MC_BNB_Cosmic/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
    }
    else if(GeneratorName == "prodgenie_bnb_nu_cosmic_uboone_field")
    {
        std::ifstream FileNames("/pnfs/uboone/persistent/users/sowjanya/v05_08_00/bnbpluscosmics_nominal_scOn_100k/ana/filesana.list");

        std::string FileName;

        while(std::getline(FileNames,FileName))
        {
            std::cout << FileName << std::endl;
            treenc -> Add((FileName).c_str());
        }
    }
    else
    {
        treenc -> Add( ("/lheppc46/data/uBData/anatrees/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
    }

    //maximum array sizes
    const int maxtracks = 10000;
    const int maxmc = 10;

    //MC truth
    Int_t           mcevts_truth; //neutrino interactions per event
    Float_t         nuvtxx_truth[maxmc]; //true vertex x (in cm)
    Float_t         nuvtxy_truth[maxmc];
    Float_t         nuvtxz_truth[maxmc];
    Int_t           ccnc_truth[maxmc]; //CC = 0, NC = 1
    Int_t           nuPDG_truth[maxmc]; //true neutrino pdg code. numu = 14
    Float_t         NuEnergyTruth[maxmc];
    Float_t         TrueLeptonMomentum[maxmc];
    Int_t           mode_truth[maxmc]; //QE = 0, RES = 1, DIS = 2
    Int_t	    PDG_truth[maxtracks];

    Int_t	   NumberOfMCTracks;

    Float_t	   XMCTrackStart[maxtracks];
    Float_t	   YMCTrackStart[maxtracks];
    Float_t	   ZMCTrackStart[maxtracks];

    Float_t	   XMCTrackEnd[maxtracks];
    Float_t	   YMCTrackEnd[maxtracks];
    Float_t	   ZMCTrackEnd[maxtracks];

    Float_t        MCTheta[maxtracks];
    Float_t        MCPhi[maxtracks];
    Float_t        MCEnergy[maxtracks];
    Int_t          MCTrueIndex[maxtracks];
    
    // Open output file
    TFile* OutputFile = new TFile(("rootfiles/Hist_MC_Truth_"+GeneratorName+"_"+Version+FileNumberStr+".root").c_str(),"RECREATE");

    treenc -> SetBranchStatus("*",0);
    treenc -> SetBranchStatus("mcevts_truth", 1);
    treenc -> SetBranchStatus("geant_list_size", 1);
    treenc -> SetBranchStatus("nuvtxx_truth", 1);
    treenc -> SetBranchStatus("nuvtxy_truth", 1);
    treenc -> SetBranchStatus("nuvtxz_truth", 1);
    treenc -> SetBranchStatus("ccnc_truth", 1);
    treenc -> SetBranchStatus("nuPDG_truth", 1);
    treenc -> SetBranchStatus("pdg", 1);
    treenc -> SetBranchStatus("mode_truth", 1);
    treenc -> SetBranchStatus("enu_truth", 1);
    treenc -> SetBranchStatus("lep_mom_truth", 1);
    treenc -> SetBranchStatus("StartPointx", 1);
    treenc -> SetBranchStatus("StartPointy", 1);
    treenc -> SetBranchStatus("StartPointz", 1);
    treenc -> SetBranchStatus("EndPointx", 1);
    treenc -> SetBranchStatus("EndPointy", 1);
    treenc -> SetBranchStatus("EndPointz", 1);
    treenc -> SetBranchStatus("theta", 1);
    treenc -> SetBranchStatus("phi", 1);
    treenc -> SetBranchStatus("Eng", 1);
    treenc -> SetBranchStatus("MCTruthIndex", 1);

    treenc -> SetBranchAddress("mcevts_truth", &mcevts_truth);
    treenc -> SetBranchAddress("geant_list_size", &NumberOfMCTracks);
    treenc -> SetBranchAddress("nuvtxx_truth", nuvtxx_truth);
    treenc -> SetBranchAddress("nuvtxy_truth", nuvtxy_truth);
    treenc -> SetBranchAddress("nuvtxz_truth", nuvtxz_truth);
    treenc -> SetBranchAddress("ccnc_truth", ccnc_truth);
    treenc -> SetBranchAddress("nuPDG_truth", nuPDG_truth);
    treenc -> SetBranchAddress("pdg", PDG_truth);
    treenc -> SetBranchAddress("mode_truth", mode_truth);
    treenc -> SetBranchAddress("enu_truth", NuEnergyTruth);
    treenc -> SetBranchAddress("lep_mom_truth", TrueLeptonMomentum);
    treenc -> SetBranchAddress("StartPointx", XMCTrackStart);
    treenc -> SetBranchAddress("StartPointy", YMCTrackStart);
    treenc -> SetBranchAddress("StartPointz", ZMCTrackStart);
    treenc -> SetBranchAddress("EndPointx", XMCTrackEnd);
    treenc -> SetBranchAddress("EndPointy", YMCTrackEnd);
    treenc -> SetBranchAddress("EndPointz", ZMCTrackEnd);
    treenc -> SetBranchAddress("theta", MCPhi);
    treenc -> SetBranchAddress("phi", MCTheta);
    treenc -> SetBranchAddress("Eng", MCEnergy);
    treenc -> SetBranchAddress("MCTruthIndex", MCTrueIndex);
    
    // Make a clone tree which gets filled
    TTree *SelectionTree = treenc->CloneTree(0);

    int ntrue = 0;
    int MCTrackCandidate;
    int MCVertexCandidate;

    double MCTrackToMCVtxDist = 0.5; //cm. distance between mc track start and mc vertex

    unsigned int NumberOfSignalTruth;

    TBranch* BrMCTrackCand = SelectionTree->Branch("MCTrackCand",&MCTrackCandidate,"MCTrackCand/I");
    TBranch* BrMCVtxCand = SelectionTree->Branch("MCVertexCand",&MCVertexCandidate,"MCVertexCand/I");

    unsigned long int Size = treenc -> GetEntries();

    // Set start and end event number for multiple threads
    unsigned long int StartEvent = Size*(ThreadNumber - 1)/NumberOfThreads;
    unsigned long int EndEvent = Size*ThreadNumber/NumberOfThreads;

    std::cout << "number of events used is: " << EndEvent-StartEvent << " of " << Size << std::endl;
    //Event Loop
    for(int i = StartEvent; i < EndEvent; i++)
    {
        if((i == 2633 || i == 22955) && GeneratorName == "prodgenie_bnb_nu_cosmic_uboone_field" ) continue;
        if(i%1 == 0) std::cout << "\t... " << i << std::endl;

        // Get tree entries
        treenc -> GetEntry(i);

        MCTrackCandidate = -1;
        MCVertexCandidate = -1;

        // Loop over all MC neutrino vertices
        for(unsigned vertex_no = 0; vertex_no < mcevts_truth; vertex_no++)
        {
            // Check if there is a numuCC vertex in the FV
            if( nuPDG_truth[vertex_no] == 14 && ccnc_truth[vertex_no] == 0 && inFV(nuvtxx_truth[vertex_no],nuvtxy_truth[vertex_no],nuvtxz_truth[vertex_no]) )
            {
                // Increase truth count
                NumberOfSignalTruth++;

                // Loop over all MC particles
                for(unsigned track_no = 0; track_no < NumberOfMCTracks; track_no++)
                {
                    // If the a muon is not contained in a singel neutrino event, set mc-track contained flag to false
                    if( PDG_truth[track_no] == 13 && MCTrueIndex[track_no] == vertex_no
                        && sqrt(pow(XMCTrackStart[track_no] - nuvtxx_truth[vertex_no],2) + pow(YMCTrackStart[track_no] - nuvtxy_truth[vertex_no],2) + pow(ZMCTrackStart[track_no] - nuvtxz_truth[vertex_no],2)) < MCTrackToMCVtxDist )
                    {
                        // Fill track candidate index
                        MCTrackCandidate = track_no;
                        MCVertexCandidate = vertex_no;
                        
                        // Fill selection Tree
                        SelectionTree->Fill();
                    }
                } // MC particle loop
            } // If numuCC in FV
        } // MC vertex loop
    }//loop over all events

    OutputFile->cd();

    SelectionTree->Write();

    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "number of events with true MC tracks : " << NumberOfSignalTruth << std::endl;
    std::cout << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;



    delete SelectionTree;

    OutputFile->Close();

    // Erase all branch addresses for the next iteration
    treenc -> ResetBranchAddresses();

    return 0;

} // end main function
