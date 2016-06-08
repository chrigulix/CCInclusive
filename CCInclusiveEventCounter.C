#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TLine.h>
#include <TTree.h>

using namespace std;

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

//This function returns the distance between a flash and a track (in one dimension, here used only for z direction)
double FlashTrackDist(double flash, double start, double end) {
    if(end >= start) {
        if(flash < end && flash > start) return 0;
        else return TMath::Min(fabs(flash-start), fabs(flash-end));
    }
    else {
        if(flash > end && flash < start) return 0;
        else return TMath::Min(fabs(flash-start), fabs(flash-end));
    }
}

// Main function
int CCInclusiveEventCounter(std::string GeneratorName, unsigned int ThreadNumber, unsigned int NumberOfThreads)
{

    string Version = "v05_08_00";

//     string GeneratorName = "prodgenie_bnb_nu_cosmic";
//     string GeneratorName = "prodgenie_bnb_nu";
//     string GeneratorName = "prodcosmics_corsika_inTime";
//     string GeneratorName = "data_onbeam_bnb";
//     string GeneratorName = "data_offbeam_bnbext";
//     string GeneratorName = "prodgenie_bnb_nu_cosmic_uboone";
//     string GeneratorName = "prodgenie_bnb_nu_cosmic_sc_uboone";

    // Initialize and fill track reco product names
    std::vector<string> TrackProdNameVec;

//     TrackProdNameVec.push_back("pandoraNuKHit");
//     TrackProdNameVec.push_back("pandoraCosmic");
    TrackProdNameVec.push_back("pandoraNu");
//     TrackProdNameVec.push_back("pmtrack");
//     TrackProdNameVec.push_back("pandoraNuPMA");
//     TrackProdNameVec.push_back("trackkalmanhit");

    // Initialize and fill vertex reco product names
    std::vector<string> VertexProdNameVec;

//     VertexProdNameVec.push_back("nuvtx");
//     VertexProdNameVec.push_back("pandoraCosmic");
    VertexProdNameVec.push_back("pandoraNu");
//     VertexProdNameVec.push_back("pmtrack");

    std::string FileNumberStr;

    if(NumberOfThreads > 1)
    {
        FileNumberStr = "_" + std::to_string(ThreadNumber);
    }
    else
    {
        FileNumberStr = "";
    }

    std::cout << "Data Sample : " << GeneratorName << std::endl;


    TChain *treenc = new TChain("analysistree/anatree");

    if(GeneratorName == "data_onbeam_bnb")
    {
//         treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/onbeam_data_bnbSWtrigger/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
        std::ifstream FileNames("/pnfs/uboone/persistent/users/aschu/devel/v05_11_01/hadd/GOODBNB/filesana.list");

        std::string FileName;

        while(std::getline(FileNames,FileName)) 
        {
            std::cout << FileName << std::endl;
            treenc -> Add((FileName).c_str());
        }
    }
    else if(GeneratorName == "data_offbeam_bnbext")
    {
//         treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/offbeam_data_bnbSWtrigger/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
        std::ifstream FileNames("/pnfs/uboone/persistent/users/aschu/devel/v05_11_01/hadd/GOODEXTBNB/filesana.list");

        std::string FileName;

        while(std::getline(FileNames,FileName)) 
        {
            std::cout << FileName << std::endl;
            treenc -> Add((FileName).c_str());
        }
    }
    else if(GeneratorName == "prodgenie_bnb_nu_cosmic_uboone")
    {
        treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/MC_BNB_Cosmic/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
    }
    else if(GeneratorName == "TEM")
    {
        treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/"+GeneratorName+"/TEMmerge.root").c_str() );
    }
    else if(GeneratorName == "MEC")
    {
        treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/"+GeneratorName+"/MECmerge.root").c_str() );
    }
    else
    {
        treenc -> Add( ("/lheppc46/data/uBData/anatrees/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
    }

    //maximum array sizes
    const int maxentries = 35000;
    const int maxtracks = 10000;
    const int maxvtx = 500;
    const int maxnu = 10;
    const int maxmc = 10;
    const int kMaxFlashes = 5000;

    //Define variables to read from Tree
    Int_t           event;
    Int_t           run;
    Int_t           subrun;
    Int_t           triggerbits; //this is only important for data. 2048 = BNB stream
    Double_t        potbnb; //this is only important for data

    Short_t         ntracks_reco; //number of reco tracks
    Short_t         trkbestplane[maxtracks]; //plane that has most hits for a given track
    Float_t         trklen[maxtracks]; //track length
    Float_t         trkstartx[maxtracks];
    Float_t         trkstarty[maxtracks];
    Float_t         trkstartz[maxtracks];
    Float_t         trkendx[maxtracks];
    Float_t         trkendy[maxtracks];
    Float_t         trkendz[maxtracks];
    Float_t         trktheta[maxtracks];
    Float_t         trkphi[maxtracks];
    Float_t         trkmomrange[maxtracks]; //track momentum calculated from track range
    Short_t         trkId[maxtracks];
    Short_t         trkorigin[maxtracks][3]; //for MC only: which true particle contributes most hits to the reco track: 2 = cosmic. 1 = neutrino
    Int_t           TrackIDTruth[maxtracks][3]; // MC id matched with reco track
    bool            vertexatstart[maxtracks]; //for analysis: is the vertex at start of the track?
    bool            vertexatend[maxtracks]; //for analysis: ist the vertex at end of track?

    Int_t           no_flashes; //number of flashes in the event
    Float_t         flash_time[kMaxFlashes]; //flash time (in microseconds)
    Float_t         flash_pe[kMaxFlashes]; //total number of photoelectrons corresponding to the flash
    Float_t         flash_zcenter[kMaxFlashes]; //z center of flash (in cm)
    Float_t         flash_ycenter[kMaxFlashes]; //y center of flash (in cm)

    Short_t         nvtx;
    Float_t         vtxx[maxvtx];
    Float_t         vtxy[maxvtx];
    Float_t         vtxz[maxvtx];

    //finding candidate nu interaction vertex in event
    bool            candvertex[maxvtx];
    Short_t         candtrack[maxvtx];
    Float_t         candlength[maxvtx];
    Short_t         numuvertex = -1;
    Short_t         mutrack = -1;
    Float_t         mutracklength = 0;

    //MC truth
    Int_t           mcevts_truth; //neutrino interactions per event
    Float_t         nuvtxx_truth[maxmc]; //true vertex x (in cm)
    Float_t         nuvtxy_truth[maxmc];
    Float_t         nuvtxz_truth[maxmc];
    Int_t           ccnc_truth[maxmc]; //CC = 0, NC = 1
    Int_t           nuPDG_truth[maxmc]; //true neutrino pdg code. numu = 14
    Int_t           mode_truth[maxmc]; //QE = 0, RES = 1, DIS = 2
    Int_t	    PDG_truth[maxtracks];

    Int_t	   NumberOfMCTracks;

    Float_t	   XMCTrackStart[maxtracks];
    Float_t	   YMCTrackStart[maxtracks];
    Float_t	   ZMCTrackStart[maxtracks];

    Float_t	   XMCTrackEnd[maxtracks];
    Float_t	   YMCTrackEnd[maxtracks];
    Float_t	   ZMCTrackEnd[maxtracks];

    Int_t          MCTrackID[maxtracks];


    //define cut variables
    double flashwidth = 80; //cm. Distance flash-track
    double distcut = 5; //cm. Distance track start/end to vertex
    double lengthcut = 75; //cm. Length of longest track
    double beammin = 3.55/*-0.36*/; //us. Beam window start
    double beammax = 5.15/*-0.36*/; //us. Beam window end
    double PEthresh = 50; //PE
    double MCTrackToMCVtxDist = 1; //cm. distance between mc track start and mc vertex
    double TrackToMCDist = 5; //cm. Distance track start/end to mcvertex

    if(GeneratorName == "data_bnb" || GeneratorName == "data_onbeam_bnb")
    {
        std::cout << "Changed beam gate window for on-beam data" << std::endl;
        beammin = 3.3;
        beammax = 4.9;
    }
    if(GeneratorName == "data_offbeam_bnbext")
    {
        std::cout << "Changed beam gate window for off-beam data" << std::endl;
        beammin = 3.65;
        beammax = 5.25;
    }
    if(GeneratorName == "prodcosmics_corsika_inTime")
    {
        std::cout << "Changed beam gate window for Corsika sample" << std::endl;
        beammin = 3.2;
        beammax = 4.8;
    }

    // Loop over all product names
    for(const auto& TrackingName : TrackProdNameVec)
    {
        for(const auto& VertexingName : VertexProdNameVec)
        {
            // Open output file
//             TFile* OutputFile = new TFile(("rootfiles/Hist_Track_"+TrackingName+ "_Vertex_" + VertexingName + "_"+GeneratorName+"_"+Version+FileNumberStr+"_Old.root").c_str(),"RECREATE");

            treenc -> SetBranchAddress("run", &run);
            treenc -> SetBranchAddress("event", &event);
            treenc -> SetBranchAddress("potbnb", &potbnb);
            treenc -> SetBranchAddress("no_flashes", &no_flashes);
            treenc -> SetBranchAddress("flash_time", flash_time);
            treenc -> SetBranchAddress("flash_pe", flash_pe);
            treenc -> SetBranchAddress("flash_zcenter", flash_zcenter);
            treenc -> SetBranchAddress("flash_ycenter", flash_ycenter);
            treenc -> SetBranchAddress("mcevts_truth", &mcevts_truth);
            treenc -> SetBranchAddress("nuvtxx_truth", nuvtxx_truth);
            treenc -> SetBranchAddress("nuvtxy_truth", nuvtxy_truth);
            treenc -> SetBranchAddress("nuvtxz_truth", nuvtxz_truth);
            treenc -> SetBranchAddress("ccnc_truth", ccnc_truth);
            treenc -> SetBranchAddress("nuPDG_truth", nuPDG_truth);
            treenc -> SetBranchAddress("pdg", PDG_truth);
            treenc -> SetBranchAddress("mode_truth", mode_truth);
            treenc -> SetBranchAddress("geant_list_size", &NumberOfMCTracks);
            treenc -> SetBranchAddress("StartPointx", XMCTrackStart);
            treenc -> SetBranchAddress("StartPointy", YMCTrackStart);
            treenc -> SetBranchAddress("StartPointz", ZMCTrackStart);
            treenc -> SetBranchAddress("EndPointx", XMCTrackEnd);
            treenc -> SetBranchAddress("EndPointy", YMCTrackEnd);
            treenc -> SetBranchAddress("EndPointz", ZMCTrackEnd);
            treenc -> SetBranchAddress("TrackId", MCTrackID);

            // Product specific stuff
            treenc -> SetBranchAddress(("ntracks_"+TrackingName).c_str(),&ntracks_reco);
            treenc -> SetBranchAddress(("trkstartx_"+TrackingName).c_str(),trkstartx);
            treenc -> SetBranchAddress(("trkstarty_"+TrackingName).c_str(),trkstarty);
            treenc -> SetBranchAddress(("trkstartz_"+TrackingName).c_str(),trkstartz);
            treenc -> SetBranchAddress(("trkendx_"+TrackingName).c_str(),trkendx);
            treenc -> SetBranchAddress(("trkendy_"+TrackingName).c_str(),trkendy);
            treenc -> SetBranchAddress(("trkendz_"+TrackingName).c_str(),trkendz);
            treenc -> SetBranchAddress(("trktheta_"+TrackingName).c_str(),trktheta);
            treenc -> SetBranchAddress(("trkphi_"+TrackingName).c_str(),trkphi);
            treenc -> SetBranchAddress(("trkorigin_"+TrackingName).c_str(),trkorigin);
            treenc -> SetBranchAddress(("trkidtruth_"+TrackingName).c_str(),TrackIDTruth);
            treenc -> SetBranchAddress(("trkpidbestplane_"+TrackingName).c_str(), trkbestplane);

            // Program hack to apply for non uniform naming of nuvtx
            if(VertexingName != "nuvtx")
            {
                treenc -> SetBranchAddress(("nvtx_"+VertexingName).c_str(), &nvtx);
                treenc -> SetBranchAddress(("vtxx_"+VertexingName).c_str(), vtxx);
                treenc -> SetBranchAddress(("vtxy_"+VertexingName).c_str(), vtxy);
                treenc -> SetBranchAddress(("vtxz_"+VertexingName).c_str(), vtxz);
            }
            else
            {
                treenc -> SetBranchAddress("nnuvtx", &nvtx);
                treenc -> SetBranchAddress("nuvtxx", vtxx);
                treenc -> SetBranchAddress("nuvtxy", vtxy);
                treenc -> SetBranchAddress("nuvtxz", vtxz);
            }

            int theflash = -1;

            double diststart = 0;
            double distend = 0;
            double length = 0;

            int ntrue = 0;

            int TrackCandidate;
            int VertexCandidate;

            int MCTrackCandidate;
            int NuMuCCTrackCandidate;

            unsigned int BadRunCount = 0;

            double TotalPOT = 0.0;

            int Size = treenc -> GetEntries();

            // Set start and end event number for multiple threads
            unsigned long int StartEvent = Size*(ThreadNumber - 1)/NumberOfThreads;
            unsigned long int EndEvent = Size*ThreadNumber/NumberOfThreads;

            std::cout << "number of events used is: " << EndEvent-StartEvent << " of " << Size << std::endl;
            //Event Loop
            for(int i = StartEvent; i < EndEvent; i++)
            {
                if(i%1000 == 0) cout << "\t... " << i << endl;

                // Get tree entries
                treenc -> GetEntry(i);
                
                if(run < 5750 && run > 5650)
                {
                    BadRunCount++;
                }

                
            }//loop over all events

//             OutputFile->cd();

            std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
            std::cout << std::endl;
            std::cout << "Track Reco Product Name : " << TrackingName << "  Vertex Reco Product Name : " << VertexingName << std::endl;
            std::cout << "Bad purity runs : " << BadRunCount << std::endl;

            std::cout << std::endl;
            std::cout << "--------------------------------------------------------------------------------------------" << std::endl;

//             OutputFile->Close();

            // Erase all branch addresses for the next iteration
            treenc -> ResetBranchAddresses();

        } // Loop over all vertexing data products
    } // Loop over all tracking data products

    return 0;

} // end main function
