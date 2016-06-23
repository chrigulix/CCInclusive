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
double cryoradius = 191.61;
double cryoz = 1086.49 + 2*67.63;

bool inFV(double x, double y, double z);
double FlashTrackDist(double flash, double start, double end);
bool inCryostat(double x, double y, double z);

// Main function
int CCInclusiveEventSelection(std::string GeneratorName, unsigned int ThreadNumber, unsigned int NumberOfThreads)
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
        std::ifstream FileNames("/pnfs/uboone/persistent/users/aschu/devel/v05_11_01/hadd/GOODBNBBEAM/filesana.list");

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
    else if(GeneratorName == "TEM" || GeneratorName == "MEC" || GeneratorName == "MA" || GeneratorName == "BITE")
    {
        treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/"+GeneratorName+"/"+GeneratorName+"merge.root").c_str() );
    }
    else
    {
        treenc -> Add( ("/lheppc46/data/uBData/anatrees/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
    }

//     treenc -> Add( ("/lheppc46/data/uBData/anatrees/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
//     treenc -> Add( ("/media/christoph/200EFBDA63AA160B/anatrees/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
//     treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/onbeam_data_bnbSWtrigger/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
//     treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/offbeam_data_bnbSWtrigger/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
//     treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/MC_BNB_Cosmic/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );

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
    Int_t          MCTrueIndex[maxtracks];


    //define cut variables
    double flashwidth = 80; //cm. Distance flash-track
    double distcut = 5; //cm. Distance track start/end to vertex
    double lengthcut = 75; //cm. Length of longest track
    double beammin = 3.55/*-0.36*/; //us. Beam window start
    double beammax = 5.15/*-0.36*/; //us. Beam window end
    double PEthresh = 50; //PE
    double MCTrackToMCVtxDist = 0.5; //cm. distance between mc track start and mc vertex
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
            TFile* OutputFile = new TFile(("rootfiles/Hist_Track_"+TrackingName+ "_Vertex_" + VertexingName + "_"+GeneratorName+"_"+Version+FileNumberStr+"_Old.root").c_str(),"RECREATE");
            // Make a clone tree which gets filled
            TTree *SelectionTree = treenc->CloneTree(0);

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
            treenc -> SetBranchAddress("MCTruthIndex", MCTrueIndex);

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
            int MCVertexCandidate;
            int NuMuCCTrackCandidate;

            unsigned int EventsWithFlash = 0;
            unsigned int EventsVtxInFV = 0;
            unsigned int EventsTrackNearVertex = 0;
            unsigned int EventsFlashMatched = 0;
            unsigned int EventsTracksInFV = 0;
            unsigned int EventsNearVtx = 0;
            unsigned int EventsTrackLong = 0;
            unsigned int EventsTruelyReco = 0;

            unsigned int MCEventsWithFlash = 0;
            unsigned int MCEventsVtxInFV = 0;
            unsigned int MCEventsTrackNearVertex = 0;
            unsigned int MCEventsFlashMatched = 0;
            unsigned int MCEventsTracksInFV = 0;
            unsigned int MCEventsNearVtx = 0;
            unsigned int MCEventsTrackLong = 0;

            unsigned int NumberOfSignalTruth = 0;
            unsigned int NumberOfSignalTruthSel = 0;
            unsigned int NumberOfBgrNCTruthSel = 0;
            unsigned int NumberOfBgrNumuBarTruthSel = 0;
            unsigned int NumberOfBgrNueTruthSel = 0;
            unsigned int NumberOfBgrCosmicSel = 0;
            unsigned int NumberOfBgrNuOutFVSel = 0;
            unsigned int NumberOfUnknownCCBgr = 0;
            unsigned int NumberOfUnknownBgr = 0;

            TBranch* BrTrackCand = SelectionTree->Branch("TrackCand",&TrackCandidate,"TrackCand/I");
            TBranch* BrVtxCand = SelectionTree->Branch("VertexCand",&VertexCandidate,"VertexCand/I");
            TBranch* BrMCTrackCand = SelectionTree->Branch("MCTrackCand",&MCTrackCandidate,"MCTrackCand/I");
            TBranch* BrMCVtxCand = SelectionTree->Branch("MCVertexCand",&MCVertexCandidate,"MCVertexCand/I");

            double TotalPOT = 0.0;

            int Size = treenc -> GetEntries();

            // Set start and end event number for multiple threads
            unsigned long int StartEvent = Size*(ThreadNumber - 1)/NumberOfThreads;
            unsigned long int EndEvent = Size*ThreadNumber/NumberOfThreads;

            std::cout << "number of events used is: " << EndEvent-StartEvent << " of " << Size << std::endl;
            
            //Event Loop
            for(int i = StartEvent; i < EndEvent; i++)
            {
                if(i == 121558 && GeneratorName == "MEC") continue;
                if(i%1000 == 0) cout << "\t... " << i << endl;

                // Get tree entries
                treenc -> GetEntry(i);

                bool flashtag = false;
                float flashmax = 0;

                // Loop over all flashes
                for(int f = 0; f < no_flashes; f++)
                {
                    // If the flash is in the beam window and above threshold set flashtag to true
                    if( (flash_time[f] > beammin && flash_time[f] < beammax) && flash_pe[f] > PEthresh )
                    {
                        flashtag = true; //the event does have a flash inside the beam window

                        // If the new flash has more PE than the current maximum, replace the maximum
                        if(flash_pe[f] > flashmax)
                        {
                            theflash = f;
                            flashmax = flash_pe[f];
                        }
                    }
                } // flash loop

                // Preset MC candidate variables
                MCTrackCandidate = -1;
                MCVertexCandidate = -1;
                NuMuCCTrackCandidate = -1;

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
                            if( PDG_truth[track_no] == 13 
                                && MCTrueIndex[track_no] == vertex_no&& sqrt(pow(XMCTrackStart[track_no] - nuvtxx_truth[vertex_no],2) + pow(YMCTrackStart[track_no] - nuvtxy_truth[vertex_no],2) + pow(ZMCTrackStart[track_no] - nuvtxz_truth[vertex_no],2)) < MCTrackToMCVtxDist )
                            {
                                // Fill track candidate index
                                NuMuCCTrackCandidate = track_no;
                            }
                        } // MC particle loop
                    } // If numuCC in FV
                } // MC vertex loop

                // If the flash tag is ture and we have POT
                if(flashtag)
                {
                    // Prepare flags
                    bool VertexInFVFlag = true;
                    bool TrackDistanceFlag = true;
                    bool FlashMatchFlag = true;
                    bool TrackContainedFlag = false;

                    // Increase events with flash > 50 PE and within beam window
                    EventsWithFlash++;
                    if(NuMuCCTrackCandidate > -1)
                        MCEventsWithFlash++;

                    // Initialize Track Candidate properties
                    TrackCandidate = -1;
                    float TrackCandLength = 0;

                    VertexCandidate = -1;

                    // Loop over all vertices
                    for(int v = 0; v < nvtx; v++)
                    {
                        // If the vertex is contained
                        if(inFV(vtxx[v], vtxy[v], vtxz[v]))
                        {
                            // Increase count of events with a vertex in FV
                            if(VertexInFVFlag)
                            {
                                EventsVtxInFV++;
                                if(NuMuCCTrackCandidate > -1)
                                    MCEventsVtxInFV++;

                                VertexInFVFlag = false;
                            }

                            unsigned int TrackCountAtVertex = 0;

                            // Loop over all reconstructed tracks
                            for(int j = 0; j < ntracks_reco; j++)
                            {
                                // Calculate distances from track start/end to vertex and calculate track lenth
                                diststart = sqrt((vtxx[v] - trkstartx[j])*(vtxx[v] - trkstartx[j]) + (vtxy[v] - trkstarty[j])*(vtxy[v] - trkstarty[j]) + (vtxz[v] - trkstartz[j])*(vtxz[v] - trkstartz[j]));
                                distend = sqrt((vtxx[v] - trkendx[j])*(vtxx[v] - trkendx[j]) + (vtxy[v] - trkendy[j])*(vtxy[v] - trkendy[j]) + (vtxz[v] - trkendz[j])*(vtxz[v] - trkendz[j]));
                                length = sqrt(pow(trkstartx[j] - trkendx[j],2) + pow(trkstarty[j] - trkendy[j],2) + pow(trkstartz[j] - trkendz[j],2));

                                // If the track vertex distance is within cut, increase track count
                                if(diststart < distcut || distend < distcut)
                                {
                                    // Increase count of events with a distance to vertex < 5 cm
                                    if(TrackDistanceFlag)
                                    {
                                        EventsTrackNearVertex++;
                                        if(NuMuCCTrackCandidate > -1)
                                            MCEventsTrackNearVertex++;

                                        TrackDistanceFlag = false;
                                    }

                                    // Increase track at vertex count
                                    TrackCountAtVertex++;

                                    // Find the longest track from this vertex
                                    if(length > TrackCandLength)
                                    {
                                        TrackCandLength = length;
                                        TrackCandidate = j;
                                        VertexCandidate = v;
                                    }
                                } // if track vertex distance is within cut distance
                            } // reco track loop
                        } // if vertex is contained
                    } // vertex loop

                    // If the longest track length is filled
                    if(TrackCandidate > -1)
                    {
                        // If the longest track is flash matched
                        if( FlashTrackDist(flash_zcenter[theflash], trkstartz[TrackCandidate], trkendz[TrackCandidate]) < flashwidth )
                        {
                            if(FlashMatchFlag)
                            {
                                EventsFlashMatched++;
                                if(NuMuCCTrackCandidate > -1)
                                    MCEventsFlashMatched++;

                                FlashMatchFlag = false;
                                // Set track contained flag true, so the other cuts can be applied on this vertex
                                TrackContainedFlag = true;
                            }

                            // If the longest track is fully contained
                            if( inFV(trkstartx[TrackCandidate], trkstarty[TrackCandidate], trkstartz[TrackCandidate])
                                    && inFV(trkendx[TrackCandidate], trkendy[TrackCandidate], trkendz[TrackCandidate])
                                    && TrackContainedFlag )
                            {
                                EventsTracksInFV++;
                                if(NuMuCCTrackCandidate > -1)
                                    MCEventsTracksInFV++;

                                // If longest track is longer than 75 cm
                                if(TrackCandLength > lengthcut)
                                {
                                    // If track origin is neutrino
                                    if(trkorigin[TrackCandidate][trkbestplane[TrackCandidate]] == 1)
                                    {
                                        // Loop over all MCTracks
                                        for(unsigned track_no = 0; track_no < NumberOfMCTracks; track_no++)
                                        {
                                            // If MCTrackID is the same as the back tracked TruthID
                                            if(MCTrackID[track_no] == TrackIDTruth[TrackCandidate][trkbestplane[TrackCandidate]])
                                            {
                                                // Store new MCTrackCandidate and MCVertexCandidate
                                                MCTrackCandidate = track_no;
                                                MCVertexCandidate = MCTrueIndex[track_no];
                                            }
                                        } // end MCtrack loop
                                    } // if neutrino origin

                                    EventsTrackLong++;
                                    if(NuMuCCTrackCandidate > -1)
                                        MCEventsTrackLong++;

                                    // If the event is a CC interaction and the selected track is of neutrino origin
                                    if(MCVertexCandidate > -1 && ccnc_truth[MCVertexCandidate] == 0 && trkorigin[TrackCandidate][trkbestplane[TrackCandidate]] == 1)
                                    {
                                        // If there is a track candidate
                                        if(nuPDG_truth[MCVertexCandidate] == 14 && inFV(nuvtxx_truth[MCVertexCandidate],nuvtxy_truth[MCVertexCandidate],nuvtxz_truth[MCVertexCandidate]))
                                        {
                                            NumberOfSignalTruthSel++;
                                        }
                                        else if(nuPDG_truth[MCVertexCandidate] == -14) // if anti-neutrino
                                        {
                                            NumberOfBgrNumuBarTruthSel++;
                                        }
                                        else if(abs(nuPDG_truth[MCVertexCandidate]) == 12) // if electron like neutrino
                                        {
                                            NumberOfBgrNueTruthSel++;
                                        }
                                        else if(!inFV(nuvtxx_truth[MCVertexCandidate],nuvtxy_truth[MCVertexCandidate],nuvtxz_truth[MCVertexCandidate])) // if not in fiducial volume
                                        {
                                            NumberOfBgrNuOutFVSel++;
                                        }
                                        else
                                        {
                                            NumberOfUnknownCCBgr++;
                                        }
                                    } // if CC interaction
                                    else if(MCVertexCandidate > -1 && ccnc_truth[MCVertexCandidate] == 1 && trkorigin[TrackCandidate][trkbestplane[TrackCandidate]] == 1) // else if NC interaction
                                    {
                                        NumberOfBgrNCTruthSel++;
                                    }
                                    else if(trkorigin[TrackCandidate][trkbestplane[TrackCandidate]] != 1) // If selected track is not associated to a neutrino
                                    {
                                        NumberOfBgrCosmicSel++;
                                    }
                                    else
                                    {
                                        NumberOfUnknownBgr++;
                                    }

                                    SelectionTree -> Fill();

                                    double TrkStartMCStartDist = sqrt(pow(XMCTrackStart[MCTrackCandidate] - trkstartx[TrackCandidate],2) + pow(YMCTrackStart[MCTrackCandidate] - trkstarty[TrackCandidate],2) + pow(ZMCTrackStart[MCTrackCandidate] - trkstartz[TrackCandidate],2));
                                    double TrkEndMCEndDist = sqrt(pow(XMCTrackEnd[MCTrackCandidate] - trkendx[TrackCandidate],2) + pow(YMCTrackEnd[MCTrackCandidate] - trkendy[TrackCandidate],2) + pow(ZMCTrackEnd[MCTrackCandidate] - trkendz[TrackCandidate],2));
                                    double TrkStartMCEndDist = sqrt(pow(XMCTrackEnd[MCTrackCandidate] - trkstartx[TrackCandidate],2) + pow(YMCTrackEnd[MCTrackCandidate] - trkstarty[TrackCandidate],2) + pow(ZMCTrackEnd[MCTrackCandidate] - trkstartz[TrackCandidate],2));
                                    double TrkEndMCStartDist = sqrt(pow(XMCTrackStart[MCTrackCandidate] - trkendx[TrackCandidate],2) + pow(YMCTrackStart[MCTrackCandidate] - trkendy[TrackCandidate],2) + pow(ZMCTrackStart[MCTrackCandidate] - trkendz[TrackCandidate],2));

                                    // if the muon track of the NuMuCC interaction is matched correctly
                                    if( MCVertexCandidate > -1 && inFV(nuvtxx_truth[MCVertexCandidate],nuvtxy_truth[MCVertexCandidate],nuvtxz_truth[MCVertexCandidate])
                                            && PDG_truth[MCTrackCandidate] == 13 && MCTrackID[MCTrackCandidate] == TrackIDTruth[TrackCandidate][trkbestplane[TrackCandidate]] )
                                    {
                                        EventsTruelyReco++;
                                    }
                                } // if track is longer than 75 cm
                            } // if track is contained
                            // Set track contained flag false
                            TrackContainedFlag = false;
                        } // If flash matched
                    } // If there is a longest track
                } // if flashtag

                // Increase the neutrino count
                ntrue++;

                // Add POT count
                TotalPOT += potbnb;
            }//loop over all events

            OutputFile->cd();

//             SelectionTree->Print();
            SelectionTree->Write();

            std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
            std::cout << std::endl;
            std::cout << "Track Reco Product Name : " << TrackingName << "  Vertex Reco Product Name : " << VertexingName << std::endl;
            std::cout << "Total POT : " << TotalPOT*1e12 << std::endl;
            std::cout << "number of CC events with vertex in FV : " << ntrue << std::endl;
            std::cout << "number of events with flash > 50 PE : " << EventsWithFlash << " " << MCEventsWithFlash << std::endl;
            std::cout << "number of events with vtx in FV : " << EventsVtxInFV << " " << MCEventsVtxInFV << std::endl;
            std::cout << "number of events with track start/end within 5cm to vtx : " << EventsTrackNearVertex << " " << MCEventsTrackNearVertex << std::endl;
            std::cout << "number of events with tracks matched within 80cm to flash : " << EventsFlashMatched << " " << MCEventsFlashMatched << std::endl;
            std::cout << "number of events with contained tracks : " << EventsTracksInFV << " " << MCEventsTracksInFV << std::endl;
            std::cout << "number of events with longest track > 75cm : " << EventsTrackLong << " " << MCEventsTrackLong << std::endl;
            std::cout << "number of events with track start end within 5cm to mc-vtx : " << EventsTruelyReco << std::endl;
            std::cout << "number of events with contained MC tracks : " << NumberOfSignalTruth << std::endl;
            std::cout << "number of well selected events : " << NumberOfSignalTruthSel << std::endl;
            std::cout << "number of NC events selected : " << NumberOfBgrNCTruthSel << std::endl;
            std::cout << "number of anti-Neutrino events selected : " << NumberOfBgrNumuBarTruthSel << std::endl;
            std::cout << "number of Nu_e events selected : " << NumberOfBgrNueTruthSel << std::endl;
            std::cout << "number of events selected cosmic : " << NumberOfBgrCosmicSel << std::endl;
            std::cout << "number of nu_mu events out of FV : " << NumberOfBgrNuOutFVSel <<std::endl;
            std::cout << "event selection efficiency : " <<  (float)NumberOfSignalTruthSel/(float)NumberOfSignalTruth << std::endl;
//             std::cout << "event selection purity : " << (float)NumberOfSignalTruthSel/(float)(NumberOfBgrNCTruthSel+NumberOfBgrNumuBarTruthSel+NumberOfBgrNueTruthSel)
            std::cout << "event selection correctness : " <<  (float)EventsTruelyReco/(float)EventsTrackLong << std::endl;
//             std::cout << "event selection missid rate : " <<  fabs((float)EventsTruelyReco-(float)NumberOfSignalTruth)/(float)NumberOfSignalTruth << std::endl;
            std::cout << std::endl;
            std::cout << "number of unknown CC bgr : " <<  NumberOfUnknownCCBgr << std::endl;
            std::cout << "number of unknown bgr : " <<  NumberOfUnknownBgr << std::endl;
            std::cout << std::endl;
            std::cout << "--------------------------------------------------------------------------------------------" << std::endl;

            delete BrMCTrackCand;
            delete BrTrackCand;
            delete BrVtxCand;
            delete BrMCVtxCand;

            delete SelectionTree;

            OutputFile->Close();

            // Erase all branch addresses for the next iteration
            treenc -> ResetBranchAddresses();

        } // Loop over all vertexing data products
    } // Loop over all tracking data products

    return 0;

} // end main function

//This function returns if a 3D point is within the fiducial volume
bool inFV(double x, double y, double z) 
{
    if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)) return true;
    else return false;
}

//This function returns the distance between a flash and a track (in one dimension, here used only for z direction)
double FlashTrackDist(double flash, double start, double end) 
{
    if(end >= start) {
        if(flash < end && flash > start) return 0;
        else return TMath::Min(fabs(flash-start), fabs(flash-end));
    }
    else {
        if(flash > end && flash < start) return 0;
        else return TMath::Min(fabs(flash-start), fabs(flash-end));
    }
}

bool inCryostat(double x, double y, double z)
{
    // If out of cylinder axis set false
    if(z < FVz/2-cryoz/2 || z > FVz/2+cryoz/2)
    {
        return false;
    }
    else if(sqrt( pow(x-FVx/2,2) + pow(y,2)) > cryoradius)
    {
        return false;
    }
    else
    {
        return true;
    }
}
