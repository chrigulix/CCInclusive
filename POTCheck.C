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
double cryoradius = 191.61;
double cryoz = 1086.49 + 2*67.63;

//This function returns if a 3D point is within the fiducial volume
bool inFV(double x, double y, double z);
double FlashTrackDist(double flash, double start, double end);
bool inCryostat(double x, double y, double z);

// Main function
int POTCheck(std::string GeneratorName, unsigned int ThreadNumber, unsigned int NumberOfThreads)
{

    std::string Version = "v05_08_00";

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

    TChain *treenc = new TChain("analysistree/pottree");

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
    else if(GeneratorName == "prodgenie_bnb_nu_cosmic_uboone_field")
    {
//         treenc -> Add( ("/pnfs/uboone/persistent/users/aschu/offbeam_data_bnbSWtrigger/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
        std::ifstream FileNames("/pnfs/uboone/persistent/users/sowjanya/v05_08_00/bnbpluscosmics_nominal_scOn_100k/ana/filesana.list");

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
    else if(GeneratorName == "test")
    {
        for(unsigned int index = 1; index < 8; index++)
        {
            treenc -> Add( ("/uboone/data/users/aschu/FilterTestFiles/ana_hist_"+std::to_string(index)+".root").c_str() );
        }
    }
    else
    {
        treenc -> Add( ("/lheppc46/data/uBData/anatrees/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );
    }
//     treenc -> Add( ("/media/christoph/200EFBDA63AA160B/anatrees/"+GeneratorName+"_"+Version+"_anatree.root").c_str() );    

    Double_t pot;
    Double_t potbnbETOR860;
    Double_t potbnbETOR875;
    Double_t potnumiETORTGT;

    //define cut variables
    double TotPOT = 0; //cm. Distance flash-track
    double TotPOTBNB860 = 0; //cm. Distance track start/end to vertex
    double TotPOTBNB875 = 0; //cm. Length of longest track
    double TotPOTBNBTGT = 0;
    
    unsigned long int Size = treenc -> GetEntries();

    // Event Loop
    for(int i = 0; i < Size; i++)
    {   
        // Asign Branches
        treenc -> SetBranchAddress("pot", &pot);
        treenc -> SetBranchAddress("potbnbETOR860", &potbnbETOR860);
        treenc -> SetBranchAddress("potbnbETOR875", &potbnbETOR875);
        treenc -> SetBranchAddress("potnumiETORTGT", &potnumiETORTGT);
        
        if(pot > .1) std::cout << i << " " << pot << std::endl;
        
        // Calculate total POTs
        TotPOT += pot;
        TotPOTBNB860 += potbnbETOR860;
        TotPOTBNB875 += potbnbETOR875;
        TotPOTBNBTGT += potnumiETORTGT;
        
    }//loop over all events
    
    std::cout << "TotPOT " << TotPOT << std::endl;
    std::cout << "TotPOTBNB860 " << TotPOTBNB860 << std::endl;
    std::cout << "TotPOTBNB875 " << TotPOTBNB875 << std::endl;
    std::cout << "TotPOTBNBTGT " << TotPOTBNBTGT << std::endl;

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