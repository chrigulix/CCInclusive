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

// Main function
void StoreCosmicDist()
{
    // 
    std::string InputFolder = "/home/christoph/anatrees/CCInclusiveNote";
    std::string OutputFolder = ".";
    
    // Reconstructed data product used
    std::string ProductName = "pandoraNu";

    // Data input file vector
    std::vector<TChain*> ChainVec;
    
    // Prepare track length -> momentum function
    MomentumSplinePreparation();

    std::vector<unsigned int> NumberOfEvents;

    // sizes
    const unsigned int maxtracks = 10000;
    const unsigned int maxvtx = 500;

    // Readout variables
    short ntracks_reco; //number of reco tracks
    short trkbestplane[maxtracks]; //plane that has most hits for a given track
    float trklen[maxtracks]; //track length
    float trkstartx[maxtracks];
    float trkstarty[maxtracks];
    float trkstartz[maxtracks];
    float trkendx[maxtracks];
    float trkendy[maxtracks];
    float trkendz[maxtracks];
    float trktheta[maxtracks];
    float trkphi[maxtracks];
    float trkmomrange[maxtracks]; //track momentum calculated from track range
    short trkId[maxtracks];
    short trkorigin[maxtracks][3]; //for MC only: which true particle contributes most hits to the reco track: 2 = cosmic. 1 = neutrino
    int TrackIDTruth[maxtracks][3]; // MC id matched with reco track
    bool vertexatstart[maxtracks]; //for analysis: is the vertex at start of the track?
    bool vertexatend[maxtracks]; //for analysis: ist the vertex at end of track?
    short nvtx;
    float vtxx[maxvtx];
    float vtxy[maxvtx];
    float vtxz[maxvtx];

    // Data containers
    std::vector<TH1F*> SelectionTrackRange;
    std::vector<TH1F*> SelectionCosTheta;
    std::vector<TH1F*> SelectionTheta;
    std::vector<TH1F*> SelectionPhi;
    std::vector<TH1F*> SelectionMomentum;
    std::vector<TH1F*> SelectionTrackLength;
    std::vector<TH1F*> SelXVtxPosition;
    std::vector<TH1F*> SelYVtxPosition;
    std::vector<TH1F*> SelZVtxPosition;

    // Number of bins
    size_t NumberOfBins = 20;

    // Data set label
    std::vector<std::string> GenLabel;
    std::vector<std::pair<double,double>> BeamWindow;

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add((InputFolder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_1_Mod.root").c_str());
    ChainVec.back() -> Add((InputFolder+"/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_2_Mod.root").c_str());
    NumberOfEvents.push_back(ChainVec.back() -> GetEntries());
    GenLabel.push_back("Data Off-Beam BNBEXT All");
    BeamWindow.push_back(std::make_pair(3.3,4.9));

    ChainVec.push_back(new TChain("analysistree/anatree"));
    ChainVec.back() -> Add((InputFolder+"/prodcosmics_corsika_inTime_v05_08_00_anatree.root").c_str());
    NumberOfEvents.push_back(ChainVec.back() -> GetEntries());
    GenLabel.push_back("In Time Corsika All");
    BeamWindow.push_back(std::make_pair(3.2,4.8));
    
    GenLabel.push_back("Cosmic Systematics All");

    for(const auto& Label : GenLabel)
    {
        SelectionTrackRange.push_back(new TH1F(("Track Range "+Label).c_str(),"Muon Track Range",NumberOfBins,0,700));
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

        SelectionTrackLength.push_back(new TH1F(("Track Length "+Label).c_str(),"Candidate Track Length",NumberOfBins,0,800));
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
    }
    
    std::vector<unsigned long int> NumberOfTracks;
    std::vector<unsigned long int> NumberOfVertices;

    // Loop over all files
    for(unsigned int file_no = 0; file_no < ChainVec.size(); file_no++)
    {
        // Print file information
        std::cout << "-------------File Progress--------------" << std::endl;
        std::cout << "File \t \t" << file_no+1 << " of " << ChainVec.size() << std::endl;
        std::cout << "----------------------------------------" << std::endl;

        // Assign branch addresses to data
        ChainVec.at(file_no) -> SetBranchAddress(("ntracks_"+ProductName).c_str(),&ntracks_reco);
        ChainVec.at(file_no) -> SetBranchAddress(("trklen_"+ProductName).c_str(), trklen);
        ChainVec.at(file_no) -> SetBranchAddress(("trkstartx_"+ProductName).c_str(),trkstartx);
        ChainVec.at(file_no) -> SetBranchAddress(("trkstarty_"+ProductName).c_str(),trkstarty);
        ChainVec.at(file_no) -> SetBranchAddress(("trkstartz_"+ProductName).c_str(),trkstartz);
        ChainVec.at(file_no) -> SetBranchAddress(("trkendx_"+ProductName).c_str(),trkendx);
        ChainVec.at(file_no) -> SetBranchAddress(("trkendy_"+ProductName).c_str(),trkendy);
        ChainVec.at(file_no) -> SetBranchAddress(("trkendz_"+ProductName).c_str(),trkendz);
        ChainVec.at(file_no) -> SetBranchAddress(("trktheta_"+ProductName).c_str(),trktheta);
        ChainVec.at(file_no) -> SetBranchAddress(("trkphi_"+ProductName).c_str(),trkphi);
        ChainVec.at(file_no) -> SetBranchAddress(("nvtx_"+ProductName).c_str(), &nvtx);
        ChainVec.at(file_no) -> SetBranchAddress(("vtxx_"+ProductName).c_str(), vtxx);
        ChainVec.at(file_no) -> SetBranchAddress(("vtxy_"+ProductName).c_str(), vtxy);
        ChainVec.at(file_no) -> SetBranchAddress(("vtxz_"+ProductName).c_str(), vtxz);
        
        NumberOfTracks.push_back(0);
        NumberOfVertices.push_back(0);

        // Loop over events
        for(unsigned int tree_index = 0; tree_index < ChainVec.at(file_no) -> GetEntries(); tree_index++)
        {
            // Progress indicator
            if(!(tree_index % 1000)) std::cout << "Event\t" << tree_index << "\t of \t" << ChainVec.at(file_no) -> GetEntries() << std::endl;

            // Get tree entry for this event
            ChainVec.at(file_no) -> GetEntry(tree_index);

            // Loop over all tracks in an event
            for(unsigned int track_no = 0; track_no < ntracks_reco; track_no++)
            {
                SelectionTrackRange.at(file_no) -> Fill(CalcRange(trkstartx[track_no],trkstarty[track_no],trkstartz[track_no],trkendx[track_no],trkendy[track_no],trkendz[track_no]));
                SelectionCosTheta.at(file_no) -> Fill(std::cos(trktheta[track_no]));
                SelectionTheta.at(file_no) -> Fill(trktheta[track_no]/Pi*180);
                SelectionPhi.at(file_no) -> Fill(trkphi[track_no]/Pi*180);
                SelectionMomentum.at(file_no) -> Fill(GetMomentum(trklen[track_no]));
                SelectionTrackLength.at(file_no) -> Fill(trklen[track_no]);
            } // end track loop

            // Loop over all vertices
            for(unsigned int vtx_no = 0; vtx_no < nvtx; vtx_no++)
            {
                SelXVtxPosition.at(file_no) -> Fill(vtxx[vtx_no]);
                SelYVtxPosition.at(file_no) -> Fill(vtxy[vtx_no]);
                SelZVtxPosition.at(file_no) -> Fill(vtxz[vtx_no]);
            }
            
            // Fill vertex and track counters
            NumberOfTracks.back() += ntracks_reco;
            NumberOfVertices.back() += nvtx;
        } // end event loop

        // Scale histograms to number of events
        SelectionTrackRange.at(file_no) -> Scale(1/(double)ChainVec.at(file_no)->GetEntries());
        SelectionCosTheta.at(file_no) -> Scale(1/(double)ChainVec.at(file_no)->GetEntries());
        SelectionTheta.at(file_no) -> Scale(1/(double)ChainVec.at(file_no)->GetEntries());
        SelectionPhi.at(file_no) -> Scale(1/(double)ChainVec.at(file_no)->GetEntries());
        SelectionMomentum.at(file_no) -> Scale(1/(double)ChainVec.at(file_no)->GetEntries());
        SelectionTrackLength.at(file_no) -> Scale(1/(double)ChainVec.at(file_no)->GetEntries());
        SelXVtxPosition.at(file_no) -> Scale(1/(double)ChainVec.at(file_no)->GetEntries());
        SelYVtxPosition.at(file_no) -> Scale(1/(double)ChainVec.at(file_no)->GetEntries());
        SelZVtxPosition.at(file_no) -> Scale(1/(double)ChainVec.at(file_no)->GetEntries());
    } // end files loop
    
    std::cout << "-------Number of Tracks per Event-------" << std::endl;
    std::cout << "Off-Beam BNBEXT \t" << (double)NumberOfTracks.at(0)/ChainVec.at(0)->GetEntries()  <<  std::endl;
    std::cout << "InTime Corsika  \t" << (double)NumberOfTracks.at(1)/ChainVec.at(1)->GetEntries() << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "------Number of Vertices per Event------" << std::endl;
    std::cout << "Off-Beam BNBEXT \t" << (double)NumberOfTracks.at(0)/(double)ChainVec.at(0)->GetEntries()  <<  std::endl;
    std::cout << "InTime Corsika  \t" << (double)NumberOfTracks.at(1)/(double)ChainVec.at(1)->GetEntries() << std::endl;    
    std::cout << "----------------------------------------" << std::endl;
    
    // Fill the relative variance into the bins of the last histogram
    for(unsigned int bin_no = 1; bin_no <= NumberOfBins; bin_no++)
    {
        SelectionTrackRange.back()->SetBinContent(bin_no,std::pow(SelectionTrackRange.at(0)->GetBinContent(bin_no)/SelectionTrackRange.at(1)->GetBinContent(bin_no) - 1,2));
        SelectionCosTheta.back()->SetBinContent(bin_no,std::pow(SelectionCosTheta.at(0)->GetBinContent(bin_no)/SelectionCosTheta.at(1)->GetBinContent(bin_no) - 1,2));
        SelectionTheta.back()->SetBinContent(bin_no,std::pow(SelectionTheta.at(0)->GetBinContent(bin_no)/SelectionTheta.at(1)->GetBinContent(bin_no) - 1,2));
        SelectionPhi.back()->SetBinContent(bin_no,std::pow(SelectionPhi.at(0)->GetBinContent(bin_no)/SelectionPhi.at(1)->GetBinContent(bin_no) - 1,2));
        
        // Avoid 0 bins
        if(SelectionMomentum.at(1)->GetBinContent(bin_no) > 0.0)
        {
            SelectionMomentum.back()->SetBinContent(bin_no,std::pow(SelectionMomentum.at(0)->GetBinContent(bin_no)/SelectionMomentum.at(1)->GetBinContent(bin_no) - 1,2));
        }
        else
        {
            SelectionMomentum.back()->SetBinContent(bin_no,1.0);
        }
        if(SelectionTrackLength.at(1)->GetBinContent(bin_no) > 0.0)
        {
            SelectionTrackLength.back()->SetBinContent(bin_no,std::pow(SelectionTrackLength.at(0)->GetBinContent(bin_no)/SelectionTrackLength.at(1)->GetBinContent(bin_no) - 1,2));
        }
        else
        {
            SelectionTrackLength.back()->SetBinContent(bin_no,1.0);
        }
        
        SelXVtxPosition.back()->SetBinContent(bin_no,std::pow(SelXVtxPosition.at(0)->GetBinContent(bin_no)/SelXVtxPosition.at(1)->GetBinContent(bin_no) - 1,2));
        SelYVtxPosition.back()->SetBinContent(bin_no,std::pow(SelYVtxPosition.at(0)->GetBinContent(bin_no)/SelYVtxPosition.at(1)->GetBinContent(bin_no) - 1,2));
        SelZVtxPosition.back()->SetBinContent(bin_no,std::pow(SelZVtxPosition.at(0)->GetBinContent(bin_no)/SelZVtxPosition.at(1)->GetBinContent(bin_no) - 1,2));
    }
    
    // Open output file
    TFile* OutputFile = new TFile((OutputFolder+"/Cosmic_Distributions_Histograms_Mod.root").c_str(),"RECREATE");

    // switch to output file
    OutputFile->cd();
    
    // Loop over all labels/files
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
    
    // Properly close file
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

    InputFile.close();

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
