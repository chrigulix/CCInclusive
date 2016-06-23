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
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TAxis.h>
#include <TSpline.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TGraphAsymmErrors.h>

//This defines our current settings for the fiducial volume
double FVx = 256.35;
double FVy = 233;
double FVz = 1036.8;
double borderx = 10.;
double bordery = 20.;
double borderz = 10.;

float GetMaximum ( const std::vector<TH1F*>& HistVector );
void AddFirstTwoHistograms ( std::vector<TH1F*>& HistVector, float Weight );
void AddFirstTwoHistograms2D ( std::vector<TH2F*>& HistVector, float Weight );
float CalcLength ( const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2 );
double FlashTrackDist ( double flash, double start, double end );
bool inDeadRegion ( double y, double z );
std::vector<TSpline5> Systematics();
void AdjustSysError ( std::vector<TH1F*>& HistVector );
bool inFV ( double x, double y, double z );

void HistoProducerMC()
{
    TGaxis::SetMaxDigits ( 4 );


    //     std::string TrackProdName="pandoraNuKHit";
//     std::string TrackProdName = "pandoraCosmic";
    std::string TrackProdName="pandoraNu";
//     std::string TrackProdName="pmtrack";
//     std::string TrackProdName="pandoraNuPMA";
//     std::string TrackProdName="trackkalmanhit";

//     std::string  VertexProdName="nuvtx";
//     std::string VertexProdName="pandoraCosmic";
    std::string VertexProdName = "pandoraNu";
//     std::string VertexProdName = "pmtrack";

//     std::string SelectionLabel = "_Old";
    std::string SelectionLabel = "_Mod";
//     std::string SelectionLabel = "_New";

//     std::string FileType = "png";
    std::string FileType = "pdf";

    std::vector<TChain*> ChainVec;

    std::vector<std::string> EfficiencyLabel;
    std::vector<std::string> MCLabel;
    std::vector<std::string> GenLabel;
    std::vector<std::string> IntLabel;
    std::vector<std::string> EffIntLabel;
    std::vector<std::string> IntTrueLabel;

    std::vector<float> ScalingFactors;
    ScalingFactors.push_back ( 1 );
    ScalingFactors.push_back ( 1 );
    ScalingFactors.push_back ( 1 );
    ScalingFactors.push_back ( 1 );

    // Binning
    unsigned int NumberOfBins = 20;
    unsigned int NumberOf2DBins = 10;

    std::vector<TH1F*> SelectionTrackRange;
    std::vector<TH1F*> SelectionEnergy;
    std::vector<TH1F*> SelectionMomentum;
    std::vector<TH1F*> SelectionTheta;
    std::vector<TH1F*> SelectionCosTheta;
    std::vector<TH1F*> SelectionPhi;

    std::vector<TH1F*> SelXTrackStartEnd;
    std::vector<TH1F*> SelYTrackStartEnd;
    std::vector<TH1F*> SelZTrackStartEnd;

    std::vector<TH1F*> SelXVtxPosition;
    std::vector<TH1F*> SelYVtxPosition;
    std::vector<TH1F*> SelZVtxPosition;

    std::vector<TEfficiency*> EffTrackRange;
    std::vector<TEfficiency*> EffEnergy;
    std::vector<TEfficiency*> EffMomentum;
    std::vector<TEfficiency*> EffTheta;
    std::vector<TEfficiency*> EffCosTheta;
    std::vector<TEfficiency*> EffPhi;

    std::vector<TEfficiency*> EffXTrackStartEnd;
    std::vector<TEfficiency*> EffYTrackStartEnd;
    std::vector<TEfficiency*> EffZTrackStartEnd;

    std::vector<TEfficiency*> EffXVtxPosition;
    std::vector<TEfficiency*> EffYVtxPosition;
    std::vector<TEfficiency*> EffZVtxPosition;

    TF1* SinTheta = new TF1 ( "const","sin(x)",0,3.142 );
    
    TPaveText TextSimulation(0.5,0.92,0.9,0.96,"nbNDC");
    TextSimulation.AddText("MicroBooNE Simulation, Preliminary");
    TextSimulation.SetTextSize(0.04);
    TextSimulation.SetTextColor(12);
    TextSimulation.SetLineColorAlpha(0,0);
    TextSimulation.SetFillColorAlpha(0,0);
    TextSimulation.SetTextAlign(33);
    
    TPaveText TextPreliminary(0.6,0.92,0.9,0.96,"nbNDC");
    TextPreliminary.AddText("MicroBooNE Preliminary");
    TextPreliminary.SetTextSize(0.04);
    TextPreliminary.SetTextColor(12);
    TextPreliminary.SetLineColorAlpha(0,0);
    TextPreliminary.SetFillColorAlpha(0,0);
    TextPreliminary.SetTextAlign(33);
    
    TPaveText TextSelection(0.1,0.92,0.3,0.96,"nbNDC");
    TextSelection.AddText("Selection I");
    TextSelection.SetTextSize(0.04);
    TextSelection.SetTextColor(1);
    TextSelection.SetLineColorAlpha(0,0);
    TextSelection.SetFillColorAlpha(0,0);
    TextSelection.SetTextAlign(13);

    TLegend* LegendEfficiency = new TLegend ( 0.15,0.7,0.45,0.85 );
    LegendEfficiency->SetLineColorAlpha ( 0,0 );
    LegendEfficiency->SetLineStyle ( 0 );
    LegendEfficiency->SetFillStyle ( 0 );
    LegendEfficiency->SetMargin ( 0.2 );
    LegendEfficiency->SetTextFont ( 43 );
    LegendEfficiency->SetTextSize ( 35 );

    EfficiencyLabel.push_back ( "Vertex in FV Eff." );
    EfficiencyLabel.push_back ( "Contained Track Eff." );
    
    TLegend* LegendMCSel = new TLegend ( 0.6,0.7,0.8,0.75 );
    LegendMCSel->SetLineColorAlpha ( 0,0 );
    LegendMCSel->SetLineStyle ( 0 );
    LegendMCSel->SetFillStyle ( 0 );
    LegendMCSel->SetMargin ( 0.2 );
    LegendMCSel->SetTextFont ( 43 );
    LegendMCSel->SetTextSize ( 35 );

    TLegend* LegendMC = new TLegend ( 0.5,0.6,0.8,0.8 );
    LegendMC->SetLineColorAlpha ( 0,0 );
    LegendMC->SetLineStyle ( 0 );
    LegendMC->SetFillStyle ( 0 );
    LegendMC->SetMargin ( 0.2 );
    LegendMC->SetTextFont ( 43 );
    LegendMC->SetTextSize ( 35 );

    MCLabel.push_back ( "MC True Vertex in FV" );
    MCLabel.push_back ( "MC True Contained Tracks" );
    MCLabel.push_back ( "Selection on MC BNB+Cosmic with Stat. Error" );

    TLegend* LegendInt = new TLegend ( 0.15,0.65,0.55,0.85 );
    LegendInt->SetLineColorAlpha ( 0,0 );
    LegendInt->SetLineStyle ( 0 );
    LegendInt->SetFillStyle ( 0 );
    LegendInt->SetMargin ( 0.2 );
    LegendInt->SetTextFont ( 43 );
    LegendInt->SetTextSize ( 35 );

    IntLabel.push_back ( "CCQE events after selection" );
    IntLabel.push_back ( "CCRES events after selection" );
    IntLabel.push_back ( "CCDIS events after selection" );

    TLegend* LegendIntTrue = new TLegend ( 0.15,0.65,0.55,0.85 );
    LegendIntTrue->SetLineColorAlpha ( 0,0 );
    LegendIntTrue->SetLineStyle ( 0 );
    LegendIntTrue->SetFillStyle ( 0 );
    LegendIntTrue->SetMargin ( 0.2 );
    LegendIntTrue->SetTextFont ( 43 );
    LegendIntTrue->SetTextSize ( 35 );

    IntTrueLabel.push_back ( "CCQE events before selection" );
    IntTrueLabel.push_back ( "CCRES events before selection" );
    IntTrueLabel.push_back ( "CCDIS events before selection" );

    TLegend* LegendEffInt = new TLegend ( 0.55,0.65,0.75,0.85 );
    LegendEffInt->SetLineColorAlpha ( 0,0 );
    LegendEffInt->SetLineStyle ( 0 );
    LegendEffInt->SetFillStyle ( 0 );
    LegendEffInt->SetMargin ( 0.2 );
    LegendEffInt->SetTextFont ( 43 );
    LegendEffInt->SetTextSize ( 35 );

    EffIntLabel.push_back ( "CCQE Efficiency" );
    EffIntLabel.push_back ( "CCRES Efficiency" );
    EffIntLabel.push_back ( "CCDIS Efficiency" );
//     EffIntLabel.push_back ( "QE Contained Eff." );
//     EffIntLabel.push_back ( "RES Contained Eff." );
//     EffIntLabel.push_back ( "DIS Contained Eff." );

//     MCLabel.push_back ( "Selection MC BNB+Cosmic Sys. Error" );

    TLegend* FlashLabel = new TLegend ( 0.7,0.7,0.9,0.9 );

    GenLabel.push_back ( "Truth BNB Nu Cosmic in FV" );
    GenLabel.push_back ( "Truth BNB Nu Cosmic in FV Muon contained " );
    GenLabel.push_back ( "Selection BNB Nu Cosmic" );
    GenLabel.push_back ( "Contained Selection BNB Nu Cosmic" );
    GenLabel.push_back ( "QE Truth BNB Nu Cosmic in FV" );
    GenLabel.push_back ( "RES Truth BNB Nu Cosmic in FV" );
    GenLabel.push_back ( "DIS Truth BNB Nu Cosmic in FV" );
    GenLabel.push_back ( "QE Truth BNB Nu Cosmic Contained" );
    GenLabel.push_back ( "RES Truth BNB Nu Cosmic Contained" );
    GenLabel.push_back ( "DIS Truth BNB Nu Cosmic Contained" );
    GenLabel.push_back ( "QE Selection BNB Nu Cosmic in FV" );
    GenLabel.push_back ( "RES Selection BNB Nu Cosmic in FV" );
    GenLabel.push_back ( "DIS Selection BNB Nu Cosmic in FV" );
    GenLabel.push_back ( "MC Systematic Errors" );

    std::vector<TSpline5> SystematicErrors = Systematics();

    std::vector<unsigned int> ColorMap = {30,38,42};

    ChainVec.push_back ( new TChain ( "anatree" ) );
    ChainVec.back() -> Add ( "/lheppc46/data/uBData/anatrees/Hist_MC_Truth_prodgenie_bnb_nu_cosmic_uboone_v05_08_00.root" );

    ChainVec.push_back ( new TChain ( "anatree" ) );
    ChainVec.back() -> Add ( ( "/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_prodgenie_bnb_nu_cosmic_uboone_v05_08_00"+ SelectionLabel +".root" ).c_str() );

    for ( const auto& Label : GenLabel )
    {
        SelectionTrackRange.push_back ( new TH1F ( ( "Track Range"+Label ).c_str(),"Track Range of Selected Track",NumberOfBins,0,1000 ) );
        SelectionTrackRange.back()->SetStats ( 0 );
        SelectionTrackRange.back()->GetXaxis()->SetTitle ( "Track Range [cm]" );
        SelectionTrackRange.back()->GetYaxis()->SetTitle ( "No. of events" );

        SelectionTheta.push_back ( new TH1F ( ( "#theta-Angle"+Label ).c_str(),"#theta-Angle of Selected Track",NumberOfBins,0,3.142 ) );
        SelectionTheta.back()->SetStats ( 0 );
        SelectionTheta.back()->GetXaxis()->SetTitle ( "#thetaangle [rad]" );
        SelectionTheta.back()->GetYaxis()->SetTitle ( "No. of events" );

        SelectionCosTheta.push_back ( new TH1F ( ( "cos#theta-Angle"+Label ).c_str(),"cos#theta of Selected Track",NumberOfBins,-1,1 ) );
        SelectionCosTheta.back()->SetStats ( 0 );
        SelectionCosTheta.back()->GetXaxis()->SetTitle ( "cos(#theta)" );
        SelectionCosTheta.back()->GetYaxis()->SetTitle ( "No. of events" );

        SelectionPhi.push_back ( new TH1F ( ( "#phi-Angle"+Label ).c_str(),"#phi-Angle of Selected Track",NumberOfBins,-3.142,3.142 ) );
        SelectionPhi.back()->SetStats ( 0 );
        SelectionPhi.back()->GetXaxis()->SetTitle ( "#phi angle [rad]" );
        SelectionPhi.back()->GetYaxis()->SetTitle ( "No. of events" );

        SelectionEnergy.push_back ( new TH1F ( ( "Energy"+Label ).c_str(),"Energy of Selected Track",NumberOfBins,0,3 ) );
        SelectionEnergy.back()->SetStats ( 0 );
        SelectionEnergy.back()->GetXaxis()->SetTitle ( "Muon Kinetic Energy [MeV]" );
        SelectionEnergy.back()->GetYaxis()->SetTitle ( "No. of events" );

        SelectionMomentum.push_back ( new TH1F ( ( "Momentum"+Label ).c_str(),"Momentum of Selected Track",NumberOfBins,0,3 ) );
        SelectionMomentum.back()->SetStats ( 0 );
        SelectionMomentum.back()->GetXaxis()->SetTitle ( "Muon Momentum [GeV/c]" );
        SelectionMomentum.back()->GetYaxis()->SetTitle ( "No. of events" );

        SelXTrackStartEnd.push_back ( new TH1F ( ( "XTrack"+Label ).c_str(),"X Track Start & End Positions",NumberOfBins,0,256 ) );
        SelXTrackStartEnd.back()->SetStats ( 0 );
        SelXTrackStartEnd.back()->GetXaxis()->SetTitle ( "Track start and end x [cm]" );
        SelXTrackStartEnd.back()->GetYaxis()->SetTitle ( "No. of events" );

        SelYTrackStartEnd.push_back ( new TH1F ( ( "YTrack"+Label ).c_str(),"Y Track Start & End Positions",NumberOfBins,-233/2,233/2 ) );
        SelYTrackStartEnd.back()->SetStats ( 0 );
        SelYTrackStartEnd.back()->GetXaxis()->SetTitle ( "Track start and end y [cm]" );
        SelYTrackStartEnd.back()->GetYaxis()->SetTitle ( "No. of events" );

        SelZTrackStartEnd.push_back ( new TH1F ( ( "ZTrack"+Label ).c_str(),"Z Track Start & End Positions",NumberOfBins,0,1036.8 ) );
        SelZTrackStartEnd.back()->SetStats ( 0 );
        SelZTrackStartEnd.back()->GetXaxis()->SetTitle ( "Track start and end z [cm]" );
        SelZTrackStartEnd.back()->GetYaxis()->SetTitle ( "No. of events" );

        SelXVtxPosition.push_back ( new TH1F ( ( "XVertex"+Label ).c_str(),"X Vertex Position",NumberOfBins,0,256 ) );
        SelXVtxPosition.back()->SetStats ( 0 );
        SelXVtxPosition.back()->GetXaxis()->SetTitle ( "Vertex x [cm]" );
        SelXVtxPosition.back()->GetYaxis()->SetTitle ( "No. of events" );

        SelYVtxPosition.push_back ( new TH1F ( ( "YVertex"+Label ).c_str(),"Y Vertex Position",NumberOfBins,-233/2,233/2 ) );
        SelYVtxPosition.back()->SetStats ( 0 );
        SelYVtxPosition.back()->GetXaxis()->SetTitle ( "Vertex y [cm]" );
        SelYVtxPosition.back()->GetYaxis()->SetTitle ( "No. of events" );

        SelZVtxPosition.push_back ( new TH1F ( ( "ZVertex"+Label ).c_str(),"Z Vertex Position",NumberOfBins,0,1036.8 ) );
        SelZVtxPosition.back()->SetStats ( 0 );
        SelZVtxPosition.back()->GetXaxis()->SetTitle ( "Vertex z [cm]" );
        SelZVtxPosition.back()->GetYaxis()->SetTitle ( "No. of events" );
    }

    int Run;
    int Subrun;
    int Event;

    int TrkID;
    int VtxID;

    float FlashPE[5000];
    int NumberOfFlashes;
    float FlashTime[5000];
    float ZFlashCenter[5000];

    int MCTrkID;
    int MCVtxID;
    int CCNCFlag[10];
    int TruthMode[10];
    int NuPDGTruth[10];
    int PDGTruth[5000];
    float NuEnergyTruth[10];
    float TrueLeptonMomentum[10];

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

    float KineticEnergy[5000][3];

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

    double beammin;
    double beammax;

    for ( unsigned int file_no = 0; file_no < ChainVec.size(); file_no++ )
    {
        beammin = 3.55;
        beammax = 5.15;

        ChainVec.at ( file_no ) -> SetBranchAddress ( "run", &Run );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "subrun", &Subrun );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "event", &Event );

        ChainVec.at ( file_no ) -> SetBranchAddress ( "TrackCand", &TrkID );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "VertexCand", &VtxID );

        ChainVec.at ( file_no ) -> SetBranchAddress ( "flash_pe", FlashPE );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "no_flashes", &NumberOfFlashes );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "flash_time", FlashTime );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "flash_zcenter", ZFlashCenter );

        ChainVec.at ( file_no ) -> SetBranchAddress ( ( "trkorigin_"+TrackProdName ).c_str(), TrkOrigin );
        ChainVec.at ( file_no ) -> SetBranchAddress ( ( "trkpidbestplane_"+TrackProdName ).c_str(), TrkBestPlane );

        ChainVec.at ( file_no ) -> SetBranchAddress ( ( "trkke_"+TrackProdName ).c_str(), KineticEnergy );
        ChainVec.at ( file_no ) -> SetBranchAddress ( ( "trkmomrange_"+TrackProdName ).c_str(), TrackMomentum );
        ChainVec.at ( file_no ) -> SetBranchAddress ( ( "trktheta_"+TrackProdName ).c_str(), TrackTheta );
        ChainVec.at ( file_no ) -> SetBranchAddress ( ( "trkphi_"+TrackProdName ).c_str(),TrackPhi );

        ChainVec.at ( file_no ) -> SetBranchAddress ( ( "trkstartx_"+TrackProdName ).c_str(),XTrackStart );
        ChainVec.at ( file_no ) -> SetBranchAddress ( ( "trkstarty_"+TrackProdName ).c_str(),YTrackStart );
        ChainVec.at ( file_no ) -> SetBranchAddress ( ( "trkstartz_"+TrackProdName ).c_str(),ZTrackStart );

        ChainVec.at ( file_no ) -> SetBranchAddress ( ( "trkendx_"+TrackProdName ).c_str(),XTrackEnd );
        ChainVec.at ( file_no ) -> SetBranchAddress ( ( "trkendy_"+TrackProdName ).c_str(),YTrackEnd );
        ChainVec.at ( file_no ) -> SetBranchAddress ( ( "trkendz_"+TrackProdName ).c_str(),ZTrackEnd );

        if ( VertexProdName != "nuvtx" )
        {
            ChainVec.at ( file_no ) -> SetBranchAddress ( ( "vtxx_"+VertexProdName ).c_str(), XVertexPosition );
            ChainVec.at ( file_no ) -> SetBranchAddress ( ( "vtxy_"+VertexProdName ).c_str(), YVertexPosition );
            ChainVec.at ( file_no ) -> SetBranchAddress ( ( "vtxz_"+VertexProdName ).c_str(), ZVertexPosition );
        }
        else
        {
            ChainVec.at ( file_no ) -> SetBranchAddress ( "nuvtxx", XVertexPosition );
            ChainVec.at ( file_no ) -> SetBranchAddress ( "nuvtxy", YVertexPosition );
            ChainVec.at ( file_no ) -> SetBranchAddress ( "nuvtxz", ZVertexPosition );
        }

        ChainVec.at ( file_no ) -> SetBranchAddress ( "MCTrackCand", &MCTrkID );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "MCVertexCand", &MCVtxID );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "ccnc_truth", CCNCFlag );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "mode_truth", TruthMode );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "nuPDG_truth", NuPDGTruth );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "pdg", PDGTruth );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "enu_truth", NuEnergyTruth );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "lep_mom_truth", TrueLeptonMomentum );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "mcevts_truth", &mcevts_truth );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "nuvtxx_truth", XnuVtxTruth );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "nuvtxy_truth", YnuVtxTruth );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "nuvtxz_truth", ZnuVtxTruth );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "nuPDG_truth", nuPDGTruth );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "geant_list_size", &NumberOfMCTracks );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "StartPointx", XMCTrackStart );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "StartPointy", YMCTrackStart );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "StartPointz", ZMCTrackStart );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "EndPointx", XMCTrackEnd );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "EndPointy", YMCTrackEnd );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "EndPointz", ZMCTrackEnd );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "theta", MCTheta );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "enu_truth", NuEnergyTruth );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "phi", MCPhi );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "Eng", MCEnergy );

        unsigned int nubar = 0;
        unsigned int nue = 0;
        unsigned int NCnu = 0;
        unsigned int Cosmic = 0;
        unsigned int UnknownOrigin = 0;
        unsigned int Signal = 0;

        unsigned int nuQE = 0;
        unsigned int nuRES = 0;
        unsigned int nuDIS = 0;
        unsigned int nuCOH = 0;

        unsigned int negPhi = 0;
        unsigned int posPhi = 0;

        float XFVCutValue = 10; //10
        float YFVCutValue = 20; //20
        float ZFVCutValue = 10; //10
        float FlashTrackCut = 80; //80

        unsigned int ContainedTracks = 0;

        for ( unsigned int tree_index = 0; tree_index < ChainVec.at ( file_no )->GetEntries(); tree_index++ )
        {
            if ( ! ( tree_index % 1000 ) )
            {
                std::cout << "Event\t" << tree_index << "\t of \t" << ChainVec.at ( file_no )->GetEntries() << std::endl;
            }

            ChainVec.at ( file_no )->GetEntry ( tree_index );

//             if(tree_index < 100)
//             {
//                 std::cout << tree_index << " " << MCTrkID << " " << inFV( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID] ) <<  " " << XMCTrackStart[MCTrkID] << " " << YMCTrackStart[MCTrkID] << " " << ZMCTrackStart[MCTrkID] << std::endl;
//             }

//             std::cout << MCTrkID << " " << inFV( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID] ) << " " << XMCTrackStart[MCTrkID] << " " << YMCTrackStart[MCTrkID] << " " << ZMCTrackStart[MCTrkID] << std::endl;
//             if(MCTrkID > NumberOfMCTracks) std::cout << MCTrkID << " " << NumberOfMCTracks << std::endl;

            if ( file_no == 0 && MCTrkID > -1 && NuPDGTruth[MCVtxID] == 14 )
            {
                SelectionTrackRange.at ( 0 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                SelectionTheta.at ( 0 )->Fill ( MCTheta[MCTrkID] );
                SelectionCosTheta.at ( 0 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                SelectionPhi.at ( 0 )->Fill ( MCPhi[MCTrkID] );
                SelectionEnergy.at ( 0 )->Fill ( MCEnergy[MCTrkID] );
                SelectionMomentum.at ( 0 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                SelXTrackStartEnd.at ( 0 )->Fill ( XMCTrackStart[MCTrkID] );
                SelXTrackStartEnd.at ( 0 )->Fill ( XMCTrackEnd[MCTrkID] );
                SelYTrackStartEnd.at ( 0 )->Fill ( YMCTrackStart[MCTrkID] );
                SelYTrackStartEnd.at ( 0 )->Fill ( YMCTrackEnd[MCTrkID] );
                SelZTrackStartEnd.at ( 0 )->Fill ( ZMCTrackStart[MCTrkID] );
                SelZTrackStartEnd.at ( 0 )->Fill ( ZMCTrackEnd[MCTrkID] );
                SelXVtxPosition.at ( 0 )->Fill ( XnuVtxTruth[MCVtxID] );
                SelYVtxPosition.at ( 0 )->Fill ( YnuVtxTruth[MCVtxID] );
                SelZVtxPosition.at ( 0 )->Fill ( ZnuVtxTruth[MCVtxID] );

                // True QE RES DIS
                if ( TruthMode[MCVtxID] == 0 )
                {
                    nuQE++;
                    SelectionTrackRange.at ( 4 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                    SelectionTheta.at ( 4 )->Fill ( MCTheta[MCTrkID] );
                    SelectionCosTheta.at ( 4 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                    SelectionPhi.at ( 4 )->Fill ( MCPhi[MCTrkID] );
                    SelectionEnergy.at ( 4 )->Fill ( MCEnergy[MCTrkID] );
                    SelectionMomentum.at ( 4 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                    SelXTrackStartEnd.at ( 4 )->Fill ( XMCTrackStart[MCTrkID] );
                    SelXTrackStartEnd.at ( 4 )->Fill ( XMCTrackEnd[MCTrkID] );
                    SelYTrackStartEnd.at ( 4 )->Fill ( YMCTrackStart[MCTrkID] );
                    SelYTrackStartEnd.at ( 4 )->Fill ( YMCTrackEnd[MCTrkID] );
                    SelZTrackStartEnd.at ( 4 )->Fill ( ZMCTrackStart[MCTrkID] );
                    SelZTrackStartEnd.at ( 4 )->Fill ( ZMCTrackEnd[MCTrkID] );
                    SelXVtxPosition.at ( 4 )->Fill ( XnuVtxTruth[MCVtxID] );
                    SelYVtxPosition.at ( 4 )->Fill ( YnuVtxTruth[MCVtxID] );
                    SelZVtxPosition.at ( 4 )->Fill ( ZnuVtxTruth[MCVtxID] );

                }
                else if ( TruthMode[MCVtxID] == 1 )
                {
                    nuRES++;
                    SelectionTrackRange.at ( 5 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                    SelectionTheta.at ( 5 )->Fill ( MCTheta[MCTrkID] );
                    SelectionCosTheta.at ( 5 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                    SelectionPhi.at ( 5 )->Fill ( MCPhi[MCTrkID] );
                    SelectionEnergy.at ( 5 )->Fill ( MCEnergy[MCTrkID] );
                    SelectionMomentum.at ( 5 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                    SelXTrackStartEnd.at ( 5 )->Fill ( XMCTrackStart[MCTrkID] );
                    SelXTrackStartEnd.at ( 5 )->Fill ( XMCTrackEnd[MCTrkID] );
                    SelYTrackStartEnd.at ( 5 )->Fill ( YMCTrackStart[MCTrkID] );
                    SelYTrackStartEnd.at ( 5 )->Fill ( YMCTrackEnd[MCTrkID] );
                    SelZTrackStartEnd.at ( 5 )->Fill ( ZMCTrackStart[MCTrkID] );
                    SelZTrackStartEnd.at ( 5 )->Fill ( ZMCTrackEnd[MCTrkID] );
                    SelXVtxPosition.at ( 5 )->Fill ( XnuVtxTruth[MCVtxID] );
                    SelYVtxPosition.at ( 5 )->Fill ( YnuVtxTruth[MCVtxID] );
                    SelZVtxPosition.at ( 5 )->Fill ( ZnuVtxTruth[MCVtxID] );

                }
                else if ( TruthMode[MCVtxID] == 2 )
                {
                    nuDIS++;
                    SelectionTrackRange.at ( 6 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                    SelectionTheta.at ( 6 )->Fill ( MCTheta[MCTrkID] );
                    SelectionCosTheta.at ( 6 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                    SelectionPhi.at ( 6 )->Fill ( MCPhi[MCTrkID] );
                    SelectionEnergy.at ( 6 )->Fill ( MCEnergy[MCTrkID] );
                    SelectionMomentum.at ( 6 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                    SelXTrackStartEnd.at ( 6 )->Fill ( XMCTrackStart[MCTrkID] );
                    SelXTrackStartEnd.at ( 6 )->Fill ( XMCTrackEnd[MCTrkID] );
                    SelYTrackStartEnd.at ( 6 )->Fill ( YMCTrackStart[MCTrkID] );
                    SelYTrackStartEnd.at ( 6 )->Fill ( YMCTrackEnd[MCTrkID] );
                    SelZTrackStartEnd.at ( 6 )->Fill ( ZMCTrackStart[MCTrkID] );
                    SelZTrackStartEnd.at ( 6 )->Fill ( ZMCTrackEnd[MCTrkID] );
                    SelXVtxPosition.at ( 6 )->Fill ( XnuVtxTruth[MCVtxID] );
                    SelYVtxPosition.at ( 6 )->Fill ( YnuVtxTruth[MCVtxID] );
                    SelZVtxPosition.at ( 6 )->Fill ( ZnuVtxTruth[MCVtxID] );
                }


                if ( inFV ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID] ) && inFV ( XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) )
                {
                    ContainedTracks++;

                    SelectionTrackRange.at ( 1 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                    SelectionTheta.at ( 1 )->Fill ( MCTheta[MCTrkID] );
                    SelectionCosTheta.at ( 1 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                    SelectionPhi.at ( 1 )->Fill ( MCPhi[MCTrkID] );
                    SelectionEnergy.at ( 1 )->Fill ( MCEnergy[MCTrkID] );
                    SelectionMomentum.at ( 1 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                    SelXTrackStartEnd.at ( 1 )->Fill ( XMCTrackStart[MCTrkID] );
                    SelXTrackStartEnd.at ( 1 )->Fill ( XMCTrackEnd[MCTrkID] );
                    SelYTrackStartEnd.at ( 1 )->Fill ( YMCTrackStart[MCTrkID] );
                    SelYTrackStartEnd.at ( 1 )->Fill ( YMCTrackEnd[MCTrkID] );
                    SelZTrackStartEnd.at ( 1 )->Fill ( ZMCTrackStart[MCTrkID] );
                    SelZTrackStartEnd.at ( 1 )->Fill ( ZMCTrackEnd[MCTrkID] );
                    SelXVtxPosition.at ( 1 )->Fill ( XnuVtxTruth[MCVtxID] );
                    SelYVtxPosition.at ( 1 )->Fill ( YnuVtxTruth[MCVtxID] );
                    SelZVtxPosition.at ( 1 )->Fill ( ZnuVtxTruth[MCVtxID] );

                    // True Contained QE RES DIS
                    if ( TruthMode[MCVtxID] == 0 )
                    {
                        nuQE++;
                        SelectionTrackRange.at ( 7 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                        SelectionTheta.at ( 7 )->Fill ( MCTheta[MCTrkID] );
                        SelectionCosTheta.at ( 7 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                        SelectionPhi.at ( 7 )->Fill ( MCPhi[MCTrkID] );
                        SelectionEnergy.at ( 7 )->Fill ( MCEnergy[MCTrkID] );
                        SelectionMomentum.at ( 7 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                        SelXTrackStartEnd.at ( 7 )->Fill ( XMCTrackStart[MCTrkID] );
                        SelXTrackStartEnd.at ( 7 )->Fill ( XMCTrackEnd[MCTrkID] );
                        SelYTrackStartEnd.at ( 7 )->Fill ( YMCTrackStart[MCTrkID] );
                        SelYTrackStartEnd.at ( 7 )->Fill ( YMCTrackEnd[MCTrkID] );
                        SelZTrackStartEnd.at ( 7 )->Fill ( ZMCTrackStart[MCTrkID] );
                        SelZTrackStartEnd.at ( 7 )->Fill ( ZMCTrackEnd[MCTrkID] );
                        SelXVtxPosition.at ( 7 )->Fill ( XnuVtxTruth[MCVtxID] );
                        SelYVtxPosition.at ( 7 )->Fill ( YnuVtxTruth[MCVtxID] );
                        SelZVtxPosition.at ( 7 )->Fill ( ZnuVtxTruth[MCVtxID] );

                    }
                    else if ( TruthMode[MCVtxID] == 1 )
                    {
                        nuRES++;
                        SelectionTrackRange.at ( 8 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                        SelectionTheta.at ( 8 )->Fill ( MCTheta[MCTrkID] );
                        SelectionCosTheta.at ( 8 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                        SelectionPhi.at ( 8 )->Fill ( MCPhi[MCTrkID] );
                        SelectionEnergy.at ( 8 )->Fill ( MCEnergy[MCTrkID] );
                        SelectionMomentum.at ( 8 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                        SelXTrackStartEnd.at ( 8 )->Fill ( XMCTrackStart[MCTrkID] );
                        SelXTrackStartEnd.at ( 8 )->Fill ( XMCTrackEnd[MCTrkID] );
                        SelYTrackStartEnd.at ( 8 )->Fill ( YMCTrackStart[MCTrkID] );
                        SelYTrackStartEnd.at ( 8 )->Fill ( YMCTrackEnd[MCTrkID] );
                        SelZTrackStartEnd.at ( 8 )->Fill ( ZMCTrackStart[MCTrkID] );
                        SelZTrackStartEnd.at ( 8 )->Fill ( ZMCTrackEnd[MCTrkID] );
                        SelXVtxPosition.at ( 8 )->Fill ( XnuVtxTruth[MCVtxID] );
                        SelYVtxPosition.at ( 8 )->Fill ( YnuVtxTruth[MCVtxID] );
                        SelZVtxPosition.at ( 8 )->Fill ( ZnuVtxTruth[MCVtxID] );

                    }
                    else if ( TruthMode[MCVtxID] == 2 )
                    {
                        nuDIS++;
                        SelectionTrackRange.at ( 9 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                        SelectionTheta.at ( 9 )->Fill ( MCTheta[MCTrkID] );
                        SelectionCosTheta.at ( 9 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                        SelectionPhi.at ( 9 )->Fill ( MCPhi[MCTrkID] );
                        SelectionEnergy.at ( 9 )->Fill ( MCEnergy[MCTrkID] );
                        SelectionMomentum.at ( 9 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                        SelXTrackStartEnd.at ( 9 )->Fill ( XMCTrackStart[MCTrkID] );
                        SelXTrackStartEnd.at ( 9 )->Fill ( XMCTrackEnd[MCTrkID] );
                        SelYTrackStartEnd.at ( 9 )->Fill ( YMCTrackStart[MCTrkID] );
                        SelYTrackStartEnd.at ( 9 )->Fill ( YMCTrackEnd[MCTrkID] );
                        SelZTrackStartEnd.at ( 9 )->Fill ( ZMCTrackStart[MCTrkID] );
                        SelZTrackStartEnd.at ( 9 )->Fill ( ZMCTrackEnd[MCTrkID] );
                        SelXVtxPosition.at ( 9 )->Fill ( XnuVtxTruth[MCVtxID] );
                        SelYVtxPosition.at ( 9 )->Fill ( YnuVtxTruth[MCVtxID] );
                        SelZVtxPosition.at ( 9 )->Fill ( ZnuVtxTruth[MCVtxID] );
                    }
                }
            }
            else if ( file_no == 1 && TrkID > -1 && TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1 && NuPDGTruth[MCVtxID] == 14 )
            {
                SelectionTrackRange.at ( 2 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                SelectionTheta.at ( 2 )->Fill ( MCTheta[MCTrkID] );
                SelectionCosTheta.at ( 2 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                SelectionPhi.at ( 2 )->Fill ( MCPhi[MCTrkID] );
                SelectionEnergy.at ( 2 )->Fill ( MCEnergy[MCTrkID] );
                SelectionMomentum.at ( 2 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                SelXTrackStartEnd.at ( 2 )->Fill ( XMCTrackStart[MCTrkID] );
                SelXTrackStartEnd.at ( 2 )->Fill ( XMCTrackEnd[MCTrkID] );
                SelYTrackStartEnd.at ( 2 )->Fill ( YMCTrackStart[MCTrkID] );
                SelYTrackStartEnd.at ( 2 )->Fill ( YMCTrackEnd[MCTrkID] );
                SelZTrackStartEnd.at ( 2 )->Fill ( ZMCTrackStart[MCTrkID] );
                SelZTrackStartEnd.at ( 2 )->Fill ( ZMCTrackEnd[MCTrkID] );
                SelXVtxPosition.at ( 2 )->Fill ( XnuVtxTruth[MCVtxID] );
                SelYVtxPosition.at ( 2 )->Fill ( YnuVtxTruth[MCVtxID] );
                SelZVtxPosition.at ( 2 )->Fill ( ZnuVtxTruth[MCVtxID] );
                
                if ( inFV ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID] ) && inFV ( XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) )
                {
                    SelectionTrackRange.at ( 3 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                    SelectionTheta.at ( 3 )->Fill ( MCTheta[MCTrkID] );
                    SelectionCosTheta.at ( 3 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                    SelectionPhi.at ( 3 )->Fill ( MCPhi[MCTrkID] );
                    SelectionEnergy.at ( 3 )->Fill ( MCEnergy[MCTrkID] );
                    SelectionMomentum.at ( 3 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                    SelXTrackStartEnd.at ( 3 )->Fill ( XMCTrackStart[MCTrkID] );
                    SelXTrackStartEnd.at ( 3 )->Fill ( XMCTrackEnd[MCTrkID] );
                    SelYTrackStartEnd.at ( 3 )->Fill ( YMCTrackStart[MCTrkID] );
                    SelYTrackStartEnd.at ( 3 )->Fill ( YMCTrackEnd[MCTrkID] );
                    SelZTrackStartEnd.at ( 3 )->Fill ( ZMCTrackStart[MCTrkID] );
                    SelZTrackStartEnd.at ( 3 )->Fill ( ZMCTrackEnd[MCTrkID] );
                    SelXVtxPosition.at ( 3 )->Fill ( XnuVtxTruth[MCVtxID] );
                    SelYVtxPosition.at ( 3 )->Fill ( YnuVtxTruth[MCVtxID] );
                    SelZVtxPosition.at ( 3 )->Fill ( ZnuVtxTruth[MCVtxID] );
                }


                // Selection QE RES DIS
                if ( TruthMode[MCVtxID] == 0 )
                {
                    nuQE++;
                    SelectionTrackRange.at ( 10 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                    SelectionTheta.at ( 10 )->Fill ( MCTheta[MCTrkID] );
                    SelectionCosTheta.at ( 10 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                    SelectionPhi.at ( 10 )->Fill ( MCPhi[MCTrkID] );
                    SelectionEnergy.at ( 10 )->Fill ( MCEnergy[MCTrkID] );
                    SelectionMomentum.at ( 10 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                    SelXTrackStartEnd.at ( 10 )->Fill ( XMCTrackStart[MCTrkID] );
                    SelXTrackStartEnd.at ( 10 )->Fill ( XMCTrackEnd[MCTrkID] );
                    SelYTrackStartEnd.at ( 10 )->Fill ( YMCTrackStart[MCTrkID] );
                    SelYTrackStartEnd.at ( 10 )->Fill ( YMCTrackEnd[MCTrkID] );
                    SelZTrackStartEnd.at ( 10 )->Fill ( ZMCTrackStart[MCTrkID] );
                    SelZTrackStartEnd.at ( 10 )->Fill ( ZMCTrackEnd[MCTrkID] );
                    SelXVtxPosition.at ( 10 )->Fill ( XnuVtxTruth[MCVtxID] );
                    SelYVtxPosition.at ( 10 )->Fill ( YnuVtxTruth[MCVtxID] );
                    SelZVtxPosition.at ( 10 )->Fill ( ZnuVtxTruth[MCVtxID] );

                }
                else if ( TruthMode[MCVtxID] == 1 )
                {
                    nuRES++;
                    SelectionTrackRange.at ( 11 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                    SelectionTheta.at ( 11 )->Fill ( MCTheta[MCTrkID] );
                    SelectionCosTheta.at ( 11 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                    SelectionPhi.at ( 11 )->Fill ( MCPhi[MCTrkID] );
                    SelectionEnergy.at ( 11 )->Fill ( MCEnergy[MCTrkID] );
                    SelectionMomentum.at ( 11 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                    SelXTrackStartEnd.at ( 11 )->Fill ( XMCTrackStart[MCTrkID] );
                    SelXTrackStartEnd.at ( 11 )->Fill ( XMCTrackEnd[MCTrkID] );
                    SelYTrackStartEnd.at ( 11 )->Fill ( YMCTrackStart[MCTrkID] );
                    SelYTrackStartEnd.at ( 11 )->Fill ( YMCTrackEnd[MCTrkID] );
                    SelZTrackStartEnd.at ( 11 )->Fill ( ZMCTrackStart[MCTrkID] );
                    SelZTrackStartEnd.at ( 11 )->Fill ( ZMCTrackEnd[MCTrkID] );
                    SelXVtxPosition.at ( 11 )->Fill ( XnuVtxTruth[MCVtxID] );
                    SelYVtxPosition.at ( 11 )->Fill ( YnuVtxTruth[MCVtxID] );
                    SelZVtxPosition.at ( 11 )->Fill ( ZnuVtxTruth[MCVtxID] );

                }
                else if ( TruthMode[MCVtxID] == 2 )
                {
                    nuDIS++;
                    SelectionTrackRange.at ( 12 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                    SelectionTheta.at ( 12 )->Fill ( MCTheta[MCTrkID] );
                    SelectionCosTheta.at ( 12 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                    SelectionPhi.at ( 12 )->Fill ( MCPhi[MCTrkID] );
                    SelectionEnergy.at ( 12 )->Fill ( MCEnergy[MCTrkID] );
                    SelectionMomentum.at ( 12 )->Fill ( TrueLeptonMomentum[MCVtxID] );

                    SelXTrackStartEnd.at ( 12 )->Fill ( XMCTrackStart[MCTrkID] );
                    SelXTrackStartEnd.at ( 12 )->Fill ( XMCTrackEnd[MCTrkID] );
                    SelYTrackStartEnd.at ( 12 )->Fill ( YMCTrackStart[MCTrkID] );
                    SelYTrackStartEnd.at ( 12 )->Fill ( YMCTrackEnd[MCTrkID] );
                    SelZTrackStartEnd.at ( 12 )->Fill ( ZMCTrackStart[MCTrkID] );
                    SelZTrackStartEnd.at ( 12 )->Fill ( ZMCTrackEnd[MCTrkID] );
                    SelXVtxPosition.at ( 12 )->Fill ( XnuVtxTruth[MCVtxID] );
                    SelYVtxPosition.at ( 12 )->Fill ( YnuVtxTruth[MCVtxID] );
                    SelZVtxPosition.at ( 12 )->Fill ( ZnuVtxTruth[MCVtxID] );
                }

                // Fill systematic errors independet of CC or NC
                SelectionTrackRange.back()->Fill ( CalcLength ( XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID] ),1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelectionTheta.back()->Fill ( TrackTheta[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelectionCosTheta.back()->Fill ( cos ( TrackTheta[TrkID] ),1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelectionPhi.back()->Fill ( TrackPhi[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelectionEnergy.back()->Fill ( KineticEnergy[TrkID][2]/1000,1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelectionMomentum.back()->Fill ( TrackMomentum[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );

                SelXTrackStartEnd.back()->Fill ( XTrackStart[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelXTrackStartEnd.back()->Fill ( XTrackEnd[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelYTrackStartEnd.back()->Fill ( YTrackStart[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelYTrackStartEnd.back()->Fill ( YTrackEnd[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelZTrackStartEnd.back()->Fill ( ZTrackStart[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelZTrackStartEnd.back()->Fill ( ZTrackEnd[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelXVtxPosition.back()->Fill ( XVertexPosition[VtxID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelYVtxPosition.back()->Fill ( YVertexPosition[VtxID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
                SelZVtxPosition.back()->Fill ( ZVertexPosition[VtxID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[MCVtxID] ) );
            }
        }
        std::cout << Signal << " " << nubar << " " << nue << " " << NCnu << " " << Cosmic << " " << UnknownOrigin << std::endl;
        std::cout << nuQE << " " << nuRES << " " << nuDIS << " " << nuCOH << std::endl;

        std::cout << "Number of negative phi in " << GenLabel.at ( file_no ) << " : " << negPhi << std::endl;
        std::cout << "Number of positive phi in " << GenLabel.at ( file_no ) << " : " << posPhi << std::endl;

        std::cout << "MC truth with contained tracks : " << ContainedTracks << std::endl;

        ChainVec.at ( file_no )->ResetBranchAddresses();
    }
    
    for ( unsigned int eff_no = 0; eff_no < EfficiencyLabel.size(); eff_no++ )
    {
        std::cout << TEfficiency::CheckConsistency ( *SelectionMomentum.at ( 2+eff_no ), *SelectionMomentum.at ( eff_no ),"w" ) << " " << TEfficiency::CheckBinning ( *SelectionMomentum.at ( 2+eff_no ), *SelectionMomentum.at ( eff_no ) ) << std::endl;
        EffTrackRange.push_back ( new TEfficiency ( *SelectionTrackRange.at ( 2+eff_no ),*SelectionTrackRange.at ( eff_no ) ) );
        EffEnergy.push_back ( new TEfficiency ( *SelectionEnergy.at ( 2+eff_no ),*SelectionEnergy.at ( eff_no ) ) );
        EffMomentum.push_back ( new TEfficiency ( *SelectionMomentum.at ( 2+eff_no ),*SelectionMomentum.at ( eff_no ) ) );
        EffTheta.push_back ( new TEfficiency ( *SelectionTheta.at ( 2+eff_no ),*SelectionTheta.at ( eff_no ) ) );
        EffCosTheta.push_back ( new TEfficiency ( *SelectionCosTheta.at ( 2+eff_no ),*SelectionCosTheta.at ( eff_no ) ) );
        EffPhi.push_back ( new TEfficiency ( *SelectionPhi.at ( 2+eff_no ),*SelectionPhi.at ( eff_no ) ) );

        EffXTrackStartEnd.push_back ( new TEfficiency ( *SelXTrackStartEnd.at ( 2+eff_no ),*SelXTrackStartEnd.at ( eff_no ) ) );
        EffYTrackStartEnd.push_back ( new TEfficiency ( *SelYTrackStartEnd.at ( 2+eff_no ),*SelYTrackStartEnd.at ( eff_no ) ) );
        EffZTrackStartEnd.push_back ( new TEfficiency ( *SelZTrackStartEnd.at ( 2+eff_no ),*SelZTrackStartEnd.at ( eff_no ) ) );

        EffXVtxPosition.push_back ( new TEfficiency ( *SelXVtxPosition.at ( 2+eff_no ),*SelXVtxPosition.at ( eff_no ) ) );
        EffYVtxPosition.push_back ( new TEfficiency ( *SelYVtxPosition.at ( 2+eff_no ),*SelYVtxPosition.at ( eff_no ) ) );
        EffZVtxPosition.push_back ( new TEfficiency ( *SelZVtxPosition.at ( 2+eff_no ),*SelZVtxPosition.at ( eff_no ) ) );
    }

    for ( unsigned int eff_no = 0; eff_no < EffIntLabel.size(); eff_no++ )
    {
        EffTrackRange.push_back ( new TEfficiency ( *SelectionTrackRange.at ( 10+eff_no ),*SelectionTrackRange.at ( 4+eff_no ) ) );
        EffEnergy.push_back ( new TEfficiency ( *SelectionEnergy.at ( 10+eff_no ),*SelectionEnergy.at ( 4+eff_no ) ) );
        EffMomentum.push_back ( new TEfficiency ( *SelectionMomentum.at ( 10+eff_no ),*SelectionMomentum.at ( 4+eff_no ) ) );
        EffTheta.push_back ( new TEfficiency ( *SelectionTheta.at ( 10+eff_no ),*SelectionTheta.at ( 4+eff_no ) ) );
        EffCosTheta.push_back ( new TEfficiency ( *SelectionCosTheta.at ( 10+eff_no ),*SelectionCosTheta.at ( 4+eff_no ) ) );
        EffPhi.push_back ( new TEfficiency ( *SelectionPhi.at ( 10+eff_no ),*SelectionPhi.at ( 4+eff_no ) ) );

        EffXTrackStartEnd.push_back ( new TEfficiency ( *SelXTrackStartEnd.at ( 10+eff_no ),*SelXTrackStartEnd.at ( 4+eff_no ) ) );
        EffYTrackStartEnd.push_back ( new TEfficiency ( *SelYTrackStartEnd.at ( 10+eff_no ),*SelYTrackStartEnd.at ( 4+eff_no ) ) );
        EffZTrackStartEnd.push_back ( new TEfficiency ( *SelZTrackStartEnd.at ( 10+eff_no ),*SelZTrackStartEnd.at ( 4+eff_no ) ) );

        EffXVtxPosition.push_back ( new TEfficiency ( *SelXVtxPosition.at ( 10+eff_no ),*SelXVtxPosition.at ( 4+eff_no ) ) );
        EffYVtxPosition.push_back ( new TEfficiency ( *SelYVtxPosition.at ( 10+eff_no ),*SelYVtxPosition.at ( 4+eff_no ) ) );
        EffZVtxPosition.push_back ( new TEfficiency ( *SelZVtxPosition.at ( 10+eff_no ),*SelZVtxPosition.at ( 4+eff_no ) ) );
    }

    LegendEffInt->AddEntry ( EffTrackRange.at ( 2 ), ( EffIntLabel.at ( 0 ) ).c_str(),"lf" );
    LegendEffInt->AddEntry ( EffTrackRange.at ( 3 ), ( EffIntLabel.at ( 1 ) ).c_str(),"lf" );
    LegendEffInt->AddEntry ( EffTrackRange.at ( 4 ), ( EffIntLabel.at ( 2 ) ).c_str(),"lf" );

    LegendEfficiency->AddEntry ( EffTrackRange.at ( 1 ), ( EfficiencyLabel.at ( 1 ) ).c_str(),"lf" );
    LegendEfficiency->AddEntry ( EffTrackRange.at ( 0 ), ( EfficiencyLabel.at ( 0 ) ).c_str(),"lf" );

    TCanvas *Canvas1 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Track Range", "Efficiency OnBeam Minus OffBeam Track Range", 1400, 1000 );
    TMultiGraph *MGTrackRange = new TMultiGraph();
    std::vector<TGraphAsymmErrors*> TrackRangeGraphs;
    EffTrackRange.at ( 0 )->SetLineWidth ( 2 );
    EffTrackRange.at ( 0 )->SetLineColor ( 8 );
    EffTrackRange.at ( 0 )->SetFillColorAlpha ( 8,0.5 );
    EffTrackRange.at ( 1 )->SetLineWidth ( 2 );
    EffTrackRange.at ( 1 )->SetLineColor ( 9 );
    EffTrackRange.at ( 1 )->SetFillColorAlpha ( 9,0.5 );
    TrackRangeGraphs.push_back(EffCosTheta.at ( 0 )->CreateGraph());
    TrackRangeGraphs.push_back(EffCosTheta.at ( 1 )->CreateGraph());
    for(const auto& Graph : TrackRangeGraphs) 
    {
        TGraphAsymmErrors* TempGraph = (TGraphAsymmErrors*) Graph->Clone();
        MGTrackRange->Add ( TempGraph );
    }
    Canvas1->cd();
    MGTrackRange->Draw ( "2AP" );
    for(const auto& Graph : TrackRangeGraphs)
    {
        for(int n = 0; n<Graph->GetN(); n++) Graph->SetPointError(n,Graph->GetErrorXlow(n),Graph->GetErrorXhigh(n),0,0);
        Graph->Draw("sameP");
    }
    MGTrackRange->GetXaxis()->SetTitle ( "Track Range [cm]" );
    MGTrackRange->GetYaxis()->SetTitle ( "Efficiency" );
    LegendEfficiency->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas1->SaveAs ( ( "EffMCRange"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas1Int = new TCanvas ( "Interaction Efficiency OnBeam Minus OffBeam Track Range", "Interaction Efficiency OnBeam Minus OffBeam Track Range", 1400, 1000 );
    TMultiGraph *MGTrackRangeInt = new TMultiGraph();
    std::vector<TGraphAsymmErrors*> TrackRangeGraphsInt;
    EffTrackRange.at ( 2 )->SetLineWidth ( 2 );
    EffTrackRange.at ( 2 )->SetFillColorAlpha (  ColorMap.at ( 0 ),0.5 );
    EffTrackRange.at ( 2 )->SetLineColor (  ColorMap.at ( 0 ) );
    EffTrackRange.at ( 3 )->SetLineWidth ( 2 );
    EffTrackRange.at ( 3 )->SetFillColorAlpha (  ColorMap.at ( 1 ),0.5 );
    EffTrackRange.at ( 3 )->SetLineColor (  ColorMap.at ( 1 ) );
    EffTrackRange.at ( 4 )->SetLineWidth ( 2 );
    EffTrackRange.at ( 4 )->SetFillColorAlpha (  ColorMap.at ( 2 ),0.5 );
    EffTrackRange.at ( 4 )->SetLineColor (  ColorMap.at ( 2 ) );
    TrackRangeGraphsInt.push_back(EffTrackRange.at ( 2 )->CreateGraph());
    TrackRangeGraphsInt.push_back(EffTrackRange.at ( 3 )->CreateGraph());
    TrackRangeGraphsInt.push_back(EffTrackRange.at ( 4 )->CreateGraph());
    for(const auto& Graph : TrackRangeGraphsInt) 
    {
        TGraphAsymmErrors* TempGraph = (TGraphAsymmErrors*) Graph->Clone();
        MGTrackRangeInt->Add ( TempGraph );
    }
    Canvas1Int->cd();
    MGTrackRangeInt->Draw ( "2AP" );
    for(const auto& Graph : TrackRangeGraphsInt)
    {
        for(int n = 0; n<Graph->GetN(); n++) Graph->SetPointError(n,Graph->GetErrorXlow(n),Graph->GetErrorXhigh(n),0,0);
        Graph->Draw("sameP");
    }
    MGTrackRangeInt->GetXaxis()->SetTitle ( "Track Range [cm]" );
    MGTrackRangeInt->GetYaxis()->SetTitle ( "Efficiency" );
    LegendEffInt->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas1Int->SaveAs ( ( "EffIntMCRange"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas2 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Theta-Angle", "Efficiency OnBeam Minus OffBeam Theta-Angle", 1400, 1000 );
    TMultiGraph *MGTheta = new TMultiGraph();
    EffTheta.at ( 0 )->SetLineWidth ( 2 );
    EffTheta.at ( 0 )->SetLineColor ( 8 );
    EffTheta.at ( 1 )->SetLineWidth ( 2 );
    EffTheta.at ( 1 )->SetLineColor ( 9 );
    MGTheta->Add ( EffTheta.at ( 0 )->CreateGraph() );
    MGTheta->Add ( EffTheta.at ( 1 )->CreateGraph() );
    Canvas2->cd();
    MGTheta->Draw ( "AP" );
    MGTheta->GetXaxis()->SetTitle ( "#theta angle [rad]" );
    MGTheta->GetYaxis()->SetTitle ( "Efficiency" );
    LegendEfficiency->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas2->SaveAs ( ( "EffMCTheta"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas2b = new TCanvas ( "Efficiency OnBeam Minus OffBeam Cos Theta-Angle", "Efficiency OnBeam Minus OffBeam Cos Theta-Angle", 1400, 1000 );
    TMultiGraph *MGCosTheta = new TMultiGraph();
    std::vector<TGraphAsymmErrors*> CosThetaGraphs;
    EffCosTheta.at ( 0 )->SetLineWidth ( 2 );
    EffCosTheta.at ( 0 )->SetLineColor ( 8 );
    EffCosTheta.at ( 0 )->SetFillColorAlpha ( 8,0.5 );
    EffCosTheta.at ( 1 )->SetLineWidth ( 2 );
    EffCosTheta.at ( 1 )->SetLineColor ( 9 );
    EffCosTheta.at ( 1 )->SetFillColorAlpha ( 9,0.5 );
    CosThetaGraphs.push_back(EffCosTheta.at ( 0 )->CreateGraph());
    CosThetaGraphs.push_back(EffCosTheta.at ( 1 )->CreateGraph());
    for(const auto& Graph : CosThetaGraphs) 
    {
        TGraphAsymmErrors* TempGraph = (TGraphAsymmErrors*) Graph->Clone();
        MGCosTheta->Add ( TempGraph );
    }
    Canvas2b->cd();
    MGCosTheta->Draw ( "2AP" );
    for(const auto& Graph : CosThetaGraphs)
    {
        for(int n = 0; n<Graph->GetN(); n++) Graph->SetPointError(n,Graph->GetErrorXlow(n),Graph->GetErrorXhigh(n),0,0);
        Graph->Draw("sameP");
    }
    MGCosTheta->GetXaxis()->SetTitle ( "cos(#theta)" );
    MGCosTheta->GetYaxis()->SetTitle ( "Efficiency" );
    LegendEfficiency->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas2b->SaveAs ( ( "EffMCCosTheta"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas2bInt = new TCanvas ( "Interaction Efficiency OnBeam Minus OffBeam Cos Theta-Angle", "Interaction Efficiency OnBeam Minus OffBeam Cos Theta-Angle", 1400, 1000 );
    TMultiGraph *MGCosThetaInt = new TMultiGraph();
    std::vector<TGraphAsymmErrors*> CosThetaGraphsInt;
    EffCosTheta.at ( 2 )->SetLineWidth ( 2 );
    EffCosTheta.at ( 2 )->SetFillColorAlpha (  ColorMap.at ( 0 ),0.5 );
    EffCosTheta.at ( 2 )->SetLineColor (  ColorMap.at ( 0 ) );
    EffCosTheta.at ( 3 )->SetLineWidth ( 2 );
    EffCosTheta.at ( 3 )->SetFillColorAlpha (  ColorMap.at ( 1 ),0.5 );
    EffCosTheta.at ( 3 )->SetLineColor (  ColorMap.at ( 1 ) );
    EffCosTheta.at ( 4 )->SetLineWidth ( 2 );
    EffCosTheta.at ( 4 )->SetFillColorAlpha (  ColorMap.at ( 2 ),0.5 );
    EffCosTheta.at ( 4 )->SetLineColor (  ColorMap.at ( 2 ) );
    CosThetaGraphsInt.push_back(EffCosTheta.at ( 2 )->CreateGraph());
    CosThetaGraphsInt.push_back(EffCosTheta.at ( 3 )->CreateGraph());
    CosThetaGraphsInt.push_back(EffCosTheta.at ( 4 )->CreateGraph());
    for(const auto& Graph : CosThetaGraphsInt) 
    {
        TGraphAsymmErrors* TempGraph = (TGraphAsymmErrors*) Graph->Clone();
        MGCosThetaInt->Add ( TempGraph );
    }
    Canvas2bInt->cd();
    MGCosThetaInt->Draw ( "2AP" );
    for(const auto& Graph : CosThetaGraphsInt)
    {
        for(int n = 0; n<Graph->GetN(); n++) Graph->SetPointError(n,Graph->GetErrorXlow(n),Graph->GetErrorXhigh(n),0,0);
        Graph->Draw("sameP");
    }
    MGCosThetaInt->GetXaxis()->SetTitle ( "cos(#theta)" );
    MGCosThetaInt->GetYaxis()->SetTitle ( "Efficiency" );
    LegendEffInt->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas2bInt->SaveAs ( ( "EffIntMCCosTheta"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas4 = new TCanvas ( "Efficiency Energy", "Efficiency Energy", 1400, 1000 );
    TMultiGraph *MGEnergy = new TMultiGraph();
    EffEnergy.at ( 0 )->SetLineWidth ( 2 );
    EffEnergy.at ( 0 )->SetLineColor ( 8 );
    EffEnergy.at ( 1 )->SetLineWidth ( 2 );
    EffEnergy.at ( 1 )->SetLineColor ( 9 );
    MGEnergy->Add ( EffEnergy.at ( 0 )->CreateGraph() );
    MGEnergy->Add ( EffEnergy.at ( 1 )->CreateGraph() );
    Canvas4->cd();
    MGEnergy->Draw ( "AP" );
    MGEnergy->GetXaxis()->SetTitle ( "Muon Energy [GeV]" );
    MGEnergy->GetYaxis()->SetTitle ( "Efficiency" );
    TextSimulation.Draw();
    TextSelection.Draw();
    LegendEfficiency->Draw();
//     Canvas4->SaveAs ( ( "EffMCEnergy"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas4a = new TCanvas ( "Efficiency Momentum", "Efficiency Momentum", 1400, 1000 );
    TMultiGraph *MGMomentum = new TMultiGraph();
    std::vector<TGraphAsymmErrors*> MomentumGraphs;
    EffMomentum.at ( 0 )->SetLineWidth ( 2 );
    EffMomentum.at ( 0 )->SetLineColor ( 8 );
    EffMomentum.at ( 0 )->SetFillColorAlpha ( 8,0.5 );
    EffMomentum.at ( 1 )->SetLineWidth ( 2 );
    EffMomentum.at ( 1 )->SetLineColor ( 9 );
    EffMomentum.at ( 1 )->SetFillColorAlpha ( 9,0.5 );
    MomentumGraphs.push_back(EffMomentum.at ( 0 )->CreateGraph());
    MomentumGraphs.push_back(EffMomentum.at ( 1 )->CreateGraph());
    for(const auto& Graph : MomentumGraphs) 
    {
        TGraphAsymmErrors* TempGraph = (TGraphAsymmErrors*) Graph->Clone();
        MGMomentum->Add ( TempGraph );
    }
    Canvas4a->cd();
    MGMomentum->Draw ( "2AP" );
    for(const auto& Graph : MomentumGraphs)
    {
        for(int n = 0; n<Graph->GetN(); n++) Graph->SetPointError(n,Graph->GetErrorXlow(n),Graph->GetErrorXhigh(n),0,0);
        Graph->Draw("sameP");
    }
    Canvas4a->cd();
    MGMomentum->GetXaxis()->SetTitle ( "Muon Momentum [GeV/c]" );
    MGMomentum->GetYaxis()->SetTitle ( "Efficiency" );
    MGMomentum->SetMaximum(1.1);
    LegendEfficiency->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas4a->SaveAs ( ( "EffMCMomentum"+SelectionLabel+"."+FileType ).c_str() );

    LegendEffInt->SetX1NDC ( 0.65 );
    LegendEffInt->SetY1NDC ( 0.65 );
    LegendEffInt->SetX2NDC ( 0.85 );
    LegendEffInt->SetY2NDC ( 0.85 );

    TCanvas *Canvas4aInt = new TCanvas ( "Interaction Efficiency Momentum", "Interaction Efficiency Momentum", 1400, 1000 );
    TMultiGraph *MGMomentumInt = new TMultiGraph();
    std::vector<TGraphAsymmErrors*> MomentumGraphsInt;
    EffMomentum.at ( 2 )->SetLineWidth ( 2 );
    EffMomentum.at ( 2 )->SetFillColorAlpha (  ColorMap.at ( 0 ),0.5 );
    EffMomentum.at ( 2 )->SetLineColor (  ColorMap.at ( 0 ) );
    EffMomentum.at ( 3 )->SetLineWidth ( 2 );
    EffMomentum.at ( 3 )->SetFillColorAlpha (  ColorMap.at ( 1 ),0.5 );
    EffMomentum.at ( 3 )->SetLineColor (  ColorMap.at ( 1 ) );
    EffMomentum.at ( 4 )->SetLineWidth ( 2 );
    EffMomentum.at ( 4 )->SetFillColorAlpha (  ColorMap.at ( 2 ),0.5 );
    EffMomentum.at ( 4 )->SetLineColor (  ColorMap.at ( 2 ) );
    MomentumGraphsInt.push_back(EffMomentum.at ( 2 )->CreateGraph());
    MomentumGraphsInt.push_back(EffMomentum.at ( 3 )->CreateGraph());
    MomentumGraphsInt.push_back(EffMomentum.at ( 4 )->CreateGraph());
    for(const auto& Graph : MomentumGraphsInt) 
    {
        TGraphAsymmErrors* TempGraph = (TGraphAsymmErrors*) Graph->Clone();
        MGMomentumInt->Add ( TempGraph );
    }
    Canvas4aInt->cd();
    MGMomentumInt->Draw ( "2AP" );
    for(const auto& Graph : MomentumGraphsInt)
    {
        for(int n = 0; n<Graph->GetN(); n++) Graph->SetPointError(n,Graph->GetErrorXlow(n),Graph->GetErrorXhigh(n),0,0);
        Graph->Draw("sameP");
    }
    MGMomentumInt->GetXaxis()->SetTitle ( "Muon Momentum [GeV/c]" );
    MGMomentumInt->GetYaxis()->SetTitle ( "Efficiency" );
    MGMomentumInt->SetMaximum(0.25);
    LegendEffInt->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas4aInt->SaveAs ( ( "EffIntMCMomentum"+SelectionLabel+"."+FileType ).c_str() );

    LegendEfficiency->SetX1NDC ( 0.57 );
    LegendEfficiency->SetY1NDC ( 0.33 );
    LegendEfficiency->SetX2NDC ( 0.87 );
    LegendEfficiency->SetY2NDC ( 0.53 );

    TCanvas *Canvas3 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Phi-Angle", "Efficiency OnBeam Minus OffBeam Phi-Angle", 1400, 1000 );
    TMultiGraph *MGPhi = new TMultiGraph();
    std::vector<TGraphAsymmErrors*> PhiGraphs;
    EffPhi.at ( 0 )->SetLineWidth ( 2 );
    EffPhi.at ( 0 )->SetLineColor ( 8 );
    EffPhi.at ( 0 )->SetFillColorAlpha ( 8,0.5 );
    EffPhi.at ( 1 )->SetLineWidth ( 2 );
    EffPhi.at ( 1 )->SetLineColor ( 9 );
    EffPhi.at ( 1 )->SetFillColorAlpha ( 9,0.5 );
    PhiGraphs.push_back(EffPhi.at ( 0 )->CreateGraph());
    PhiGraphs.push_back(EffPhi.at ( 1 )->CreateGraph());
    for(const auto& Graph : PhiGraphs) 
    {
        TGraphAsymmErrors* TempGraph = (TGraphAsymmErrors*) Graph->Clone();
        MGPhi->Add ( TempGraph );
    }
    Canvas3->cd();
    MGPhi->Draw ( "2AP" );
    for(const auto& Graph : PhiGraphs)
    {
        for(int n = 0; n<Graph->GetN(); n++) Graph->SetPointError(n,Graph->GetErrorXlow(n),Graph->GetErrorXhigh(n),0,0);
        Graph->Draw("sameP");
    }
    MGPhi->GetXaxis()->SetTitle ( "#phi angle [rad]" );
    MGPhi->GetYaxis()->SetTitle ( "Efficiency" );
    LegendEfficiency->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas3->SaveAs ( ( "EffMCPhi"+SelectionLabel+"."+FileType ).c_str() );

    LegendEffInt->SetX1NDC ( 0.12 );
    LegendEffInt->SetY1NDC ( 0.12 );
    LegendEffInt->SetX2NDC ( 0.32 );
    LegendEffInt->SetY2NDC ( 0.32 );

    TCanvas *Canvas3Int = new TCanvas ( "Interaction Efficiency OnBeam Minus OffBeam Phi-Angle", "Interaction Efficiency OnBeam Minus OffBeam Phi-Angle", 1400, 1000 );
    TMultiGraph *MGPhiInt = new TMultiGraph();
    std::vector<TGraphAsymmErrors*> PhiGraphsInt;
    EffPhi.at ( 2 )->SetLineWidth ( 2 );
    EffPhi.at ( 2 )->SetFillColorAlpha (  ColorMap.at ( 0 ),0.5 );
    EffPhi.at ( 2 )->SetLineColor (  ColorMap.at ( 0 ) );
    EffPhi.at ( 3 )->SetLineWidth ( 2 );
    EffPhi.at ( 3 )->SetFillColorAlpha (  ColorMap.at ( 1 ),0.5 );
    EffPhi.at ( 3 )->SetLineColor (  ColorMap.at ( 1 ) );
    EffPhi.at ( 4 )->SetLineWidth ( 2 );
    EffPhi.at ( 4 )->SetFillColorAlpha (  ColorMap.at ( 2 ),0.5 );
    EffPhi.at ( 4 )->SetLineColor (  ColorMap.at ( 2 ) );
    PhiGraphsInt.push_back(EffPhi.at ( 2 )->CreateGraph());
    PhiGraphsInt.push_back(EffPhi.at ( 3 )->CreateGraph());
    PhiGraphsInt.push_back(EffPhi.at ( 4 )->CreateGraph());
    for(const auto& Graph : PhiGraphsInt) 
    {
        TGraphAsymmErrors* TempGraph = (TGraphAsymmErrors*) Graph->Clone();
        MGPhiInt->Add ( TempGraph );
    }
    Canvas3Int->cd();
    MGPhiInt->Draw ( "2AP" );
    for(const auto& Graph : PhiGraphsInt)
    {
        for(int n = 0; n<Graph->GetN(); n++) Graph->SetPointError(n,Graph->GetErrorXlow(n),Graph->GetErrorXhigh(n),0,0);
        Graph->Draw("sameP");
    }
    MGPhiInt->GetXaxis()->SetTitle ( "#phi angle [rad]" );
    MGPhiInt->GetYaxis()->SetTitle ( "Efficiency" );
//     MGPhiInt->SetMaximum(0.22);
    MGPhiInt->SetMinimum(0.0);
    LegendEffInt->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas3Int->SaveAs ( ( "EffIntMCPhi"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas5 = new TCanvas ( "Efficiency OnBeam Minus OffBeam X Start & End Point ", "Efficiency OnBeam Minus OffBeam X Start & End Point ", 1400, 1000 );
    TMultiGraph *MGXStart = new TMultiGraph();
    Canvas5->cd();
//     EffXTrackStartEnd.at ( 1 )->SetMaximum ( 1.5*GetMaximum(EffXTrackStartEnd) );
//     EffXTrackStartEnd.at ( 1 )->SetMinimum ( 0.0 );
    EffXTrackStartEnd.at ( 1 )->SetLineColor ( 9 );
    EffXTrackStartEnd.at ( 1 )->Draw ( "A" );
    EffXTrackStartEnd.at ( 0 )->SetLineColor ( 8 );
    EffXTrackStartEnd.at ( 0 )->Draw ( "SAME" );
    LegendEfficiency->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas5->SaveAs ( ( "EffMCXTrack"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas6 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Y Start & End Point ", "Efficiency OnBeam Minus OffBeam Y Start & End Point ", 1400, 1000 );
    TMultiGraph *MGYStart = new TMultiGraph();
    Canvas6->cd();
//     EffYTrackStartEnd.at ( 1 )->SetMaximum ( 1.8*GetMaximum(EffYTrackStartEnd) );
//     EffYTrackStartEnd.at ( 1 )->SetMinimum ( 0.0 );
    EffYTrackStartEnd.at ( 1 )->SetLineColor ( 9 );
    EffYTrackStartEnd.at ( 1 )->Draw ( "A" );
    EffYTrackStartEnd.at ( 0 )->SetLineColor ( 8 );
    EffYTrackStartEnd.at ( 0 )->Draw ( "SAME" );
    LegendEfficiency->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas6->SaveAs ( ( "EffMCYTrack"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas7 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Z Start & End Point ", "Efficiency OnBeam Minus OffBeam Z Start & End Point ", 1400, 1000 );
    TMultiGraph *MGZStart = new TMultiGraph();
    Canvas7->cd();
//     EffZTrackStartEnd.at ( 1 )->SetMaximum ( 1.5*GetMaximum(EffZTrackStartEnd) );
//     EffZTrackStartEnd.at ( 1 )->SetMinimum ( 0.0 );
    EffZTrackStartEnd.at ( 1 )->SetLineColor ( 9 );
    EffZTrackStartEnd.at ( 1 )->Draw ( "A" );
    EffZTrackStartEnd.at ( 0 )->SetLineColor ( 8 );
    EffZTrackStartEnd.at ( 0 )->Draw ( "SAME" );
    TextSimulation.Draw();
//     LegendMC->Draw();
//     Canvas7->SaveAs ( ( "EffMCZTrack"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas8 = new TCanvas ( "Efficiency OnBeam Minus OffBeam X Vertex Postion", "Efficiency OnBeam Minus OffBeam X Vertex Postion", 1400, 1000 );
    TMultiGraph *MGXVertex = new TMultiGraph();
    Canvas8->cd();
//     EffXVtxPosition.at ( 1 )->SetMaximum ( 1.5*GetMaximum(EffXVtxPosition) );
//     EffXVtxPosition.at ( 1 )->SetMinimum ( 0.0 );
    EffXVtxPosition.at ( 1 )->SetLineColor ( 9 );
    EffXVtxPosition.at ( 1 )->Draw ( "A" );
    EffXVtxPosition.at ( 0 )->SetLineColor ( 8 );
    EffXVtxPosition.at ( 0 )->Draw ( "SAME" );
    TextSimulation.Draw();
    TextSelection.Draw();
//     LegendMC->Draw();
//     Canvas8->SaveAs ( ( "EffMCXVertex"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas9 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Y Vertex Postion", "Efficiency OnBeam Minus OffBeam Y Vertex Postion", 1400, 1000 );
    TMultiGraph *MGYVertex = new TMultiGraph();
    Canvas9->cd();
//     EffYVtxPosition.at ( 1 )->SetMaximum ( 1.8*GetMaximum(EffYVtxPosition) );
//     EffYVtxPosition.at ( 1 )->SetMinimum ( 0.0 );
    EffYVtxPosition.at ( 1 )->SetLineColor ( 9 );
    EffYVtxPosition.at ( 1 )->Draw ( "A" );
    EffYVtxPosition.at ( 0 )->SetLineColor ( 8 );
    EffYVtxPosition.at ( 0 )->Draw ( "SAME" );
    TextSimulation.Draw();
    TextSelection.Draw();
//     LegendMC->Draw();
//     Canvas9->SaveAs ( ( "EffMCYVertex"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas10 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Z Vertex Postion", "Efficiency OnBeam Minus OffBeam Z Vertex Postion", 1400, 1000 );
    TMultiGraph *MGZVertex = new TMultiGraph();
    Canvas10->cd();
//     EffZVtxPosition.at ( 1 )->SetMaximum ( 1.5*GetMaximum(EffZVtxPosition) );
//     EffZVtxPosition.at ( 1 )->SetMinimum ( 0.0 );
    EffZVtxPosition.at ( 1 )->SetLineColor ( 9 );
    EffZVtxPosition.at ( 1 )->Draw ( "A" );
    EffZVtxPosition.at ( 0 )->SetLineColor ( 8 );
    EffZVtxPosition.at ( 0 )->Draw ( "SAME" );
    TextSimulation.Draw();
    TextSelection.Draw();
//     LegendMC->Draw();
//     Canvas10->SaveAs ( ( "EffMCZVertex"+SelectionLabel+"."+FileType ).c_str() );


    for ( unsigned int file_no = 0; file_no < ScalingFactors.size(); file_no++ )
    {
        SelectionTrackRange.at ( file_no )->Sumw2();
        SelectionTheta.at ( file_no )->Sumw2();
        SelectionCosTheta.at ( file_no )->Sumw2();
        SelectionPhi.at ( file_no )->Sumw2();
        SelectionEnergy.at ( file_no )->Sumw2();
        SelectionMomentum.at ( file_no )->Sumw2();
        SelXTrackStartEnd.at ( file_no )->Sumw2();
        SelYTrackStartEnd.at ( file_no )->Sumw2();
        SelZTrackStartEnd.at ( file_no )->Sumw2();
        SelXVtxPosition.at ( file_no )->Sumw2();
        SelYVtxPosition.at ( file_no )->Sumw2();
        SelZVtxPosition.at ( file_no )->Sumw2();
    }
    
    LegendMCSel->AddEntry(SelectionTrackRange.at ( 2 ), "After Selection","f");

    LegendMC->AddEntry ( SelectionTrackRange.at ( 0 ), ( MCLabel.at ( 0 ) ).c_str(),"f" );
    LegendMC->AddEntry ( SelectionTrackRange.at ( 1 ), ( MCLabel.at ( 1 ) ).c_str(),"f" );
//     LegendMC->AddEntry ( SelectionTrackRange.at ( 2 ), ( MCLabel.at ( 2 ) ).c_str(),"f" );
//     for ( unsigned int bgrhist_no = 0; bgrhist_no < BgrLabel.size(); bgrhist_no++ )
//     {
//         LegendMC->AddEntry( BgrTrackRange.at(bgrhist_no), (BgrLabel.at(bgrhist_no)).c_str(),"f" );
//     }

    LegendInt->AddEntry ( SelectionTrackRange.at ( 10 ), ( IntLabel.at ( 0 ) ).c_str(),"f" );
    LegendInt->AddEntry ( SelectionTrackRange.at ( 11 ), ( IntLabel.at ( 1 ) ).c_str(),"f" );
    LegendInt->AddEntry ( SelectionTrackRange.at ( 12 ), ( IntLabel.at ( 2 ) ).c_str(),"f" );

    LegendIntTrue->AddEntry ( SelectionTrackRange.at ( 4 ), ( IntTrueLabel.at ( 0 ) ).c_str(),"f" );
    LegendIntTrue->AddEntry ( SelectionTrackRange.at ( 5 ), ( IntTrueLabel.at ( 1 ) ).c_str(),"f" );
    LegendIntTrue->AddEntry ( SelectionTrackRange.at ( 6 ), ( IntTrueLabel.at ( 2 ) ).c_str(),"f" );


    TCanvas *Canvas11 = new TCanvas ( "OnBeam Minus OffBeam Track Range", "OnBeam Minus OffBeam Track Range", 1400, 1000 );
    Canvas11->cd();
    SelectionTrackRange.at ( 1 )->SetMaximum ( 1.1*GetMaximum ( SelectionTrackRange ) );
    SelectionTrackRange.at ( 1 )->SetMinimum ( 0.0 );
    SelectionTrackRange.at ( 1 )->SetFillColor ( 9 );
    SelectionTrackRange.at ( 1 )->Draw ( "E2" );
    SelectionTrackRange.at ( 0 )->SetFillColor ( 8 );
    SelectionTrackRange.at ( 0 )->Draw ( "E2SAME" );
    TextSimulation.Draw();
    TextSelection.Draw();
    LegendMC->Draw();
//     Canvas11->SaveAs ( ( "MCRange"+SelectionLabel+"."+FileType ).c_str() );
    
    TCanvas *Canvas11Sel = new TCanvas ( "OnBeam Minus OffBeam Sel Track Range", "OnBeam Minus OffBeam Sel Track Range", 1400, 1000 );
    Canvas11Sel->cd();
    SelectionTrackRange.at ( 2 )->SetMaximum ( 1.2*SelectionTrackRange.at ( 2 )->GetBinContent(SelectionTrackRange.at ( 2 )->GetMaximumBin()) );
    SelectionTrackRange.at ( 2 )->SetMinimum ( 0.0 );
    SelectionTrackRange.at ( 2 )->SetFillColor ( 46 );
    SelectionTrackRange.at ( 2 )->Draw ( "E2" );
    LegendMCSel->Draw();
    TextSelection.Draw();
    TextSimulation.Draw();
//     Canvas11Sel->SaveAs ( ( "MCRangeSel"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas11Int = new TCanvas ( "OnBeam Minus OffBeam Track Range Int", "OnBeam Minus OffBeam Track Range Int", 1400, 1000 );
    Canvas11Int->cd();
    SelectionTrackRange.at ( 10 )->SetMaximum ( 1.1*SelectionTrackRange.at ( 10 )->GetBinContent ( SelectionTrackRange.at ( 10 )->GetMaximumBin() ) );
    SelectionTrackRange.at ( 10 )->SetMinimum ( 0.0 );
    SelectionTrackRange.at ( 10 )->SetFillColor ( ColorMap.at ( 0 ) );
    SelectionTrackRange.at ( 10 )->Draw ( "E2" );
    for ( unsigned int iter = 11; iter < SelectionTrackRange.size()-1; iter++ )
    {
        SelectionTrackRange.at ( iter )->SetFillColor ( ColorMap.at ( iter-10 ) );
        SelectionTrackRange.at ( iter )->Draw ( "E2SAME" );
    }
    TextSimulation.Draw();
    TextSelection.Draw();
    LegendInt->Draw();
//     Canvas11Int->SaveAs ( ( "MCRange_Int"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas11IntTrue = new TCanvas ( "OnBeam Minus OffBeam Track Range Int True", "OnBeam Minus OffBeam Track Range Int Ture", 1400, 1000 );
    Canvas11IntTrue->cd();
    SelectionTrackRange.at ( 4 )->SetMaximum ( 1.1*SelectionTrackRange.at ( 4 )->GetBinContent ( SelectionTrackRange.at ( 4 )->GetMaximumBin() ) );
    SelectionTrackRange.at ( 4 )->SetMinimum ( 0.0 );
    SelectionTrackRange.at ( 4 )->SetFillColor ( ColorMap.at ( 0 ) );
    SelectionTrackRange.at ( 4 )->Draw ( "E2" );
    for ( unsigned int iter = 5; iter < 7; iter++ )
    {
        SelectionTrackRange.at ( iter )->SetFillColor ( ColorMap.at ( iter-4 ) );
        SelectionTrackRange.at ( iter )->Draw ( "E2SAME" );
    }
    TextSimulation.Draw();
    TextSelection.Draw();
    LegendIntTrue->Draw();
//     Canvas11Int->SaveAs ( ( "MCRange_Int_True"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas12 = new TCanvas ( "OnBeam Minus OffBeam Theta-Angle", "OnBeam Minus OffBeam Theta-Angle", 1400, 1000 );
    Canvas12->cd();
    SelectionTheta.at ( 1 )->SetMaximum ( 1.2*GetMaximum ( SelectionTheta ) );
    SelectionTheta.at ( 1 )->SetMinimum ( 0.0 );
    SelectionTheta.at ( 1 )->SetFillColor ( 9 );
    SelectionTheta.at ( 1 )->Draw ( "E2" );
    SelectionTheta.at ( 0 )->SetFillColor ( 8 );
    SelectionTheta.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    TextSelection.Draw();
    TextSimulation.Draw();
//     Canvas12->SaveAs ( ( "MCTheta"+SelectionLabel+"."+FileType ).c_str() );

    for ( auto& ThetaHistogram : SelectionTheta )
    {
        ThetaHistogram->Divide ( SinTheta,1. );
    }

    TCanvas *Canvas12a = new TCanvas ( "OnBeam Minus OffBeam Theta-Angle Omega", "OnBeam Minus OffBeam Theta-Angle Omega", 1400, 1000 );
    Canvas12a->cd();
    SelectionTheta.at ( 1 )->SetMaximum ( 1.2*GetMaximum ( SelectionTheta ) );
    SelectionTheta.at ( 1 )->SetMinimum ( 0.0 );
    SelectionTheta.at ( 1 )->GetYaxis()->SetTitle ( "No. of events" );
    SelectionTheta.at ( 1 )->SetFillColor ( 9 );
    SelectionTheta.at ( 1 )->Draw ( "E2" );
    SelectionTheta.at ( 0 )->SetFillColor ( 8 );
    SelectionTheta.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas12a->SaveAs ( ( "MCThetaOmega"+SelectionLabel+"."+FileType ).c_str() );

//     LegendMC->SetX1NDC ( 0.2 );
//     LegendMC->SetY1NDC ( 0.6 );
//     LegendMC->SetX2NDC ( 0.5 );
//     LegendMC->SetY2NDC ( 0.8 );

    TCanvas *Canvas12b = new TCanvas ( "OnBeam Minus OffBeam Cos Theta-Angle", "OnBeam Minus OffBeam Cos Theta-Angle", 1400, 1000 );
    Canvas12b->cd();
    SelectionCosTheta.at ( 1 )->SetMaximum ( 1.2*GetMaximum ( SelectionCosTheta ) );
    SelectionCosTheta.at ( 1 )->SetMinimum ( 0.0 );
    SelectionCosTheta.at ( 1 )->SetFillColor ( 9 );
    SelectionCosTheta.at ( 1 )->Draw ( "E2" );
    SelectionCosTheta.at ( 0 )->SetFillColor ( 8 );
    SelectionCosTheta.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas12b->SaveAs ( ( "MCCosTheta"+SelectionLabel+"."+FileType ).c_str() );
    
    TCanvas *Canvas12bSel = new TCanvas ( "OnBeam Minus OffBeam Sel Cos Theta-Angle", "OnBeam Minus OffBeam Sel Cos Theta-Angle", 1400, 1000 );
    Canvas12bSel->cd();
    SelectionCosTheta.at ( 2 )->SetMaximum ( 1.2*SelectionCosTheta.at ( 2 )->GetBinContent(SelectionCosTheta.at ( 2 )->GetMaximumBin()) );
    SelectionCosTheta.at ( 2 )->SetMinimum ( 0.0 );
    SelectionCosTheta.at ( 2 )->SetFillColor ( 46 );
    SelectionCosTheta.at ( 2 )->Draw ( "E2" );
    LegendMCSel->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas12bSel->SaveAs ( ( "MCCosThetaSel"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas12bInt = new TCanvas ( "OnBeam Minus OffBeam Cos Theta-Angle Int", "OnBeam Minus OffBeam Cos Theta-Angle Int", 1400, 1000 );
    Canvas12bInt->cd();
    SelectionCosTheta.at ( 10 )->SetMaximum ( 1.1*SelectionCosTheta.at ( 10 )->GetBinContent ( SelectionCosTheta.at ( 10 )->GetMaximumBin() ) );
    SelectionCosTheta.at ( 10 )->SetMinimum ( 0.0 );
    SelectionCosTheta.at ( 10 )->SetFillColor ( ColorMap.at ( 0 ) );
    SelectionCosTheta.at ( 10 )->Draw ( "E2" );
    for ( unsigned int iter = 11; iter < SelectionCosTheta.size()-1; iter++ )
    {
        SelectionCosTheta.at ( iter )->SetFillColor ( ColorMap.at ( iter-10 ) );
        SelectionCosTheta.at ( iter )->Draw ( "E2SAME" );
    }
    LegendInt->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas12bInt->SaveAs ( ( "MCCosTheta_Int"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas12bIntTrue = new TCanvas ( "OnBeam Minus OffBeam Cos Theta-Angle Int Ture", "OnBeam Minus OffBeam Cos Theta-Angle Int Ture", 1400, 1000 );
    Canvas12bIntTrue->cd();
    SelectionCosTheta.at ( 4 )->SetMaximum ( 1.1*SelectionCosTheta.at ( 4 )->GetBinContent ( SelectionCosTheta.at ( 4 )->GetMaximumBin() ) );
    SelectionCosTheta.at ( 4 )->SetMinimum ( 0.0 );
    SelectionCosTheta.at ( 4 )->SetFillColor ( ColorMap.at ( 0 ) );
    SelectionCosTheta.at ( 4 )->Draw ( "E2" );
    for ( unsigned int iter = 5; iter < 7; iter++ )
    {
        SelectionCosTheta.at ( iter )->SetFillColor ( ColorMap.at ( iter-4 ) );
        SelectionCosTheta.at ( iter )->Draw ( "E2SAME" );
    }
    LegendIntTrue->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas12bIntTrue->SaveAs ( ( "MCCosTheta_Int_True"+SelectionLabel+"."+FileType ).c_str() );

    LegendMC->SetX1NDC ( 0.5 );
    LegendMC->SetY1NDC ( 0.4 );
    LegendMC->SetX2NDC ( 0.8 );
    LegendMC->SetY2NDC ( 0.6 );

    TCanvas *Canvas13 = new TCanvas ( "OnBeam Minus OffBeam Phi-Angle", "OnBeam Minus OffBeam Phi-Angle", 1400, 1000 );
    Canvas13->cd();
    SelectionPhi.at ( 1 )->SetMaximum ( 1.2*GetMaximum ( SelectionPhi ) );
    SelectionPhi.at ( 1 )->SetMinimum ( 0.0 );
    SelectionPhi.at ( 1 )->SetFillColor ( 9 );
    SelectionPhi.at ( 1 )->Draw ( "E2" );
    SelectionPhi.at ( 0 )->SetFillColor ( 8 );
    SelectionPhi.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas13->SaveAs ( ( "MCPhi"+SelectionLabel+"."+FileType ).c_str() );
    
    LegendMCSel->SetX1NDC ( 0.6 );
    LegendMCSel->SetY1NDC ( 0.25 );
    LegendMCSel->SetX2NDC ( 0.8 );
    LegendMCSel->SetY2NDC ( 0.3 );
    
    TCanvas *Canvas13Sel = new TCanvas ( "OnBeam Minus OffBeam Sel Phi-Angle", "OnBeam Minus OffBeam Sel Phi-Angle", 1400, 1000 );
    Canvas13Sel->cd();
    SelectionPhi.at ( 2 )->SetMaximum ( 1.2*SelectionPhi.at ( 2 )->GetBinContent(SelectionPhi.at ( 2 )->GetMaximumBin()) );
    SelectionPhi.at ( 2 )->SetMinimum ( 0.0 );
    SelectionPhi.at ( 2 )->SetFillColor ( 46 );
    SelectionPhi.at ( 2 )->Draw ( "E2" );
    LegendMCSel->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas13Sel->SaveAs ( ( "MCPhiSel"+SelectionLabel+"."+FileType ).c_str() );

    LegendInt->SetX1NDC ( 0.45 );
    LegendInt->SetY1NDC ( 0.65 );
    LegendInt->SetX2NDC ( 0.85 );
    LegendInt->SetY2NDC ( 0.85 );

    TCanvas *Canvas13Int = new TCanvas ( "OnBeam Minus OffBeam Phi-Angle Int", "OnBeam Minus OffBeam Phi-Angle Int", 1400, 1000 );
    Canvas13Int->cd();
    SelectionPhi.at ( 10 )->SetMaximum ( 1.9*SelectionPhi.at ( 10 )->GetBinContent ( SelectionPhi.at ( 10 )->GetMaximumBin() ) );
    SelectionPhi.at ( 10 )->SetMinimum ( 0.0 );
    SelectionPhi.at ( 10 )->SetFillColor ( ColorMap.at ( 0 ) );
    SelectionPhi.at ( 10 )->Draw ( "E2" );
    for ( unsigned int iter = 11; iter < SelectionPhi.size()-1; iter++ )
    {
        SelectionPhi.at ( iter )->SetFillColor ( ColorMap.at ( iter-10 ) );
        SelectionPhi.at ( iter )->Draw ( "E2SAME" );
    }
    LegendInt->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas13Int->SaveAs ( ( "MCPhi_Int"+SelectionLabel+"."+FileType ).c_str() );

    LegendIntTrue->SetX1NDC ( 0.45 );
    LegendIntTrue->SetY1NDC ( 0.65 );
    LegendIntTrue->SetX2NDC ( 0.85 );
    LegendIntTrue->SetY2NDC ( 0.85 );

    TCanvas *Canvas13IntTrue = new TCanvas ( "OnBeam Minus OffBeam Phi-Angle Int Ture", "OnBeam Minus OffBeam Phi-Angle Int Ture", 1400, 1000 );
    Canvas13IntTrue->cd();
    SelectionPhi.at ( 4 )->SetMaximum ( 1.9*SelectionPhi.at ( 4 )->GetBinContent ( SelectionPhi.at ( 4 )->GetMaximumBin() ) );
    SelectionPhi.at ( 4 )->SetMinimum ( 0.0 );
    SelectionPhi.at ( 4 )->SetFillColor ( ColorMap.at ( 0 ) );
    SelectionPhi.at ( 4 )->Draw ( "E2" );
    for ( unsigned int iter = 5; iter < 7; iter++ )
    {
        SelectionPhi.at ( iter )->SetFillColor ( ColorMap.at ( iter-4 ) );
        SelectionPhi.at ( iter )->Draw ( "E2SAME" );
    }
    LegendIntTrue->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas13IntTrue->SaveAs ( ( "MCPhi_Int_True"+SelectionLabel+"."+FileType ).c_str() );

    LegendMC->SetX1NDC ( 0.5 );
    LegendMC->SetY1NDC ( 0.6 );
    LegendMC->SetX2NDC ( 0.8 );
    LegendMC->SetY2NDC ( 0.8 );

    TCanvas *Canvas14 = new TCanvas ( "Energy", "Energy", 1400, 1000 );
    Canvas14->cd();
    SelectionEnergy.at ( 1 )->SetMaximum ( 1.2*GetMaximum ( SelectionEnergy ) );
    SelectionEnergy.at ( 1 )->SetMinimum ( 0.0 );
    SelectionEnergy.at ( 1 )->SetFillColor ( 9 );
    SelectionEnergy.at ( 1 )->Draw ( "E2" );
    SelectionEnergy.at ( 0 )->SetFillColor ( 8 );
    SelectionEnergy.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas14->SaveAs ( ( "MCEnergy"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas14a = new TCanvas ( "Momentum", "Momentum", 1400, 1000 );
    Canvas14a->cd();
    SelectionMomentum.at ( 1 )->SetMaximum ( 1.2*GetMaximum ( SelectionMomentum ) );
    SelectionMomentum.at ( 1 )->SetMinimum ( 0.0 );
    SelectionMomentum.at ( 1 )->SetFillColor ( 9 );
    SelectionMomentum.at ( 1 )->Draw ( "E2" );
    SelectionMomentum.at ( 0 )->SetFillColor ( 8 );
    SelectionMomentum.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas14a->SaveAs ( ( "MCMomentum"+SelectionLabel+"."+FileType ).c_str() );
    
    LegendMCSel->SetX1NDC ( 0.6 );
    LegendMCSel->SetY1NDC ( 0.7 );
    LegendMCSel->SetX2NDC ( 0.8 );
    LegendMCSel->SetY2NDC ( 0.75 );
    
    TCanvas *Canvas14aSel = new TCanvas ( "Momentum Sel", "Momentum Sel", 1400, 1000 );
    Canvas14aSel->cd();
    SelectionMomentum.at ( 2 )->SetMaximum ( 1.2*SelectionMomentum.at ( 2 )->GetBinContent(SelectionMomentum.at ( 2 )->GetMaximumBin()) );
    SelectionMomentum.at ( 2 )->SetMinimum ( 0.0 );
    SelectionMomentum.at ( 2 )->SetFillColor ( 46 );
    SelectionMomentum.at ( 2 )->Draw ( "E2" );
    LegendMCSel->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas14aSel->SaveAs ( ( "MCMomentumSel"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas14aInt = new TCanvas ( "Momentum Int", "Momentum Int", 1400, 1000 );
    Canvas14aInt->cd();
    SelectionMomentum.at ( 10 )->SetMaximum ( 1.1*SelectionMomentum.at ( 10 )->GetBinContent ( SelectionMomentum.at ( 10 )->GetMaximumBin() ) );
    SelectionMomentum.at ( 10 )->SetMinimum ( 0.0 );
    SelectionMomentum.at ( 10 )->SetFillColor ( ColorMap.at ( 0 ) );
    SelectionMomentum.at ( 10 )->Draw ( "E2" );
    for ( unsigned int iter = 11; iter < SelectionMomentum.size()-1; iter++ )
    {
        SelectionMomentum.at ( iter )->SetFillColor ( ColorMap.at ( iter-10 ) );
        SelectionMomentum.at ( iter )->Draw ( "E2SAME" );
    }
    LegendInt->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas14aInt->SaveAs ( ( "MCMomentum_Int"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas14aIntTrue = new TCanvas ( "Momentum Int Ture", "Momentum Int Ture", 1400, 1000 );
    Canvas14aIntTrue->cd();
    SelectionMomentum.at ( 4 )->SetMaximum ( 1.1*SelectionMomentum.at ( 4 )->GetBinContent ( SelectionMomentum.at ( 4 )->GetMaximumBin() ) );
    SelectionMomentum.at ( 4 )->SetMinimum ( 0.0 );
    SelectionMomentum.at ( 4 )->SetFillColor ( ColorMap.at ( 0 ) );
    SelectionMomentum.at ( 4 )->Draw ( "E2" );
    for ( unsigned int iter = 5; iter < 7; iter++ )
    {
        SelectionMomentum.at ( iter )->SetFillColor ( ColorMap.at ( iter-4 ) );
        SelectionMomentum.at ( iter )->Draw ( "E2SAME" );
    }
    LegendIntTrue->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
    Canvas14aIntTrue->SaveAs ( ( "MCMomentum_Int_True"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas15 = new TCanvas ( "OnBeam Minus OffBeam X Start & End Point ", "OnBeam Minus OffBeam X Start & End Point ", 1400, 1000 );
    Canvas15->cd();
    SelXTrackStartEnd.at ( 1 )->SetMaximum ( 1.5*GetMaximum ( SelXTrackStartEnd ) );
    SelXTrackStartEnd.at ( 1 )->SetMinimum ( 0.0 );
    SelXTrackStartEnd.at ( 1 )->SetFillColor ( 9 );
    SelXTrackStartEnd.at ( 1 )->Draw ( "E2" );
    SelXTrackStartEnd.at ( 0 )->SetFillColor ( 8 );
    SelXTrackStartEnd.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas15->SaveAs ( ( "MCXTrack"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas16 = new TCanvas ( "OnBeam Minus OffBeam Y Start & End Point ", "OnBeam Minus OffBeam Y Start & End Point ", 1400, 1000 );
    Canvas16->cd();
    SelYTrackStartEnd.at ( 1 )->SetMaximum ( 1.8*GetMaximum ( SelYTrackStartEnd ) );
    SelYTrackStartEnd.at ( 1 )->SetMinimum ( 0.0 );
    SelYTrackStartEnd.at ( 1 )->SetFillColor ( 9 );
    SelYTrackStartEnd.at ( 1 )->Draw ( "E2" );
    SelYTrackStartEnd.at ( 0 )->SetFillColor ( 8 );
    SelYTrackStartEnd.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas16->SaveAs ( ( "MCYTrack"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas17 = new TCanvas ( "OnBeam Minus OffBeam Z Start & End Point ", "OnBeam Minus OffBeam Z Start & End Point ", 1400, 1000 );
    Canvas17->cd();
    SelZTrackStartEnd.at ( 1 )->SetMaximum ( 1.5*GetMaximum ( SelZTrackStartEnd ) );
    SelZTrackStartEnd.at ( 1 )->SetMinimum ( 0.0 );
    SelZTrackStartEnd.at ( 1 )->SetFillColor ( 9 );
    SelZTrackStartEnd.at ( 1 )->Draw ( "E2" );
    SelZTrackStartEnd.at ( 0 )->SetFillColor ( 8 );
    SelZTrackStartEnd.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas17->SaveAs ( ( "MCZTrack"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas18 = new TCanvas ( "OnBeam Minus OffBeam X Vertex Postion", "OnBeam Minus OffBeam X Vertex Postion", 1400, 1000 );
    Canvas18->cd();
    SelXVtxPosition.at ( 1 )->SetMaximum ( 1.5*GetMaximum ( SelXVtxPosition ) );
    SelXVtxPosition.at ( 1 )->SetMinimum ( 0.0 );
    SelXVtxPosition.at ( 1 )->SetFillColor ( 9 );
    SelXVtxPosition.at ( 1 )->Draw ( "E2" );
    SelXVtxPosition.at ( 0 )->SetFillColor ( 8 );
    SelXVtxPosition.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas18->SaveAs ( ( "MCXVertex"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas19 = new TCanvas ( "OnBeam Minus OffBeam Y Vertex Postion", "OnBeam Minus OffBeam Y Vertex Postion", 1400, 1000 );
    Canvas19->cd();
    SelYVtxPosition.at ( 1 )->SetMaximum ( 1.8*GetMaximum ( SelYVtxPosition ) );
    SelYVtxPosition.at ( 1 )->SetMinimum ( 0.0 );
    SelYVtxPosition.at ( 1 )->SetFillColor ( 9 );
    SelYVtxPosition.at ( 1 )->Draw ( "E2" );
    SelYVtxPosition.at ( 0 )->SetFillColor ( 8 );
    SelYVtxPosition.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas19->SaveAs ( ( "MCYVertex"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas20 = new TCanvas ( "OnBeam Minus OffBeam Z Vertex Postion", "OnBeam Minus OffBeam Z Vertex Postion", 1400, 1000 );
    Canvas20->cd();
    SelZVtxPosition.at ( 1 )->SetMaximum ( 1.5*GetMaximum ( SelZVtxPosition ) );
    SelZVtxPosition.at ( 1 )->SetMinimum ( 0.0 );
    SelZVtxPosition.at ( 1 )->SetFillColor ( 9 );
    SelZVtxPosition.at ( 1 )->Draw ( "E2" );
    SelZVtxPosition.at ( 0 )->SetFillColor ( 8 );
    SelZVtxPosition.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    TextSimulation.Draw();
    TextSelection.Draw();
//     Canvas20->SaveAs ( ( "MCZVertex"+SelectionLabel+"."+FileType ).c_str() );
}

float GetMaximum ( const std::vector<TH1F*>& HistVector )
{
    float Maximum = 0;

    for ( unsigned int hist_no = 0; hist_no < 2; hist_no++ )
    {
        float TempMax = HistVector.at ( hist_no )->GetBinContent ( HistVector.at ( hist_no )->GetMaximumBin() );

        if ( TempMax > Maximum )
        {
            Maximum = TempMax;
        }
    }

    return Maximum;
}

void AddFirstTwoHistograms ( std::vector<TH1F*>& HistVector, float Weight )
{
    if ( HistVector.size() > 1 )
    {
        HistVector.at ( 0 )->Add ( HistVector.at ( 1 ), Weight );
        delete HistVector.at ( 1 );
        HistVector.erase ( HistVector.begin() +1 );
    }
    else
    {
        std::cout << "Histograms not added!" << std::endl;
    }
}

void AddFirstTwoHistograms2D ( std::vector<TH2F*>& HistVector, float Weight )
{
    if ( HistVector.size() > 1 )
    {
        HistVector.at ( 0 )->Add ( HistVector.at ( 1 ), Weight );
        delete HistVector.at ( 1 );
        HistVector.erase ( HistVector.begin() +1 );
    }
    else
    {
        std::cout << "Histograms not added!" << std::endl;
    }
}

float CalcLength ( const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2 )
{
    return sqrt ( pow ( x_1-x_2, 2 ) + pow ( y_1-y_2, 2 ) + pow ( z_1-z_2, 2 ) );
}

double FlashTrackDist ( double flash, double start, double end )
{
    if ( end >= start )
    {
        if ( flash < end && flash > start )
        {
            return 0;
        }
        else
        {
            return TMath::Min ( fabs ( flash-start ), fabs ( flash-end ) );
        }
    }
    else
    {
        if ( flash > end && flash < start )
        {
            return 0;
        }
        else
        {
            return TMath::Min ( fabs ( flash-start ), fabs ( flash-end ) );
        }
    }
}

bool inDeadRegion ( double y, double z )
{
    if ( ( y < ( 0.63*z+20 ) ) && ( y > ( 0.63*z-130 ) ) )
    {
        return true;
    }
    else if ( ( y < ( 0.63*z-185 ) ) && ( y > ( 0.63*z-232.3 ) ) )
    {
        return true;
    }
    else if ( ( y > ( -0.63*z+429.3 ) ) && ( y < ( -0.63*z+476.5 ) ) )
    {
        return true;
    }
    else if ( z > 700 && z < 750 )
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::vector<TSpline5> Systematics()
{
    // Number of columns in file
    unsigned short NumberOfColumns = 5;

    // Initialize data structure
    std::vector <std::vector<float>> BeamSystematics;
    BeamSystematics.resize ( NumberOfColumns );

    // Line and cell string for ifstream
    std::string FileLine;
    std::string Cell;

    // Open field systematic error file
    std::ifstream SysFile ( "bnb_sys_error_uboone.txt" );

    // check file
    if ( SysFile.bad() )
    {
        std::cout << "No such file or directory: " << "bnb_sys_error_uboone.txt" << std::endl;
        exit ( -1 );
    }

    // Loop over lines until files end
    while ( std::getline ( SysFile,FileLine ) )
    {
        // If not a header line
        if ( FileLine[0] != 'E' )
        {
            // Loop over all columns
            for ( unsigned column_no = 0; column_no < NumberOfColumns; column_no++ )
            {
                // Only read data if data stream works
                if ( SysFile >> Cell )
                {
                    // Fill systematic error data into data structure
                    BeamSystematics.at ( column_no ).push_back ( std::stof ( Cell ) );
                }
            } // End of column loop
        } // if not header
    } // line loop

    // Initialize Graph vector
    std::vector<TGraph*> GraphVector;

    // Fill graphs with beam systematic data
    for ( unsigned int entry_no = 1; entry_no < NumberOfColumns; entry_no++ )
    {
        GraphVector.push_back ( new TGraph ( BeamSystematics.at ( 0 ).size(),BeamSystematics.at ( 0 ).data(),BeamSystematics.at ( entry_no ).data() ) );
    }

    // Initialize spline vector
    std::vector<TSpline5> SplineVector;

    // Produce spline vector by fitting all graphs
    for ( const auto& Graph : GraphVector )
    {
        SplineVector.push_back ( TSpline5 ( "",Graph ) );
        delete Graph;
    }

    return SplineVector;
}

void AdjustSysError ( std::vector<TH1F*>& HistVector )
{
    for ( unsigned int bin_no = 1; bin_no < HistVector.back()->GetNbinsX() +1; bin_no++ )
    {
        HistVector.back()->SetBinError ( bin_no, HistVector.back()->GetBinContent ( bin_no ) - HistVector.at ( 1 )->GetBinContent ( bin_no ) + HistVector.at ( 1 )->GetBinError ( bin_no ) );
        HistVector.back()->SetBinContent ( bin_no, HistVector.at ( 1 )->GetBinContent ( bin_no ) );
    }
}

bool inFV ( double x, double y, double z )
{
    if ( x < ( FVx - borderx ) && ( x > borderx ) && ( y < ( FVy/2. - bordery ) ) && ( y > ( -FVy/2. + bordery ) ) && ( z < ( FVz - borderz ) ) && ( z > borderz ) )
    {
        return true;
    }
    else
    {
        return false;
    }
}
