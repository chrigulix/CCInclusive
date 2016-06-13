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
    
    std::vector<int> ColorVec = {3,4,5};

    TF1* SinTheta = new TF1 ( "const","sin(x)",0,3.142 );

    TLegend* LegendEfficiency = new TLegend ( 0.15,0.7,0.45,0.85 );
    LegendEfficiency->SetLineColorAlpha ( 0,0 );
    LegendEfficiency->SetLineStyle ( 0 );
    LegendEfficiency->SetFillStyle ( 0 );
    LegendEfficiency->SetMargin ( 0.2 );
    LegendEfficiency->SetTextFont ( 43 );
    LegendEfficiency->SetTextSize ( 30 );

    EfficiencyLabel.push_back ( "FV Efficiency" );
    EfficiencyLabel.push_back ( "Contained Efficiency" );

    TLegend* LegendMC = new TLegend ( 0.5,0.6,0.8,0.8 );
    LegendMC->SetLineColorAlpha ( 0,0 );
    LegendMC->SetLineStyle ( 0 );
    LegendMC->SetFillStyle ( 0 );
    LegendMC->SetMargin ( 0.2 );
    LegendMC->SetTextFont ( 43 );
    LegendMC->SetTextSize ( 30 );

    MCLabel.push_back ( "True Track in FV" );
    MCLabel.push_back ( "True Contained Tracks" );
    MCLabel.push_back ( "Selection on MC BNB+Cosmic with Stat. Error" );

    TLegend* LegendInt = new TLegend ( 0.6,0.69,0.89,0.89 );
    LegendInt->SetLineColorAlpha ( 0,0 );
    LegendInt->SetLineStyle ( 0 );
    LegendInt->SetFillStyle ( 0 );
    LegendInt->SetMargin ( 0.2 );
    LegendInt->SetTextFont ( 43 );
    LegendInt->SetTextSize ( 30 );

    IntLabel.push_back ( "QE Intercations" );
    IntLabel.push_back ( "RES Interactions" );
    IntLabel.push_back ( "DIS Interactions" );

    TLegend* LegendEffInt = new TLegend ( 0.65,0.65,0.85,0.85 );
    LegendEffInt->SetLineColorAlpha ( 0,0 );
    LegendEffInt->SetLineStyle ( 0 );
    LegendEffInt->SetFillStyle ( 0 );
    LegendEffInt->SetMargin ( 0.2 );
    LegendEffInt->SetTextFont ( 43 );
    LegendEffInt->SetTextSize ( 30 );

    EffIntLabel.push_back ( "QE Efficiency" );
    EffIntLabel.push_back ( "RES Efficiency" );
    EffIntLabel.push_back ( "DIS Efficiency" );
//     EffIntLabel.push_back ( "QE Contained Eff." );
//     EffIntLabel.push_back ( "RES Contained Eff." );
//     EffIntLabel.push_back ( "DIS Contained Eff." );

//     MCLabel.push_back ( "Selection MC BNB+Cosmic Sys. Error" );

    TLegend* FlashLabel = new TLegend ( 0.7,0.7,0.9,0.9 );

    GenLabel.push_back ( "Truth BNB Nu Cosmic in FV" );
    GenLabel.push_back ( "Truth BNB Nu Cosmic in FV Muon contained " );
    GenLabel.push_back ( "MC Prodgenie BNB Nu Cosmic" );
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

    std::vector<unsigned int> ColorMap = {28,42,30,38};

    ChainVec.push_back ( new TChain ( "anatree" ) );
    ChainVec.back() -> Add ( "/lheppc46/data/uBData/anatrees/Hist_MC_Truth_prodgenie_bnb_nu_cosmic_uboone_v05_08_00.root" );

    ChainVec.push_back ( new TChain ( "anatree" ) );
    ChainVec.back() -> Add ( ( "/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_prodgenie_bnb_nu_cosmic_uboone_v05_08_00"+ SelectionLabel +".root" ).c_str() );

    for ( const auto& Label : GenLabel )
    {
        SelectionTrackRange.push_back ( new TH1F ( ( "Track Range"+Label ).c_str(),"Track Range of Selected Track",NumberOfBins,0,1000 ) );
        SelectionTrackRange.back()->SetStats ( 0 );
        SelectionTrackRange.back()->GetXaxis()->SetTitle ( "Track Range [cm]" );
        SelectionTrackRange.back()->GetYaxis()->SetTitle ( "# Events" );

        SelectionTheta.push_back ( new TH1F ( ( "#theta-Angle"+Label ).c_str(),"#theta-Angle of Selected Track",NumberOfBins,0,3.142 ) );
        SelectionTheta.back()->SetStats ( 0 );
        SelectionTheta.back()->GetXaxis()->SetTitle ( "#theta [rad]" );
        SelectionTheta.back()->GetYaxis()->SetTitle ( "# Events" );
        SelectionTheta.back()->GetYaxis()->SetTitleOffset ( 1.3 );

        SelectionCosTheta.push_back ( new TH1F ( ( "cos#theta-Angle"+Label ).c_str(),"cos#theta of Selected Track",NumberOfBins,-1,1 ) );
        SelectionCosTheta.back()->SetStats ( 0 );
        SelectionCosTheta.back()->GetXaxis()->SetTitle ( "cos#theta [ ]" );
        SelectionCosTheta.back()->GetYaxis()->SetTitle ( "# Events}" );
        SelectionCosTheta.back()->GetYaxis()->SetTitleOffset ( 1.3 );

        SelectionPhi.push_back ( new TH1F ( ( "#phi-Angle"+Label ).c_str(),"#phi-Angle of Selected Track",NumberOfBins,-3.142,3.142 ) );
        SelectionPhi.back()->SetStats ( 0 );
        SelectionPhi.back()->GetXaxis()->SetTitle ( "#phi angle [rad]" );
        SelectionPhi.back()->GetYaxis()->SetTitle ( "# Events" );
        SelectionPhi.back()->GetYaxis()->SetTitleOffset ( 1.3 );

        SelectionEnergy.push_back ( new TH1F ( ( "Energy"+Label ).c_str(),"Energy of Selected Track",NumberOfBins,0,3 ) );
        SelectionEnergy.back()->SetStats ( 0 );
        SelectionEnergy.back()->GetXaxis()->SetTitle ( "Muon Kinetic Energy [MeV]" );
        SelectionEnergy.back()->GetYaxis()->SetTitle ( "# Events" );
        SelectionEnergy.back()->GetYaxis()->SetTitleOffset ( 1.3 );

        SelectionMomentum.push_back ( new TH1F ( ( "Momentum"+Label ).c_str(),"Momentum of Selected Track",NumberOfBins,0,3 ) );
        SelectionMomentum.back()->SetStats ( 0 );
        SelectionMomentum.back()->GetXaxis()->SetTitle ( "Muon Momentum [GeV/c]" );
        SelectionMomentum.back()->GetYaxis()->SetTitle ( "# Events" );
        SelectionMomentum.back()->GetYaxis()->SetTitleOffset ( 1.3 );

        SelXTrackStartEnd.push_back ( new TH1F ( ( "XTrack"+Label ).c_str(),"X Track Start & End Positions",NumberOfBins,0,256 ) );
        SelXTrackStartEnd.back()->SetStats ( 0 );
        SelXTrackStartEnd.back()->GetXaxis()->SetTitle ( "x [cm]" );
        SelXTrackStartEnd.back()->GetYaxis()->SetTitle ( "# Events" );
        SelXTrackStartEnd.back()->GetYaxis()->SetTitleOffset ( 1.3 );

        SelYTrackStartEnd.push_back ( new TH1F ( ( "YTrack"+Label ).c_str(),"Y Track Start & End Positions",NumberOfBins,-233/2,233/2 ) );
        SelYTrackStartEnd.back()->SetStats ( 0 );
        SelYTrackStartEnd.back()->GetXaxis()->SetTitle ( "y [cm]" );
        SelYTrackStartEnd.back()->GetYaxis()->SetTitle ( "# Events" );
        SelYTrackStartEnd.back()->GetYaxis()->SetTitleOffset ( 1.3 );

        SelZTrackStartEnd.push_back ( new TH1F ( ( "ZTrack"+Label ).c_str(),"Z Track Start & End Positions",NumberOfBins,0,1036.8 ) );
        SelZTrackStartEnd.back()->SetStats ( 0 );
        SelZTrackStartEnd.back()->GetXaxis()->SetTitle ( "z [cm]" );
        SelZTrackStartEnd.back()->GetYaxis()->SetTitle ( "# Events" );
        SelZTrackStartEnd.back()->GetYaxis()->SetTitleOffset ( 1.3 );

        SelXVtxPosition.push_back ( new TH1F ( ( "XVertex"+Label ).c_str(),"X Vertex Position",NumberOfBins,0,256 ) );
        SelXVtxPosition.back()->SetStats ( 0 );
        SelXVtxPosition.back()->GetXaxis()->SetTitle ( "x [cm]" );
        SelXVtxPosition.back()->GetYaxis()->SetTitle ( "# Events" );
        SelXVtxPosition.back()->GetYaxis()->SetTitleOffset ( 1.3 );

        SelYVtxPosition.push_back ( new TH1F ( ( "YVertex"+Label ).c_str(),"Y Vertex Position",NumberOfBins,-233/2,233/2 ) );
        SelYVtxPosition.back()->SetStats ( 0 );
        SelYVtxPosition.back()->GetXaxis()->SetTitle ( "y [cm]" );
        SelYVtxPosition.back()->GetYaxis()->SetTitle ( "# Events" );
        SelYVtxPosition.back()->GetYaxis()->SetTitleOffset ( 1.3 );

        SelZVtxPosition.push_back ( new TH1F ( ( "ZVertex"+Label ).c_str(),"Z Vertex Position",NumberOfBins,0,1036.8 ) );
        SelZVtxPosition.back()->SetStats ( 0 );
        SelZVtxPosition.back()->GetXaxis()->SetTitle ( "z [cm]" );
        SelZVtxPosition.back()->GetYaxis()->SetTitle ( "# Events" );
        SelZVtxPosition.back()->GetYaxis()->SetTitleOffset ( 1.3 );
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
    int CCNCFlag[10];
    int TruthMode[10];
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
        ChainVec.at ( file_no ) -> SetBranchAddress ( "ccnc_truth", CCNCFlag );
        ChainVec.at ( file_no ) -> SetBranchAddress ( "mode_truth", TruthMode );
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

            if ( file_no == 0 && MCTrkID > -1 && PDGTruth[MCTrkID] == 13 )
            {
                SelectionTrackRange.at ( 0 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                SelectionTheta.at ( 0 )->Fill ( MCTheta[MCTrkID] );
                SelectionCosTheta.at ( 0 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                SelectionPhi.at ( 0 )->Fill ( MCPhi[MCTrkID] );
                SelectionEnergy.at ( 0 )->Fill ( MCEnergy[MCTrkID] );
                SelectionMomentum.at ( 0 )->Fill ( TrueLeptonMomentum[0] );

                SelXTrackStartEnd.at ( 0 )->Fill ( XMCTrackStart[MCTrkID] );
                SelXTrackStartEnd.at ( 0 )->Fill ( XMCTrackEnd[MCTrkID] );
                SelYTrackStartEnd.at ( 0 )->Fill ( YMCTrackStart[MCTrkID] );
                SelYTrackStartEnd.at ( 0 )->Fill ( YMCTrackEnd[MCTrkID] );
                SelZTrackStartEnd.at ( 0 )->Fill ( ZMCTrackStart[MCTrkID] );
                SelZTrackStartEnd.at ( 0 )->Fill ( ZMCTrackEnd[MCTrkID] );
                SelXVtxPosition.at ( 0 )->Fill ( XnuVtxTruth[0] );
                SelYVtxPosition.at ( 0 )->Fill ( YnuVtxTruth[0] );
                SelZVtxPosition.at ( 0 )->Fill ( ZnuVtxTruth[0] );

                // True QE RES DIS
                if ( TruthMode[0] == 0 )
                {
                    nuQE++;
                    SelectionTrackRange.at ( 3 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                    SelectionTheta.at ( 3 )->Fill ( MCTheta[MCTrkID] );
                    SelectionCosTheta.at ( 3 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                    SelectionPhi.at ( 3 )->Fill ( MCPhi[MCTrkID] );
                    SelectionEnergy.at ( 3 )->Fill ( MCEnergy[MCTrkID] );
                    SelectionMomentum.at ( 3 )->Fill ( TrueLeptonMomentum[0] );

                    SelXTrackStartEnd.at ( 3 )->Fill ( XMCTrackStart[MCTrkID] );
                    SelXTrackStartEnd.at ( 3 )->Fill ( XMCTrackEnd[MCTrkID] );
                    SelYTrackStartEnd.at ( 3 )->Fill ( YMCTrackStart[MCTrkID] );
                    SelYTrackStartEnd.at ( 3 )->Fill ( YMCTrackEnd[MCTrkID] );
                    SelZTrackStartEnd.at ( 3 )->Fill ( ZMCTrackStart[MCTrkID] );
                    SelZTrackStartEnd.at ( 3 )->Fill ( ZMCTrackEnd[MCTrkID] );
                    SelXVtxPosition.at ( 3 )->Fill ( XnuVtxTruth[0] );
                    SelYVtxPosition.at ( 3 )->Fill ( YnuVtxTruth[0] );
                    SelZVtxPosition.at ( 3 )->Fill ( ZnuVtxTruth[0] );

                }
                else if ( TruthMode[0] == 1 )
                {
                    nuRES++;
                    SelectionTrackRange.at ( 4 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                    SelectionTheta.at ( 4 )->Fill ( MCTheta[MCTrkID] );
                    SelectionCosTheta.at ( 4 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                    SelectionPhi.at ( 4 )->Fill ( MCPhi[MCTrkID] );
                    SelectionEnergy.at ( 4 )->Fill ( MCEnergy[MCTrkID] );
                    SelectionMomentum.at ( 4 )->Fill ( TrueLeptonMomentum[0] );

                    SelXTrackStartEnd.at ( 4 )->Fill ( XMCTrackStart[MCTrkID] );
                    SelXTrackStartEnd.at ( 4 )->Fill ( XMCTrackEnd[MCTrkID] );
                    SelYTrackStartEnd.at ( 4 )->Fill ( YMCTrackStart[MCTrkID] );
                    SelYTrackStartEnd.at ( 4 )->Fill ( YMCTrackEnd[MCTrkID] );
                    SelZTrackStartEnd.at ( 4 )->Fill ( ZMCTrackStart[MCTrkID] );
                    SelZTrackStartEnd.at ( 4 )->Fill ( ZMCTrackEnd[MCTrkID] );
                    SelXVtxPosition.at ( 4 )->Fill ( XnuVtxTruth[0] );
                    SelYVtxPosition.at ( 4 )->Fill ( YnuVtxTruth[0] );
                    SelZVtxPosition.at ( 4 )->Fill ( ZnuVtxTruth[0] );

                }
                else if ( TruthMode[0] == 2 )
                {
                    nuDIS++;
                    SelectionTrackRange.at ( 5 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                    SelectionTheta.at ( 5 )->Fill ( MCTheta[MCTrkID] );
                    SelectionCosTheta.at ( 5 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                    SelectionPhi.at ( 5 )->Fill ( MCPhi[MCTrkID] );
                    SelectionEnergy.at ( 5 )->Fill ( MCEnergy[MCTrkID] );
                    SelectionMomentum.at ( 5 )->Fill ( TrueLeptonMomentum[0] );

                    SelXTrackStartEnd.at ( 5 )->Fill ( XMCTrackStart[MCTrkID] );
                    SelXTrackStartEnd.at ( 5 )->Fill ( XMCTrackEnd[MCTrkID] );
                    SelYTrackStartEnd.at ( 5 )->Fill ( YMCTrackStart[MCTrkID] );
                    SelYTrackStartEnd.at ( 5 )->Fill ( YMCTrackEnd[MCTrkID] );
                    SelZTrackStartEnd.at ( 5 )->Fill ( ZMCTrackStart[MCTrkID] );
                    SelZTrackStartEnd.at ( 5 )->Fill ( ZMCTrackEnd[MCTrkID] );
                    SelXVtxPosition.at ( 5 )->Fill ( XnuVtxTruth[0] );
                    SelYVtxPosition.at ( 5 )->Fill ( YnuVtxTruth[0] );
                    SelZVtxPosition.at ( 5 )->Fill ( ZnuVtxTruth[0] );
                }


                if ( inFV ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID] ) && inFV ( XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) )
                {
                    ContainedTracks++;

                    SelectionTrackRange.at ( 1 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                    SelectionTheta.at ( 1 )->Fill ( MCTheta[MCTrkID] );
                    SelectionCosTheta.at ( 1 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                    SelectionPhi.at ( 1 )->Fill ( MCPhi[MCTrkID] );
                    SelectionEnergy.at ( 1 )->Fill ( MCEnergy[MCTrkID] );
                    SelectionMomentum.at ( 1 )->Fill ( TrueLeptonMomentum[0] );

                    SelXTrackStartEnd.at ( 1 )->Fill ( XMCTrackStart[MCTrkID] );
                    SelXTrackStartEnd.at ( 1 )->Fill ( XMCTrackEnd[MCTrkID] );
                    SelYTrackStartEnd.at ( 1 )->Fill ( YMCTrackStart[MCTrkID] );
                    SelYTrackStartEnd.at ( 1 )->Fill ( YMCTrackEnd[MCTrkID] );
                    SelZTrackStartEnd.at ( 1 )->Fill ( ZMCTrackStart[MCTrkID] );
                    SelZTrackStartEnd.at ( 1 )->Fill ( ZMCTrackEnd[MCTrkID] );
                    SelXVtxPosition.at ( 1 )->Fill ( XnuVtxTruth[0] );
                    SelYVtxPosition.at ( 1 )->Fill ( YnuVtxTruth[0] );
                    SelZVtxPosition.at ( 1 )->Fill ( ZnuVtxTruth[0] );

                    // True Contained QE RES DIS
                    if ( TruthMode[0] == 0 )
                    {
                        nuQE++;
                        SelectionTrackRange.at ( 6 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                        SelectionTheta.at ( 6 )->Fill ( MCTheta[MCTrkID] );
                        SelectionCosTheta.at ( 6 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                        SelectionPhi.at ( 6 )->Fill ( MCPhi[MCTrkID] );
                        SelectionEnergy.at ( 6 )->Fill ( MCEnergy[MCTrkID] );
                        SelectionMomentum.at ( 6 )->Fill ( TrueLeptonMomentum[0] );

                        SelXTrackStartEnd.at ( 6 )->Fill ( XMCTrackStart[MCTrkID] );
                        SelXTrackStartEnd.at ( 6 )->Fill ( XMCTrackEnd[MCTrkID] );
                        SelYTrackStartEnd.at ( 6 )->Fill ( YMCTrackStart[MCTrkID] );
                        SelYTrackStartEnd.at ( 6 )->Fill ( YMCTrackEnd[MCTrkID] );
                        SelZTrackStartEnd.at ( 6 )->Fill ( ZMCTrackStart[MCTrkID] );
                        SelZTrackStartEnd.at ( 6 )->Fill ( ZMCTrackEnd[MCTrkID] );
                        SelXVtxPosition.at ( 6 )->Fill ( XnuVtxTruth[0] );
                        SelYVtxPosition.at ( 6 )->Fill ( YnuVtxTruth[0] );
                        SelZVtxPosition.at ( 6 )->Fill ( ZnuVtxTruth[0] );

                    }
                    else if ( TruthMode[0] == 1 )
                    {
                        nuRES++;
                        SelectionTrackRange.at ( 7 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                        SelectionTheta.at ( 7 )->Fill ( MCTheta[MCTrkID] );
                        SelectionCosTheta.at ( 7 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                        SelectionPhi.at ( 7 )->Fill ( MCPhi[MCTrkID] );
                        SelectionEnergy.at ( 7 )->Fill ( MCEnergy[MCTrkID] );
                        SelectionMomentum.at ( 7 )->Fill ( TrueLeptonMomentum[0] );

                        SelXTrackStartEnd.at ( 7 )->Fill ( XMCTrackStart[MCTrkID] );
                        SelXTrackStartEnd.at ( 7 )->Fill ( XMCTrackEnd[MCTrkID] );
                        SelYTrackStartEnd.at ( 7 )->Fill ( YMCTrackStart[MCTrkID] );
                        SelYTrackStartEnd.at ( 7 )->Fill ( YMCTrackEnd[MCTrkID] );
                        SelZTrackStartEnd.at ( 7 )->Fill ( ZMCTrackStart[MCTrkID] );
                        SelZTrackStartEnd.at ( 7 )->Fill ( ZMCTrackEnd[MCTrkID] );
                        SelXVtxPosition.at ( 7 )->Fill ( XnuVtxTruth[0] );
                        SelYVtxPosition.at ( 7 )->Fill ( YnuVtxTruth[0] );
                        SelZVtxPosition.at ( 7 )->Fill ( ZnuVtxTruth[0] );

                    }
                    else if ( TruthMode[0] == 2 )
                    {
                        nuDIS++;
                        SelectionTrackRange.at ( 8 )->Fill ( CalcLength ( XMCTrackStart[MCTrkID],YMCTrackStart[MCTrkID],ZMCTrackStart[MCTrkID],XMCTrackEnd[MCTrkID],YMCTrackEnd[MCTrkID],ZMCTrackEnd[MCTrkID] ) );
                        SelectionTheta.at ( 8 )->Fill ( MCTheta[MCTrkID] );
                        SelectionCosTheta.at ( 8 )->Fill ( cos ( MCTheta[MCTrkID] ) );
                        SelectionPhi.at ( 8 )->Fill ( MCPhi[MCTrkID] );
                        SelectionEnergy.at ( 8 )->Fill ( MCEnergy[MCTrkID] );
                        SelectionMomentum.at ( 8 )->Fill ( TrueLeptonMomentum[0] );

                        SelXTrackStartEnd.at ( 8 )->Fill ( XMCTrackStart[MCTrkID] );
                        SelXTrackStartEnd.at ( 8 )->Fill ( XMCTrackEnd[MCTrkID] );
                        SelYTrackStartEnd.at ( 8 )->Fill ( YMCTrackStart[MCTrkID] );
                        SelYTrackStartEnd.at ( 8 )->Fill ( YMCTrackEnd[MCTrkID] );
                        SelZTrackStartEnd.at ( 8 )->Fill ( ZMCTrackStart[MCTrkID] );
                        SelZTrackStartEnd.at ( 8 )->Fill ( ZMCTrackEnd[MCTrkID] );
                        SelXVtxPosition.at ( 8 )->Fill ( XnuVtxTruth[0] );
                        SelYVtxPosition.at ( 8 )->Fill ( YnuVtxTruth[0] );
                        SelZVtxPosition.at ( 8 )->Fill ( ZnuVtxTruth[0] );
                    }
                }
            }
            else if ( file_no == 1 && TrkID > -1 && TrkOrigin[TrkID][TrkBestPlane[TrkID]] == 1 && PDGTruth[MCTrkID] == 13 )
            {
                SelectionTrackRange.at ( 2 )->Fill ( CalcLength ( XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID] ) );
                SelectionTheta.at ( 2 )->Fill ( TrackTheta[TrkID] );
                SelectionCosTheta.at ( 2 )->Fill ( cos ( TrackTheta[TrkID] ) );
                SelectionPhi.at ( 2 )->Fill ( TrackPhi[TrkID] );
                SelectionEnergy.at ( 2 )->Fill ( KineticEnergy[TrkID][2]/1000 );
                SelectionMomentum.at ( 2 )->Fill ( TrackMomentum[TrkID] );

                SelXTrackStartEnd.at ( 2 )->Fill ( XTrackStart[TrkID] );
                SelXTrackStartEnd.at ( 2 )->Fill ( XTrackEnd[TrkID] );
                SelYTrackStartEnd.at ( 2 )->Fill ( YTrackStart[TrkID] );
                SelYTrackStartEnd.at ( 2 )->Fill ( YTrackEnd[TrkID] );
                SelZTrackStartEnd.at ( 2 )->Fill ( ZTrackStart[TrkID] );
                SelZTrackStartEnd.at ( 2 )->Fill ( ZTrackEnd[TrkID] );
                SelXVtxPosition.at ( 2 )->Fill ( XVertexPosition[VtxID] );
                SelYVtxPosition.at ( 2 )->Fill ( YVertexPosition[VtxID] );
                SelZVtxPosition.at ( 2 )->Fill ( ZVertexPosition[VtxID] );
                

                // Selection QE RES DIS
                if ( TruthMode[0] == 0 )
                {
                    nuQE++;
                    SelectionTrackRange.at ( 9 )->Fill ( CalcLength ( XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID] ) );
                    SelectionTheta.at ( 9 )->Fill ( TrackTheta[TrkID] );
                    SelectionCosTheta.at ( 9 )->Fill ( cos ( TrackTheta[TrkID] ) );
                    SelectionPhi.at ( 9 )->Fill ( TrackPhi[TrkID] );
                    SelectionEnergy.at ( 9 )->Fill ( KineticEnergy[TrkID][2]/1000 );
                    SelectionMomentum.at ( 9 )->Fill ( TrackMomentum[TrkID] );

                    SelXTrackStartEnd.at ( 9 )->Fill ( XTrackStart[TrkID] );
                    SelXTrackStartEnd.at ( 9 )->Fill ( XTrackEnd[TrkID] );
                    SelYTrackStartEnd.at ( 9 )->Fill ( YTrackStart[TrkID] );
                    SelYTrackStartEnd.at ( 9 )->Fill ( YTrackEnd[TrkID] );
                    SelZTrackStartEnd.at ( 9 )->Fill ( ZTrackStart[TrkID] );
                    SelZTrackStartEnd.at ( 9 )->Fill ( ZTrackEnd[TrkID] );
                    SelXVtxPosition.at ( 9 )->Fill ( XVertexPosition[VtxID] );
                    SelYVtxPosition.at ( 9 )->Fill ( YVertexPosition[VtxID] );
                    SelZVtxPosition.at ( 9 )->Fill ( ZVertexPosition[VtxID] );

                }
                else if ( TruthMode[0] == 1 )
                {
                    nuRES++;
                    SelectionTrackRange.at ( 10 )->Fill ( CalcLength ( XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID] ) );
                    SelectionTheta.at ( 10 )->Fill ( TrackTheta[TrkID] );
                    SelectionCosTheta.at ( 10 )->Fill ( cos ( TrackTheta[TrkID] ) );
                    SelectionPhi.at ( 10 )->Fill ( TrackPhi[TrkID] );
                    SelectionEnergy.at ( 10 )->Fill ( KineticEnergy[TrkID][2]/1000 );
                    SelectionMomentum.at ( 10 )->Fill ( TrackMomentum[TrkID] );

                    SelXTrackStartEnd.at ( 10 )->Fill ( XTrackStart[TrkID] );
                    SelXTrackStartEnd.at ( 10 )->Fill ( XTrackEnd[TrkID] );
                    SelYTrackStartEnd.at ( 10 )->Fill ( YTrackStart[TrkID] );
                    SelYTrackStartEnd.at ( 10 )->Fill ( YTrackEnd[TrkID] );
                    SelZTrackStartEnd.at ( 10 )->Fill ( ZTrackStart[TrkID] );
                    SelZTrackStartEnd.at ( 10 )->Fill ( ZTrackEnd[TrkID] );
                    SelXVtxPosition.at ( 10 )->Fill ( XVertexPosition[VtxID] );
                    SelYVtxPosition.at ( 10 )->Fill ( YVertexPosition[VtxID] );
                    SelZVtxPosition.at ( 10 )->Fill ( ZVertexPosition[VtxID] );

                }
                else if ( TruthMode[0] == 2 )
                {
                    nuDIS++;
                    SelectionTrackRange.at ( 11 )->Fill ( CalcLength ( XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID] ) );
                    SelectionTheta.at ( 11 )->Fill ( TrackTheta[TrkID] );
                    SelectionCosTheta.at ( 11 )->Fill ( cos ( TrackTheta[TrkID] ) );
                    SelectionPhi.at ( 11 )->Fill ( TrackPhi[TrkID] );
                    SelectionEnergy.at ( 11 )->Fill ( KineticEnergy[TrkID][2]/1000 );
                    SelectionMomentum.at ( 11 )->Fill ( TrackMomentum[TrkID] );

                    SelXTrackStartEnd.at ( 11 )->Fill ( XTrackStart[TrkID] );
                    SelXTrackStartEnd.at ( 11 )->Fill ( XTrackEnd[TrkID] );
                    SelYTrackStartEnd.at ( 11 )->Fill ( YTrackStart[TrkID] );
                    SelYTrackStartEnd.at ( 11 )->Fill ( YTrackEnd[TrkID] );
                    SelZTrackStartEnd.at ( 11 )->Fill ( ZTrackStart[TrkID] );
                    SelZTrackStartEnd.at ( 11 )->Fill ( ZTrackEnd[TrkID] );
                    SelXVtxPosition.at ( 11 )->Fill ( XVertexPosition[VtxID] );
                    SelYVtxPosition.at ( 11 )->Fill ( YVertexPosition[VtxID] );
                    SelZVtxPosition.at ( 11 )->Fill ( ZVertexPosition[VtxID] );
                }


                // Fill systematic errors independet of CC or NC
                SelectionTrackRange.back()->Fill ( CalcLength ( XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID] ),1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelectionTheta.back()->Fill ( TrackTheta[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelectionCosTheta.back()->Fill ( cos ( TrackTheta[TrkID] ),1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelectionPhi.back()->Fill ( TrackPhi[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelectionEnergy.back()->Fill ( KineticEnergy[TrkID][2]/1000,1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelectionMomentum.back()->Fill ( TrackMomentum[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );

                SelXTrackStartEnd.back()->Fill ( XTrackStart[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelXTrackStartEnd.back()->Fill ( XTrackEnd[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelYTrackStartEnd.back()->Fill ( YTrackStart[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelYTrackStartEnd.back()->Fill ( YTrackEnd[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelZTrackStartEnd.back()->Fill ( ZTrackStart[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelZTrackStartEnd.back()->Fill ( ZTrackEnd[TrkID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelXVtxPosition.back()->Fill ( XVertexPosition[VtxID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelYVtxPosition.back()->Fill ( YVertexPosition[VtxID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
                SelZVtxPosition.back()->Fill ( ZVertexPosition[VtxID],1+SystematicErrors.at ( 0 ).Eval ( NuEnergyTruth[0] ) );
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
        std::cout << TEfficiency::CheckConsistency(*SelectionEnergy.at(2), *SelectionEnergy.at(eff_no),"w") << " " << TEfficiency::CheckBinning(*SelectionEnergy.at(2), *SelectionEnergy.at(eff_no)) << std::endl; 
        EffTrackRange.push_back ( new TEfficiency ( *SelectionTrackRange.at ( 2 ),*SelectionTrackRange.at ( eff_no ) ) );
        EffEnergy.push_back ( new TEfficiency ( *SelectionEnergy.at ( 2 ),*SelectionEnergy.at ( eff_no ) ) );
        EffMomentum.push_back ( new TEfficiency ( *SelectionMomentum.at ( 2 ),*SelectionMomentum.at ( eff_no ) ) );
        EffTheta.push_back ( new TEfficiency ( *SelectionTheta.at ( 2 ),*SelectionTheta.at ( eff_no ) ) );
        EffCosTheta.push_back ( new TEfficiency ( *SelectionCosTheta.at ( 2 ),*SelectionCosTheta.at ( eff_no ) ) );
        EffPhi.push_back ( new TEfficiency ( *SelectionPhi.at ( 2 ),*SelectionPhi.at ( eff_no ) ) );

        EffXTrackStartEnd.push_back ( new TEfficiency ( *SelXTrackStartEnd.at ( 2 ),*SelXTrackStartEnd.at ( eff_no ) ) );
        EffYTrackStartEnd.push_back ( new TEfficiency ( *SelYTrackStartEnd.at ( 2 ),*SelYTrackStartEnd.at ( eff_no ) ) );
        EffZTrackStartEnd.push_back ( new TEfficiency ( *SelZTrackStartEnd.at ( 2 ),*SelZTrackStartEnd.at ( eff_no ) ) );

        EffXVtxPosition.push_back ( new TEfficiency ( *SelXVtxPosition.at ( 2 ),*SelXVtxPosition.at ( eff_no ) ) );
        EffYVtxPosition.push_back ( new TEfficiency ( *SelYVtxPosition.at ( 2 ),*SelYVtxPosition.at ( eff_no ) ) );
        EffZVtxPosition.push_back ( new TEfficiency ( *SelZVtxPosition.at ( 2 ),*SelZVtxPosition.at ( eff_no ) ) );
    }
    
    for ( unsigned int eff_no = 0; eff_no < EffIntLabel.size(); eff_no++ )
    {
        EffTrackRange.push_back ( new TEfficiency ( *SelectionTrackRange.at ( 9+eff_no ),*SelectionTrackRange.at ( 2+eff_no ) ) );
        EffEnergy.push_back ( new TEfficiency ( *SelectionEnergy.at ( 9+eff_no ),*SelectionEnergy.at ( 2+eff_no ) ) );
        EffMomentum.push_back ( new TEfficiency ( *SelectionMomentum.at ( 9+eff_no ),*SelectionMomentum.at ( 2+eff_no ) ) );
        EffTheta.push_back ( new TEfficiency ( *SelectionTheta.at ( 9+eff_no ),*SelectionTheta.at ( 2+eff_no ) ) );
        EffCosTheta.push_back ( new TEfficiency ( *SelectionCosTheta.at ( 9+eff_no ),*SelectionCosTheta.at ( 2+eff_no ) ) );
        EffPhi.push_back ( new TEfficiency ( *SelectionPhi.at ( 9+eff_no ),*SelectionPhi.at ( 2+eff_no ) ) );

        EffXTrackStartEnd.push_back ( new TEfficiency ( *SelXTrackStartEnd.at ( 9+eff_no ),*SelXTrackStartEnd.at ( 2+eff_no ) ) );
        EffYTrackStartEnd.push_back ( new TEfficiency ( *SelYTrackStartEnd.at ( 9+eff_no ),*SelYTrackStartEnd.at ( 2+eff_no ) ) );
        EffZTrackStartEnd.push_back ( new TEfficiency ( *SelZTrackStartEnd.at ( 9+eff_no ),*SelZTrackStartEnd.at ( 2+eff_no ) ) );

        EffXVtxPosition.push_back ( new TEfficiency ( *SelXVtxPosition.at ( 9+eff_no ),*SelXVtxPosition.at ( 2+eff_no ) ) );
        EffYVtxPosition.push_back ( new TEfficiency ( *SelYVtxPosition.at ( 9+eff_no ),*SelYVtxPosition.at ( 2+eff_no ) ) );
        EffZVtxPosition.push_back ( new TEfficiency ( *SelZVtxPosition.at ( 9+eff_no ),*SelZVtxPosition.at ( 2+eff_no ) ) );
    }
    
    LegendEffInt->AddEntry ( EffTrackRange.at ( 2 ), ( EffIntLabel.at ( 0 ) ).c_str(),"lep" );
    LegendEffInt->AddEntry ( EffTrackRange.at ( 3 ), ( EffIntLabel.at ( 1 ) ).c_str(),"lep" );
    LegendEffInt->AddEntry ( EffTrackRange.at ( 4 ), ( EffIntLabel.at ( 2 ) ).c_str(),"lep" );

    LegendEfficiency->AddEntry ( EffTrackRange.at ( 1 ), ( EfficiencyLabel.at ( 1 ) ).c_str(),"lep" );
    LegendEfficiency->AddEntry ( EffTrackRange.at ( 0 ), ( EfficiencyLabel.at ( 0 ) ).c_str(),"lep" );

    TCanvas *Canvas1 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Track Range", "Efficiency OnBeam Minus OffBeam Track Range", 1400, 1000 );
    TMultiGraph *MGTrackRange = new TMultiGraph();
    EffTrackRange.at ( 0 )->SetLineWidth ( 2 );
    EffTrackRange.at ( 0 )->SetLineColor ( 8 );
    EffTrackRange.at ( 1 )->SetLineWidth ( 2 );
    EffTrackRange.at ( 1 )->SetLineColor (9);
    MGTrackRange->Add ( EffTrackRange.at ( 0 )->CreateGraph() );
    MGTrackRange->Add ( EffTrackRange.at ( 1 )->CreateGraph() );
    Canvas1->cd();
    MGTrackRange->Draw ( "AP" );
    MGTrackRange->GetXaxis()->SetTitle("Track Range [cm]");
    MGTrackRange->GetYaxis()->SetTitle("Efficiency [ ]");
    LegendEfficiency->Draw();
    Canvas1->SaveAs ( ( "EffMCRange"+SelectionLabel+"."+FileType ).c_str() );
    
    TCanvas *Canvas1Int = new TCanvas ( "Interaction Efficiency OnBeam Minus OffBeam Track Range", "Interaction Efficiency OnBeam Minus OffBeam Track Range", 1400, 1000 );
    TMultiGraph *MGTrackRangeInt = new TMultiGraph();
    EffTrackRange.at ( 2 )->SetLineWidth ( 2 );
    EffTrackRange.at ( 2 )->SetLineColor ( 46 );
    EffTrackRange.at ( 3 )->SetLineWidth ( 2 );
    EffTrackRange.at ( 3 )->SetLineColor ( 42 );
    EffTrackRange.at ( 4 )->SetLineWidth ( 2 );
    EffTrackRange.at ( 4 )->SetLineColor ( 30 );
    MGTrackRangeInt->Add ( EffTrackRange.at ( 2 )->CreateGraph() );
    MGTrackRangeInt->Add ( EffTrackRange.at ( 3 )->CreateGraph() );
    MGTrackRangeInt->Add ( EffTrackRange.at ( 4 )->CreateGraph() );
    Canvas1Int->cd();
    MGTrackRangeInt->Draw ( "AP" );
    MGTrackRangeInt->GetXaxis()->SetTitle("Track Range [cm]");
    MGTrackRangeInt->GetYaxis()->SetTitle("Efficiency [ ]");
    LegendEffInt->Draw();
    Canvas1Int->SaveAs ( ( "EffIntMCRange"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas2 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Theta-Angle", "Efficiency OnBeam Minus OffBeam Theta-Angle", 1400, 1000 );
    TMultiGraph *MGTheta = new TMultiGraph();
    EffTheta.at ( 0 )->SetLineWidth ( 2 );
    EffTheta.at ( 0 )->SetLineColor ( 8 );
    EffTheta.at ( 1 )->SetLineWidth ( 2 );
    EffTheta.at ( 1 )->SetLineColor (9);
    MGTheta->Add ( EffTheta.at ( 0 )->CreateGraph() );
    MGTheta->Add ( EffTheta.at ( 1 )->CreateGraph() );
    Canvas2->cd();
    MGTheta->Draw ( "AP" );
    MGTheta->GetXaxis()->SetTitle("#theta-Angle [rad]");
    MGTheta->GetYaxis()->SetTitle("Efficiency [ ]");
    LegendEfficiency->Draw();
    Canvas2->SaveAs ( ( "EffMCTheta"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas2b = new TCanvas ( "Efficiency OnBeam Minus OffBeam Cos Theta-Angle", "Efficiency OnBeam Minus OffBeam Cos Theta-Angle", 1400, 1000 );
    TMultiGraph *MGCosTheta = new TMultiGraph();
    EffCosTheta.at ( 0 )->SetLineWidth ( 2 );
    EffCosTheta.at ( 0 )->SetLineColor ( 8 );
    EffCosTheta.at ( 1 )->SetLineWidth ( 2 );
    EffCosTheta.at ( 1 )->SetLineColor (9);
    MGCosTheta->Add ( EffCosTheta.at ( 0 )->CreateGraph() );
    MGCosTheta->Add ( EffCosTheta.at ( 1 )->CreateGraph() );
    Canvas2b->cd();
    MGCosTheta->Draw ( "AP" );
    MGCosTheta->GetXaxis()->SetTitle("cos #theta [ ]");
    MGCosTheta->GetYaxis()->SetTitle("Efficiency [ ]");
    LegendEfficiency->Draw();
    Canvas2b->SaveAs ( ( "EffMCCosTheta"+SelectionLabel+"."+FileType ).c_str() );
    
    TCanvas *Canvas2bInt = new TCanvas ( "Interaction Efficiency OnBeam Minus OffBeam Cos Theta-Angle", "Interaction Efficiency OnBeam Minus OffBeam Cos Theta-Angle", 1400, 1000 );
    TMultiGraph *MGCosThetaInt = new TMultiGraph();
    EffCosTheta.at ( 2 )->SetLineWidth ( 2 );
    EffCosTheta.at ( 2 )->SetLineColor ( 46 );
    EffCosTheta.at ( 3 )->SetLineWidth ( 2 );
    EffCosTheta.at ( 3 )->SetLineColor ( 42 );
    EffCosTheta.at ( 4 )->SetLineWidth ( 2 );
    EffCosTheta.at ( 4 )->SetLineColor ( 30 );
    MGCosThetaInt->Add ( EffCosTheta.at ( 2 )->CreateGraph() );
    MGCosThetaInt->Add ( EffCosTheta.at ( 3 )->CreateGraph() );
    MGCosThetaInt->Add ( EffCosTheta.at ( 4 )->CreateGraph() );
    Canvas2bInt->cd();
    MGCosThetaInt->Draw ( "AP" );
    MGCosThetaInt->GetXaxis()->SetTitle("cos #theta [ ]");
    MGCosThetaInt->GetYaxis()->SetTitle("Efficiency [ ]");
    LegendEffInt->Draw();
    Canvas2bInt->SaveAs ( ( "EffIntMCCosTheta"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas4 = new TCanvas ( "Efficiency Energy", "Efficiency Energy", 1400, 1000 );
    TMultiGraph *MGEnergy = new TMultiGraph();
    EffEnergy.at ( 0 )->SetLineWidth ( 2 );
    EffEnergy.at ( 0 )->SetLineColor ( 8 );
    EffEnergy.at ( 1 )->SetLineWidth ( 2 );
    EffEnergy.at ( 1 )->SetLineColor (9);
    MGEnergy->Add ( EffEnergy.at ( 0 )->CreateGraph() );
    MGEnergy->Add ( EffEnergy.at ( 1 )->CreateGraph() );
    Canvas4->cd();
    MGEnergy->Draw ( "AP" );
    MGEnergy->GetXaxis()->SetTitle("Muon Energy [GeV]");
    MGEnergy->GetYaxis()->SetTitle("Efficiency [ ]");
    LegendEfficiency->Draw();
    Canvas4->SaveAs ( ( "EffMCEnergy"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas4a = new TCanvas ( "Efficiency Momentum", "Efficiency Momentum", 1400, 1000 );
    TMultiGraph *MGMomentum = new TMultiGraph();
    EffMomentum.at ( 0 )->SetLineWidth ( 2 );
    EffMomentum.at ( 0 )->SetLineColor ( 8 );
    EffMomentum.at ( 1 )->SetLineWidth ( 2 );
    EffMomentum.at ( 1 )->SetLineColor (9);
    MGMomentum->Add ( EffMomentum.at ( 0 )->CreateGraph() );
    MGMomentum->Add ( EffMomentum.at ( 1 )->CreateGraph() );
    Canvas4a->cd();
    MGMomentum->Draw ( "AP" );
    MGMomentum->GetXaxis()->SetTitle("Muon Momentum [GeV/c]");
    MGMomentum->GetYaxis()->SetTitle("Efficiency [ ]");
    LegendEfficiency->Draw();
    Canvas4a->SaveAs ( ( "EffMCMomentum"+SelectionLabel+"."+FileType ).c_str() );
    
    TCanvas *Canvas4aInt = new TCanvas ( "Interaction Efficiency Momentum", "Interaction Efficiency Momentum", 1400, 1000 );
    TMultiGraph *MGMomentumInt = new TMultiGraph();
    EffMomentum.at ( 2 )->SetLineWidth ( 2 );
    EffMomentum.at ( 2 )->SetLineColor ( 46 );
    EffMomentum.at ( 3 )->SetLineWidth ( 2 );
    EffMomentum.at ( 3 )->SetLineColor ( 42 );
    EffMomentum.at ( 4 )->SetLineWidth ( 2 );
    EffMomentum.at ( 4 )->SetLineColor ( 30 );
    MGMomentumInt->Add ( EffMomentum.at ( 2 )->CreateGraph() );
    MGMomentumInt->Add ( EffMomentum.at ( 3 )->CreateGraph() );
    MGMomentumInt->Add ( EffMomentum.at ( 4 )->CreateGraph() );
    Canvas4aInt->cd();
    MGMomentumInt->Draw ( "AP" );
    MGMomentumInt->GetXaxis()->SetTitle("Muon Momentum [GeV/c]");
    MGMomentumInt->GetYaxis()->SetTitle("Efficiency [ ]");
    LegendEffInt->Draw();
    Canvas4aInt->SaveAs ( ( "EffIntMCMomentum"+SelectionLabel+"."+FileType ).c_str() );

    LegendEfficiency->SetX1NDC ( 0.55 );
    LegendEfficiency->SetY1NDC ( 0.35 );
    LegendEfficiency->SetX2NDC ( 0.85 );
    LegendEfficiency->SetY2NDC ( 0.55 );

    TCanvas *Canvas3 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Phi-Angle", "Efficiency OnBeam Minus OffBeam Phi-Angle", 1400, 1000 );
    TMultiGraph *MGPhi = new TMultiGraph();
    EffPhi.at ( 0 )->SetLineWidth ( 2 );
    EffPhi.at ( 0 )->SetLineColor ( 8 );
    EffPhi.at ( 1 )->SetLineWidth ( 2 );
    EffPhi.at ( 1 )->SetLineColor (9);
    MGPhi->Add ( EffPhi.at ( 0 )->CreateGraph() );
    MGPhi->Add ( EffPhi.at ( 1 )->CreateGraph() );
    Canvas3->cd();
    MGPhi->Draw ( "AP" );
    MGPhi->GetXaxis()->SetTitle("#phi-Angle [rad]");
    MGPhi->GetYaxis()->SetTitle("Efficiency [ ]");
    LegendEfficiency->Draw();
    Canvas3->SaveAs ( ( "EffMCPhi"+SelectionLabel+"."+FileType ).c_str() );
    
    LegendEffInt->SetX1NDC ( 0.65 );
    LegendEffInt->SetY1NDC ( 0.35 );
    LegendEffInt->SetX2NDC ( 0.85 );
    LegendEffInt->SetY2NDC ( 0.55 );
    
    TCanvas *Canvas3Int = new TCanvas ( "Interaction Efficiency OnBeam Minus OffBeam Phi-Angle", "Interaction Efficiency OnBeam Minus OffBeam Phi-Angle", 1400, 1000 );
    TMultiGraph *MGPhiInt = new TMultiGraph();
    EffPhi.at ( 2 )->SetLineWidth ( 2 );
    EffPhi.at ( 2 )->SetLineColor ( 46 );
    EffPhi.at ( 3 )->SetLineWidth ( 2 );
    EffPhi.at ( 3 )->SetLineColor ( 42 );
    EffPhi.at ( 4 )->SetLineWidth ( 2 );
    EffPhi.at ( 4 )->SetLineColor ( 30 );
    MGPhiInt->Add ( EffPhi.at ( 2 )->CreateGraph() );
    MGPhiInt->Add ( EffPhi.at ( 3 )->CreateGraph() );
    MGPhiInt->Add ( EffPhi.at ( 4 )->CreateGraph() );
    Canvas3Int->cd();
    MGPhiInt->Draw ( "AP" );
    MGPhiInt->GetXaxis()->SetTitle("#phi-Angle [rad]");
    MGPhiInt->GetYaxis()->SetTitle("Efficiency [ ]");
    LegendEffInt->Draw();
    Canvas3Int->SaveAs ( ( "EffIntMCPhi"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas5 = new TCanvas ( "Efficiency OnBeam Minus OffBeam X Start & End Point ", "Efficiency OnBeam Minus OffBeam X Start & End Point ", 1400, 1000 );
    TMultiGraph *MGXStart = new TMultiGraph();
    Canvas5->cd();
//     EffXTrackStartEnd.at ( 1 )->SetMaximum ( 1.5*GetMaximum(EffXTrackStartEnd) );
//     EffXTrackStartEnd.at ( 1 )->SetMinimum ( 0.0 );
    EffXTrackStartEnd.at ( 1 )->SetLineColor (9);
    EffXTrackStartEnd.at ( 1 )->Draw ( "A" );
    EffXTrackStartEnd.at ( 0 )->SetLineColor ( 8 );
    EffXTrackStartEnd.at ( 0 )->Draw ( "SAME" );
    LegendEfficiency->Draw();
    Canvas5->SaveAs ( ( "EffMCXTrack"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas6 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Y Start & End Point ", "Efficiency OnBeam Minus OffBeam Y Start & End Point ", 1400, 1000 );
    TMultiGraph *MGYStart = new TMultiGraph();
    Canvas6->cd();
//     EffYTrackStartEnd.at ( 1 )->SetMaximum ( 1.8*GetMaximum(EffYTrackStartEnd) );
//     EffYTrackStartEnd.at ( 1 )->SetMinimum ( 0.0 );
    EffYTrackStartEnd.at ( 1 )->SetLineColor (9);
    EffYTrackStartEnd.at ( 1 )->Draw ( "A" );
    EffYTrackStartEnd.at ( 0 )->SetLineColor ( 8 );
    EffYTrackStartEnd.at ( 0 )->Draw ( "SAME" );
    LegendEfficiency->Draw();
    Canvas6->SaveAs ( ( "EffMCYTrack"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas7 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Z Start & End Point ", "Efficiency OnBeam Minus OffBeam Z Start & End Point ", 1400, 1000 );
    TMultiGraph *MGZStart = new TMultiGraph();
    Canvas7->cd();
//     EffZTrackStartEnd.at ( 1 )->SetMaximum ( 1.5*GetMaximum(EffZTrackStartEnd) );
//     EffZTrackStartEnd.at ( 1 )->SetMinimum ( 0.0 );
    EffZTrackStartEnd.at ( 1 )->SetLineColor (9);
    EffZTrackStartEnd.at ( 1 )->Draw ( "A" );
    EffZTrackStartEnd.at ( 0 )->SetLineColor ( 8 );
    EffZTrackStartEnd.at ( 0 )->Draw ( "SAME" );
//     LegendMC->Draw();
    Canvas7->SaveAs ( ( "EffMCZTrack"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas8 = new TCanvas ( "Efficiency OnBeam Minus OffBeam X Vertex Postion", "Efficiency OnBeam Minus OffBeam X Vertex Postion", 1400, 1000 );
    TMultiGraph *MGXVertex = new TMultiGraph();
    Canvas8->cd();
//     EffXVtxPosition.at ( 1 )->SetMaximum ( 1.5*GetMaximum(EffXVtxPosition) );
//     EffXVtxPosition.at ( 1 )->SetMinimum ( 0.0 );
    EffXVtxPosition.at ( 1 )->SetLineColor (9);
    EffXVtxPosition.at ( 1 )->Draw ( "A" );
    EffXVtxPosition.at ( 0 )->SetLineColor ( 8 );
    EffXVtxPosition.at ( 0 )->Draw ( "SAME" );
//     LegendMC->Draw();
    Canvas8->SaveAs ( ( "EffMCXVertex"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas9 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Y Vertex Postion", "Efficiency OnBeam Minus OffBeam Y Vertex Postion", 1400, 1000 );
    TMultiGraph *MGYVertex = new TMultiGraph();
    Canvas9->cd();
//     EffYVtxPosition.at ( 1 )->SetMaximum ( 1.8*GetMaximum(EffYVtxPosition) );
//     EffYVtxPosition.at ( 1 )->SetMinimum ( 0.0 );
    EffYVtxPosition.at ( 1 )->SetLineColor (9);
    EffYVtxPosition.at ( 1 )->Draw ( "A" );
    EffYVtxPosition.at ( 0 )->SetLineColor ( 8 );
    EffYVtxPosition.at ( 0 )->Draw ( "SAME" );
//     LegendMC->Draw();
    Canvas9->SaveAs ( ( "EffMCYVertex"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas10 = new TCanvas ( "Efficiency OnBeam Minus OffBeam Z Vertex Postion", "Efficiency OnBeam Minus OffBeam Z Vertex Postion", 1400, 1000 );
    TMultiGraph *MGZVertex = new TMultiGraph();
    Canvas10->cd();
//     EffZVtxPosition.at ( 1 )->SetMaximum ( 1.5*GetMaximum(EffZVtxPosition) );
//     EffZVtxPosition.at ( 1 )->SetMinimum ( 0.0 );
    EffZVtxPosition.at ( 1 )->SetLineColor (9);
    EffZVtxPosition.at ( 1 )->Draw ( "A" );
    EffZVtxPosition.at ( 0 )->SetLineColor ( 8 );
    EffZVtxPosition.at ( 0 )->Draw ( "SAME" );
//     LegendMC->Draw();
    Canvas10->SaveAs ( ( "EffMCZVertex"+SelectionLabel+"."+FileType ).c_str() );


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

    LegendMC->AddEntry ( SelectionTrackRange.at ( 0 ), ( MCLabel.at ( 0 ) ).c_str(),"f" );
    LegendMC->AddEntry ( SelectionTrackRange.at ( 1 ), ( MCLabel.at ( 1 ) ).c_str(),"f" );
//     LegendMC->AddEntry ( SelectionTrackRange.at ( 2 ), ( MCLabel.at ( 2 ) ).c_str(),"f" );
//     for ( unsigned int bgrhist_no = 0; bgrhist_no < BgrLabel.size(); bgrhist_no++ )
//     {
//         LegendMC->AddEntry( BgrTrackRange.at(bgrhist_no), (BgrLabel.at(bgrhist_no)).c_str(),"f" );
//     }

    TCanvas *Canvas11 = new TCanvas ( "OnBeam Minus OffBeam Track Range", "OnBeam Minus OffBeam Track Range", 1400, 1000 );
    Canvas11->cd();
    SelectionTrackRange.at ( 1 )->SetMaximum ( 1.1*GetMaximum ( SelectionTrackRange ) );
    SelectionTrackRange.at ( 1 )->SetMinimum ( 0.0 );
    SelectionTrackRange.at ( 1 )->SetFillColor (9);
    SelectionTrackRange.at ( 1 )->Draw ( "E2" );
    SelectionTrackRange.at ( 0 )->SetFillColor ( 8 );
    SelectionTrackRange.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas11->SaveAs ( ( "MCRange"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas11Int = new TCanvas ( "OnBeam Minus OffBeam Track Range Int", "OnBeam Minus OffBeam Track Range Int", 1400, 1000 );
    Canvas11Int->cd();
    SelectionTrackRange.at ( 9 )->SetMaximum ( 1.1*SelectionTrackRange.at(9)->GetBinContent(SelectionTrackRange.at(9)->GetMaximumBin()) );
    SelectionTrackRange.at ( 9 )->SetMinimum ( 0.0 );
    SelectionTrackRange.at ( 9 )->SetFillColor ( ColorVec.at(0) );
    SelectionTrackRange.at ( 9 )->Draw ( "E2" );
    for ( unsigned int iter = 10; iter < SelectionTrackRange.size()-1; iter++ )
    {
        SelectionTrackRange.at ( iter )->SetFillColor ( ColorVec.at(iter-10) );
        SelectionTrackRange.at ( iter )->Draw ( "E2SAME" );
    }
//     LegendMC->Draw();
    Canvas11Int->SaveAs ( ( "MCRange_Int"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas11IntTrue = new TCanvas ( "OnBeam Minus OffBeam Track Range Int True", "OnBeam Minus OffBeam Track Range Int Ture", 1400, 1000 );
    Canvas11IntTrue->cd();
    SelectionTrackRange.at ( 3 )->SetMaximum ( 1.1*SelectionTrackRange.at(3)->GetBinContent(SelectionTrackRange.at(3)->GetMaximumBin()) );
    SelectionTrackRange.at ( 3 )->SetMinimum ( 0.0 );
    SelectionTrackRange.at ( 3 )->SetFillColor ( ColorVec.at(0) );
    SelectionTrackRange.at ( 3 )->Draw ( "E2" );
    for ( unsigned int iter = 4; iter < 6; iter++ )
    {
        SelectionTrackRange.at ( iter )->SetFillColor ( ColorVec.at(iter-5) );
        SelectionTrackRange.at ( iter )->Draw ( "E2SAME" );
    }
//     LegendMC->Draw();
    Canvas11Int->SaveAs ( ( "MCRange_Int_True"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas12 = new TCanvas ( "OnBeam Minus OffBeam Theta-Angle", "OnBeam Minus OffBeam Theta-Angle", 1400, 1000 );
    Canvas12->cd();
    SelectionTheta.at ( 1 )->SetMaximum ( 1.2*GetMaximum ( SelectionTheta ) );
    SelectionTheta.at ( 1 )->SetMinimum ( 0.0 );
    SelectionTheta.at ( 1 )->SetFillColor (9);
    SelectionTheta.at ( 1 )->Draw ( "E2" );
    SelectionTheta.at ( 0 )->SetFillColor ( 8 );
    SelectionTheta.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas12->SaveAs ( ( "MCTheta"+SelectionLabel+"."+FileType ).c_str() );

    for ( auto& ThetaHistogram : SelectionTheta )
    {
        ThetaHistogram->Divide ( SinTheta,1. );
    }

    TCanvas *Canvas12a = new TCanvas ( "OnBeam Minus OffBeam Theta-Angle Omega", "OnBeam Minus OffBeam Theta-Angle Omega", 1400, 1000 );
    Canvas12a->cd();
    SelectionTheta.at ( 1 )->SetMaximum ( 1.2*GetMaximum ( SelectionTheta ) );
    SelectionTheta.at ( 1 )->SetMinimum ( 0.0 );
    SelectionTheta.at ( 1 )->GetYaxis()->SetTitle ( "# Events" );
    SelectionTheta.at ( 1 )->SetFillColor (9);
    SelectionTheta.at ( 1 )->Draw ( "E2" );
    SelectionTheta.at ( 0 )->SetFillColor ( 8 );
    SelectionTheta.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas12a->SaveAs ( ( "MCThetaOmega"+SelectionLabel+"."+FileType ).c_str() );

    LegendMC->SetX1NDC ( 0.2 );
    LegendMC->SetY1NDC ( 0.6 );
    LegendMC->SetX2NDC ( 0.5 );
    LegendMC->SetY2NDC ( 0.8 );

    TCanvas *Canvas12b = new TCanvas ( "OnBeam Minus OffBeam Cos Theta-Angle", "OnBeam Minus OffBeam Cos Theta-Angle", 1400, 1000 );
    Canvas12b->cd();
    SelectionCosTheta.at ( 1 )->SetMaximum ( 1.2*GetMaximum ( SelectionCosTheta ) );
    SelectionCosTheta.at ( 1 )->SetMinimum ( 0.0 );
    SelectionCosTheta.at ( 1 )->SetFillColor (9);
    SelectionCosTheta.at ( 1 )->Draw ( "E2" );
    SelectionCosTheta.at ( 0 )->SetFillColor ( 8 );
    SelectionCosTheta.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas12b->SaveAs ( ( "MCCosTheta"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas12bInt = new TCanvas ( "OnBeam Minus OffBeam Cos Theta-Angle Int", "OnBeam Minus OffBeam Cos Theta-Angle Int", 1400, 1000 );
    Canvas12bInt->cd();
    SelectionCosTheta.at ( 9 )->SetMaximum ( 1.1*SelectionCosTheta.at(9)->GetBinContent(SelectionCosTheta.at(9)->GetMaximumBin()) );
    SelectionCosTheta.at ( 9 )->SetMinimum ( 0.0 );
    SelectionCosTheta.at ( 9 )->SetFillColor ( ColorVec.at(0) );
    SelectionCosTheta.at ( 9 )->Draw ( "E2" );
    for ( unsigned int iter = 10; iter < SelectionCosTheta.size()-1; iter++ )
    {
        SelectionCosTheta.at ( iter )->SetFillColor ( ColorVec.at(iter-11) );
        SelectionCosTheta.at ( iter )->Draw ( "E2SAME" );
    }
//     LegendMC->Draw();
    Canvas12bInt->SaveAs ( ( "MCCosTheta_Int"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas12bIntTrue = new TCanvas ( "OnBeam Minus OffBeam Cos Theta-Angle Int Ture", "OnBeam Minus OffBeam Cos Theta-Angle Int Ture", 1400, 1000 );
    Canvas12bIntTrue->cd();
    SelectionCosTheta.at ( 3 )->SetMaximum ( 1.1*SelectionCosTheta.at(3)->GetBinContent(SelectionCosTheta.at(3)->GetMaximumBin()) );
    SelectionCosTheta.at ( 3 )->SetMinimum ( 0.0 );
    SelectionCosTheta.at ( 3 )->SetFillColor ( ColorVec.at(0) );
    SelectionCosTheta.at ( 3 )->Draw ( "E2" );
    for ( unsigned int iter = 4; iter < 6; iter++ )
    {
        SelectionCosTheta.at ( iter )->SetFillColor ( ColorVec.at(iter-5) );
        SelectionCosTheta.at ( iter )->Draw ( "E2SAME" );
    }
    Canvas12bInt->SaveAs ( ( "MCCosTheta_Int_True"+SelectionLabel+"."+FileType ).c_str() );

    LegendMC->SetX1NDC ( 0.5 );
    LegendMC->SetY1NDC ( 0.4 );
    LegendMC->SetX2NDC ( 0.8 );
    LegendMC->SetY2NDC ( 0.6 );

    TCanvas *Canvas13 = new TCanvas ( "OnBeam Minus OffBeam Phi-Angle", "OnBeam Minus OffBeam Phi-Angle", 1400, 1000 );
    Canvas13->cd();
    SelectionPhi.at ( 1 )->SetMaximum ( 1.2*GetMaximum ( SelectionPhi ) );
    SelectionPhi.at ( 1 )->SetMinimum ( 0.0 );
    SelectionPhi.at ( 1 )->SetFillColor (9);
    SelectionPhi.at ( 1 )->Draw ( "E2" );
    SelectionPhi.at ( 0 )->SetFillColor ( 8 );
    SelectionPhi.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas13->SaveAs ( ( "MCPhi"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas13Int = new TCanvas ( "OnBeam Minus OffBeam Phi-Angle Int", "OnBeam Minus OffBeam Phi-Angle Int", 1400, 1000 );
    Canvas13Int->cd();
    SelectionPhi.at ( 9 )->SetMaximum ( 1.1*SelectionPhi.at(9)->GetBinContent(SelectionPhi.at(9)->GetMaximumBin()) );
    SelectionPhi.at ( 9 )->SetMinimum ( 0.0 );
    SelectionPhi.at ( 9 )->SetFillColor ( ColorVec.at(0) );
    SelectionPhi.at ( 9 )->Draw ( "E2" );
    for ( unsigned int iter = 10; iter < SelectionPhi.size()-1; iter++ )
    {
        SelectionPhi.at ( iter )->SetFillColor ( ColorVec.at(iter-11) );
        SelectionPhi.at ( iter )->Draw ( "E2SAME" );
    }
//     LegendMC->Draw();
    Canvas13Int->SaveAs ( ( "MCPhi_Int"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas13IntTrue = new TCanvas ( "OnBeam Minus OffBeam Phi-Angle Int Ture", "OnBeam Minus OffBeam Phi-Angle Int Ture", 1400, 1000 );
    Canvas13IntTrue->cd();
    SelectionPhi.at ( 3 )->SetMaximum ( 1.1*SelectionPhi.at(3)->GetBinContent(SelectionPhi.at(3)->GetMaximumBin()) );
    SelectionPhi.at ( 3 )->SetMinimum ( 0.0 );
    SelectionPhi.at ( 3 )->SetFillColor ( ColorVec.at(0) );
    SelectionPhi.at ( 3 )->Draw ( "E2" );
    for ( unsigned int iter = 4; iter < 6; iter++ )
    {
        SelectionPhi.at ( iter )->SetFillColor ( ColorVec.at(iter-5) );
        SelectionPhi.at ( iter )->Draw ( "E2SAME" );
    }
    Canvas13Int->SaveAs ( ( "MCPhi_Int_True"+SelectionLabel+"."+FileType ).c_str() );

    LegendMC->SetX1NDC ( 0.5 );
    LegendMC->SetY1NDC ( 0.6 );
    LegendMC->SetX2NDC ( 0.8 );
    LegendMC->SetY2NDC ( 0.8 );

    TCanvas *Canvas14 = new TCanvas ( "Energy", "Energy", 1400, 1000 );
    Canvas14->cd();
    SelectionEnergy.at ( 1 )->SetMaximum ( 1.2*GetMaximum ( SelectionEnergy ) );
    SelectionEnergy.at ( 1 )->SetMinimum ( 0.0 );
    SelectionEnergy.at ( 1 )->SetFillColor (9);
    SelectionEnergy.at ( 1 )->Draw ( "E2" );
    SelectionEnergy.at ( 0 )->SetFillColor ( 8 );
    SelectionEnergy.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas14->SaveAs ( ( "MCEnergy"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas14a = new TCanvas ( "Momentum", "Momentum", 1400, 1000 );
    Canvas14a->cd();
    SelectionMomentum.at ( 1 )->SetMaximum ( 1.2*GetMaximum ( SelectionMomentum ) );
    SelectionMomentum.at ( 1 )->SetMinimum ( 0.0 );
    SelectionMomentum.at ( 1 )->SetFillColor (9);
    SelectionMomentum.at ( 1 )->Draw ( "E2" );
    SelectionMomentum.at ( 0 )->SetFillColor ( 8 );
    SelectionMomentum.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas14a->SaveAs ( ( "MCMomentum"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas14aInt = new TCanvas ( "Momentum Int", "Momentum Int", 1400, 1000 );
    Canvas14aInt->cd();
    SelectionMomentum.at ( 9 )->SetMaximum ( 1.1*SelectionMomentum.at(9)->GetBinContent(SelectionMomentum.at(9)->GetMaximumBin()) );
    SelectionMomentum.at ( 9 )->SetMinimum ( 0.0 );
    SelectionMomentum.at ( 9 )->SetFillColor ( ColorVec.at(0) );
    SelectionMomentum.at ( 9 )->Draw ( "E2" );
    for ( unsigned int iter = 10; iter < SelectionMomentum.size()-1; iter++ )
    {
        SelectionMomentum.at ( iter )->SetFillColor ( ColorVec.at(iter-11) );
        SelectionMomentum.at ( iter )->Draw ( "E2SAME" );
    }
//     LegendMC->Draw();
    Canvas14aInt->SaveAs ( ( "MCPhi_Int"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas14aIntTrue = new TCanvas ( "Momentum Int Ture", "Momentum Int Ture", 1400, 1000 );
    Canvas14aIntTrue->cd();
    SelectionMomentum.at ( 3 )->SetMaximum ( 1.1*SelectionMomentum.at(3)->GetBinContent(SelectionMomentum.at(3)->GetMaximumBin()) );
    SelectionMomentum.at ( 3 )->SetMinimum ( 0.0 );
    SelectionMomentum.at ( 3 )->SetFillColor ( ColorVec.at(0) );
    SelectionMomentum.at ( 3 )->Draw ( "E2" );
    for ( unsigned int iter = 4; iter < 6; iter++ )
    {
        SelectionMomentum.at ( iter )->SetFillColor ( ColorVec.at(iter-5) );
        SelectionMomentum.at ( iter )->Draw ( "E2SAME" );
    }
    Canvas14aInt->SaveAs ( ( "MCPhi_Int_True"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas15 = new TCanvas ( "OnBeam Minus OffBeam X Start & End Point ", "OnBeam Minus OffBeam X Start & End Point ", 1400, 1000 );
    Canvas15->cd();
    SelXTrackStartEnd.at ( 1 )->SetMaximum ( 1.5*GetMaximum ( SelXTrackStartEnd ) );
    SelXTrackStartEnd.at ( 1 )->SetMinimum ( 0.0 );
    SelXTrackStartEnd.at ( 1 )->SetFillColor (9);
    SelXTrackStartEnd.at ( 1 )->Draw ( "E2" );
    SelXTrackStartEnd.at ( 0 )->SetFillColor ( 8 );
    SelXTrackStartEnd.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas15->SaveAs ( ( "MCXTrack"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas16 = new TCanvas ( "OnBeam Minus OffBeam Y Start & End Point ", "OnBeam Minus OffBeam Y Start & End Point ", 1400, 1000 );
    Canvas16->cd();
    SelYTrackStartEnd.at ( 1 )->SetMaximum ( 1.8*GetMaximum ( SelYTrackStartEnd ) );
    SelYTrackStartEnd.at ( 1 )->SetMinimum ( 0.0 );
    SelYTrackStartEnd.at ( 1 )->SetFillColor (9);
    SelYTrackStartEnd.at ( 1 )->Draw ( "E2" );
    SelYTrackStartEnd.at ( 0 )->SetFillColor ( 8 );
    SelYTrackStartEnd.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas16->SaveAs ( ( "MCYTrack"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas17 = new TCanvas ( "OnBeam Minus OffBeam Z Start & End Point ", "OnBeam Minus OffBeam Z Start & End Point ", 1400, 1000 );
    Canvas17->cd();
    SelZTrackStartEnd.at ( 1 )->SetMaximum ( 1.5*GetMaximum ( SelZTrackStartEnd ) );
    SelZTrackStartEnd.at ( 1 )->SetMinimum ( 0.0 );
    SelZTrackStartEnd.at ( 1 )->SetFillColor (9);
    SelZTrackStartEnd.at ( 1 )->Draw ( "E2" );
    SelZTrackStartEnd.at ( 0 )->SetFillColor ( 8 );
    SelZTrackStartEnd.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas17->SaveAs ( ( "MCZTrack"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas18 = new TCanvas ( "OnBeam Minus OffBeam X Vertex Postion", "OnBeam Minus OffBeam X Vertex Postion", 1400, 1000 );
    Canvas18->cd();
    SelXVtxPosition.at ( 1 )->SetMaximum ( 1.5*GetMaximum ( SelXVtxPosition ) );
    SelXVtxPosition.at ( 1 )->SetMinimum ( 0.0 );
    SelXVtxPosition.at ( 1 )->SetFillColor (9);
    SelXVtxPosition.at ( 1 )->Draw ( "E2" );
    SelXVtxPosition.at ( 0 )->SetFillColor ( 8 );
    SelXVtxPosition.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas18->SaveAs ( ( "MCXVertex"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas19 = new TCanvas ( "OnBeam Minus OffBeam Y Vertex Postion", "OnBeam Minus OffBeam Y Vertex Postion", 1400, 1000 );
    Canvas19->cd();
    SelYVtxPosition.at ( 1 )->SetMaximum ( 1.8*GetMaximum ( SelYVtxPosition ) );
    SelYVtxPosition.at ( 1 )->SetMinimum ( 0.0 );
    SelYVtxPosition.at ( 1 )->SetFillColor (9);
    SelYVtxPosition.at ( 1 )->Draw ( "E2" );
    SelYVtxPosition.at ( 0 )->SetFillColor ( 8 );
    SelYVtxPosition.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas19->SaveAs ( ( "MCYVertex"+SelectionLabel+"."+FileType ).c_str() );

    TCanvas *Canvas20 = new TCanvas ( "OnBeam Minus OffBeam Z Vertex Postion", "OnBeam Minus OffBeam Z Vertex Postion", 1400, 1000 );
    Canvas20->cd();
    SelZVtxPosition.at ( 1 )->SetMaximum ( 1.5*GetMaximum ( SelZVtxPosition ) );
    SelZVtxPosition.at ( 1 )->SetMinimum ( 0.0 );
    SelZVtxPosition.at ( 1 )->SetFillColor (9);
    SelZVtxPosition.at ( 1 )->Draw ( "E2" );
    SelZVtxPosition.at ( 0 )->SetFillColor ( 8 );
    SelZVtxPosition.at ( 0 )->Draw ( "E2SAME" );
    LegendMC->Draw();
    Canvas20->SaveAs ( ( "MCZVertex"+SelectionLabel+"."+FileType ).c_str() );
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
