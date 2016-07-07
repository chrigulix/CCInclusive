// ROOT logon file with microboone style by Christoph Rudolf von Rohr
// Download this file to your root working folder. Also check the default 
// "MicroBooNE Preliminary" or "MicroBooNE Simulation, Preliminary" text box
// in the last line of this code.

// Some headers
#include <TStyle.h>
#include <TPad.h>
#include <TAxis.h>
#include <TPavesText.h>

// Root logon function
void rootlogon()
{
    int Font = 42;
    double LabelSize = 0.04;
    double TitleSize = 0.05;
    
    // Set text font and text size
    gStyle->SetTextFont( Font );
    gStyle->SetTextSize(LabelSize);
    
    // Set axis label title fonts
    gStyle->SetLabelFont( Font,"x");
    gStyle->SetTitleFont( Font,"x");
    gStyle->SetLabelFont( Font,"y");
    gStyle->SetTitleFont( Font,"y");
    gStyle->SetLabelFont( Font,"z");
    gStyle->SetTitleFont( Font,"z");

    // Set axis label title fonts
    gStyle->SetLabelSize(LabelSize,"x");
    gStyle->SetTitleSize(TitleSize,"x");
    gStyle->SetLabelSize(LabelSize,"y");
    gStyle->SetTitleSize(TitleSize,"y");
    gStyle->SetLabelSize(LabelSize,"z");
    gStyle->SetTitleSize(TitleSize,"z");
    
    // Set axis title offsets
    gStyle->SetTitleOffset(1.,"x");
    gStyle->SetTitleOffset(1.04,"y");

    // Supress title and statistics in canvas
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    // Set default color
    int Color = 0; // WHITE
    gStyle->SetFrameBorderMode(Color);
    gStyle->SetFrameFillColor(Color);
    gStyle->SetCanvasBorderMode(Color);
    gStyle->SetCanvasColor(Color);
    gStyle->SetPadBorderMode(Color);
    gStyle->SetPadColor(Color);
    gStyle->SetStatColor(Color);

    // Create ticks on all four sides of a plot
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // Set data point marker style and size
    gStyle->SetMarkerStyle(1);
    gStyle->SetMarkerSize(1.2);
    
    // Set line width for different cases
    gStyle->SetHistLineWidth(2.);
    gStyle->SetLineWidth(2.);
    gStyle->SetLineStyleString(2,"[12 12]");// postscript dashes
    
    // The following specifies the "MicroBooNE Simulation, Preliminary" text box
    TPaveText* TextSimulation = new TPaveText(0.5,0.92,0.9,0.96,"nbNDC");
    TextSimulation->AddText("MicroBooNE Simulation, Preliminary");
    TextSimulation->SetTextSize(0.05);
    TextSimulation->SetTextColor(12);
    TextSimulation->SetLineColorAlpha(0,0);
    TextSimulation->SetFillColorAlpha(0,0);
    TextSimulation->SetTextAlign(33);
    
    // The following specifies the "MicroBooNE Preliminary" text box
    TPaveText* TextPreliminary = new TPaveText(0.6,0.92,0.9,0.96,"nbNDC");
    TextPreliminary->AddText("MicroBooNE Preliminary");
    TextPreliminary->SetTextSize(0.05);
    TextPreliminary->SetTextColor(12);
    TextPreliminary->SetLineColorAlpha(0,0);
    TextPreliminary->SetFillColorAlpha(0,0);
    TextPreliminary->SetTextAlign(33);
    
    // The following specifies the "MicroBooNE" text box
    TPaveText* TextUBooNE = new TPaveText(0.6,0.92,0.9,0.96,"nbNDC");
    TextUBooNE->AddText("MicroBooNE");
    TextUBooNE->SetTextSize(0.05);
    TextUBooNE->SetTextColor(12);
    TextUBooNE->SetLineColorAlpha(0,0);
    TextUBooNE->SetFillColorAlpha(0,0);
    TextUBooNE->SetTextAlign(33);
    
//     TPad* Pad = (TPad*) gPad;
    
    // Choose wich text box you want to have by activating either of these
//     gPad->GetListOfPrimitives()->Add(TextSimulation->Clone());
//     gPad->GetListOfPrimitives()->AddLast(TextPreliminary);

}
//by CRvR