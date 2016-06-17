// ROOT logon file with microboone style

#include <TStyle.h>
#include <TPad.h>
#include <TAxis.h>


void rootlogon()
{
// gROOT->LoadMacro("AtlasStyle.C");
// AtlasStyle();
    int font = 42; // Serif
    double tsize = 0.04;
    gStyle->SetTextFont(font);

    gStyle->SetTextSize(tsize);
    gStyle->SetLabelFont(font,"x");
    gStyle->SetTitleFont(font,"x");
    gStyle->SetLabelFont(font,"y");
    gStyle->SetTitleFont(font,"y");
    gStyle->SetLabelFont(font,"z");
    gStyle->SetTitleFont(font,"z");

    gStyle->SetLabelSize(tsize,"x");
    gStyle->SetTitleSize(tsize,"x");
    gStyle->SetLabelSize(tsize,"y");
    gStyle->SetTitleSize(tsize,"y");
    gStyle->SetLabelSize(tsize,"z");
    gStyle->SetTitleSize(tsize,"z");
    gStyle->SetTitleOffset(1.3,"x");
    gStyle->SetTitleOffset(1.1,"y");

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    int icol = 0; // WHITE
    gStyle->SetFrameBorderMode(icol);
    gStyle->SetFrameFillColor(icol);
    gStyle->SetCanvasBorderMode(icol);
    gStyle->SetCanvasColor(icol);
    gStyle->SetPadBorderMode(icol);
    gStyle->SetPadColor(icol);
    gStyle->SetStatColor(icol);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetMarkerStyle(1);
    gStyle->SetMarkerSize(1.2);
    gStyle->SetHistLineWidth(2.);
    gStyle->SetLineWidth(2.);
    gStyle->SetLineStyleString(2,"[12 12]");// postscript dashes
}
