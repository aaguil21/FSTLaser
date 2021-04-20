#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>

#include <string>

void plotTrackRes()
{
  gStyle->SetOptStat(1101);
  gStyle->SetOptFit(1001);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.2); gStyle->SetStatH(0.3);

  string inputfile = "./IstTracking.root";
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mXResidual = (TH1F*)File_InPut->Get("h_mXResidual");
  TH1F *h_mYResidual = (TH1F*)File_InPut->Get("h_mYResidual");

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,1000,500);
  c_play->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_play->cd(i_pad+1)->SetLeftMargin(0.15);
    c_play->cd(i_pad+1)->SetBottomMargin(0.15);
    c_play->cd(i_pad+1)->SetTicks(1,1);
    c_play->cd(i_pad+1)->SetGrid(0,0);
  }

  c_play->cd(1);
  h_mXResidual->GetXaxis()->SetTitle("x-residual (mm)");
  h_mXResidual->GetXaxis()->SetTitleSize(0.06);
  h_mXResidual->GetYaxis()->SetTitle("No. Tracks");
  h_mXResidual->GetYaxis()->SetTitleSize(0.06);
  h_mXResidual->Draw();
  h_mXResidual->Fit("gaus");

  c_play->cd(2);
  h_mYResidual->GetXaxis()->SetTitle("y-residual (mm)");
  h_mYResidual->GetXaxis()->SetTitleSize(0.06);
  h_mYResidual->GetYaxis()->SetTitle("No. Tracks");
  h_mYResidual->GetYaxis()->SetTitleSize(0.06);
  h_mYResidual->Draw();
  h_mYResidual->Fit("gaus");

  c_play->SaveAs("./TrackRes.eps");
}
