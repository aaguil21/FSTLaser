#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>

#include <string>

void plotAngleCorr()
{
  gStyle->SetOptStat(1101);
  gStyle->SetOptFit(1001);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.2); gStyle->SetStatH(0.3);

  string inputfile = "./IstTracking.root";
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mAdc_Layer1 = (TH1F*)File_InPut->Get("h_mAdc_Layer1");
  TH1F *h_mAdcAngleCorr_Layer1 = (TH1F*)File_InPut->Get("h_mAdcAngleCorr_Layer1");

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,500,1000);
  c_play->Divide(1,2);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_play->cd(i_pad+1)->SetLeftMargin(0.15);
    c_play->cd(i_pad+1)->SetBottomMargin(0.15);
    c_play->cd(i_pad+1)->SetTicks(1,1);
    c_play->cd(i_pad+1)->SetGrid(0,0);
  }

  c_play->cd(1);
  h_mAdc_Layer1->SetTitle("ADC spectrum");
  h_mAdc_Layer1->GetXaxis()->SetTitle("adc");
  h_mAdc_Layer1->GetXaxis()->SetTitleSize(0.06);
  h_mAdc_Layer1->GetYaxis()->SetTitle("No. Hits");
  h_mAdc_Layer1->GetYaxis()->SetTitleSize(0.06);
  h_mAdc_Layer1->Draw();
  h_mAdc_Layer1->Fit("landau");

  c_play->cd(2);
  h_mAdcAngleCorr_Layer1->SetTitle("angle corrected ADC spectrum");
  h_mAdcAngleCorr_Layer1->GetXaxis()->SetTitle("adc");
  h_mAdcAngleCorr_Layer1->GetXaxis()->SetTitleSize(0.06);
  h_mAdcAngleCorr_Layer1->GetYaxis()->SetTitle("No. Hits");
  h_mAdcAngleCorr_Layer1->GetYaxis()->SetTitleSize(0.06);
  h_mAdcAngleCorr_Layer1->Draw();
  h_mAdcAngleCorr_Layer1->Fit("landau");

  c_play->SaveAs("./ADCspectrum.eps");
}
