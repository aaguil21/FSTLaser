#include <TCanvas.h>
#include <TH1F.h>

void plotCosmicAngle(int sensorId = 0)
{
  string inputname = Form("../../output/simulation/FstMcProjection_Sensor%d.root",sensorId);
  TFile *File_InPut = TFile::Open(inputname.c_str());
  TH1F *h_mPolarAngle = (TH1F*)File_InPut->Get("h_mPolarAngle");
  TH1F *h_mAzimuthalAngle = (TH1F*)File_InPut->Get("h_mAzimuthalAngle");

  TCanvas *c_AngleDisplay = new TCanvas("c_AngleDisplay","c_AngleDisplay",10,10,400,800);
  c_AngleDisplay->Divide(1,2);
  for(int i_pad = 0; i_pad < 2; ++ i_pad)
  {
    c_AngleDisplay->cd(i_pad+1);
    c_AngleDisplay->cd(i_pad+1)->SetLeftMargin(0.15);
    c_AngleDisplay->cd(i_pad+1)->SetBottomMargin(0.15);
    c_AngleDisplay->cd(i_pad+1)->SetTicks(1,1);
    c_AngleDisplay->cd(i_pad+1)->SetGrid(0,0);
  }

  c_AngleDisplay->cd(1);
  h_mPolarAngle->GetXaxis()->SetTitle("#theta (rad)");
  h_mPolarAngle->GetXaxis()->SetTitleSize(0.06);
  h_mPolarAngle->GetXaxis()->CenterTitle();
  h_mPolarAngle->GetYaxis()->SetTitle("Counts");
  h_mPolarAngle->GetYaxis()->SetTitleSize(0.06);
  h_mPolarAngle->GetYaxis()->CenterTitle();
  h_mPolarAngle->GetYaxis()->SetRangeUser(0.0,h_mPolarAngle->GetMaximum()*1.2);
  h_mPolarAngle->Draw("h");

  c_AngleDisplay->cd(2);
  h_mAzimuthalAngle->GetXaxis()->SetTitle("#phi (rad)");
  h_mAzimuthalAngle->GetXaxis()->SetTitleSize(0.06);
  h_mAzimuthalAngle->GetXaxis()->CenterTitle();
  h_mAzimuthalAngle->GetYaxis()->SetTitle("Counts");
  h_mAzimuthalAngle->GetYaxis()->SetTitleSize(0.06);
  h_mAzimuthalAngle->GetYaxis()->CenterTitle();
  h_mAzimuthalAngle->GetYaxis()->SetRangeUser(0.0,h_mAzimuthalAngle->GetMaximum()*1.2);
  h_mAzimuthalAngle->Draw("h");

  c_AngleDisplay->SaveAs("./figures/c_AngleDisplay.eps");
}
