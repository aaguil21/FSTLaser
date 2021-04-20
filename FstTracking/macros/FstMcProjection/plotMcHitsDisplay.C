#include <TCanvas.h>
#include <TH1F.h>

void plotMcHitsDisplay(int sensorId = 0)
{
  string inputname = Form("../../output/simulation/FstMcProjection_Sensor%d.root",sensorId);
  TFile *File_InPut = TFile::Open(inputname.c_str());
  TH1F *h_mIst1Display = (TH1F*)File_InPut->Get("h_mIst1Display");
  TH1F *h_mIst3Display = (TH1F*)File_InPut->Get("h_mIst3Display");
  TH1F *h_mFstDisplay  = (TH1F*)File_InPut->Get("h_mFstDisplay");

  TCanvas *c_HitsDisplay = new TCanvas("c_HitsDisplay","c_HitsDisplay",10,10,1200,400);
  c_HitsDisplay->Divide(3,1);
  for(int i_pad = 0; i_pad < 3; ++ i_pad)
  {
    c_HitsDisplay->cd(i_pad+1);
    c_HitsDisplay->cd(i_pad+1)->SetLeftMargin(0.15);
    c_HitsDisplay->cd(i_pad+1)->SetBottomMargin(0.15);
    c_HitsDisplay->cd(i_pad+1)->SetTicks(1,1);
    c_HitsDisplay->cd(i_pad+1)->SetGrid(0,0);
  }

  c_HitsDisplay->cd(1);
  h_mIst1Display->GetXaxis()->SetTitle("x (mm)");
  h_mIst1Display->GetXaxis()->SetTitleSize(0.06);
  h_mIst1Display->GetXaxis()->CenterTitle();
  h_mIst1Display->GetYaxis()->SetTitle("y (mm)");
  h_mIst1Display->GetYaxis()->SetTitleSize(0.06);
  h_mIst1Display->GetYaxis()->CenterTitle();
  h_mIst1Display->Draw("colz");

  c_HitsDisplay->cd(2);
  h_mIst3Display->GetXaxis()->SetTitle("x (mm)");
  h_mIst3Display->GetXaxis()->SetTitleSize(0.06);
  h_mIst3Display->GetXaxis()->CenterTitle();
  h_mIst3Display->GetYaxis()->SetTitle("y (mm)");
  h_mIst3Display->GetYaxis()->SetTitleSize(0.06);
  h_mIst3Display->GetYaxis()->CenterTitle();
  h_mIst3Display->Draw("colz");

  c_HitsDisplay->cd(3);
  h_mFstDisplay->GetXaxis()->SetTitle("r (mm)");
  h_mFstDisplay->GetXaxis()->SetTitleSize(0.06);
  h_mFstDisplay->GetXaxis()->CenterTitle();
  h_mFstDisplay->GetYaxis()->SetTitle("#phi (rad)");
  h_mFstDisplay->GetYaxis()->SetTitleSize(0.06);
  h_mFstDisplay->GetYaxis()->CenterTitle();
  h_mFstDisplay->Draw("colz");

  string FigName = Form("./figures/c_HitsDisplay_Sensor%d.eps",sensorId);
  c_HitsDisplay->SaveAs(FigName.c_str());
}
