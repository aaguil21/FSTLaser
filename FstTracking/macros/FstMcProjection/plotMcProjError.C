#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>

using namespace std;

void plotMcProjError(int sensorId = 0)
{
  string mode = "Scan";
  // string mode = "Simple";
  string inputname = Form("../../output/simulation/FstMcProjection_Sensor%d.root",sensorId);
  TFile *File_InPut = TFile::Open(inputname.c_str());
  TH1F *h_mFstSimResR_2Layer     = (TH1F*)File_InPut->Get("h_mFstSimResR_2Layer");
  TH1F *h_mFstSimResPhi_2Layer   = (TH1F*)File_InPut->Get("h_mFstSimResPhi_2Layer");
  TH2F *h_mFstSimResRPhi_2Layer  = (TH2F*)File_InPut->Get("h_mFstSimResRPhi_2Layer");

  TCanvas *c_ProjError = new TCanvas("c_ProjError","c_ProjError",10,10,1200,400);
  c_ProjError->Divide(3,1);
  for(int i_pad = 0; i_pad < 3; ++ i_pad)
  {
    c_ProjError->cd(i_pad+1);
    c_ProjError->cd(i_pad+1)->SetLeftMargin(0.15);
    c_ProjError->cd(i_pad+1)->SetBottomMargin(0.15);
    c_ProjError->cd(i_pad+1)->SetTicks(1,1);
    c_ProjError->cd(i_pad+1)->SetGrid(0,0);
  }

  c_ProjError->cd(1);
  h_mFstSimResRPhi_2Layer->SetTitle("IST");
  h_mFstSimResRPhi_2Layer->GetXaxis()->SetTitle("r^{IST}_{proj} - r^{IST}_{gen} (mm)");
  h_mFstSimResRPhi_2Layer->GetXaxis()->SetTitleSize(0.06);
  h_mFstSimResRPhi_2Layer->GetXaxis()->CenterTitle();
  h_mFstSimResRPhi_2Layer->GetYaxis()->SetTitle("#phi^{IST}_{proj} - #phi^{IST}_{gen} (rad)");
  h_mFstSimResRPhi_2Layer->GetYaxis()->SetTitleSize(0.06);
  h_mFstSimResRPhi_2Layer->GetYaxis()->CenterTitle();
  h_mFstSimResRPhi_2Layer->Draw("colz");

  c_ProjError->cd(2);
  h_mFstSimResR_2Layer->SetTitle("IST");
  h_mFstSimResR_2Layer->GetXaxis()->SetTitle("r^{IST}_{proj} - r^{IST}_{gen} (mm)");
  h_mFstSimResR_2Layer->GetXaxis()->SetTitleSize(0.06);
  h_mFstSimResR_2Layer->GetXaxis()->CenterTitle();
  h_mFstSimResR_2Layer->SetLineStyle(1);
  h_mFstSimResR_2Layer->SetLineWidth(2);
  h_mFstSimResR_2Layer->SetLineColor(1);
  h_mFstSimResR_2Layer->Draw("hE");

  c_ProjError->cd(3);
  h_mFstSimResPhi_2Layer->SetTitle("IST");
  h_mFstSimResPhi_2Layer->GetXaxis()->SetTitle("#phi^{IST}_{proj} - #phi^{IST}_{gen} (rad)");
  h_mFstSimResPhi_2Layer->GetXaxis()->SetTitleSize(0.06);
  h_mFstSimResPhi_2Layer->GetXaxis()->CenterTitle();
  h_mFstSimResPhi_2Layer->SetLineStyle(1);
  h_mFstSimResPhi_2Layer->SetLineWidth(2);
  h_mFstSimResPhi_2Layer->SetLineColor(1);
  h_mFstSimResPhi_2Layer->Draw("hE");

  string FigName = Form("./figures/c_ProjError_Sensor%d.eps",sensorId);
  c_ProjError->SaveAs(FigName.c_str());
}

