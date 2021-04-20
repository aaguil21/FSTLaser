#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include "../../src/FstUtil/FstCons.h"

using namespace std;

void plotMcProjection(int sensorId = 0, string mod = "Mod03", string hv = "HV70V")
{
  string mode = "Scan";
  // string mode = "Simple";
  string inputname = Form("../../output/simulation/FstMcProjection_Sensor%d.root",sensorId);
  TFile *File_InPut = TFile::Open(inputname.c_str());
  TH1F *h_mIstProjResX_2Layer    = (TH1F*)File_InPut->Get("h_mIstProjResX_2Layer");
  TH1F *h_mIstProjResY_2Layer    = (TH1F*)File_InPut->Get("h_mIstProjResY_2Layer");
  TH2F *h_mIstProjResXY_2Layer   = (TH2F*)File_InPut->Get("h_mIstProjResXY_2Layer");

  TH1F *h_mFstProjResX_2Layer    = (TH1F*)File_InPut->Get("h_mFstProjResX_2Layer");
  TH1F *h_mFstProjResY_2Layer    = (TH1F*)File_InPut->Get("h_mFstProjResY_2Layer");
  TH2F *h_mFstProjResXY_2Layer   = (TH2F*)File_InPut->Get("h_mFstProjResXY_2Layer");

  TH1F *h_mFstProjResR_2Layer    = (TH1F*)File_InPut->Get("h_mFstProjResR_2Layer");
  TH1F *h_mFstProjResPhi_2Layer  = (TH1F*)File_InPut->Get("h_mFstProjResPhi_2Layer");
  TH2F *h_mFstProjResRPhi_2Layer = (TH2F*)File_InPut->Get("h_mFstProjResRPhi_2Layer");

  TH1F *h_mFstProjResX_2Layer_Rstrips[FST::mFstNumRstripPerSensor];
  TH1F *h_mFstProjResY_2Layer_Rstrips[FST::mFstNumRstripPerSensor];
  TH2F *h_mFstProjResXY_2Layer_Rstrips[FST::mFstNumRstripPerSensor]; // 2D distribution

  TH1F *h_mFstProjResR_2Layer_Rstrips[FST::mFstNumRstripPerSensor];
  TH1F *h_mFstProjResPhi_2Layer_Rstrips[FST::mFstNumRstripPerSensor];
  TH2F *h_mFstProjResRPhi_2Layer_Rstrips[FST::mFstNumRstripPerSensor]; // 2D distribution
  for(int i_rstrip = 0; i_rstrip < FST::mFstNumRstripPerSensor; ++i_rstrip)
  {
    string HistName;
    HistName = Form("h_mFstProjResX_2Layer_Sensor%d_Rstrip%d",sensorId,i_rstrip);
    if(sensorId > 0) HistName = Form("h_mFstProjResX_2Layer_Sensor%d_Rstrip%d",sensorId,i_rstrip+4);
    h_mFstProjResX_2Layer_Rstrips[i_rstrip]  = (TH1F*)File_InPut->Get(HistName.c_str());
    HistName = Form("h_mFstProjResY_2Layer_Sensor%d_Rstrip%d",sensorId,i_rstrip);
    if(sensorId > 0) HistName = Form("h_mFstProjResY_2Layer_Sensor%d_Rstrip%d",sensorId,i_rstrip+4);
    h_mFstProjResY_2Layer_Rstrips[i_rstrip]  = (TH1F*)File_InPut->Get(HistName.c_str());
    HistName = Form("h_mFstProjResXY_2Layer_Sensor%d_Rstrip%d",sensorId,i_rstrip);
    if(sensorId > 0) HistName = Form("h_mFstProjResXY_2Layer_Sensor%d_Rstrip%d",sensorId,i_rstrip+4);
    h_mFstProjResXY_2Layer_Rstrips[i_rstrip] = (TH2F*)File_InPut->Get(HistName.c_str());
    
    HistName = Form("h_mFstProjResR_2Layer_Sensor%d_Rstrip%d",sensorId,i_rstrip);
    if(sensorId > 0) HistName = Form("h_mFstProjResR_2Layer_Sensor%d_Rstrip%d",sensorId,i_rstrip+4);
    h_mFstProjResR_2Layer_Rstrips[i_rstrip]    = (TH1F*)File_InPut->Get(HistName.c_str());
    HistName = Form("h_mFstProjResPhi_2Layer_Sensor%d_Rstrip%d",sensorId,i_rstrip);
    if(sensorId > 0) HistName = Form("h_mFstProjResPhi_2Layer_Sensor%d_Rstrip%d",sensorId,i_rstrip+4);
    h_mFstProjResPhi_2Layer_Rstrips[i_rstrip]  = (TH1F*)File_InPut->Get(HistName.c_str());
    HistName = Form("h_mFstProjResRPhi_2Layer_Sensor%d_Rstrip%d",sensorId,i_rstrip);
    if(sensorId > 0) HistName = Form("h_mFstProjResRPhi_2Layer_Sensor%d_Rstrip%d",sensorId,i_rstrip+4);
    h_mFstProjResRPhi_2Layer_Rstrips[i_rstrip] = (TH2F*)File_InPut->Get(HistName.c_str());
  }

  TH1F *h_mFstSimResR_2Layer     = (TH1F*)File_InPut->Get("h_mFstSimResR_2Layer");
  TH1F *h_mFstSimResPhi_2Layer   = (TH1F*)File_InPut->Get("h_mFstSimResPhi_2Layer");
  TH2F *h_mFstSimResRPhi_2Layer  = (TH2F*)File_InPut->Get("h_mFstSimResRPhi_2Layer");

  string inputData = Form("../../output/configuration/FstTracking_%s_%s_Th4.0Tb2Ped2.5Ped3.5_withPed_withCMNCorr.root",mod.c_str(),hv.c_str());
  TFile *File_InPutData = TFile::Open(inputData.c_str());
  TH1F *h_mClustersTrackFstResR_2Layer;
  TH1F *h_mClustersTrackFstResPhi_2Layer;
  {
    string HistName;
    HistName = Form("h_m%sClustersTrackFstResR_2Layer_Sensor%d",mode.c_str(),sensorId);
    h_mClustersTrackFstResR_2Layer = (TH1F*)File_InPutData->Get(HistName.c_str());
    HistName = Form("h_m%sClustersTrackFstResPhi_2Layer_Sensor%d",mode.c_str(),sensorId);
    h_mClustersTrackFstResPhi_2Layer = (TH1F*)File_InPutData->Get(HistName.c_str());
  }
  TH1F *h_mClustersTrackFstResR_2Layer_Rstrips[FST::mFstNumRstripPerSensor];
  TH1F *h_mClustersTrackFstResPhi_2Layer_Rstrips[FST::mFstNumRstripPerSensor];
  for(int i_rstrp = 0; i_rstrp < FST::mFstNumRstripPerSensor; ++i_rstrp)
  {
    string HistName;
    HistName = Form("h_m%sClustersTrackFstResR_2Layer_Sensor%d_Rstrip%d",mode.c_str(),sensorId,i_rstrp); // R Residual
    if(sensorId > 0) HistName = Form("h_m%sClustersTrackFstResR_2Layer_Sensor%d_Rstrip%d",mode.c_str(),sensorId,i_rstrp+4);
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp] = (TH1F*)File_InPutData->Get(HistName.c_str());
    HistName = Form("h_m%sClustersTrackFstResPhi_2Layer_Sensor%d_Rstrip%d",mode.c_str(),sensorId,i_rstrp); // R Residual
    if(sensorId > 0) HistName = Form("h_m%sClustersTrackFstResPhi_2Layer_Sensor%d_Rstrip%d",mode.c_str(),sensorId,i_rstrp+4);
    h_mClustersTrackFstResPhi_2Layer_Rstrips[i_rstrp] = (TH1F*)File_InPutData->Get(HistName.c_str());
  }

  // calculate scaling factor
  double NormR_InteMC = h_mFstProjResR_2Layer->GetMaximum();
  double NormR_InteData = h_mClustersTrackFstResR_2Layer->GetMaximum();
  double NormR_InteIST = h_mFstSimResR_2Layer->GetMaximum();
  double NormR_Data[4];
  double NormR_MC[4];
  for(int i_rstrp = 0; i_rstrp < 4; ++i_rstrp)
  {
    NormR_Data[i_rstrp] = h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->GetMaximum();
    NormR_MC[i_rstrp] = h_mFstProjResR_2Layer_Rstrips[i_rstrp]->GetMaximum();
  }
  double NormPhi_InteMC = h_mFstProjResPhi_2Layer->GetMaximum();
  double NormPhi_InteData = h_mClustersTrackFstResPhi_2Layer->GetMaximum();
  double NormPhi_InteIST = h_mFstSimResPhi_2Layer->GetMaximum();
  double NormPhi_Data[4];
  double NormPhi_MC[4];
  for(int i_rstrp = 0; i_rstrp < 4; ++i_rstrp)
  {
    NormPhi_Data[i_rstrp] = h_mClustersTrackFstResPhi_2Layer_Rstrips[i_rstrp]->GetMaximum();
    NormPhi_MC[i_rstrp] = h_mFstProjResPhi_2Layer_Rstrips[i_rstrp]->GetMaximum();
  }

  TCanvas *c_Residual = new TCanvas("c_Residual","c_Residual",10,10,900,900);
  c_Residual->Divide(3,3);
  for(int i_pad = 0; i_pad < 9; ++ i_pad)
  {
    c_Residual->cd(i_pad+1);
    c_Residual->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Residual->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Residual->cd(i_pad+1)->SetTicks(1,1);
    c_Residual->cd(i_pad+1)->SetGrid(0,0);
  }

  // IST
  c_Residual->cd(1);
  h_mIstProjResXY_2Layer->SetTitle("IST");
  h_mIstProjResXY_2Layer->GetXaxis()->SetTitle("x^{IST}_{proj} - x^{IST}_{gen} (mm)");
  h_mIstProjResXY_2Layer->GetXaxis()->SetTitleSize(0.06);
  h_mIstProjResXY_2Layer->GetXaxis()->CenterTitle();
  h_mIstProjResXY_2Layer->GetYaxis()->SetTitle("y^{IST}_{proj} - y^{IST}_{gen} (mm)");
  h_mIstProjResXY_2Layer->GetYaxis()->SetTitleSize(0.06);
  h_mIstProjResXY_2Layer->GetYaxis()->CenterTitle();
  h_mIstProjResXY_2Layer->Draw("colz");

  c_Residual->cd(2);
  h_mIstProjResX_2Layer->SetTitle("IST");
  h_mIstProjResX_2Layer->GetXaxis()->SetTitle("x^{IST}_{proj} - x^{IST}_{gen} (mm)");
  h_mIstProjResX_2Layer->GetXaxis()->SetTitleSize(0.06);
  h_mIstProjResX_2Layer->GetXaxis()->CenterTitle();
  h_mIstProjResX_2Layer->Draw("hE");

  c_Residual->cd(3);
  h_mIstProjResY_2Layer->SetTitle("IST");
  h_mIstProjResY_2Layer->GetXaxis()->SetTitle("y^{IST}_{proj} - y^{IST}_{gen} (mm)");
  h_mIstProjResY_2Layer->GetXaxis()->CenterTitle();
  h_mIstProjResY_2Layer->GetXaxis()->SetTitleSize(0.06);
  h_mIstProjResY_2Layer->Draw("hE");

  //FST X-Y
  c_Residual->cd(4);
  h_mFstProjResXY_2Layer->SetTitle("FST");
  h_mFstProjResXY_2Layer->GetXaxis()->SetTitle("x_{ro} - x_{proj} (mm)");
  h_mFstProjResXY_2Layer->GetXaxis()->SetTitleSize(0.06);
  h_mFstProjResXY_2Layer->GetXaxis()->CenterTitle();
  h_mFstProjResXY_2Layer->GetYaxis()->SetTitle("y_{ro} - y_{proj} (mm)");
  h_mFstProjResXY_2Layer->GetYaxis()->SetTitleSize(0.06);
  h_mFstProjResXY_2Layer->GetYaxis()->CenterTitle();
  h_mFstProjResXY_2Layer->Draw("colz");

  c_Residual->cd(5);
  h_mFstProjResX_2Layer->SetTitle("FST");
  h_mFstProjResX_2Layer->GetXaxis()->SetTitle("x_{ro} - x_{proj} (mm)");
  h_mFstProjResX_2Layer->GetXaxis()->SetTitleSize(0.06);
  h_mFstProjResX_2Layer->GetXaxis()->CenterTitle();
  h_mFstProjResX_2Layer->SetLineStyle(1);
  h_mFstProjResX_2Layer->SetLineWidth(2);
  h_mFstProjResX_2Layer->SetLineColor(1);
  h_mFstProjResX_2Layer->Draw("hE");

  c_Residual->cd(6);
  h_mFstProjResY_2Layer->SetTitle("FST");
  h_mFstProjResY_2Layer->GetXaxis()->SetTitle("y_{ro} - y_{proj} (mm)");
  h_mFstProjResY_2Layer->GetXaxis()->SetTitleSize(0.06);
  h_mFstProjResY_2Layer->GetXaxis()->CenterTitle();
  h_mFstProjResY_2Layer->SetLineStyle(1);
  h_mFstProjResY_2Layer->SetLineWidth(2);
  h_mFstProjResY_2Layer->SetLineColor(1);
  h_mFstProjResY_2Layer->Draw("hE");

  //FST R-Phi
  c_Residual->cd(7);
  h_mFstProjResRPhi_2Layer->SetTitle("FST");
  h_mFstProjResRPhi_2Layer->GetXaxis()->SetTitle("r_{ro} - r_{proj} (mm)");
  h_mFstProjResRPhi_2Layer->GetXaxis()->SetTitleSize(0.06);
  h_mFstProjResRPhi_2Layer->GetXaxis()->CenterTitle();
  h_mFstProjResRPhi_2Layer->GetYaxis()->SetTitle("#phi_{ro} - #phi_{proj} (rad)");
  h_mFstProjResRPhi_2Layer->GetYaxis()->SetTitleSize(0.06);
  h_mFstProjResRPhi_2Layer->GetYaxis()->CenterTitle();
  h_mFstProjResRPhi_2Layer->Draw("colz");

  c_Residual->cd(8);
  h_mFstProjResR_2Layer->SetTitle("FST");
  h_mFstProjResR_2Layer->GetXaxis()->SetTitle("r_{ro} - r_{proj} (mm)");
  h_mFstProjResR_2Layer->GetXaxis()->SetTitleSize(0.06);
  h_mFstProjResR_2Layer->GetXaxis()->CenterTitle();
  h_mFstProjResR_2Layer->GetXaxis()->SetRangeUser(-100.0,100.0);
  h_mFstProjResR_2Layer->SetLineStyle(1);
  h_mFstProjResR_2Layer->SetLineWidth(2);
  h_mFstProjResR_2Layer->SetLineColor(1);
  h_mFstProjResR_2Layer->Scale(NormR_InteData/NormR_InteMC);
  h_mFstProjResR_2Layer->Draw("h");
  h_mFstSimResR_2Layer->SetLineStyle(2);
  h_mFstSimResR_2Layer->SetLineWidth(2);
  h_mFstSimResR_2Layer->SetLineColor(2);
  h_mFstSimResR_2Layer->Scale(NormR_InteData/NormR_InteIST);
  h_mFstSimResR_2Layer->Draw("h same");
  h_mClustersTrackFstResR_2Layer->SetMarkerStyle(24);
  h_mClustersTrackFstResR_2Layer->SetMarkerSize(0.8);
  h_mClustersTrackFstResR_2Layer->SetMarkerColor(1);
  h_mClustersTrackFstResR_2Layer->SetLineColor(1);
  h_mClustersTrackFstResR_2Layer->Draw("pE same");
  TLegend *leg = new TLegend(0.60,0.5,0.85,0.7);
  leg->AddEntry(h_mClustersTrackFstResR_2Layer,"Data","p");
  leg->AddEntry(h_mFstProjResR_2Layer,"MC","l");
  leg->AddEntry(h_mFstSimResR_2Layer,"IST Proj Err.","l");
  leg->Draw("same");

  c_Residual->cd(9);
  h_mFstProjResPhi_2Layer->SetTitle("FST");
  h_mFstProjResPhi_2Layer->GetXaxis()->SetTitle("#phi_{ro} - #phi_{proj} (rad)");
  h_mFstProjResPhi_2Layer->GetXaxis()->SetTitleSize(0.06);
  h_mFstProjResPhi_2Layer->GetXaxis()->CenterTitle();
  h_mFstProjResPhi_2Layer->GetXaxis()->SetRangeUser(-0.05,0.05);
  h_mFstProjResPhi_2Layer->SetLineStyle(1);
  h_mFstProjResPhi_2Layer->SetLineWidth(2);
  h_mFstProjResPhi_2Layer->SetLineColor(1);
  h_mFstProjResPhi_2Layer->Scale(NormPhi_InteData/NormPhi_InteMC);
  h_mFstProjResPhi_2Layer->Draw("h");
  h_mFstSimResPhi_2Layer->SetLineStyle(2);
  h_mFstSimResPhi_2Layer->SetLineWidth(2);
  h_mFstSimResPhi_2Layer->SetLineColor(2);
  h_mFstSimResPhi_2Layer->Scale(NormPhi_InteData/NormPhi_InteIST);
  h_mFstSimResPhi_2Layer->Draw("h same");
  h_mClustersTrackFstResPhi_2Layer->SetMarkerStyle(24);
  h_mClustersTrackFstResPhi_2Layer->SetMarkerSize(0.8);
  h_mClustersTrackFstResPhi_2Layer->SetMarkerColor(1);
  h_mClustersTrackFstResPhi_2Layer->SetLineColor(1);
  h_mClustersTrackFstResPhi_2Layer->Draw("pE same");
  leg->Draw("same");

  string FigName = Form("./figures/c_Residual_%s_%s_Sensor%d.eps",mode.c_str(),hv.c_str(),sensorId);
  c_Residual->SaveAs(FigName.c_str());

  TCanvas *c_ResDiff = new TCanvas("c_ResDiff","c_ResDiff",10,10,900,1200);
  c_ResDiff->Divide(3,4);
  for(int i_pad = 0; i_pad < 12; ++ i_pad)
  {
    c_ResDiff->cd(i_pad+1);
    c_ResDiff->cd(i_pad+1)->SetLeftMargin(0.15);
    c_ResDiff->cd(i_pad+1)->SetBottomMargin(0.15);
    c_ResDiff->cd(i_pad+1)->SetTicks(1,1);
    c_ResDiff->cd(i_pad+1)->SetGrid(0,0);
  }

  for(int i_rstrp = 0; i_rstrp < 4; ++i_rstrp)
  {
    c_ResDiff->cd(i_rstrp*3+1);
    h_mFstProjResRPhi_2Layer_Rstrips[i_rstrp]->GetXaxis()->SetTitle("r_{ro} - r_{proj} (mm)");
    h_mFstProjResRPhi_2Layer_Rstrips[i_rstrp]->GetXaxis()->SetTitleSize(0.06);
    h_mFstProjResRPhi_2Layer_Rstrips[i_rstrp]->GetXaxis()->CenterTitle();
    h_mFstProjResRPhi_2Layer_Rstrips[i_rstrp]->GetYaxis()->SetTitle("#phi_{ro} - #phi_{proj} (rad)");
    h_mFstProjResRPhi_2Layer_Rstrips[i_rstrp]->GetYaxis()->SetTitleSize(0.06);
    h_mFstProjResRPhi_2Layer_Rstrips[i_rstrp]->GetYaxis()->CenterTitle();
    h_mFstProjResRPhi_2Layer_Rstrips[i_rstrp]->Draw("colz");

    c_ResDiff->cd(i_rstrp*3+2);
    h_mFstProjResR_2Layer_Rstrips[i_rstrp]->GetXaxis()->SetTitle("r_{ro} - r_{proj} (mm)");
    h_mFstProjResR_2Layer_Rstrips[i_rstrp]->GetXaxis()->SetTitleSize(0.06);
    h_mFstProjResR_2Layer_Rstrips[i_rstrp]->GetXaxis()->CenterTitle();
    h_mFstProjResR_2Layer_Rstrips[i_rstrp]->GetXaxis()->SetRangeUser(-100.0,100.0);
    if(NormR_Data[i_rstrp] > 0 && NormR_MC[i_rstrp] > 0) h_mFstProjResR_2Layer_Rstrips[i_rstrp]->Scale(NormR_Data[i_rstrp]/NormR_MC[i_rstrp]);
    h_mFstProjResR_2Layer_Rstrips[i_rstrp]->GetYaxis()->SetRangeUser(0.0,h_mFstProjResR_2Layer_Rstrips[i_rstrp]->GetMaximum()*1.1);
    h_mFstProjResR_2Layer_Rstrips[i_rstrp]->SetLineStyle(1);
    h_mFstProjResR_2Layer_Rstrips[i_rstrp]->SetLineWidth(2);
    h_mFstProjResR_2Layer_Rstrips[i_rstrp]->SetLineColor(i_rstrp+1);
    h_mFstProjResR_2Layer_Rstrips[i_rstrp]->Draw("h");
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->SetMarkerStyle(24);
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->SetMarkerSize(0.8);
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->SetMarkerColor(i_rstrp+1);
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->SetLineColor(i_rstrp+1);
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->Draw("pE same");
    TLegend *leg = new TLegend(0.60,0.4,0.85,0.7);
    string legend = Form("Rstrip %d",i_rstrp); 
    leg->AddEntry((TObject*)0,legend.c_str(),"");
    legend = Form("%s Cluster",mode.c_str());
    leg->AddEntry(h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp],legend.c_str(),"p");
    leg->AddEntry(h_mFstProjResR_2Layer_Rstrips[i_rstrp],"MC","l");
    leg->Draw("same");

    c_ResDiff->cd(i_rstrp*3+3);
    h_mFstProjResPhi_2Layer_Rstrips[i_rstrp]->GetXaxis()->SetTitle("#phi_{ro} - #phi_{proj} (rad)");
    h_mFstProjResPhi_2Layer_Rstrips[i_rstrp]->GetXaxis()->SetTitleSize(0.06);
    h_mFstProjResPhi_2Layer_Rstrips[i_rstrp]->GetXaxis()->CenterTitle();
    h_mFstProjResPhi_2Layer_Rstrips[i_rstrp]->GetXaxis()->SetRangeUser(-0.05,0.05);
    if(NormPhi_Data[i_rstrp] > 0 && NormPhi_MC[i_rstrp] > 0) h_mFstProjResPhi_2Layer_Rstrips[i_rstrp]->Scale(NormPhi_Data[i_rstrp]/NormPhi_MC[i_rstrp]);
    h_mFstProjResPhi_2Layer_Rstrips[i_rstrp]->GetYaxis()->SetRangeUser(0.0,h_mFstProjResPhi_2Layer_Rstrips[i_rstrp]->GetMaximum()*1.1);
    h_mFstProjResPhi_2Layer_Rstrips[i_rstrp]->SetLineStyle(1);
    h_mFstProjResPhi_2Layer_Rstrips[i_rstrp]->SetLineWidth(2);
    h_mFstProjResPhi_2Layer_Rstrips[i_rstrp]->SetLineColor(i_rstrp+1);
    h_mFstProjResPhi_2Layer_Rstrips[i_rstrp]->Draw("h");
    h_mClustersTrackFstResPhi_2Layer_Rstrips[i_rstrp]->SetMarkerStyle(24);
    h_mClustersTrackFstResPhi_2Layer_Rstrips[i_rstrp]->SetMarkerSize(0.8);
    h_mClustersTrackFstResPhi_2Layer_Rstrips[i_rstrp]->SetMarkerColor(i_rstrp+1);
    h_mClustersTrackFstResPhi_2Layer_Rstrips[i_rstrp]->SetLineColor(i_rstrp+1);
    h_mClustersTrackFstResPhi_2Layer_Rstrips[i_rstrp]->Draw("pE same");
    leg->Draw("same");
  }

  FigName = Form("./figures/c_ResDiff_%s_%s_Sensor%d.eps",mode.c_str(),hv.c_str(),sensorId);
  c_ResDiff->SaveAs(FigName.c_str());
}
