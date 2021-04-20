#include <iostream>
#include <string>

#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include <TF1.h>
#include <TProfile2D.h>
#include <TGraph.h>
#include <TStyle.h>
#include "./draw.h"
#include "../../src/FstUtil/FstCons.h"

using namespace std;

void plotEfficiencyHV(string mod = "Mod01")
{
  bool isSavePed          = true;
  bool isApplyCMNCorr     = true;
  float nFstHitsCut       = 4.0;
  int numOfUsedTimeBins   = 2;
  float nFstThresholdCut2 = 2.5;
  float nFstThresholdCut1 = 3.5;

  std::string pedMode = "withPed";
  if(!isSavePed) pedMode = "woPed";
  std::string cmnMode = "withCMNCorr";
  if(!isApplyCMNCorr) cmnMode = "woCMNCorr";

  const int numOfHV = 4;
  string hv[numOfHV] = {"HV70V","HV100V","HV120V","HV140V"};
  double highVolt[numOfHV] = {70.0,100.0,120.0,140.0};
  TFile *File_InPut[numOfHV];
  TH1F *h_mSensorScanClustersTrackFstEffR_2Layer[numOfHV][FST::mFstNumSensorsPerModule][FST::mFstNumMatching];
  TH1F *h_mSensorScanClustersTrackFstEffPhi_2Layer[numOfHV][FST::mFstNumSensorsPerModule][FST::mFstNumMatching];

  for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
  {
    string inputfile = Form("../../output/signal/%s/FstEfficiency_%s_%s_Th%1.1fTb%dPed%1.1fPed%1.1f_%s_%s.root",mod.c_str(),mod.c_str(),hv[i_hv].c_str(),nFstHitsCut,numOfUsedTimeBins,nFstThresholdCut2,nFstThresholdCut1,pedMode.c_str(),cmnMode.c_str());
    File_InPut[i_hv] = TFile::Open(inputfile.c_str());
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
      {
	string HistName;
	HistName = Form("h_mScanClustersTrackFstEffR_2Layer_Sensor%d_SF%d",i_sensor,i_match);
	h_mSensorScanClustersTrackFstEffR_2Layer[i_hv][i_sensor][i_match] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());

	HistName = Form("h_mScanClustersTrackFstEffPhi_2Layer_Sensor%d_SF%d",i_sensor,i_match);
	h_mSensorScanClustersTrackFstEffPhi_2Layer[i_hv][i_sensor][i_match] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());
      }
    }
  }

  TH1F *h_play = new TH1F("h_play","h_play",400,-9.5,390.5);
  for(int i_bin = 0; i_bin < 400; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("R (mm)");
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetNdivisions(505);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetTitle("Efficiency");
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetRangeUser(0.0,1.1);

  const int markerColor[4] = {2, 4, 6,kGray+2};
  const int markerStyle[4] = {24, 24, 24, 20};
  TCanvas *c_EfficiencyR = new TCanvas("c_EfficiencyR","c_EfficiencyR",10,10,1200,400);
  c_EfficiencyR->Divide(3,1);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_EfficiencyR->cd(i_pad+1)->SetLeftMargin(0.15);
    c_EfficiencyR->cd(i_pad+1)->SetBottomMargin(0.15);
    c_EfficiencyR->cd(i_pad+1)->SetTicks(1,1);
    c_EfficiencyR->cd(i_pad+1)->SetGrid(0,0);
  }

  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    c_EfficiencyR->cd(i_sensor+1);
    string title = Form("Sensor%d",i_sensor);
    h_play->SetTitle(title.c_str());
    h_play->GetXaxis()->SetRangeUser(FST::rMin[i_sensor],FST::rMax[i_sensor]);
    h_play->DrawCopy("hE");
    PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMin[i_sensor], 0.0, 1.05, 1, 2, 2);
    PlotLine(FST::mFstRMax[i_sensor], FST::mFstRMax[i_sensor], 0.0, 1.05, 1, 2, 2);
    PlotLine(FST::rMin[i_sensor], FST::rMax[i_sensor], 0.90, 0.90, 4, 1, 2);
    PlotLine(FST::rMin[i_sensor], FST::rMax[i_sensor], 0.95, 0.95, 4, 1, 2);
    for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
    {
      h_mSensorScanClustersTrackFstEffR_2Layer[i_hv][i_sensor][0]->SetMarkerStyle(markerStyle[i_hv]);
      h_mSensorScanClustersTrackFstEffR_2Layer[i_hv][i_sensor][0]->SetMarkerColor(markerColor[i_hv]);
      h_mSensorScanClustersTrackFstEffR_2Layer[i_hv][i_sensor][0]->SetMarkerSize(0.8);
      h_mSensorScanClustersTrackFstEffR_2Layer[i_hv][i_sensor][0]->SetLineColor(markerColor[i_hv]);
      h_mSensorScanClustersTrackFstEffR_2Layer[i_hv][i_sensor][0]->Draw("PE same");
    }
    TLegend *leg = new TLegend(0.45,0.25,0.75,0.55);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetMargin(0.1);
    // leg->SetNColumns(2);
    leg->AddEntry(h_mSensorScanClustersTrackFstEffR_2Layer[0][i_sensor][0],hv[0].c_str(),"P");
    leg->AddEntry(h_mSensorScanClustersTrackFstEffR_2Layer[1][i_sensor][0],hv[1].c_str(),"P");
    leg->AddEntry(h_mSensorScanClustersTrackFstEffR_2Layer[2][i_sensor][0],hv[2].c_str(),"P");
    leg->AddEntry(h_mSensorScanClustersTrackFstEffR_2Layer[3][i_sensor][0],hv[3].c_str(),"P");
    leg->Draw("same");
  }
  c_EfficiencyR->SaveAs("./figures/c_EfficiencyRHVScan.eps");

  TH1F *h_playPhi = new TH1F("h_playPhi","h_playPhi",100,-0.5,0.5);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_playPhi->SetBinContent(i_bin+1,-10.0);
    h_playPhi->SetBinError(i_bin+1,1.0);
  }
  h_playPhi->SetStats(0);
  h_playPhi->GetXaxis()->SetTitle("#phi (rad)");
  h_playPhi->GetXaxis()->SetTitleSize(0.06);
  h_playPhi->GetXaxis()->SetNdivisions(505);
  h_playPhi->GetXaxis()->SetLabelSize(0.04);
  h_playPhi->GetYaxis()->SetTitle("Efficiency");
  h_playPhi->GetYaxis()->SetNdivisions(505);
  h_playPhi->GetYaxis()->SetTitleSize(0.06);
  h_playPhi->GetYaxis()->SetLabelSize(0.04);
  h_playPhi->GetYaxis()->SetRangeUser(0.0,1.1);

  TCanvas *c_EfficiencyPhi = new TCanvas("c_EfficiencyPhi","c_EfficiencyPhi",10,10,1200,400);
  c_EfficiencyPhi->Divide(3,1);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_EfficiencyPhi->cd(i_pad+1)->SetLeftMargin(0.15);
    c_EfficiencyPhi->cd(i_pad+1)->SetBottomMargin(0.15);
    c_EfficiencyPhi->cd(i_pad+1)->SetTicks(1,1);
    c_EfficiencyPhi->cd(i_pad+1)->SetGrid(0,0);
  }

  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    c_EfficiencyPhi->cd(i_sensor+1);
    string title = Form("Sensor%d",i_sensor);
    h_playPhi->SetTitle(title.c_str());
    h_playPhi->GetXaxis()->SetRangeUser(FST::phiMin[i_sensor],FST::phiMax[i_sensor]);
    h_playPhi->DrawCopy("hE");
    PlotLine(FST::mFstPhiMin[i_sensor], FST::mFstPhiMin[i_sensor], 0.0, 1.05, 1, 2, 2);
    PlotLine(FST::mFstPhiMax[i_sensor], FST::mFstPhiMax[i_sensor], 0.0, 1.05, 1, 2, 2);
    PlotLine(FST::phiMin[i_sensor], FST::phiMax[i_sensor], 0.90, 0.90, 4, 1, 2);
    PlotLine(FST::phiMin[i_sensor], FST::phiMax[i_sensor], 0.95, 0.95, 4, 1, 2);
    for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
    {
      h_mSensorScanClustersTrackFstEffPhi_2Layer[i_hv][i_sensor][0]->SetMarkerStyle(markerStyle[i_hv]);
      h_mSensorScanClustersTrackFstEffPhi_2Layer[i_hv][i_sensor][0]->SetMarkerColor(markerColor[i_hv]);
      h_mSensorScanClustersTrackFstEffPhi_2Layer[i_hv][i_sensor][0]->SetMarkerSize(0.8);
      h_mSensorScanClustersTrackFstEffPhi_2Layer[i_hv][i_sensor][0]->SetLineColor(markerColor[i_hv]);
      h_mSensorScanClustersTrackFstEffPhi_2Layer[i_hv][i_sensor][0]->Draw("PE same");
    }
    TLegend *leg = new TLegend(0.45,0.25,0.75,0.55);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetMargin(0.1);
    // leg->SetNColumns(2);
    leg->AddEntry(h_mSensorScanClustersTrackFstEffPhi_2Layer[0][i_sensor][0],hv[0].c_str(),"P");
    leg->AddEntry(h_mSensorScanClustersTrackFstEffPhi_2Layer[1][i_sensor][0],hv[1].c_str(),"P");
    leg->AddEntry(h_mSensorScanClustersTrackFstEffPhi_2Layer[2][i_sensor][0],hv[2].c_str(),"P");
    leg->AddEntry(h_mSensorScanClustersTrackFstEffPhi_2Layer[3][i_sensor][0],hv[3].c_str(),"P");
    leg->Draw("same");
  }
  c_EfficiencyPhi->SaveAs("./figures/c_EfficiencyPhiHVScan.eps");
}
