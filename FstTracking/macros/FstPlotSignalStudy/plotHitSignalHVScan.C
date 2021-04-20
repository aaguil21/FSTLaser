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

void plotHitSignalHVScan(string mod = "Mod01")
{
  const int numOfHV = 11;
  string hv[numOfHV] = {"HV40V","HV50V","HV60V","HV70V","HV80V","HV100V","HV120V","HV140V","HV160V","HV180V","HV200V"};
  double highVolt[numOfHV] = {40.0,50.0,60.0,70.0,80.0,100.0,120.0,140.0,160.0,180.0,200.0};
  // const int numOfHV = 4;
  // string hv[numOfHV] = {"HV70V","HV100V","HV120V","HV140V"};
  // double highVolt[numOfHV] = {70.0,100.0,120.0,140.0};
  TFile *File_InPut[numOfHV];
  TH1F *h_mFstHitsSignal[numOfHV][FST::mFstNumSensorsPerModule];
  TH1F *h_mFstHitsNoise[numOfHV][FST::mFstNumSensorsPerModule];
  TH1F *h_mFstHitsSNRatio[numOfHV][FST::mFstNumSensorsPerModule];

  // read in signal histograms
  for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
  {
    string inputfile = Form("../../output/signal/%s/HVScan/FstQAStudy_%s_%s_Th4.0Tb2Ped2.5Ped3.5_withPed_withCMNCorr.root",mod.c_str(),mod.c_str(),hv[i_hv].c_str());
    File_InPut[i_hv] = TFile::Open(inputfile.c_str());
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      string HistName;
      HistName = Form("h_mFstHitsSignal_Sensor%d",i_sensor);
      h_mFstHitsSignal[i_hv][i_sensor] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());
      HistName = Form("h_mFstHitsNoise_Sensor%d",i_sensor);
      h_mFstHitsNoise[i_hv][i_sensor] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());
      HistName = Form("h_mFstHitsSNRatio_Sensor%d",i_sensor);
      h_mFstHitsSNRatio[i_hv][i_sensor] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());
    }
  }

  // calculate mean signal & S/N
  TGraphAsymmErrors *g_mMeanFstHitsSignal[FST::mFstNumSensorsPerModule];
  TGraphAsymmErrors *g_mMeanFstHitsNoise[FST::mFstNumSensorsPerModule];
  TGraphAsymmErrors *g_mMeanFstHitsSNRatio[FST::mFstNumSensorsPerModule];
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    string GraphName;

    GraphName = Form("g_mMeanFstHitsSignal_Sensor%d",i_sensor);
    g_mMeanFstHitsSignal[i_sensor] = new TGraphAsymmErrors();
    g_mMeanFstHitsSignal[i_sensor]->SetName(GraphName.c_str());

    GraphName = Form("g_mMeanFstHitsNoise_Sensor%d",i_sensor);
    g_mMeanFstHitsNoise[i_sensor] = new TGraphAsymmErrors();
    g_mMeanFstHitsNoise[i_sensor]->SetName(GraphName.c_str());

    GraphName = Form("g_mMeanFstHitsSNRatio_Sensor%d",i_sensor);
    g_mMeanFstHitsSNRatio[i_sensor] = new TGraphAsymmErrors();
    g_mMeanFstHitsSNRatio[i_sensor]->SetName(GraphName.c_str());
  }

  for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
  {
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      double meanSignal  = h_mFstHitsSignal[i_hv][i_sensor]->GetMean();
      double errSignal  = h_mFstHitsSignal[i_hv][i_sensor]->GetMeanError();
      double meanNoise   = h_mFstHitsNoise[i_hv][i_sensor]->GetMean();
      double meanSNRatio = h_mFstHitsSNRatio[i_hv][i_sensor]->GetMean();
      double errSNRatio = h_mFstHitsSNRatio[i_hv][i_sensor]->GetMeanError();
      g_mMeanFstHitsSignal[i_sensor]->SetPoint(i_hv,highVolt[i_hv],meanSignal);
      g_mMeanFstHitsSignal[i_sensor]->SetPointError(i_hv,0.0,0.0,errSignal,errSignal);
      g_mMeanFstHitsNoise[i_sensor]->SetPoint(i_hv,highVolt[i_hv],meanNoise);
      g_mMeanFstHitsSNRatio[i_sensor]->SetPoint(i_hv,highVolt[i_hv],meanSNRatio);
      g_mMeanFstHitsSNRatio[i_sensor]->SetPointError(i_hv,0.0,0.0,errSNRatio,errSNRatio);
    }
  }

  TH1F *h_play = new TH1F("h_play","h_play",200,-9.5,190.5);
  for(int i_bin = 0; i_bin < 200; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("High Voltage (V)");
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetNdivisions(505);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetTitle("ADC");
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetLabelSize(0.04);

  // mean signal
  TCanvas *c_SignalDiff = new TCanvas("c_SignalDiff","c_SignalDiff",10,10,400,1200);
  c_SignalDiff->Divide(1,3);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_SignalDiff->cd(i_pad+1)->SetLeftMargin(0.15);
    c_SignalDiff->cd(i_pad+1)->SetBottomMargin(0.15);
    c_SignalDiff->cd(i_pad+1)->SetTicks(1,1);
    c_SignalDiff->cd(i_pad+1)->SetGrid(0,0);
  }
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    c_SignalDiff->cd(i_sensor+1);
    string title = Form("Signal vs. HV: Sensor%d",i_sensor);
    h_play->SetTitle(title.c_str());
    h_play->SetTitleSize(0.06);
    h_play->GetXaxis()->SetRangeUser(-1.5,160.5);
    h_play->GetYaxis()->SetRangeUser(350.0,500.0);
    h_play->DrawCopy("hE");
    g_mMeanFstHitsSignal[i_sensor]->SetMarkerStyle(20);
    g_mMeanFstHitsSignal[i_sensor]->SetMarkerSize(1.4);
    g_mMeanFstHitsSignal[i_sensor]->SetMarkerColor(kGray+2);
    g_mMeanFstHitsSignal[i_sensor]->Draw("p");
  }
  c_SignalDiff->SaveAs("./figures/c_SignalDiffHVScan.eps");

  // mean noise
  TCanvas *c_NoiseDiff = new TCanvas("c_NoiseDiff","c_NoiseDiff",10,10,400,1200);
  c_NoiseDiff->Divide(1,3);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_NoiseDiff->cd(i_pad+1)->SetLeftMargin(0.15);
    c_NoiseDiff->cd(i_pad+1)->SetBottomMargin(0.15);
    c_NoiseDiff->cd(i_pad+1)->SetTicks(1,1);
    c_NoiseDiff->cd(i_pad+1)->SetGrid(0,0);
  }
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    c_NoiseDiff->cd(i_sensor+1);
    string title = Form("Noise vs. HV: Sensor%d",i_sensor);
    h_play->SetTitle(title.c_str());
    h_play->SetTitleSize(0.06);
    h_play->GetXaxis()->SetRangeUser(-1.5,160.5);
    h_play->GetYaxis()->SetRangeUser(5.0,30.0);
    h_play->DrawCopy("hE");
    g_mMeanFstHitsNoise[i_sensor]->SetMarkerStyle(20);
    g_mMeanFstHitsNoise[i_sensor]->SetMarkerSize(1.4);
    g_mMeanFstHitsNoise[i_sensor]->SetMarkerColor(kGray+2);
    g_mMeanFstHitsNoise[i_sensor]->Draw("p");
  }
  c_NoiseDiff->SaveAs("./figures/c_NoiseDiffHVScan.eps");

  // mean S/N ratio
  TCanvas *c_SNRatioDiff = new TCanvas("c_SNRatioDiff","c_SNRatioDiff",10,10,400,1200);
  c_SNRatioDiff->Divide(1,3);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_SNRatioDiff->cd(i_pad+1)->SetLeftMargin(0.15);
    c_SNRatioDiff->cd(i_pad+1)->SetBottomMargin(0.15);
    c_SNRatioDiff->cd(i_pad+1)->SetTicks(1,1);
    c_SNRatioDiff->cd(i_pad+1)->SetGrid(0,0);
  }
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    c_SNRatioDiff->cd(i_sensor+1);
    string title = Form("S/N vs. HV: Sensor%d",i_sensor);
    h_play->SetTitle(title.c_str());
    h_play->SetTitleSize(0.06);
    h_play->GetXaxis()->SetRangeUser(-1.5,160.5);
    h_play->GetYaxis()->SetTitle("S/N");
    h_play->GetYaxis()->SetRangeUser(15.0,35.0);
    h_play->DrawCopy("hE");
    g_mMeanFstHitsSNRatio[i_sensor]->SetMarkerStyle(20);
    g_mMeanFstHitsSNRatio[i_sensor]->SetMarkerSize(1.4);
    g_mMeanFstHitsSNRatio[i_sensor]->SetMarkerColor(kGray+2);
    g_mMeanFstHitsSNRatio[i_sensor]->Draw("p");
  }
  c_SNRatioDiff->SaveAs("./figures/c_SNRatioDiffHVScan.eps");

  const int markerColor[3] = {kGray+2, 2, 4};
  const int markerStyle[3] = {20, 24, 24};
  // mean signal
  TCanvas *c_Signal = new TCanvas("c_Signal","c_Signal",10,10,800,800);
  c_Signal->cd()->SetLeftMargin(0.15);
  c_Signal->cd()->SetBottomMargin(0.15);
  c_Signal->cd()->SetTicks(1,1);
  c_Signal->cd()->SetGrid(0,0);
  h_play->SetTitle("Hit Signal vs. HV");
  // h_play->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(-4.5,145.5);
  h_play->GetYaxis()->SetTitle("Hit Signal (ADC)");
  h_play->GetYaxis()->SetRangeUser(0.0,500.0);
  h_play->DrawCopy("hE");
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    g_mMeanFstHitsSignal[i_sensor]->SetMarkerStyle(markerStyle[i_sensor]);
    g_mMeanFstHitsSignal[i_sensor]->SetMarkerSize(2.0);
    g_mMeanFstHitsSignal[i_sensor]->SetMarkerColor(markerColor[i_sensor]);
    g_mMeanFstHitsSignal[i_sensor]->Draw("pE same");
  }
  TLegend *leg = new TLegend(0.20,0.65,0.65,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetMargin(0.1);
  // leg->SetNColumns(2);
  leg->AddEntry(g_mMeanFstHitsSignal[0],"Sensor0 (Inner Sector)","P");
  leg->AddEntry(g_mMeanFstHitsSignal[1],"Sensor1 (Outer Sector)","P");
  leg->AddEntry(g_mMeanFstHitsSignal[2],"Sensor2 (Outer Sector)","P");
  leg->Draw("same");
  c_Signal->SaveAs("./figures/c_SignalHVScan.eps");

  // mean noise
  TCanvas *c_Noise = new TCanvas("c_Noise","c_Noise",10,10,800,800);
  c_Noise->cd()->SetLeftMargin(0.15);
  c_Noise->cd()->SetBottomMargin(0.15);
  c_Noise->cd()->SetTicks(1,1);
  c_Noise->cd()->SetGrid(0,0);
  h_play->SetTitle("Hit Noise vs. HV");
  // h_play->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(-4.5,145.5);
  h_play->GetYaxis()->SetTitle("Hit Noise (ADC)");
  h_play->GetYaxis()->SetRangeUser(5.0,30.0);
  h_play->DrawCopy("hE");
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    g_mMeanFstHitsNoise[i_sensor]->SetMarkerStyle(markerStyle[i_sensor]);
    g_mMeanFstHitsNoise[i_sensor]->SetMarkerSize(2.0);
    g_mMeanFstHitsNoise[i_sensor]->SetMarkerColor(markerColor[i_sensor]);
    g_mMeanFstHitsNoise[i_sensor]->Draw("p same");
  }
  leg->Draw("same");
  c_Noise->SaveAs("./figures/c_NoiseHVScan.eps");

  // mean S/N ratio
  TCanvas *c_SNRatio = new TCanvas("c_SNRatio","c_SNRatio",10,10,800,800);
  c_SNRatio->cd()->SetLeftMargin(0.15);
  c_SNRatio->cd()->SetBottomMargin(0.15);
  c_SNRatio->cd()->SetTicks(1,1);
  c_SNRatio->cd()->SetGrid(0,0);
  h_play->SetTitle("Hit S/N vs. HV");
  // h_play->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(-1.5,160.5);
  h_play->GetYaxis()->SetTitle("S/N");
  h_play->GetYaxis()->SetRangeUser(0.0,35.0);
  h_play->DrawCopy("hE");
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    g_mMeanFstHitsSNRatio[i_sensor]->SetMarkerStyle(markerStyle[i_sensor]);
    g_mMeanFstHitsSNRatio[i_sensor]->SetMarkerSize(2.0);
    g_mMeanFstHitsSNRatio[i_sensor]->SetMarkerColor(markerColor[i_sensor]);
    g_mMeanFstHitsSNRatio[i_sensor]->Draw("p same");
  }
  leg->Draw("same");
  c_SNRatio->SaveAs("./figures/c_SNRatioHVScan.eps");
}
