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

void plotSignalHVScan(string mod = "Mod01")
{
  const int numOfHV = 11;
  string hv[numOfHV] = {"HV40V","HV50V","HV60V","HV70V","HV80V","HV100V","HV120V","HV140V","HV160V","HV180V","HV200V"};
  double highVolt[numOfHV] = {40.0,50.0,60.0,70.0,80.0,100.0,120.0,140.0,160.0,180.0,200.0};
  TFile *File_InPut[numOfHV];
  TH1F *h_mFstScanClustersSignal[numOfHV][FST::mFstNumSensorsPerModule];
  TH1F *h_mFstScanClustersNoise[numOfHV][FST::mFstNumSensorsPerModule];
  TH1F *h_mFstScanClustersSNRatio[numOfHV][FST::mFstNumSensorsPerModule];

  // read in signal histograms
  for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
  {
    string inputfile = Form("../../output/signal/%s/HVScan/FstQAStudy_%s_%s_Th4.0Tb2Ped2.5Ped3.5_withPed_withCMNCorr.root",mod.c_str(),mod.c_str(),hv[i_hv].c_str());
    File_InPut[i_hv] = TFile::Open(inputfile.c_str());
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      string HistName;
      HistName = Form("h_mFstScanClustersSignal_Sensor%d",i_sensor);
      h_mFstScanClustersSignal[i_hv][i_sensor] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());
      HistName = Form("h_mFstScanClustersNoise_Sensor%d",i_sensor);
      h_mFstScanClustersNoise[i_hv][i_sensor] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());
      HistName = Form("h_mFstScanClustersSNRatio_Sensor%d",i_sensor);
      h_mFstScanClustersSNRatio[i_hv][i_sensor] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());
    }
  }

  // calculate mean signal & S/N
  TGraphAsymmErrors *g_mMeanFstScanClustersSignal[FST::mFstNumSensorsPerModule];
  TGraphAsymmErrors *g_mMeanFstScanClustersNoise[FST::mFstNumSensorsPerModule];
  TGraphAsymmErrors *g_mMeanFstScanClustersSNRatio[FST::mFstNumSensorsPerModule];
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    string GraphName;

    GraphName = Form("g_mMeanFstScanClustersSignal_Sensor%d",i_sensor);
    g_mMeanFstScanClustersSignal[i_sensor] = new TGraphAsymmErrors();
    g_mMeanFstScanClustersSignal[i_sensor]->SetName(GraphName.c_str());

    GraphName = Form("g_mMeanFstScanClustersNoise_Sensor%d",i_sensor);
    g_mMeanFstScanClustersNoise[i_sensor] = new TGraphAsymmErrors();
    g_mMeanFstScanClustersNoise[i_sensor]->SetName(GraphName.c_str());

    GraphName = Form("g_mMeanFstScanClustersSNRatio_Sensor%d",i_sensor);
    g_mMeanFstScanClustersSNRatio[i_sensor] = new TGraphAsymmErrors();
    g_mMeanFstScanClustersSNRatio[i_sensor]->SetName(GraphName.c_str());
  }

  for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
  {
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      double meanSignal  = h_mFstScanClustersSignal[i_hv][i_sensor]->GetMean();
      double errSignal  = h_mFstScanClustersSignal[i_hv][i_sensor]->GetMeanError();
      double meanNoise   = h_mFstScanClustersNoise[i_hv][i_sensor]->GetMean();
      double meanSNRatio = h_mFstScanClustersSNRatio[i_hv][i_sensor]->GetMean();
      double errSNRatio = h_mFstScanClustersSNRatio[i_hv][i_sensor]->GetMeanError();
      g_mMeanFstScanClustersSignal[i_sensor]->SetPoint(i_hv,highVolt[i_hv],meanSignal);
      g_mMeanFstScanClustersSignal[i_sensor]->SetPointError(i_hv,0.0,0.0,errSignal,errSignal);
      g_mMeanFstScanClustersNoise[i_sensor]->SetPoint(i_hv,highVolt[i_hv],meanNoise);
      g_mMeanFstScanClustersSNRatio[i_sensor]->SetPoint(i_hv,highVolt[i_hv],meanSNRatio);
      g_mMeanFstScanClustersSNRatio[i_sensor]->SetPointError(i_hv,0.0,0.0,errSNRatio,errSNRatio);
    }
  }

  TH1F *h_play = new TH1F("h_play","h_play",250,-9.5,240.5);
  for(int i_bin = 0; i_bin < 250; ++i_bin)
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

  const int markerColor[3] = {kGray+2, 2, 4};
  const int markerStyle[3] = {20, 24, 24};
  // mean signal
  TCanvas *c_Signal = new TCanvas("c_Signal","c_Signal",10,10,800,800);
  c_Signal->cd()->SetLeftMargin(0.15);
  c_Signal->cd()->SetBottomMargin(0.15);
  c_Signal->cd()->SetTicks(1,1);
  c_Signal->cd()->SetGrid(0,0);
  h_play->SetTitle("Cluster Signal vs. HV");
  h_play->GetXaxis()->SetRangeUser(-4.5,220.5);
  h_play->GetYaxis()->SetTitle("Cluster Signal (ADC)");
  h_play->GetYaxis()->SetRangeUser(200.0,700.0);
  h_play->DrawCopy("hE");
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    g_mMeanFstScanClustersSignal[i_sensor]->SetMarkerStyle(markerStyle[i_sensor]);
    g_mMeanFstScanClustersSignal[i_sensor]->SetMarkerSize(2.0);
    g_mMeanFstScanClustersSignal[i_sensor]->SetMarkerColor(markerColor[i_sensor]);
    g_mMeanFstScanClustersSignal[i_sensor]->Draw("p same");
  }
  TLegend *leg = new TLegend(0.20,0.70,0.65,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetMargin(0.1);
  // leg->SetNColumns(2);
  leg->AddEntry(g_mMeanFstScanClustersSignal[0],"Sensor0 (Inner Sector)","P");
  leg->AddEntry(g_mMeanFstScanClustersSignal[1],"Sensor1 (Outer Sector)","P");
  leg->AddEntry(g_mMeanFstScanClustersSignal[2],"Sensor2 (Outer Sector)","P");
  leg->Draw("same");

  string FigName = Form("./figures/%s/HVScan/SignalHVScan_%s_ScanCluster.eps",mod.c_str(),mod.c_str());
  c_Signal->SaveAs(FigName.c_str());

  // mean noise
  TCanvas *c_Noise = new TCanvas("c_Noise","c_Noise",10,10,800,800);
  c_Noise->cd()->SetLeftMargin(0.15);
  c_Noise->cd()->SetBottomMargin(0.15);
  c_Noise->cd()->SetTicks(1,1);
  c_Noise->cd()->SetGrid(0,0);
  h_play->SetTitle("Cluster Noise vs. HV");
  // h_play->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(-4.5,200.5);
  h_play->GetYaxis()->SetTitle("Cluster Noise (ADC)");
  h_play->GetYaxis()->SetRangeUser(5.0,30.0);
  h_play->DrawCopy("hE");
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    g_mMeanFstScanClustersNoise[i_sensor]->SetMarkerStyle(markerStyle[i_sensor]);
    g_mMeanFstScanClustersNoise[i_sensor]->SetMarkerSize(2.0);
    g_mMeanFstScanClustersNoise[i_sensor]->SetMarkerColor(markerColor[i_sensor]);
    g_mMeanFstScanClustersNoise[i_sensor]->Draw("p same");
  }
  leg->Draw("same");

  FigName = Form("./figures/%s/HVScan/NoiseHVScan_%s_ScanCluster.eps",mod.c_str(),mod.c_str());
  c_Noise->SaveAs(FigName.c_str());

  // mean S/N ratio
  TCanvas *c_SNRatio = new TCanvas("c_SNRatio","c_SNRatio",10,10,800,800);
  c_SNRatio->cd()->SetLeftMargin(0.15);
  c_SNRatio->cd()->SetBottomMargin(0.15);
  c_SNRatio->cd()->SetTicks(1,1);
  c_SNRatio->cd()->SetGrid(0,0);
  h_play->SetTitle("Cluster S/N vs. HV");
  // h_play->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(-1.5,220.5);
  h_play->GetYaxis()->SetTitle("S/N");
  h_play->GetYaxis()->SetRangeUser(5.0,50.0);
  h_play->DrawCopy("hE");
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    g_mMeanFstScanClustersSNRatio[i_sensor]->SetMarkerStyle(markerStyle[i_sensor]);
    g_mMeanFstScanClustersSNRatio[i_sensor]->SetMarkerSize(2.0);
    g_mMeanFstScanClustersSNRatio[i_sensor]->SetMarkerColor(markerColor[i_sensor]);
    g_mMeanFstScanClustersSNRatio[i_sensor]->Draw("p same");
  }
  leg->Draw("same");

  FigName = Form("./figures/%s/HVScan/SNRatioHVScan_%s_ScanCluster.eps",mod.c_str(),mod.c_str());
  c_SNRatio->SaveAs(FigName.c_str());
}
