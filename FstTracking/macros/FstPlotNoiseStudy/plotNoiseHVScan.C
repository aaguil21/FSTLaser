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

void plotNoiseHVScan(string mod = "Mod04")
{
  const int defTimeBin = 0;

  // const int numOfHV = 12; // Mod01
  // string hv[numOfHV] = {"HVOff","HV40V","HV50V","HV60V","HV70V","HV80V","HV100V","HV120V","HV140V","HV160V","HV180V","HV200V"};
  // int highVolt[numOfHV] = {0,40,50,60,70,80,100,120,140,160,180,200};

  const int numOfHV = 16; // Mod03
  string hv[numOfHV] = {"HV0V","HV5V","HV10V","HV20V","HV30V","HV40V","HV50V","HV60V","HV70V","HV80V","HV90V","HV100V","HV110V","HV120V","HV130V","HV140V"};
  int highVolt[numOfHV] = {0,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140};

  TFile *File_InPut[numOfHV];
  TH1F *h_mPedSigma[numOfHV][FST::numRStrip][FST::numTBins];
  TH1F *h_mRanSigma[numOfHV][FST::numRStrip][FST::numTBins];

  // read in noise histograms
  for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
  {
    string input_Sensor = Form("../../output/noise/%s/HVScan/FstPedNoise_%s_%s.root",mod.c_str(),mod.c_str(),hv[i_hv].c_str());
    File_InPut[i_hv] = TFile::Open(input_Sensor.c_str());
    for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
    {
      for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
      {
	std::string HistName;
	HistName = Form("h_mPedSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
	h_mPedSigma[i_hv][i_rstrip][i_tb] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());
	HistName = Form("h_mRanSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
	h_mRanSigma[i_hv][i_rstrip][i_tb] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());
      }
    }
  }

  // calculate mean noise
  // TGraphAsymmErrors *g_mMeanPedSigma[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
  // TGraphAsymmErrors *g_mMeanRanSigma[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
  TGraphAsymmErrors *g_mMeanPedSigma[FST::mFstNumSensorsPerModule];
  TGraphAsymmErrors *g_mMeanRanSigma[FST::mFstNumSensorsPerModule];
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    string GraphName;

    GraphName = Form("g_mMeanPedSigma_Sensor%d",i_sensor);
    g_mMeanPedSigma[i_sensor] = new TGraphAsymmErrors();
    g_mMeanPedSigma[i_sensor]->SetName(GraphName.c_str());

    GraphName = Form("g_mMeanRanSigma_Sensor%d",i_sensor);
    g_mMeanRanSigma[i_sensor] = new TGraphAsymmErrors();
    g_mMeanRanSigma[i_sensor]->SetName(GraphName.c_str());
  }

  for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
  {
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      int pedCounts = 0;
      double meanPedSigam = 0.0;
      int ranCounts = 0;
      double meanRanSigam = 0.0;
      for(int i_rstrip = 0; i_rstrip < FST::mFstNumRstripPerSensor; ++i_rstrip)
      {
	if(i_sensor == 0)
	{
	  for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	  { // calculate mean total and random noise
	    if(h_mPedSigma[i_hv][i_rstrip][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanPedSigam += h_mPedSigma[i_hv][i_rstrip][defTimeBin]->GetBinContent(i_ch+1);
	      pedCounts++;
	    }
	    if(h_mRanSigma[i_hv][i_rstrip][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanRanSigam += h_mRanSigma[i_hv][i_rstrip][defTimeBin]->GetBinContent(i_ch+1);
	      ranCounts++;
	    }
	  }
	}
	if(i_sensor == 1)
	{
	  for(int i_ch = 0; i_ch < 64; ++i_ch)
	  { // calculate mean total and random noise
	    if(h_mPedSigma[i_hv][i_rstrip+4][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanPedSigam += h_mPedSigma[i_hv][i_rstrip+4][defTimeBin]->GetBinContent(i_ch+1);
	      pedCounts++;
	    }
	    if(h_mRanSigma[i_hv][i_rstrip+4][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanRanSigam += h_mRanSigma[i_hv][i_rstrip+4][defTimeBin]->GetBinContent(i_ch+1);
	      ranCounts++;
	    }
	  }
	}
	if(i_sensor == 2)
	{
	  for(int i_ch = 64; i_ch < 128; ++i_ch)
	  { // calculate mean total and random noise
	    if(h_mPedSigma[i_hv][i_rstrip+4][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanPedSigam += h_mPedSigma[i_hv][i_rstrip+4][defTimeBin]->GetBinContent(i_ch+1);
	      pedCounts++;
	    }
	    if(h_mRanSigma[i_hv][i_rstrip+4][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanRanSigam += h_mRanSigma[i_hv][i_rstrip+4][defTimeBin]->GetBinContent(i_ch+1);
	      ranCounts++;
	    }
	  }
	}
      }
      // cout << "i_sensor = " << i_sensor << ", pedCounts = " << pedCounts << ", ranCounts = " << ranCounts << endl;
      g_mMeanPedSigma[i_sensor]->SetPoint(i_hv,highVolt[i_hv],meanPedSigam/pedCounts);
      g_mMeanRanSigma[i_sensor]->SetPoint(i_hv,highVolt[i_hv],meanRanSigam/ranCounts);
    }
  }

  TH1F *h_play = new TH1F("h_play","h_play",250,-9.5,240.5);
  for(int i_bin = 0; i_bin < 250; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetStats(0);
  h_play->GetXaxis()->SetRangeUser(-5.5,220.5);
  h_play->GetXaxis()->SetTitle("High Voltage (V)");
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetNdivisions(505);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetRangeUser(-1.5,50.0);
  h_play->GetYaxis()->SetTitle("ADC");
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetLabelSize(0.04);

  const int markerColor[3] = {kGray+2, 2, 4};
  const int markerStyle[3] = {20, 24, 24};
  TCanvas *c_Noise = new TCanvas("c_Noise","c_Noise",10,10,800,800);
  c_Noise->cd()->SetLeftMargin(0.15);
  c_Noise->cd()->SetBottomMargin(0.15);
  c_Noise->cd()->SetTicks(1,1);
  c_Noise->cd()->SetGrid(0,0);
  string title = Form("Single Channel Random Noise vs. HV @ %s",mod.c_str());
  h_play->SetTitle(title.c_str());
  h_play->GetXaxis()->SetRangeUser(-5.5,220.5);
  h_play->GetYaxis()->SetTitle("Singal Channel Noise (ADC)");
  h_play->DrawCopy("hE");
  PlotLine(-5.5, 220.5,0.0,0.0,1,1,2);
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    g_mMeanRanSigma[i_sensor]->SetMarkerStyle(markerStyle[i_sensor]);
    g_mMeanRanSigma[i_sensor]->SetMarkerSize(2.0);
    g_mMeanRanSigma[i_sensor]->SetMarkerColor(markerColor[i_sensor]);
    g_mMeanRanSigma[i_sensor]->Draw("P same");
  }
  TLegend *leg = new TLegend(0.40,0.65,0.85,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetMargin(0.1);
  // leg->SetNColumns(2);
  leg->AddEntry(g_mMeanRanSigma[0],"Sensor0 (Inner Sector)","P");
  leg->AddEntry(g_mMeanRanSigma[1],"Sensor1 (Outer Sector)","P");
  leg->AddEntry(g_mMeanRanSigma[2],"Sensor2 (Outer Sector)","P");
  leg->Draw("same");

  string FigName = Form("./figures/%s/NoiseHVScan_%s_SingelChannel.eps",mod.c_str(),mod.c_str());
  c_Noise->SaveAs(FigName.c_str());
}
