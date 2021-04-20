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

void plotSensorDataNoiseDate(string mod = "Mod01", string hv = "HV70V")
{
  const int defTimeBin = 0;
  std::vector<string> date;

  cout << "Read in run list->" << endl;
  string inputfile = Form("../../list/FST/data/%s/Date_%s_%s.list",mod.c_str(),mod.c_str(),hv.c_str());
  cout << "Open run list: " << inputfile.c_str() << endl;
  ifstream file_date (inputfile.c_str());
  if( !file_date.is_open() )
  {
    cout << "Abort. Fail to read in channel map: " << inputfile.c_str() << endl;
    return false;
  }

  cout << "reading in run list:" << endl;
  string runlist;
  while (file_date >> runlist)
  {
    cout << "date = " << runlist.c_str() << endl;
    date.push_back(runlist);
  }
  file_date.close();
  
  const int numOfDate = date.size();

  TFile *File_InPut[numOfDate];
  TH1F *h_mPedSigma[numOfDate][FST::numRStrip][FST::numTBins];
  TH1F *h_mCMNSigma[numOfDate][FST::numRStrip][FST::numTBins];
  TH1F *h_mRanSigma[numOfDate][FST::numRStrip][FST::numTBins];

  // read in noise histograms
  for(int i_date = 0; i_date < numOfDate; ++i_date)
  {
    string input_Sensor = Form("../../output/noise/%s/FstDataNoise_%s_%s_%s.root",mod.c_str(),mod.c_str(),hv.c_str(),date[i_date].c_str());
    File_InPut[i_date] = TFile::Open(input_Sensor.c_str());
    for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
    {
      for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
      {
	std::string HistName;
	HistName = Form("h_mPedSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
	h_mPedSigma[i_date][i_rstrip][i_tb] = (TH1F*)File_InPut[i_date]->Get(HistName.c_str());
	HistName = Form("h_mCMNSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
	h_mCMNSigma[i_date][i_rstrip][i_tb] = (TH1F*)File_InPut[i_date]->Get(HistName.c_str());
	HistName = Form("h_mRanSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
	h_mRanSigma[i_date][i_rstrip][i_tb] = (TH1F*)File_InPut[i_date]->Get(HistName.c_str());
      }
    }
  }

  // calculate mean noise
  // TGraphAsymmErrors *g_mMeanPedSigma[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
  // TGraphAsymmErrors *g_mMeanRanSigma[FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
  TH1F *h_mMeanPedSigma[FST::mFstNumSensorsPerModule];
  TH1F *h_mMeanCMNSigma[FST::mFstNumSensorsPerModule];
  TH1F *h_mMeanRanSigma[FST::mFstNumSensorsPerModule];
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    std::string HistName;
    HistName = Form("h_mMeanPedSigma_Sensor%d",i_sensor);
    h_mMeanPedSigma[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),numOfDate,-0.5,numOfDate-0.5);

    HistName = Form("h_mMeanCMNSigma_Sensor%d",i_sensor);
    h_mMeanCMNSigma[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),numOfDate,-0.5,numOfDate-0.5);

    HistName = Form("h_mMeanRanSigma_Sensor%d",i_sensor);
    h_mMeanRanSigma[i_sensor] = new TH1F(HistName.c_str(),HistName.c_str(),numOfDate,-0.5,numOfDate-0.5);
  }

  for(int i_date = 0; i_date < numOfDate; ++i_date)
  {
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      int pedCounts = 0;
      double meanPedSigam = 0.0;
      int cmnCounts = 0;
      double meanCMNSigam = 0.0;
      int ranCounts = 0;
      double meanRanSigam = 0.0;
      for(int i_rstrip = 0; i_rstrip < FST::mFstNumRstripPerSensor; ++i_rstrip)
      {
	if(i_sensor == 0)
	{
	  for(int i_ch = 0; i_ch < FST::numChannels; ++i_ch)
	  { // calculate mean total and random noise
	    if(h_mPedSigma[i_date][i_rstrip][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanPedSigam += h_mPedSigma[i_date][i_rstrip][defTimeBin]->GetBinContent(i_ch+1);
	      pedCounts++;
	    }
	    if(h_mCMNSigma[i_date][i_rstrip][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanCMNSigam += h_mCMNSigma[i_date][i_rstrip][defTimeBin]->GetBinContent(i_ch+1);
	      cmnCounts++;
	    }
	    if(h_mRanSigma[i_date][i_rstrip][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanRanSigam += h_mRanSigma[i_date][i_rstrip][defTimeBin]->GetBinContent(i_ch+1);
	      ranCounts++;
	    }
	  }
	}
	if(i_sensor == 1)
	{
	  for(int i_ch = 0; i_ch < 64; ++i_ch)
	  { // calculate mean total and random noise
	    if(h_mPedSigma[i_date][i_rstrip+4][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanPedSigam += h_mPedSigma[i_date][i_rstrip+4][defTimeBin]->GetBinContent(i_ch+1);
	      pedCounts++;
	    }
	    if(h_mCMNSigma[i_date][i_rstrip+4][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanCMNSigam += h_mCMNSigma[i_date][i_rstrip+4][defTimeBin]->GetBinContent(i_ch+1);
	      cmnCounts++;
	    }
	    if(h_mRanSigma[i_date][i_rstrip+4][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanRanSigam += h_mRanSigma[i_date][i_rstrip+4][defTimeBin]->GetBinContent(i_ch+1);
	      ranCounts++;
	    }
	  }
	}
	if(i_sensor == 2)
	{
	  for(int i_ch = 64; i_ch < 128; ++i_ch)
	  { // calculate mean total and random noise
	    if(h_mPedSigma[i_date][i_rstrip+4][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanPedSigam += h_mPedSigma[i_date][i_rstrip+4][defTimeBin]->GetBinContent(i_ch+1);
	      pedCounts++;
	    }
	    if(h_mCMNSigma[i_date][i_rstrip+4][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanCMNSigam += h_mCMNSigma[i_date][i_rstrip+4][defTimeBin]->GetBinContent(i_ch+1);
	      cmnCounts++;
	    }
	    if(h_mRanSigma[i_date][i_rstrip+4][defTimeBin]->GetBinError(i_ch+1) > 0)
	    {
	      meanRanSigam += h_mRanSigma[i_date][i_rstrip+4][defTimeBin]->GetBinContent(i_ch+1);
	      ranCounts++;
	    }
	  }
	}
      }
      // cout << "i_sensor = " << i_sensor << ", pedCounts = " << pedCounts << ", cmnCounts = " << cmnCounts << ", ranCounts = " << ranCounts << endl;
      // cout << "i_sensor = " << i_sensor << ", meanCMNSigam/cmnCounts = " << meanCMNSigam/cmnCounts <<endl;
      h_mMeanPedSigma[i_sensor]->SetBinContent(i_date+1,meanPedSigam/pedCounts);
      h_mMeanCMNSigma[i_sensor]->SetBinContent(i_date+1,meanCMNSigam/cmnCounts);
      h_mMeanRanSigma[i_sensor]->SetBinContent(i_date+1,meanRanSigam/ranCounts);
    }
  }

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
    string title = Form("Noise vs. Date: Sensor%d",i_sensor);
    h_mMeanPedSigma[i_sensor]->SetTitle(title.c_str());
    h_mMeanPedSigma[i_sensor]->SetStats(0);
    for(int i_date = 0; i_date < numOfDate; ++i_date)
    {
      h_mMeanPedSigma[i_sensor]->GetXaxis()->SetBinLabel(i_date+1,date[i_date].c_str());
    }
    h_mMeanPedSigma[i_sensor]->GetXaxis()->SetLabelSize(0.06);
    h_mMeanPedSigma[i_sensor]->GetYaxis()->SetRangeUser(0.0,60.0);
    h_mMeanPedSigma[i_sensor]->SetMarkerStyle(20);
    h_mMeanPedSigma[i_sensor]->SetMarkerSize(1.0);
    h_mMeanPedSigma[i_sensor]->SetMarkerColor(kGray+2);
    h_mMeanPedSigma[i_sensor]->Draw("p");
    h_mMeanCMNSigma[i_sensor]->SetMarkerStyle(24);
    h_mMeanCMNSigma[i_sensor]->SetMarkerSize(1.0);
    h_mMeanCMNSigma[i_sensor]->SetMarkerColor(4);
    h_mMeanCMNSigma[i_sensor]->Draw("P same");
    h_mMeanRanSigma[i_sensor]->SetMarkerStyle(24);
    h_mMeanRanSigma[i_sensor]->SetMarkerSize(1.0);
    h_mMeanRanSigma[i_sensor]->SetMarkerColor(2);
    h_mMeanRanSigma[i_sensor]->Draw("P same");
    TLegend *leg = new TLegend(0.20,0.65,0.75,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry(h_mMeanPedSigma[i_sensor],"Total Noise","P");
    leg->AddEntry(h_mMeanCMNSigma[i_sensor],"Common Mode Noise","P");
    leg->AddEntry(h_mMeanRanSigma[i_sensor],"Random Noise","P");
    leg->Draw("same");
  }
  string FigName = Form("./figures/%s/DataNoiseDate_%s_%s.pdf",mod.c_str(),mod.c_str(),hv.c_str());
  c_NoiseDiff->SaveAs(FigName.c_str());
}
