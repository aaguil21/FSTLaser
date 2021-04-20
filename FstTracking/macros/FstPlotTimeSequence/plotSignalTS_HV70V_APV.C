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

float ErrorAdd(float x, float y)
{
  return sqrt(x*x+y*y);
}

float ErrTimes(float x, float y, float dx, float dy)
{
  return x*y*ErrorAdd(dx/x,dy/y);
}

float ErrDiv(float x, float y, float dx, float dy)
{
  return x/y*ErrorAdd(dx/x,dy/y);
}

void plotSignalTS_HV70V_APV(int apv = 1, int phibin = 1)
{
  const int numOfTS = 6;
  string date[numOfTS] = {"0318","0406","0410","0411","0414","0416"};

  TH1F *h_mMeanSignalHits_Apv[4];
  TH1F *h_mMeanNoiseHits_Apv[4];
  TH1F *h_mMeanSNRatioHits_Apv[4];
  TH1F *h_mMeanMaxTbHits_Apv[4];
  TH1F *h_mMeanSignalClusters_Apv[4];
  for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
  {
    string HistName;
    HistName = Form("h_mMeanSignalHits_Apv%d_Rstrip%d_Phi%d",apv+4,i_rstrip,phibin);
    h_mMeanSignalHits_Apv[i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
    HistName = Form("h_mMeanNoiseHits_Apv%d_Rstrip%d_Phi%d",apv+4,i_rstrip,phibin);
    h_mMeanNoiseHits_Apv[i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
    HistName = Form("h_mMeanSNRatioHits_Apv%d_Rstrip%d_Phi%d",apv+4,i_rstrip,phibin);
    h_mMeanSNRatioHits_Apv[i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
    HistName = Form("h_mMeanMaxTbHits_Apv%d_Rstrip%d_Phi%d",apv+4,i_rstrip,phibin);
    h_mMeanMaxTbHits_Apv[i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
    HistName = Form("h_mMeanSignalClusters_Apv%d_Rstrip%d_Phi%d",apv+4,i_rstrip,phibin);
    h_mMeanSignalClusters_Apv[i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
  }

  TFile *File_InPut[numOfTS];
  TH1F *h_mSignalHits_Apv[numOfTS][4];
  TH1F *h_mNoiseHits_Apv[numOfTS][4];
  TH1F *h_mSNRatioHits_Apv[numOfTS][4];
  TH1F *h_mMaxTbHits_Apv[numOfTS][4];
  TH1F *h_mFstSimpleClustersSignal_Apv[numOfTS][4];
  for(int i_ts = 0; i_ts < numOfTS; ++i_ts)
  {
    string inputfile = Form("../../output/timesequence/FstQAStudy_HV70V_Th4.5Tb3_withPed_woCMNCorr_%s2020.root",date[i_ts].c_str());
    cout << "open inputfile: " << inputfile.c_str() << endl;
    File_InPut[i_ts] = TFile::Open(inputfile.c_str());

    for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
    {
      string HistName;
      HistName = Form("h_mSignalHits_Apv%d_Rstrip%d_Phi%d",apv+4,i_rstrip,phibin);
      h_mSignalHits_Apv[i_ts][i_rstrip] = (TH1F*)File_InPut[i_ts]->Get(HistName.c_str());
      HistName = Form("h_mNoiseHits_Apv%d_Rstrip%d_Phi%d",apv+4,i_rstrip,phibin);
      h_mNoiseHits_Apv[i_ts][i_rstrip] = (TH1F*)File_InPut[i_ts]->Get(HistName.c_str());
      HistName = Form("h_mSNRatioHits_Apv%d_Rstrip%d_Phi%d",apv+4,i_rstrip,phibin);
      h_mSNRatioHits_Apv[i_ts][i_rstrip] = (TH1F*)File_InPut[i_ts]->Get(HistName.c_str());
      HistName = Form("h_mMaxTbHits_Apv%d_Rstrip%d_Phi%d",apv+4,i_rstrip,phibin);
      h_mMaxTbHits_Apv[i_ts][i_rstrip] = (TH1F*)File_InPut[i_ts]->Get(HistName.c_str());
      HistName = Form("h_mFstSimpleClustersSignal_Apv%d_Rstrip%d_Phi%d",apv+4,i_rstrip,phibin);
      h_mFstSimpleClustersSignal_Apv[i_ts][i_rstrip] = (TH1F*)File_InPut[i_ts]->Get(HistName.c_str());
    }
  }

  string outputname = Form("./figures/SignalTS_70V_Apv%d_Phi%d.pdf",apv+5,phibin);

  TCanvas *c_Signal = new TCanvas("c_Signal","c_Signal",10,10,1000,1000);
  c_Signal->Divide(5,5);
  for(int i_pad = 0;i_pad < 25; ++i_pad)
  {
    c_Signal->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Signal->cd(i_pad+1)->SetRightMargin(0.15);
    c_Signal->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Signal->cd(i_pad+1)->SetTicks(1,1);
    c_Signal->cd(i_pad+1)->SetGrid(0,0);
  }

  string output_start = Form("./figures/SignalTS_70V_Apv%d_Phi%d.pdf[",apv+5,phibin);
  c_Signal->Print(output_start.c_str()); // open pdf file

  int numOfUsedTS = 0;
  while(numOfUsedTS < numOfTS)
  {
    for(int i_pad = 0;i_pad < 25; ++i_pad)
    {
      c_Signal->cd(i_pad+1)->Clear();
    }
    for(int i_ts = 0; i_ts < 5; ++i_ts)
    {
      if(numOfUsedTS == numOfTS) continue;
      c_Signal->cd(5*i_ts+1);
      {
	TLegend *leg_FST = new TLegend(0.4,0.6,0.85,0.8);
	leg_FST->SetBorderSize(0);
	leg_FST->SetFillColor(10);
	for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
	{
	  string title = Form("FST Hits @ Date %s",date[numOfUsedTS].c_str());
	  h_mSignalHits_Apv[numOfUsedTS][i_rstrip]->SetStats(0);
	  h_mSignalHits_Apv[numOfUsedTS][i_rstrip]->SetTitle(title.c_str());
	  h_mSignalHits_Apv[numOfUsedTS][i_rstrip]->GetXaxis()->SetTitle("ADC");
	  h_mSignalHits_Apv[numOfUsedTS][i_rstrip]->GetXaxis()->SetTitleSize(0.06);
	  h_mSignalHits_Apv[numOfUsedTS][i_rstrip]->SetLineColor(i_rstrip+1);

	  if(i_rstrip == 0) h_mSignalHits_Apv[numOfUsedTS][i_rstrip]->Draw();
	  else h_mSignalHits_Apv[numOfUsedTS][i_rstrip]->Draw("same");

	  string LegName = Form("R_strip %d: Mean %1.2f",i_rstrip,h_mSignalHits_Apv[numOfUsedTS][i_rstrip]->GetMean());
	  leg_FST->AddEntry(h_mSignalHits_Apv[numOfUsedTS][i_rstrip],LegName.c_str(),"L");

	  double mean_orig    = h_mSignalHits_Apv[0][i_rstrip]->GetMean();
	  double err_orig     = h_mSignalHits_Apv[0][i_rstrip]->GetMeanError();
	  double mean_current = h_mSignalHits_Apv[numOfUsedTS][i_rstrip]->GetMean();
	  double err_current  = h_mSignalHits_Apv[numOfUsedTS][i_rstrip]->GetMeanError();
	  double ratio = mean_current/mean_orig;
	  double error = ErrDiv(mean_current, mean_orig, err_current, err_orig);
	  h_mMeanSignalHits_Apv[i_rstrip]->SetBinContent(numOfUsedTS+1,ratio);
	  h_mMeanSignalHits_Apv[i_rstrip]->SetBinError(numOfUsedTS+1,error);
	  // h_mMeanSignalHits_Apv[i_rstrip]->SetBinContent(numOfUsedTS+1,h_mSignalHits_Apv[numOfUsedTS][i_rstrip]->GetMean());
	}
	leg_FST->Draw("same");
      }

      c_Signal->cd(5*i_ts+2);
      {
	TLegend *leg_FST = new TLegend(0.4,0.6,0.85,0.8);
	leg_FST->SetBorderSize(0);
	leg_FST->SetFillColor(10);
	for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
	{
	  string title = Form("FST Clusters @ Date %s",date[numOfUsedTS].c_str());
	  h_mFstSimpleClustersSignal_Apv[numOfUsedTS][i_rstrip]->SetStats(0);
	  h_mFstSimpleClustersSignal_Apv[numOfUsedTS][i_rstrip]->SetTitle(title.c_str());
	  h_mFstSimpleClustersSignal_Apv[numOfUsedTS][i_rstrip]->GetXaxis()->SetTitle("ADC");
	  h_mFstSimpleClustersSignal_Apv[numOfUsedTS][i_rstrip]->GetXaxis()->SetTitleSize(0.06);
	  h_mFstSimpleClustersSignal_Apv[numOfUsedTS][i_rstrip]->SetLineColor(i_rstrip+1);

	  if(i_rstrip == 0) h_mFstSimpleClustersSignal_Apv[numOfUsedTS][i_rstrip]->Draw();
	  else h_mFstSimpleClustersSignal_Apv[numOfUsedTS][i_rstrip]->Draw("same");

	  string LegName = Form("R_strip %d: Mean %1.2f",i_rstrip,h_mFstSimpleClustersSignal_Apv[numOfUsedTS][i_rstrip]->GetMean());
	  leg_FST->AddEntry(h_mFstSimpleClustersSignal_Apv[numOfUsedTS][i_rstrip],LegName.c_str(),"L");

	  double mean_orig    = h_mFstSimpleClustersSignal_Apv[0][i_rstrip]->GetMean();
	  double err_orig     = h_mFstSimpleClustersSignal_Apv[0][i_rstrip]->GetMeanError();
	  double mean_current = h_mFstSimpleClustersSignal_Apv[numOfUsedTS][i_rstrip]->GetMean();
	  double err_current  = h_mFstSimpleClustersSignal_Apv[numOfUsedTS][i_rstrip]->GetMeanError();
	  double ratio = mean_current/mean_orig;
	  double error = ErrDiv(mean_current, mean_orig, err_current, err_orig);
	  h_mMeanSignalClusters_Apv[i_rstrip]->SetBinContent(numOfUsedTS+1,ratio);
	  h_mMeanSignalClusters_Apv[i_rstrip]->SetBinError(numOfUsedTS+1,error);
	  // h_mMeanSignalClusters_Apv[i_rstrip]->SetBinContent(numOfUsedTS+1,h_mFstSimpleClustersSignal_Apv[numOfUsedTS][i_rstrip]->GetMean());
	}
	leg_FST->Draw("same");
      }

      c_Signal->cd(5*i_ts+3);
      {
	TLegend *leg_FST = new TLegend(0.4,0.6,0.85,0.8);
	leg_FST->SetBorderSize(0);
	leg_FST->SetFillColor(10);
	for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
	{
	  string title = Form("FST Noise @ Date %s",date[numOfUsedTS].c_str());
	  h_mNoiseHits_Apv[numOfUsedTS][i_rstrip]->SetStats(0);
	  h_mNoiseHits_Apv[numOfUsedTS][i_rstrip]->SetTitle(title.c_str());
	  h_mNoiseHits_Apv[numOfUsedTS][i_rstrip]->GetXaxis()->SetTitle("ADC");
	  h_mNoiseHits_Apv[numOfUsedTS][i_rstrip]->GetXaxis()->SetTitleSize(0.06);
	  h_mNoiseHits_Apv[numOfUsedTS][i_rstrip]->SetLineColor(i_rstrip+1);

	  if(i_rstrip == 0) h_mNoiseHits_Apv[numOfUsedTS][i_rstrip]->Draw();
	  else h_mNoiseHits_Apv[numOfUsedTS][i_rstrip]->Draw("same");

	  string LegName = Form("R_strip %d: Mean %1.2f",i_rstrip,h_mNoiseHits_Apv[numOfUsedTS][i_rstrip]->GetMean());
	  leg_FST->AddEntry(h_mNoiseHits_Apv[numOfUsedTS][i_rstrip],LegName.c_str(),"L");

	  double mean_orig    = h_mNoiseHits_Apv[0][i_rstrip]->GetMean();
	  double err_orig     = h_mNoiseHits_Apv[0][i_rstrip]->GetMeanError();
	  double mean_current = h_mNoiseHits_Apv[numOfUsedTS][i_rstrip]->GetMean();
	  double err_current  = h_mNoiseHits_Apv[numOfUsedTS][i_rstrip]->GetMeanError();
	  double ratio = mean_current/mean_orig;
	  double error = ErrDiv(mean_current, mean_orig, err_current, err_orig);
	  h_mMeanNoiseHits_Apv[i_rstrip]->SetBinContent(numOfUsedTS+1,ratio);
	  h_mMeanNoiseHits_Apv[i_rstrip]->SetBinError(numOfUsedTS+1,error);
	  // h_mMeanNoiseHits_Apv[i_rstrip]->SetBinContent(numOfUsedTS+1,h_mNoiseHits_Apv[numOfUsedTS][i_rstrip]->GetMean());
	}
	leg_FST->Draw("same");
      }

      c_Signal->cd(5*i_ts+4);
      {
	TLegend *leg_FST = new TLegend(0.50,0.6,0.90,0.8);
	leg_FST->SetBorderSize(0);
	leg_FST->SetFillColor(10);
	for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
	{
	  string title = Form("FST Signal/Noise @ Date %s",date[numOfUsedTS].c_str());
	  h_mSNRatioHits_Apv[numOfUsedTS][i_rstrip]->SetStats(0);
	  h_mSNRatioHits_Apv[numOfUsedTS][i_rstrip]->SetTitle(title.c_str());
	  h_mSNRatioHits_Apv[numOfUsedTS][i_rstrip]->GetXaxis()->SetTitle("Signal/Noise");
	  h_mSNRatioHits_Apv[numOfUsedTS][i_rstrip]->GetXaxis()->SetTitleSize(0.06);
	  h_mSNRatioHits_Apv[numOfUsedTS][i_rstrip]->SetLineColor(i_rstrip+1);

	  if(i_rstrip == 0) h_mSNRatioHits_Apv[numOfUsedTS][i_rstrip]->Draw();
	  else h_mSNRatioHits_Apv[numOfUsedTS][i_rstrip]->Draw("same");

	  string LegName = Form("R_strip %d: Mean %1.2f",i_rstrip,h_mSNRatioHits_Apv[numOfUsedTS][i_rstrip]->GetMean());
	  leg_FST->AddEntry(h_mSNRatioHits_Apv[numOfUsedTS][i_rstrip],LegName.c_str(),"L");

	  double mean_orig    = h_mSNRatioHits_Apv[0][i_rstrip]->GetMean();
	  double err_orig     = h_mSNRatioHits_Apv[0][i_rstrip]->GetMeanError();
	  double mean_current = h_mSNRatioHits_Apv[numOfUsedTS][i_rstrip]->GetMean();
	  double err_current  = h_mSNRatioHits_Apv[numOfUsedTS][i_rstrip]->GetMeanError();
	  double ratio = mean_current/mean_orig;
	  double error = ErrDiv(mean_current, mean_orig, err_current, err_orig);
	  h_mMeanSNRatioHits_Apv[i_rstrip]->SetBinContent(numOfUsedTS+1,ratio);
	  h_mMeanSNRatioHits_Apv[i_rstrip]->SetBinError(numOfUsedTS+1,error);
	  // h_mMeanSNRatioHits_Apv[i_rstrip]->SetBinContent(numOfUsedTS+1,h_mSNRatioHits_Apv[numOfUsedTS][i_rstrip]->GetMean());
	}
	leg_FST->Draw("same");
      }

      c_Signal->cd(5*i_ts+5);
      {
	TLegend *leg_FST = new TLegend(0.55,0.6,0.8,0.8);
	leg_FST->SetBorderSize(0);
	leg_FST->SetFillColor(10);
	for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
	{
	  h_mMaxTbHits_Apv[numOfUsedTS][i_rstrip]->SetStats(0);
	  h_mMaxTbHits_Apv[numOfUsedTS][i_rstrip]->SetTitle("Max Time Bin Hits");
	  h_mMaxTbHits_Apv[numOfUsedTS][i_rstrip]->GetXaxis()->SetTitle("Maxt Time Bin");
	  h_mMaxTbHits_Apv[numOfUsedTS][i_rstrip]->GetXaxis()->SetTitleSize(0.06);
	  h_mMaxTbHits_Apv[numOfUsedTS][i_rstrip]->SetLineColor(i_rstrip+1);

	  if(i_rstrip == 0) h_mMaxTbHits_Apv[numOfUsedTS][i_rstrip]->Draw();
	  else h_mMaxTbHits_Apv[numOfUsedTS][i_rstrip]->Draw("same");

	  string LegName = Form("R_strip %d: Mean %1.2f",i_rstrip,h_mMaxTbHits_Apv[numOfUsedTS][i_rstrip]->GetMean());
	  leg_FST->AddEntry(h_mMaxTbHits_Apv[numOfUsedTS][i_rstrip],LegName.c_str(),"L");

	  double mean_orig    = h_mMaxTbHits_Apv[0][i_rstrip]->GetMean();
	  double err_orig     = h_mMaxTbHits_Apv[0][i_rstrip]->GetMeanError();
	  double mean_current = h_mMaxTbHits_Apv[numOfUsedTS][i_rstrip]->GetMean();
	  double err_current  = h_mMaxTbHits_Apv[numOfUsedTS][i_rstrip]->GetMeanError();
	  double ratio = mean_current/mean_orig;
	  double error = ErrDiv(mean_current, mean_orig, err_current, err_orig);
	  h_mMeanMaxTbHits_Apv[i_rstrip]->SetBinContent(numOfUsedTS+1,ratio);
	  h_mMeanMaxTbHits_Apv[i_rstrip]->SetBinError(numOfUsedTS+1,error);
	  // h_mMeanMaxTbHits_Apv[i_rstrip]->SetBinContent(numOfUsedTS+1,h_mMaxTbHits_Apv[numOfUsedTS][i_rstrip]->GetMean());
	}
	leg_FST->Draw("same");
      }
      numOfUsedTS++;
    }
    c_Signal->Update();
    c_Signal->Print(outputname.c_str());
  }


  TCanvas *c_Signal_TS = new TCanvas("c_Signal_TS","c_Signal_TS",10,10,1000,1000);
  c_Signal_TS->Divide(1,5);
  for(int i_pad = 0;i_pad < 5; ++i_pad)
  {
    c_Signal_TS->cd(i_pad+1)->SetLeftMargin(0.10);
    c_Signal_TS->cd(i_pad+1)->SetRightMargin(0.0);
    c_Signal_TS->cd(i_pad+1)->SetTopMargin(0.05);
    c_Signal_TS->cd(i_pad+1)->SetTicks(1,1);
    c_Signal_TS->cd(i_pad+1)->SetGrid(0,0);
  }

  c_Signal_TS->cd(1);
  {
    TLegend *leg_FST = new TLegend(0.7,0.6,0.8,0.9);
    leg_FST->SetBorderSize(0);
    leg_FST->SetFillColor(10);
    for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
    {
      h_mMeanSignalHits_Apv[i_rstrip]->SetTitle("");
      h_mMeanSignalHits_Apv[i_rstrip]->SetStats(0);
      for(int i_ts = 0; i_ts < numOfTS; ++i_ts)
      {
	h_mMeanSignalHits_Apv[i_rstrip]->GetXaxis()->SetBinLabel(i_ts+1,date[i_ts].c_str());
      }
      h_mMeanSignalHits_Apv[i_rstrip]->GetXaxis()->SetLabelSize(0.14);
      h_mMeanSignalHits_Apv[i_rstrip]->GetYaxis()->SetTitle("#frac{<ADC>_{Hits}}{<ADC>_{Hits @ 0320}}");
      h_mMeanSignalHits_Apv[i_rstrip]->GetYaxis()->SetTitleSize(0.10);
      h_mMeanSignalHits_Apv[i_rstrip]->GetYaxis()->SetTitleOffset(0.5);
      // h_mMeanSignalHits_Apv[i_rstrip]->GetYaxis()->SetRangeUser(200.0,650.0);
      h_mMeanSignalHits_Apv[i_rstrip]->GetYaxis()->SetRangeUser(0.3,1.5);
      h_mMeanSignalHits_Apv[i_rstrip]->GetYaxis()->SetLabelSize(0.08);
      h_mMeanSignalHits_Apv[i_rstrip]->GetYaxis()->SetNdivisions(505);
      h_mMeanSignalHits_Apv[i_rstrip]->SetLineColor(i_rstrip+1);

      if(i_rstrip == 0) h_mMeanSignalHits_Apv[i_rstrip]->Draw();
      else h_mMeanSignalHits_Apv[i_rstrip]->Draw("same");

      string LegName = Form("R_strip %d",i_rstrip);
      leg_FST->AddEntry(h_mMeanSignalHits_Apv[i_rstrip],LegName.c_str(),"L");
    }
    leg_FST->Draw("same");
    PlotLine(-0.5, numOfTS-0.5, 1.0, 1.0, 1, 2, 2);
  }

  c_Signal_TS->cd(2);
  {
    for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
    {
      h_mMeanSignalClusters_Apv[i_rstrip]->SetTitle("");
      h_mMeanSignalClusters_Apv[i_rstrip]->SetStats(0);
      for(int i_ts = 0; i_ts < numOfTS; ++i_ts)
      {
	h_mMeanSignalClusters_Apv[i_rstrip]->GetXaxis()->SetBinLabel(i_ts+1,date[i_ts].c_str());
      }
      h_mMeanSignalClusters_Apv[i_rstrip]->GetXaxis()->SetLabelSize(0.14);
      h_mMeanSignalClusters_Apv[i_rstrip]->GetYaxis()->SetTitle("#frac{<ADC>_{Clusters}}{<ADC>_{Clusters @ 0320}}");
      h_mMeanSignalClusters_Apv[i_rstrip]->GetYaxis()->SetTitleSize(0.10);
      h_mMeanSignalClusters_Apv[i_rstrip]->GetYaxis()->SetTitleOffset(0.5);
      // h_mMeanSignalClusters_Apv[i_rstrip]->GetYaxis()->SetRangeUser(200.0,650.0);
      h_mMeanSignalClusters_Apv[i_rstrip]->GetYaxis()->SetRangeUser(0.3,1.5);
      h_mMeanSignalClusters_Apv[i_rstrip]->GetYaxis()->SetLabelSize(0.08);
      h_mMeanSignalClusters_Apv[i_rstrip]->GetYaxis()->SetNdivisions(505);
      h_mMeanSignalClusters_Apv[i_rstrip]->SetLineColor(i_rstrip+1);

      if(i_rstrip == 0) h_mMeanSignalClusters_Apv[i_rstrip]->Draw();
      else h_mMeanSignalClusters_Apv[i_rstrip]->Draw("same");
    }
    PlotLine(-0.5, numOfTS-0.5, 1.0, 1.0, 1, 2, 2);
  }

  c_Signal_TS->cd(3);
  {
    for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
    {
      h_mMeanNoiseHits_Apv[i_rstrip]->SetTitle("");
      h_mMeanNoiseHits_Apv[i_rstrip]->SetStats(0);
      for(int i_ts = 0; i_ts < numOfTS; ++i_ts)
      {
	h_mMeanNoiseHits_Apv[i_rstrip]->GetXaxis()->SetBinLabel(i_ts+1,date[i_ts].c_str());
      }
      h_mMeanNoiseHits_Apv[i_rstrip]->GetXaxis()->SetLabelSize(0.14);
      h_mMeanNoiseHits_Apv[i_rstrip]->GetYaxis()->SetTitle("Noise #frac{<ADC>_{Hits}}{<ADC>_{Hits @ 0320}}");
      h_mMeanNoiseHits_Apv[i_rstrip]->GetYaxis()->SetTitleSize(0.10);
      h_mMeanNoiseHits_Apv[i_rstrip]->GetYaxis()->SetTitleOffset(0.5);
      // h_mMeanNoiseHits_Apv[i_rstrip]->GetYaxis()->SetRangeUser(10.0,35.0);
      h_mMeanNoiseHits_Apv[i_rstrip]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_mMeanNoiseHits_Apv[i_rstrip]->GetYaxis()->SetLabelSize(0.08);
      h_mMeanNoiseHits_Apv[i_rstrip]->GetYaxis()->SetNdivisions(505);
      h_mMeanNoiseHits_Apv[i_rstrip]->SetLineColor(i_rstrip+1);

      if(i_rstrip == 0) h_mMeanNoiseHits_Apv[i_rstrip]->Draw();
      else h_mMeanNoiseHits_Apv[i_rstrip]->Draw("same");
    }
    PlotLine(-0.5, numOfTS-0.5, 1.0, 1.0, 1, 2, 2);
  }

  c_Signal_TS->cd(4);
  {
    for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
    {
      h_mMeanSNRatioHits_Apv[i_rstrip]->SetTitle("");
      h_mMeanSNRatioHits_Apv[i_rstrip]->SetStats(0);
      for(int i_ts = 0; i_ts < numOfTS; ++i_ts)
      {
	h_mMeanSNRatioHits_Apv[i_rstrip]->GetXaxis()->SetBinLabel(i_ts+1,date[i_ts].c_str());
      }
      h_mMeanSNRatioHits_Apv[i_rstrip]->GetXaxis()->SetLabelSize(0.14);
      h_mMeanSNRatioHits_Apv[i_rstrip]->GetYaxis()->SetTitle("#frac{<Signal/Noise>_{Hits}}{<Signal/Noise>_{Hits @ 0320}}");
      h_mMeanSNRatioHits_Apv[i_rstrip]->GetYaxis()->SetTitleSize(0.10);
      h_mMeanSNRatioHits_Apv[i_rstrip]->GetYaxis()->SetTitleOffset(0.5);
      // h_mMeanSNRatioHits_Apv[i_rstrip]->GetYaxis()->SetRangeUser(5.0,30.0);
      h_mMeanSNRatioHits_Apv[i_rstrip]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_mMeanSNRatioHits_Apv[i_rstrip]->GetYaxis()->SetLabelSize(0.08);
      h_mMeanSNRatioHits_Apv[i_rstrip]->GetYaxis()->SetNdivisions(505);
      h_mMeanSNRatioHits_Apv[i_rstrip]->SetLineColor(i_rstrip+1);

      if(i_rstrip == 0) h_mMeanSNRatioHits_Apv[i_rstrip]->Draw();
      else h_mMeanSNRatioHits_Apv[i_rstrip]->Draw("same");
    }
    PlotLine(-0.5, numOfTS-0.5, 1.0, 1.0, 1, 2, 2);
  }

  c_Signal_TS->cd(5);
  {
    for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
    {
      h_mMeanMaxTbHits_Apv[i_rstrip]->SetTitle("");
      h_mMeanMaxTbHits_Apv[i_rstrip]->SetStats(0);
      for(int i_ts = 0; i_ts < numOfTS; ++i_ts)
      {
	h_mMeanMaxTbHits_Apv[i_rstrip]->GetXaxis()->SetBinLabel(i_ts+1,date[i_ts].c_str());
      }
      h_mMeanMaxTbHits_Apv[i_rstrip]->GetXaxis()->SetLabelSize(0.14);
      h_mMeanMaxTbHits_Apv[i_rstrip]->GetYaxis()->SetTitle("#frac{<Max Time Bin>_{Hits}}{<Max Time Bin>_{Hits @ 0320}}");
      h_mMeanMaxTbHits_Apv[i_rstrip]->GetYaxis()->SetTitleSize(0.10);
      h_mMeanMaxTbHits_Apv[i_rstrip]->GetYaxis()->SetTitleOffset(0.5);
      // h_mMeanMaxTbHits_Apv[i_rstrip]->GetYaxis()->SetRangeUser(2.0,5.0);
      h_mMeanMaxTbHits_Apv[i_rstrip]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_mMeanMaxTbHits_Apv[i_rstrip]->GetYaxis()->SetLabelSize(0.08);
      h_mMeanMaxTbHits_Apv[i_rstrip]->GetYaxis()->SetNdivisions(505);
      h_mMeanMaxTbHits_Apv[i_rstrip]->SetLineColor(i_rstrip+1);

      if(i_rstrip == 0) h_mMeanMaxTbHits_Apv[i_rstrip]->Draw();
      else h_mMeanMaxTbHits_Apv[i_rstrip]->Draw("same");
    }
    PlotLine(-0.5, numOfTS-0.5, 1.0, 1.0, 1, 2, 2);
  }

  c_Signal_TS->Update();
  c_Signal_TS->Print(outputname.c_str());

  string output_stop = Form("./figures/SignalTS_70V_Apv%d_Phi%d.pdf]",apv+5,phibin);
  c_Signal->Print(output_stop.c_str()); // open pdf file
}

