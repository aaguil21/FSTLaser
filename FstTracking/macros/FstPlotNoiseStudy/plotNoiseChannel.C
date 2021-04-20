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

void plotNoiseChannel(string mod = "Mod03", string hv = "HV70V")
{
  const int defTimeBin = 0;

  string input_Sensor = Form("../../output/noise/%s/HVScan/FstPedNoise_%s_%s.root",mod.c_str(),mod.c_str(),hv.c_str());
  TFile *File_Sensor = TFile::Open(input_Sensor.c_str());
  TH1F *h_mPedSigma_Sensor[FST::numRStrip][FST::numTBins];
  TH1F *h_mRanSigma_Sensor[FST::numRStrip][FST::numTBins];
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      std::string HistName;
      HistName = Form("h_mPedSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mPedSigma_Sensor[i_rstrip][i_tb] = (TH1F*)File_Sensor->Get(HistName.c_str());
      HistName = Form("h_mRanSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mRanSigma_Sensor[i_rstrip][i_tb] = (TH1F*)File_Sensor->Get(HistName.c_str());
    }
  }

  string input_Sensor0V = Form("../../output/noise/%s/HVScan/FstPedNoise_%s_HV0V.root",mod.c_str(),mod.c_str());
  TFile *File_Sensor0V = TFile::Open(input_Sensor0V.c_str());
  TH1F *h_mPedSigma_Sensor0V[FST::numRStrip][FST::numTBins];
  TH1F *h_mRanSigma_Sensor0V[FST::numRStrip][FST::numTBins];
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      std::string HistName;
      HistName = Form("h_mPedSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mPedSigma_Sensor0V[i_rstrip][i_tb] = (TH1F*)File_Sensor0V->Get(HistName.c_str());
      HistName = Form("h_mRanSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mRanSigma_Sensor0V[i_rstrip][i_tb] = (TH1F*)File_Sensor0V->Get(HistName.c_str());
    }
  }
  
  string input_InnerChip = Form("../../output/noise/%s/FstChipNoise_%s_Inner.root",mod.c_str(),mod.c_str());
  TFile *File_InnerChip = TFile::Open(input_InnerChip.c_str());
  TH1F *h_mPedSigma_InnerChip[FST::numRStrip][FST::numTBins];
  TH1F *h_mRanSigma_InnerChip[FST::numRStrip][FST::numTBins];
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      std::string HistName;
      HistName = Form("h_mPedSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mPedSigma_InnerChip[i_rstrip][i_tb] = (TH1F*)File_InnerChip->Get(HistName.c_str());
      HistName = Form("h_mRanSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mRanSigma_InnerChip[i_rstrip][i_tb] = (TH1F*)File_InnerChip->Get(HistName.c_str());
    }
  }
  
  string input_OuterChip = Form("../../output/noise/%s/FstChipNoise_%s_Outer.root",mod.c_str(),mod.c_str());
  TFile *File_OuterChip = TFile::Open(input_OuterChip.c_str());
  TH1F *h_mPedSigma_OuterChip[FST::numRStrip][FST::numTBins];
  TH1F *h_mRanSigma_OuterChip[FST::numRStrip][FST::numTBins];
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      std::string HistName;
      HistName = Form("h_mPedSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mPedSigma_OuterChip[i_rstrip][i_tb] = (TH1F*)File_OuterChip->Get(HistName.c_str());
      HistName = Form("h_mRanSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mRanSigma_OuterChip[i_rstrip][i_tb] = (TH1F*)File_OuterChip->Get(HistName.c_str());
    }
  }

  TCanvas *c_PedSigma = new TCanvas("c_PedSigma","c_PedSigma",10,10,1600,800);
  c_PedSigma->Divide(4,2);
  for(int i_pad = 0; i_pad < 8; ++i_pad)
  {
    c_PedSigma->cd(i_pad+1)->SetLeftMargin(0.15);
    c_PedSigma->cd(i_pad+1)->SetBottomMargin(0.15);
    c_PedSigma->cd(i_pad+1)->SetTicks(1,1);
    c_PedSigma->cd(i_pad+1)->SetGrid(0,0);
  }

  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    c_PedSigma->cd(i_rstrip+1);
    string title = Form("Total Noise: R-Strip %d",i_rstrip);
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->SetTitle(title.c_str());
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->SetStats(0);
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->GetXaxis()->SetTitle("Channel");
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->GetXaxis()->SetTitleSize(0.06);
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->GetXaxis()->SetLabelSize(0.06);
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->GetYaxis()->SetRangeUser(-1.5,90.0);
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->GetYaxis()->SetTitle("ADC");
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->GetYaxis()->SetTitleSize(0.06);
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->GetYaxis()->SetLabelSize(0.06);
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->SetMarkerStyle(24);
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->SetMarkerSize(1.0);
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->SetMarkerColor(1);
    h_mPedSigma_Sensor[i_rstrip][defTimeBin]->Draw("p");
    h_mPedSigma_Sensor0V[i_rstrip][defTimeBin]->SetMarkerStyle(24);
    h_mPedSigma_Sensor0V[i_rstrip][defTimeBin]->SetMarkerSize(1.0);
    h_mPedSigma_Sensor0V[i_rstrip][defTimeBin]->SetMarkerColor(4);
    h_mPedSigma_Sensor0V[i_rstrip][defTimeBin]->Draw("p same");
    h_mPedSigma_InnerChip[i_rstrip][defTimeBin]->SetMarkerStyle(24);
    h_mPedSigma_InnerChip[i_rstrip][defTimeBin]->SetMarkerSize(1.0);
    h_mPedSigma_InnerChip[i_rstrip][defTimeBin]->SetMarkerColor(2);
    if(i_rstrip < 4) 
    {
      h_mPedSigma_InnerChip[i_rstrip][defTimeBin]->DrawCopy("p same");
    }
    else
    {
      h_mPedSigma_InnerChip[i_rstrip][defTimeBin]->GetXaxis()->SetRangeUser(-0.5,31.5);
      h_mPedSigma_InnerChip[i_rstrip][defTimeBin]->DrawCopy("p same");
      h_mPedSigma_InnerChip[i_rstrip][defTimeBin]->GetXaxis()->SetRangeUser(63.5,127.5);
      h_mPedSigma_InnerChip[i_rstrip][defTimeBin]->DrawCopy("p same");
      h_mPedSigma_OuterChip[i_rstrip][defTimeBin]->SetMarkerStyle(24);
      h_mPedSigma_OuterChip[i_rstrip][defTimeBin]->SetMarkerSize(1.0);
      h_mPedSigma_OuterChip[i_rstrip][defTimeBin]->SetMarkerColor(2);
      h_mPedSigma_OuterChip[i_rstrip][defTimeBin]->GetXaxis()->SetRangeUser(31.5,63.5);
      h_mPedSigma_OuterChip[i_rstrip][defTimeBin]->DrawCopy("p same");
    }
    PlotLine(-0.5, 127.5,0.0,0.0,1,1,2);
    string leg_HV = Form("With Sensors (%s)",hv.c_str());
    TLegend *leg = new TLegend(0.15,0.65,0.85,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry(h_mPedSigma_InnerChip[i_rstrip][defTimeBin],"W/O Sensors","P");
    leg->AddEntry(h_mPedSigma_Sensor0V[i_rstrip][defTimeBin],"With Sensors (HV0V)","P");
    leg->AddEntry(h_mPedSigma_Sensor[i_rstrip][defTimeBin],leg_HV.c_str(),"P");
    leg->Draw("same");
  }
  string FigName_Ped = Form("./figures/%s/c_PedSigma_%s_%s.eps",mod.c_str(),mod.c_str(),hv.c_str());
  c_PedSigma->SaveAs(FigName_Ped.c_str());

  TCanvas *c_RanSigma = new TCanvas("c_RanSigma","c_RanSigma",10,10,1600,800);
  c_RanSigma->Divide(4,2);
  for(int i_pad = 0; i_pad < 8; ++i_pad)
  {
    c_RanSigma->cd(i_pad+1)->SetLeftMargin(0.15);
    c_RanSigma->cd(i_pad+1)->SetBottomMargin(0.15);
    c_RanSigma->cd(i_pad+1)->SetTicks(1,1);
    c_RanSigma->cd(i_pad+1)->SetGrid(0,0);
  }

  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    c_RanSigma->cd(i_rstrip+1);
    string title = Form("Random Noise: R-Strip %d",i_rstrip);
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->SetTitle(title.c_str());
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->SetStats(0);
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->GetXaxis()->SetTitle("Channel");
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->GetXaxis()->SetTitleSize(0.06);
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->GetXaxis()->SetLabelSize(0.06);
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->GetYaxis()->SetRangeUser(-1.5,90.0);
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->GetYaxis()->SetTitle("ADC");
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->GetYaxis()->SetTitleSize(0.06);
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->GetYaxis()->SetLabelSize(0.06);
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->SetMarkerStyle(24);
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->SetMarkerSize(1.0);
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->SetMarkerColor(1);
    h_mRanSigma_Sensor[i_rstrip][defTimeBin]->Draw("p");
    h_mRanSigma_Sensor0V[i_rstrip][defTimeBin]->SetMarkerStyle(24);
    h_mRanSigma_Sensor0V[i_rstrip][defTimeBin]->SetMarkerSize(1.0);
    h_mRanSigma_Sensor0V[i_rstrip][defTimeBin]->SetMarkerColor(4);
    h_mRanSigma_Sensor0V[i_rstrip][defTimeBin]->Draw("p same");
    h_mRanSigma_InnerChip[i_rstrip][defTimeBin]->SetMarkerStyle(24);
    h_mRanSigma_InnerChip[i_rstrip][defTimeBin]->SetMarkerSize(1.0);
    h_mRanSigma_InnerChip[i_rstrip][defTimeBin]->SetMarkerColor(2);
    if(i_rstrip < 4) 
    {
      h_mRanSigma_InnerChip[i_rstrip][defTimeBin]->DrawCopy("p same");
    }
    else
    {
      h_mRanSigma_InnerChip[i_rstrip][defTimeBin]->GetXaxis()->SetRangeUser(-0.5,31.5);
      h_mRanSigma_InnerChip[i_rstrip][defTimeBin]->DrawCopy("p same");
      h_mRanSigma_InnerChip[i_rstrip][defTimeBin]->GetXaxis()->SetRangeUser(63.5,127.5);
      h_mRanSigma_InnerChip[i_rstrip][defTimeBin]->DrawCopy("p same");
      h_mRanSigma_OuterChip[i_rstrip][defTimeBin]->SetMarkerStyle(24);
      h_mRanSigma_OuterChip[i_rstrip][defTimeBin]->SetMarkerSize(1.0);
      h_mRanSigma_OuterChip[i_rstrip][defTimeBin]->SetMarkerColor(2);
      h_mRanSigma_OuterChip[i_rstrip][defTimeBin]->GetXaxis()->SetRangeUser(31.5,63.5);
      h_mRanSigma_OuterChip[i_rstrip][defTimeBin]->DrawCopy("p same");
    }
    PlotLine(-0.5, 127.5,0.0,0.0,1,1,2);
    string leg_HV = Form("With Seonsrs (%s)",hv.c_str());
    TLegend *leg = new TLegend(0.15,0.65,0.85,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetMargin(0.1);
    leg->AddEntry(h_mRanSigma_InnerChip[i_rstrip][defTimeBin],"W/O Sensors","P");
    leg->AddEntry(h_mRanSigma_Sensor0V[i_rstrip][defTimeBin],"With Sensors (HV0V)","P");
    leg->AddEntry(h_mRanSigma_Sensor[i_rstrip][defTimeBin],leg_HV.c_str(),"P");
    leg->Draw("same");
  }
  string FigName_Ran = Form("./figures/%s/c_RanSigma_%s_%s.eps",mod.c_str(),mod.c_str(),hv.c_str());
  c_RanSigma->SaveAs(FigName_Ran.c_str());
}
