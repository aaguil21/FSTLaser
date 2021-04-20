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

void plotIstChipNoiseQA(string module = "Mod36")
{
  gStyle->SetOptStat(111111);
  gStyle->SetStatX(0.95); gStyle->SetStatY(0.95);
  gStyle->SetStatW(0.35); gStyle->SetStatH(0.35);

  string inputfile = Form("../../output/noise/Ist%s/IstChipNoise_%s.root",module.c_str(),module.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH2F *h_mPedDisplay[4][FST::numTBins];
  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      std::string HistName = Form("h_mPedDisplay_Layer%d_TimeBin%d",i_layer,i_tb);
      h_mPedDisplay[i_layer][i_tb] = (TH2F*)File_InPut->Get(HistName.c_str());
    }
  }

  TH1F *h_mPedMean_FST[FST::numRStrip][FST::numTBins]; // for RStrip
  TH1F *h_mPedSigma_FST[FST::numRStrip][FST::numTBins];
  TH1F *h_mCMNSigma_FST[FST::numRStrip][FST::numTBins];
  TH1F *h_mRanSigma_FST[FST::numRStrip][FST::numTBins];
  TH1F *h_mDisplayPedSigma_FST[FST::numRStrip][FST::numTBins];
  TH1F *h_mDisplayCMNSigma_FST[FST::numRStrip][FST::numTBins];
  TH1F *h_mDisplayRanSigma_FST[FST::numRStrip][FST::numTBins];
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      string HistName;
      HistName = Form("h_mPedMean_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mPedMean_FST[i_rstrip][i_tb] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mPedSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mPedSigma_FST[i_rstrip][i_tb] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mCMNSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mCMNSigma_FST[i_rstrip][i_tb] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mRanSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mRanSigma_FST[i_rstrip][i_tb] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mDisplayPedSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mDisplayPedSigma_FST[i_rstrip][i_tb] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mDisplayCMNSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mDisplayCMNSigma_FST[i_rstrip][i_tb] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mDisplayRanSigma_FST_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_mDisplayRanSigma_FST[i_rstrip][i_tb] = (TH1F*)File_InPut->Get(HistName.c_str());
    }
  }

  string outputname = Form("./figures/Ist%s/IstChipNoiseQA_%s.pdf",module.c_str(),module.c_str());

  TCanvas *c_Noise = new TCanvas("c_Noise","c_Noise",10,10,1800,800);
  c_Noise->Divide(9,4);
  for(int i_pad = 0; i_pad < 36; ++i_pad)
  {
    c_Noise->cd(i_pad+1)->SetLeftMargin(0.1);
    // c_Noise->cd(i_pad+1)->SetRightMargin(0.15);
    c_Noise->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Noise->cd(i_pad+1)->SetTicks(1,1);
    c_Noise->cd(i_pad+1)->SetGrid(0,0);
  }

  string output_start = Form("./figures/Ist%s/IstChipNoiseQA_%s.pdf[",module.c_str(),module.c_str());
  c_Noise->Print(output_start.c_str()); // open pdf file

  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      int i_pad = i_layer*9 + i_tb;
      c_Noise->cd(i_pad+1);
      string title = Form("Layer%d & TB%d",i_layer,i_tb);
      h_mPedDisplay[i_layer][i_tb]->SetTitle(title.c_str());
      h_mPedDisplay[i_layer][i_tb]->SetStats(0);
      h_mPedDisplay[i_layer][i_tb]->Draw("colz");
    }
  }
  c_Noise->Update();
  c_Noise->Print(outputname.c_str());

  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      int i_pad = i_rstrip*9 + i_tb;
      if(i_rstrip > 3) i_pad = (i_rstrip-4)*9 + i_tb;
      c_Noise->cd(i_pad+1);
      string title = Form("RStrip%d & TB%d",i_rstrip,i_tb);
      h_mPedMean_FST[i_rstrip][i_tb]->SetTitle(title.c_str());
      h_mPedMean_FST[i_rstrip][i_tb]->SetStats(0);
      h_mPedMean_FST[i_rstrip][i_tb]->GetXaxis()->SetRangeUser(-0.5,127.5);
      h_mPedMean_FST[i_rstrip][i_tb]->GetXaxis()->SetLabelSize(0.06);
      h_mPedMean_FST[i_rstrip][i_tb]->GetXaxis()->SetTitle("channel");
      h_mPedMean_FST[i_rstrip][i_tb]->GetYaxis()->SetRangeUser(200,1000);
      h_mPedMean_FST[i_rstrip][i_tb]->GetYaxis()->SetLabelSize(0.06);
      h_mPedMean_FST[i_rstrip][i_tb]->GetYaxis()->SetTitle("ADC");
      h_mPedMean_FST[i_rstrip][i_tb]->Draw();
      PlotLine(31.5, 31.5,200,1000,2,1,2);
      PlotLine(63.5, 63.5,200,1000,2,1,2);
      PlotLine(95.5, 95.5,200,1000,2,1,2);
    }
    if(i_rstrip == 3 || i_rstrip == FST::numRStrip-1)
    {
      c_Noise->cd();
      string TitleName = "Pedestal";
      TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.055,1.0,TitleName.c_str());
      canvas_title->SetBorderSize(1);
      canvas_title->Draw();
      c_Noise->Update();
      c_Noise->Print(outputname.c_str());
    }
  }

  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      int i_pad = i_rstrip*9 + i_tb;
      if(i_rstrip > 3) i_pad = (i_rstrip-4)*9 + i_tb;
      c_Noise->cd(i_pad+1);
      string title = Form("RStrip%d & TB%d",i_rstrip,i_tb);
      h_mPedSigma_FST[i_rstrip][i_tb]->SetTitle(title.c_str());
      h_mPedSigma_FST[i_rstrip][i_tb]->SetStats(0);
      h_mPedSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetRangeUser(-0.5,127.5);
      h_mPedSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetLabelSize(0.06);
      h_mPedSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetTitle("channel");
      h_mPedSigma_FST[i_rstrip][i_tb]->GetYaxis()->SetRangeUser(0,80);
      h_mPedSigma_FST[i_rstrip][i_tb]->GetYaxis()->SetLabelSize(0.06);
      h_mPedSigma_FST[i_rstrip][i_tb]->GetYaxis()->SetTitle("ADC");
      h_mPedSigma_FST[i_rstrip][i_tb]->Draw();
      PlotLine(31.5, 31.5,0,80,2,1,2);
      PlotLine(63.5, 63.5,0,80,2,1,2);
      PlotLine(95.5, 95.5,0,80,2,1,2);
    }
    if(i_rstrip == 3 || i_rstrip == FST::numRStrip-1)
    {
      c_Noise->cd();
      string TitleName = "Total Noise";
      TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.055,1.0,TitleName.c_str());
      canvas_title->SetBorderSize(1);
      canvas_title->Draw();
      c_Noise->Update();
      c_Noise->Print(outputname.c_str());
    }
  }

  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      int i_pad = i_rstrip*9 + i_tb;
      if(i_rstrip > 3) i_pad = (i_rstrip-4)*9 + i_tb;
      c_Noise->cd(i_pad+1);
      string title = Form("RStrip%d & TB%d",i_rstrip,i_tb);
      h_mCMNSigma_FST[i_rstrip][i_tb]->SetTitle(title.c_str());
      h_mCMNSigma_FST[i_rstrip][i_tb]->SetStats(0);
      h_mCMNSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetRangeUser(-0.5,127.5);
      h_mCMNSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetLabelSize(0.06);
      h_mCMNSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetTitle("channel");
      h_mCMNSigma_FST[i_rstrip][i_tb]->GetYaxis()->SetRangeUser(0,80);
      h_mCMNSigma_FST[i_rstrip][i_tb]->GetYaxis()->SetLabelSize(0.06);
      h_mCMNSigma_FST[i_rstrip][i_tb]->GetYaxis()->SetTitle("ADC");
      h_mCMNSigma_FST[i_rstrip][i_tb]->Draw();
      PlotLine(31.5, 31.5,0,80,2,1,2);
      PlotLine(63.5, 63.5,0,80,2,1,2);
      PlotLine(95.5, 95.5,0,80,2,1,2);
    }
    if(i_rstrip == 3 || i_rstrip == FST::numRStrip-1)
    {
      c_Noise->cd();
      string TitleName = "Common Mode Noise";
      TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.055,1.0,TitleName.c_str());
      canvas_title->SetBorderSize(1);
      canvas_title->Draw();
      c_Noise->Update();
      c_Noise->Print(outputname.c_str());
    }
  }

  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      int i_pad = i_rstrip*9 + i_tb;
      if(i_rstrip > 3) i_pad = (i_rstrip-4)*9 + i_tb;
      c_Noise->cd(i_pad+1);
      string title = Form("RStrip%d & TB%d",i_rstrip,i_tb);
      h_mRanSigma_FST[i_rstrip][i_tb]->SetTitle(title.c_str());
      h_mRanSigma_FST[i_rstrip][i_tb]->SetStats(0);
      h_mRanSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetRangeUser(-0.5,127.5);
      h_mRanSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetLabelSize(0.06);
      h_mRanSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetTitle("channel");
      h_mRanSigma_FST[i_rstrip][i_tb]->GetYaxis()->SetRangeUser(0,80);
      h_mRanSigma_FST[i_rstrip][i_tb]->GetYaxis()->SetLabelSize(0.06);
      h_mRanSigma_FST[i_rstrip][i_tb]->GetYaxis()->SetTitle("ADC");
      h_mRanSigma_FST[i_rstrip][i_tb]->Draw();
      PlotLine(31.5, 31.5,0,80,2,1,2);
      PlotLine(63.5, 63.5,0,80,2,1,2);
      PlotLine(95.5, 95.5,0,80,2,1,2);
    }
    if(i_rstrip == 3 || i_rstrip == FST::numRStrip-1)
    {
      c_Noise->cd();
      string TitleName = "Random Noise";
      TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.055,1.0,TitleName.c_str());
      canvas_title->SetBorderSize(1);
      canvas_title->Draw();
      c_Noise->Update();
      c_Noise->Print(outputname.c_str());
    }
  }

  // noise distribution
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      int i_pad = i_rstrip*9 + i_tb;
      if(i_rstrip > 3) i_pad = (i_rstrip-4)*9 + i_tb;
      c_Noise->cd(i_pad+1);
      string title = Form("RStrip%d & TB%d",i_rstrip,i_tb);
      h_mDisplayPedSigma_FST[i_rstrip][i_tb]->SetTitle(title.c_str());
      // h_mDisplayPedSigma_FST[i_rstrip][i_tb]->SetStats(0);
      h_mDisplayPedSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetLabelSize(0.06);
      h_mDisplayPedSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetTitle("ADC");
      h_mDisplayPedSigma_FST[i_rstrip][i_tb]->GetYaxis()->SetLabelSize(0.06);
      h_mDisplayPedSigma_FST[i_rstrip][i_tb]->Draw();
      h_mDisplayPedSigma_FST[i_rstrip][i_tb]->Fit("gaus");
    }
    if(i_rstrip == 3 || i_rstrip == FST::numRStrip-1)
    {
      c_Noise->cd();
      string TitleName = "Total Noise";
      TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.055,1.0,TitleName.c_str());
      canvas_title->SetBorderSize(1);
      canvas_title->Draw();
      c_Noise->Update();
      c_Noise->Print(outputname.c_str());
    }
  }

  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      int i_pad = i_rstrip*9 + i_tb;
      if(i_rstrip > 3) i_pad = (i_rstrip-4)*9 + i_tb;
      c_Noise->cd(i_pad+1);
      string title = Form("RStrip%d & TB%d",i_rstrip,i_tb);
      h_mDisplayCMNSigma_FST[i_rstrip][i_tb]->SetTitle(title.c_str());
      // h_mDisplayCMNSigma_FST[i_rstrip][i_tb]->SetStats(0);
      h_mDisplayCMNSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetLabelSize(0.06);
      h_mDisplayCMNSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetTitle("ADC");
      h_mDisplayCMNSigma_FST[i_rstrip][i_tb]->GetYaxis()->SetLabelSize(0.06);
      h_mDisplayCMNSigma_FST[i_rstrip][i_tb]->Draw();
      h_mDisplayCMNSigma_FST[i_rstrip][i_tb]->Fit("gaus");
    }
    if(i_rstrip == 3 || i_rstrip == FST::numRStrip-1)
    {
      c_Noise->cd();
      string TitleName = "Common Mode Noise";
      TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.055,1.0,TitleName.c_str());
      canvas_title->SetBorderSize(1);
      canvas_title->Draw();
      c_Noise->Update();
      c_Noise->Print(outputname.c_str());
    }
  }

  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      int i_pad = i_rstrip*9 + i_tb;
      if(i_rstrip > 3) i_pad = (i_rstrip-4)*9 + i_tb;
      c_Noise->cd(i_pad+1);
      string title = Form("RStrip%d & TB%d",i_rstrip,i_tb);
      h_mDisplayRanSigma_FST[i_rstrip][i_tb]->SetTitle(title.c_str());
      // h_mDisplayRanSigma_FST[i_rstrip][i_tb]->SetStats(0);
      h_mDisplayRanSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetLabelSize(0.06);
      h_mDisplayRanSigma_FST[i_rstrip][i_tb]->GetXaxis()->SetTitle("ADC");
      h_mDisplayRanSigma_FST[i_rstrip][i_tb]->GetYaxis()->SetLabelSize(0.06);
      h_mDisplayRanSigma_FST[i_rstrip][i_tb]->Draw();
      h_mDisplayRanSigma_FST[i_rstrip][i_tb]->Fit("gaus");
    }
    if(i_rstrip == 3 || i_rstrip == FST::numRStrip-1)
    {
      c_Noise->cd();
      string TitleName = "Random Noise";
      TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.055,1.0,TitleName.c_str());
      canvas_title->SetBorderSize(1);
      canvas_title->Draw();
      c_Noise->Update();
      c_Noise->Print(outputname.c_str());
    }
  }

  TGraph *g_mPedSigma[4][FST::numTBins];
  TGraph *g_mCMNSigma[4][FST::numTBins];
  TGraph *g_mRanSigma[4][FST::numTBins];
  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      std::string gName;
      gName = Form("g_mPedSigma_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mPedSigma[i_layer][i_tb] = (TGraph*)File_InPut->Get(gName.c_str());

      gName = Form("g_mCMNSigma_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mCMNSigma[i_layer][i_tb] = (TGraph*)File_InPut->Get(gName.c_str());

      gName = Form("g_mRanSigma_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mRanSigma[i_layer][i_tb] = (TGraph*)File_InPut->Get(gName.c_str());
    }
  }

  TH1F *h_play = new TH1F("h_play","h_play",2500,-500.5,1999.5);
  for(int i_bin = 0; i_bin < 2000; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-1000.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("ch");
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(0*128-50,8*128+50);
  h_play->GetYaxis()->SetTitle("noise adc");
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetRangeUser(0.0,50);

  TCanvas *c_NoiseSum = new TCanvas("c_NoiseSum","c_NoiseSum",10,10,900,900);
  c_NoiseSum->Divide(3,3);
  for(int i_pad = 0; i_pad < 9; ++i_pad)
  {
    c_NoiseSum->cd(i_pad+1)->SetLeftMargin(0.15);
    c_NoiseSum->cd(i_pad+1)->SetBottomMargin(0.15);
    c_NoiseSum->cd(i_pad+1)->SetTicks(1,1);
    c_NoiseSum->cd(i_pad+1)->SetGrid(0,0);
  }

  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);

    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      c_NoiseSum->cd(i_tb+1);
      string title = Form("FST Noise @ TimeBin%d",i_tb);
      if(i_layer > 0) title = Form("IST%d Noise @ TimeBin%d",i_layer,i_tb);
      if(i_layer > 0) h_play->GetXaxis()->SetRangeUser(0,128*12);
      h_play->SetTitle(title.c_str());
      h_play->DrawCopy("pE");

      g_mPedSigma[i_layer][i_tb]->SetMarkerStyle(20);
      g_mPedSigma[i_layer][i_tb]->SetMarkerColor(kGray+1);
      g_mPedSigma[i_layer][i_tb]->SetMarkerSize(1.0);
      g_mPedSigma[i_layer][i_tb]->Draw("P same");

      g_mCMNSigma[i_layer][i_tb]->SetMarkerStyle(24);
      g_mCMNSigma[i_layer][i_tb]->SetMarkerColor(2);
      g_mCMNSigma[i_layer][i_tb]->SetMarkerSize(1.0);
      g_mCMNSigma[i_layer][i_tb]->Draw("P same");

      g_mRanSigma[i_layer][i_tb]->SetMarkerStyle(24);
      g_mRanSigma[i_layer][i_tb]->SetMarkerColor(4);
      g_mRanSigma[i_layer][i_tb]->SetMarkerSize(1.0);
      g_mRanSigma[i_layer][i_tb]->Draw("P same");
      if(i_tb == 0)
      {
	leg->AddEntry(g_mPedSigma[i_layer][i_tb],"Total Noise","P");
	leg->AddEntry(g_mCMNSigma[i_layer][i_tb],"CMN","P");
	leg->AddEntry(g_mRanSigma[i_layer][i_tb],"Differential Noise","P");
	leg->Draw("same");
      }
    }
    c_NoiseSum->Update();
    c_NoiseSum->Print(outputname.c_str());
  }

  TH1F *h_mMeanPedSigma_RStrip[FST::numRStrip];
  TH1F *h_mMeanCMNSigma_RStrip[FST::numRStrip];
  TH1F *h_mMeanRanSigma_RStrip[FST::numRStrip];
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    string HistName; 
    HistName = Form("h_mMeanPedSigma_RStrip%d",i_rstrip);
    h_mMeanPedSigma_RStrip[i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numTBins,-0.5,FST::numTBins-0.5);

    HistName = Form("h_mMeanCMNSigma_RStrip%d",i_rstrip);
    h_mMeanCMNSigma_RStrip[i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numTBins,-0.5,FST::numTBins-0.5);

    HistName = Form("h_mMeanRanSigma_RStrip%d",i_rstrip);
    h_mMeanRanSigma_RStrip[i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numTBins,-0.5,FST::numTBins-0.5);
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      int counter_Ped = 0;
      double sum_Ped = 0.0;
      int counter_CMN = 0;
      double sum_CMN = 0.0;
      int counter_Ran = 0;
      double sum_Ran = 0.0;
      for(int i_bin = 0; i_bin < FST::numPhiSeg; ++i_bin)
      // for(int i_bin = 0; i_bin < 64; ++i_bin)
      {
	sum_Ped += h_mPedSigma_FST[i_rstrip][i_tb]->GetBinContent(i_bin+1);
	counter_Ped++;
	sum_CMN += h_mCMNSigma_FST[i_rstrip][i_tb]->GetBinContent(i_bin+1);
	counter_CMN++;
	sum_Ran += h_mRanSigma_FST[i_rstrip][i_tb]->GetBinContent(i_bin+1);
	counter_Ran++;
      }
      if(counter_Ped > 0) h_mMeanPedSigma_RStrip[i_rstrip]->SetBinContent(i_tb+1,sum_Ped/counter_Ped);
      if(counter_CMN > 0) h_mMeanCMNSigma_RStrip[i_rstrip]->SetBinContent(i_tb+1,sum_CMN/counter_CMN);
      if(counter_Ran > 0) h_mMeanRanSigma_RStrip[i_rstrip]->SetBinContent(i_tb+1,sum_Ran/counter_Ran);
    }
  }
  TCanvas *c_NoiseMean = new TCanvas("c_NoiseMean","c_NoiseMean",10,10,800,800);
  c_NoiseMean->Divide(1,4);
  for(int i_pad = 0; i_pad < 4; ++i_pad)
  {
    c_NoiseMean->cd(i_pad+1)->SetLeftMargin(0.15);
    c_NoiseMean->cd(i_pad+1)->SetBottomMargin(0.15);
    c_NoiseMean->cd(i_pad+1)->SetTicks(1,1);
    c_NoiseMean->cd(i_pad+1)->SetGrid(0,0);
  }

  c_NoiseMean->cd(1);
  TLegend *leg_mean = new TLegend(0.7,0.2,0.8,0.5);
  leg_mean->SetBorderSize(0);
  leg_mean->SetFillColor(10);
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    h_mMeanPedSigma_RStrip[i_rstrip]->SetTitle("FST Total Noise");
    h_mMeanPedSigma_RStrip[i_rstrip]->SetStats(0);
    h_mMeanPedSigma_RStrip[i_rstrip]->GetXaxis()->SetTitle("Time Bin");
    h_mMeanPedSigma_RStrip[i_rstrip]->GetXaxis()->SetTitleSize(0.06);
    h_mMeanPedSigma_RStrip[i_rstrip]->GetXaxis()->SetLabelSize(0.06);
    h_mMeanPedSigma_RStrip[i_rstrip]->GetYaxis()->SetTitle("<Noise_{Total}>");
    h_mMeanPedSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleSize(0.10);
    h_mMeanPedSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleOffset(0.5);
    h_mMeanPedSigma_RStrip[i_rstrip]->GetYaxis()->SetRangeUser(0.0,20.0);
    h_mMeanPedSigma_RStrip[i_rstrip]->GetYaxis()->SetLabelSize(0.08);
    h_mMeanPedSigma_RStrip[i_rstrip]->SetLineColor(i_rstrip+1);

    if(i_rstrip == 0) h_mMeanPedSigma_RStrip[i_rstrip]->Draw();
    else h_mMeanPedSigma_RStrip[i_rstrip]->Draw("same");

    string LegName = Form("R_strip %d",i_rstrip);
    leg_mean->AddEntry(h_mMeanRanSigma_RStrip[i_rstrip],LegName.c_str(),"L");
  }
  leg_mean->Draw("same");

  c_NoiseMean->cd(2);
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    h_mMeanCMNSigma_RStrip[i_rstrip]->SetTitle("FST Commen Mode Noise");
    h_mMeanCMNSigma_RStrip[i_rstrip]->SetStats(0);
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetXaxis()->SetTitle("Time Bin");
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetXaxis()->SetTitleSize(0.06);
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetXaxis()->SetLabelSize(0.06);
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetTitle("<Noise_{CMN}>");
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleSize(0.10);
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleOffset(0.5);
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetRangeUser(0.0,20.0);
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetLabelSize(0.08);
    h_mMeanCMNSigma_RStrip[i_rstrip]->SetLineColor(i_rstrip+1);

    if(i_rstrip == 0) h_mMeanCMNSigma_RStrip[i_rstrip]->Draw();
    else h_mMeanCMNSigma_RStrip[i_rstrip]->Draw("same");
  }

  c_NoiseMean->cd(3);
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    h_mMeanRanSigma_RStrip[i_rstrip]->SetTitle("FST Differential Noise");
    h_mMeanRanSigma_RStrip[i_rstrip]->SetStats(0);
    h_mMeanRanSigma_RStrip[i_rstrip]->GetXaxis()->SetTitle("Time Bin");
    h_mMeanRanSigma_RStrip[i_rstrip]->GetXaxis()->SetTitleSize(0.06);
    h_mMeanRanSigma_RStrip[i_rstrip]->GetXaxis()->SetLabelSize(0.06);
    h_mMeanRanSigma_RStrip[i_rstrip]->GetYaxis()->SetTitle("<Noise_{Ran}>");
    h_mMeanRanSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleSize(0.10);
    h_mMeanRanSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleOffset(0.5);
    h_mMeanRanSigma_RStrip[i_rstrip]->GetYaxis()->SetRangeUser(0.0,15.0);
    h_mMeanRanSigma_RStrip[i_rstrip]->GetYaxis()->SetLabelSize(0.08);
    h_mMeanRanSigma_RStrip[i_rstrip]->SetLineColor(i_rstrip+1);

    if(i_rstrip == 0) h_mMeanRanSigma_RStrip[i_rstrip]->Draw();
    else h_mMeanRanSigma_RStrip[i_rstrip]->Draw("same");
  }
  c_NoiseMean->Update();
  c_NoiseMean->Print(outputname.c_str());

  // mean IST
  TH1F *h_mMeanPedSigma_IST[3];
  TH1F *h_mMeanCMNSigma_IST[3];
  TH1F *h_mMeanRanSigma_IST[3];
  for(int i_layer = 1; i_layer < 4; ++i_layer)
  {
    string HistName = Form("h_mMeanPedSigma_IST%d",i_layer);
    h_mMeanPedSigma_IST[i_layer-1] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numTBins,-0.5,FST::numTBins-0.5);

    HistName = Form("h_mMeanCMNSigma_IST%d",i_layer);
    h_mMeanCMNSigma_IST[i_layer-1] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numTBins,-0.5,FST::numTBins-0.5);

    HistName = Form("h_mMeanRanSigma_IST%d",i_layer);
    h_mMeanRanSigma_IST[i_layer-1] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numTBins,-0.5,FST::numTBins-0.5);
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      double ch = -1.0;
      double noise_ran = -1.0;
      double noise_cmn = -1.0;
      double noise_tot = -1.0;
      int counts = 0;
      double sum_ran = 0.0; double sqr_ran = 0.0;
      double sum_cmn = 0.0; double sqr_cmn = 0.0;
      double sum_tot = 0.0; double sqr_tot = 0.0;
      for(int i_point = 0; i_point < g_mRanSigma[i_layer][i_tb]->GetN(); ++i_point)
      {
	g_mRanSigma[i_layer][i_tb]->GetPoint(i_point,ch,noise_ran);
	g_mCMNSigma[i_layer][i_tb]->GetPoint(i_point,ch,noise_cmn);
	g_mPedSigma[i_layer][i_tb]->GetPoint(i_point,ch,noise_tot);
	if(noise_tot > 0)
	{
	  sum_ran += noise_ran;
	  sqr_ran += noise_ran*noise_ran;
	  sum_cmn += noise_cmn;
	  sqr_cmn += noise_cmn*noise_cmn;
	  sum_tot += noise_tot;
	  sqr_tot += noise_tot*noise_tot;
	  counts++;
	}
      }
      // cout << "i_layer = " << i_layer << ", i_tb = " << i_tb << ", count_even = " << count_even << ", count_odd = " << count_odd << endl;
      if(counts > 0) 
      {
	h_mMeanRanSigma_IST[i_layer-1]->SetBinContent(i_tb+1,sum_ran/counts);
	// double std_ran = sqrt((sqr_ran-sum_ran*sum_ran/(double)counts)/(double)(counts-1));
	// h_mMeanRanSigma_IST[i_layer]->SetBinError(i_tb+1,std_ran/sqrt(counts));

	h_mMeanCMNSigma_IST[i_layer-1]->SetBinContent(i_tb+1,sum_cmn/counts);
	// double std_cmn = sqrt((sqr_cmn-sum_cmn*sum_cmn/(double)counts)/(double)(counts-1));
	// h_mMeanCMNSigma_IST[i_layer]->SetBinError(i_tb+1,std_cmn/sqrt(counts));

	h_mMeanPedSigma_IST[i_layer-1]->SetBinContent(i_tb+1,sum_tot/counts);
	// double std_tot = sqrt((sqr_tot-sum_tot*sum_tot/(double)counts)/(double)(counts-1));
	// h_mMeanPedSigma_IST[i_layer]->SetBinError(i_tb+1,std_tot/sqrt(counts));
      }
    }
  }
  
  for(int i_layer = 1; i_layer < 4; ++i_layer)
  {
    c_NoiseMean->cd(1);
    string title = Form("IST%d Total Noise",i_layer);
    h_mMeanPedSigma_IST[i_layer-1]->SetTitle(title.c_str());
    h_mMeanPedSigma_IST[i_layer-1]->SetStats(0);
    h_mMeanPedSigma_IST[i_layer-1]->GetXaxis()->SetTitle("Time Bin");
    h_mMeanPedSigma_IST[i_layer-1]->GetXaxis()->SetTitleSize(0.06);
    h_mMeanPedSigma_IST[i_layer-1]->GetXaxis()->SetLabelSize(0.06);
    h_mMeanPedSigma_IST[i_layer-1]->GetYaxis()->SetTitle("<Noise_{Total}>");
    h_mMeanPedSigma_IST[i_layer-1]->GetYaxis()->SetTitleSize(0.10);
    h_mMeanPedSigma_IST[i_layer-1]->GetYaxis()->SetTitleOffset(0.5);
    h_mMeanPedSigma_IST[i_layer-1]->GetYaxis()->SetRangeUser(0.0,30.0);
    h_mMeanPedSigma_IST[i_layer-1]->GetYaxis()->SetNdivisions(505);
    h_mMeanPedSigma_IST[i_layer-1]->GetYaxis()->SetLabelSize(0.08);
    h_mMeanPedSigma_IST[i_layer-1]->SetLineColor(1);
    h_mMeanPedSigma_IST[i_layer-1]->Draw("h");

    c_NoiseMean->cd(2);
    title = Form("IST%d Common Mode Noise",i_layer);
    h_mMeanCMNSigma_IST[i_layer-1]->SetTitle(title.c_str());
    h_mMeanCMNSigma_IST[i_layer-1]->SetStats(0);
    h_mMeanCMNSigma_IST[i_layer-1]->GetXaxis()->SetTitle("Time Bin");
    h_mMeanCMNSigma_IST[i_layer-1]->GetXaxis()->SetTitleSize(0.06);
    h_mMeanCMNSigma_IST[i_layer-1]->GetXaxis()->SetLabelSize(0.06);
    h_mMeanCMNSigma_IST[i_layer-1]->GetYaxis()->SetTitle("<Noise_{CMN}>");
    h_mMeanCMNSigma_IST[i_layer-1]->GetYaxis()->SetTitleSize(0.10);
    h_mMeanCMNSigma_IST[i_layer-1]->GetYaxis()->SetTitleOffset(0.5);
    h_mMeanCMNSigma_IST[i_layer-1]->GetYaxis()->SetRangeUser(0.0,30.0);
    h_mMeanCMNSigma_IST[i_layer-1]->GetYaxis()->SetNdivisions(505);
    h_mMeanCMNSigma_IST[i_layer-1]->GetYaxis()->SetLabelSize(0.08);
    h_mMeanCMNSigma_IST[i_layer-1]->SetLineColor(1);
    h_mMeanCMNSigma_IST[i_layer-1]->Draw("h");

    c_NoiseMean->cd(3);
    title = Form("IST%d Differential Noise",i_layer);
    h_mMeanRanSigma_IST[i_layer-1]->SetTitle(title.c_str());
    h_mMeanRanSigma_IST[i_layer-1]->SetStats(0);
    h_mMeanRanSigma_IST[i_layer-1]->GetXaxis()->SetTitle("Time Bin");
    h_mMeanRanSigma_IST[i_layer-1]->GetXaxis()->SetTitleSize(0.06);
    h_mMeanRanSigma_IST[i_layer-1]->GetXaxis()->SetLabelSize(0.06);
    h_mMeanRanSigma_IST[i_layer-1]->GetYaxis()->SetTitle("<Noise_{Ran}>");
    h_mMeanRanSigma_IST[i_layer-1]->GetYaxis()->SetTitleSize(0.10);
    h_mMeanRanSigma_IST[i_layer-1]->GetYaxis()->SetTitleOffset(0.5);
    h_mMeanRanSigma_IST[i_layer-1]->GetYaxis()->SetRangeUser(0.0,30.0);
    h_mMeanRanSigma_IST[i_layer-1]->GetYaxis()->SetNdivisions(505);
    h_mMeanRanSigma_IST[i_layer-1]->GetYaxis()->SetLabelSize(0.08);
    h_mMeanRanSigma_IST[i_layer-1]->SetLineColor(1);
    h_mMeanRanSigma_IST[i_layer-1]->Draw("h");

    c_NoiseMean->Update();
    c_NoiseMean->Print(outputname.c_str());
  }

  TH1F *h_ratio[FST::numRStrip][FST::numTBins];
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      string HistName = Form("h_ratio_RStrip%d_TimeBin%d",i_rstrip,i_tb);
      h_ratio[i_rstrip][i_tb] = (TH1F*)h_mRanSigma_FST[i_rstrip][i_tb]->Clone(HistName.c_str());
      h_ratio[i_rstrip][i_tb]->SetTitle(HistName.c_str());
      h_ratio[i_rstrip][i_tb]->Reset();
      h_ratio[i_rstrip][i_tb]->Divide(h_mRanSigma_FST[i_rstrip][i_tb],h_mPedSigma_FST[i_rstrip][i_tb],1,1,"B");
      int counter = 0;
      double sum = 0.0;
      for(int i_bin = 0; i_bin < h_ratio[i_rstrip][i_tb]->GetNbinsX(); ++i_bin)
      {
	sum += h_ratio[i_rstrip][i_tb]->GetBinContent(i_bin+1);
	counter++;
      }
    }
  }
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      int i_pad = i_rstrip*9 + i_tb;
      if(i_rstrip > 3) i_pad = (i_rstrip-4)*9 + i_tb;
      c_Noise->cd(i_pad+1);
      string title = Form("Noise_{Diff}/Noise_{total} RStrip%d & TB%d",i_rstrip,i_tb);
      h_ratio[i_rstrip][i_tb]->SetTitle(title.c_str());
      h_ratio[i_rstrip][i_tb]->SetStats(0);
      h_ratio[i_rstrip][i_tb]->GetXaxis()->SetLabelSize(0.06);
      h_ratio[i_rstrip][i_tb]->GetYaxis()->SetRangeUser(0,1.0);
      h_ratio[i_rstrip][i_tb]->GetYaxis()->SetLabelSize(0.06);
      h_ratio[i_rstrip][i_tb]->Draw();
      PlotLine(31.5, 31.5,0,1.0,2,1,2);
      PlotLine(63.5, 63.5,0,1.0,2,1,2);
      PlotLine(95.5, 95.5,0,1.0,2,1,2);
    }
    if(i_rstrip == 3 || i_rstrip == FST::numRStrip-1)
    {
      c_Noise->Update();
      c_Noise->Print(outputname.c_str());
    }
  }

  // get Ran/Total ratio
  TH1F *h_meanRatio_Rstrip[FST::numRStrip];
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    string HistName = Form("h_meanRatio_Rstrip%d",i_rstrip);
    h_meanRatio_Rstrip[i_rstrip] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numTBins,-0.5,FST::numTBins-0.5);
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      int counter = 0;
      double sum = 0.0;
      for(int i_bin = 0; i_bin < FST::numPhiSeg; ++i_bin)
      // for(int i_bin = 0; i_bin < 64; ++i_bin)
      {
	sum += h_ratio[i_rstrip][i_tb]->GetBinContent(i_bin+1);
	counter++;
      }
      h_meanRatio_Rstrip[i_rstrip]->SetBinContent(i_tb+1,sum/counter);
    }
  }
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    c_NoiseMean->cd(1);
    h_meanRatio_Rstrip[i_rstrip]->SetTitle("FST Noise");
    h_meanRatio_Rstrip[i_rstrip]->SetStats(0);
    h_meanRatio_Rstrip[i_rstrip]->GetXaxis()->SetTitle("Time Bin");
    h_meanRatio_Rstrip[i_rstrip]->GetXaxis()->SetTitleSize(0.06);
    h_meanRatio_Rstrip[i_rstrip]->GetXaxis()->SetLabelSize(0.06);
    h_meanRatio_Rstrip[i_rstrip]->GetYaxis()->SetTitle("Noise_{Diff}/Noise_{Total}");
    h_meanRatio_Rstrip[i_rstrip]->GetYaxis()->SetTitleSize(0.10);
    h_meanRatio_Rstrip[i_rstrip]->GetYaxis()->SetTitleOffset(0.5);
    h_meanRatio_Rstrip[i_rstrip]->GetYaxis()->SetRangeUser(0.0,1.0);
    h_meanRatio_Rstrip[i_rstrip]->GetYaxis()->SetNdivisions(505);
    h_meanRatio_Rstrip[i_rstrip]->GetYaxis()->SetLabelSize(0.08);
    h_meanRatio_Rstrip[i_rstrip]->SetLineColor(i_rstrip+1);

    if(i_rstrip == 0) h_meanRatio_Rstrip[i_rstrip]->Draw();
    else h_meanRatio_Rstrip[i_rstrip]->Draw("same");
  }
  leg_mean->Draw("same");

  TH1F *h_meanRatio_IST[3]; // even col and odd col
  for(int i_layer = 1; i_layer < 4; ++i_layer)
  {
    string HistName = Form("h_meanRatio_IST%d",i_layer);
    h_meanRatio_IST[i_layer-1] = new TH1F(HistName.c_str(),HistName.c_str(),FST::numTBins,-0.5,FST::numTBins-0.5);
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      double ch = -1.0;
      double noise_ran = -1.0;
      double noise_tot = -1.0;
      int counts = 0;
      double sum = 0.0;
      for(int i_point = 0; i_point < g_mRanSigma[i_layer][i_tb]->GetN(); ++i_point)
      {
	g_mRanSigma[i_layer][i_tb]->GetPoint(i_point,ch,noise_ran);
	g_mPedSigma[i_layer][i_tb]->GetPoint(i_point,ch,noise_tot);
	if(noise_tot > 0)
	{
	  sum += noise_ran/noise_tot;
	  counts++;
	}
      }
      if(counts > 0) h_meanRatio_IST[i_layer-1]->SetBinContent(i_tb+1,sum/counts);
    }
  }
  
  for(int i_layer = 1; i_layer < 4; ++i_layer)
  {
    c_NoiseMean->cd(i_layer+1);
    string title = Form("IST%d Noise",i_layer);
    h_meanRatio_IST[i_layer-1]->SetTitle(title.c_str());
    h_meanRatio_IST[i_layer-1]->SetStats(0);
    h_meanRatio_IST[i_layer-1]->GetXaxis()->SetTitle("Time Bin");
    h_meanRatio_IST[i_layer-1]->GetXaxis()->SetTitleSize(0.06);
    h_meanRatio_IST[i_layer-1]->GetXaxis()->SetLabelSize(0.06);
    h_meanRatio_IST[i_layer-1]->GetYaxis()->SetTitle("Noise_{Diff}/Noise_{Total}");
    h_meanRatio_IST[i_layer-1]->GetYaxis()->SetTitleSize(0.10);
    h_meanRatio_IST[i_layer-1]->GetYaxis()->SetTitleOffset(0.5);
    h_meanRatio_IST[i_layer-1]->GetYaxis()->SetRangeUser(0.0,1.0);
    h_meanRatio_IST[i_layer-1]->GetYaxis()->SetNdivisions(505);
    h_meanRatio_IST[i_layer-1]->GetYaxis()->SetLabelSize(0.08);
    h_meanRatio_IST[i_layer-1]->SetLineColor(1);
    h_meanRatio_IST[i_layer-1]->Draw();
  }
  c_NoiseMean->Update();
  c_NoiseMean->Print(outputname.c_str());
  
#if 0
  // noise in RO channel
  TGraph *g_mRoPedMean[4][FST::numTBins];
  TGraph *g_mRoPedSigma[4][FST::numTBins];
  TGraph *g_mRoCMNSigma[4][FST::numTBins];
  TGraph *g_mRoRanSigma[4][FST::numTBins];
  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      std::string gName = Form("g_mRoPedMean_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mRoPedMean[i_layer][i_tb] = (TGraph*)File_InPut->Get(gName.c_str());

      gName = Form("g_mRoPedSigma_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mRoPedSigma[i_layer][i_tb] = (TGraph*)File_InPut->Get(gName.c_str()); 

      gName = Form("g_mRoCMNSigma_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mRoCMNSigma[i_layer][i_tb] = (TGraph*)File_InPut->Get(gName.c_str());

      gName = Form("g_mRoRanSigma_Layer%d_TimeBin%d",i_layer,i_tb);
      g_mRoRanSigma[i_layer][i_tb] = (TGraph*)File_InPut->Get(gName.c_str());
    }
  }

  TH1F *h_play_RO = new TH1F("h_play_RO","h_play_RO",15000,-500.5,14499.5);
  for(int i_bin = 0; i_bin < 15000; ++i_bin)
  {
    h_play_RO->SetBinContent(i_bin+1,-1000.0);
    h_play_RO->SetBinError(i_bin+1,1.0);
  }
  h_play_RO->SetStats(0);
  h_play_RO->GetXaxis()->SetTitle("RO Order");
  h_play_RO->GetXaxis()->SetTitleSize(0.06);
  // h_play_RO->GetXaxis()->SetRangeUser(4*1152-50,6*1152+50);
  h_play_RO->GetYaxis()->SetTitleSize(0.06);
  h_play_RO->GetYaxis()->SetRangeUser(0.0,50);

  TCanvas *c_NoiseRO = new TCanvas("c_NoiseRO","c_NoiseRO",10,10,900,900);
  c_NoiseRO->Divide(1,4);
  for(int i_pad = 0; i_pad < 4; ++i_pad)
  {
    c_NoiseRO->cd(i_pad+1)->SetLeftMargin(0.10);
    c_NoiseRO->cd(i_pad+1)->SetBottomMargin(0.15);
    c_NoiseRO->cd(i_pad+1)->SetTicks(1,1);
    c_NoiseRO->cd(i_pad+1)->SetGrid(0,0);
  }

  const string mDetName[4] = {"FST","IST1","IST2","IST3"};
  for(int i_layer = 0; i_layer < 4; ++i_layer)
  {
    if(i_layer == 2) continue;

    if(i_layer == 0)h_play_RO->GetXaxis()->SetRangeUser(4*1152-500,6*1152+500);
    if(i_layer != 0)h_play_RO->GetXaxis()->SetRangeUser(0*1152-500,12*1152+500);

    c_NoiseRO->cd(1);
    string titleRO = Form("%s Pedestal vs. RO Order",mDetName[i_layer].c_str());
    h_play_RO->SetTitle(titleRO.c_str());
    if(i_layer == 0) h_play_RO->GetYaxis()->SetRangeUser(0.0,1200.0);
    if(i_layer != 0) h_play_RO->GetYaxis()->SetRangeUser(0.0,2400.0);
    h_play_RO->GetYaxis()->SetTitle("Pedestal ADC");
    h_play_RO->DrawCopy("pE");

    TLegend *legRO = new TLegend(0.4,0.6,0.85,0.85);
    legRO->SetBorderSize(0);
    legRO->SetFillColor(10);
    legRO->SetNColumns(3);
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      g_mRoPedMean[i_layer][i_tb]->SetMarkerStyle(24);
      g_mRoPedMean[i_layer][i_tb]->SetMarkerColor(i_tb+1);
      g_mRoPedMean[i_layer][i_tb]->SetMarkerSize(1.0);
      g_mRoPedMean[i_layer][i_tb]->Draw("P same");
      string leg_ro = Form("TimeBin %d",i_tb);
      legRO->AddEntry(g_mRoPedMean[i_layer][i_tb],leg_ro.c_str(),"P");
    }
    legRO->Draw("same");
    if(i_layer == 0) PlotLine(5*1152,5*1152,0.0,1200,2,1,2);
    if(i_layer != 0)
    {
      for(int i_apv = 0; i_apv < 13; ++i_apv)
      {
	PlotLine(i_apv*1152,i_apv*1152,0.0,1200,2,1,2);
      }
    }

    c_NoiseRO->cd(2);
    titleRO = Form("%s Total Noise_ vs. RO Order",mDetName[i_layer].c_str());
    h_play_RO->SetTitle(titleRO.c_str());
    h_play_RO->GetYaxis()->SetRangeUser(0.0,50.0);
    h_play_RO->GetYaxis()->SetTitle("Total Noise ADC");
    h_play_RO->DrawCopy("pE");
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      g_mRoPedSigma[i_layer][i_tb]->SetMarkerStyle(24);
      g_mRoPedSigma[i_layer][i_tb]->SetMarkerColor(i_tb+1);
      g_mRoPedSigma[i_layer][i_tb]->SetMarkerSize(1.0);
      g_mRoPedSigma[i_layer][i_tb]->Draw("P same");
    }
    if(i_layer == 0) PlotLine(5*1152,5*1152,0.0,50,2,1,2);
    if(i_layer != 0)
    {
      for(int i_apv = 0; i_apv < 13; ++i_apv)
      {
	PlotLine(i_apv*1152,i_apv*1152,0.0,50,2,1,2);
      }
    }

    c_NoiseRO->cd(3);
    titleRO = Form("%s Common Mode Noise_ vs. RO Order",mDetName[i_layer].c_str());
    h_play_RO->SetTitle(titleRO.c_str());
    h_play_RO->GetYaxis()->SetRangeUser(0.0,50.0);
    h_play_RO->GetYaxis()->SetTitle("Common Mode Noise ADC");
    h_play_RO->DrawCopy("pE");
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      g_mRoCMNSigma[i_layer][i_tb]->SetMarkerStyle(24);
      g_mRoCMNSigma[i_layer][i_tb]->SetMarkerColor(i_tb+1);
      g_mRoCMNSigma[i_layer][i_tb]->SetMarkerSize(1.0);
      g_mRoCMNSigma[i_layer][i_tb]->Draw("P same");
    }
    if(i_layer == 0) PlotLine(5*1152,5*1152,0.0,50,2,1,2);
    if(i_layer != 0)
    {
      for(int i_apv = 0; i_apv < 13; ++i_apv)
      {
	PlotLine(i_apv*1152,i_apv*1152,0.0,50,2,1,2);
      }
    }

    c_NoiseRO->cd(4);
    titleRO = Form("%s Random Noise_ vs. RO Order",mDetName[i_layer].c_str());
    h_play_RO->SetTitle(titleRO.c_str());
    h_play_RO->GetYaxis()->SetRangeUser(0.0,50.0);
    h_play_RO->GetYaxis()->SetTitle("Differential Noise ADC");
    h_play_RO->DrawCopy("pE");
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      g_mRoRanSigma[i_layer][i_tb]->SetMarkerStyle(24);
      g_mRoRanSigma[i_layer][i_tb]->SetMarkerColor(i_tb+1);
      g_mRoRanSigma[i_layer][i_tb]->SetMarkerSize(1.0);
      g_mRoRanSigma[i_layer][i_tb]->Draw("P same");
    }
    if(i_layer == 0) PlotLine(5*1152,5*1152,0.0,50,2,1,2);
    if(i_layer != 0)
    {
      for(int i_apv = 0; i_apv < 13; ++i_apv)
      {
	PlotLine(i_apv*1152,i_apv*1152,0.0,50,2,1,2);
      }
    }

    c_NoiseRO->Update();
    c_NoiseRO->Print(outputname.c_str());

    // for each time bin
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      c_NoiseRO->cd(1);
      titleRO = Form("%s Pedestal vs. RO Order at TB %d",mDetName[i_layer].c_str(),i_tb);
      h_play_RO->SetTitle(titleRO.c_str());
      if(i_layer == 0) h_play_RO->GetYaxis()->SetRangeUser(0.0,1200.0);
      if(i_layer != 0) h_play_RO->GetYaxis()->SetRangeUser(0.0,2400.0);
      h_play_RO->GetYaxis()->SetTitle("Pedestal ADC");
      h_play_RO->DrawCopy("pE");
      g_mRoPedMean[i_layer][i_tb]->SetMarkerStyle(24);
      g_mRoPedMean[i_layer][i_tb]->SetMarkerColor(i_tb+1);
      g_mRoPedMean[i_layer][i_tb]->SetMarkerSize(1.0);
      g_mRoPedMean[i_layer][i_tb]->Draw("P same");
      if(i_layer == 0) PlotLine(5*1152,5*1152,0.0,1200,2,1,2);
      if(i_layer != 0)
      {
	for(int i_apv = 0; i_apv < 13; ++i_apv)
	{
	  PlotLine(i_apv*1152,i_apv*1152,0.0,1200,2,1,2);
	}
      }

      c_NoiseRO->cd(2);
      titleRO = Form("%s Total Noise vs. RO Order at TB %d",mDetName[i_layer].c_str(),i_tb);
      h_play_RO->SetTitle(titleRO.c_str());
      h_play_RO->GetYaxis()->SetRangeUser(0.0,50.0);
      h_play_RO->GetYaxis()->SetTitle("Total Noise ADC");
      h_play_RO->DrawCopy("pE");
      g_mRoPedSigma[i_layer][i_tb]->SetMarkerStyle(24);
      g_mRoPedSigma[i_layer][i_tb]->SetMarkerColor(i_tb+1);
      g_mRoPedSigma[i_layer][i_tb]->SetMarkerSize(1.0);
      g_mRoPedSigma[i_layer][i_tb]->Draw("P same");
      if(i_layer == 0) PlotLine(5*1152,5*1152,0.0,50,2,1,2);
      if(i_layer != 0)
      {
	for(int i_apv = 0; i_apv < 13; ++i_apv)
	{
	  PlotLine(i_apv*1152,i_apv*1152,0.0,50,2,1,2);
	}
      }

      c_NoiseRO->cd(3);
      titleRO = Form("%s Common Mode Noise vs. RO Order at TB %d",mDetName[i_layer].c_str(),i_tb);
      h_play_RO->SetTitle(titleRO.c_str());
      h_play_RO->GetYaxis()->SetRangeUser(0.0,50.0);
      h_play_RO->GetYaxis()->SetTitle("Common Mode Noise ADC");
      h_play_RO->DrawCopy("pE");
      g_mRoCMNSigma[i_layer][i_tb]->SetMarkerStyle(24);
      g_mRoCMNSigma[i_layer][i_tb]->SetMarkerColor(i_tb+1);
      g_mRoCMNSigma[i_layer][i_tb]->SetMarkerSize(1.0);
      g_mRoCMNSigma[i_layer][i_tb]->Draw("P same");
      if(i_layer == 0) PlotLine(5*1152,5*1152,0.0,50,2,1,2);
      if(i_layer != 0)
      {
	for(int i_apv = 0; i_apv < 13; ++i_apv)
	{
	  PlotLine(i_apv*1152,i_apv*1152,0.0,50,2,1,2);
	}
      }

      c_NoiseRO->cd(4);
      titleRO = Form("%s Differential Noise vs. RO Order at TB %d",mDetName[i_layer].c_str(),i_tb);
      h_play_RO->SetTitle(titleRO.c_str());
      h_play_RO->GetYaxis()->SetRangeUser(0.0,50.0);
      h_play_RO->GetYaxis()->SetTitle("Differential Noise ADC");
      h_play_RO->DrawCopy("pE");
      g_mRoRanSigma[i_layer][i_tb]->SetMarkerStyle(24);
      g_mRoRanSigma[i_layer][i_tb]->SetMarkerColor(i_tb+1);
      g_mRoRanSigma[i_layer][i_tb]->SetMarkerSize(1.0);
      g_mRoRanSigma[i_layer][i_tb]->Draw("P same");
      if(i_layer == 0) PlotLine(5*1152,5*1152,0.0,50,2,1,2);
      if(i_layer != 0)
      {
	for(int i_apv = 0; i_apv < 13; ++i_apv)
	{
	  PlotLine(i_apv*1152,i_apv*1152,0.0,50,2,1,2);
	}
      }

      c_NoiseRO->Update();
      c_NoiseRO->Print(outputname.c_str());
    }
  }

  // Zoomed in View for all time bins
  {
    h_play_RO->GetXaxis()->SetRangeUser(4*1152,4*1152+9*16);

    c_NoiseRO->cd(1);
    h_play_RO->SetTitle("FST Pedestal vs. RO Order (First 16 channels)");
    h_play_RO->GetYaxis()->SetRangeUser(0.0,1200.0);
    h_play_RO->GetYaxis()->SetTitle("Pedestal ADC");
    h_play_RO->DrawCopy("pE");
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      g_mRoPedMean_FST[i_tb]->SetMarkerStyle(24);
      g_mRoPedMean_FST[i_tb]->SetMarkerColor(i_tb+1);
      g_mRoPedMean_FST[i_tb]->SetMarkerSize(1.0);
      g_mRoPedMean_FST[i_tb]->Draw("P same");
      string leg_ro = Form("TimeBin %d",i_tb);
    }
    legRO->Draw("same");

    c_NoiseRO->cd(2);
    h_play_RO->SetTitle("FST Total Noise vs. RO Order (First 16 channels)");
    h_play_RO->GetYaxis()->SetRangeUser(0.0,50.0);
    h_play_RO->GetYaxis()->SetTitle("Total Noise ADC");
    h_play_RO->DrawCopy("pE");
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      g_mRoPedSigma_FST[i_tb]->SetMarkerStyle(24);
      g_mRoPedSigma_FST[i_tb]->SetMarkerColor(i_tb+1);
      g_mRoPedSigma_FST[i_tb]->SetMarkerSize(1.0);
      g_mRoPedSigma_FST[i_tb]->Draw("P same");
    }

    c_NoiseRO->cd(3);
    h_play_RO->SetTitle("FST Common Mode Noise vs. RO Order (First 16 channels)");
    h_play_RO->GetYaxis()->SetRangeUser(0.0,50.0);
    h_play_RO->GetYaxis()->SetTitle("Common Mode Noise ADC");
    h_play_RO->DrawCopy("pE");
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      g_mRoCMNSigma_FST[i_tb]->SetMarkerStyle(24);
      g_mRoCMNSigma_FST[i_tb]->SetMarkerColor(i_tb+1);
      g_mRoCMNSigma_FST[i_tb]->SetMarkerSize(1.0);
      g_mRoCMNSigma_FST[i_tb]->Draw("P same");
    }

    c_NoiseRO->cd(4);
    h_play_RO->SetTitle("FST Random Noise vs. RO Order (First 16 channels)");
    h_play_RO->GetYaxis()->SetRangeUser(0.0,50.0);
    h_play_RO->GetYaxis()->SetTitle("Differential Noise ADC");
    h_play_RO->DrawCopy("pE");
    for(int i_tb = 0; i_tb < FST::numTBins; ++i_tb)
    {
      g_mRoRanSigma_FST[i_tb]->SetMarkerStyle(24);
      g_mRoRanSigma_FST[i_tb]->SetMarkerColor(i_tb+1);
      g_mRoRanSigma_FST[i_tb]->SetMarkerSize(1.0);
      g_mRoRanSigma_FST[i_tb]->Draw("P same");
    }

    c_NoiseRO->Update();
    c_NoiseRO->Print(outputname.c_str());
  }
#endif


  string output_stop = Form("./figures/Ist%s/IstChipNoiseQA_%s.pdf]",module.c_str(),module.c_str());
  c_Noise->Print(output_stop.c_str()); // open pdf file

  /*
  string outputfile = Form("../../output/noise/Ist%s/IstChipMeanNoise_%s.root",module.c_str(),module.c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_layer = 1; i_layer < 4; ++i_layer)
  {
    h_mMeanPedSigma_IST[i_layer-1]->Write();
    h_mMeanCMNSigma_IST[i_layer-1]->Write();
    h_mMeanRanSigma_IST[i_layer-1]->Write();
  }
  File_OutPut->Close();
  */
}

