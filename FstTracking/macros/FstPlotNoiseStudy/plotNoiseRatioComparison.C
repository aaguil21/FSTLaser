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

void plotNoiseRatioComparison(string module = "Mod03", string istMod = "Mod36")
{
  TH1F *h_mChipMeanPedSigma_RStrip[FST::numRStrip];
  TH1F *h_mChipMeanCMNSigma_RStrip[FST::numRStrip];
  TH1F *h_mChipMeanRanSigma_RStrip[FST::numRStrip];

  string inputInnerFst = Form("../../output/noise/%s/FstChipMeanNoise_%s_Inner.root",module.c_str(),module.c_str());
  TFile *File_InPutInnerFst = TFile::Open(inputInnerFst.c_str());
  for(int i_rstrip = 0; i_rstrip < 4; ++i_rstrip)
  {
    string HistName; 
    HistName = Form("h_mMeanPedSigma_RStrip%d",i_rstrip);
    h_mChipMeanPedSigma_RStrip[i_rstrip] = (TH1F*)File_InPutInnerFst->Get(HistName.c_str());

    HistName = Form("h_mMeanCMNSigma_RStrip%d",i_rstrip);
    h_mChipMeanCMNSigma_RStrip[i_rstrip] = (TH1F*)File_InPutInnerFst->Get(HistName.c_str());

    HistName = Form("h_mMeanRanSigma_RStrip%d",i_rstrip);
    h_mChipMeanRanSigma_RStrip[i_rstrip] = (TH1F*)File_InPutInnerFst->Get(HistName.c_str());
  }

  string inputOuterFst = Form("../../output/noise/%s/FstChipMeanNoise_%s_Outer.root",module.c_str(),module.c_str());
  TFile *File_InPutOuterFst = TFile::Open(inputOuterFst.c_str());
  for(int i_rstrip = 4; i_rstrip < 8; ++i_rstrip)
  {
    string HistName; 
    HistName = Form("h_mMeanPedSigma_RStrip%d",i_rstrip);
    h_mChipMeanPedSigma_RStrip[i_rstrip] = (TH1F*)File_InPutOuterFst->Get(HistName.c_str());

    HistName = Form("h_mMeanCMNSigma_RStrip%d",i_rstrip);
    h_mChipMeanCMNSigma_RStrip[i_rstrip] = (TH1F*)File_InPutOuterFst->Get(HistName.c_str());

    HistName = Form("h_mMeanRanSigma_RStrip%d",i_rstrip);
    h_mChipMeanRanSigma_RStrip[i_rstrip] = (TH1F*)File_InPutOuterFst->Get(HistName.c_str());
  }

  TH1F *h_mMeanPedSigma_RStrip[FST::numRStrip];
  TH1F *h_mMeanCMNSigma_RStrip[FST::numRStrip];
  TH1F *h_mMeanRanSigma_RStrip[FST::numRStrip];

  TH1F *h_mRatioPedSigma_RStrip[FST::numRStrip];
  TH1F *h_mRatioCMNSigma_RStrip[FST::numRStrip];
  TH1F *h_mRatioRanSigma_RStrip[FST::numRStrip];

  string inputFst = Form("../../output/noise/%s/FstDataMeanNoise_%s_HV70V_07122020.root",module.c_str(),module.c_str());
  TFile *File_InPutFst = TFile::Open(inputFst.c_str());
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    string HistName; 
    HistName = Form("h_mMeanPedSigma_RStrip%d",i_rstrip);
    h_mMeanPedSigma_RStrip[i_rstrip] = (TH1F*)File_InPutFst->Get(HistName.c_str());
    h_mRatioPedSigma_RStrip[i_rstrip] = (TH1F*)h_mMeanPedSigma_RStrip[i_rstrip]->Clone();
    h_mRatioPedSigma_RStrip[i_rstrip]->Divide(h_mChipMeanPedSigma_RStrip[i_rstrip]);

    HistName = Form("h_mMeanCMNSigma_RStrip%d",i_rstrip);
    h_mMeanCMNSigma_RStrip[i_rstrip] = (TH1F*)File_InPutFst->Get(HistName.c_str());
    h_mRatioCMNSigma_RStrip[i_rstrip] = (TH1F*)h_mMeanCMNSigma_RStrip[i_rstrip]->Clone();
    h_mRatioCMNSigma_RStrip[i_rstrip]->Divide(h_mChipMeanCMNSigma_RStrip[i_rstrip]);

    HistName = Form("h_mMeanRanSigma_RStrip%d",i_rstrip);
    h_mMeanRanSigma_RStrip[i_rstrip] = (TH1F*)File_InPutFst->Get(HistName.c_str());
    h_mRatioRanSigma_RStrip[i_rstrip] = (TH1F*)h_mMeanRanSigma_RStrip[i_rstrip]->Clone();
    h_mRatioRanSigma_RStrip[i_rstrip]->Divide(h_mChipMeanRanSigma_RStrip[i_rstrip]);
  }

  TH1F *h_mChipMeanPedSigma_IST[3];
  TH1F *h_mChipMeanCMNSigma_IST[3];
  TH1F *h_mChipMeanRanSigma_IST[3];
  string inputIst = Form("../../output/noise/Ist%s/IstChipMeanNoise_%s.root",istMod.c_str(),istMod.c_str());
  TFile *File_InPutIst = TFile::Open(inputIst.c_str());
  for(int i_layer = 1; i_layer < 4; ++i_layer)
  {
    string HistName = Form("h_mMeanPedSigma_IST%d",i_layer);
    h_mChipMeanPedSigma_IST[i_layer-1] = (TH1F*)File_InPutIst->Get(HistName.c_str());

    HistName = Form("h_mMeanCMNSigma_IST%d",i_layer);
    h_mChipMeanCMNSigma_IST[i_layer-1] = (TH1F*)File_InPutIst->Get(HistName.c_str());

    HistName = Form("h_mMeanRanSigma_IST%d",i_layer);
    h_mChipMeanRanSigma_IST[i_layer-1] = (TH1F*)File_InPutIst->Get(HistName.c_str());
  }

  TH1F *h_mMeanPedSigma_IST[3];
  TH1F *h_mMeanCMNSigma_IST[3];
  TH1F *h_mMeanRanSigma_IST[3];
  TH1F *h_mRatioPedSigma_IST[3];
  TH1F *h_mRatioCMNSigma_IST[3];
  TH1F *h_mRatioRanSigma_IST[3];
  string inputIstSensor = Form("../../output/noise/Ist%s/IstPedMeanNoise_%s_HV60V.root",istMod.c_str(),istMod.c_str());
  TFile *File_InPutIstSensor = TFile::Open(inputIstSensor.c_str());
  for(int i_layer = 1; i_layer < 4; ++i_layer)
  {
    string HistName = Form("h_mMeanPedSigma_IST%d",i_layer);
    h_mMeanPedSigma_IST[i_layer-1] = (TH1F*)File_InPutIstSensor->Get(HistName.c_str());
    h_mRatioPedSigma_IST[i_layer-1] = (TH1F*)h_mMeanPedSigma_IST[i_layer-1]->Clone();
    h_mRatioPedSigma_IST[i_layer-1]->Divide(h_mChipMeanPedSigma_IST[i_layer-1]);

    HistName = Form("h_mMeanCMNSigma_IST%d",i_layer);
    h_mMeanCMNSigma_IST[i_layer-1] = (TH1F*)File_InPutIstSensor->Get(HistName.c_str());
    h_mRatioCMNSigma_IST[i_layer-1] = (TH1F*)h_mMeanCMNSigma_IST[i_layer-1]->Clone();
    h_mRatioCMNSigma_IST[i_layer-1]->Divide(h_mChipMeanCMNSigma_IST[i_layer-1]);

    HistName = Form("h_mMeanRanSigma_IST%d",i_layer);
    h_mMeanRanSigma_IST[i_layer-1] = (TH1F*)File_InPutIstSensor->Get(HistName.c_str());
    h_mRatioRanSigma_IST[i_layer-1] = (TH1F*)h_mMeanRanSigma_IST[i_layer-1]->Clone();
    h_mRatioRanSigma_IST[i_layer-1]->Divide(h_mChipMeanRanSigma_IST[i_layer-1]);
  }

  TCanvas *c_NoiseRatio = new TCanvas("c_NoiseRatio","c_NoiseRatio",10,10,800,400);
  c_NoiseRatio->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_NoiseRatio->cd(i_pad+1)->SetLeftMargin(0.15);
    c_NoiseRatio->cd(i_pad+1)->SetBottomMargin(0.15);
    c_NoiseRatio->cd(i_pad+1)->SetTicks(1,1);
    c_NoiseRatio->cd(i_pad+1)->SetGrid(0,0);
  }

  TLegend *leg_mean = new TLegend(0.2,0.65,0.8,0.85);
  leg_mean->SetBorderSize(0);
  leg_mean->SetFillColor(10);
  leg_mean->SetNColumns(2);
  c_NoiseRatio->cd(1);
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    h_mRatioCMNSigma_RStrip[i_rstrip]->SetTitle("FST Commen Mode Noise");
    h_mRatioCMNSigma_RStrip[i_rstrip]->SetStats(0);
    h_mRatioCMNSigma_RStrip[i_rstrip]->GetXaxis()->SetTitle("Time Bin");
    h_mRatioCMNSigma_RStrip[i_rstrip]->GetXaxis()->SetTitleSize(0.04);
    h_mRatioCMNSigma_RStrip[i_rstrip]->GetXaxis()->SetLabelSize(0.04);
    h_mRatioCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetTitle("<Noise^{Sensor}_{CMN}>/<Noise^{Chip}_{CMN}>");
    h_mRatioCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleSize(0.06);
    h_mRatioCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleOffset(1.0);
    h_mRatioCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetLabelSize(0.04);
    h_mRatioCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetRangeUser(0.0,8.0);
    h_mRatioCMNSigma_RStrip[i_rstrip]->SetLineColor(i_rstrip+1);
    h_mRatioCMNSigma_RStrip[i_rstrip]->SetLineWidth(3);

    if(i_rstrip == 0) h_mRatioCMNSigma_RStrip[i_rstrip]->Draw();
    else h_mRatioCMNSigma_RStrip[i_rstrip]->Draw("same");

    string LegName = Form("R_strip %d",i_rstrip);
    leg_mean->AddEntry(h_mRatioCMNSigma_RStrip[i_rstrip],LegName.c_str(),"L");
  }

  h_mRatioCMNSigma_IST[1]->SetMarkerStyle(20);
  h_mRatioCMNSigma_IST[1]->SetMarkerColor(kGray+2);
  h_mRatioCMNSigma_IST[1]->SetMarkerSize(1.0);
  h_mRatioCMNSigma_IST[1]->Draw("p same");
  leg_mean->AddEntry(h_mRatioCMNSigma_IST[1],"IST 2013","p");
  leg_mean->Draw("same");

  c_NoiseRatio->cd(2);
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    h_mRatioRanSigma_RStrip[i_rstrip]->SetTitle("FST Differential Noise");
    h_mRatioRanSigma_RStrip[i_rstrip]->SetStats(0);
    h_mRatioRanSigma_RStrip[i_rstrip]->GetXaxis()->SetTitle("Time Bin");
    h_mRatioRanSigma_RStrip[i_rstrip]->GetXaxis()->SetTitleSize(0.04);
    h_mRatioRanSigma_RStrip[i_rstrip]->GetXaxis()->SetLabelSize(0.04);
    h_mRatioRanSigma_RStrip[i_rstrip]->GetYaxis()->SetTitle("<Noise^{Sensor}_{Ran}>/<Noise^{Chip}_{Ran}>");
    h_mRatioRanSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleSize(0.06);
    h_mRatioRanSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleOffset(1.0);
    h_mRatioRanSigma_RStrip[i_rstrip]->GetYaxis()->SetLabelSize(0.04);
    h_mRatioRanSigma_RStrip[i_rstrip]->GetYaxis()->SetRangeUser(0.0,8.0);
    h_mRatioRanSigma_RStrip[i_rstrip]->SetLineColor(i_rstrip+1);
    h_mRatioRanSigma_RStrip[i_rstrip]->SetLineWidth(3);

    if(i_rstrip == 0) h_mRatioRanSigma_RStrip[i_rstrip]->Draw();
    else h_mRatioRanSigma_RStrip[i_rstrip]->Draw("same");
  }

  h_mRatioRanSigma_IST[1]->SetMarkerStyle(20);
  h_mRatioRanSigma_IST[1]->SetMarkerColor(kGray+2);
  h_mRatioRanSigma_IST[1]->SetMarkerSize(1.0);
  h_mRatioRanSigma_IST[1]->Draw("p same");
  leg_mean->Draw("same");

  c_NoiseRatio->SaveAs("./figures/Comparison/NoiseRatioComparison.eps");
}
