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

void plotIstSensorPedNoiseComparison(string module = "Mod03", string istMod = "Mod36")
{
  TH1F *h_mMeanPedSigma_RStrip[FST::numRStrip];
  TH1F *h_mMeanCMNSigma_RStrip[FST::numRStrip];
  TH1F *h_mMeanRanSigma_RStrip[FST::numRStrip];
  TH1F *h_mMeanPedSigma_ist[3];
  TH1F *h_mMeanCMNSigma_ist[3];
  TH1F *h_mMeanRanSigma_ist[3];

  // string inputFst = Form("../../output/noise/%s/FstPedMeanNoise_%s_HV70V.root",module.c_str(),module.c_str());
  string inputFst = Form("../../output/noise/%s/FstDataMeanNoise_%s_HV70V_07122020.root",module.c_str(),module.c_str());
  TFile *File_InPutFst = TFile::Open(inputFst.c_str());
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    string HistName; 
    HistName = Form("h_mMeanPedSigma_RStrip%d",i_rstrip);
    h_mMeanPedSigma_RStrip[i_rstrip] = (TH1F*)File_InPutFst->Get(HistName.c_str());

    HistName = Form("h_mMeanCMNSigma_RStrip%d",i_rstrip);
    h_mMeanCMNSigma_RStrip[i_rstrip] = (TH1F*)File_InPutFst->Get(HistName.c_str());

    HistName = Form("h_mMeanRanSigma_RStrip%d",i_rstrip);
    h_mMeanRanSigma_RStrip[i_rstrip] = (TH1F*)File_InPutFst->Get(HistName.c_str());
  }
  for(int i_layer = 1; i_layer < 4; ++i_layer)
  {
    string HistName = Form("h_mMeanPedSigma_IST%d",i_layer);
    h_mMeanPedSigma_ist[i_layer-1] = (TH1F*)File_InPutFst->Get(HistName.c_str());

    HistName = Form("h_mMeanCMNSigma_IST%d",i_layer);
    h_mMeanCMNSigma_ist[i_layer-1] = (TH1F*)File_InPutFst->Get(HistName.c_str());

    HistName = Form("h_mMeanRanSigma_IST%d",i_layer);
    h_mMeanRanSigma_ist[i_layer-1] = (TH1F*)File_InPutFst->Get(HistName.c_str());
  }

  TH1F *h_mMeanPedSigma_IST[3];
  TH1F *h_mMeanCMNSigma_IST[3];
  TH1F *h_mMeanRanSigma_IST[3];
  string inputIst = Form("../../output/noise/Ist%s/IstPedMeanNoise_%s_HV60V.root",istMod.c_str(),istMod.c_str());
  TFile *File_InPutIst = TFile::Open(inputIst.c_str());
  for(int i_layer = 1; i_layer < 4; ++i_layer)
  {
    string HistName = Form("h_mMeanPedSigma_IST%d",i_layer);
    h_mMeanPedSigma_IST[i_layer-1] = (TH1F*)File_InPutIst->Get(HistName.c_str());

    HistName = Form("h_mMeanCMNSigma_IST%d",i_layer);
    h_mMeanCMNSigma_IST[i_layer-1] = (TH1F*)File_InPutIst->Get(HistName.c_str());

    HistName = Form("h_mMeanRanSigma_IST%d",i_layer);
    h_mMeanRanSigma_IST[i_layer-1] = (TH1F*)File_InPutIst->Get(HistName.c_str());
  }

  TCanvas *c_NoiseMean = new TCanvas("c_NoiseMean","c_NoiseMean",10,10,800,400);
  c_NoiseMean->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_NoiseMean->cd(i_pad+1)->SetLeftMargin(0.15);
    c_NoiseMean->cd(i_pad+1)->SetBottomMargin(0.15);
    c_NoiseMean->cd(i_pad+1)->SetTicks(1,1);
    c_NoiseMean->cd(i_pad+1)->SetGrid(0,0);
  }

  TLegend *leg_mean = new TLegend(0.2,0.65,0.8,0.85);
  leg_mean->SetBorderSize(0);
  leg_mean->SetFillColor(10);
  leg_mean->SetNColumns(2);
  c_NoiseMean->cd(1);
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    h_mMeanCMNSigma_RStrip[i_rstrip]->SetTitle("FST Commen Mode Noise");
    h_mMeanCMNSigma_RStrip[i_rstrip]->SetStats(0);
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetXaxis()->SetTitle("Time Bin");
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetXaxis()->SetTitleSize(0.04);
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetXaxis()->SetLabelSize(0.04);
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetTitle("<Noise_{CMN}>");
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleSize(0.06);
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleOffset(1.0);
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetLabelSize(0.04);
    h_mMeanCMNSigma_RStrip[i_rstrip]->GetYaxis()->SetRangeUser(5.0,30.0);
    h_mMeanCMNSigma_RStrip[i_rstrip]->SetLineColor(i_rstrip+1);
    h_mMeanCMNSigma_RStrip[i_rstrip]->SetLineWidth(3);

    if(i_rstrip == 0) h_mMeanCMNSigma_RStrip[i_rstrip]->Draw();
    else h_mMeanCMNSigma_RStrip[i_rstrip]->Draw("same");

    string LegName = Form("R_strip %d",i_rstrip);
    leg_mean->AddEntry(h_mMeanCMNSigma_RStrip[i_rstrip],LegName.c_str(),"L");
  }
  h_mMeanCMNSigma_IST[1]->SetMarkerStyle(20);
  h_mMeanCMNSigma_IST[1]->SetMarkerColor(kGray+2);
  h_mMeanCMNSigma_IST[1]->SetMarkerSize(1.0);
  h_mMeanCMNSigma_IST[1]->Draw("p same");
  leg_mean->AddEntry(h_mMeanCMNSigma_IST[1],"IST 2013","p");

  h_mMeanCMNSigma_ist[0]->SetMarkerStyle(20);
  h_mMeanCMNSigma_ist[0]->SetMarkerColor(2);
  h_mMeanCMNSigma_ist[0]->SetMarkerSize(1.0);
  h_mMeanCMNSigma_ist[0]->Draw("p same");
  leg_mean->AddEntry(h_mMeanCMNSigma_ist[0],"IST Current","p");
  leg_mean->Draw("same");

  c_NoiseMean->cd(2);
  for(int i_rstrip = 0; i_rstrip < FST::numRStrip; ++i_rstrip)
  {
    h_mMeanRanSigma_RStrip[i_rstrip]->SetTitle("FST Differential Noise");
    h_mMeanRanSigma_RStrip[i_rstrip]->SetStats(0);
    h_mMeanRanSigma_RStrip[i_rstrip]->GetXaxis()->SetTitle("Time Bin");
    h_mMeanRanSigma_RStrip[i_rstrip]->GetXaxis()->SetTitleSize(0.04);
    h_mMeanRanSigma_RStrip[i_rstrip]->GetXaxis()->SetLabelSize(0.04);
    h_mMeanRanSigma_RStrip[i_rstrip]->GetYaxis()->SetTitle("<Noise_{Ran}>");
    h_mMeanRanSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleSize(0.06);
    h_mMeanRanSigma_RStrip[i_rstrip]->GetYaxis()->SetTitleOffset(1.0);
    h_mMeanRanSigma_RStrip[i_rstrip]->GetYaxis()->SetLabelSize(0.04);
    h_mMeanRanSigma_RStrip[i_rstrip]->GetYaxis()->SetRangeUser(5.0,30.0);
    h_mMeanRanSigma_RStrip[i_rstrip]->SetLineColor(i_rstrip+1);
    h_mMeanRanSigma_RStrip[i_rstrip]->SetLineWidth(3);

    if(i_rstrip == 0) h_mMeanRanSigma_RStrip[i_rstrip]->Draw();
    else h_mMeanRanSigma_RStrip[i_rstrip]->Draw("same");
  }

  h_mMeanRanSigma_ist[0]->SetMarkerStyle(20);
  h_mMeanRanSigma_ist[0]->SetMarkerColor(2);
  h_mMeanRanSigma_ist[0]->SetMarkerSize(1.0);
  h_mMeanRanSigma_ist[0]->Draw("p same");

  h_mMeanRanSigma_IST[1]->SetMarkerStyle(20);
  h_mMeanRanSigma_IST[1]->SetMarkerColor(kGray+2);
  h_mMeanRanSigma_IST[1]->SetMarkerSize(1.0);
  h_mMeanRanSigma_IST[1]->Draw("p same");

  leg_mean->Draw("same");
  
  c_NoiseMean->SaveAs("./figures/Comparison/SensorPedNoiseComparison.eps");
}
