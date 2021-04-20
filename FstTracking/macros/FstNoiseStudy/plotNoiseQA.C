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
// #include "./draw.h"
#include "../../src/FstUtil/FstCons.h"

using namespace std;
void plotNoiseQA(string module = "Mod03", string mode = "Ped")
{
  const int numOfHV = 16;
  string HV[numOfHV] = {"HV0V","HV5V","HV10V","HV20V","HV30V","HV40V","HV50V","HV60V","HV70V","HV80V","HV90V","HV100V","HV110V","HV120V","HV130V","HV140V"};
  int hv[numOfHV] = {0,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140};

  TFile *File_InPut[numOfHV];
  TH1F *h_mPedSigma_RStrip0[numOfHV];
  TH1F *h_mPedSigma_RStrip3[numOfHV];
  TH1F *h_mPedSigma_RStrip4[numOfHV];
  TH1F *h_mPedSigma_RStrip7[numOfHV];
  TH1F *h_mRanSigma_RStrip0[numOfHV];
  TH1F *h_mRanSigma_RStrip3[numOfHV];
  TH1F *h_mRanSigma_RStrip4[numOfHV];
  TH1F *h_mRanSigma_RStrip7[numOfHV];
  for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
  {
    string inputfile = Form("../../output/noise/Fst%sNoise_%s_%s.root",mode.c_str(),module.c_str(),HV[i_hv].c_str());
    cout << "Open InPut File: " << inputfile.c_str() << endl;
    File_InPut[i_hv] = TFile::Open(inputfile.c_str());
    h_mPedSigma_RStrip0[i_hv] = (TH1F*)File_InPut[i_hv]->Get("h_mPedSigma_FST_RStrip0_TimeBin0");
    h_mPedSigma_RStrip3[i_hv] = (TH1F*)File_InPut[i_hv]->Get("h_mPedSigma_FST_RStrip3_TimeBin0");
    h_mPedSigma_RStrip4[i_hv] = (TH1F*)File_InPut[i_hv]->Get("h_mPedSigma_FST_RStrip4_TimeBin0");
    h_mPedSigma_RStrip7[i_hv] = (TH1F*)File_InPut[i_hv]->Get("h_mPedSigma_FST_RStrip7_TimeBin0");
    h_mRanSigma_RStrip0[i_hv] = (TH1F*)File_InPut[i_hv]->Get("h_mRanSigma_FST_RStrip0_TimeBin0");
    h_mRanSigma_RStrip3[i_hv] = (TH1F*)File_InPut[i_hv]->Get("h_mRanSigma_FST_RStrip3_TimeBin0");
    h_mRanSigma_RStrip4[i_hv] = (TH1F*)File_InPut[i_hv]->Get("h_mRanSigma_FST_RStrip4_TimeBin0");
    h_mRanSigma_RStrip7[i_hv] = (TH1F*)File_InPut[i_hv]->Get("h_mRanSigma_FST_RStrip7_TimeBin0");
  }

  int apv0_RStrip0[4] = {0,   8, 16, 24}; 
  int apv1_RStrip0[4] = {32, 40, 48, 56}; 
  int apv2_RStrip3[4] = {64, 72, 80, 88}; 
  int apv3_RStrip3[4] = {96,104,112,120}; 
  int apv4_RStrip4[4] = {0,   8, 16, 24}; 
  int apv5_RStrip7[4] = {32, 40, 48, 56}; 
  int apv6_RStrip4[4] = {64, 72, 80, 88}; 
  int apv7_RStrip7[4] = {96,104,112,120}; 

  const int channel = 0;

  TH1F *h_mPedSigma_Apv0Rstrip0 = new TH1F("h_mPedSigma_Apv0Rstrip0","h_mPedSigma_Apv0Rstrip0",200,-20.5,179.5);
  TH1F *h_mPedSigma_Apv3Rstrip3 = new TH1F("h_mPedSigma_Apv3Rstrip3","h_mPedSigma_Apv3Rstrip3",200,-20.5,179.5);
  TH1F *h_mPedSigma_Apv4Rstrip4 = new TH1F("h_mPedSigma_Apv4Rstrip4","h_mPedSigma_Apv4Rstrip4",200,-20.5,179.5);
  TH1F *h_mPedSigma_Apv5Rstrip7 = new TH1F("h_mPedSigma_Apv5Rstrip7","h_mPedSigma_Apv5Rstrip7",200,-20.5,179.5);
  TH1F *h_mRanSigma_Apv0Rstrip0 = new TH1F("h_mRanSigma_Apv0Rstrip0","h_mRanSigma_Apv0Rstrip0",200,-20.5,179.5);
  TH1F *h_mRanSigma_Apv3Rstrip3 = new TH1F("h_mRanSigma_Apv3Rstrip3","h_mRanSigma_Apv3Rstrip3",200,-20.5,179.5);
  TH1F *h_mRanSigma_Apv4Rstrip4 = new TH1F("h_mRanSigma_Apv4Rstrip4","h_mRanSigma_Apv4Rstrip4",200,-20.5,179.5);
  TH1F *h_mRanSigma_Apv5Rstrip7 = new TH1F("h_mRanSigma_Apv5Rstrip7","h_mRanSigma_Apv5Rstrip7",200,-20.5,179.5);

  for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
  {
    // APV 0 & RStrip 0
    int PedBin = h_mPedSigma_RStrip0[i_hv]->FindBin(apv0_RStrip0[channel]);
    double PedNoise = h_mPedSigma_RStrip0[i_hv]->GetBinContent(PedBin);
    PedBin = h_mPedSigma_Apv0Rstrip0->FindBin(hv[i_hv]);
    h_mPedSigma_Apv0Rstrip0->SetBinContent(PedBin,PedNoise);
    // APV 0 & RStrip 0
    int RanBin = h_mRanSigma_RStrip0[i_hv]->FindBin(apv0_RStrip0[channel]);
    double RanNoise = h_mRanSigma_RStrip0[i_hv]->GetBinContent(RanBin);
    RanBin = h_mRanSigma_Apv0Rstrip0->FindBin(hv[i_hv]);
    h_mRanSigma_Apv0Rstrip0->SetBinContent(RanBin,RanNoise);
  }

  const int normalCh = 14;
  TH1F *h_mPedSigma = new TH1F("h_mPedSigma","h_mPedSigma",200,-20.5,179.5);
  TH1F *h_mRanSigma = new TH1F("h_mRanSigma","h_mRanSigma",200,-20.5,179.5);
  for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
  {
    // APV 0 & RStrip 0
    int PedBin = h_mPedSigma_RStrip0[i_hv]->FindBin(normalCh);
    double PedNoise = h_mPedSigma_RStrip0[i_hv]->GetBinContent(PedBin);
    PedBin = h_mPedSigma->FindBin(hv[i_hv]);
    h_mPedSigma->SetBinContent(PedBin,PedNoise);
    // APV 0 & RStrip 0
    int RanBin = h_mRanSigma_RStrip0[i_hv]->FindBin(normalCh);
    double RanNoise = h_mRanSigma_RStrip0[i_hv]->GetBinContent(RanBin);
    RanBin = h_mRanSigma->FindBin(hv[i_hv]);
    h_mRanSigma->SetBinContent(RanBin,RanNoise);
  }

  TCanvas *c_noise = new TCanvas("c_noise","c_noise",10,10,800,400);
  c_noise->Divide(2,1);
  for(int i_pad = 0; i_pad < 1; ++i_pad)
  {
    c_noise->cd(i_pad+1)->SetLeftMargin(0.15);
    c_noise->cd(i_pad+1)->SetRightMargin(0.15);
    c_noise->cd(i_pad+1)->SetBottomMargin(0.15);
    c_noise->cd(i_pad+1)->SetTicks(1,1);
    c_noise->cd(i_pad+1)->SetGrid(0,0);
  }

  c_noise->cd(1);
  h_mPedSigma_Apv0Rstrip0->SetTitle("Total Noise");
  h_mPedSigma_Apv0Rstrip0->SetMarkerSize(1.2);
  h_mPedSigma_Apv0Rstrip0->SetMarkerStyle(24);
  h_mPedSigma_Apv0Rstrip0->SetMarkerColor(2);
  h_mPedSigma_Apv0Rstrip0->GetXaxis()->SetRangeUser(-9.5,150);
  h_mPedSigma_Apv0Rstrip0->GetYaxis()->SetRangeUser(5,60);
  h_mPedSigma_Apv0Rstrip0->Draw("P");
  h_mPedSigma->SetMarkerSize(1.2);
  h_mPedSigma->SetMarkerStyle(20);
  h_mPedSigma->SetMarkerColor(kGray);
  h_mPedSigma->Draw("P same");

  TLegend *leg = new TLegend(0.3,0.6,0.8,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(h_mPedSigma_Apv0Rstrip0,"Suspicious Channel","P");
  leg->AddEntry(h_mPedSigma,"Normal Channel","P");
  leg->Draw("same");

  c_noise->cd(2);
  h_mRanSigma_Apv0Rstrip0->SetTitle("Random Noise");
  h_mRanSigma_Apv0Rstrip0->SetMarkerSize(1.2);
  h_mRanSigma_Apv0Rstrip0->SetMarkerStyle(24);
  h_mRanSigma_Apv0Rstrip0->SetMarkerColor(2);
  h_mRanSigma_Apv0Rstrip0->GetXaxis()->SetRangeUser(-9.5,150);
  h_mRanSigma_Apv0Rstrip0->GetYaxis()->SetRangeUser(5,60);
  h_mRanSigma_Apv0Rstrip0->Draw("P");
  h_mRanSigma->SetMarkerSize(1.2);
  h_mRanSigma->SetMarkerStyle(20);
  h_mRanSigma->SetMarkerColor(kGray);
  h_mRanSigma->Draw("P same");
  leg->Draw("same");

  c_noise->SaveAs("./NoiseScan.pdf");

  /*
  TH1F *h_mPedSigma_Apv0Rstrip0[4];
  TH1F *h_mPedSigma_Apv1Rstrip0[4];
  TH1F *h_mPedSigma_Apv2Rstrip3[4];
  TH1F *h_mPedSigma_Apv3Rstrip3[4];
  TH1F *h_mPedSigma_Apv4Rstrip4[4];
  TH1F *h_mPedSigma_Apv5Rstrip7[4];
  TH1F *h_mPedSigma_Apv6Rstrip4[4];
  TH1F *h_mPedSigma_Apv7Rstrip7[4];
  for(int i_ch = 0; i_ch < 4; ++i_ch)
  {
    string HistName = Form("h_mPedSigma_Apv0Rstrip0_ch%d",i_ch);
    h_mPedSigma_Apv0Rstrip0[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mPedSigma_Apv1Rstrip0_ch%d",i_ch);
    h_mPedSigma_Apv1Rstrip0[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mPedSigma_Apv2Rstrip3_ch%d",i_ch);
    h_mPedSigma_Apv2Rstrip3[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mPedSigma_Apv3Rstrip3_ch%d",i_ch);
    h_mPedSigma_Apv3Rstrip3[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mPedSigma_Apv4Rstrip4_ch%d",i_ch);
    h_mPedSigma_Apv4Rstrip4[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mPedSigma_Apv5Rstrip7_ch%d",i_ch);
    h_mPedSigma_Apv5Rstrip7[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mPedSigma_Apv6Rstrip4_ch%d",i_ch);
    h_mPedSigma_Apv6Rstrip4[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mPedSigma_Apv7Rstrip7_ch%d",i_ch);
    h_mPedSigma_Apv7Rstrip7[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
  }

  TH1F *h_mRanSigma_Apv0Rstrip0[4];
  TH1F *h_mRanSigma_Apv1Rstrip0[4];
  TH1F *h_mRanSigma_Apv2Rstrip3[4];
  TH1F *h_mRanSigma_Apv3Rstrip3[4];
  TH1F *h_mRanSigma_Apv4Rstrip4[4];
  TH1F *h_mRanSigma_Apv5Rstrip7[4];
  TH1F *h_mRanSigma_Apv6Rstrip4[4];
  TH1F *h_mRanSigma_Apv7Rstrip7[4];
  for(int i_ch = 0; i_ch < 4; ++i_ch)
  {
    string HistName = Form("h_mRanSigma_Apv0Rstrip0_ch%d",i_ch);
    h_mRanSigma_Apv0Rstrip0[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mRanSigma_Apv1Rstrip0_ch%d",i_ch);
    h_mRanSigma_Apv1Rstrip0[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mRanSigma_Apv2Rstrip3_ch%d",i_ch);
    h_mRanSigma_Apv2Rstrip3[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mRanSigma_Apv3Rstrip3_ch%d",i_ch);
    h_mRanSigma_Apv3Rstrip3[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mRanSigma_Apv4Rstrip4_ch%d",i_ch);
    h_mRanSigma_Apv4Rstrip4[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mRanSigma_Apv5Rstrip7_ch%d",i_ch);
    h_mRanSigma_Apv5Rstrip7[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mRanSigma_Apv6Rstrip4_ch%d",i_ch);
    h_mRanSigma_Apv6Rstrip4[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
    string HistName = Form("h_mRanSigma_Apv7Rstrip7_ch%d",i_ch);
    h_mRanSigma_Apv7Rstrip7[i_ch] = new TH1F(HistName.c_str(),HistName.c_str(),200,-0.5,199.5);
  }
  */
}
