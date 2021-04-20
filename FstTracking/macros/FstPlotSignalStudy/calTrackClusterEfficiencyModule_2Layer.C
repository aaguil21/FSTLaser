#include <iostream>
#include <string>

#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include "./draw.h"
#include "../../src/FstUtil/FstCons.h"

using namespace std;

void calTrackClusterEfficiencyModule_2Layer(string mod = "Mod04", string hv = "HV140V", bool isSavePed = true, bool isApplyCMNCorr = true, float nFstHitsCut = 4.0, int numOfUsedTimeBins = 2, float nFstThresholdCut2 = 2.5, float nFstThresholdCut1 = 3.5)
{
  std::string pedMode = "withPed";
  if(!isSavePed) pedMode = "woPed";
  std::string cmnMode = "withCMNCorr";
  if(!isApplyCMNCorr) cmnMode = "woCMNCorr";

  string inputfile = Form("../../output/signal/%s/FstTracking_%s_%s_Th%1.1fTb%dPed%1.1fPed%1.1f_%s_%s.root",mod.c_str(),mod.c_str(),hv.c_str(),nFstHitsCut,numOfUsedTimeBins,nFstThresholdCut2,nFstThresholdCut1,pedMode.c_str(),cmnMode.c_str());

  TFile *File_InPut = TFile::Open(inputfile.c_str());

  // simple cluster
  TH2F *h_mSimpleClustersTrackIstCounts_2Layer[FST::mFstNumMatching];
  TH1F *h_mSimpleClustersTrackIstCountsR_2Layer[FST::mFstNumMatching];
  TH1F *h_mSimpleClustersTrackIstCountsPhi_2Layer[FST::mFstNumMatching];

  TH2F *h_mSimpleClustersTrackFstCounts_2Layer[FST::mFstNumMatching];
  TH1F *h_mSimpleClustersTrackFstCountsR_2Layer[FST::mFstNumMatching];
  TH1F *h_mSimpleClustersTrackFstCountsPhi_2Layer[FST::mFstNumMatching];
  for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
  {
    string HistName;
    HistName = Form("h_mSimpleClustersTrackIstCountsModule_2Layer_SF%d",i_match);
    h_mSimpleClustersTrackIstCounts_2Layer[i_match] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mSimpleClustersTrackIstCounts_2Layer[i_match]->Sumw2();

    // IST R projection
    HistName = Form("h_mSimpleClustersTrackIstCountsR_2Layer_SF%d",i_match);
    int projIstBinX0 = h_mSimpleClustersTrackIstCounts_2Layer[i_match]->GetYaxis()->FindBin(FST::mFstPhiMin[1]+4*FST::pitchPhi);
    int projIstBinX1 = h_mSimpleClustersTrackIstCounts_2Layer[i_match]->GetYaxis()->FindBin(FST::mFstPhiMax[0]-4*FST::pitchPhi);
    h_mSimpleClustersTrackIstCountsR_2Layer[i_match] = (TH1F*)h_mSimpleClustersTrackIstCounts_2Layer[i_match]->ProjectionX(HistName.c_str(),projIstBinX0,projIstBinX1);

    // IST Phi projection
    HistName = Form("h_mSimpleClustersTrackIstCountsPhi_2Layer_SF%d",i_match);
    int projIstBinY0 = h_mSimpleClustersTrackIstCounts_2Layer[i_match]->GetXaxis()->FindBin(FST::mFstRMin[1]+0.4*FST::pitchR);
    int projIstBinY1 = h_mSimpleClustersTrackIstCounts_2Layer[i_match]->GetXaxis()->FindBin(FST::mFstRMax[1]-0.4*FST::pitchR);
    h_mSimpleClustersTrackIstCountsPhi_2Layer[i_match] = (TH1F*)h_mSimpleClustersTrackIstCounts_2Layer[i_match]->ProjectionY(HistName.c_str(),projIstBinY0,projIstBinY1);

    HistName = Form("h_mSimpleClustersTrackFstCountsModule_2Layer_SF%d",i_match);
    h_mSimpleClustersTrackFstCounts_2Layer[i_match] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mSimpleClustersTrackFstCounts_2Layer[i_match]->Sumw2();

    // FST R projection
    HistName = Form("h_mSimpleClustersTrackFstCountsR_2Layer_SF%d",i_match);
    int projFstBinX0 = h_mSimpleClustersTrackFstCounts_2Layer[i_match]->GetYaxis()->FindBin(FST::mFstPhiMin[1]+4*FST::pitchPhi);
    int projFstBinX1 = h_mSimpleClustersTrackFstCounts_2Layer[i_match]->GetYaxis()->FindBin(FST::mFstPhiMax[0]-4*FST::pitchPhi);
    h_mSimpleClustersTrackFstCountsR_2Layer[i_match] = (TH1F*)h_mSimpleClustersTrackFstCounts_2Layer[i_match]->ProjectionX(HistName.c_str(),projFstBinX0,projFstBinX1);

    // FST Phi projection
    HistName = Form("h_mSimpleClustersTrackFstCountsPhi_2Layer_SF%d",i_match);
    int projFstBinY0 = h_mSimpleClustersTrackFstCounts_2Layer[i_match]->GetXaxis()->FindBin(FST::mFstRMin[1]+0.4*FST::pitchR);
    int projFstBinY1 = h_mSimpleClustersTrackFstCounts_2Layer[i_match]->GetXaxis()->FindBin(FST::mFstRMax[1]-0.4*FST::pitchR);
    h_mSimpleClustersTrackFstCountsPhi_2Layer[i_match] = (TH1F*)h_mSimpleClustersTrackFstCounts_2Layer[i_match]->ProjectionY(HistName.c_str(),projFstBinY0,projFstBinY1);
  }

  TH2F *h_mSimpleClustersTrackFstEff_2Layer[FST::mFstNumMatching];
  TH1F *h_mSimpleClustersTrackFstEffR_2Layer[FST::mFstNumMatching];
  TH1F *h_mSimpleClustersTrackFstEffPhi_2Layer[FST::mFstNumMatching];
  for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
  {
    string HistName;
    HistName = Form("h_mSimpleClustersTrackFstEff_2Layer_SF%d",i_match);
    h_mSimpleClustersTrackFstEff_2Layer[i_match] = (TH2F*)h_mSimpleClustersTrackIstCounts_2Layer[i_match]->Clone(HistName.c_str());
    h_mSimpleClustersTrackFstEff_2Layer[i_match]->SetTitle(HistName.c_str());
    h_mSimpleClustersTrackFstEff_2Layer[i_match]->Reset();
    if(h_mSimpleClustersTrackIstCounts_2Layer[i_match]->GetEntries() > 0)
    {
      h_mSimpleClustersTrackFstEff_2Layer[i_match]->Divide(h_mSimpleClustersTrackFstCounts_2Layer[i_match],h_mSimpleClustersTrackIstCounts_2Layer[i_match],1,1,"B");
    }

    HistName = Form("h_mSimpleClustersTrackFstEffR_2Layer_SF%d",i_match);
    h_mSimpleClustersTrackFstEffR_2Layer[i_match] = (TH1F*)h_mSimpleClustersTrackIstCountsR_2Layer[i_match]->Clone(HistName.c_str());
    h_mSimpleClustersTrackFstEffR_2Layer[i_match]->SetTitle(HistName.c_str());
    h_mSimpleClustersTrackFstEffR_2Layer[i_match]->Reset();
    if(h_mSimpleClustersTrackIstCountsR_2Layer[i_match]->GetEntries() > 0)
    {
      h_mSimpleClustersTrackFstEffR_2Layer[i_match]->Divide(h_mSimpleClustersTrackFstCountsR_2Layer[i_match],h_mSimpleClustersTrackIstCountsR_2Layer[i_match],1,1,"B");
    }

    HistName = Form("h_mSimpleClustersTrackFstEffPhi_2Layer_SF%d",i_match);
    h_mSimpleClustersTrackFstEffPhi_2Layer[i_match] = (TH1F*)h_mSimpleClustersTrackIstCountsPhi_2Layer[i_match]->Clone(HistName.c_str());
    h_mSimpleClustersTrackFstEffPhi_2Layer[i_match]->SetTitle(HistName.c_str());
    h_mSimpleClustersTrackFstEffPhi_2Layer[i_match]->Reset();
    if(h_mSimpleClustersTrackIstCountsPhi_2Layer[i_match]->GetEntries() > 0)
    {
      h_mSimpleClustersTrackFstEffPhi_2Layer[i_match]->Divide(h_mSimpleClustersTrackFstCountsPhi_2Layer[i_match],h_mSimpleClustersTrackIstCountsPhi_2Layer[i_match],1,1,"B");
    }
  }

  // scan cluster
  TH2F *h_mScanClustersTrackIstCounts_2Layer[FST::mFstNumMatching];
  TH1F *h_mScanClustersTrackIstCountsR_2Layer[FST::mFstNumMatching];
  TH1F *h_mScanClustersTrackIstCountsPhi_2Layer[FST::mFstNumMatching];

  TH2F *h_mScanClustersTrackFstCounts_2Layer[FST::mFstNumMatching];
  TH1F *h_mScanClustersTrackFstCountsR_2Layer[FST::mFstNumMatching];
  TH1F *h_mScanClustersTrackFstCountsPhi_2Layer[FST::mFstNumMatching];
  for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
  {
    string HistName;
    HistName = Form("h_mScanClustersTrackIstCountsModule_2Layer_SF%d",i_match);
    h_mScanClustersTrackIstCounts_2Layer[i_match] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mScanClustersTrackIstCounts_2Layer[i_match]->Sumw2();

    // IST R projection
    HistName = Form("h_mScanClustersTrackIstCountsR_2Layer_SF%d",i_match);
    int projIstBinX0 = h_mScanClustersTrackIstCounts_2Layer[i_match]->GetYaxis()->FindBin(FST::mFstPhiMin[1]+4*FST::pitchPhi);
    int projIstBinX1 = h_mScanClustersTrackIstCounts_2Layer[i_match]->GetYaxis()->FindBin(FST::mFstPhiMax[0]-4*FST::pitchPhi);
    h_mScanClustersTrackIstCountsR_2Layer[i_match] = (TH1F*)h_mScanClustersTrackIstCounts_2Layer[i_match]->ProjectionX(HistName.c_str(),projIstBinX0,projIstBinX1);

    // IST Phi projection
    HistName = Form("h_mScanClustersTrackIstCountsPhi_2Layer_SF%d",i_match);
    int projIstBinY0 = h_mScanClustersTrackIstCounts_2Layer[i_match]->GetXaxis()->FindBin(FST::mFstRMin[1]+0.4*FST::pitchR);
    int projIstBinY1 = h_mScanClustersTrackIstCounts_2Layer[i_match]->GetXaxis()->FindBin(FST::mFstRMax[1]-0.4*FST::pitchR);
    h_mScanClustersTrackIstCountsPhi_2Layer[i_match] = (TH1F*)h_mScanClustersTrackIstCounts_2Layer[i_match]->ProjectionY(HistName.c_str(),projIstBinY0,projIstBinY1);

    HistName = Form("h_mScanClustersTrackFstCountsModule_2Layer_SF%d",i_match);
    h_mScanClustersTrackFstCounts_2Layer[i_match] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mScanClustersTrackFstCounts_2Layer[i_match]->Sumw2();

    // FST R projection
    HistName = Form("h_mScanClustersTrackFstCountsR_2Layer_SF%d",i_match);
    int projFstBinX0 = h_mScanClustersTrackFstCounts_2Layer[i_match]->GetYaxis()->FindBin(FST::mFstPhiMin[1]+4*FST::pitchPhi);
    int projFstBinX1 = h_mScanClustersTrackFstCounts_2Layer[i_match]->GetYaxis()->FindBin(FST::mFstPhiMax[0]-4*FST::pitchPhi);
    h_mScanClustersTrackFstCountsR_2Layer[i_match] = (TH1F*)h_mScanClustersTrackFstCounts_2Layer[i_match]->ProjectionX(HistName.c_str(),projFstBinX0,projFstBinX1);

    // FST Phi projection
    HistName = Form("h_mScanClustersTrackFstCountsPhi_2Layer_SF%d",i_match);
    int projFstBinY0 = h_mScanClustersTrackFstCounts_2Layer[i_match]->GetXaxis()->FindBin(FST::mFstRMin[1]+0.4*FST::pitchR);
    int projFstBinY1 = h_mScanClustersTrackFstCounts_2Layer[i_match]->GetXaxis()->FindBin(FST::mFstRMax[1]-0.4*FST::pitchR);
    h_mScanClustersTrackFstCountsPhi_2Layer[i_match] = (TH1F*)h_mScanClustersTrackFstCounts_2Layer[i_match]->ProjectionY(HistName.c_str(),projFstBinY0,projFstBinY1);
  }

  TH2F *h_mScanClustersTrackFstEff_2Layer[FST::mFstNumMatching];
  TH1F *h_mScanClustersTrackFstEffR_2Layer[FST::mFstNumMatching];
  TH1F *h_mScanClustersTrackFstEffPhi_2Layer[FST::mFstNumMatching];
  for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
  {
    string HistName;
    HistName = Form("h_mScanClustersTrackFstEff_2Layer_SF%d",i_match);
    h_mScanClustersTrackFstEff_2Layer[i_match] = (TH2F*)h_mScanClustersTrackIstCounts_2Layer[i_match]->Clone(HistName.c_str());
    h_mScanClustersTrackFstEff_2Layer[i_match]->SetTitle(HistName.c_str());
    h_mScanClustersTrackFstEff_2Layer[i_match]->Reset();
    if(h_mScanClustersTrackIstCounts_2Layer[i_match]->GetEntries() > 0)
    {
      h_mScanClustersTrackFstEff_2Layer[i_match]->Divide(h_mScanClustersTrackFstCounts_2Layer[i_match],h_mScanClustersTrackIstCounts_2Layer[i_match],1,1,"B");
    }

    HistName = Form("h_mScanClustersTrackFstEffR_2Layer_SF%d",i_match);
    h_mScanClustersTrackFstEffR_2Layer[i_match] = (TH1F*)h_mScanClustersTrackIstCountsR_2Layer[i_match]->Clone(HistName.c_str());
    h_mScanClustersTrackFstEffR_2Layer[i_match]->SetTitle(HistName.c_str());
    h_mScanClustersTrackFstEffR_2Layer[i_match]->Reset();
    if(h_mScanClustersTrackIstCountsR_2Layer[i_match]->GetEntries() > 0)
    {
      h_mScanClustersTrackFstEffR_2Layer[i_match]->Divide(h_mScanClustersTrackFstCountsR_2Layer[i_match],h_mScanClustersTrackIstCountsR_2Layer[i_match],1,1,"B");
    }

    HistName = Form("h_mScanClustersTrackFstEffPhi_2Layer_SF%d",i_match);
    h_mScanClustersTrackFstEffPhi_2Layer[i_match] = (TH1F*)h_mScanClustersTrackIstCountsPhi_2Layer[i_match]->Clone(HistName.c_str());
    h_mScanClustersTrackFstEffPhi_2Layer[i_match]->SetTitle(HistName.c_str());
    h_mScanClustersTrackFstEffPhi_2Layer[i_match]->Reset();
    if(h_mScanClustersTrackIstCountsPhi_2Layer[i_match]->GetEntries() > 0)
    {
      h_mScanClustersTrackFstEffPhi_2Layer[i_match]->Divide(h_mScanClustersTrackFstCountsPhi_2Layer[i_match],h_mScanClustersTrackIstCountsPhi_2Layer[i_match],1,1,"B");
    }
  }

  string outputname = "./figures/ModuleEfficiency_TrackCluster_2Layer.pdf";
  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,900,900);
  c_play->Divide(3,3);
  for(int i_pad = 0; i_pad < 9; ++i_pad)
  {
    c_play->cd(i_pad+1)->SetLeftMargin(0.15);
    c_play->cd(i_pad+1)->SetBottomMargin(0.15);
    c_play->cd(i_pad+1)->SetTicks(1,1);
    c_play->cd(i_pad+1)->SetGrid(0,0);
  }

  string output_start = "./figures/ModuleEfficiency_TrackCluster_2Layer.pdf[";
  c_play->Print(output_start.c_str()); // open pdf file

  for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
  {
    c_play->cd();
    string TitleName = Form("Match %d", i_match);
    TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.085,1.0,TitleName.c_str());
    canvas_title->SetBorderSize(1);
    canvas_title->Draw();
    // projection position from IST
    c_play->cd(1);
    h_mScanClustersTrackIstCounts_2Layer[i_match]->GetXaxis()->SetTitle("r_{proj} (mm)");
    h_mScanClustersTrackIstCounts_2Layer[i_match]->GetYaxis()->SetTitle("phi_{proj} (rad)");
    h_mScanClustersTrackIstCounts_2Layer[i_match]->Draw("colz");
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMax[i_sensor], FST::mFstPhiMin[i_sensor], FST::mFstPhiMin[i_sensor], 1, 2, 2);
      PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMax[i_sensor], FST::mFstPhiMax[i_sensor], FST::mFstPhiMax[i_sensor], 1, 2, 2);
      PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMin[i_sensor], FST::mFstPhiMin[i_sensor], FST::mFstPhiMax[i_sensor], 1, 2, 2);
      PlotLine(FST::mFstRMax[i_sensor], FST::mFstRMax[i_sensor], FST::mFstPhiMin[i_sensor], FST::mFstPhiMax[i_sensor], 1, 2, 2);
    }

    c_play->cd(2);
    h_mScanClustersTrackIstCountsR_2Layer[i_match]->GetXaxis()->SetTitle("r_{proj} (mm)");
    h_mScanClustersTrackIstCountsR_2Layer[i_match]->GetYaxis()->SetRangeUser(0.0,1.1*h_mScanClustersTrackIstCountsR_2Layer[i_match]->GetMaximum());
    h_mScanClustersTrackIstCountsR_2Layer[i_match]->SetLineColor(1);
    h_mScanClustersTrackIstCountsR_2Layer[i_match]->Draw("HIST");
    h_mSimpleClustersTrackIstCountsR_2Layer[i_match]->SetLineColor(1);
    h_mSimpleClustersTrackIstCountsR_2Layer[i_match]->Draw("HIST same");
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule-1; ++i_sensor)
    {
      PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMin[i_sensor], 0.5, h_mScanClustersTrackIstCountsR_2Layer[i_match]->GetMaximum(), 1, 2, 2);
      PlotLine(FST::mFstRMax[i_sensor], FST::mFstRMax[i_sensor], 0.5, h_mScanClustersTrackIstCountsR_2Layer[i_match]->GetMaximum(), 1, 2, 2);
    }

    c_play->cd(3);
    h_mScanClustersTrackIstCountsPhi_2Layer[i_match]->GetXaxis()->SetTitle("phi_{proj} (mm)");
    h_mScanClustersTrackIstCountsPhi_2Layer[i_match]->SetLineColor(2);
    h_mScanClustersTrackIstCountsPhi_2Layer[i_match]->Draw("HIST");
    h_mSimpleClustersTrackIstCountsPhi_2Layer[i_match]->SetLineColor(1);
    h_mSimpleClustersTrackIstCountsPhi_2Layer[i_match]->Draw("HIST same");
    for(int i_sensor = 1; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      PlotLine(FST::mFstPhiMin[i_sensor], FST::mFstPhiMin[i_sensor], 0.5, h_mScanClustersTrackIstCountsPhi_2Layer[i_match]->GetMaximum(), 1, 2, 2);
      PlotLine(FST::mFstPhiMax[i_sensor], FST::mFstPhiMax[i_sensor], 0.5, h_mScanClustersTrackIstCountsPhi_2Layer[i_match]->GetMaximum(), 1, 2, 2);
    }

    // measured position from FST
    c_play->cd(4);
    h_mScanClustersTrackFstCounts_2Layer[i_match]->GetXaxis()->SetTitle("r_{proj} (mm)");
    h_mScanClustersTrackFstCounts_2Layer[i_match]->GetYaxis()->SetTitle("phi_{proj} (rad)");
    h_mScanClustersTrackFstCounts_2Layer[i_match]->Draw("colz");
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMax[i_sensor], FST::mFstPhiMin[i_sensor], FST::mFstPhiMin[i_sensor], 1, 2, 2);
      PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMax[i_sensor], FST::mFstPhiMax[i_sensor], FST::mFstPhiMax[i_sensor], 1, 2, 2);
      PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMin[i_sensor], FST::mFstPhiMin[i_sensor], FST::mFstPhiMax[i_sensor], 1, 2, 2);
      PlotLine(FST::mFstRMax[i_sensor], FST::mFstRMax[i_sensor], FST::mFstPhiMin[i_sensor], FST::mFstPhiMax[i_sensor], 1, 2, 2);
    }

    c_play->cd(5);
    h_mScanClustersTrackFstCountsR_2Layer[i_match]->GetXaxis()->SetTitle("r_{proj} (mm)");
    h_mScanClustersTrackFstCountsR_2Layer[i_match]->GetYaxis()->SetRangeUser(0.0,1.1*h_mScanClustersTrackFstCountsR_2Layer[i_match]->GetMaximum());
    h_mScanClustersTrackFstCountsR_2Layer[i_match]->SetLineColor(2);
    h_mScanClustersTrackFstCountsR_2Layer[i_match]->Draw("HIST");
    h_mSimpleClustersTrackFstCountsR_2Layer[i_match]->SetLineColor(1);
    h_mSimpleClustersTrackFstCountsR_2Layer[i_match]->Draw("HIST same");
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule-1; ++i_sensor)
    {
      PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMin[i_sensor], 0.5, h_mScanClustersTrackFstCountsR_2Layer[i_match]->GetMaximum(), 1, 2, 2);
      PlotLine(FST::mFstRMax[i_sensor], FST::mFstRMax[i_sensor], 0.5, h_mScanClustersTrackFstCountsR_2Layer[i_match]->GetMaximum(), 1, 2, 2);
    }

    c_play->cd(6);
    h_mScanClustersTrackFstCountsPhi_2Layer[i_match]->GetXaxis()->SetTitle("phi_{proj} (rad)");
    h_mScanClustersTrackFstCountsPhi_2Layer[i_match]->SetLineColor(2);
    h_mScanClustersTrackFstCountsPhi_2Layer[i_match]->Draw("HIST");
    h_mSimpleClustersTrackFstCountsPhi_2Layer[i_match]->SetLineColor(1);
    h_mSimpleClustersTrackFstCountsPhi_2Layer[i_match]->Draw("HIST same");
    for(int i_sensor = 1; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      PlotLine(FST::mFstPhiMin[i_sensor], FST::mFstPhiMin[i_sensor], 0.5, h_mScanClustersTrackFstCountsPhi_2Layer[i_match]->GetMaximum(), 1, 2, 2);
      PlotLine(FST::mFstPhiMax[i_sensor], FST::mFstPhiMax[i_sensor], 0.5, h_mScanClustersTrackFstCountsPhi_2Layer[i_match]->GetMaximum(), 1, 2, 2);
    }

    c_play->cd(7);
    h_mScanClustersTrackFstEff_2Layer[i_match]->SetStats(0);
    h_mScanClustersTrackFstEff_2Layer[i_match]->GetXaxis()->SetTitle("r_{proj} (mm)");
    h_mScanClustersTrackFstEff_2Layer[i_match]->GetYaxis()->SetTitle("phi_{proj} (rad)");
    h_mScanClustersTrackFstEff_2Layer[i_match]->Draw("colz");
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMax[i_sensor], FST::mFstPhiMin[i_sensor], FST::mFstPhiMin[i_sensor], 1, 2, 2);
      PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMax[i_sensor], FST::mFstPhiMax[i_sensor], FST::mFstPhiMax[i_sensor], 1, 2, 2);
      PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMin[i_sensor], FST::mFstPhiMin[i_sensor], FST::mFstPhiMax[i_sensor], 1, 2, 2);
      PlotLine(FST::mFstRMax[i_sensor], FST::mFstRMax[i_sensor], FST::mFstPhiMin[i_sensor], FST::mFstPhiMax[i_sensor], 1, 2, 2);
    }

    c_play->cd(8);
    h_mScanClustersTrackFstEffR_2Layer[i_match]->SetStats(0);
    h_mScanClustersTrackFstEffR_2Layer[i_match]->GetXaxis()->SetTitle("r_{proj} (mm)");
    h_mScanClustersTrackFstEffR_2Layer[i_match]->GetYaxis()->SetRangeUser(0.0,1.05);
    h_mScanClustersTrackFstEffR_2Layer[i_match]->SetLineColor(2);
    h_mScanClustersTrackFstEffR_2Layer[i_match]->Draw("pE");
    h_mSimpleClustersTrackFstEffR_2Layer[i_match]->SetLineColor(1);
    h_mSimpleClustersTrackFstEffR_2Layer[i_match]->Draw("pE same");
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule-1; ++i_sensor)
    {
      PlotLine(FST::mFstRMin[i_sensor], FST::mFstRMin[i_sensor], 0.0, 1.05, 1, 2, 2);
      PlotLine(FST::mFstRMax[i_sensor], FST::mFstRMax[i_sensor], 0.0, 1.05, 1, 2, 2);
      PlotLine(FST::rMin[i_sensor], FST::rMax[i_sensor], 0.90, 0.90, 4, 1, 2);
      PlotLine(FST::rMin[i_sensor], FST::rMax[i_sensor], 0.95, 0.95, 4, 1, 2);
    }

    c_play->cd(9);
    h_mScanClustersTrackFstEffPhi_2Layer[i_match]->SetStats(0);
    h_mScanClustersTrackFstEffPhi_2Layer[i_match]->GetXaxis()->SetTitle("phi_{proj} (rad)");
    h_mScanClustersTrackFstEffPhi_2Layer[i_match]->GetYaxis()->SetRangeUser(0.0,1.05);
    h_mScanClustersTrackFstEffPhi_2Layer[i_match]->SetLineColor(2);
    h_mScanClustersTrackFstEffPhi_2Layer[i_match]->Draw("pE");
    h_mSimpleClustersTrackFstEffPhi_2Layer[i_match]->SetLineColor(1);
    h_mSimpleClustersTrackFstEffPhi_2Layer[i_match]->Draw("pE same");
    for(int i_sensor = 1; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      PlotLine(FST::mFstPhiMin[i_sensor], FST::mFstPhiMin[i_sensor], 0.0, 1.05, 1, 2, 2);
      PlotLine(FST::mFstPhiMax[i_sensor], FST::mFstPhiMax[i_sensor], 0.0, 1.05, 1, 2, 2);
      PlotLine(FST::phiMin[i_sensor], FST::phiMax[i_sensor], 0.90, 0.90, 4, 1, 2);
      PlotLine(FST::phiMin[i_sensor], FST::phiMax[i_sensor], 0.95, 0.95, 4, 1, 2);
    }
    if(i_match == 0)
    {
      TLegend *leg = new TLegend(0.3,0.4,0.7,0.6);
      leg->SetBorderSize(0);
      leg->SetFillColor(10);
      leg->AddEntry(h_mSimpleClustersTrackFstEffPhi_2Layer[i_match],"Simple Cluster","L");
      leg->AddEntry(h_mScanClustersTrackFstEffPhi_2Layer[i_match],"Scan Cluster","L");
      leg->Draw("same");
    }

    c_play->Update();
    c_play->Print(outputname.c_str());
  }

  string output_stop = "./figures/ModuleEfficiency_TrackCluster_2Layer.pdf]";
  c_play->Print(output_stop.c_str()); // open pdf file

  string outputfile = Form("../../output/signal/%s/FstModuleEff_%s_Th%1.1fTb%dPed%1.1fPed%1.1f_%s_%s.root",mod.c_str(),hv.c_str(),nFstHitsCut,numOfUsedTimeBins,nFstThresholdCut2,nFstThresholdCut1,pedMode.c_str(),cmnMode.c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
  {
    h_mSimpleClustersTrackFstEff_2Layer[i_match]->Write();
    h_mSimpleClustersTrackFstEffR_2Layer[i_match]->Write();
    h_mSimpleClustersTrackFstEffPhi_2Layer[i_match]->Write();

    h_mScanClustersTrackFstEff_2Layer[i_match]->Write();
    h_mScanClustersTrackFstEffR_2Layer[i_match]->Write();
    h_mScanClustersTrackFstEffPhi_2Layer[i_match]->Write();
  }
  File_OutPut->Close();
}
