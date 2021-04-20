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

void plotNumOfRawHits(string mod = "Mod04", string hv = "HV70V", bool isSavePed = true, bool isApplyCMNCorr = true, float nFstHitsCut = 3.5, int numOfUsedTimeBins = 2, float nFstThresholdCut2 = 2.5, float nFstThresholdCut1 = 3.5)
{
  std::string pedMode = "withPed";
  if(!isSavePed) pedMode = "woPed";
  std::string cmnMode = "withCMNCorr";
  if(!isApplyCMNCorr) cmnMode = "woCMNCorr";

  string inputfile = Form("../../output/signal/%s/FstQAStudy_%s_%s_Th%1.1fTb%dPed%1.1fPed%1.1f_%s_%s.root",mod.c_str(),mod.c_str(),hv.c_str(),nFstHitsCut,numOfUsedTimeBins,nFstThresholdCut2,nFstThresholdCut1,pedMode.c_str(),cmnMode.c_str());

  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mNumFstSimpleClusters = (TH1F*)File_InPut->Get("h_mNumFstSimpleClusters");
  TH1F *h_mNumFstRawHitsSimpleCluster = (TH1F*)File_InPut->Get("h_mNumFstRawHitsSimpleCluster");
  TH1F *h_mNumFstScanClusters = (TH1F*)File_InPut->Get("h_mNumFstScanClusters");
  TH1F *h_mNumFstRawHitsScanCluster = (TH1F*)File_InPut->Get("h_mNumFstRawHitsScanCluster");

  TCanvas *c_EventDisplay = new TCanvas("c_EventDisplay","c_EventDisplay",10,10,800,400);
  c_EventDisplay->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_EventDisplay->cd(i_pad+1)->SetLeftMargin(0.15);
    c_EventDisplay->cd(i_pad+1)->SetRightMargin(0.15);
    c_EventDisplay->cd(i_pad+1)->SetBottomMargin(0.15);
    c_EventDisplay->cd(i_pad+1)->SetTicks(1,1);
    c_EventDisplay->cd(i_pad+1)->SetGrid(0,0);
    c_EventDisplay->cd(i_pad+1)->SetLogy();
  }

  c_EventDisplay->cd(1);
  h_mNumFstSimpleClusters->SetTitle("Number of Clusters");
  h_mNumFstSimpleClusters->GetXaxis()->SetTitle("# Clusters");
  h_mNumFstSimpleClusters->GetXaxis()->SetTitleSize(0.06);
  h_mNumFstSimpleClusters->GetXaxis()->SetRangeUser(-0.5,44.5);
  h_mNumFstSimpleClusters->SetLineColor(1);
  h_mNumFstSimpleClusters->Draw();
  h_mNumFstScanClusters->SetLineColor(2);
  h_mNumFstScanClusters->Draw("same");

  TLegend *leg = new TLegend(0.4,0.6,0.85,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(h_mNumFstSimpleClusters,"Simple Cluster","L");
  leg->AddEntry(h_mNumFstScanClusters,"Scan Cluster","L");
  leg->Draw("same");

  c_EventDisplay->cd(2);
  h_mNumFstRawHitsSimpleCluster->SetTitle("Number of RawHits per Cluster");
  h_mNumFstRawHitsSimpleCluster->GetXaxis()->SetTitle("# RawHits per Cluster");
  h_mNumFstRawHitsSimpleCluster->GetXaxis()->SetTitleSize(0.06);
  h_mNumFstRawHitsSimpleCluster->GetXaxis()->SetRangeUser(-0.5,44.5);
  h_mNumFstRawHitsSimpleCluster->SetLineColor(1);
  h_mNumFstRawHitsSimpleCluster->Draw();
  h_mNumFstRawHitsScanCluster->SetLineColor(2);
  h_mNumFstRawHitsScanCluster->Draw("same");

  c_EventDisplay->SaveAs("./figures/NumOfRawHits.pdf");
}

