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
#include "../FstPlotConfiguration/draw.h"
#include "../../src/FstUtil/FstCons.h"

using namespace std;

void plotMcEffComparison(string mod = "Mod03", string hv = "HV70V", int sensorId = 0)
{
  // read in Data efficiency 
  string inputData = Form("../../output/configuration/FstEfficiency_%s_%s_Th4.0Tb2Ped2.5Ped3.5_withPed_withCMNCorr.root",mod.c_str(),hv.c_str());
  TFile *File_InPutData = TFile::Open(inputData.c_str());
  TH1F *h_mSensorScanClustersTrackFstEffR_2Layer[FST::mFstNumSensorsPerModule][FST::mFstNumMatching];
  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
    {
      string HistName;
      HistName = Form("h_mScanClustersTrackFstEffR_2Layer_Sensor%d_SF%d",i_sensor,i_match);
      h_mSensorScanClustersTrackFstEffR_2Layer[i_sensor][i_match] = (TH1F*)File_InPutData->Get(HistName.c_str());
    }
  }

  // read in MC efficiency
  string inputMc = Form("../../output/simulation/FstMcEfficiency_Sensor%d.root",sensorId);
  TFile *File_InPutMc = TFile::Open(inputMc.c_str());
  TH1F *h_mFstMcEffR_2Layer[FST::mFstNumMatching];
  for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
  {
    string HistName;
    HistName = Form("h_mFstEffR_2Layer_SF%d",i_match);
    h_mFstMcEffR_2Layer[i_match] = (TH1F*)File_InPutMc->Get(HistName.c_str());
  }

  const int defMatch = 0;

  TCanvas *c_McEffComparison = new TCanvas("c_McEffComparison","c_McEffComparison",10,10,800,800);
  c_McEffComparison->cd()->SetLeftMargin(0.15);
  c_McEffComparison->cd()->SetBottomMargin(0.15);
  c_McEffComparison->cd()->SetTicks(1,1);
  c_McEffComparison->cd()->SetGrid(0,0);

  h_mFstMcEffR_2Layer[defMatch]->SetStats(0);
  h_mFstMcEffR_2Layer[defMatch]->GetXaxis()->SetTitle("r_{proj} (mm)");
  h_mFstMcEffR_2Layer[defMatch]->GetYaxis()->SetRangeUser(0.0,1.05);
  h_mFstMcEffR_2Layer[defMatch]->SetLineColor(1);
  h_mFstMcEffR_2Layer[defMatch]->Draw("HIST");

  h_mSensorScanClustersTrackFstEffR_2Layer[sensorId][defMatch]->SetMarkerStyle(24);
  h_mSensorScanClustersTrackFstEffR_2Layer[sensorId][defMatch]->SetMarkerSize(1.4);
  h_mSensorScanClustersTrackFstEffR_2Layer[sensorId][defMatch]->SetMarkerColor(1);
  h_mSensorScanClustersTrackFstEffR_2Layer[sensorId][defMatch]->SetLineColor(1);
  h_mSensorScanClustersTrackFstEffR_2Layer[sensorId][defMatch]->Draw("pE Same");

  PlotLine(FST::mFstRMin[sensorId], FST::mFstRMin[sensorId], 0.0, 1.05, 1, 2, 2);
  PlotLine(FST::mFstRMax[sensorId], FST::mFstRMax[sensorId], 0.0, 1.05, 1, 2, 2);
  PlotLine(FST::rMin[sensorId], FST::rMax[sensorId], 0.90, 0.90, 4, 1, 2);
  PlotLine(FST::rMin[sensorId], FST::rMax[sensorId], 0.95, 0.95, 4, 1, 2);

  TLegend *leg = new TLegend(0.40,0.2,0.55,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(h_mSensorScanClustersTrackFstEffR_2Layer[sensorId][defMatch],"Data","p");
  leg->AddEntry(h_mFstMcEffR_2Layer[defMatch],"MC","l");
  leg->Draw("same");

  string FigName = Form("./figures/c_McEffComparison_Sensor%d.png",sensorId);
  c_McEffComparison->SaveAs(FigName.c_str());
}
