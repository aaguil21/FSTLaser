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
#include <TPaveStats.h>
#include "./draw.h"
#include "../../src/FstUtil/FstCons.h"

using namespace std;

double gaussian(double *var, double *par)
{
  double residual = var[0];

  double norm  = par[0];
  double mean  = par[1];
  double sigma = par[2];
  double binwidth = par[3];

  double sigmaSquare = sigma*sigma;
  double Norm = norm/(sqrt(2.0*TMath::Pi())*sigma);
  double Power = -0.5*(residual-mean)*(residual-mean)/sigmaSquare;

  double y = Norm*exp(Power)*binwidth;

  return y;
}

void plotResidual_FSTClusterTracksHV(string mod = "Mod03")
{
  bool isSavePed = true;
  bool isApplyCMNCorr = true;
  float nFstHitsCut = 4.0;
  int numOfUsedTimeBins = 2;
  float nFstThresholdCut2 = 2.5; 
  float nFstThresholdCut1 = 3.5;

  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1001);
  gStyle->SetStatX(0.95); gStyle->SetStatY(0.90);
  gStyle->SetStatW(0.20); gStyle->SetStatH(0.20);

  std::string pedMode = "withPed";
  if(!isSavePed) pedMode = "woPed";
  std::string cmnMode = "withCMNCorr";
  if(!isApplyCMNCorr) cmnMode = "woCMNCorr";

  const int numOfHV = 4;
  string hv[numOfHV] = {"HV70V","HV100V","HV120V","HV140V"};
  double highVolt[numOfHV] = {70.0,100.0,120.0,140.0};
  TFile *File_InPut[numOfHV];

  // FST scan clusters
  TH1F *h_mScanClustersTrackFstResR_2Layer[numOfHV][FST::mFstNumSensorsPerModule];
  TH1F *h_mScanClustersTrackFstResPhi_2Layer[numOfHV][FST::mFstNumSensorsPerModule];
  TH1F *h_mScanClustersTrackFstResR_2Layer_Rstrips[numOfHV][FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];
  TH1F *h_mScanClustersTrackFstResPhi_2Layer_Rstrips[numOfHV][FST::mFstNumSensorsPerModule][FST::mFstNumRstripPerSensor];

  for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
  {
    string inputfile = Form("../../output/configuration/FstTracking_%s_%s_Th%1.1fTb%dPed%1.1fPed%1.1f_%s_%s.root",mod.c_str(),hv[i_hv].c_str(),nFstHitsCut,numOfUsedTimeBins,nFstThresholdCut2,nFstThresholdCut1,pedMode.c_str(),cmnMode.c_str());
    File_InPut[i_hv] = TFile::Open(inputfile.c_str());
    for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
    {
      for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
      {
	string HistName;
	HistName = Form("h_mScanClustersTrackFstResR_2Layer_Sensor%d",i_sensor);
	h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());
	HistName = Form("h_mScanClustersTrackFstResPhi_2Layer_Sensor%d",i_sensor);
	h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());

	for(int i_rstrip = 0; i_rstrip < FST::mFstNumRstripPerSensor; ++i_rstrip)
	{
	  HistName = Form("h_mScanClustersTrackFstResR_2Layer_Sensor%d_Rstrip%d",i_sensor,i_rstrip); // R Residual
	  if(i_sensor > 0) HistName = Form("h_mScanClustersTrackFstResR_2Layer_Sensor%d_Rstrip%d",i_sensor,i_rstrip+4);
	  h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());
	  HistName = Form("h_mScanClustersTrackFstResPhi_2Layer_Sensor%d_Rstrip%d",i_sensor,i_rstrip); // Phi Residual
	  if(i_sensor > 0) HistName = Form("h_mScanClustersTrackFstResPhi_2Layer_Sensor%d_Rstrip%d",i_sensor,i_rstrip+4);
	  h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip] = (TH1F*)File_InPut[i_hv]->Get(HistName.c_str());
	}
      }
    }
  }

  string outputname = Form("./figures/Residual_FSTClusterTracksHV_%s.pdf",mod.c_str());
  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,400);
  c_play->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_play->cd(i_pad+1)->SetLeftMargin(0.15);
    c_play->cd(i_pad+1)->SetRightMargin(0.15);
    c_play->cd(i_pad+1)->SetBottomMargin(0.15);
    c_play->cd(i_pad+1)->SetTicks(1,1);
    c_play->cd(i_pad+1)->SetGrid(0,0);
  }

  string output_start = Form("./figures/Residual_FSTClusterTracksHV_%s.pdf[",mod.c_str());
  c_play->Print(output_start.c_str()); // open pdf file

  for(int i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; ++i_sensor)
  {
    c_play->cd();
    string TitleName = Form("Sensor %d",i_sensor);
    TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.055,1.0,TitleName.c_str());
    canvas_title->SetBorderSize(1);
    canvas_title->Draw();
    for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
    {
      c_play->cd(1);
      h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->SetTitle("Corrected R-residual FST Scan Clusters");
      h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->GetXaxis()->SetTitle("r-residual (mm)");
      h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->GetXaxis()->SetTitleSize(0.06);
      h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->GetYaxis()->SetTitle("No. Tracks");
      h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->GetYaxis()->SetTitleSize(0.06);
      if(i_hv == 0) 
      {
	h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->SetMarkerStyle(24);
	h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->SetMarkerColor(1);
	h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->SetMarkerSize(1.0);
	h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->Draw("pE");
      }
      else 
      {
	double maxNorm   = h_mScanClustersTrackFstResR_2Layer[0][i_sensor]->GetMaximum();
	double maxRstrip = h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->GetMaximum();
	h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->Scale(maxNorm/maxRstrip);
	h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->SetLineColor(i_hv+1);
	h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->SetLineStyle(1);
	h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->SetLineWidth(2);
	h_mScanClustersTrackFstResR_2Layer[i_hv][i_sensor]->Draw("HIST same");
      }
    }
    TLegend *leg = new TLegend(0.25,0.55,0.45,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    // leg->SetMargin(0.1);
    leg->AddEntry(h_mScanClustersTrackFstResR_2Layer[0][i_sensor],hv[0].c_str(),"p");
    leg->AddEntry(h_mScanClustersTrackFstResR_2Layer[1][i_sensor],hv[1].c_str(),"l");
    leg->AddEntry(h_mScanClustersTrackFstResR_2Layer[2][i_sensor],hv[2].c_str(),"l");
    leg->AddEntry(h_mScanClustersTrackFstResR_2Layer[3][i_sensor],hv[3].c_str(),"l");
    leg->Draw("same");

    for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
    {
      c_play->cd(2);
      h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->SetTitle("Corrected #phi-residual FST Scan Clusters");
      h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->GetXaxis()->SetTitle("#phi-residual (rad)");
      h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->GetXaxis()->SetTitleSize(0.06);
      h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->GetYaxis()->SetTitle("No. Tracks");
      h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->GetYaxis()->SetTitleSize(0.06);
      if(i_hv == 0) 
      {
	h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->SetMarkerStyle(24);
	h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->SetMarkerColor(1);
	h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->SetMarkerSize(1.0);
	h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->GetXaxis()->SetRangeUser(-0.05,0.05);
	h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->Draw("PE");
      }
      else 
      {
	double maxNorm   = h_mScanClustersTrackFstResPhi_2Layer[0][i_sensor]->GetMaximum();
	double maxRstrip = h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->GetMaximum();
	h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->Scale(maxNorm/maxRstrip);
	h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->SetLineColor(i_hv+1);
	h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->SetLineStyle(1);
	h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->SetLineWidth(2);
	h_mScanClustersTrackFstResPhi_2Layer[i_hv][i_sensor]->Draw("HIST same");
      }
    }
    c_play->Update();
    c_play->Print(outputname.c_str());

    for(int i_rstrip = 0; i_rstrip < FST::mFstNumRstripPerSensor; ++i_rstrip)
    { 
      c_play->cd();
      string TitleName = Form("Sensor %d & R Strip %d",i_sensor, i_rstrip);
      TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.085,1.0,TitleName.c_str());
      canvas_title->SetBorderSize(1);
      canvas_title->Draw();
      string HistName;

      c_play->cd(1)->Clear();
      c_play->cd(1);
      for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
      {
	HistName = Form("Corrected R-residual FST Scan Clusters Sensor %d Rstrip%d",i_sensor, i_rstrip);
	h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetTitle(HistName.c_str());
	h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->GetXaxis()->SetTitle("r-residual (mm)");
	h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->GetXaxis()->SetTitleSize(0.06);
	h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->GetYaxis()->SetTitle("No. Tracks");
	h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->GetYaxis()->SetTitleSize(0.06);
	if(i_hv == 0) 
	{
	  h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetMarkerStyle(24);
	  h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetMarkerColor(1);
	  h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetMarkerSize(1.0);
	  h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->Draw("pE");
	}
	else 
	{
	  double maxNorm   = h_mScanClustersTrackFstResR_2Layer_Rstrips[0][i_sensor][i_rstrip]->GetMaximum();
	  double maxRstrip = h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->GetMaximum();
	  h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->Scale(maxNorm/maxRstrip);
	  h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetLineColor(i_hv+1);
	  h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetLineStyle(1);
	  h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetLineWidth(2);
	  h_mScanClustersTrackFstResR_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->Draw("HIST same");
	}
      }
      leg->Draw("same");

      c_play->cd(2)->Clear();
      c_play->cd(2);
      for(int i_hv = 0; i_hv < numOfHV; ++i_hv)
      {
	HistName = Form("Corrected #phi-residual FST Scan Clusters Sensor %d Rstrip%d",i_sensor, i_rstrip);
	h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetTitle(HistName.c_str());
	h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->GetXaxis()->SetTitle("#phi-residual (rad)");
	h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->GetXaxis()->SetTitleSize(0.06);
	h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->GetYaxis()->SetTitle("No. Tracks");
	h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->GetYaxis()->SetTitleSize(0.06);
	if(i_hv == 0) 
	{
	  h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetMarkerStyle(24);
	  h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetMarkerColor(1);
	  h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetMarkerSize(1.0);
	  h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->GetXaxis()->SetRangeUser(-0.05,0.05);
	  h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->Draw("pE");
	}
	else 
	{
	  double maxNorm   = h_mScanClustersTrackFstResPhi_2Layer_Rstrips[0][i_sensor][i_rstrip]->GetMaximum();
	  double maxRstrip = h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->GetMaximum();
	  h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->Scale(maxNorm/maxRstrip);
	  h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetLineColor(i_hv+1);
	  h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetLineStyle(1);
	  h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->SetLineWidth(2);
	  h_mScanClustersTrackFstResPhi_2Layer_Rstrips[i_hv][i_sensor][i_rstrip]->Draw("HIST same");
	}
      }
      c_play->Update();
      c_play->Print(outputname.c_str());
    }
  }

  string output_stop = Form("./figures/Residual_FSTClusterTracksHV_%s.pdf]",mod.c_str()); 
  c_play->Print(output_stop.c_str()); // close pdf file
}
