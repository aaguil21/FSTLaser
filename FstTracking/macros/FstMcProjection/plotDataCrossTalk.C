#include <iostream>
#include <fstream>
#include <string>

#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>

#include "../../src/FstUtil/FstCons.h"
#include "../FstPlotConfiguration/draw.h"

using namespace std;

void plotDataCrossTalk(int sensorId = 1, string mod = "Mod03", string hv = "HV70V")
{
  string mode = "Scan";
  // string mode = "Simple";


  string inputfile = Form("/Users/xusun/WorkSpace/STAR/ForwardSiliconTracker/FstTracking/output/configuration/FstTracking_%s_%s_Th4.0Tb2Ped2.5Ped3.5_withPed_withCMNCorr.root",mod.c_str(),hv.c_str());

  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH1F *h_mClustersTrackFstResR_2Layer_Rstrips[FST::mFstNumRstripPerSensor];
  for(int i_rstrp = 0; i_rstrp < FST::mFstNumRstripPerSensor; ++i_rstrp)
  {
    string HistName = Form("h_m%sClustersTrackFstResR_2Layer_Sensor%d_Rstrip%d",mode.c_str(),sensorId,i_rstrp); // R Residual
    if(sensorId > 0) HistName = Form("h_m%sClustersTrackFstResR_2Layer_Sensor%d_Rstrip%d",mode.c_str(),sensorId,i_rstrp+4);
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp] = (TH1F*)File_InPut->Get(HistName.c_str());
  }

  TH1F *h_mFstRecoResR_2Layer_Rstrips[4];
  for(int i_rstrp = 0; i_rstrp < 4; ++i_rstrp)
  {
    h_mFstRecoResR_2Layer_Rstrips[i_rstrp] = (TH1F*)h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->Clone();
    h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->Reset();
  }

  TCanvas *c_Residual = new TCanvas("c_Residual","c_Residual",10,10,1600,400);
  c_Residual->Divide(4,1);
  for(int i_pad = 0; i_pad < 4; ++ i_pad)
  {
    c_Residual->cd(i_pad+1);
    c_Residual->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Residual->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Residual->cd(i_pad+1)->SetTicks(1,1);
    c_Residual->cd(i_pad+1)->SetGrid(0,0);
  }

  // fill the reco histogram
  for(int i_track = 0; i_track < 10000; ++i_track)
  {
    for(int i_rstrp = 0; i_rstrp < 4; ++i_rstrp)
    {
      h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->Fill(h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->GetRandom());
    }
  }


  double inteR[4];
  double ctRate[4][7];
  for(int i_rstrp = 0; i_rstrp < 4; ++i_rstrp)
  {
    c_Residual->cd(i_rstrp+1);
    string title = Form("R Residual @ Rstrip%d",i_rstrp);
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->SetTitle(title.c_str());
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->GetXaxis()->SetTitle("r_{ro} - r_{reco} (mm)");
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->GetXaxis()->SetTitleSize(0.06);
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->GetXaxis()->CenterTitle();
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->GetYaxis()->SetTitle("Counts");
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->GetYaxis()->SetTitleSize(0.06);
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->GetYaxis()->CenterTitle();
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->SetLineColor(1);
    h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->Draw("hE");
    int binR_start = h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->FindBin(-150.0);
    int binR_stop  = h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->FindBin(150.0);
    inteR[i_rstrp] = h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->Integral(binR_start,binR_stop);
    double temptRate = 0.0;
    // cout << "i_rstrp = " << i_rstrp << ", binR_start = " << binR_start << ", binR_stop = " << binR_stop << endl;
    cout << "i_rstrp = " << i_rstrp <<  ", inteR = " << inteR[i_rstrp] << endl;
    double offset[7] = {0.0,-0.1,0.0,-0.1,0.0,-0.1,0.0};
    for(int i_delta = 0; i_delta < 7; ++i_delta)
    {
      double deltaStart = (i_delta-3.5)*FST::pitchR;
      double deltaStop  = (i_delta-2.5)*FST::pitchR;
      int binDelta_start = h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->FindBin(deltaStart);
      int binDelta_stop  = h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->FindBin(deltaStop);
      double inteDelta = h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->Integral(binDelta_start,binDelta_stop-1);
      double ratio = inteDelta/inteR[i_rstrp];
      // cout << "i_delta = " << i_delta << ", deltaStart = " << deltaStart << ", deltaStop = " << deltaStop << endl;
      // cout << "binDelta_start = " << binDelta_start << ", binDelta_stop = " << binDelta_stop << endl;
      cout << "i_delta = " << i_delta << ", inteDelta = " << inteDelta << ", ratio = " << inteDelta/inteR[i_rstrp] << endl;
      temptRate += ratio;
      ctRate[i_rstrp][i_delta] = temptRate;
      // cout << "temptRate = " << temptRate << endl;
      PlotLine(deltaStart, deltaStart, 0.5, 0.5*h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->GetMaximum(), 1, 2, 2);
      string percent = Form("%1.1f%%",ratio*100.0);
      double max = h_mClustersTrackFstResR_2Layer_Rstrips[i_rstrp]->GetMaximum();
      plotTopLegend((char*)percent.c_str(),deltaStart,(0.6+offset[i_delta])*max,0.04,1,0.0,42,0,1);
    }
    cout << endl;
  }

  c_Residual->SaveAs("./figures/c_CrossTalkRate.eps");

  // save ctRate to a text file
  string outputfile = "./crosstalkRate.txt";
  ofstream file_crosstalkRate;
  file_crosstalkRate.open(outputfile.c_str());

  if (!file_crosstalkRate.is_open())
  {
    std::cout << "failed to open " << outputfile.c_str() << endl;
    return -1;
  }
  else
  {

    for(int i_rstrp = 0; i_rstrp < 4; ++i_rstrp)
    {
      file_crosstalkRate << std::fixed << std::setprecision(4) << i_rstrp << "    " << ctRate[i_rstrp][0] << "    " << ctRate[i_rstrp][1] << "    " << ctRate[i_rstrp][2] << "    " << ctRate[i_rstrp][3] << "    " << ctRate[i_rstrp][4] << "    " << ctRate[i_rstrp][5] << "    " << ctRate[i_rstrp][6] << endl;
    }
  }

  file_crosstalkRate.close();


  TCanvas *c_Reco = new TCanvas("c_Reco","c_Reco",10,500,1600,400);
  c_Reco->Divide(4,1);
  for(int i_pad = 0; i_pad < 4; ++ i_pad)
  {
    c_Reco->cd(i_pad+1);
    c_Reco->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Reco->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Reco->cd(i_pad+1)->SetTicks(1,1);
    c_Reco->cd(i_pad+1)->SetGrid(0,0);
  }

  for(int i_rstrp = 0; i_rstrp < 4; ++i_rstrp)
  {
    c_Reco->cd(i_rstrp+1);
    string title = Form("R Reco @ Rstrip%d",i_rstrp);
    h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->SetTitle(title.c_str());
    h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetXaxis()->SetTitle("r_{ro} - r_{reco} (mm)");
    h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetXaxis()->SetTitleSize(0.06);
    h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetXaxis()->CenterTitle();
    h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetYaxis()->SetTitle("Counts");
    h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetYaxis()->SetTitleSize(0.06);
    h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetYaxis()->CenterTitle();
    h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->SetLineColor(1);
    h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->Draw("hE");
    PlotLine(0.5*FST::pitchR, 0.5*FST::pitchR, 0.5, h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetMaximum(), 1, 2, 2);
    PlotLine(1.5*FST::pitchR, 1.5*FST::pitchR, 0.5, h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetMaximum(), 1, 2, 2);
    PlotLine(2.5*FST::pitchR, 2.5*FST::pitchR, 0.5, h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetMaximum(), 1, 2, 2);
    PlotLine(3.5*FST::pitchR, 3.5*FST::pitchR, 0.5, h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetMaximum(), 1, 2, 2);
    PlotLine(-0.5*FST::pitchR, -0.5*FST::pitchR, 0.5, h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetMaximum(), 1, 2, 2);
    PlotLine(-1.5*FST::pitchR, -1.5*FST::pitchR, 0.5, h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetMaximum(), 1, 2, 2);
    PlotLine(-2.5*FST::pitchR, -2.5*FST::pitchR, 0.5, h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetMaximum(), 1, 2, 2);
    PlotLine(-3.5*FST::pitchR, -3.5*FST::pitchR, 0.5, h_mFstRecoResR_2Layer_Rstrips[i_rstrp]->GetMaximum(), 1, 2, 2);
  }
}
