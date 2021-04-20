#include <iostream>
#include <string>

#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include "../FstPlotConfiguration/draw.h"
#include "../../src/FstUtil/FstCons.h"

using namespace std;

void calMcEfficiency(int sensorId = 0)
{
  string inputname = Form("../../output/simulation/FstMcProjection_Sensor%d.root",sensorId);
  TFile *File_InPut = TFile::Open(inputname.c_str());
  TH2F *h_mIstCounts_2Layer[FST::mFstNumMatching];
  TH1F *h_mIstCountsR_2Layer[FST::mFstNumMatching];
  TH1F *h_mIstCountsPhi_2Layer[FST::mFstNumMatching];

  TH2F *h_mFstCounts_2Layer[FST::mFstNumMatching];
  TH1F *h_mFstCountsR_2Layer[FST::mFstNumMatching];
  TH1F *h_mFstCountsPhi_2Layer[FST::mFstNumMatching];
  for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
  {
    string HistName;
    // simple clusters
    HistName = Form("h_mIstCounts_2Layer_Sensor%d_SF%d",sensorId,i_match);
    h_mIstCounts_2Layer[i_match] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mIstCounts_2Layer[i_match]->Sumw2();
    int projIstBinX0 = h_mIstCounts_2Layer[i_match]->GetYaxis()->FindBin(FST::mFstPhiMin[sensorId]+4*FST::pitchPhi);
    int projIstBinX1 = h_mIstCounts_2Layer[i_match]->GetYaxis()->FindBin(FST::mFstPhiMax[sensorId]-4*FST::pitchPhi);
    HistName = Form("h_mIstCountsR_2Layer_SF%d",i_match);
    h_mIstCountsR_2Layer[i_match] = (TH1F*)h_mIstCounts_2Layer[i_match]->ProjectionX(HistName.c_str(),projIstBinX0,projIstBinX1);

    int projIstBinY0 = h_mIstCounts_2Layer[i_match]->GetXaxis()->FindBin(FST::mFstRMin[sensorId]+0.4*FST::pitchR);
    int projIstBinY1 = h_mIstCounts_2Layer[i_match]->GetXaxis()->FindBin(FST::mFstRMax[sensorId]-0.4*FST::pitchR);
    HistName = Form("h_mIstCountsPhi_2Layer_SF%d",i_match);
    h_mIstCountsPhi_2Layer[i_match] = (TH1F*)h_mIstCounts_2Layer[i_match]->ProjectionY(HistName.c_str(),projIstBinY0,projIstBinY1);

    HistName = Form("h_mFstCounts_2Layer_Sensor%d_SF%d",sensorId,i_match);
    h_mFstCounts_2Layer[i_match] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mFstCounts_2Layer[i_match]->Sumw2();
    int projFstBinX0 = h_mFstCounts_2Layer[i_match]->GetYaxis()->FindBin(FST::mFstPhiMin[sensorId]+4*FST::pitchPhi);
    int projFstBinX1 = h_mFstCounts_2Layer[i_match]->GetYaxis()->FindBin(FST::mFstPhiMax[sensorId]-4*FST::pitchPhi);
    HistName = Form("h_mFstCountsR_2Layer_SF%d",i_match);
    h_mFstCountsR_2Layer[i_match] = (TH1F*)h_mFstCounts_2Layer[i_match]->ProjectionX(HistName.c_str(),projFstBinX0,projFstBinX1);

    int projFstBinY0 = h_mFstCounts_2Layer[i_match]->GetXaxis()->FindBin(FST::mFstRMin[sensorId]+0.4*FST::pitchR);
    int projFstBinY1 = h_mFstCounts_2Layer[i_match]->GetXaxis()->FindBin(FST::mFstRMax[sensorId]-0.4*FST::pitchR);
    HistName = Form("h_mFstCountsPhi_2Layer_SF%d",i_match);
    h_mFstCountsPhi_2Layer[i_match] = (TH1F*)h_mFstCounts_2Layer[i_match]->ProjectionY(HistName.c_str(),projFstBinY0,projFstBinY1);
  }

  TH2F *h_mFstEff_2Layer[FST::mFstNumMatching];
  TH1F *h_mFstEffR_2Layer[FST::mFstNumMatching];
  TH1F *h_mFstEffPhi_2Layer[FST::mFstNumMatching];
  for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
  {
    string HistName;
    HistName = Form("h_mFstEff_2Layer_SF%d",i_match);
    h_mFstEff_2Layer[i_match] = (TH2F*)h_mIstCounts_2Layer[i_match]->Clone(HistName.c_str());
    h_mFstEff_2Layer[i_match]->SetTitle(HistName.c_str());
    h_mFstEff_2Layer[i_match]->Reset();
    h_mFstEff_2Layer[i_match]->Divide(h_mFstCounts_2Layer[i_match],h_mIstCounts_2Layer[i_match],1,1,"B");

    HistName = Form("h_mFstEffR_2Layer_SF%d",i_match);
    h_mFstEffR_2Layer[i_match] = (TH1F*)h_mIstCountsR_2Layer[i_match]->Clone(HistName.c_str());
    h_mFstEffR_2Layer[i_match]->SetTitle(HistName.c_str());
    h_mFstEffR_2Layer[i_match]->Reset();
    h_mFstEffR_2Layer[i_match]->Divide(h_mFstCountsR_2Layer[i_match],h_mIstCountsR_2Layer[i_match],1,1,"B");

    HistName = Form("h_mFstEffPhi_2Layer_SF%d",i_match);
    h_mFstEffPhi_2Layer[i_match] = (TH1F*)h_mIstCountsPhi_2Layer[i_match]->Clone(HistName.c_str());
    h_mFstEffPhi_2Layer[i_match]->SetTitle(HistName.c_str());
    h_mFstEffPhi_2Layer[i_match]->Reset();
    h_mFstEffPhi_2Layer[i_match]->Divide(h_mFstCountsPhi_2Layer[i_match],h_mIstCountsPhi_2Layer[i_match],1,1,"B");
  }


  string outputname = Form("./figures/McEfficiency_Sensor%d.pdf",sensorId);
  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,900,900);
  c_play->Divide(3,3);
  for(int i_pad = 0; i_pad < 9; ++i_pad)
  {
    c_play->cd(i_pad+1)->SetLeftMargin(0.15);
    c_play->cd(i_pad+1)->SetBottomMargin(0.15);
    c_play->cd(i_pad+1)->SetTicks(1,1);
    c_play->cd(i_pad+1)->SetGrid(0,0);
  }

  string output_start= Form("./figures/McEfficiency_Sensor%d.pdf[",sensorId);
  c_play->Print(output_start.c_str()); // open pdf file

  for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
  {
    // projection position from IST
    c_play->cd(1);
    h_mIstCounts_2Layer[i_match]->GetXaxis()->SetTitle("r_{proj} (mm)");
    h_mIstCounts_2Layer[i_match]->GetYaxis()->SetTitle("phi_{proj} (rad)");
    h_mIstCounts_2Layer[i_match]->Draw("colz");
    PlotLine(FST::mFstRMin[sensorId], FST::mFstRMax[sensorId], FST::mFstPhiMin[sensorId], FST::mFstPhiMin[sensorId], 1, 2, 2);
    PlotLine(FST::mFstRMin[sensorId], FST::mFstRMax[sensorId], FST::mFstPhiMax[sensorId], FST::mFstPhiMax[sensorId], 1, 2, 2);
    PlotLine(FST::mFstRMin[sensorId], FST::mFstRMin[sensorId], FST::mFstPhiMin[sensorId], FST::mFstPhiMax[sensorId], 1, 2, 2);
    PlotLine(FST::mFstRMax[sensorId], FST::mFstRMax[sensorId], FST::mFstPhiMin[sensorId], FST::mFstPhiMax[sensorId], 1, 2, 2);

    c_play->cd(2);
    h_mIstCountsR_2Layer[i_match]->GetXaxis()->SetTitle("r_{proj} (mm)");
    h_mIstCountsR_2Layer[i_match]->GetYaxis()->SetRangeUser(0.0,1.1*h_mIstCountsR_2Layer[i_match]->GetMaximum());
    h_mIstCountsR_2Layer[i_match]->Draw("HIST");
    PlotLine(FST::mFstRMin[sensorId], FST::mFstRMin[sensorId], 0.5, h_mFstCountsR_2Layer[i_match]->GetMaximum(), 1, 2, 2);
    PlotLine(FST::mFstRMax[sensorId], FST::mFstRMax[sensorId], 0.5, h_mFstCountsR_2Layer[i_match]->GetMaximum(), 1, 2, 2);

    c_play->cd(3);
    h_mIstCountsPhi_2Layer[i_match]->GetXaxis()->SetTitle("phi_{proj} (mm)");
    // h_mIstCountsPhi_2Layer[i_match]->GetXaxis()->SetRangeUser(-1.2*phiMaxFst,1.2*phiMaxFst);
    h_mIstCountsPhi_2Layer[i_match]->Draw("HIST");
    PlotLine(FST::mFstPhiMin[sensorId], FST::mFstPhiMin[sensorId], 0.5, h_mFstCountsPhi_2Layer[i_match]->GetMaximum(), 1, 2, 2);
    PlotLine(FST::mFstPhiMax[sensorId], FST::mFstPhiMax[sensorId], 0.5, h_mFstCountsPhi_2Layer[i_match]->GetMaximum(), 1, 2, 2);

    // measured position from FST
    c_play->cd(4);
    h_mFstCounts_2Layer[i_match]->GetXaxis()->SetTitle("r_{proj} (mm)");
    h_mFstCounts_2Layer[i_match]->GetYaxis()->SetTitle("phi_{proj} (rad)");
    h_mFstCounts_2Layer[i_match]->Draw("colz");
    PlotLine(FST::mFstRMin[sensorId], FST::mFstRMax[sensorId], FST::mFstPhiMin[sensorId], FST::mFstPhiMin[sensorId], 1, 2, 2);
    PlotLine(FST::mFstRMin[sensorId], FST::mFstRMax[sensorId], FST::mFstPhiMax[sensorId], FST::mFstPhiMax[sensorId], 1, 2, 2);
    PlotLine(FST::mFstRMin[sensorId], FST::mFstRMin[sensorId], FST::mFstPhiMin[sensorId], FST::mFstPhiMax[sensorId], 1, 2, 2);
    PlotLine(FST::mFstRMax[sensorId], FST::mFstRMax[sensorId], FST::mFstPhiMin[sensorId], FST::mFstPhiMax[sensorId], 1, 2, 2);

    c_play->cd(5);
    h_mFstCountsR_2Layer[i_match]->GetXaxis()->SetTitle("r_{proj} (mm)");
    h_mFstCountsR_2Layer[i_match]->GetYaxis()->SetRangeUser(0.0,1.1*h_mFstCountsR_2Layer[i_match]->GetMaximum());
    h_mFstCountsR_2Layer[i_match]->Draw("HIST");
    PlotLine(FST::mFstRMin[sensorId], FST::mFstRMin[sensorId], 0.5, h_mFstCountsR_2Layer[i_match]->GetMaximum(), 1, 2, 2);
    PlotLine(FST::mFstRMax[sensorId], FST::mFstRMax[sensorId], 0.5, h_mFstCountsR_2Layer[i_match]->GetMaximum(), 1, 2, 2);

    c_play->cd(6);
    h_mFstCountsPhi_2Layer[i_match]->GetXaxis()->SetTitle("phi_{proj} (mm)");
    // h_mFstCountsPhi_2Layer[i_match]->GetXaxis()->SetRangeUser(-1.2*phiMaxFst,1.2*phiMaxFst);
    h_mFstCountsPhi_2Layer[i_match]->Draw("HIST");
    PlotLine(FST::mFstPhiMin[sensorId], FST::mFstPhiMin[sensorId], 0.5, h_mFstCountsPhi_2Layer[i_match]->GetMaximum(), 1, 2, 2);
    PlotLine(FST::mFstPhiMax[sensorId], FST::mFstPhiMax[sensorId], 0.5, h_mFstCountsPhi_2Layer[i_match]->GetMaximum(), 1, 2, 2);

    c_play->cd(7);
    h_mFstEff_2Layer[i_match]->SetStats(0);
    h_mFstEff_2Layer[i_match]->GetXaxis()->SetTitle("r_{proj} (mm)");
    h_mFstEff_2Layer[i_match]->GetYaxis()->SetTitle("phi_{proj} (rad)");
    h_mFstEff_2Layer[i_match]->Draw("colz");
    PlotLine(FST::mFstRMin[sensorId], FST::mFstRMax[sensorId], FST::mFstPhiMin[sensorId], FST::mFstPhiMin[sensorId], 1, 2, 2);
    PlotLine(FST::mFstRMin[sensorId], FST::mFstRMax[sensorId], FST::mFstPhiMax[sensorId], FST::mFstPhiMax[sensorId], 1, 2, 2);
    PlotLine(FST::mFstRMin[sensorId], FST::mFstRMin[sensorId], FST::mFstPhiMin[sensorId], FST::mFstPhiMax[sensorId], 1, 2, 2);
    PlotLine(FST::mFstRMax[sensorId], FST::mFstRMax[sensorId], FST::mFstPhiMin[sensorId], FST::mFstPhiMax[sensorId], 1, 2, 2);

    c_play->cd(8);
    h_mFstEffR_2Layer[i_match]->SetStats(0);
    h_mFstEffR_2Layer[i_match]->GetXaxis()->SetTitle("r_{proj} (mm)");
    h_mFstEffR_2Layer[i_match]->GetYaxis()->SetRangeUser(0.0,1.05);
    h_mFstEffR_2Layer[i_match]->SetLineColor(1);
    h_mFstEffR_2Layer[i_match]->Draw("pE");
    PlotLine(FST::mFstRMin[sensorId], FST::mFstRMin[sensorId], 0.0, 1.05, 1, 2, 2);
    PlotLine(FST::mFstRMax[sensorId], FST::mFstRMax[sensorId], 0.0, 1.05, 1, 2, 2);
    PlotLine(FST::rMin[sensorId], FST::rMax[sensorId], 0.90, 0.90, 4, 1, 2);
    PlotLine(FST::rMin[sensorId], FST::rMax[sensorId], 0.95, 0.95, 4, 1, 2);

    c_play->cd(9);
    h_mFstEffPhi_2Layer[i_match]->SetStats(0);
    h_mFstEffPhi_2Layer[i_match]->GetXaxis()->SetTitle("phi_{proj} (rad)");
    // h_mFstEffPhi_2Layer[i_match]->GetXaxis()->SetRangeUser(-1.2*phiMaxFst,1.2*phiMaxFst);
    h_mFstEffPhi_2Layer[i_match]->GetYaxis()->SetRangeUser(0.0,1.05);
    h_mFstEffPhi_2Layer[i_match]->SetLineColor(1);
    h_mFstEffPhi_2Layer[i_match]->Draw("pE");
    PlotLine(FST::mFstPhiMin[sensorId], FST::mFstPhiMin[sensorId], 0.0, 1.05, 1, 2, 2);
    PlotLine(FST::mFstPhiMax[sensorId], FST::mFstPhiMax[sensorId], 0.0, 1.05, 1, 2, 2);
    PlotLine(FST::phiMin[sensorId], FST::phiMax[sensorId], 0.90, 0.90, 4, 1, 2);
    PlotLine(FST::phiMin[sensorId], FST::phiMax[sensorId], 0.95, 0.95, 4, 1, 2);

    c_play->Update();
    c_play->Print(outputname.c_str());
  }

  string output_stop= Form("./figures/McEfficiency_Sensor%d.pdf]",sensorId);
  c_play->Print(output_stop.c_str()); // open pdf file

  string outputfile = Form("../../output/simulation/FstMcEfficiency_Sensor%d.root",sensorId);
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_match = 0; i_match < FST::mFstNumMatching; ++i_match)
  {
    h_mFstEff_2Layer[i_match]->Write();
    h_mFstEffR_2Layer[i_match]->Write();
    h_mFstEffPhi_2Layer[i_match]->Write();
  }
  File_OutPut->Close();
}
