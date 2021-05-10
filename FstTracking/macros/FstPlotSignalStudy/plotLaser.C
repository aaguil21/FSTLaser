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
#include <deque>
#include <fstream>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <vector>
 


void plotLaserScan(string date = "01272021_1", int numEvents_Laser = 21, int data  = 0, string mod = "Mod03", string hv = "FST70", bool isSavePed = true, bool isApplyCMNCorr = true, float nFstHitsCut = 4.0, int numOfUsedTimeBins = 2, float nFstThresholdCut2 = 2.5, float nFstThresholdCut1 = 3.5)
{
  gROOT->Reset();

  // gStyle->SetPalette(kRainBow);
  gStyle->SetPalette(kBlackBody);
  gStyle->SetPaintTextFormat("1.0f");

  std::string pedMode = "withPed";
  if(!isSavePed) pedMode = "woPed";
  std::string cmnMode = "withCMNCorr";
  if(!isApplyCMNCorr) cmnMode = "woCMNCorr";

  string inputfile = Form("../../output/signal/%s/FstQAStudy_%s_%s_Th%1.1fTb%dPed%1.1fPed%1.1f_%s_%s_%s_data_%d.root",mod.c_str(),mod.c_str(),hv.c_str(),nFstHitsCut,numOfUsedTimeBins,nFstThresholdCut2,nFstThresholdCut1,pedMode.c_str(),cmnMode.c_str(), date.c_str(), data);
  TFile *File_InPut = TFile::Open(inputfile.c_str());
    

//laser signal, noise, and ratio distributions 
  TH3F *signalLaser[FST::mFstNumSensorsPerModule][numEvents_Laser];
  TH3F *noiseLaser[FST::mFstNumSensorsPerModule][numEvents_Laser];
  TH3F *ratioLaser[FST::mFstNumSensorsPerModule][numEvents_Laser];
  TH3F *cmnNoiseLaser[FST::mFstNumSensorsPerModule][numEvents_Laser];
  TH3F *timebinLaser[FST::mFstNumSensorsPerModule][numEvents_Laser];

  TH1F *signalLaser_hits[FST::mFstNumSensorsPerModule][numEvents_Laser][FST::numRStrip][FST::numPhiSeg];
  TH1F *noiseLaser_hits[FST::mFstNumSensorsPerModule][numEvents_Laser][FST::numRStrip][FST::numPhiSeg];
  TH1F *ratioLaser_hits[FST::mFstNumSensorsPerModule][numEvents_Laser][FST::numRStrip][FST::numPhiSeg];
  TH1F *cmnNoiseLaser_hits[FST::mFstNumSensorsPerModule][numEvents_Laser][FST::numRStrip][FST::numPhiSeg];
  TH1F *timebinLaser_hits[FST::mFstNumSensorsPerModule][numEvents_Laser][FST::numRStrip][FST::numPhiSeg];

  TH2F *signalMap[FST::mFstNumSensorsPerModule][numEvents_Laser];
  TH2F *sigTime[FST::mFstNumSensorsPerModule][numEvents_Laser];
  TH2F *noiseTime[FST::mFstNumSensorsPerModule][numEvents_Laser];


int numEntries = 10;

  for (Int_t i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; i_sensor++) {
    for (Int_t i_event = 0; i_event < numEvents_Laser; i_event++) { 
      string HistName = Form("h_mFstLaser_HitsSignal_Sensor%d_Event%d", i_sensor, i_event);
      signalLaser[i_sensor][i_event] = (TH3F*)File_InPut->Get(HistName.c_str());
      HistName = Form("h_mFstLaser_NoiseSignal_Sensor%d_Event%d", i_sensor, i_event);
      noiseLaser[i_sensor][i_event] = (TH3F*)File_InPut->Get(HistName.c_str());
      HistName = Form("h_mFstLaser_SignalRatio_Sensor%d_Event%d", i_sensor, i_event);
      ratioLaser[i_sensor][i_event] = (TH3F*)File_InPut->Get(HistName.c_str());      
      HistName = Form("h_mFstLaser_CmnNoiseSignal_Sensor%d_Event%d", i_sensor, i_event);
      cmnNoiseLaser[i_sensor][i_event] = (TH3F*)File_InPut->Get(HistName.c_str());
      HistName = Form("h_mFstLaser_Timebin_Sensor%d_Event%d", i_sensor, i_event);
      timebinLaser[i_sensor][i_event] = (TH3F*)File_InPut->Get(HistName.c_str());
 
      HistName = Form("h_mFstLaser_SigTime_Sensor%d_Event%d", i_sensor, i_event);
      sigTime[i_sensor][i_event] = (TH2F*)File_InPut->Get(HistName.c_str());
      HistName = Form("h_mFstLaser_NoiseTime_Sensor%d_Event%d", i_sensor, i_event);
      noiseTime[i_sensor][i_event] = (TH2F*)File_InPut->Get(HistName.c_str());


      string canvasName = Form("laserPlot_event%d_sensor%d", i_event, i_sensor);
      TCanvas *laserPlot = new TCanvas(canvasName.c_str(),canvasName.c_str(), 800, 800);
      string label = Form("Sensor %d, %s", i_sensor, positionText.at(i_event).c_str());
      TPaveLabel* title = new TPaveLabel(0.1,0.96,0.9,0.99,label.c_str());
      title->Draw();
      string padName = Form("laserPad_event%d_sensor%d", i_event, i_sensor);
      TPad* laserPad = new TPad(padName.c_str(), padName.c_str(),0.01,0.05,0.95,0.95);
      laserPad->Draw();
      laserPad->cd();
      laserPad->Divide(3,3);      

      int counterSig = 0;
      int counterRatio = 0;
      int counterNoise = 0;
      int counterCmn = 0;
      int counterTime = 0;


      laserPad->cd(1);
	TLegend *lasLeg1 = new TLegend(0.1,0.7,0.48,0.9);
      laserPad->cd(2);
	TLegend *lasLeg2 = new TLegend(0.1,0.7,0.48,0.9);
      laserPad->cd(3);
	TLegend *lasLeg3 = new TLegend(0.1,0.7,0.48,0.9);
      laserPad->cd(4);
	TLegend *lasLeg4 = new TLegend(0.1,0.7,0.48,0.9);
      laserPad->cd(5);
	TLegend *lasLeg5 = new TLegend(0.1,0.7,0.48,0.9);

      laserPad->cd(6);
      signalMap[i_sensor][i_event] = (TH2F*)signalLaser[i_sensor][i_event]->Project3D("yx");
      signalMap[i_sensor][i_event]->GetXaxis()->SetTitle("R-Strip");
      signalMap[i_sensor][i_event]->GetYaxis()->SetTitle("Phi-Segment");
      signalMap[i_sensor][i_event]->Draw("colz");
      signalMap[i_sensor][i_event]->GetYaxis()->SetRange(64,100);
      signalMap[i_sensor][i_event]->GetXaxis()->SetRange(4,8);
      laserPad->Update();

      laserPad->cd(7);

      sigTime[i_sensor][i_event]->GetXaxis()->SetTitle("Max TimeBin");
      sigTime[i_sensor][i_event]->GetYaxis()->SetTitle("Signal");
      sigTime[i_sensor][i_event]->GetYaxis()->SetRangeUser(0,2000);
      sigTime[i_sensor][i_event]->Draw("colz");
      laserPad->Update();
      laserPad->cd(8);
      noiseTime[i_sensor][i_event]->GetXaxis()->SetTitle("Max TimeBin");
      noiseTime[i_sensor][i_event]->GetYaxis()->SetTitle("Noise");
      noiseTime[i_sensor][i_event]->GetYaxis()->SetRangeUser(0,40);
      noiseTime[i_sensor][i_event]->Draw("colz");

//      gStyle->SetOptStat(111100);
  //    gStyle->SetStatX(44);
    //  gStyle->SetStatY(8.3);
      int color[7] = {kRed, kBlue, kGreen, kBlack, kMagenta, kCyan, kGray}; 


      deque<int> indexRowSig = {};
      deque<int> indexColSig = {};
      int maxSig = 0;
      for(int laserRow = 0; laserRow < 128; laserRow++) {
        for(int laserCol = 0; laserCol < 8; laserCol++) {
 
	    string histoName = Form("signal_event%d_sensor%d_row%d_col%d", i_event, i_sensor, laserRow, laserCol);
	    signalLaser_hits[i_sensor][i_event][laserCol][laserRow] = (TH1F*)signalLaser[i_sensor][i_event]->ProjectionZ(histoName.c_str(), laserCol+1, laserCol+1, laserRow+1, laserRow+1);
	 
        histoName = Form("noise_event%d_sensor%d_row%d_col%d", i_event, i_sensor, laserRow-1, laserCol-1);
        noiseLaser_hits[i_sensor][i_event][laserCol][laserRow] = (TH1F*)noiseLaser[i_sensor][i_event]->ProjectionZ(histoName.c_str(), laserCol+1, laserCol+1, laserRow+1, laserRow+1);

        histoName = Form("ratio_event%d_sensor%d_row%d_col%d", i_event, i_sensor, laserRow-1, laserCol-1);
        ratioLaser_hits[i_sensor][i_event][laserCol][laserRow] = (TH1F*)ratioLaser[i_sensor][i_event]->ProjectionZ(histoName.c_str(), laserCol+1, laserCol+1, laserRow+1, laserRow+1);
          
        histoName = Form("cmnNoise_event%d_sensor%d_row%d_col%d", i_event, i_sensor, laserRow-1, laserCol-1);
        cmnNoiseLaser_hits[i_sensor][i_event][laserCol][laserRow] = (TH1F*)cmnNoiseLaser[i_sensor][i_event]->ProjectionZ(histoName.c_str(), laserCol+1, laserCol+1, laserRow+1, laserRow+1);

        histoName = Form("timebin_event%d_sensor%d_row%d_col%d", i_event, i_sensor, laserRow-1, laserCol-1);
        timebinLaser_hits[i_sensor][i_event][laserCol][laserRow] = (TH1F*)timebinLaser[i_sensor][i_event]->ProjectionZ(histoName.c_str(), laserCol+1, laserCol+1, laserRow+1, laserRow+1);

	      if( signalLaser_hits[i_sensor][i_event][laserCol][laserRow]->GetEntries() >= numEntries) {
            if (signalLaser_hits[i_sensor][i_event][laserCol][laserRow]->GetMaximum() > maxSig ) {
            maxSig = signalLaser_hits[i_sensor][i_event][laserCol][laserRow]->GetMaximum();
            indexRowSig.push_front(laserRow);
	        indexColSig.push_front( laserCol);
            }
	        else {
	        indexRowSig.push_back(laserRow);
            indexColSig.push_back(laserCol);

	        }
	      }
        } 
      }
      
      
      for( int plotIndex = 0; plotIndex < indexRowSig.size(); plotIndex++) {
	    laserPad->cd(1);
	    signalLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetStats(0);
	    signalLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetTitle("Signal");
	    signalLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetXaxis()->SetTitle("Signal");
	    signalLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetYaxis()->SetRangeUser(0,signalLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetMaximum() * 1.2);
	    signalLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetXaxis()->SetRangeUser(0,2000);
	    signalLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetMarkerStyle(2);
	    signalLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetLineColor(color[counterSig]);
	    signalLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->Draw("same");
	    string legName = Form("Ent: %d, row: %d, col: %d",int( signalLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetEntries()), indexRowSig.at(plotIndex), indexColSig.at(plotIndex));
	    lasLeg1->AddEntry(signalLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)], legName.c_str(), "l");
	    lasLeg1->Draw("same");
	    counterSig++;	
	    laserPad->Update();	  


            laserPad->cd(2);
            noiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetStats(0);
            noiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetTitle("Noise");
	    noiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetXaxis()->SetRangeUser(10,30);
	    noiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetXaxis()->SetTitle("Noise");
            noiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetMarkerStyle(2);
            noiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetLineColor(color[counterNoise]);
            noiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->Draw("same");
            string legName2 = Form("Ent: %d, row: %d, col: %d",int( noiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetEntries()), indexRowSig.at(plotIndex), indexColSig.at(plotIndex));
            lasLeg2->AddEntry(noiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)], legName2.c_str(), "l");
            lasLeg2->Draw("same");
            counterNoise++;
	    laserPad->Update();


            laserPad->cd(3);
            ratioLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetStats(0);
            ratioLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetTitle("Signal/Noise");
	    ratioLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetXaxis()->SetTitle("S/N");
            ratioLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetMarkerStyle(2);
            ratioLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetLineColor(color[counterRatio]);
            ratioLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->Draw("same");
            string legName3 = Form("Ent: %d, row: %d, col: %d",int( ratioLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetEntries()), indexRowSig.at(plotIndex), indexColSig.at(plotIndex));
            lasLeg3->AddEntry(ratioLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)], legName3.c_str(), "l");
            lasLeg3->Draw("same");
	    counterRatio++;
	    laserPad->Update();

            laserPad->cd(4);
            cmnNoiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetStats(0);
            cmnNoiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetTitle("Common Mode Noise");
	    cmnNoiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetXaxis()->SetTitle("CMN Noise");
            cmnNoiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetMarkerStyle(2);
            cmnNoiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetLineColor(color[counterCmn]);
            cmnNoiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->Draw("same");
            string legName4 = Form("Ent: %d, row: %d, col: %d",int( cmnNoiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetEntries()), indexRowSig.at(plotIndex), indexColSig.at(plotIndex));
            lasLeg4->AddEntry(cmnNoiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)], legName4.c_str(), "l");
            lasLeg4->Draw("same");
            counterCmn++;
	    laserPad->Update();


            laserPad->cd(5);
            timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetStats(0);
            timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetTitle("Max Timebin");
            timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetMarkerStyle(2);
	    timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetXaxis()->SetTitle("Timebin");
            timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetLineColor(color[counterTime]);
            timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->Draw("same");
            string legName5 = Form("Ent: %d, row: %d, col: %d",int( timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetEntries()), indexRowSig.at(plotIndex), indexColSig.at(plotIndex));
            lasLeg5->AddEntry(timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)], legName5.c_str(), "l");
            lasLeg5->Draw("same");
            counterTime++;
	    laserPad->Update();
       
      }

      laserPlot->Update();
     string pdfName = Form("./figures/LaserEvent_%s.pdf", date.c_str());
      if (i_event == 0 && i_sensor == 0) pdfName = Form("./figures/LaserEvent_%s.pdf(", date.c_str());
      if (i_event == (numEvents_Laser-1) && i_sensor == 2) pdfName = Form("./figures/LaserEvent_%s.pdf)", date.c_str());
      laserPlot->SaveAs(pdfName.c_str());
    }
  }

}

