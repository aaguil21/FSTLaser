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
 


void plotLaserScan(string date = "03152021_4", int numEvents_Laser = 30, int data  = 0, string mod = "Mod03", string hv = "FST70", bool isSavePed = true, bool isApplyCMNCorr = true, float nFstHitsCut = 4.0, int numOfUsedTimeBins = 2, float nFstThresholdCut2 = 2.5, float nFstThresholdCut1 = 3.5)
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


//For laser scan use histograms 
  for (int j_sensor = 0; j_sensor < FST::mFstNumSensorsPerModule; j_sensor++) {
    for (int j_event = 0; j_event < numEvents_Laser; j_event++) {
      for(int j_laserRow = 0; j_laserRow < 128; j_laserRow++) {
        for(int j_laserCol = 0; j_laserCol < 8; j_laserCol++) {
	  string histScan = Form("h_signal_sen%d_ev%d_col%d_row%d", j_sensor, j_event, j_laserCol, j_laserRow);
	  signalLaser_hits[j_sensor][j_event][j_laserCol][j_laserRow] = new TH1F(histScan.c_str(), histScan.c_str(), 200, -0.5, 4800.5);
	  histScan = Form("h_noise_sen%d_ev%d_col%d_row%d", j_sensor, j_event, j_laserCol, j_laserRow);
	  noiseLaser_hits[j_sensor][j_event][j_laserCol][j_laserRow] = new TH1F(histScan.c_str(), histScan.c_str(), 200, -0.5, 40.5);
	  histScan = Form("h_ratio_sen%d_ev%d_col%d_row%d", j_sensor, j_event, j_laserCol, j_laserRow);
	  ratioLaser_hits[j_sensor][j_event][j_laserCol][j_laserRow] = new TH1F(histScan.c_str(), histScan.c_str(), 200, -0.5, 400.5);
	  histScan = Form("h_cmnNoise_sen%d_ev%d_col%d_row%d", j_sensor, j_event, j_laserCol, j_laserRow);
	  cmnNoiseLaser_hits[j_sensor][j_event][j_laserCol][j_laserRow] = new TH1F(histScan.c_str(), histScan.c_str(), 200, -100.5, 100.5);
	  histScan = Form("h_time_sen%d_ev%d_col%d_row%d", j_sensor, j_event, j_laserCol, j_laserRow);
	  timebinLaser_hits[j_sensor][j_event][j_laserCol][j_laserRow] = new TH1F(histScan.c_str(), histScan.c_str(), FST::numTBins,0,FST::numTBins);
	}
      }  
      string hist2DScan = Form("h_sigMap_sen%d_ev%d", j_sensor, j_event);
      signalMap[j_sensor][j_event] = new TH2F(hist2DScan.c_str(),hist2DScan.c_str(), 8, 0, FST::numRStrip, 128, 0, FST::numPhiSeg);
      hist2DScan = Form("h_sigTime_sen%d_ev%d", j_sensor, j_event);
      sigTime[j_sensor][j_event] = new TH2F(hist2DScan.c_str(),hist2DScan.c_str(), FST::numTBins,0,FST::numTBins, 200, -0.5,4800.5);
      hist2DScan = Form("h_noiseTime_sen%d_ev%d", j_sensor, j_event);
      noiseTime[j_sensor][j_event] = new TH2F(hist2DScan.c_str(),hist2DScan.c_str(), FST::numTBins,0,FST::numTBins, 200, -0.5,40.5);
    }
  }


// adding code so that the files with notes include the position for each event range
ifstream logFile;
string EventFile = Form("../../../logs/%s_log.txt", date.c_str());
logFile.open(EventFile.c_str());
string positionHolder;
vector<string> positionText;

while (getline(logFile,positionHolder)) {

  std::string::size_type pos = positionHolder.find("y");
  positionHolder.erase(0, pos);
  std::string tune = Form("tune = 51.50");
  std::string::size_type pos2 = positionHolder.find(tune.c_str());
  if (pos2 != std::string::npos){
    positionHolder.erase(pos2, tune.length()+3);
  }
  positionText.push_back(positionHolder);
} 


//TTree loading  for laser scan
   TTree *evScan; File_InPut->GetObject("Event_Signal", evScan);
   std::vector<Double_t> *signalScan   = 0;
   std::vector<Double_t> *noiseScan    = 0;
   std::vector<Double_t> *cmnNoiseScan = 0;
   std::vector<Int_t> *eventIdScan  = 0;
   std::vector<Int_t> *maxTbScan    = 0;
   std::vector<Int_t> *sensorIdScan = 0;
   std::vector<Int_t> *colNumScan   = 0;
   std::vector<Int_t> *rowNumScan   = 0;

   TBranch *bsignalScan   = 0;
   TBranch *bnoiseScan    = 0;
   TBranch *bcmnNoiseScan = 0;
   TBranch *beventIdScan  = 0;
   TBranch *bsensorIdScan = 0;
   TBranch *browNumScan   = 0;
   TBranch *bmaxTbScan    = 0;
   TBranch *bcolNumScan   = 0;

   evScan->SetBranchAddress("Branch1", &signalScan, &bsignalScan);
   evScan->SetBranchAddress("Branch2", &noiseScan, &bnoiseScan);   
   evScan->SetBranchAddress("Branch3", &maxTbScan, &bmaxTbScan);
   evScan->SetBranchAddress("Branch4", &cmnNoiseScan, &bcmnNoiseScan);
   evScan->SetBranchAddress("Branch5", &eventIdScan, &beventIdScan);
   evScan->SetBranchAddress("Branch6", &sensorIdScan, &bsensorIdScan);
   evScan->SetBranchAddress("Branch7", &rowNumScan, &browNumScan);
   evScan->SetBranchAddress("Branch8", &colNumScan, &bcolNumScan);

   Int_t numScanEntries = evScan->GetEntries();
   Long64_t firstEntry = evScan->LoadTree(numScanEntries-1);
   beventIdScan->GetEntry(firstEntry);
   bsignalScan->GetEntry(firstEntry);
   bnoiseScan->GetEntry(firstEntry);
   bcmnNoiseScan->GetEntry(firstEntry);
   bmaxTbScan->GetEntry(firstEntry);
   bsensorIdScan->GetEntry(firstEntry);
   bcolNumScan->GetEntry(firstEntry);
   browNumScan->GetEntry(firstEntry);
   Int_t firstEvent = 0;
   Int_t lastEvent  = 0;

   for(int i_vec = 1; i_vec < eventIdScan->size(); i_vec++) {
    if ( signalScan->at(i_vec) > 120 ) {
       firstEvent = eventIdScan->at(i_vec);
       break;
    }
   }
   Int_t size = eventIdScan->size();
   for(int i_vec = 1; i_vec < eventIdScan->size(); i_vec++) {
    if ( signalScan->at(size - i_vec) > 120 ) {
       lastEvent = eventIdScan->at(size - i_vec);
       break;
    }
   }
   Int_t evScanRange = floor((lastEvent - firstEvent)/(numEvents_Laser)  );


cout << "Event Range " << evScanRange << endl;
Double_t eventRatio = 0.0;

//Filling histograms from TTree loading event data from signal study

for(int i_vec = 1; i_vec < signalScan->size(); i_vec++) {
  for (int k_event = 0; k_event < numEvents_Laser; k_event++) {
    if( (eventIdScan->at(i_vec) >= (firstEvent + evScanRange*k_event + evScanRange*eventRatio)) && (eventIdScan->at(i_vec) < (firstEvent + evScanRange*(k_event+1) - evScanRange*eventRatio))) {
      for (int j_sensor = 0; j_sensor < FST::mFstNumSensorsPerModule; j_sensor++) {
	    if(j_sensor == int(sensorIdScan->at(i_vec))) {
          for(int j_laserRow = 0; j_laserRow < 128; j_laserRow++) {
	        if(j_laserRow == int(rowNumScan->at(i_vec))) {
              for(int j_laserCol = 0; j_laserCol < 8; j_laserCol++) {
		        if(j_laserCol == int(colNumScan->at(i_vec))) {
		          if ( signalScan->at(i_vec) > 80 ) {     // Signal threshold for the event to be plotted
        
	                signalLaser_hits[j_sensor][k_event][j_laserCol][j_laserRow]->Fill(signalScan->at(i_vec));
                    noiseLaser_hits[j_sensor][k_event][j_laserCol][j_laserRow]->Fill(noiseScan->at(i_vec));
                    cmnNoiseLaser_hits[j_sensor][k_event][j_laserCol][j_laserRow]->Fill(cmnNoiseScan->at(i_vec));
                    timebinLaser_hits[j_sensor][k_event][j_laserCol][j_laserRow]->Fill(maxTbScan->at(i_vec));
                    ratioLaser_hits[j_sensor][k_event][j_laserCol][j_laserRow]->Fill(signalScan->at(i_vec)/noiseScan->at(i_vec));
                    signalMap[j_sensor][k_event]->Fill(colNumScan->at(i_vec), rowNumScan->at(i_vec));
                    sigTime[j_sensor][k_event]->Fill( maxTbScan->at(i_vec), signalScan->at(i_vec));
                    noiseTime[j_sensor][k_event]->Fill( maxTbScan->at(i_vec), noiseScan->at(i_vec));
                      
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


int numEntries = 10;

  for (Int_t i_sensor = 0; i_sensor < FST::mFstNumSensorsPerModule; i_sensor++) {
    for (Int_t i_event = 0; i_event < numEvents_Laser; i_event++) { 

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

      int color[7] = {kRed, kBlue, kGreen, kBlack, kMagenta, kCyan, kGray}; 
      deque<int> indexRowSig = {};
      deque<int> indexColSig = {};
      int maxSig = 0;
      for(int laserRow = 0; laserRow < 128; laserRow++) {
        for(int laserCol = 0; laserCol < 8; laserCol++) {
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
        noiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetYaxis()->SetRangeUser(0,noiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetMaximum() * 1.2);
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
        ratioLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetYaxis()->SetRangeUser(0,ratioLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetMaximum() * 1.2);
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
        cmnNoiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetYaxis()->SetRangeUser(0,cmnNoiseLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetMaximum() * 1.2);
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
        timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetYaxis()->SetRangeUser(0,timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetMaximum() * 1.2);
        timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->SetLineColor(color[counterTime]);
        timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->Draw("same");
        string legName5 = Form("Ent: %d, row: %d, col: %d",int( timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)]->GetEntries()), indexRowSig.at(plotIndex), indexColSig.at(plotIndex));
        lasLeg5->AddEntry(timebinLaser_hits[i_sensor][i_event][indexColSig.at(plotIndex)][indexRowSig.at(plotIndex)], legName5.c_str(), "l");
        lasLeg5->Draw("same");
        counterTime++;
	    laserPad->Update();
       
      }

      laserPlot->Update();
      string pdfName = Form("./figures/LaserEventScanTest_%s.pdf", date.c_str());
      if (i_event == 0 && i_sensor == 0) pdfName = Form("./figures/LaserEventScanTest_%s.pdf(", date.c_str());
      if (i_event == (numEvents_Laser-1) && i_sensor == 2) pdfName = Form("./figures/LaserEventScanTest_%s.pdf)", date.c_str());
      laserPlot->SaveAs(pdfName.c_str());
    }
  }
}

