#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>

#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include "./draw.h"
#include "../../src/FstUtil/FstCons.h"

using namespace std;

/*
	This macro produces a file comparing:
	1. Efficiency results with and without CMN correction for various configs;
	2. Efficiency results for various configs with CMN correction;
	3. Efficiency results for various configs without CMN correction.
*/

void plotEfficiencyTS_HV140V()
{
  const int numOfTS = 32;
  string date[numOfTS] = {"0329","0330","0331","0402","0403","0404","0406","0407","0408","0413","0417","0423","0427","0429","0502","0504","0507","0509","0511","0512","0514","0516","0518","0519","0520","0521","0525","0526","0527","0528","0529","0531"};
  
  const int numOfConfigs = 8;
  string configs[numOfConfigs] = {"../../output/timesequence/FstTracking_HV140V_Th4.5Tb3_withPed_withCMNCorr_",
  				   "../../output/timesequence/FstTracking_HV140V_Th4.5Tb2_withPed_withCMNCorr_",
  				   "../../output/timesequence/FstTracking_HV140V_Th4.0Tb3_withPed_withCMNCorr_",
  				   "../../output/timesequence/FstTracking_HV140V_Th4.0Tb2_withPed_withCMNCorr_",
  				   "../../output/timesequence/FstTracking_HV140V_Th4.5Tb3_withPed_woCMNCorr_",
  				   "../../output/timesequence/FstTracking_HV140V_Th4.5Tb2_withPed_woCMNCorr_",
  				   "../../output/timesequence/FstTracking_HV140V_Th4.0Tb3_withPed_woCMNCorr_",
  				   "../../output/timesequence/FstTracking_HV140V_Th4.0Tb2_withPed_woCMNCorr_"
  				   };
  
  const int numOfThTb = 4;
  string ThTbAry[numOfThTb] = {"Th4.5Tb3", "Th4.5Tb2", "Th4.0Tb3", "Th4.0Tb2"};
  				  
  TFile *File_InPut[numOfConfigs][numOfTS];
  
  TH1F *h_mMatchingEfficiencyTS_Ratio_config[2][4][5][8]; // w/wo CMN | ThTb | Which algorithm | Which Scale Factor
  TH1F *h_mMatchingEfficiencyTS_raw_config[2][4][5][8];

  // initializing histograms
  for(int i_match = 0; i_match < 8; ++i_match)
  { 
    string HistName;
    // configs with and without CMN histograms
    for(int i_CMN = 0; i_CMN < 2; ++i_CMN)
    {
      string CMNstring;
      if(i_CMN == 0) CMNstring = "withCMN";
      else if(i_CMN == 1) CMNstring = "woCMN";
      for(int i_ThTb = 0; i_ThTb < 4; ++i_ThTb)
      {
        string ThTb;
        if(i_ThTb == 0) ThTb = "Th4.5Tb3";
        else if(i_ThTb == 1) ThTb = "Th4.5Tb2";
        else if(i_ThTb == 2) ThTb = "Th4.0Tb3";
        else if(i_ThTb == 3) ThTb = "Th4.0Tb2";
        
        // raw histograms with CMN
        HistName = Form("h_mMatchedHitEfficiencyRaw_SF%d_%s_%s",i_match,CMNstring.c_str(),ThTb.c_str());
        h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][0][i_match] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
        HistName = Form("h_mMatchedSimpleEfficiencyRaw_2Layer_SF%d_%s_%s",i_match,CMNstring.c_str(),ThTb.c_str());
        h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][1][i_match] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
        HistName = Form("h_mMatchedSimpleEfficiencyRaw_3Layer_SF%d_%s_%s",i_match,CMNstring.c_str(),ThTb.c_str());
        h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][2][i_match] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
        HistName = Form("h_mMatchedScanEfficiencyRaw_2Layer_SF%d_%s_%s",i_match,CMNstring.c_str(),ThTb.c_str());
        h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][3][i_match] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
        HistName = Form("h_mMatchedScanEfficiencyRaw_3Layer_SF%d_%s_%s",i_match,CMNstring.c_str(),ThTb.c_str());
        h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][4][i_match] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
         
        // Ratio Histograms
        HistName = Form("h_mMatchedHitEfficiency_SF%d_%s_%s",i_match,CMNstring.c_str(),ThTb.c_str());
        h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][0][i_match] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
        HistName = Form("h_mMatchedSimpleEfficiency_2Layer_SF%d_%s_%s",i_match,CMNstring.c_str(),ThTb.c_str());
        h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][1][i_match] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
        HistName = Form("h_mMatchedSimpleEfficiency_3Layer_SF%d_%s_%s",i_match,CMNstring.c_str(),ThTb.c_str());
        h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][2][i_match] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
        HistName = Form("h_mMatchedScanEfficiency_2Layer_SF%d_%s_%s",i_match,CMNstring.c_str(),ThTb.c_str());
        h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][3][i_match] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
        HistName = Form("h_mMatchedScanEfficiency_3Layer_SF%d_%s_%s",i_match,CMNstring.c_str(),ThTb.c_str());
        h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][4][i_match] = new TH1F(HistName.c_str(),HistName.c_str(),numOfTS,-0.5,numOfTS-0.5);
      }
    }
  }  
  
  for(int i_config = 0; i_config < numOfConfigs; ++i_config)
  {
    int numOfUsedTS = 0; 
    double *firstVals[5];
    for(int i_ts = 0; i_ts < numOfTS; ++i_ts)
    { 
      string inputfile = Form("../../output/timesequence/%s%s2020.root",configs[i_config].c_str(),date[i_ts].c_str());
      //cout << "open inputfile: " << inputfile.c_str() << endl;
      File_InPut[i_config][i_ts] = TFile::Open(inputfile.c_str());
     
      // retrieve efficiencies
      TVectorT<double> *v_mMatchedHitEfficiency = (TVectorT<double>*)File_InPut[i_config][i_ts]->Get("v_mMatchedHitEfficiency");
      TVectorT<double> *v_mMatchedSimpleEfficiency_2Layer = (TVectorT<double>*)File_InPut[i_config][i_ts]->Get("v_mMatchedSimpleEfficiency_2Layer");
      TVectorT<double> *v_mMatchedSimpleEfficiency_3Layer = (TVectorT<double>*)File_InPut[i_config][i_ts]->Get("v_mMatchedSimpleEfficiency_3Layer");
      TVectorT<double> *v_mMatchedScanEfficiency_2Layer = (TVectorT<double>*)File_InPut[i_config][i_ts]->Get("v_mMatchedScanEfficiency_2Layer");
      TVectorT<double> *v_mMatchedScanEfficiency_3Layer = (TVectorT<double>*)File_InPut[i_config][i_ts]->Get("v_mMatchedScanEfficiency_3Layer");
     
      // save efficiencies to a matrix structure
      double *effMatrix[5];
      effMatrix[0] = v_mMatchedHitEfficiency->GetMatrixArray();
      effMatrix[1] = v_mMatchedSimpleEfficiency_2Layer->GetMatrixArray();
      effMatrix[2] = v_mMatchedSimpleEfficiency_3Layer->GetMatrixArray();
      effMatrix[3] = v_mMatchedScanEfficiency_2Layer->GetMatrixArray();
      effMatrix[4] = v_mMatchedScanEfficiency_3Layer->GetMatrixArray();
      
      // copy first values
      if(i_ts == 0) std::copy(std::begin(effMatrix), std::end(effMatrix), std::begin(firstVals));
      
      // filling histograms
      for(int i_match = 0; i_match < 8; ++i_match)
      {
        for(int i_type = 0; i_type < 5; ++i_type)
        {
          int i_CMN = -1;
    	  int i_ThTb = -1;
          if(i_config >= 0 && i_config <= 3) 
	  {
            i_CMN = 0;
            i_ThTb = i_config;
          }
          else 
	  {
	    i_CMN = 1;
	    i_ThTb = i_config - 4;
          }
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_match]->SetBinContent(numOfUsedTS+1, effMatrix[i_type][i_match]);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_match]->SetBinError(numOfUsedTS+1, 0.0);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_match]->SetBinContent(numOfUsedTS+1, effMatrix[i_type][i_match]/firstVals[i_type][i_match]);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_match]->SetBinError(numOfUsedTS+1, 0.0);
        }
      }
      numOfUsedTS++;
    }
  }

  std::string outputname;
  std::string output_start;
  std::string output_stop;

  ////////////////Start Plotting///////////////////////////////////

  TCanvas *c_EffRaw_CMN = new TCanvas("c_EffRaw_CMN","c_EffRaw_CMN",10,10,1000,1600);
  c_EffRaw_CMN->Divide(1,8);
  for(int i_pad = 0; i_pad < 8; ++i_pad)
  {
    c_EffRaw_CMN->cd(i_pad+1)->SetLeftMargin(0.125);
    c_EffRaw_CMN->cd(i_pad+1)->SetRightMargin(0.0);
    c_EffRaw_CMN->cd(i_pad+1)->SetTopMargin(0.10);
    c_EffRaw_CMN->cd(i_pad+1)->SetTicks(1,1);
    c_EffRaw_CMN->cd(i_pad+1)->SetGrid(0,0);
  }

  TCanvas *c_EffRatio_CMN = new TCanvas("c_EffRatio_CMN","c_EffRatio_CMN",10,10,1000,1600);
  c_EffRatio_CMN->Divide(1,8);
  for(int i_pad = 0; i_pad < 8; ++i_pad)
  {
    c_EffRatio_CMN->cd(i_pad+1)->SetLeftMargin(0.125);
    c_EffRatio_CMN->cd(i_pad+1)->SetRightMargin(0.0);
    c_EffRatio_CMN->cd(i_pad+1)->SetTopMargin(0.10);
    c_EffRatio_CMN->cd(i_pad+1)->SetTicks(1,1);
    c_EffRatio_CMN->cd(i_pad+1)->SetGrid(0,0);
  }

  for(int i_ThTb = 0; i_ThTb < 4; ++i_ThTb)
  {
    // set an output name and initialize a canvas for CMN plots
    outputname = Form("./figures/EfficiencyCMNComparison_HV140V_%s.pdf",ThTbAry[i_ThTb].c_str());

    output_start = Form("./figures/EfficiencyCMNComparison_HV140V_%s.pdf[",ThTbAry[i_ThTb].c_str());
    c_EffRaw_CMN->Print(output_start.c_str()); // open pdf file

    for(int i_type = 0; i_type < 4;)
    {
      for(int i_SF = 0; i_SF < 8; ++i_SF)
      {
        c_EffRaw_CMN->cd(i_SF+1)->Clear();
        c_EffRaw_CMN->cd(i_SF+1);
        
        TLegend *leg_FST = new TLegend(0.1,0.2,0.2,0.4);
        leg_FST->SetBorderSize(0);
        leg_FST->SetFillColor(10);
        leg_FST->SetFillStyle(0);
        for(int i_CMN = 0; i_CMN < 2; i_CMN++)
        {
          string title;
          if(i_type == 0)
          {
            title = Form("FST Hit Efficiency SF%d Raw %s",i_SF,ThTbAry[i_ThTb].c_str());
          }
          if(i_type == 1)
          {
            title = Form("FST Simple Cluster Efficiency (2 layer) SF%d Raw %s",i_SF,ThTbAry[i_ThTb].c_str());
          }
          if(i_type == 2)
          {
            title = Form("FST Simple Cluster Efficiency (3 layer) SF%d Raw %s",i_SF,ThTbAry[i_ThTb].c_str());
          }
          if(i_type == 3)
          {
            title = Form("FST Scan Cluster Efficiency (2 layer) SF%d Raw %s",i_SF,ThTbAry[i_ThTb].c_str());
          }
          if(i_type == 4)
          {
            title = Form("FST Scan Cluster Efficiency (3 layer) SF%d Raw %s",i_SF,ThTbAry[i_ThTb].c_str());
          }
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->SetTitle(title.c_str());
          gStyle->SetTitleFontSize(0.1);
          gStyle->SetTitleY(1.0);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->SetStats(0);
          for(int i_ts = 0; i_ts < numOfTS; ++i_ts)
          {
            h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->GetXaxis()->SetBinLabel(i_ts+1,date[i_ts].c_str());
          }
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->GetXaxis()->SetLabelSize(0.13);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->GetYaxis()->SetTitle("Efficiency");
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->GetYaxis()->SetTitleSize(0.10);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->GetYaxis()->SetTitleOffset(0.4);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->GetYaxis()->SetRangeUser(0.0,1.1);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->GetYaxis()->SetLabelSize(0.08);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->GetYaxis()->SetNdivisions(505);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->SetLineColor(i_CMN+1);
          
          if(i_CMN == 0) h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->Draw();
          else h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF]->Draw("same");
          
          string LegName;
          if(i_CMN == 0) LegName = "With CMN Corr";
          else if(i_CMN == 1) LegName = "Without CMN Corr";
          leg_FST->AddEntry(h_mMatchingEfficiencyTS_raw_config[i_CMN][i_ThTb][i_type][i_SF],LegName.c_str(),"L");
        }
        leg_FST->Draw("same");
      }
        
      c_EffRaw_CMN->Update();
      c_EffRaw_CMN->Print(outputname.c_str());
      // Now for ratio histograms
      //c_EffRatio_CMN->Print(output_start.c_str()); // open pdf file
      for(int i_SF = 0; i_SF < 8; ++i_SF)
      {
        c_EffRatio_CMN->cd(i_SF+1)->Clear();
        c_EffRatio_CMN->cd(i_SF+1);
      
        TLegend *leg_FST = new TLegend(0.1,0.2,0.2,0.4);
        leg_FST->SetBorderSize(0);
        leg_FST->SetFillColor(10);
        leg_FST->SetFillStyle(0);
        for(int i_CMN = 0; i_CMN < 2; i_CMN++)
        {
          string title;
          if(i_type == 0)
          {
            title = Form("FST Hit Efficiency SF%d Ratio %s",i_SF,ThTbAry[i_ThTb].c_str());
          }
          if(i_type == 1)
          {
            title = Form("FST Simple Cluster Efficiency (2 layer) SF%d Ratio %s",i_SF,ThTbAry[i_ThTb].c_str());
          }
          if(i_type == 2)
          {
            title = Form("FST Simple Cluster Efficiency (3 layer) SF%d Ratio %s",i_SF,ThTbAry[i_ThTb].c_str());
          }
          if(i_type == 3)
          {
            title = Form("FST Scan Cluster Efficiency (2 layer) SF%d Ratio %s",i_SF,ThTbAry[i_ThTb].c_str());
          }
          if(i_type == 4)
          {
            title = Form("FST Scan Cluster Efficiency (3 layer) SF%d Ratio %s",i_SF,ThTbAry[i_ThTb].c_str());
          }
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->SetTitle(title.c_str());
          gStyle->SetTitleFontSize(0.1);
          gStyle->SetTitleY(1.0);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->SetStats(0);
          for(int i_ts = 0; i_ts < numOfTS; ++i_ts)
          {
            h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->GetXaxis()->SetBinLabel(i_ts+1,date[i_ts].c_str());
          }
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->GetXaxis()->SetLabelSize(0.13);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->GetYaxis()->SetTitle("#frac{Efficiency}{Efficiency_{@0329}}");
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->GetYaxis()->SetTitleSize(0.10);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->GetYaxis()->SetTitleOffset(0.4);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->GetYaxis()->SetRangeUser(0.4,1.2);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->GetYaxis()->SetLabelSize(0.08);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->GetYaxis()->SetNdivisions(505);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->SetLineColor(i_CMN+1);
          
          if(i_CMN == 0) h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->Draw();
          else h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF]->Draw("same");
          
          string LegName;
          if(i_CMN == 0) LegName = "With CMN Corr";
          else if(i_CMN == 1) LegName = "Without CMN Corr";
          leg_FST->AddEntry(h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_ThTb][i_type][i_SF],LegName.c_str(),"L");
        }
        leg_FST->Draw("same");
      }
      c_EffRatio_CMN->Update();
      c_EffRatio_CMN->Print(outputname.c_str());
      if(i_type == 0) i_type++;
      else if(i_type == 1) i_type+=2;
      else i_type++;
    }

    output_stop = Form("./figures/EfficiencyCMNComparison_HV140V_%s.pdf]",ThTbAry[i_ThTb].c_str());  
    c_EffRaw_CMN->Print(output_stop.c_str()); // open pdf file
  }

  TCanvas *c_EffRaw_Configs = new TCanvas("c_EffRaw_Configs","c_EffRaw_Configs",10,10,1000,1600);
  c_EffRaw_Configs->Divide(1,8);
  for(int i_pad = 0; i_pad < 8; ++i_pad)
  {
    c_EffRaw_Configs->cd(i_pad+1)->SetLeftMargin(0.125);
    c_EffRaw_Configs->cd(i_pad+1)->SetRightMargin(0.0);
    c_EffRaw_Configs->cd(i_pad+1)->SetTopMargin(0.10);
    c_EffRaw_Configs->cd(i_pad+1)->SetTicks(1,1);
    c_EffRaw_Configs->cd(i_pad+1)->SetGrid(0,0);
  }

  TCanvas *c_EffRatio_Configs = new TCanvas("c_EffRatio_Configs","c_EffRatio_Configs",10,10,1000,1600);
  c_EffRatio_Configs->Divide(1,8);
  for(int i_pad = 0; i_pad < 8; ++i_pad)
  {
    c_EffRatio_Configs->cd(i_pad+1)->SetLeftMargin(0.125);
    c_EffRatio_Configs->cd(i_pad+1)->SetRightMargin(0.0);
    c_EffRatio_Configs->cd(i_pad+1)->SetTopMargin(0.10);
    c_EffRatio_Configs->cd(i_pad+1)->SetTicks(1,1);
    c_EffRatio_Configs->cd(i_pad+1)->SetGrid(0,0);
  }
  // set an output name and initialize a canvas for config plots with CMN
  for(int i_CMN = 0; i_CMN < 2; ++i_CMN)
  {
    string CMNstring;
    if(i_CMN == 0) CMNstring = "withCMNCorr";
    else if(i_CMN == 1) CMNstring = "woCMNCorr";
    outputname = Form("./figures/EfficiencyConfigComparison_HV140V_%s.pdf",CMNstring.c_str());
    
    output_start = Form("./figures/EfficiencyConfigComparison_HV140V_%s.pdf[",CMNstring.c_str());

    c_EffRaw_Configs->Print(output_start.c_str()); // open pdf file

    for(int i_type = 0; i_type < 4;)
    {
      for(int i_SF = 0; i_SF < 8; ++i_SF)
      {
        c_EffRaw_Configs->cd(i_SF+1)->Clear();
        c_EffRaw_Configs->cd(i_SF+1);
        
        TLegend *leg_FST = new TLegend(0.1,0.2,0.2,0.4);
        leg_FST->SetBorderSize(0);
        leg_FST->SetFillColor(10);
        leg_FST->SetFillStyle(0);
        for(int i_config = 0; i_config < 4; i_config++)
        {
          string title;
          if(i_type == 0)
          {
            title = Form("FST Hit Efficiency SF%d Raw %s",i_SF,CMNstring.c_str());
          }
          if(i_type == 1)
          {
            title = Form("FST Simple Cluster Efficiency (2 layer) SF%d Raw %s",i_SF,CMNstring.c_str());
          }
          if(i_type == 2)
          {
            title = Form("FST Simple Cluster Efficiency (3 layer) SF%d Raw %s",i_SF,CMNstring.c_str());
          }
          if(i_type == 3)
          {
            title = Form("FST Scan Cluster Efficiency (2 layer) SF%d Raw %s",i_SF,CMNstring.c_str());
          }
          if(i_type == 4)
          {
            title = Form("FST Scan Cluster Efficiency (3 layer) SF%d Raw %s",i_SF,CMNstring.c_str());
          }
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->SetTitle(title.c_str());
          gStyle->SetTitleFontSize(0.1);
          gStyle->SetTitleY(1.0);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->SetStats(0);
          for(int i_ts = 0; i_ts < numOfTS; ++i_ts)
          {
            h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->GetXaxis()->SetBinLabel(i_ts+1,date[i_ts].c_str());
          }
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->GetXaxis()->SetLabelSize(0.13);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->GetYaxis()->SetTitle("Efficiency");
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->GetYaxis()->SetTitleSize(0.10);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->GetYaxis()->SetTitleOffset(0.4);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->GetYaxis()->SetRangeUser(0.2,1.05);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->GetYaxis()->SetLabelSize(0.08);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->GetYaxis()->SetNdivisions(505);
          h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->SetLineColor(i_config+1);
          
          if(i_config == 0) h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->Draw();
          else h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF]->Draw("same");
          
          string LegName;
          LegName = ThTbAry[i_config];
          leg_FST->AddEntry(h_mMatchingEfficiencyTS_raw_config[i_CMN][i_config][i_type][i_SF],LegName.c_str(),"L");
        }
        leg_FST->Draw("same");
      }
        
      c_EffRaw_Configs->Update();
      c_EffRaw_Configs->Print(outputname.c_str());
     
      for(int i_SF = 0; i_SF < 8; ++i_SF)
      {
        c_EffRatio_Configs->cd(i_SF+1)->Clear();
        c_EffRatio_Configs->cd(i_SF+1);
      
        TLegend *leg_FST = new TLegend(0.1,0.2,0.2,0.4);
        leg_FST->SetBorderSize(0);
        leg_FST->SetFillColor(10);
        leg_FST->SetFillStyle(0);
        for(int i_config = 0; i_config < 4; i_config++)
        {
          string title;
          if(i_type == 0)
          {
            title = Form("FST Hit Efficiency SF%d Ratio %s",i_SF,CMNstring.c_str());
          }
          if(i_type == 1)
          {
            title = Form("FST Simple Cluster Efficiency (2 layer) SF%d Ratio %s",i_SF,CMNstring.c_str());
          }
          if(i_type == 2)
          {
            title = Form("FST Simple Cluster Efficiency (3 layer) SF%d Ratio %s",i_SF,CMNstring.c_str());
          }
          if(i_type == 3)
          {
            title = Form("FST Scan Cluster Efficiency (2 layer) SF%d Ratio %s",i_SF,CMNstring.c_str());
          }
          if(i_type == 4)
          {
            title = Form("FST Scan Cluster Efficiency (3 layer) SF%d Ratio %s",i_SF,CMNstring.c_str());
          }
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->SetTitle(title.c_str());
          gStyle->SetTitleFontSize(0.1);
          gStyle->SetTitleY(1.0);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->SetStats(0);
          for(int i_ts = 0; i_ts < numOfTS; ++i_ts)
          {
            h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->GetXaxis()->SetBinLabel(i_ts+1,date[i_ts].c_str());
          }
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->GetXaxis()->SetLabelSize(0.13);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->GetYaxis()->SetTitle("#frac{Efficiency}{Efficiency_{@0329}}");
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->GetYaxis()->SetTitleSize(0.10);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->GetYaxis()->SetTitleOffset(0.4);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->GetYaxis()->SetRangeUser(0.5,1.25);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->GetYaxis()->SetLabelSize(0.08);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->GetYaxis()->SetNdivisions(505);
          h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->SetLineColor(i_config+1);
          
          if(i_config == 0) h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->Draw();
          else h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF]->Draw("same");
          
          string LegName;
          LegName = ThTbAry[i_config];
          leg_FST->AddEntry(h_mMatchingEfficiencyTS_Ratio_config[i_CMN][i_config][i_type][i_SF],LegName.c_str(),"L");
        }
        leg_FST->Draw("same");
      }
      c_EffRatio_Configs->Update();
      c_EffRatio_Configs->Print(outputname.c_str());
      
      if(i_type == 0) i_type++;
      else if(i_type == 1) i_type+=2;
      else i_type++;
    }
    output_stop = Form("./figures/EfficiencyConfigComparison_HV140V_%s.pdf]",CMNstring.c_str());  
    c_EffRatio_Configs->Print(output_stop.c_str()); // open pdf file
  }
}
