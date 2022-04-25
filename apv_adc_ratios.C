#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
using namespace std::chrono;
#include "TGraph.h"
#include "TMinuit.h"

#include "./include/include_files.h"
#include "./include/APV_strips.h"
#include "./include/search_file.C"

Double_t poly2_fit(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t gaus_fit(Double_t *x, Double_t *par){
	return par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
}

//Fit functions for returning parameters:
Double_t fitFunction(Double_t *x, Double_t *par){
	return poly2_fit(x, par) + gaus_fit(x, &par[3]);
}

//Define some global variables
int apv_chan_adc[128][2];
int apv_strip_adc[128][2];

int nAPVs = 30;
int occupancy_bin_num = 10000;
int first_apv = 13;
int last_apv = 14;
int ratio_bins = 300;

//Set defaults for ADC on channel n and n+1
double chan_n_adc = 0.0;
double chan_n_plus_1_ADC = 0.0;
double ADC_ratio = 0.0;

int nstrips = 3840;
Int_t max_strips = 8000;
Int_t nAPV_strips = 128;
Int_t gaus_min_bin;
Int_t gaus_max_bin;


void apv_adc_ratios(int runnum = 13444, int const ADC_cut = 500){
	auto start = high_resolution_clock::now();
	TChain *TC = new TChain("T");

	const char * DATA_DIR =  "/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk";
	// const char * DATA_DIR = "/lustre19/expphy/volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles";
	const char * protorootfile = Form("/%i/e1209019_replayed_%i.root", runnum, runnum);

	TString rootfile = Form("%s%s", DATA_DIR, protorootfile);
	cout << "Input root file is: " << rootfile << endl;
	TC->Add( rootfile );

	//Loop for TChain
	// for(int seg=0; seg < 10; seg++){
	// 	protorootfile = Form("/%i/e1209019_replayed_%i.root", runnum, runnum);
	// 	TString rootfile = Form("%s%s", DATA_DIR, protorootfile);
	// 	TC->Add( rootfile );
		
	// }

	//CANVASES
	TCanvas *c_APV_ratio_ADCmax_chan_U[nAPVs];

	//Turn off all Branches to save memory/time
	TC->SetBranchStatus("*", false);

	//Turn on desired branches
	
	//Hit branches:
	TC->SetBranchStatus("bb.gem.hit.nstripu", true);
	TC->SetBranchStatus("bb.gem.hit.nstripv", true);
	TC->SetBranchStatus("bb.gem.hit.ADCmaxsampU", true);
	TC->SetBranchStatus("bb.gem.hit.ADCmaxsampV", true);
	TC->SetBranchStatus("singletrack", true);

	//Strip branches
	// ****Ndata branches:
	TC->SetBranchStatus("Ndata.bb.gem.m0.strip.ADCsamples", true);
	TC->SetBranchStatus("Ndata.bb.gem.m0.strip.IsU", true);
	TC->SetBranchStatus("Ndata.bb.gem.m0.strip.IsV", true);
	TC->SetBranchStatus("bb.gem.m0.strip.nstrips_keepU", true);
	//****Variable branches:
	TC->SetBranchStatus("bb.gem.m0.strip.IsU", true);
	TC->SetBranchStatus("bb.gem.m0.strip.IsV", true);
	TC->SetBranchStatus("bb.gem.m0.strip.istrip", true);
	TC->SetBranchStatus("bb.gem.m0.strip.ADCsamples", true);
	TC->SetBranchStatus("bb.gem.m0.strip.ADCmax", true);
	TC->SetBranchStatus("bb.gem.m0.strip.ADCsum", true);
	TC->SetBranchStatus("bb.gem.m0.strip.nstrips_keep", true);
	TC->SetBranchStatus("bb.gem.m0.strip.isampmax", true);


	//Define variables that will hold the branch values/variables
	//Ndata variables:
	double singletrack;
	Int_t Ndata_strip_ADC;
	Int_t Ndata_IsU;
	Int_t Ndata_IsV;
	double nstrips_keepU;
	double nstripu[max_strips];
	double nstripv[max_strips];
	double nstrips_keep;
	
	//Strip, ADC, etc, variable values:
	Int_t isampmax[max_strips];
	double IsU[max_strips];
	double IsV[max_strips];
	double istrip[max_strips];
	double ADC_samples[max_strips];
	double ADC_max[max_strips];
	double ADC_sum[max_strips];


	//Assign branches to these variables:
	//Ndata
	TC->SetBranchAddress("Ndata.bb.gem.m0.strip.ADCsamples", &Ndata_strip_ADC);
	TC->SetBranchAddress("Ndata.bb.gem.m0.strip.IsU", &Ndata_IsU);
	TC->SetBranchAddress("Ndata.bb.gem.m0.strip.IsV", &Ndata_IsV);

	//Hits
	TC->SetBranchAddress("bb.gem.hit.nstripu", &nstripu);
	TC->SetBranchAddress("bb.gem.hit.nstripv", &nstripv);

	//Strips
	TC->SetBranchAddress("bb.gem.m0.strip.ADCsamples", &ADC_samples);
	TC->SetBranchAddress("bb.gem.m0.strip.ADCmax", &ADC_max);
	TC->SetBranchAddress("bb.gem.m0.strip.ADCsum", &ADC_sum);
	TC->SetBranchAddress("bb.gem.m0.strip.IsU", &(IsU));
	TC->SetBranchAddress("bb.gem.m0.strip.IsV", &(IsV));
	TC->SetBranchAddress("bb.gem.m0.strip.istrip", &istrip);
	TC->SetBranchAddress("bb.gem.m0.strip.isampmax", &isampmax);
	TC->SetBranchAddress("bb.gem.m0.strip.nstrips_keep", &nstrips_keep);
	TC->SetBranchAddress("bb.gem.m0.strip.nstrips_keepU", &nstrips_keepU);
	TC->SetBranchAddress("singletrack", &singletrack);


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------     MAIN LOOP -------------------------------
//                            Where the fun happens!!etc********************
//This fist loop is to step through APVs. We will look at individual APVs and this is
//where that main loop starts.

	//Some general variables
	double max_U_ADC = TC->GetMaximum("bb.gem.hit.ADCmaxsampU");
	double max_V_ADC = TC->GetMaximum("bb.gem.hit.ADCmaxsampV");

	//Define histograms:
	//1-D histograms
	TH1D *h_APV_ratio_ADCmax_chan_U[nAPVs];

	//2-D histograms
	TH2D *h2_ADCmax_U = new TH2D("h2_ADCmax_samples_U", "" , nstrips, 0, nstrips, nstrips, 0, nstrips);

	cout << endl << "Run number: " << runnum << endl;
	cout << endl << "Analyzing data for APVs " << first_apv << " up to " << last_apv << "." << endl << endl;

//************************************************
//-----------------APV LOOP-----------------------
//************************************************

	for(int apv_cnt = first_apv; apv_cnt < last_apv; apv_cnt++){
		cout <<  "Processing data for apv: " << apv_cnt << endl << endl;
		
		//Max and min strips for the specific APV
		int min_APV_strip = APV_strip_nums(apv_cnt, "min");
		int max_APV_strip = APV_strip_nums(apv_cnt, "max");

	//Create ratio histogram entry for each APV

		h_APV_ratio_ADCmax_chan_U[apv_cnt] = new TH1D(Form("h_APV%i_ratio_ADCmax_chan_U", apv_cnt), "", ratio_bins, 0, nAPVs);
	
//EVENT //Loop all selected events:		
		//We can either loop over all events or select some events:
		//Loop over all events:

		int apv_event_num = 0;
		while(TC->GetEntry(apv_event_num++)){
			if(apv_event_num%10000 == 0){cout << "Analyzing event " << apv_event_num << " of " << TC->GetEntries() << " total events." << endl;}
			TC->GetEntry(apv_event_num);

		
//EVENT //Loop over selected events:

		// for(int evt = 500; evt <= 500; evt++){
		
		// 	// if(evt == 706 || evt == 1186) {continue;}
		// 	cout << "evt: " << evt << endl;
		// 	TC->GetEntry(evt);
			

			
			//------------U-strips-------------
			for(int i = 0; i < Ndata_IsU; i++){
				
				if(IsU[i]){
					h2_ADCmax_U->Fill(istrip[i], ADC_max[i]);

					//Only select strips that are in the selected range for the specific APV
					//Need to cross reference the global strip numbers and relate them to those of the APV_strip_nums

					if(istrip[i] >= min_APV_strip && istrip[i] <= max_APV_strip){
						//Fill an array with the maximum ADC on each channel.
						//Storing it in a 2D array to grab the maximum entry in each bin

						//-----CHANNELS-------
						apv_chan_adc[UV_APV_strip_to_channel(int(istrip[i])%128)][0] = UV_APV_strip_to_channel(int(istrip[i])%128);
						apv_chan_adc[UV_APV_strip_to_channel(int(istrip[i])%128)][1] = ADC_max[i];

						//-----STRIPS---------
						apv_strip_adc[int(istrip[i])%128][0] = istrip[i];
						apv_strip_adc[int(istrip[i])%128][1] = ADC_max[i];

					}
				
				}

			}
			//End of looping through U-Strips

			//----------RATIOS------------
			for(int j = 0; j < 127; j++){
				chan_n_adc = apv_chan_adc[j][1];
				chan_n_plus_1_ADC = apv_chan_adc[j+1][1];
				ADC_ratio = 0.0;
				
				//Always take larger ADC divided by smaller ADC for neighboring channels
				if(chan_n_adc > chan_n_plus_1_ADC && chan_n_plus_1_ADC != 0 && chan_n_adc > ADC_cut){
					double ADC_ratio = chan_n_adc/chan_n_plus_1_ADC;
					h_APV_ratio_ADCmax_chan_U[apv_cnt]->Fill(ADC_ratio);
				}
				else if(chan_n_plus_1_ADC > chan_n_adc && chan_n_adc != 0 && chan_n_plus_1_ADC > ADC_cut){
					double ADC_ratio = chan_n_plus_1_ADC/chan_n_adc;
					h_APV_ratio_ADCmax_chan_U[apv_cnt]->Fill(ADC_ratio);
				}

	
			}
		
		
		}
		//END of EVENTS

		//----------FITS------------//
		// Double_t par[6];
		Int_t lastXbin = (0.1)*h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetNbinsX();
		cout << "last x bin: " << lastXbin << endl;

		int gaus_last_point = 13;
		h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetXaxis()->SetRangeUser(1, 9);
		Double_t gaus_min_bin = (0.9)*((0.1)*(h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetMinimumBin()));
		Double_t landau_max_bin = (1.25)*((0.1)*(h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetMinimumBin()));
		h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetXaxis()->SetRangeUser(gaus_min_bin, gaus_last_point);

		Double_t gaus_max_bin = (0.1)*(h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetMaximumBin());
		h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetXaxis()->UnZoom();

		cout << "***********************************" << endl << endl;
		cout << "min bin: " << gaus_min_bin << "      maxb: " << gaus_max_bin << endl << endl;
		cout << "***********************************" << endl << endl;
		

		c_APV_ratio_ADCmax_chan_U[apv_cnt] = new TCanvas(Form("APV%i Ratio of Channels - Ustrips", apv_cnt), Form("c_ratio_chan_apv_%i", apv_cnt), 700, 500);
		h_APV_ratio_ADCmax_chan_U[apv_cnt]->Draw();
		h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetXaxis()->UnZoom();
		h_APV_ratio_ADCmax_chan_U[apv_cnt]->SetTitle(Form("APV Channel Ratios - Run: %d, APV: %i)", runnum, apv_cnt));
		h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetXaxis()->SetTitle("ADC Ratio");
		h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetYaxis()->SetTitle("Entries");
		h_APV_ratio_ADCmax_chan_U[apv_cnt]->SetMarkerStyle(2);
		h_APV_ratio_ADCmax_chan_U[apv_cnt]->SetMarkerColor(06);
		h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetYaxis()->SetTitleOffset(1.5);
		// h_APV_ratio_ADCmax_chan_U[apv_cnt]->SetStats(0);
		
		TF1 *fitFcn = new TF1("fitFcn", fitFunction, 1, 15, 6);
		fitFcn->SetNpx(500);
		fitFcn->SetParameters(1, 1, 1, 1, 1, 1);
		h_APV_ratio_ADCmax_chan_U[apv_cnt]->Fit("fitFcn", "R");

		// fitFcn->SetParameters(13789.1, -4465.81, 523.945, 7684.37, 8.95072, 3.56006);

		c_APV_ratio_ADCmax_chan_U[apv_cnt]->Update();
		

		if(apv_cnt == first_apv){
			c_APV_ratio_ADCmax_chan_U[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/APV_chan_ratios_U_strips_all_%i_ADCcut_%i.pdf(", runnum, ADC_cut));
		}
		else if (apv_cnt > first_apv && apv_cnt < (last_apv-1)){
			c_APV_ratio_ADCmax_chan_U[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/APV_ratios_U_strips_all_%i_ADCcut_%i.pdf", runnum, ADC_cut));
		}

		if(apv_cnt == (last_apv - 1)){
			c_APV_ratio_ADCmax_chan_U[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/APV_ratios_U_strips_all_%i_ADCcut_%i.pdf)", runnum, ADC_cut));
		}
		// c_APV_ratio_ADCmax_chan_U[apv_cnt]->Close();
		// gSystem->ProcessEvents();

	}
	//END OF APVs

auto stop = high_resolution_clock::now();
auto duration = duration_cast<minutes>(stop - start);	
cout << "Time elapsed: " << duration.count() << " minutes." << endl;
}


	