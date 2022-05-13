#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
using namespace std::chrono;


#include "./include/include_files.h"
#include "./include/APV_strips.h"
#include "./include/search_file.C"
#include "./include/draw_plots.h"
#include "./include/save_outputs.h"

Double_t fitFunc(Double_t * x, Double_t *par){

	Double_t total_fit = 0.0;
	Double_t g1 = 0.0;
	Double_t g2 = 0.0;
	Double_t g3 = 0.0;
	
	g1 = exp(par[0] +(par[1]*x[0]));
	g2 = par[2]*exp((-0.5)*pow(((x[0] -  par[3])/par[4]),2));
	g3 = par[5]*exp((-0.5)*pow(((x[0] -  par[6])/par[7]),2));

	total_fit = g1 + g2 + g3;
	
	return total_fit;
}

TFile *outputRootFile = new TFile("apv_adc_ratios_histogram_99.root", "RECREATE");

//Define some global variables
int apv_chan_adc[128][2];
int apv_chan_adc_crosstalk_reject[128][2];
int apv_strip_adc[128][2];

double apv_adc_cut_thresh = 50.0;

int nAPVs = 30;
int occupancy_bin_num = 10000;
int first_apv = 0;
int last_apv = 30;
int ratio_bins = 300;
int nstrips_keep_cut;


//Set defaults for ADC on channel n and n+1
double chan_n_adc = 0.0;
double chan_n_plus_1_ADC = 0.0;
double ADC_ratio = 0.0;
double xtalk_mean;
double xtalk_sigma;

int nstrips = 3840;
Int_t max_strips = 8000;
Int_t nAPV_strips = 128;
Int_t gaus_min_bin;
Int_t gaus_max_bin;

double crosstalk_mean;
double crosstalk_integral;
double crosstalk_occupancy;
double crosstalk_occupancy_U;
Int_t Ndata_strip_ADC_samples;
Int_t Ndata_strip_ADC_max;
int Ndata_strip_ADC_max_keep;
int Ndata_strip_ADC_max_keepU;
int count_in_isu;


//BOOLEANS
bool build_all = false;
bool build_occupancies = true;
bool build_ratios = false;

void apv_adc_ratios(int const ADC_cut = 500){

	ofstream myfile;
	myfile.open("crosstalk_v_occupancy_U.txt");
	myfile << Form("ADCcut; %i", ADC_cut) << Form("; Noise_cut; %i", int(apv_adc_cut_thresh)) << endl;
	myfile << "run	occupancy	occupancy_U		crosstalk	gaus_mean	gaus_sigma	Ndata_ADCmax	Ndata_ADCmax_cut	Ndata_IsU	Ndata_IsU_cut		\n \n";

	auto start = high_resolution_clock::now();
	int runs[] = {11562};
	// int runs[] = {11449, 11451, 11456, 11494, 11562, 11580, 11595, 11997, 12001, 12013, 12030, 12050, 12060, 12073, 12342, 12423, 12424, 12425};
	// int runs[] = {11449, 11451, 11456, 11494, 11562, 11580, 11595, 11997, 12001, 12013, 12030, 12050, 12060, 12073, 12342, 12423, 12424, 12425, 12550, 12620, 12662, 12728, 13060, 13309, 13325, 13344, 13370, 13400, 13454, 13474, 13505, 13554, 13560, 13615, 13620, 13660, 13661, 13664, 13666, 13680, 13685, 13732, 13770, 13799};
	int num_runs = (sizeof(runs)/sizeof(runs[0]));

	double arr_xtalk_v_occupancies[num_runs][2];
	double arr_xtalk_v_occupancies_U[num_runs][2];

	double arr_occupancy[num_runs];
	double arr_occupancy_U[num_runs];
	double arr_crosstalk[num_runs];
	double arr_gaus_mean[num_runs];
	double arr_gaus_sigma[num_runs];
	double arr_Ndata_ADCmax[num_runs];
	double arr_Ndata_ADCmax_cut[num_runs];
	double arr_Ndata_IsU[num_runs];
	double arr_Ndata_IsU_cut[num_runs];

	//-------HISTOGRAMS FOR STORING INFO ACROSS RUNS Occupancies---------------
	TH1D *h_xtalk_v_occupancies_runs = new TH1D(Form("h_xtalk_v_occupancies_runs_%i_%i", runs[0], runs[num_runs-1]), "", occupancy_bin_num, 0, 1);

	//Beginning of RUNS Loop
	for(int irun = 0; irun < num_runs; irun++){
		int runnum = runs[irun];
		
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
		TCanvas *c_APV_ratio_source_mean4_ADCmax_chan_U[nAPVs];
		TCanvas *c_APV_ratio_source_mean9_ADCmax_chan_U[nAPVs];
		TCanvas *c_APV_ratio_source_upper_ratios_ADCmax_chan_U[nAPVs];

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
		TC->SetBranchAddress("Ndata.bb.gem.m0.strip.ADCsamples", &Ndata_strip_ADC_samples);
		TC->SetBranchAddress("Ndata.bb.gem.m0.strip.ADCmax", &Ndata_strip_ADC_max);
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
		
		
		//---------Total Occupancies--------------
		TH1D *h_total_occupancy[num_runs];
		h_total_occupancy[irun] = new TH1D(Form("h_total_occupancy_all_strips_%i", runs[irun]), "", occupancy_bin_num, 0, 1);
		TH1D *h_total_occupancy_U[num_runs];
		h_total_occupancy_U[irun] = new TH1D(Form("h_total_occupancy_U_%i", runs[irun]), "", occupancy_bin_num, 0, 1);

		//2-D histograms
		TH2D *h2_ADCmax_U = new TH2D("h2_ADCmax_samples_U", "" , nstrips, 0, nstrips, nstrips, 0, nstrips);

		TH1D *h_APV_ratio_ADCmax_chan_U[nAPVs];
		TH1D *h_APV_ratio_larger_ADC[nAPVs];
		TH1D *h_APV_ratio_smaller_ADC[nAPVs];
		TH2D *h2_APV_ratio_source_mean4_ADCmax_chan_U[nAPVs];
		TH2D *h2_APV_ratio_source_mean9_ADCmax_chan_U[nAPVs];
		TH2D *h2_APV_ratio_source_upper_ratios_ADCmax_chan_U[nAPVs];

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
			h_APV_ratio_larger_ADC[apv_cnt] = new TH1D(Form("h_APV%i_ratio_larger_adc", apv_cnt), "", ratio_bins, 0, nstrips);
			h_APV_ratio_smaller_ADC[apv_cnt] = new TH1D(Form("h_APV%i_ratio_smaller_adc", apv_cnt), "", 750, 0, 750);

		//Hisogram for Analyzing ratio contributions at different ratio values

			h2_APV_ratio_source_mean4_ADCmax_chan_U[apv_cnt] = new TH2D(Form("h2_APV%i_ratio_source_mean4_ADCmax_chan_U", apv_cnt), "", 750, 0, 750, ratio_bins, 0, nstrips);
			h2_APV_ratio_source_mean9_ADCmax_chan_U[apv_cnt] = new TH2D(Form("h2_APV%i_ratio_source_mean9_ADCmax_chan_U", apv_cnt), "", 750, 0, 750, ratio_bins, 0, nstrips);
			h2_APV_ratio_source_upper_ratios_ADCmax_chan_U[apv_cnt] = new TH2D(Form("h2_APV%i_ratio_source_upper_ratios_ADCmax_chan_U", apv_cnt), "", 750, 0, 750, ratio_bins, 0, nstrips);

		
	//EVENT //Loop all selected events:		
			// We can either loop over all events or select some events:
			// Loop over all events:

			int apv_event_num = 0;
			while(TC->GetEntry(apv_event_num++)){
				if(apv_event_num%10000 == 0){cout << "Analyzing event " << apv_event_num << " of " << TC->GetEntries() << " total events." << endl;}
				TC->GetEntry(apv_event_num);

	// // EVENT //Loop over selected events:

			// for(int evt = 0; evt <= 10; evt++){
			
			// 	// if(evt == 706 || evt == 1186) {continue;}
			// 	cout << "evt: " << evt << endl;
			// 	TC->GetEntry(evt);

		
	//----------Calculate Occupancies---------------
				int total_event_num = 0;
				double total_occupancy = 0.0;
				double total_occupancy_U = 0.0;
				

	//Scan through the strips and apply an extra ADC cut
	//We want only those that are passing apv_adc_cut_thresh
	//We do this for all strips and for U/V only
				// cout << "*****************************" << endl;
				// cout << "ndata_strip befoe cuts: " << Ndata_strip_ADC_max << endl;
				Ndata_strip_ADC_max_keep = 0;
				for(int istrip = 0; istrip < Ndata_strip_ADC_max ; istrip++){
					if(ADC_max[istrip] > apv_adc_cut_thresh){
						Ndata_strip_ADC_max_keep++;
					}
				}
				// cout << "ndata_strip after cuts: " << Ndata_strip_ADC_max_keep << endl;
				// cout << "*****************************" << endl;

	//We do this for all U strips now		
				// cout << "*****************************" << endl;
				// cout << "ndata_strip_u befoe cuts: " << Ndata_IsU << endl;

				Ndata_strip_ADC_max_keepU = 0;
				count_in_isu = 0;
				//------------U-strips-----------
				for(int i = 0; i < Ndata_IsU; i++){
					
					if(IsU[i]){
					count_in_isu++;
	//------------Threshold ADC cut on all apvs (Can maybe implement this on the pedestal RMS CUT???)
						if(ADC_max[i] > apv_adc_cut_thresh){
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

							Ndata_strip_ADC_max_keepU++;
						}
						
					}

				}
				// cout << "ndata_strip_u after cuts: " << Ndata_strip_ADC_max_keepU << endl;
				// cout << "count_in_isu: " << count_in_isu << endl;
				// cout << "*****************************" << endl;

	//Now calculated OCCUPANCIES
				// total_occupancy = nstrips_keep/(2.0*nstrips);
				// total_occupancy_U = nstrips_keepU/nstrips;
				total_occupancy = Ndata_strip_ADC_max_keep/(2.0*nstrips);
				total_occupancy_U = Ndata_strip_ADC_max_keepU/(1.0*nstrips);
				h_total_occupancy[irun]->Fill(total_occupancy);
				h_total_occupancy_U[irun]->Fill(total_occupancy_U);


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

						h_APV_ratio_larger_ADC[apv_cnt]->Fill(chan_n_adc);
						h_APV_ratio_smaller_ADC[apv_cnt]->Fill(chan_n_plus_1_ADC);
						
						//Ratio Sources --> Larger value on Y, smaller value on X
						//For ratios around a value of 4
						if( (ADC_ratio > 3) && (ADC_ratio < 5) ){
							h2_APV_ratio_source_mean4_ADCmax_chan_U[apv_cnt]->Fill(chan_n_plus_1_ADC, chan_n_adc);
						}

						//For ratios around a value of 10
						if( (ADC_ratio > 7.5) && (ADC_ratio < 11.5) ){
							h2_APV_ratio_source_mean9_ADCmax_chan_U[apv_cnt]->Fill(chan_n_plus_1_ADC, chan_n_adc);
						}

						if( (ADC_ratio > 17) ){
							h2_APV_ratio_source_upper_ratios_ADCmax_chan_U[apv_cnt]->Fill(chan_n_plus_1_ADC, chan_n_adc);
						}

					}
					else if(chan_n_plus_1_ADC > chan_n_adc && chan_n_adc != 0 && chan_n_plus_1_ADC > ADC_cut){
						double ADC_ratio = chan_n_plus_1_ADC/chan_n_adc;
						h_APV_ratio_ADCmax_chan_U[apv_cnt]->Fill(ADC_ratio);

						h_APV_ratio_larger_ADC[apv_cnt]->Fill(chan_n_plus_1_ADC);
						h_APV_ratio_smaller_ADC[apv_cnt]->Fill(chan_n_adc);
						
						//Ratio Sources --> Larger value on Y, smaller value on X
						//For ratios around a value of 4
						if( (ADC_ratio > 3) && (ADC_ratio < 5) ){
							h2_APV_ratio_source_mean4_ADCmax_chan_U[apv_cnt]->Fill(chan_n_adc, chan_n_plus_1_ADC);
						}
						
						//For ratios around a value of 10
						if( (ADC_ratio > 7.5) && (ADC_ratio < 11.5) ){
							h2_APV_ratio_source_mean9_ADCmax_chan_U[apv_cnt]->Fill(chan_n_adc, chan_n_plus_1_ADC);
						}

						if( (ADC_ratio > 17) ){
							h2_APV_ratio_source_upper_ratios_ADCmax_chan_U[apv_cnt]->Fill(chan_n_adc, chan_n_plus_1_ADC);
						}
					}
				}

			}
			//END of EVENTS
		
			//----------FITS------------//
			Double_t par[8];
			Double_t indep_par[8];
			Int_t lastXbin = (0.1)*h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetNbinsX();
			cout << "last x bin: " << lastXbin << endl;

			int gaus_last_point = 13;
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetXaxis()->SetRangeUser(1, 9);
			Double_t gaus_min_bin = (0.1)*(h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetMinimumBin());
			Double_t landau_max_bin = (1.25)*((0.1)*(h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetMinimumBin()));
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetXaxis()->SetRangeUser(gaus_min_bin, gaus_last_point);
			Double_t gaus1_max_val = h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetMaximum();
			Double_t gaus_max_bin = (0.1)*(h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetMaximumBin());
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetXaxis()->UnZoom();

			cout << "***********************************" << endl << endl;
			cout << "min bin: " << gaus_min_bin << "      maxb: " << gaus_max_bin << "   mean ADC: " << h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetRMS() << endl << endl;
			cout << "***********************************" << endl << endl;
			

			c_APV_ratio_ADCmax_chan_U[apv_cnt] = new TCanvas(Form("APV%i Ratio of Channels - Ustrips", apv_cnt), Form("c_ratio_chan_apv_%i", apv_cnt), 700, 500);
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->Draw();
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetXaxis()->UnZoom();
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->SetTitle(Form("APV Channel Ratios - Run: %d, APV: %i, ADCcut: %i, Noise cut: %i)", runnum, apv_cnt, ADC_cut, int(apv_adc_cut_thresh)));
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetXaxis()->SetTitle("ADC Ratio");
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetYaxis()->SetTitle("Entries");
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->SetMarkerStyle(2);
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->SetMarkerColor(06);
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetYaxis()->SetTitleOffset(1.5);
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->SetStats(0);


			TF1 *myFitFunc = new TF1("myFitFunc", fitFunc, 1, 30, 8);
			
			myFitFunc->SetLineWidth(3);
			myFitFunc->SetLineColor(2);
			myFitFunc->SetLineStyle(9);

			myFitFunc->SetParName(0, "Exp Const");
			myFitFunc->SetParName(1, "Exp Slope");
			myFitFunc->SetParName(2, "G1 Norm");
			myFitFunc->SetParName(3, "G1 Center");
			myFitFunc->SetParName(4, "G1 Sigma");
			myFitFunc->SetParName(5, "G2 Norm");
			myFitFunc->SetParName(6, "G2 Center");
			myFitFunc->SetParName(7, "G2 Sigma");
			
			// //Initial parameters for pol2 - gaus1 - gaus1
			// // myFitFunc->SetParameters(13789.1, -4465.81, 523.945, 8000.0, 4.0, 8.96587, 1.16e+06, -148.77, 47.73);

			// Initial parameters for expo - gaus1 - gaus2
			// Good parameters for ADCcut 300 & 500:
			// myFitFunc->SetParameters(9.30973, -2.26944e-01, 8000, 4.0, 8.96587, 1.16e+06, 5, 47.73);

			// myFitFunc->SetParLimits(0, 9.25, 10.75);
			// myFitFunc->SetParLimits(1, -12.10, -1.15);
			// myFitFunc->SetParLimits(2, 0, gaus1_max_val);
			// myFitFunc->SetParLimits(3, 5, 15);
			// myFitFunc->SetParLimits(4, 1, 3);
			// myFitFunc->SetParLimits(5, 0, 100000);
			// myFitFunc->SetParLimits(6, 2, 8);
			// myFitFunc->SetParLimits(7, 0, 100000);

			myFitFunc->SetParameters(9.30973, -2.26944e-01, 8000, 4.0, 8.96587, 1.16e+06, 5, 47.73);
			// myFitFunc->SetParLimits(0, 9.25, 10.75);
			myFitFunc->SetParLimits(0, 9.25, 13.74);
			myFitFunc->SetParLimits(1, -12.10, -0.50);
			myFitFunc->SetParLimits(2, 0, gaus1_max_val);
			myFitFunc->SetParLimits(3, 5, 12);
			myFitFunc->SetParLimits(4, .5, 2.5);
			myFitFunc->SetParLimits(5, 0, 4000);
			myFitFunc->SetParLimits(6, 18, 24);
			myFitFunc->SetParLimits(7, 0, 7);
			
			// //ADCcut 200
			// myFitFunc->SetParameters(9.30973, -2.26944e-01, 8000, 4.0, 8.96587, 1.16e+06, 5, 47.73);

			// myFitFunc->SetParLimits(0, 9.25, 12);
			// myFitFunc->SetParLimits(1, -12.10, -1.15);
			// myFitFunc->SetParLimits(2, 0, 20000);
			// myFitFunc->SetParLimits(3, 3, 4.1);
			// myFitFunc->SetParLimits(4, 0, 3);
			// myFitFunc->SetParLimits(5, 0, 1.0e+8);
			// myFitFunc->SetParLimits(6, 7, 12);
			// myFitFunc->SetParLimits(7, 4, 1000);
			
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->Fit("myFitFunc", "0");

			myFitFunc->GetParameters(&par[0]);
			Double_t fit_chi = myFitFunc->GetChisquare();

			for(int i=0; i < 10; i++){
				h_APV_ratio_ADCmax_chan_U[apv_cnt]->Fit("myFitFunc", "eQ0");
			
				myFitFunc->GetParameters(&par[0]);
				myFitFunc->SetParameters(par);
				fit_chi = myFitFunc->GetChisquare();
			}
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->Fit("myFitFunc", "R+");

			Double_t fit_integral = myFitFunc->Integral(1, 30);
			cout << "max bin content: " << gaus1_max_val << endl;
			cout << "***********************************" << endl << endl;
			cout << "Chi2 of total fit: " << fit_chi << endl << endl;
			cout << "***********************************" << endl << endl;
			cout << "Integral of fit on [1, 30] = " << fit_integral << endl << endl;
			cout << "***********************************" << endl << endl;
			myFitFunc->Draw("same");

			Double_t fit_x[300], expo_y[300], gaus1_y[300], gaus2_y[300];
			Int_t fit_n = 300;
			for(Int_t i=0; i <= fit_n; i++){
				fit_x[i]=i*0.1+1;
				expo_y[i] = exp(par[0] +(par[1]*fit_x[i]));
				gaus1_y[i] = par[2]*exp((-0.5)*pow(((fit_x[i] -  par[3])/par[4]),2));
				gaus2_y[i] = par[5]*exp((-0.5)*pow(((fit_x[i] -  par[6])/par[7]),2));

			}
		
			TGraph* gr_expo = new TGraph(fit_n, fit_x, expo_y);
			gr_expo->SetLineColor(3);
			gr_expo->SetLineStyle(10);
			TGraph* gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
			gr_gaus1->SetLineColor(6);
			gr_gaus1->SetLineStyle(7);
			TGraph* gr_gaus2 = new TGraph(fit_n, fit_x, gaus2_y);
			gr_gaus2->SetLineColor(7);
			gr_gaus2->SetLineStyle(2);
			
			
			xtalk_mean = par[3];
			xtalk_sigma = par[4];

			gr_expo->Draw("same");
			gr_gaus1->Draw("same");
			gr_gaus2->Draw("same");

			
			crosstalk_mean = gr_gaus1->GetMean();
			
			crosstalk_integral = gr_gaus1->Integral(1, 30);
			crosstalk_occupancy_U = crosstalk_integral/nstrips;
			crosstalk_occupancy = crosstalk_integral/(2.0*nstrips);

			cout << "crosstalk integral U: " << crosstalk_integral << endl;
			
			// TPaveLabel *tpl_integral = new TPaveLabel(19, 9000, 28, 11000, Form("Integral of fit: %i", int(fit_integral)));
			// tpl_integral->SetBorderSize(1);
			// tpl_integral->Draw("same");

			TLegend *legend = new TLegend(0.6,0.65,0.88,0.85);
			legend->AddEntry(h_APV_ratio_ADCmax_chan_U[apv_cnt], "ADC Ratios", "l");
			legend->AddEntry(myFitFunc, "Total Fit", "l");
			legend->AddEntry(gr_expo, "expo fit", "l");
			legend->AddEntry(gr_gaus1, "gaus fit #1", "l");
			legend->AddEntry(gr_gaus2, "gaus fit #2", "l");
			legend->AddEntry(myFitFunc, Form("Fit Chi-2: %f", fit_chi), "l");
			legend->Draw("same");
			h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetYaxis()->SetRangeUser(0, h_APV_ratio_ADCmax_chan_U[apv_cnt]->GetMaximum());
			c_APV_ratio_ADCmax_chan_U[apv_cnt]->Update();

			outputRootFile->WriteObject(h_APV_ratio_ADCmax_chan_U[apv_cnt], Form("h_APV%i_ratio_ADCmax_chan_U_run%i", apv_cnt, runs[irun]));
			
			print_single_PDF(c_APV_ratio_ADCmax_chan_U[apv_cnt], Form("multirun_APV_channel_ratios_U_ADCcut_%i_NoiseCut_%i", ADC_cut, int(apv_adc_cut_thresh)), Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/", ADC_cut), first_apv, apv_cnt, last_apv);

			// if(irun == 0){
			// 	c_APV_ratio_ADCmax_chan_U[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/multirun_APV%i_channel_ratios_U_ADCcut_%i_NoiseCut_%i.pdf(", ADC_cut, apv_cnt, ADC_cut, int(apv_adc_cut_thresh)));
			// }
			// else if (irun > 0 && irun < (num_runs-1)){
			// 	c_APV_ratio_ADCmax_chan_U[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/multirun_APV%i_channel_ratios_U_ADCcut_%i_NoiseCut_%i.pdf", ADC_cut, apv_cnt, ADC_cut, int(apv_adc_cut_thresh)));
			// }

			// if(irun == (num_runs - 1)){
			// 	c_APV_ratio_ADCmax_chan_U[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/multirun_APV%i_channel_ratios_U_ADCcut_%i_NoiseCut_%i.pdf)", ADC_cut, apv_cnt, ADC_cut, int(apv_adc_cut_thresh)));
			// }
			// // c_APV_ratio_ADCmax_chan_U[apv_cnt]->Close();
			// gSystem->ProcessEvents();
			
			

			// if(apv_cnt == first_apv){
			// 	c_APV_ratio_ADCmax_chan_U[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/APV_chan_ratios_U_strips_all_%i_ADCcut_%i.pdf(", runnum, ADC_cut));
			// }
			// else if (apv_cnt > first_apv && apv_cnt < (last_apv-1)){
			// 	c_APV_ratio_ADCmax_chan_U[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/APV_ratios_U_strips_all_%i_ADCcut_%i.pdf", runnum, ADC_cut));
			// }

			// if(apv_cnt == (last_apv - 1)){
			// 	c_APV_ratio_ADCmax_chan_U[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/APV_ratios_U_strips_all_%i_ADCcut_%i.pdf)", runnum, ADC_cut));
			// }
			// c_APV_ratio_ADCmax_chan_U[apv_cnt]->Close();
			// gSystem->ProcessEvents();

		//------------------------------------------------------------------
		//------------------------------------------------------------------
		//------------Ratios and ratio-related plots

			if(true){

				TH1D *h_px_ratio_mean4[nAPVs];
				TH1D *h_py_ratio_mean4[nAPVs];
				TH1D *h_px_ratio_mean9[nAPVs];
				TH1D *h_py_ratio_mean9[nAPVs];
				TH1D *h_px_ratio_upper_ratios[nAPVs];
				TH1D *h_py_ratio_upper_ratios[nAPVs];

				h_px_ratio_mean4[apv_cnt] = h2_APV_ratio_source_mean4_ADCmax_chan_U[apv_cnt]->ProjectionX("h_px_ratio_mean4", 0, -1);
				h_py_ratio_mean4[apv_cnt] = h2_APV_ratio_source_mean4_ADCmax_chan_U[apv_cnt]->ProjectionY("h_py_ratio_mean4", 0, -1);

				h_px_ratio_mean9[apv_cnt] = h2_APV_ratio_source_mean9_ADCmax_chan_U[apv_cnt]->ProjectionX("h_px_ratio_mean9",0, -1);
				h_px_ratio_mean9[apv_cnt]->GetXaxis()->SetRange(0, 500);
				h_py_ratio_mean9[apv_cnt] = h2_APV_ratio_source_mean9_ADCmax_chan_U[apv_cnt]->ProjectionY("h_py_ratio_mean9",0, -1);

				h_px_ratio_upper_ratios[apv_cnt] = h2_APV_ratio_source_upper_ratios_ADCmax_chan_U[apv_cnt]->ProjectionX("h_px_ratio_upper_ratios", 0, -1);
				h_px_ratio_upper_ratios[apv_cnt]->GetXaxis()->SetRange(0, 300);
				h_py_ratio_upper_ratios[apv_cnt] = h2_APV_ratio_source_upper_ratios_ADCmax_chan_U[apv_cnt]->ProjectionY("h_py_ratio_upper_ratios", 0, -1);


				//Plotting Ratio sources
				bool first = false;
				bool last = false;
				if(apv_cnt == first_apv){first = true;}
				else if(apv_cnt == last_apv-1){last = true;}

				cout << "first: " << first << "  last: " << last << "  apv_cnt: " << apv_cnt << endl;

				outputRootFile->WriteObject(h2_APV_ratio_source_mean4_ADCmax_chan_U[apv_cnt], Form("h2_APV%i_ratio_source_mean4_ADCmax_chan_U", apv_cnt));

				TCanvas *c_APV_ratio_source_mean4_ADCmax_chan_U = plot_2DH(h2_APV_ratio_source_mean4_ADCmax_chan_U[apv_cnt], Form("c_APV%i_ratio_source_mean4_ADCmax_chan_U", apv_cnt), Form("APV%i Ratio - Large ADC vs Smaller ADC for Ratios = 3 thru 5, Run: %i, ADCcut = %i", apv_cnt, runnum, ADC_cut), "ADC of Denominator", "ADC of Numerator", "colz");
				print_multi_PDF(c_APV_ratio_source_mean4_ADCmax_chan_U, Form("APV_ratio_sources_various_region_plots_ADCcut_%i_NoiseCut_%i", ADC_cut, int(apv_adc_cut_thresh)), Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/", ADC_cut), first, false);

				outputRootFile->WriteObject(h2_APV_ratio_source_mean9_ADCmax_chan_U[apv_cnt], Form("h2_APV%i_ratio_source_mean9_ADCmax_chan_U", apv_cnt));

				TCanvas *c_APV_ratio_source_mean9_ADCmax_chan_U = plot_2DH(h2_APV_ratio_source_mean9_ADCmax_chan_U[apv_cnt], Form("c_APV%i_ratio_source_mean9_ADCmax_chan_U", apv_cnt), Form("APV%i Ratio - Large ADC vs Smaller ADC for Ratios = 7.5 thru 11.5, Run: %i, ADCcut = %i", apv_cnt, runnum, ADC_cut), "ADC of Denominator", "ADC of Numerator", "colz");
				print_multi_PDF(c_APV_ratio_source_mean9_ADCmax_chan_U, Form("APV_ratio_sources_various_region_plots_ADCcut_%i_NoiseCut_%i", ADC_cut, int(apv_adc_cut_thresh)), Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/", ADC_cut), false, false);

				outputRootFile->WriteObject(h2_APV_ratio_source_upper_ratios_ADCmax_chan_U[apv_cnt], Form("h2_APV%i_ratio_source_upper_ratios_ADCmax_chan_U", apv_cnt));

				TCanvas *c_APV_ratio_source_upper_ratios_ADCmax_chan_U = plot_2DH(h2_APV_ratio_source_upper_ratios_ADCmax_chan_U[apv_cnt], Form("c_APV%i_ratio_source_upper_ratios_ADCmax_chan_U", apv_cnt), Form("APV%i Ratio - Large ADC vs Smaller ADC for Ratios Greater Than 17, Run: %i, ADCcut = %i, Noise cut = %i", apv_cnt, runnum, ADC_cut, int(apv_adc_cut_thresh)), "ADC of Denominator", "ADC of Numerator", "colz");
				print_multi_PDF(c_APV_ratio_source_upper_ratios_ADCmax_chan_U, Form("APV_ratio_sources_various_region_plots_ADCcut_%i_NoiseCut_%i", ADC_cut, int(apv_adc_cut_thresh)), Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/", ADC_cut), false, false);

				outputRootFile->WriteObject(h_px_ratio_mean4[apv_cnt], Form("h2_APV%i_px_ratio_source_mean4_ADCmax_chan_U", apv_cnt));

				TCanvas *c_px_ratio_mean4 = plot_1DH(h_px_ratio_mean4[apv_cnt], "c_px_ratio_mean4", Form("APV%i Proj X: Ratio - Large ADC vs Smaller ADC for Ratios = 3 thru 5, Run %i, ADCcut = %i, Noise cut = %i", apv_cnt, runnum, ADC_cut, int(apv_adc_cut_thresh)), "ADC", "Entries", "colz");
				print_multi_PDF(c_px_ratio_mean4, Form("APV_ratio_sources_various_region_plots_ADCcut_%i_NoiseCut_%i", ADC_cut, int(apv_adc_cut_thresh)), Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/", ADC_cut), false, false);
				
				outputRootFile->WriteObject(h_py_ratio_mean4[apv_cnt], Form("h2_APV%i_py_ratio_source_mean4_ADCmax_chan_U", apv_cnt));
				
				TCanvas *c_py_ratio_mean4 = plot_1DH(h_py_ratio_mean4[apv_cnt], "c_py_ratio_mean4", Form("APV%i Proj Y: Ratio - Large ADC vs Smaller ADC for Ratios = 3 thru 5, Run: %i, ADCcut = %i, Noise cut = %i", apv_cnt, runnum, ADC_cut, int(apv_adc_cut_thresh)), "ADC", "Entries", "colz");
				print_multi_PDF(c_py_ratio_mean4, Form("APV_ratio_sources_various_region_plots_ADCcut_%i_NoiseCut_%i", ADC_cut, int(apv_adc_cut_thresh)), Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/", ADC_cut), false, false);

				outputRootFile->WriteObject(h_px_ratio_mean9[apv_cnt], Form("h2_APV%i_px_ratio_source_mean9_ADCmax_chan_U", apv_cnt));
				
				TCanvas *c_px_ratio_mean9 = plot_1DH(h_px_ratio_mean9[apv_cnt], "c_px_ratio_mean9", Form("APV%i Proj X: Ratio - Large ADC vs Smaller ADC for Ratios = 7.5 thru 11.5, Run: %i, ADCcut = %i, Noise cut = %i", apv_cnt, runnum, ADC_cut, int(apv_adc_cut_thresh)), "ADC", "Entries)", "colz");
				print_multi_PDF(c_px_ratio_mean9, Form("APV_ratio_sources_various_region_plots_ADCcut_%i_NoiseCut_%i", ADC_cut, int(apv_adc_cut_thresh)), Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/", ADC_cut), false, false);

				outputRootFile->WriteObject(h_py_ratio_mean9[apv_cnt], Form("h2_APV%i_py_ratio_source_mean9_ADCmax_chan_U", apv_cnt));

				TCanvas *c_py_ratio_mean9 = plot_1DH(h_py_ratio_mean9[apv_cnt], "c_py_ratio_mean9", Form("APV%i Proj Y: Ratio - Large ADC vs Smaller ADC for Ratios = 7.5 thru 11.5, Run: %i, ADCcut = %i, Noise cut = %i", apv_cnt, runnum, ADC_cut, int(apv_adc_cut_thresh)), "ADC", "Entries)", "colz");
				print_multi_PDF(c_py_ratio_mean9, Form("APV_ratio_sources_various_region_plots_ADCcut_%i_NoiseCut_%i", ADC_cut, int(apv_adc_cut_thresh)), Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/", ADC_cut), false, false);

				outputRootFile->WriteObject(h_px_ratio_upper_ratios[apv_cnt], Form("h2_APV%i_px_ratio_source_upper_ratios_ADCmax_chan_U", apv_cnt));

				TCanvas *c_px_ratio_upper_ratios = plot_1DH(h_px_ratio_upper_ratios[apv_cnt], "c_px_ratio_upper_ratios",Form("APV%i Proj X: Ratio - Large ADC vs Smaller ADC for Ratios Greater Than 17, Run %i, ADCcut = %i, Noise cut = %i", apv_cnt, runnum, ADC_cut, int(apv_adc_cut_thresh)), "ADC", "Entries", "colz");
				print_multi_PDF(c_px_ratio_upper_ratios, Form("APV_ratio_sources_various_region_plots_ADCcut_%i_NoiseCut_%i", ADC_cut, int(apv_adc_cut_thresh)), Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/", ADC_cut), false, false);

				outputRootFile->WriteObject(h_py_ratio_upper_ratios[apv_cnt], Form("h2_APV%i_py_ratio_source_upper_ratios_ADCmax_chan_U", apv_cnt));

				TCanvas *c_py_ratio_upper_ratios = plot_1DH(h_py_ratio_upper_ratios[apv_cnt], "c_py_ratio_upper_ratios",Form("APV%i Proj Y: Ratio - Large ADC vs Smaller ADC for Ratios Greater Than 17, Run %i, ADCcut = %i, Noise cut = %i", apv_cnt, runnum, ADC_cut, int(apv_adc_cut_thresh)), "ADC", "Entries", "colz");
				print_multi_PDF(c_py_ratio_upper_ratios, Form("APV_ratio_sources_various_region_plots_ADCcut_%i_NoiseCut_%i", ADC_cut, int(apv_adc_cut_thresh)), Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/ratio/ADCcut_%i/", ADC_cut), false, last);
			}
			outputRootFile->WriteObject(h_APV_ratio_larger_ADC[apv_cnt], Form("h_APV%i_ratio_larger_ADC", apv_cnt));
			plot_1DH(h_APV_ratio_larger_ADC[apv_cnt], "c_APV_ratio_larger", Form("Larger ADC of Ratio Calc, Run: %i, ADCcut = %i", runnum, ADC_cut), "ADC", "Entries");

			outputRootFile->WriteObject(h_APV_ratio_smaller_ADC[apv_cnt], Form("h_APV%i_ratio_smaller_ADC", apv_cnt));
			plot_1DH(h_APV_ratio_smaller_ADC[apv_cnt], "c_APV_ratio_smaller", Form("Smaller ADC of Ratio Calc, Run: %i, ADCcut = %i", runnum, ADC_cut), "ADC", "Entries");

		}
		//END OF APVs
		double final_total_occupancy_U = h_total_occupancy_U[irun]->GetMean();
		arr_xtalk_v_occupancies[irun][0] = h_total_occupancy[irun]->GetMean();
		arr_xtalk_v_occupancies[irun][1] = crosstalk_occupancy;
		arr_xtalk_v_occupancies_U[irun][0] = h_total_occupancy_U[irun]->GetMean();
		arr_xtalk_v_occupancies_U[irun][1] = crosstalk_occupancy_U;

		arr_occupancy[irun] = h_total_occupancy[irun]->GetMean();
		arr_occupancy_U[irun] = h_total_occupancy_U[irun]->GetMean();
		arr_gaus_mean[irun] = xtalk_mean;
		arr_gaus_sigma[irun] = xtalk_sigma;
		arr_Ndata_ADCmax[irun] = Ndata_strip_ADC_max;
		arr_Ndata_ADCmax_cut[irun] = Ndata_strip_ADC_max_keep;
		arr_Ndata_IsU[irun] = count_in_isu;
		arr_Ndata_IsU_cut[irun] = Ndata_strip_ADC_max_keepU;

		outputRootFile->WriteObject(h_total_occupancy[irun], Form("h_total_occupancy_%i", runs[irun]));
		plot_1DH(h_total_occupancy[irun], "c_total_occupancy", Form("Total Occupancy on All Strips, Run: %i, ADCcut = %i, Noise cut = %i", runnum, ADC_cut, int(apv_adc_cut_thresh)), "Occupancy", "Entries");

		outputRootFile->WriteObject(h_total_occupancy_U[irun], Form("h_total_occupancy_U_%i", runs[irun]));
		plot_1DH(h_total_occupancy_U[irun], "c_total_occupancy_U", Form("Total Occupancy on All U Strips, Run: %i, ADCcut = %i, Noise cut = %i", runnum, ADC_cut, int(apv_adc_cut_thresh)), "Occupancy on U Strips", "Entries");
		
		cout << "**************************************************" << endl;
		cout << "Occupancy on U strips: " << arr_xtalk_v_occupancies_U[0][0] << "   Xtalk Occupancy (aginst U): " << arr_xtalk_v_occupancies_U[0][1] << "  Xtalk occupancy (Against All strips): " << crosstalk_occupancy << endl;
		cout << "Occupancy on U mean: " << final_total_occupancy_U << endl;
		cout << "Crosstalk Gaussian mean: " << xtalk_mean << "  sigma: " << xtalk_sigma << endl;

		myfile << runs[irun] << "; " << arr_xtalk_v_occupancies[irun][0] << "; " << arr_xtalk_v_occupancies_U[irun][0] << "; " << arr_xtalk_v_occupancies_U[irun][1] << "; " << xtalk_mean << "; " << xtalk_sigma << "; " << Ndata_strip_ADC_max << "; " << Ndata_strip_ADC_max_keep << "; " << count_in_isu << "; " << Ndata_strip_ADC_max_keepU << endl;
		cout << "**************************************************" << endl;
		
	}
	//END of RUNS Loop
	
//print out the crosstalk array stuff:


TH1D *h_xtalk_v_occupancy = new TH1D("h_xtalk_v_occupancy", "", 1000, 0, 1);

for(int i = 0; i< num_runs; i++){
	h_xtalk_v_occupancy->SetBinContent(1000*arr_xtalk_v_occupancies_U[i][0], arr_xtalk_v_occupancies_U[i][1]);
}

TCanvas *c_xtalk_v_occupancy = new TCanvas("c_xtalk_v_occupancy", "", 700, 500);
h_xtalk_v_occupancy->Draw();
h_xtalk_v_occupancy->SetTitle(Form("Crosstalk Occupancy vs Total Cccupancy on U-Strips, ADCcut = %i, Noise Cut = %i", ADC_cut, int(apv_adc_cut_thresh)));
h_xtalk_v_occupancy->GetXaxis()->SetTitle("Total Occupancy on U-Strips");
h_xtalk_v_occupancy->GetYaxis()->SetTitle("Crosstalk Occupancy");
c_xtalk_v_occupancy->Update();

// plot_1DH(h_xtalk_v_occupancy, "xtalk_v_occupancy", Form("Crosstalk Occupancy vs Total Cccupancy on U-Strips, ADCcut = %i, Noise Cut = %i", ADC_cut, int(apv_adc_cut_thresh)), "Total Occupancy on U-Strips", "Crosstalk Occupancy", {0.0, 0.0}, {0.0, 0.0});

auto stop = high_resolution_clock::now();
auto duration = duration_cast<minutes>(stop - start);	
cout << "Time elapsed: " << duration.count() << " minutes." << endl;
}


	