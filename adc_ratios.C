#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>

#include "./include/include_files.h"
#include "./include/APV_strips.h"
#include "./include/search_file.C"

//Define some global variables
int const ADC_cut = 300;
int apv_chan_adc[128][2];
int apv_strip_adc[128][2];
Int_t max_strips = (3840);
Int_t max_events = 100;
Int_t nAPV_strips = 128;
Int_t gaus_min_bin;
Int_t gaus_max_bin;

//********BOOLEANS********//
bool build_APV_occupancy = false; 
bool build_APV_occupancy_U = false;
bool save_APV_occupancy_U = false;
bool build_ratio_all = true;
bool build_ratio_single = true;
bool fit_gaus = false;

TF1 *xtalk_gaus;
TF1 *xtalk_fit;


// for(int cut_cnt = 0; cut_cnt < (sizeof(db_cuts)/sizeof(db_cuts[0])); cut_cnt++){
// 	tpt_db_cuts->AddText(db_cuts[cut_cnt]);
// }

void adc_ratios(int runnum = 13560){ 
	TChain *TC = new TChain("T");
	// const char * DATA_DIR = Form("/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/%i/", runnum);
	const char * DATA_DIR = "/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/before_cuts/";
	const char * protorootfile = Form("e1209019_replayed_%i.root", runnum);
	// const char * protorootfile = Form("e1209019_fullreplay_%i_*", runnum);
	// const char * protorootfile = Form("e1209019_replayed_%i_%s.root", runnum, cut);
	// const char * protorootfile = "e1209019_fullreplay_13656_stream0_seg9_9.root";
	// const char * protorootfile = Form("*%i*%s*", runnum, cut);
	// const char * protorootfile = Form("*%i*", runnum);

	TString rootfile = Form("%s%s", DATA_DIR, protorootfile);
	cout << "Input rootfile is: " << rootfile << endl;
	TC->Add( rootfile );

	// Loop for Chain
	// const char * protorootfile;
	// for(int seg=0; seg <11; seg++){
	// 	protorootfile = Form("e1209019_fullreplay_11562_stream0_seg%i_%i.root", seg, seg);
	// 	TString rootfile = Form("%s%s", DATA_DIR, protorootfile);
	// 	cout << "Input rootfile is: " << rootfile << endl;
	// 	TC->Add( rootfile );
	// }
	

	//Cuts from db_bb.gem.dat
	// TString db_cuts[] = {
	// 	search_file("bb.gem.threshold_sample"),
	// 	search_file("bb.gem.threshold_stripsum"),
	// 	search_file("bb.gem.threshold_clustersum"),
	// 	search_file("bb.gem.ADCasym_cut"),
	// 	search_file("bb.gem.maxstrip_t0"),
	// 	search_file("bb.gem.maxstrip_tcut"),
	// 	search_file("bb.gem.addstrip_tcut"),
	// 	search_file("bb.gem.addstrip_ccor_cut"),
	// 	search_file("bb.gem.suppressfirstlast"),
	// 	search_file("bb.gem.deltat_cut"),
	// 	search_file("bb.gem.corrcoeff_cut")};
	// TPaveText *tpt_db_cuts = new TPaveText(0.4, 50.0, 0.65, 105.0);
	// TText *db_cuts_title = tpt_db_cuts->AddText("Tracking Cuts:");
	// db_cuts_title->SetTextFont(53);
	// db_cuts_title->SetTextSize(12);
	// for(int k = 0; k < (sizeof(db_cuts)/sizeof(db_cuts[0])); k++){
	// 		tpt_db_cuts->AddText(db_cuts[k]);
	// 	}
	
	// string searched_var = search_file("bb.gem.maxstrip_t0");

	//Turn off all Branches to save memory/time
	TC->SetBranchStatus("*", false);

	//Turn on desired branches
	TC->SetBranchStatus("fEvtHdr.fEvtNum", true);
	//Hit branches:
	TC->SetBranchStatus("bb.gem.hit.nstripu", true);
	TC->SetBranchStatus("bb.gem.hit.nstripv", true);
	TC->SetBranchStatus("bb.gem.hit.ADCmaxsampU", true);
	TC->SetBranchStatus("bb.gem.hit.ADCmaxsampV", true);

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


	///MAIN LOOP -- Where the fun happens!!etc
	//This fist loop is to step through APVs. We will look at individual APVs and this is
	//where that main loop starts.
	int nAPVs = 30;
	int occupancy_bin_num = 10000;
	int first_apv = 1;
	int last_apv = 10;
	

	//Grab some maximum values to use as limits and constraints
	double max_U_ADC = TC->GetMaximum("bb.gem.hit.ADCmaxsampU");
	double max_V_ADC = TC->GetMaximum("bb.gem.hit.ADCmaxsampV");

	//Define 1-D Histograms
	TH1D *h_APV_U_ADCmax_chan_ratios[nAPVs];


	//Define 2-D Histograms
	TH2D *h2_ADCmax = new TH2D("h2_ADCsamples", "", max_strips, 0, max_strips, 2, 0, max_strips);
	TH2D *h2_U_ADCmax = new TH2D("h2_U_ADCmax", "", max_strips, 0, max_strips, max_U_ADC+50, 0, max_U_ADC+50);

	//This is a general histogram that will hold all the histograms for all of the APVs.
	TH2D *h_apv_ADCmax_istrip[nAPVs];
	TH2D *h2_APV_U_ADCmax_apvstrip[nAPVs];
	TH2D *h2_APV_U_ADCmax_apvchan[nAPVs];


//************************************************************//
//******************      MAIN APV LOOP     ******************//	
//************************************************************//
	cout << endl << "Run number: " << runnum << endl;
	cout << endl << "Calculating total occupancy on layer." << endl;

	int total_event_num = 0;
	
	cout << "Analyzing data for APVs " << first_apv << " through " << last_apv << "." << endl;

	for(int apv_cnt = first_apv; apv_cnt < last_apv; apv_cnt++) {
		cout << "Processing data for apv: " << apv_cnt << endl;
		int min_APV_strip = APV_strip_nums(apv_cnt, "min");
		int max_APV_strip = APV_strip_nums(apv_cnt, "max");

		h_APV_U_ADCmax_chan_ratios[apv_cnt] = new TH1D(Form("h_APV%i_U_ADCmax_chan_ratios", apv_cnt), "", 300, 0, 30);
		h_apv_ADCmax_istrip[apv_cnt] = new TH2D(Form("h2_APV%i_U_ADCmax_istrip",apv_cnt), "",128, APV_strip_nums(apv_cnt, "min"), APV_strip_nums(apv_cnt, "max"), max_U_ADC, 0, max_U_ADC );
		h2_APV_U_ADCmax_apvstrip[apv_cnt] = new TH2D(Form("h2_APV%i_U_ADCmax_apvstrip", apv_cnt), "",128, 0, 128, max_U_ADC, 0, max_U_ADC );
		h2_APV_U_ADCmax_apvchan[apv_cnt] = new TH2D(Form("h2_APV%i_U_ADCmax_apvchan", apv_cnt), "", 128, 0, 128, max_U_ADC, 0, max_U_ADC );



		//We can either loop over all events or select some events:
		//Loop over all events:
		// int apv_event_num = 0;
		// while(TC->GetEntry(apv_event_num++)){
		// 	TC->GetEntry(apv_event_num);

		
		//Loop over selected events:
		for(int evt = 5; evt <= 60000; evt++){
			TC->GetEntry(evt);
			int occupancy_cnt = 0;

			//Only look at U-strips and loop over them. We only look up to the number of U-strips fired (Ndata_IsU)
			for(int i = 0; i < Ndata_IsU; i++){
				h2_ADCmax->Fill(istrip[i], ADC_max[i]);
				
				if(IsU[i]){
				h2_U_ADCmax->Fill(istrip[i], ADC_max[i]);

					//Only select strips that are in the selected range for the specific APV
					//Need to cross reference the global strip numbers and relate them to those of the APV_strip_nums
					
					if(istrip[i] >= min_APV_strip && istrip[i] <= max_APV_strip){
						//strip and ADC into each APV histogram
						h_apv_ADCmax_istrip[apv_cnt]->Fill(istrip[i], ADC_max[i]); 
						//Put the ADC value into the histograms for APV strips and Channels
						h2_APV_U_ADCmax_apvstrip[apv_cnt]->Fill( int(istrip[i])%128, ADC_max[i] );
						//In these we need to convert the strip nubmers to the multiplexer APV channels
						h2_APV_U_ADCmax_apvchan[apv_cnt]->Fill( UV_APV_strip_to_channel(int(istrip[i])%128), ADC_max[i]);
						//Fill an array with the maximum APV on each channel
						apv_chan_adc[UV_APV_strip_to_channel(int(istrip[i])%128)][0] = UV_APV_strip_to_channel(int(istrip[i])%128);
						apv_chan_adc[UV_APV_strip_to_channel(int(istrip[i])%128)][1] = ADC_max[i];

						apv_strip_adc[int(istrip[i])%128][0] = istrip[i];
						apv_strip_adc[int(istrip[i])%128][1] = ADC_max[i];

						occupancy_cnt++;
					}

				}
		
			}


			//The histograms for this event have been filled.
			//SCAN THROUGH FILLED HISTOGRAMS

			//First set some defaults
			double chan_n_adc = 0.0;
			double chan_n_plus_1_ADC = 0.0;
			double ADC_ratio = 0.0;
//***********************
///*******RATIOS*********
			if(build_ratio_all){
			//SCANNING THROUGH THE FILLED HISTOGRAMS
				for(int i = 0; i < 127; i++){
					double chan_n_adc = apv_chan_adc[i][1];
					double chan_n_plus_1_ADC = apv_chan_adc[i+1][1];
					double ADC_ratio = 0.0;

					//Always take larger ADC divided by smaller ADC for neighboring channels
					if(chan_n_adc > chan_n_plus_1_ADC && chan_n_plus_1_ADC != 0 && chan_n_adc > ADC_cut){
						double ADC_ratio = chan_n_adc/chan_n_plus_1_ADC;
						h_APV_U_ADCmax_chan_ratios[apv_cnt]->Fill(ADC_ratio);

					}
					else if(chan_n_plus_1_ADC > chan_n_adc && chan_n_adc != 0 && chan_n_plus_1_ADC > ADC_cut){
						double ADC_ratio = chan_n_plus_1_ADC/chan_n_adc;
						h_APV_U_ADCmax_chan_ratios[apv_cnt]->Fill(ADC_ratio);

					}

				}
			}
		//End of looping through events
		}

		/////////////////////////
		///////T-Graphs//////////
		/////////////////////////
		
		double apv_chan_adc_x[128];
		double apv_chan_adc_y[128];
		double apv_strip_adc_x[128];
		double apv_strip_adc_y[128];


		for(int i=0; i<128; i++){
			apv_chan_adc_x[i] = i;
			apv_chan_adc_y[i] = apv_chan_adc[i][1];

			apv_strip_adc_x[i] = i;
			apv_strip_adc_y[i] = apv_strip_adc[i][1];
		}
		// TCanvas *c_adcmax_apvchan_graph[nAPVs];
		// c_adcmax_apvchan_graph[apv_cnt] = new TCanvas(Form("ADCmax vs APV Channel - APV %i", apv_cnt), "c_adcmax_apvchan_13_graph", 700, 500);
		// TGraph* apv_chan_adc_graph[nAPVs];
		// apv_chan_adc_graph[apv_cnt] = new TGraph(128, apv_chan_adc_x, apv_chan_adc_y);
		// apv_chan_adc_graph[apv_cnt]->Draw("ALP");
		// apv_chan_adc_graph[apv_cnt]->SetTitle("ADCmax vs APV Channel - APV 13 - Ustrips");
		// apv_chan_adc_graph[apv_cnt]->GetXaxis()->SetTitle("APV Channel [n]");
		// apv_chan_adc_graph[apv_cnt]->GetYaxis()->SetTitle("ADC");
		// c_adcmax_apvchan_graph[apv_cnt]->Update();

		// TCanvas *c_adcmax_apvstrip_graph[nAPVs];
		// c_adcmax_apvstrip_graph[apv_cnt]  = new TCanvas(Form("ADCmax vs APV Strip - APV %i", apv_cnt), Form("c_adcmax_apvstrip_APV%i_graph", apv_cnt), 700, 500);
		// TGraph *apv_strip_adc_graph[nAPVs];
		// apv_strip_adc_graph[apv_cnt] = new TGraph(128, apv_strip_adc_x, apv_strip_adc_y); 
		// apv_strip_adc_graph[apv_cnt]->Draw("ALP");
		// apv_strip_adc_graph[apv_cnt]->SetTitle("ADCmax vs APV Strip - APV 13");
		// apv_strip_adc_graph[apv_cnt]->GetXaxis()->SetTitle("APV Strip [n]");
		// apv_strip_adc_graph[apv_cnt]->GetYaxis()->SetTitle("ADC");
		// c_adcmax_apvstrip_graph[apv_cnt]->Update();

		// cout << "size of graph: " << (sizeof(apv_chan_adc)/sizeof(apv_chan_adc[0])) << endl;


		/////////////////////////////////////////////
		//*******PLOT HISTOGRAMS**********histograms
		/////////////////////////////////////////////
		
		// TCanvas *c_adcmax_istrip_u = new TCanvas("ADCmax vs istrip [U-strips only]", "", 600, 500);
		// h2_U_ADCmax->Draw();
		// h2_U_ADCmax->SetTitle("ADCmax vs istrip [bb.gem.m0.strip.IsU only]");
		// h2_U_ADCmax->GetXaxis()->SetTitle("bb.gem.m0.strip.istrip [i]");
		// h2_U_ADCmax->GetYaxis()->SetTitle("bb.gem.m0.strip.ADCmax");
		// // h2_U_ADCmax->GetXaxis()->SetRange(2944, 3071);
		// h2_U_ADCmax->SetMarkerStyle(2);
		// h2_U_ADCmax->SetMarkerColor(06);
		// c_adcmax_istrip_u->Update();

		// TCanvas *c_h_apv_ADCmax_istrip[apv_cnt];
		// c_h_apv_ADCmax_istrip[apv_cnt] = new TCanvas("ADCmax vs strip - APV %i - Ustrips only", "c_adcmax_istrip_13", 600, 500);
		// h_apv_ADCmax_istrip[apv_cnt]->Draw();
		// h_apv_ADCmax_istrip[apv_cnt]->SetTitle("ADCmax vs strip - APV 13 - Ustrips only");
		// h_apv_ADCmax_istrip[apv_cnt]->GetXaxis()->SetTitle("istrips");
		// h_apv_ADCmax_istrip[apv_cnt]->GetYaxis()->SetTitle("ADCmax");
		// h_apv_ADCmax_istrip[apv_cnt]->SetMarkerStyle(2);
		// h_apv_ADCmax_istrip[apv_cnt]->SetMarkerColor(06);
		// c_h_apv_ADCmax_istrip[apv_cnt]->Update();

		// // Plot APV --> ADCmax vs Strips on the APV
		// TCanvas *c_adcmax_apvstrip[nAPVs];
		// c_adcmax_apvstrip[apv_cnt] = new TCanvas(Form("ADCmax vs strip - APV %i - APV strips", apv_cnt), Form("c_adcmax_apvstrip_apv%i", apv_cnt), 600, 500);

		// h2_APV_U_ADCmax_apvstrip[apv_cnt]->Draw("");
		// h2_APV_U_ADCmax_apvstrip[apv_cnt]->SetTitle("ADCmax vs APV strip - APV 13 - Ustrips only");
		// h2_APV_U_ADCmax_apvstrip[apv_cnt]->GetXaxis()->SetTitle("APV Strips");
		// h2_APV_U_ADCmax_apvstrip[apv_cnt]->GetYaxis()->SetTitle("ADCmax");
		// h2_APV_U_ADCmax_apvstrip[apv_cnt]->SetMarkerStyle(2);
		// c_adcmax_apvstrip[apv_cnt]->Update();

		// // Plot APV 13 --> ADCmax vs Channels on the APV
		// TCanvas *c_adcmax_apvchan[nAPVs];
		// c_adcmax_apvchan[apv_cnt] = new TCanvas(Form("ADCmax vs APV channel - APV %i - U strips", apv_cnt), Form("c_adcmax_apvchan_apv%i", apv_cnt), 700, 500);
		// h2_APV_U_ADCmax_apvchan[apv_cnt]->Draw("B");
		// h2_APV_U_ADCmax_apvchan[apv_cnt]->SetTitle("ADCmax vs APV channel - APV 13 - Ustrips only");
		// h2_APV_U_ADCmax_apvchan[apv_cnt]->GetXaxis()->SetTitle("APV channels");
		// h2_APV_U_ADCmax_apvchan[apv_cnt]->GetYaxis()->SetTitle("ADCmax");
		// h2_APV_U_ADCmax_apvchan[apv_cnt]->SetMarkerStyle(2);
		// c_adcmax_apvchan[apv_cnt]->Update();
		
		if(fit_gaus){
			h_APV_U_ADCmax_chan_ratios[apv_cnt]->GetXaxis()->SetRangeUser(1, 9);
			gaus_min_bin = (0.1)*(h_APV_U_ADCmax_chan_ratios[apv_cnt]->GetMinimumBin());
			h_APV_U_ADCmax_chan_ratios[apv_cnt]->GetXaxis()->SetRangeUser(gaus_min_bin, 15);
			gaus_max_bin = (0.1)*(h_APV_U_ADCmax_chan_ratios[apv_cnt]->GetMaximumBin());
			cout << "***********************************" << endl << endl;
			cout << "min bin: " << gaus_min_bin << "      maxb: " << gaus_max_bin << endl << endl;
			cout << "***********************************" << endl << endl;
			h_APV_U_ADCmax_chan_ratios[apv_cnt]->GetXaxis()->UnZoom();

			//*********FITTING CROSSTALK SIGNALS*********
			
			
		}

		if(build_ratio_all){
			// xtalk_gaus = new TF1("xtalk_gaus", "gaus", 1+gaus_min_bin, (2*gaus_max_bin) + gaus_min_bin -1);
			// h_APV_U_ADCmax_chan_ratios[apv_cnt]->Fit("xtalk_gaus", "R");
			// xtalk_fit = new TF1("xtalk_fit", "pol2", 1, gaus_min_bin);
			// h_APV_U_ADCmax_chan_ratios[apv_cnt]->Fit("xtalk_fit", "R");
	
			// TF1 total = new TF1("total_fit", "pol2(0) + gaus(3)", 1, (2*gaus_max_bin) + gaus_min_bin -1);
			// xtalk_fit->GetParameters(&par[0]);
			// xtalk_gaus->GetParameters(&par[3]);
			// total->SetParameters(par);
			// h_APV_U_ADCmax_chan_ratios[apv_cnt]->Fit(total, "R+");

			TCanvas *c_adcmax_ratios_apv[nAPVs];
			c_adcmax_ratios_apv[apv_cnt] = new TCanvas(Form("APV%i Ratio of U-Strip channel ADCmax", apv_cnt), Form("c_ratio_apv_%i", apv_cnt), 700, 500);
			h_APV_U_ADCmax_chan_ratios[apv_cnt]->Draw();
			// xtalk_gaus->Draw("same");
			// xtalk_fit->Draw("same");
			
			h_APV_U_ADCmax_chan_ratios[apv_cnt]->SetTitle(Form("Strip ratios (Run: %d, APV: %i, ADC cut: %i)", runnum, apv_cnt, ADC_cut));
			h_APV_U_ADCmax_chan_ratios[apv_cnt]->GetXaxis()->SetTitle("ADC Ratio");
			h_APV_U_ADCmax_chan_ratios[apv_cnt]->GetYaxis()->SetTitle("Entries");
			h_APV_U_ADCmax_chan_ratios[apv_cnt]->SetMarkerStyle(2);
			h_APV_U_ADCmax_chan_ratios[apv_cnt]->SetMarkerColor(06);
			c_adcmax_ratios_apv[apv_cnt]->Update();

			if(apv_cnt == first_apv){
				c_adcmax_ratios_apv[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/APV_ratios_U_strips_all_%i_ADCcut_%i.pdf(", runnum,ADC_cut));
			}
			else if (apv_cnt > first_apv && apv_cnt < (last_apv-1)){
				c_adcmax_ratios_apv[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/APV_ratios_U_strips_all_%i_ADCcut_%i.pdf", runnum, ADC_cut));
			}

			if(apv_cnt == (last_apv - 1)){
				c_adcmax_ratios_apv[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/APV_ratios_U_strips_all_%i_ADCcut_%i.pdf)", runnum, ADC_cut));
			}
			// c_adcmax_ratios_apv[apv_cnt]->Close();
			// gSystem->ProcessEvents();
		}
	
		
	//END OF APV_CNT LOOP
	}
	

///////////////////
//END OF MAIN LOOP
/////////////////////////////////////////////
}