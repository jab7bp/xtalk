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
int apv_chan_adc[128][2];
int apv_strip_adc[128][2];
int const ADC_cut = 600;
Int_t max_strips = (3840);
Int_t max_events = 100;
Int_t nAPV_strips = 128;

int single_apv = 13;

//********BOOLEANS********//
bool build_APV_occupancy = true;
bool plot_APV_occupancy = true; 
bool save_APV_occupancy = true;
bool build_APV_occupancy_U = true;
bool plot_APV_occupancy_U = true;
bool save_APV_occupancy_U = true;

bool build_total_occupancy = true;
bool plot_total_occupancy = true;
bool save_total_occupancy = true;
bool build total_occupancy_U = true;
bool plot_total_occupancy_U = true;
bool save_total_occupancy_U = true;

bool build_ratio_single = true;
bool plot_ratio_single = true;
bool save_ratio_single = true;

bool build_ratio_all = true;
bool plot_ratio_all = true;
bool save_ratio_all = true;

// for(int cut_cnt = 0; cut_cnt < (sizeof(db_cuts)/sizeof(db_cuts[0])); cut_cnt++){
// 	tpt_db_cuts->AddText(db_cuts[cut_cnt]);
// }

void xtalk(int runnum = 13444, const char *cut = "1A"){ 
	
	const char * DATA_DIR = "/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/";
	// const char * protorootfile = Form("e1209019_replayed_%i.root", runnum);
	// const char * protorootfile = Form("e1209019_fullreplay_%i_*", runnum);
	const char * protorootfile = Form("e1209019_replayed_%i_%s.root", runnum, cut);
	// const char * protorootfile = "e1209019_fullreplay_13656_stream0_seg9_9.root";
	// const char * protorootfile = Form("*%i*%s*", runnum, cut);

	TString rootfile = Form("%s%s", DATA_DIR, protorootfile);
	cout << "Input rootfile is: " << rootfile << endl;

	//Loop for Chain
	

	TChain *TC = new TChain("T");
	TC->Add( rootfile );

	//Define an output file to print things protorootfile
	ofstream myfile;
	myfile.open("output.txt");
	myfile << "occupancy: \n \n";
	myfile << "-------------------------" << endl;

	//Cuts from db_bb.gem.dat
	TString db_cuts[] = {
		search_file("bb.gem.threshold_sample"),
		search_file("bb.gem.threshold_stripsum"),
		search_file("bb.gem.threshold_clustersum"),
		search_file("bb.gem.ADCasym_cut"),
		search_file("bb.gem.maxstrip_t0"),
		search_file("bb.gem.maxstrip_tcut"),
		search_file("bb.gem.addstrip_tcut"),
		search_file("bb.gem.addstrip_ccor_cut"),
		search_file("bb.gem.suppressfirstlast"),
		search_file("bb.gem.deltat_cut"),
		search_file("bb.gem.corrcoeff_cut")};
	TPaveText *tpt_db_cuts = new TPaveText(0.4, 50.0, 0.65, 105.0);
	TText *db_cuts_title = tpt_db_cuts->AddText("Tracking Cuts:");
	db_cuts_title->SetTextFont(53);
	db_cuts_title->SetTextSize(12);
	for(int k = 0; k < (sizeof(db_cuts)/sizeof(db_cuts[0])); k++){
			tpt_db_cuts->AddText(db_cuts[k]);
		}
	
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
	int first_apv = 0;
	int last_apv = nAPVs;
	
	double apv_occupancy_mean[nAPVs];
	//LABELS and Lines
	TPaveLabel *tpl_apv_occupancy_mean_U[nAPVs];
	TLine *tl_apv_occupancy_lines[nAPVs];

	//Histogram with all individual APV occupancies
	TH1D *h_all_APV_occupancies_U = new TH1D("h_all_APV_occupancies_U", "", (occupancy_bin_num*nAPVs), 0.0, double(nAPVs));

	//Grab some maximum values to use as limits and constraints
	int first_occupancy_bin = first_apv*occupancy_bin_num;
	double max_U_ADC = TC->GetMaximum("bb.gem.hit.ADCmaxsampU");
	double max_V_ADC = TC->GetMaximum("bb.gem.hit.ADCmaxsampV");

	//Define 1-D Histograms
	TH1D *h_APV_U_ADCmax_chan_ratios[nAPVs];
	TH1D *h_APV_occupancy_U[nAPVs];

	//Define 2-D Histograms
	TH2D *h2_ADCmax = new TH2D("h2_ADCsamples", "", max_strips, 0, max_strips, 2, 0, max_strips);
	TH2D *h2_U_ADCmax = new TH2D("h2_U_ADCmax", "", max_strips, 0, max_strips, max_U_ADC+50, 0, max_U_ADC+50);

	//This is a general histogram that will hold all the histograms for all of the APVs.
	TH2D *h_apv_ADCmax_istrip[nAPVs];
	TH2D *h2_APV_U_ADCmax_apvstrip[nAPVs];
	TH2D *h2_APV_U_ADCmax_apvchan[nAPVs];
	TH2D *h2_APV_U_ADCmax_neigh_chan[nAPVs];
	TH2D *h2_APV_U_ADCmax_neigh_chan_lrg_sml[nAPVs];

//************************************************************//
//******************      MAIN APV LOOP     ******************//	
//************************************************************//
	cout << endl << "Run number: " << runnum << endl;
	cout << endl << "Calculating total occupancy on layer." << endl;
	
	TH1D *h_total_occupancy = new TH1D("h_total_occupancy_all_strips", "", occupancy_bin_num, 0, 1);
	TH1D *h_total_occupancy_U = new TH1D("h_total_occupancy_U", "", occupancy_bin_num, 0, 1);
	int total_event_num = 0;
	double total_occupancy = 0;
	double total_occupancy_U = 0;
	while(TC->GetEntry(total_event_num++)){
		TC->GetEntry(total_event_num);
		total_occupancy = nstrips_keep/max_strips;
		total_occupancy_U = nstrips_keepU/max_strips;
		h_total_occupancy->Fill(total_occupancy);
		h_total_occupancy_U->Fill(total_occupancy_U);
	}

	double last_occupancy_bin = (h_total_occupancy->FindLastBinAbove(0, 1) + 500)/double(occupancy_bin_num);
	double last_occupancy_bin_U = (h_total_occupancy_U->FindLastBinAbove(0, 1) + 500)/double(occupancy_bin_num);
	cout << "Finished calculating all total occupancies." << endl;

	TCanvas *c_total_occupancy = new TCanvas("c_total_occupancy", "", 700, 500);
	h_total_occupancy->Draw();
	h_total_occupancy->SetTitle(Form("Total Occupancy On GEM, All strips, All events - Run: %i - Cut: %s", runnum, cut));
	h_total_occupancy->GetXaxis()->SetTitle("Occupancy");
	h_total_occupancy->SetAxisRange(0, last_occupancy_bin, "X");
	h_total_occupancy->GetYaxis()->SetTitle("Entries");
	tpt_db_cuts->Draw();
	c_total_occupancy->Update();
	c_total_occupancy->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/%s/total_occupancy_%i_%s.pdf", cut, runnum, cut));


	TCanvas *c_total_occupancy_U = new TCanvas("c_total_occupancy_U", "", 700, 500);
	h_total_occupancy_U->Draw();
	h_total_occupancy_U->SetTitle(Form("Total Occupancy On GEM U Strips, ALL Events - Run: %i - Cut: %s", runnum, cut));
	h_total_occupancy_U->GetXaxis()->SetTitle("Occupancy");
	h_total_occupancy_U->SetAxisRange(0, last_occupancy_bin_U, "X");
	h_total_occupancy_U->GetYaxis()->SetTitle("Entries");
	c_total_occupancy_U->Update();
	c_total_occupancy_U->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/%s/total_occupancy_U_%i_%s.pdf", cut, runnum, cut));
	
	cout << "Analyzing data for APVs " << first_apv << " through " << last_apv << "." << endl;

	for(int apv_cnt = first_apv; apv_cnt < last_apv; apv_cnt++) {
		cout << "Processing data for apv: " << apv_cnt << endl;
		int min_APV_strip = APV_strip_nums(apv_cnt, "min");
		int max_APV_strip = APV_strip_nums(apv_cnt, "max");

		h_APV_U_ADCmax_chan_ratios[apv_cnt] = new TH1D(Form("h_APV%i_U_ADCmax_chan_ratios", apv_cnt), "", 300, 0, 30);
		h_APV_occupancy_U[apv_cnt] = new TH1D(Form("h_APV%i_occupancy_U", apv_cnt), "", occupancy_bin_num, 0.0, 1);
		h_apv_ADCmax_istrip[apv_cnt] = new TH2D(Form("h2_APV%i_U_ADCmax_istrip",apv_cnt), "",128, APV_strip_nums(apv_cnt, "min"), APV_strip_nums(apv_cnt, "max"), max_U_ADC, 0, max_U_ADC );
		h2_APV_U_ADCmax_apvstrip[apv_cnt] = new TH2D(Form("h2_APV%i_U_ADCmax_apvstrip", apv_cnt), "",128, 0, 128, max_U_ADC, 0, max_U_ADC );
		h2_APV_U_ADCmax_apvchan[apv_cnt] = new TH2D(Form("h2_APV%i_U_ADCmax_apvchan", apv_cnt), "", 128, 0, 128, max_U_ADC, 0, max_U_ADC );
		h2_APV_U_ADCmax_neigh_chan[apv_cnt] = new TH2D(Form("h2_APV%i_U_ADCmax_neigh_chan", apv_cnt), "", max_U_ADC, 0, max_U_ADC,max_U_ADC,  0, max_U_ADC);
		h2_APV_U_ADCmax_neigh_chan_lrg_sml[apv_cnt] = new TH2D(Form("h2_APV%i_U_ADCmax_neigh_chan_lrg_sml", apv_cnt), "", max_U_ADC, 0, max_U_ADC,max_U_ADC,  0, max_U_ADC);

		//We can either loop over all events or select some events:
		//Loop over all events:
		int apv_event_num = 0;
		while(TC->GetEntry(apv_event_num++)){
			TC->GetEntry(apv_event_num);

		
		// //Loop over selected events:
		// for(int evt = 0; evt <= 5000; evt++){
		// 	TC->GetEntry(evt);
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
				//Make occupancy plots for single APVs:
				// if(IsU[i] && (istrip[i] >= APV_strip_nums(apv_cnt, "min")) && (istrip[i] <= APV_strip_nums(apv_cnt, "max")) && ADC_max

				
			}
			// cout << "occupancy after event: " << evt << "   " << occupancy_cnt <<  "  & occupancy = " << double(occupancy_cnt)/nAPV_strips << endl;
			if(double(occupancy_cnt)/nAPV_strips > 0.0){h_APV_occupancy_U[apv_cnt]->Fill(double(occupancy_cnt)/nAPV_strips);}
			

			//The histograms for this event have been filled.
			//SCAN THROUGH FILLED HISTOGRAMS

			//First set some defaults
			double chan_n_adc = 0.0;
			double chan_n_plus_1_ADC = 0.0;
			double ADC_ratio = 0.0;

//*****************************/
//*********RATIOS**************/
			if(plot_ratio_all){
			//SCANNING THROUGH THE FILLED HISTOGRAMS
				for(int i = 0; i < 127; i++){
					double chan_n_adc = apv_chan_adc[i][1];
					double chan_n_plus_1_ADC = apv_chan_adc[i+1][1];
					double ADC_ratio = 0.0;

					//Always take larger ADC divided by smaller ADC for neighboring channels
					if(chan_n_adc > chan_n_plus_1_ADC && chan_n_plus_1_ADC != 0 && chan_n_adc > ADC_cut){
						double ADC_ratio = chan_n_adc/chan_n_plus_1_ADC;
						h_APV_U_ADCmax_chan_ratios[apv_cnt]->Fill(ADC_ratio);
						h2_APV_U_ADCmax_neigh_chan_lrg_sml[apv_cnt]->Fill(chan_n_adc, chan_n_plus_1_ADC);
					}
					else if(chan_n_plus_1_ADC > chan_n_adc && chan_n_adc != 0 && chan_n_plus_1_ADC > ADC_cut){
						double ADC_ratio = chan_n_plus_1_ADC/chan_n_adc;
						h_APV_U_ADCmax_chan_ratios[apv_cnt]->Fill(ADC_ratio);
						h2_APV_U_ADCmax_neigh_chan_lrg_sml[apv_cnt]->Fill(chan_n_plus_1_ADC, chan_n_adc);
					}
					//Fill channel n on the y-axis. Fill channel n+1 on the x-axis. 
					h2_APV_U_ADCmax_neigh_chan[apv_cnt]->Fill(chan_n_plus_1_ADC, chan_n_adc);
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

		apv_occupancy_mean[apv_cnt]=h_APV_occupancy_U[apv_cnt]->GetMean();

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

		// TCanvas *c_adcmax_ratios_apv[nAPVs];
		// c_adcmax_ratios_apv[apv_cnt] = new TCanvas(Form("APV%i Ratio of U-Strip channel ADCmax", apv_cnt), Form("c_ratio_apv_%i", apv_cnt), 700, 500);
		// h_APV_U_ADCmax_chan_ratios[apv_cnt]->Draw();
		// h_APV_U_ADCmax_chan_ratios[apv_cnt]->SetTitle(Form("Strip ratios (Run: %d, APV: %i, ADC cut: %i)", runnum, apv_cnt, ADC_cut));
		// h_APV_U_ADCmax_chan_ratios[apv_cnt]->GetXaxis()->SetTitle("ADC Ratio");
		// h_APV_U_ADCmax_chan_ratios[apv_cnt]->GetYaxis()->SetTitle("Entries");
		// h_APV_U_ADCmax_chan_ratios[apv_cnt]->SetMarkerStyle(2);
		// h_APV_U_ADCmax_chan_ratios[apv_cnt]->SetMarkerColor(06);
		// c_adcmax_ratios_apv[apv_cnt]->Update();
		
		// TCanvas *c_APV_U_ADCmax_neigh_chan[nAPVs];
		// c_APV_U_ADCmax_neigh_chan[apv_cnt] = new TCanvas(Form("APV%i, ADCmax on strip n vs ADCmax on strip n + 1", apv_cnt), Form("APV%i n and n+1", apv_cnt), 600, 500);
		// h2_APV_U_ADCmax_neigh_chan[apv_cnt]->GetZaxis()->SetRangeUser(10,2000);
		// h2_APV_U_ADCmax_neigh_chan[apv_cnt]->Draw("colz");
		// h2_APV_U_ADCmax_neigh_chan[apv_cnt]->SetTitle(Form("ADCmax on strip n vs ADCmax on strip n+1 (Run: %d)", runnum));
		// h2_APV_U_ADCmax_neigh_chan[apv_cnt]->GetXaxis()->SetTitle("ADCmax on strip (n+1)");
		// h2_APV_U_ADCmax_neigh_chan[apv_cnt]->GetYaxis()->SetTitle("ADCmax on strip n");
		// c_APV_U_ADCmax_neigh_chan[apv_cnt]->Update();

		// TCanvas *c_APV_U_ADCmax_neigh_chan_lrg_sml[nAPVs];
		// c_APV_U_ADCmax_neigh_chan_lrg_sml[apv_cnt] = new TCanvas(Form("APV%i, ADCmax: Larger of 2 strips on x, smaller on y", apv_cnt), "", 600, 500);
		// h2_APV_U_ADCmax_neigh_chan_lrg_sml[apv_cnt]->GetZaxis()->SetRangeUser(0,3000);
		// h2_APV_U_ADCmax_neigh_chan_lrg_sml[apv_cnt]->Draw("colz");
		// h2_APV_U_ADCmax_neigh_chan_lrg_sml[apv_cnt]->SetTitle(Form("ADCmax for n and (n+1): Larger ADC on x-axis; Smaller on y-axis. (Run: %d)", runnum));
		// h2_APV_U_ADCmax_neigh_chan_lrg_sml[apv_cnt]->GetXaxis()->SetTitle("ADCmax of larger ADC");
		// h2_APV_U_ADCmax_neigh_chan_lrg_sml[apv_cnt]->GetYaxis()->SetTitle("ADCmax of smaller ADC");
		// c_APV_U_ADCmax_neigh_chan_lrg_sml[apv_cnt]->Update();
		
		//SAVE APV OCCUPANCY PLOTS TO PDF
		double last_APV_occupancy_bin_U = (h_APV_occupancy_U[apv_cnt]->FindLastBinAbove(0, 1) + 500)/double(occupancy_bin_num);

		if(plot_APV_occupancy_U){
			TCanvas *c_APV_occupancy_U[nAPVs];
			c_APV_occupancy_U[apv_cnt] = new TCanvas(Form("c_APV%i_occupancy_U", apv_cnt), "", 700, 500);
			if(apv_cnt == first_apv){
				h_APV_occupancy_U[apv_cnt]->Draw();
				h_APV_occupancy_U[apv_cnt]->SetTitle(Form("Occupancy on APV %i - Run: %i - Cut: %s - Ustrips Only", apv_cnt, runnum, cut));
				h_APV_occupancy_U[apv_cnt]->GetXaxis()->SetTitle("Occupancy");
				h_APV_occupancy_U[apv_cnt]->GetYaxis()->SetTitle("Entries");
				h_APV_occupancy_U[apv_cnt]->SetAxisRange(0, last_APV_occupancy_bin_U, "X");
				c_APV_occupancy_U[apv_cnt]->Update();
				c_APV_occupancy_U[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/%s/APV_occupancy_U_%i_%s.pdf(", cut, runnum, cut));
				// c_APV_occupancy_U[apv_cnt]->Clear();
			}
			else if (apv_cnt > first_apv && apv_cnt < (last_apv-1)){
				h_APV_occupancy_U[apv_cnt]->Draw();
				h_APV_occupancy_U[apv_cnt]->SetTitle(Form("Occupancy on APV %i - Run: %i - Cut: %s - Ustrips Only", apv_cnt, runnum, cut));
				h_APV_occupancy_U[apv_cnt]->GetXaxis()->SetTitle("Occupancy");
				h_APV_occupancy_U[apv_cnt]->GetYaxis()->SetTitle("Entries");
				h_APV_occupancy_U[apv_cnt]->SetAxisRange(0, last_APV_occupancy_bin_U, "X");
				c_APV_occupancy_U[apv_cnt]->Update();
				c_APV_occupancy_U[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/%s/APV_occupancy_U_%i_%s.pdf", cut, runnum, cut));
				// c_APV_occupancy_U[apv_cnt]->Clear();
			}
			if(apv_cnt == (last_apv - 1)){
				h_APV_occupancy_U[apv_cnt]->Draw();
				h_APV_occupancy_U[apv_cnt]->SetTitle(Form("Occupancy on APV %i - Run: %i - Cut: %s - Ustrips Only", apv_cnt, runnum, cut));
				h_APV_occupancy_U[apv_cnt]->GetXaxis()->SetTitle("Occupancy");
				h_APV_occupancy_U[apv_cnt]->GetYaxis()->SetTitle("Entries");
				h_APV_occupancy_U[apv_cnt]->SetAxisRange(0, last_APV_occupancy_bin_U, "X");
				c_APV_occupancy_U[apv_cnt]->Update();
				c_APV_occupancy_U[apv_cnt]->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/%s/APV_occupancy_U_%i_%s.pdf)", cut, runnum, cut));
				
			}
			c_APV_occupancy_U[apv_cnt]->Close();
			gSystem->ProcessEvents();
		}

		//Put each APV occupancy histogram into the histogram for all APV h_all_APV_occupancies_U
		int apv_occupancy_bin = 0;
		for(int bin_cnt = first_occupancy_bin; bin_cnt < (first_occupancy_bin + occupancy_bin_num); bin_cnt ++){
			h_all_APV_occupancies_U->SetBinContent(bin_cnt, h_APV_occupancy_U[apv_cnt]->GetBinContent(apv_occupancy_bin));
			apv_occupancy_bin++;
		}
		
		first_occupancy_bin += occupancy_bin_num;
		
	// cout<< "mean: " << h_APV_occupancy_U[13]->GetMean();
	
	
		
	//END OF APV_CNT LOOP
	}
	
	TCanvas *c_all_APV_occupancies_U = new TCanvas("c_all_APV_occupancies_U", "", 700, 500);
	h_all_APV_occupancies_U->Draw("same");
	h_all_APV_occupancies_U->SetStats(0);
	h_all_APV_occupancies_U->SetTitle(Form("Occupancies on individual APVs - Run: %i - Cut: %s - Ustrips only", runnum, cut));
	h_all_APV_occupancies_U->GetXaxis()->SetTitle("APVs and respective occupany bins (0.0 - 0.30)");
	h_all_APV_occupancies_U->GetYaxis()->SetTitle("Entries");
	c_all_APV_occupancies_U->Update();
	cout << "h_apv_occ_u: " << (h_all_APV_occupancies_U->GetMaximum())*(2.0/3.0) << endl;
	for(int apv = first_apv; apv < last_apv; apv++){
		tpl_apv_occupancy_mean_U[apv] = new TPaveLabel(apv +.15, (h_all_APV_occupancies_U->GetMaximum())*(2.0/3.0), apv +.85, (h_all_APV_occupancies_U->GetMaximum()), Form("APV %i, Mean = %0.3f", apv, h_APV_occupancy_U[apv]->GetMean()));
		// tpl_apv_occupancy_mean_U[apv] = new TPaveLabel(apv +.15, 140, apv +.85, 230, Form("APV %i, Mean = %0.3f", apv, h_APV_occupancy_U[apv]->GetMean()));
		tpl_apv_occupancy_mean_U[apv]->SetBorderSize(1);
		tpl_apv_occupancy_mean_U[apv]->SetTextAngle(90.);
		tpl_apv_occupancy_mean_U[apv]->SetTextFont(53);
		tpl_apv_occupancy_mean_U[apv]->SetTextSize(12);
		tpl_apv_occupancy_mean_U[apv]->Draw();
		tl_apv_occupancy_lines[apv] = new TLine(apv, 0, apv, (h_all_APV_occupancies_U->GetMaximum())*(1.05));
		tl_apv_occupancy_lines[apv]->SetLineStyle(4);
		tl_apv_occupancy_lines[apv]->Draw();
		c_all_APV_occupancies_U->Update();
	}

	c_all_APV_occupancies_U->Print(Form("/work/halla/sbs/jboyd/analysis/xtalk/plots/%s/APV_all_occupancies_%i_%s.pdf", cut, runnum, cut));

///////////////////
//END OF MAIN LOOP
/////////////////////////////////////////////
}