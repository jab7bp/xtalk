#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
#include <limits>
using namespace std::chrono;

#include "/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/include/APV_strips.h"
#include "/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/include/search_file.C"
#include "/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/include/draw_plots.h"
#include "/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/include/save_outputs.h"
#include "/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/include/calc_functions.h"


TFile *outputRootFile = new TFile("xtalk_correlation.root", "RECREATE");

//Define some global variables


TH1D *h_channels_ADCs[128];
TH2D *h2_channels_ADCs = new TH2D("h2_channels_ADCs", "", 128, 0, 127, 128, 0, 127);
TH1D *h_chan1_corr_coeff = new TH1D("h_chan1_corr_coeff", "", 128, 0, 127);
TH1D *h_chan2_corr_coeff = new TH1D("h_chan2_corr_coeff", "", 128, 0, 127);
TH1D *h_chan3_corr_coeff = new TH1D("h_chan3_corr_coeff", "", 128, 0, 127);
TH2D *h2_correlation_coeff = new TH2D("h2_correlation_coeff", "", 128, 0, 127, 128, 0, 127);
TH2D *h2arr_correlation_coeff = new TH2D("h2arr_correlation_coeff", "", 128, 0, 127, 128, 0, 127);
TH2D *h2_correlation_coeff_fill = new TH2D("h2_correlation_coeff_fill", "", 128, 0, 127, 128, 0, 127);
TH2D *h2_neigh_chan_7_8_ADCs = new TH2D("h2_neigh_chan_7_8_ADCs", "", 128, 0, 127, 128, 0, 127);
TH2D *h2arr_neigh_chan_7_8_ADCs = new TH2D("h2arr_neigh_chan_7_8_ADCs", "", 128, 0, 127, 128, 0, 127);


int apv_chan_adc[128][2];
int apv_chan_adc_crosstalk_reject[128][2];
int apv_strip_adc[128][2];

double apv_adc_cut_thresh = 50.0;

int nAPVs = 30;
int first_apv = 13;
int last_apv = 14;

int first_event = 0;
int last_event = 80000;
double apv_adc_channels[128][80000];
double ADC_max_thresh = 4000;

double channel_sigma[128];
double channel_mu[128];


double correlation_coeff_1D[128];
double correlation_coeff_2D[128][128];

int nstrips_keep_cut;


//Set defaults for ADC on channel n and n+1
double chan_n_adc = 0.0;
double chan_n_plus_1_ADC = 0.0;


int nstrips = 3840;
Int_t max_strips = 8000;
Int_t nAPV_strips = 128;


Int_t Ndata_strip_ADC_samples;
Int_t Ndata_strip_ADC_max;
int Ndata_strip_ADC_max_keep;
int Ndata_strip_ADC_max_keepU;
int count_in_isu;


//BOOLEANS
bool build_all = false;
bool build_occupancies = true;
bool build_ratios = false;

void xtalk_correlation(int const ADC_cut = 500){



	ofstream myfile;
	myfile.open("crosstalk_v_occupancy_U.txt");
	myfile << Form("ADCcut; %i", ADC_cut) << Form("; Noise_cut; %i", int(apv_adc_cut_thresh)) << endl;
	myfile << "run	occupancy	occupancy_U		crosstalk	gaus_mean	gaus_sigma	Ndata_ADCmax	Ndata_ADCmax_cut	Ndata_IsU	Ndata_IsU_cut		\n \n";

	auto start = high_resolution_clock::now();
	int runs[] = {11562};
	// int runs[] = {11449, 11451, 11456, 11494, 11562, 11580, 11595, 11997, 12001, 12013, 12030, 12050, 12060, 12073, 12342, 12423, 12424, 12425};
	// int runs[] = {11449, 11451, 11456, 11494, 11562, 11580, 11595, 11997, 12001, 12013, 12030, 12050, 12060, 12073, 12342, 12423, 12424, 12425, 12550, 12620, 12662, 12728, 13060, 13309, 13325, 13344, 13370, 13400, 13454, 13474, 13505, 13554, 13560, 13615, 13620, 13660, 13661, 13664, 13666, 13680, 13685, 13732, 13770, 13799};
	int num_runs = (sizeof(runs)/sizeof(runs[0]));

	double arr_Ndata_ADCmax[num_runs];
	double arr_Ndata_ADCmax_cut[num_runs];
	double arr_Ndata_IsU[num_runs];
	double arr_Ndata_IsU_cut[num_runs];

	
	//Beginning of RUNS Loop
	for(int irun = 0; irun < num_runs; irun++){
		int runnum = runs[irun];
		
		TChain *TC = new TChain("T");

		const char * DATA_DIR =  "/lustre/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk";
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
		for( int i = 0; i < 128; i++){
			h_channels_ADCs[i] = new TH1D(Form("h_channel%i_ADCs", i), "", ADC_max_thresh, 0, ADC_max_thresh-1);
		}
		

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

		
	//EVENT //Loop all selected events:		
			// We can either loop over all events or select some events:
			// Loop over all events:

			// int apv_event_num = 0;
			// while(TC->GetEntry(apv_event_num++)){
			// 	if(apv_event_num%10000 == 0){cout << "Analyzing event " << apv_event_num << " of " << TC->GetEntries() << " total events." << endl;}
			// 	TC->GetEntry(apv_event_num);

	// EVENT //Loop over selected events:
			cout << "**************************************************" << endl;
			cout << "Number of entries: " << TC->GetEntries() << endl;
			cout << "**************************************************" << endl;
			for(int evt = 0; evt <= last_event; evt++){
				

				// cout << "evt: " << evt << endl;
				// if(evt == 706 || evt == 1186) {continue;}
				if( evt%10000 == 0) {cout << "evt: " << evt << endl;}
				TC->GetEntry(evt);

		
	//----------Calculate Occupancies---------------
				int total_event_num = 0;
		
				

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
						// if(ADC_max[i] > apv_adc_cut_thresh){
							h2_ADCmax_U->Fill(istrip[i], ADC_max[i]);

							//Only select strips that are in the selected range for the specific APV
							//Need to cross reference the global strip numbers and relate them to those of the APV_strip_nums


							if( (istrip[i] >= min_APV_strip) && (istrip[i] <= max_APV_strip) ){
								//Fill an array with the maximum ADC on each channel.
								//Storing it in a 2D array to grab the maximum entry in each bin

								//-----CHANNELS-------
								int chan = UV_APV_strip_to_channel(int(istrip[i])%128);
								// if(chan == 0 ){ cout << "Chan is 0 on evt: " << evt << " with adc: " << ADC_max[i] << endl;}
								// cout << "istrip[i]: " << istrip[i] << "     strip: " << int(istrip[i])%128 << "     Channel on APV: " << chan << endl;

								// apv_chan_adc[UV_APV_strip_to_channel(int(istrip[i])%128)][0] = chan;
								// apv_chan_adc[UV_APV_strip_to_channel(int(istrip[i])%128)][1] = ADC_max[i];
								if(ADC_max[i] < ADC_max_thresh && ADC_max[i] > 0.0){ 
									cout << "evt: " << evt << " getting: " << ADC_max[i] << endl;
									apv_adc_channels[chan][evt] = ADC_max[i];
								}
									
								//-----STRIPS---------
								apv_strip_adc[int(istrip[i])%128][0] = istrip[i];
								apv_strip_adc[int(istrip[i])%128][1] = ADC_max[i];
								

							}

							Ndata_strip_ADC_max_keepU++;
						// }

						
					}


				}
				
				// cout << "ndata_strip_u after cuts: " << Ndata_strip_ADC_max_keepU << endl;
				// cout << "count_in_isu: " << count_in_isu << endl;
				// cout << "*****************************" << endl;



				//End of looping through U-Strips
				

			}
			cout << "Filling Channel ADC histograms" << endl;
			for(int i = 0; i < last_event; i++){
				if(apv_adc_channels[8][i] > 0){
					h_channels_ADCs[8]->Fill(apv_adc_channels[8][i]);
				}
			}

			// for(int i = 0; i < 128; i++){
			// 	for(int j = 0; j < last_event; j++){
			// 		h2_channels_ADCs->Fill(apv_adc_channels[1][j], apv_adc_channels[2][j]);
			// 	}
			// }

			//END of EVENTS
			cout << "length of apv_adc_channels: " << (sizeof(apv_adc_channels[chan])/sizeof(apv_adc_channels[chan][0])) << endl;
			

			cout << "**************************************************" << endl;
			cout << "Finished going through all events." << endl;
			cout << "**************************************************" << endl;
		}
		//END OF APVs
		
		for(int i = 0; i < 128; i++){
			cout << "Chan: " << i << "   ADCs: ";
			for(int j = 20; j < 40; j++){
				cout << apv_adc_channels[i][j] << " ";
			}
			cout << endl;
		}

		TCanvas *c_chan_ADCs[128];
		
		arr_Ndata_ADCmax[irun] = Ndata_strip_ADC_max;
		arr_Ndata_ADCmax_cut[irun] = Ndata_strip_ADC_max_keep;
		arr_Ndata_IsU[irun] = count_in_isu;
		arr_Ndata_IsU_cut[irun] = Ndata_strip_ADC_max_keepU;

		// for(int i = 0; i < 128; i++){
		// 	c_chan_ADCs[i] = new TCanvas(Form("c_chan_%i_ADCs",i), "", 700, 500);
		// 	h_channels_ADCs[i]->SetTitle(Form("APV13 Channel_%i ADCmax Values", i));
		// 	h_channels_ADCs[i]->GetXaxis()->SetTitle(Form("Channel %i ADCmax", i));
		// 	h_channels_ADCs[i]->GetYaxis()->SetTitle("Entries");
		// 	h_channels_ADCs[i]->Draw();

		// 	print_single_PDF(c_chan_ADCs[i], Form("APV13_channels_ADCs_run_%i", runs[irun]), "/work/halla/sbs/jboyd/analysis/xtalk/plots/correlation/", 0, i, 128);
		// }
		
		
		
//CALCULATE THE CORRELATION COEFFICIENT


		//All VALUES Mean
		for(int chan = 0; chan < 128; chan++){
			int array_n = (sizeof(apv_adc_channels[chan])/sizeof(apv_adc_channels[chan][0]));
			double array_sum = 0.0, mean = 0.0, StdDev_sum = 0.0, StdDev = 0.0;

			for(int i = 0; i < array_n; i++){
				array_sum += apv_adc_channels[chan][i];
			}
			mean = array_sum/array_n;

			channel_mu[chan] = mean;
			
			for(int j = 0; j < array_n; j++){
				StdDev_sum += pow(apv_adc_channels[chan][j] - mean, 2);
			}

			StdDev = sqrt(StdDev_sum / array_n);

			channel_sigma[chan] = StdDev;
			cout << "chan: " << chan <<"  Mean: " << channel_mu[chan] << "   StdDev: " << channel_sigma[chan] << endl;
		}

		// for(int i = 0; i < 128; i++){
		// 	channel_mu[i] = h_channels_ADCs[i]->GetMean();
		// 	channel_sigma[i] = h_channels_ADCs[i]->GetStdDev();
		// 	cout << "chan: " << i <<"  Mu_chan: " << channel_mu[i] << "   sigma_chan: " << channel_sigma[i] << endl;
		// }
		
	

		// Calculate the coefficient using two channel arrays. Array dimensions is: 128 [chan] x N [events]
		// Sum through all events on each channel.
		// Equation: sum += (channel_i_ADC_for_event_k - channel_i_mu) * (channel_j_ADC_for_event_k - channel_j_mu)
		// correlation coefficient: (sum / N events) / (channel_i_sigma * channel_j_sigma)
	    cout << "**************************************************" << endl << endl;	
		cout << "Calculating the correlation coefficient: " << endl;

		for(int i = 0; i < 128; i++){
			for(int j = 0; j < 128; j++){
				double sum = 0.0;
				for(int k = 0; k < last_event; k++){
					sum += (apv_adc_channels[i][k] - channel_mu[i]) * (apv_adc_channels[j][k] - channel_mu[j]);
				}
				correlation_coeff_2D[i][j] = (sum) / ( (double)last_event * (channel_sigma[i] * channel_sigma[j]) ) ;
			}
		}

		cout << "Filling correlation coefficient histogram. " << endl;
		for(int i = 0; i < 128; i++){
			for(int j = 0; j < 128; j++){
				h2_correlation_coeff->SetBinContent(i, j, (correlation_coeff_2D[i][j]));
			}
		}

		TCanvas *c2_correlation_coeff = new TCanvas("c2_correlation_coeff", "", 700, 500);
		h2_correlation_coeff->Draw("colz");
		// TCanvas *c2_correlation_coeff = new TCanvas("c2_correlation_coeff", "", 700, 500);
		// h2_correlation_coeff->Draw("colz");

		// myfile << runs[irun] << "; " << arr_xtalk_v_occupancies[irun][0] << "; " << arr_xtalk_v_occupancies_U[irun][0] << "; " << arr_xtalk_v_occupancies_U[irun][1] << "; " << xtalk_mean << "; " << xtalk_sigma << "; " << Ndata_strip_ADC_max << "; " << Ndata_strip_ADC_max_keep << "; " << count_in_isu << "; " << Ndata_strip_ADC_max_keepU << endl;
		// cout << "**************************************************" << endl;
		// for(int i = 0; i < 128; i++){
		// 	cout << "channel: " << i << ", strip: " << UV_APV_channel_to_strip(i) << "  ---  ";
		// 	for(int j = first_event; j < 10; j++){
		// 		cout << correlation_coeff_2D[i][j] << "   ";
		// 	}
		// 	cout << endl;
		// }


		cout << "Filling neighboring channels ADC histograms. " << endl;


		// h2_neigh_chan_7_8_ADCs->Draw();

		// for (int i = 0; i < 128; i++)       // loops through each row of vy
		//    {   cout << "chan " << i << ": ";
		//    for (int j = 0; j < 10; j++) // loops through each element of each row 
		//           cout << " " << apv_adc_channels[i][j];           // prints the jth element of the ith row
		//       cout << endl;
  //  		}

	}
	//END of RUNS Loop
	
//print out the crosstalk array stuff:



// plot_1DH(h_xtalk_v_occupancy, "xtalk_v_occupancy", Form("Crosstalk Occupancy vs Total Cccupancy on U-Strips, ADCcut = %i, Noise Cut = %i", ADC_cut, int(apv_adc_cut_thresh)), "Total Occupancy on U-Strips", "Crosstalk Occupancy", {0.0, 0.0}, {0.0, 0.0});

auto stop = high_resolution_clock::now();
auto duration = duration_cast<minutes>(stop - start);	
cout << "Time elapsed: " << duration.count() << " minutes." << endl;
}


	