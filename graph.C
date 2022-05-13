#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
using namespace std::chrono;
#include "TStyle.h"


#include "./include/include_files.h"
#include "./include/APV_strips.h"
#include "./include/search_file.C"

typedef vector <double> record_t;
typedef vector <record_t> data_t;

//-----------------------------------------------------------------------------
// Let's overload the stream input operator to read a list of CSV fields (which a CSV record).
// Remember, a record is a list of doubles separated by commas ','.
istream& operator >> ( istream& ins, record_t& record )
  {
  // make sure that the returned record contains only the stuff we read now
  record.clear();

  // read the entire line into a string (a CSV record is terminated by a newline)
  string line;
  getline( ins, line );

  // now we'll use a stringstream to separate the fields out of the line
  stringstream ss( line );
  string field;
  while (getline( ss, field, ';' ))
    {
    // for each field we wish to convert it to a double
    // (since we require that the CSV contains nothing but floating-point values)
    stringstream fs( field );
    double f = 0.0;  // (default value is 0.0)
    fs >> f;

    // add the newly-converted field to the end of the record
    record.push_back( f );
    }

  // Now we have read a single line, converted into a list of fields, converted the fields
  // from strings to doubles, and stored the results in the argument record, so
  // we just return the argument stream as required for this kind of input overload function.
  return ins;
  }

//-----------------------------------------------------------------------------
// Let's likewise overload the stream input operator to read a list of CSV records.
// This time it is a little easier, just because we only need to worry about reading
// records, and not fields.
istream& operator >> ( istream& ins, data_t& data )
  {
  // make sure that the returned data only contains the CSV data we read here
  data.clear();

  // For every record we can read from the file, append it to our resulting data
  record_t record;
  while (ins >> record)
    {
    data.push_back( record );
    }

  // Again, return the argument stream as required for this kind of input stream overload.
  return ins;  
  }

void graph(){

ifstream myfile;
myfile.open("crosstalk_v_occupancy_U.txt");

data_t data;
myfile >> data;

if (!myfile.eof())
{
	cout << "error in file!\n";

}

int first_row = 2;
Int_t num_runs = data.size()-first_row;
myfile.close();

cout << " the file containts " << data.size() << " records. \n";

unsigned max_record_size = 0;
cout << "data.size: " << data.size() << endl;
for (unsigned n = first_row; n < data.size(); n++){
	if (max_record_size < data[ n ].size()){
		max_record_size = data [ n ].size();
	}

}

int rows[num_runs];
Int_t runnum[num_runs];
double x[num_runs];
for(int i=0; i < num_runs; i++){x[i] = i;}
double occupancy[num_runs];
double occupancy_U[num_runs];
double crosstalk[num_runs];
double gaus_mean[num_runs];
double gaus_sigma[num_runs];
double Ndata_ADCmax[num_runs];
double Ndata_ADC_max_cut[num_runs];
double Ndata_IsU[num_runs];
double Ndata_IsU_cut[num_runs];

int ADC_cut = data[0][1];
int Noise_cut = data[0][3];

for(int row = first_row; row < data.size(); row++){
	rows[row - first_row] = row;
	runnum[row - first_row] = data[row][0];
	occupancy[row - first_row] = data[row][1];
	occupancy_U[row - first_row] = data[row][2];
	crosstalk[row - first_row] = data[row][3];
	gaus_mean[row - first_row] = data[row][4];
	gaus_sigma[row - first_row] = data[row][5];
	Ndata_ADCmax[row - first_row] = data[row][6];
	Ndata_ADC_max_cut[row - first_row] = data[row][7];
	Ndata_IsU[row - first_row] = data[row][8];
	Ndata_IsU_cut[row - first_row] = data[row][9];
	
}
TCanvas *c_gaus_mean = new TCanvas("c_gaus_mean", "", 700, 500);
TGraph *gr_gaus_mean = new TGraph(num_runs, occupancy, gaus_mean);
gr_gaus_mean->SetTitle(Form("Xtalk Gaussian Mean vs Total Occupancy (All strips), ADCcut: %i, Noise cut: %i", ADC_cut, Noise_cut));
gr_gaus_mean->GetXaxis()->SetTitle("Total Occupancy on all U AND V Strips");
gr_gaus_mean->GetYaxis()->SetTitle("Gaussian Mean of Crosstalk Fit (U-strips)");
gr_gaus_mean->Draw("AP*");

TCanvas *c_gaus_mean_U = new TCanvas("c_gaus_mean_U", "", 700, 500);
TGraph *gr_gaus_mean_U = new TGraph(num_runs, occupancy_U, gaus_mean);
gr_gaus_mean_U->SetTitle(Form("Xtalk Gaussian Mean vs U Occupancy (U strips only), ADCcut: %i, Noise cut: %i", ADC_cut, Noise_cut));
gr_gaus_mean_U->GetXaxis()->SetTitle("Occupancy on U-strips ONLY");
gr_gaus_mean_U->GetYaxis()->SetTitle("Gaussian Mean of Crosstalk Fit (U-strips)");
gr_gaus_mean_U->GetYaxis()->SetRangeUser(0, 10);
gr_gaus_mean_U->Draw("AP*");

TCanvas *c_gaus_sigma = new TCanvas("c_gaus_sigma", "", 700, 500);
TGraph *gr_gaus_sigma = new TGraph(num_runs, occupancy, gaus_sigma);
gr_gaus_sigma->SetTitle(Form("Xtalk Gaussian sigma vs Total Occupancy (All strips), ADCcut: %i, Noise cut: %i", ADC_cut, Noise_cut));
gr_gaus_sigma->GetXaxis()->SetTitle("Total Occupancy on all U AND V Strips");
gr_gaus_sigma->GetYaxis()->SetTitle("Gaussian sigma of Crosstalk Fit (U-strips)");
gr_gaus_sigma->GetYaxis()->SetRangeUser(0, 4);
gr_gaus_sigma->Draw("AP*");

TCanvas *c_gaus_sigma_U = new TCanvas("c_gaus_sigma_U", "", 700, 500);
TGraph *gr_gaus_sigma_U = new TGraph(num_runs, occupancy_U, gaus_sigma);
gr_gaus_sigma_U->SetTitle(Form("Xtalk Gaussian sigma vs U Occupancy (U strips only), ADCcut: %i, Noise cut: %i", ADC_cut, Noise_cut));
gr_gaus_sigma_U->GetXaxis()->SetTitle("Occupancy on U-strips ONLY");
gr_gaus_sigma_U->GetYaxis()->SetTitle("Gaussian sigma of Crosstalk Fit (U-strips)");
gr_gaus_sigma_U->GetYaxis()->SetRangeUser(0, 4);
gr_gaus_sigma_U->Draw("AP*");

TCanvas *c_xtalk_occ = new TCanvas("c_xtalk_occ", "", 700, 500);
TGraph *gr_xtalk_occ = new TGraph(num_runs, occupancy, crosstalk);
gr_xtalk_occ->SetTitle(Form("Crosstalk Integral vs Total Occupancy (all strips), ADCcut: %i, Noise cut: %i", ADC_cut, Noise_cut));
gr_xtalk_occ->GetXaxis()->SetTitle("Occupancy on all U AND V Strips");
gr_xtalk_occ->GetYaxis()->SetTitle("Integral of Crosstalk Gaussian");
gr_xtalk_occ->GetYaxis()->SetRangeUser(0, 0.1);
gr_xtalk_occ->Draw("AP*");

TCanvas *c_xtalk_occ_U = new TCanvas("c_xtalk_occ_U", "", 700, 500);
TGraph *gr_xtalk_occ_U = new TGraph(num_runs, occupancy_U, crosstalk);
gr_xtalk_occ_U->SetTitle(Form("Crosstalk Integral vs U Occupancy (U strips only), ADCcut: %i, Noise cut: %i", ADC_cut, Noise_cut));
gr_xtalk_occ_U->GetXaxis()->SetTitle("Occupancy on U strips only");
gr_xtalk_occ_U->GetYaxis()->SetTitle("Integral of Crosstalk Gaussian");
gr_xtalk_occ_U->GetYaxis()->SetRangeUser(0, 0.1);
gr_xtalk_occ_U->Draw("AP*");



// double occupancy_U[] = {0.237205, 0.117608, 0.217408, 0.20054, 0.212072, 0.243356, 0.200603, 0.277269, 0.243893, 0.245375, 0.29481, 0.297213, 0.295859, 0.287497, 0.269897, 0.212023, 0.256669, 0.2487788, 0.251899};

// double crosstalk[] = {0.0762284, 0.0209963, 0.0557538, 0.0373592, 0.0084385, 0.0145028, 0.035471, 0.0348434, 0.269705, 0.0176466, 0.0220618, 0.0318019, 0.0399211, 0.0505838, 0.0341498, 0.0573485, 0.0568664, 0.0900653, 0.0913458};



// TGraph *graph = new TGraph(19, occupancy_U, crosstalk);

// graph->Draw("AP*");
// graph->SetTitle("Crosstalk Occupancy vs Total Cccupancy on U-Strips, ADCcut = 500");
// graph->GetXaxis()->SetTitle("Total Occupancy on U-Strips");
// graph->GetYaxis()->SetTitle("Crosstalk Occupancy");


}