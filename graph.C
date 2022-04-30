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

void graph(){

// ifstream myfile;
// myfile.open("crosstalk_v_occupancy_U.txt");

// std::string line;

// while( std::getline(myfile, line)){

// 	std::stringstream ss(line);

// 	std::string run, occupancy, occupancy_U, crosstalk, gaus_mean, gaus_sigma;
// 	std::getline(ss, run, ';'); std::cout<<"\""<<run<<"\"";
// 	std::getline(ss, occupancy, ';'); std::cout<<"\""<<occupancy<<"\"";
// }


double occupancy_U[] = {0.237205, 0.117608, 0.217408, 0.20054, 0.212072, 0.243356, 0.200603, 0.277269, 0.243893, 0.245375, 0.29481, 0.297213, 0.295859, 0.287497, 0.269897, 0.212023, 0.256669, 0.2487788, 0.251899};

double crosstalk[] = {0.0762284, 0.0209963, 0.0557538, 0.0373592, 0.0084385, 0.0145028, 0.035471, 0.0348434, 0.269705, 0.0176466, 0.0220618, 0.0318019, 0.0399211, 0.0505838, 0.0341498, 0.0573485, 0.0568664, 0.0900653, 0.0913458};

// TH1D *histo = new TH1D("histo", "", 1000, 0, 1);
// for(int i = 0; i < 19; i++){
// 	histo->SetBinContent(1000*occupancy_U[i], 100*crosstalk[i]);
// }

TGraph *graph = new TGraph(19, occupancy_U, crosstalk);

graph->Draw("AP*");
graph->SetTitle("Crosstalk Occupancy vs Total Cccupancy on U-Strips, ADCcut = 500");
graph->GetXaxis()->SetTitle("Total Occupancy on U-Strips");
graph->GetYaxis()->SetTitle("Crosstalk Occupancy");


}