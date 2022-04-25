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

void apv_multiplexer(TString conversion, int num){

	if(conversion == "strip2chan"){
		cout << "Physical strip: " << num << endl << "Correspoding channel: " << UV_APV_strip_to_channel(num) << endl;
	}

	if(conversion == "chan2strip"){
		cout << "Multiplexer channel: " << num << endl << "Corresponding physical strip: " << UV_APV_channel_to_strip(num) << endl;
	}

}