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
#include "./include/calc_functions.h"

void test(){

	int apv_cnt = 15;

	double array[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

	int array_n = (sizeof(array)/sizeof(array[0]));
	double array_sum = 0.0, mean, StdDev_sum = 0.0, StdDev;

	for(int i = 0; i < array_n; i++){
		array_sum += array[i];
	}
	mean = array_sum/array_n;
	
	for(int i = 0; i < array_n; i++){
		StdDev_sum += pow(array[i] - mean, 2);
	}

	StdDev = sqrt(StdDev_sum / array_n);
	cout << "Mean: " << mean << "   StdDev: " << StdDev << endl;

	// for(int apv_cnt = 0; apv_cnt< 30; apv_cnt++){
	// 	int min_APV_strip = APV_strip_nums(apv_cnt, "min");
	// 	int max_APV_strip = APV_strip_nums(apv_cnt, "max");
	// 	cout << "apv: " << apv_cnt << "   min: " << min_APV_strip << "  max: " << max_APV_strip << endl;

	// }
	// cout << "tesint new function: channel 11 = strip: " << UV_APV_channel_to_strip(11) << endl;
	// cout << "tesint new function: channel 12 = strip: " << UV_APV_channel_to_strip(12) << endl;
	// cout << "tesint new function: channel 15 = strip: " << UV_APV_channel_to_strip(15) << endl;
	// cout << "tesint new function: channel 27 = strip: " << UV_APV_channel_to_strip(27) << endl;
	// cout << "tesint new function: channel 28 = strip: " << UV_APV_channel_to_strip(28) << endl;
	// cout << "tesint new function: channel 29 = strip: " << UV_APV_channel_to_strip(29) << endl;
	// cout << "tesint new function: channel 39 = strip: " << UV_APV_channel_to_strip(39) << endl;
	// cout << "tesint new function: channel 43 = strip: " << UV_APV_channel_to_strip(43) << endl;
	// cout << "tesint new function: channel 55 = strip: " << UV_APV_channel_to_strip(55) << endl;
	// cout << "tesint new function: channel 59 = strip: " << UV_APV_channel_to_strip(59) << endl;

}