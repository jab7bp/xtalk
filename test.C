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

void test(){
	for(int i = 0; i< 128; i++){
		cout << "strip: " << i << "   Channel: " << UV_APV_strip_to_channel(i) << endl;

	}

	cout << "tesint new function: channel 0 = strip: " << UV_APV_channel_to_strip(0) << endl;
}