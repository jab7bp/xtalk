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
	cout << "tesint new function: channel 11 = strip: " << UV_APV_channel_to_strip(11) << endl;
	cout << "tesint new function: channel 12 = strip: " << UV_APV_channel_to_strip(12) << endl;
	cout << "tesint new function: channel 15 = strip: " << UV_APV_channel_to_strip(15) << endl;
	cout << "tesint new function: channel 27 = strip: " << UV_APV_channel_to_strip(27) << endl;
	cout << "tesint new function: channel 28 = strip: " << UV_APV_channel_to_strip(28) << endl;
	cout << "tesint new function: channel 29 = strip: " << UV_APV_channel_to_strip(29) << endl;
	cout << "tesint new function: channel 39 = strip: " << UV_APV_channel_to_strip(39) << endl;
	cout << "tesint new function: channel 43 = strip: " << UV_APV_channel_to_strip(43) << endl;
	cout << "tesint new function: channel 55 = strip: " << UV_APV_channel_to_strip(55) << endl;
	cout << "tesint new function: channel 59 = strip: " << UV_APV_channel_to_strip(59) << endl;

}