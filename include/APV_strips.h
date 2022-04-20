#ifndef APV_strips
#define APV_strips

double APV_strip_nums(int APV, TString minmax){
	int apv_min_max_strips[30][2] =	{{0, 127},
									{128, 255},
									{256, 383},
									{384, 511},
									{512, 639},
									{640, 767},
									{768, 895},
									{896, 1023},
									{1024, 1151},
									{1152, 1279},
									{1280, 1407},
									{1408, 1535},
									{1536, 1663},
									{1664, 1791},
									{1792, 1919},
									{1920, 2047},
									{2048, 2175},
									{2176, 2303},
									{2304, 2431},
									{2432, 2559},
									{2560, 2687},
									{2688, 2815},
									{2816, 2943},
									{2944, 3071},
									{3072, 3199},
									{3200, 3327},
									{3328, 3455},
									{3456, 3583},
									{3584, 3711},
									{3712, 3839}};
	
	int apv_min = apv_min_max_strips[APV][0];
	int apv_max = apv_min_max_strips[APV][1];

	if(minmax == "min"){ return apv_min; }
	else if (minmax == "max"){ return apv_max; }
	
	
}

int UV_APV_strip_to_channel(int strip) {
	
	//Mapping for APV channel to strip
	const int _mapped_strip_uva_uv[128] = {
	     31,  15, 127, 111,  27,  11, 123, 107,  23,   7,
	    119, 103,  19,   3, 115,  99,  30,  14, 126, 110,
	     26,  10, 122, 106,  22,   6, 118, 102,  18,   2,
	    114,  98,  29,  13, 125, 109,  25,   9, 121, 105,
	     21,   5, 117, 101,  17,   1, 113,  97,  28,  12,
	    124, 108,  24,   8, 120, 104,  20,   4, 116, 100,
	     16,   0, 112,  96,  32,  48,  64,  80,  36,  52,
	     68,  84,  40,  56,  72,  88,  44,  60,  76,  92,
	     33,  49,  65,  81,  37,  53,  69,  85,  41,  57,
	     73,  89,  45,  61,  77,  93,  34,  50,  66,  82,
	     38,  54,  70,  86,  42,  58,  74,  90,  46,  62,
	     78,  94,  35,  51,  67,  83,  39,  55,  71,  87,
	     43,  59,  75,  91,  47,  63,  79,  95
	};

	if(strip >= 128){
		cout << endl << " WARNING!!! APVs only have strips from 0 - 127. Cannot have strips >= 128." << endl;
		int APV_channel = 99999;
		return APV_channel;
	}
	else{
		int n = sizeof(_mapped_strip_uva_uv)/sizeof(_mapped_strip_uva_uv[0]);
		auto chan = find(_mapped_strip_uva_uv, _mapped_strip_uva_uv + n, strip);
		int APV_channel = distance(_mapped_strip_uva_uv, chan);
		return APV_channel;
	}
	
}

#endif