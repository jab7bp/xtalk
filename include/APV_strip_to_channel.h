#ifndef APV_strips
#define APV_strips

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