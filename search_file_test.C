#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
using std::cerr;
using std::cout;
using std::endl;
#include <cstdlib>
#include <limits>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"


#include "./include/search_file.h"



void search_file_test(int runnum = 10400){
		
	std::vector<std::pair<int,char>> gemped_info;
	gemped_info.push_back( std::pair<int,string>( 11436, "Thu Oct 14 01:58:43 2021"));

	// TString gemped_file = Form("/work/halla/sbs/jboyd/SBS_OFFLINE/SBS_REPLAY/SBS-replay/DB_xtalk/%s", get_gemped_file(runnum));
	// cout << "gemped file: " << gemped_file << endl;

	// string search = "bb.gem.m0.rmsu";

	// ifstream inFile;
	// fstream file(gemped_file);
	// string line;
	// int line_num = 0;
	// int row = 1;

	// inFile.open(gemped_file);

	// if(!inFile){
	// 	cout << "Unable to open file" << endl;
	// 	exit(1);
	// }

	// line_num = get_gemDB_lineno(search);

	// cout << "line num: " << line_num << endl;
}

// int ped_runnum = atoi(num_part); 



// // TObjArray *tokens = search_gemDB("bb.gem.pedfile").Tokenize( "_" );

// // TString skey = ( (TObjString*) (*tokens)[4] )->GetString();

	
// 		int search_line_num = get_gemDB_lineno("bb.gem.pedfile");

// 		cout << "search Line num: " << search_line_num << endl;

// 		// TObjArray *tokens = search_gemDB("bb.gem.pedfile", 1229).Tokenize( "_" );

// 		// TString runnum_parse = ( (TObjString*) (*tokens)[0] )->GetString();

// 		// cout << "runnum parse: " << runnum_parse << endl;

// 		// TObjArray *runnum_tokens = runnum_parse.Tokenize(".");

// 		// TString runnum_part = ( (TObjString*) (*runnum_tokens)[0] )->GetString();

// 		// TString num_part( runnum_part(3, 7) );

// 		// int ped_runnum = atoi(num_part); 

// 		// cout << "num_part: " << num_part << "   ped_runnum: " << ped_runnum << endl;
	

