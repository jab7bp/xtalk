#include <iostream>
#include <string>
#include <fstream>
using std::cerr;
using std::cout;
using std::endl;
#include <cstdlib>
#include <limits>
#include "TString.h"

// std::fstream& GotoLine(std::fstream& db_dat, int line){
// 	db_dat.seekg(std::ios::beg);
// 	for(int i=0; i < line-1; i++){
// 		db_dat.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
// 	}
// 	return db_dat;
// }

TString GetLine(string filename, int db_line){
	fstream myfile (filename);
	string y;
	for(int lineno = 0; getline (myfile, y) && lineno < db_line;lineno++){
		if(lineno == db_line-1){
			break;
		}
	}
	return y;

}
	
TString search_file(string db_cut_variable){
	using namespace std;
	// string db_cut_variable = "bb.gem.maxstrip_t0";
	string filename = "/work/halla/sbs/jboyd/SBS_OFFLINE/SBS_REPLAY/SBS-replay/DB_xtalk/db_bb.gem.dat";
	fstream db_file("/work/halla/sbs/jboyd/SBS_OFFLINE/SBS_REPLAY/SBS-replay/DB_xtalk/db_bb.gem.dat");
	string x;
	bool found_cut = false;
	int db_cut_line=1;

	if (db_file.is_open()){
		while( getline(db_file,x) ){
			if(x.find(db_cut_variable, 0) != string::npos){
				// cout << "found bb.gem.maxstrip_t0 on line: " << db_cut_line << endl;
				found_cut = true;
				if(found_cut){break;}
			}
		db_cut_line++;
		}
		
	}
	string cut_variable;
	
	TString var = GetLine(filename, db_cut_line);
	
	// GotoLine(db_file, db_cut_line) >> cut_variable;


	db_file.close();
	return var;
	

}