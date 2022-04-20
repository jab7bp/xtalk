#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <limits>
#include <stdio>
#include <cstdlib>
using std::cout;


string GetLine(string filename, int db_line){
	fstream myfile (filename);
	string y;
	for(int lineno = 0; getline (myfile, y) && lineno < db_line;lineno++){
		if(lineno == db_line-1){
			cout << y << endl;
		}
	}
	return y;
}
	
string search_file(string db_cut_variable){
	
	string db_cut_variable = "bb.gem.maxstrip_t0";
	string filename = "/work/halla/sbs/jboyd/SBS_OFFLINE/SBS_REPLAY/SBS-replay/DB_xtalk/db_bb.gem.dat";
	fstream db_file("/work/halla/sbs/jboyd/SBS_OFFLINE/SBS_REPLAY/SBS-replay/DB_xtalk/db_bb.gem.dat");
	string x;
	bool found_cut = false;
	int db_cut_line=1;

	if (db_file.is_open()){
		while( getline(db_file,x) ){
			if(x.find(db_cut_variable, 0) != string::npos){
				cout << "found bb.gem.maxstrip_t0 on line: " << db_cut_line << endl;
				found_cut = true;
				if(found_cut){break;}
			}
		db_cut_line++;
		}
	}
	string cut_variable;
	string var = GetLine(filename, 16);
	
	// GotoLine(db_file, db_cut_line) >> cut_variable;


	db_file.close();
	return var;

}