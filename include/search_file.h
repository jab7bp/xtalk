#ifndef search_file
#define search_file

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

int get_gemDB_lineno(string search_item){
	using namespace std;
	// string search_item = "bb.gem.maxstrip_t0";
	string filename = "/work/halla/sbs/jboyd/SBS_OFFLINE/SBS_REPLAY/SBS-replay/DB_xtalk/db_bb.gem.dat";
	fstream db_file(filename);

	string x;
	bool found_line = false;
	int search_item_line=1;

	if (db_file.is_open()){
		while( getline(db_file,x) ){
			if(x.find(search_item, 0) != string::npos){
				// cout << "found bb.gem.maxstrip_t0 on line: " << search_item_line << endl;
				found_line = true;
				if(found_line){break;}
			}
		search_item_line++;
		}
		
	}
	string cut_variable;
	
	TString var = GetLine(filename, search_item_line);
	
	// GotoLine(db_file, search_item_line) >> cut_variable;


	db_file.close();
	return search_item_line;
	

}

int parse_ped_line(TString db_ped_line){
	TObjArray *tokens = db_ped_line.Tokenize( "_" );
	TString runnum_parse = ( (TObjString*) (*tokens)[4] )->GetString();
	TObjArray *runnum_tokens = runnum_parse.Tokenize(".");
	TString runnum_part = ( (TObjString*) (*runnum_tokens)[0] )->GetString();
	TString num_part( runnum_part(3, 7) );
	int ped_runnum = atoi(num_part); 
	return ped_runnum;
}

// string get_pedfile(TString db_ped_line){
// 	TObjArray *tokens = db_ped_line.Tokenize( "_" );
// 	TString runnum_parse = ( (TObjString*) (*tokens)[4] )->GetString();
// 	TObjArray *runnum_tokens = runnum_parse.Tokenize(".");
// 	TString runnum_part = ( (TObjString*) (*runnum_tokens)[0] )->GetString();
// 	TString num_part( runnum_part(3, 7) );
// 	int ped_runnum = atoi(num_part); 
// 	return ped_runnum;
// }

std::fstream& GotoLine(std::fstream& file, int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

const char* get_gemped_file(int runnum){

	int counter = 1;
	int iter = 0;
	int ped_runs_line[20][2];

	string db_file = "/work/halla/sbs/jboyd/SBS_OFFLINE/SBS_REPLAY/SBS-replay/DB_xtalk/db_bb.gem.dat";

	string search = "bb.gem.pedfile";
	ifstream inFile;
	fstream file(db_file);
	string line;
	int row = 1;
	int next_row = 0;

	inFile.open(db_file);

	if(!inFile){
		cout << "Unable to open file" << endl;
		exit(1);
	}
	int first_cnt = get_gemDB_lineno(search);
	int latest_run_number = parse_ped_line(GetLine(db_file, first_cnt));
	ped_runs_line[0][0] = first_cnt;
	ped_runs_line[0][1] = latest_run_number;

	while(latest_run_number < runnum){
		GotoLine(file, ped_runs_line[counter-1][0] + 1);
		size_t pos;
		while(file.good())
		  {
		      getline(file,line); // get line from file
		      pos=line.find(search); // search

		      if(pos!=string::npos) // string::npos is returned if string is not found
		        {
		            next_row = row + first_cnt + iter;
		            ped_runs_line[counter][0] = next_row;
		            latest_run_number = parse_ped_line(GetLine(db_file, next_row));
		            ped_runs_line[counter][1] = latest_run_number;
		            break;
		        }
		        row++;
		  }
		counter++;
		iter++;

	}

	cout << "going to use this pedfile: " << ped_runs_line[counter - 2][1] << " on line: " << ped_runs_line[counter - 2][0] << endl;


	TObjArray *tokens2 = GetLine(db_file, ped_runs_line[counter - 2][0]).Tokenize( "=" );
	TString pedfile_parse = ( (TObjString*) (*tokens2)[1] )->GetString();
	TString pedfile_part_TS( pedfile_parse(1, 40) );
	std::string pedfile_part(pedfile_part_TS.Data());

	return pedfile_part.c_str();

}

#endif