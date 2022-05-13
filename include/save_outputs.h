#ifndef save_outputs
#define save_outputs

void print_single_PDF(TCanvas *canvas, const char* filename, const char* dir, int start, int current, int end){
	const char* output = Form("%s%s.pdf", dir, filename);
	if(start == current){
		canvas->Print(Form("%s(", output));
	}
	else if( (current > start) &&  (current < (end -1) ) ){
		canvas->Print(output);
	}
	if(current == (end - 1)){
		canvas->Print(Form("%s)", output));
	}
	canvas->Close();
	gSystem->ProcessEvents();

}

void print_multi_PDF(TCanvas *canvas, const char* filename, const char* dir, bool first, bool last){
	const char* output = Form("%s%s.pdf", dir, filename);

	if(first){
		canvas->Print(Form("%s(", output));
	}
	else if( !(first) && !(last) ){
		canvas->Print(output);
	}
	if( last ){
		canvas->Print(Form("%s)", output));
	}
	canvas->Close();
	gSystem->ProcessEvents();

}

#endif