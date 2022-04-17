#include<iostream>
#include<TH1D.h>
#include<TF1.h>
#include<TCanvas.h>
#include<TRandom.h>
#include<TStyle.h>
#include<TLegend.h>
#include<TROOT.h>
#include<TPaveText.h>

void lansim() {
	gROOT->Reset();
	TStyle * plain = new TStyle("plain","plain");
	plain->SetCanvasBorderMode(0);
	plain->SetPadBorderMode(0);
	plain->SetPadColor(0);
	plain->SetCanvasColor(0);
	plain->SetTitleColor(1);
	plain->SetStatColor(0);
	plain->SetTitleFillColor(0);
	gROOT->SetStyle("plain");
	gStyle->SetPalette(1);
	

	int max_entries = 20000;
	int n_bins = 500;
	int lan_cen = 300;
	int lan_wid = 50;

	//Create empty histograms
	TH1F * hlandau1 = new TH1F("hlandau1", "", n_bins, 0, 2000);
	TH1F * hlandau2 = new TH1F("hlandau2", "", n_bins, 0, 2000);
	TH1F * ratio = new TH1F("ratio", "", n_bins, 0, 30);

	//Fill histograms with random Landau distributions (a different one in each)
	for(double_t i = 0; i < max_entries; i++) {
		float ADC1 = gRandom->Landau(lan_cen, lan_wid);
		float ADC2 = gRandom->Landau(lan_cen, lan_wid);
		hlandau1->Fill(ADC1);
		hlandau2->Fill(ADC2);

		float r;

		if(ADC1 > ADC2) {
			r = ADC1/ADC2;
		}
		else {
			r = ADC2/ADC1;
		}
		ratio->Fill(r);
	}
		
	//NO CROSS TALK
	//Calculate the ratio:
	double_t BinRatio;

	

	TCanvas * c1 = new TCanvas("c1", "landau1", 200, 10, 600, 600);
	hlandau1->Draw();

	TCanvas * c2 = new TCanvas("c2", "landau2", 200, 10, 600, 600);
	hlandau2->Draw();
	


	TCanvas * cratio = new TCanvas("cratio", "ratio", 200, 10, 600, 600);
	ratio->Draw();
	
}