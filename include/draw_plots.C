#ifndef draw_plots
#define draw_plots

void plot_1DH(TH1D * histogram, char * c_title = "", TString title  = "", TString xTitle = "", TString yTitle = "", double xRange[]={0.0, 0.0}, double yRange[]={0.0, 0.0}, TString options = "" )
{
	TCanvas * canvas = new TCanvas(Form("c_%s", c_title), "", 700, 500);
	histogram->SetTitle(title);
	histogram->GetXaxis()->SetTitle(xTitle);
	histogram->GetYaxis()->SetTitle(yTitle);
	if(TMath::MaxElement(2, xRange) != 0){ histogram->GetXaxis()->SetRange(xRange[0], xRange[1]); }
	if(TMath::MaxElement(2, yRange) != 0){ histogram->GetYaxis()->SetRangeUser(yRange[0], yRange[1]); }
	histogram->Draw(options);
}


// void plot_2DH(TH2D * histo_name, TString title  = "", TString xTitle = "", TString yTitle = "", TString xRange = "", TString yRange = ""TString options = "" );

#endif