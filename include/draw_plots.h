#ifndef draw_plots
#define draw_plots

TCanvas *plot_1DH(TH1D *histogram, const char* c_title = "", const char* title  = "", const char* xTitle = "", const char* yTitle = "", const char* options = "" , double xLow = 0.0, double xHigh = 0.0, double yLow = 0.0, double yHigh = 0.0)
{
	double xRange[2] = {xLow, xHigh};
	double yRange[2] = {yLow, yHigh};
	TCanvas * canvas = new TCanvas(c_title, "", 700, 500);
	histogram->SetTitle(title);
	histogram->GetXaxis()->SetTitle(xTitle);
	histogram->GetYaxis()->SetTitle(yTitle);
	if(TMath::MaxElement(2, xRange) != 0){ histogram->GetXaxis()->SetRange(xRange[0], xRange[1]); }
	if(TMath::MaxElement(2, yRange) != 0){ histogram->GetYaxis()->SetRangeUser(yRange[0], yRange[1]); }
	histogram->Draw(options);

	return canvas;
}

TCanvas * plot_2DH(TH2D *histogram, const char* c_title = "", const char* title  = "", const char* xTitle = "", const char* yTitle = "", const char* options = "" , double xLow = 0.0, double xHigh = 0.0, double yLow = 0.0, double yHigh = 0.0)
{
	double xRange[2] = {xLow, xHigh};
	double yRange[2] = {yLow, yHigh};
	TCanvas * canvas = new TCanvas(c_title, "", 700, 500);
	histogram->SetTitle(title);
	histogram->GetXaxis()->SetTitle(xTitle);
	histogram->GetYaxis()->SetTitle(yTitle);
	if(TMath::MaxElement(2, xRange) != 0){ histogram->GetXaxis()->SetRange(xRange[0], xRange[1]); }
	if(TMath::MaxElement(2, yRange) != 0){ histogram->GetYaxis()->SetRangeUser(yRange[0], yRange[1]); }
	histogram->Draw(options);

	return canvas;
}

#endif