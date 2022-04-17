// TCanvas *c_adcmax_istrip_u = new TCanvas("ADCmax vs istrip [U-strips only]", "", 600, 500);
	
// h2_U_ADCmax->Draw();
// h2_U_ADCmax->SetTitle("ADCmax vs istrip [bb.gem.m0.strip.IsU only]");
// h2_U_ADCmax->GetXaxis()->SetTitle("bb.gem.m0.strip.istrip [i]");
// h2_U_ADCmax->GetYaxis()->SetTitle("bb.gem.m0.strip.ADCmax");
// // h2_U_ADCmax->GetXaxis()->SetRange(2944, 3071);
// h2_U_ADCmax->SetMarkerStyle(2);
// h2_U_ADCmax->SetMarkerColor(06);
// c_adcmax_istrip_u->Update();

// TCanvas *c_adcmax_istrip_13 = new TCanvas("ADCmax vs strip - APV 13 - Ustrips only", "c_adcmax_istrip_13", 600, 500);
// h2_APV13_U_ADCmax_istrip->Draw();
// h2_APV13_U_ADCmax_istrip->SetTitle("ADCmax vs strip - APV 13 - Ustrips only");
// h2_APV13_U_ADCmax_istrip->GetXaxis()->SetTitle("istrips");
// h2_APV13_U_ADCmax_istrip->GetYaxis()->SetTitle("ADCmax");
// h2_APV13_U_ADCmax_istrip->SetMarkerStyle(2);
// h2_APV13_U_ADCmax_istrip->SetMarkerColor(06);
// c_adcmax_istrip_13->Update();

//Plot APV 13 --> ADCmax vs Strips on the APV
// TCanvas *c_adcmax_apvstrip_13 = new TCanvas("ADCmax vs strip - APV 13 - APV strips", "c_adcmax_apvstrip_13", 600, 500);
// h2_APV13_U_ADCmax_apvstrip->Draw();
// h2_APV13_U_ADCmax_apvstrip->SetTitle("ADCmax vs APV strip - APV 13 - Ustrips only");
// h2_APV13_U_ADCmax_apvstrip->GetXaxis()->SetTitle("APV Strips");
// h2_APV13_U_ADCmax_apvstrip->GetYaxis()->SetTitle("ADCmax");
// h2_APV13_U_ADCmax_apvstrip->SetMarkerStyle(2);
// h2_APV13_U_ADCmax_apvstrip->SetMarkerColor(06);
// c_adcmax_apvstrip_13->Update();

//Plot APV 13 --> ADCmax vs Channels on the APV
// TCanvas *c_adcmax_apvchan_13 = new TCanvas("ADCmax vs APV channel - APV 13 - U strips", "c_adcmax_apvchan_13", 600, 500);
// h2_APV13_U_ADCmax_apvchan->Draw();
// h2_APV13_U_ADCmax_apvchan->SetTitle("ADCmax vs APV channel - APV 13 - Ustrips only");
// h2_APV13_U_ADCmax_apvchan->GetXaxis()->SetTitle("APV channels");
// h2_APV13_U_ADCmax_apvchan->GetYaxis()->SetTitle("ADCmax");
// h2_APV13_U_ADCmax_apvchan->SetMarkerStyle(2);
// h2_APV13_U_ADCmax_apvchan->SetMarkerColor(06);
// c_adcmax_apvchan_13->Update();

TCanvas *c_adcmax_ratios_apv13 = new TCanvas("Ratio of channels for ADCmax and Ustrips", "c_ratio_apv_13", 600, 500);
h_APV13_U_ADCmax_chan_ratios->Draw();
h_APV13_U_ADCmax_chan_ratios->SetTitle("Ratio of channels ADCmax for U strips");
h_APV13_U_ADCmax_chan_ratios->GetXaxis()->SetTitle("ADC Ratio");
h_APV13_U_ADCmax_chan_ratios->GetYaxis()->SetTitle("Entries");
h_APV13_U_ADCmax_chan_ratios->SetMarkerStyle(2);
h_APV13_U_ADCmax_chan_ratios->SetMarkerColor(06);
c_adcmax_ratios_apv13->Update();

TCanvas *c_APV13_U_ADCmax_neigh_chan= new TCanvas("ADCmax on strip n vs ADCmax on strip n + 1", "c_ratio_apv_13", 600, 500);
h2_APV13_U_ADCmax_neigh_chan->Draw();

TCanvas *c_APV13_U_ADCmax_neigh_chan_lrg_sml = new TCanvas("ADCmax- Larger of 2 strips on x, smaller on y", "", 600, 500);
h2_APV13_U_ADCmax_neigh_chan_lrg_sml->Draw();