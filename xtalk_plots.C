// This function makes the plots
TH1F *make_plot(TString run){
  
  
  //TString Rootfile = "/chafs2/work1/sbs/Rootfiles/gmn_replayed_" + run + "_stream0_seg0_0.root";
  //TString Rootfile = "/Users/john/UVa/SBS/inv_mass/gmn_replayed_" + run + "_stream0_seg0_0.root";
  TString Rootfile = "/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/e1209019_replayed_" + run + "_stream0_seg0_0_firstevent0_nevent500000_1.root";

  TChain *t = new TChain("T");

  t->Add(Rootfile);
  
  TH1D *h = new TH1D("m0_ADC_samples","",3500,-500,3000);
  TH1D *u = new TH1D("isU", "", 5000, -2500, 2500);


  
  //TCut cut1 = "bb.gem.track.nhits[0] > 3"; //number of hits on tracks > 3
  //TCut cut2 = "(bb.sh.e + bb.ps.e)/bb.tr.p > 0.78 && (bb.ps.e + bb.ps.e)/bb.tr.p < 1.15";  //For this cut we plot the variable and cut on the peak (see below)
  //TCut cut2 = "(bb.sh.e + bb.ps.e)/bb.tr.p > 0.6 && (bb.ps.e + bb.ps.e)/bb.tr.p < 0.8";

  //TCut cut_tot = cut1 + cut2;
  
  //This calcualtion is W = sqrt(M^2 + 2M(E - E') - 2E*E'(1 - cos(theta))

  //t->Draw(Form("sqrt(%f^2 + 2*%f*(%f - bb.tr.p[0]) - 2*%f*bb.tr.p[0]*(1 - cos(acos(bb.tr.pz[0]/bb.tr.p[0])))) >> hW_I" + current,M_p, M_p, E, E), cut_tot);
  
  t->Draw("bb.gem.m0.strip.rawADCsamples >> m0_ADC_samples");
  //u->Draw("bb.gem.m0.strip.IsU >> IsU");
  
  //t->Draw("(bb.sh.e+bb.ps.e)/bb.tr.p", "bb.ps.e > 0.1");

  // TLorentzVector *TLV = new TLorentzVector(0.0, 0.0, 0.0, 0.938);

  return 0;

}




////// This is the main function ///////
void xtalk_plots(){

  gStyle->SetOptStat(0);

  // Set lines where we are cutting on the elastic peak
  // double low_x = 1.55;
  // double high_x = 2.;


  //TString runs[np] = {"11207","11214","11218","11224","11228"};
  //TString currents[np] = {"1","2","3","4","5"};

  TCanvas *c = new TCanvas("c","",1000,800);

  make_plot("13444");

}
