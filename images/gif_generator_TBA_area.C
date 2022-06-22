double min_phd=0;
double max_phd=0;
double min_TBA=0;
double max_TBA=0;
double min_rate=0;
double max_rate=0;
double max_time_s=0;

void BinLogX(TH2*h)
{

   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);
   }
   axis->Set(bins, new_bins);
   // delete new_bins;
}



void gif_generator_TBA_area(){
  // TFile *f_data = TFile::Open("SR1_till_today_large_tree_tag.root");
  // TFile *f_data = TFile::Open("kr83.root");
  // TFile *f_data = TFile::Open("SR1_4.3.1.root");
  // TFile *f_data = TFile::Open("5_3_8_5991-6020a.root");


  // TFile *f_data = TFile::Open("SR1_5.4.1.root");
  TFile *f_data = TFile::Open("SR1_5.4.2.root");


  const Int_t fps=1;
  const Int_t gif_length_s=20;
  const Double_t time_window_size_s=5.*60.*60.;
  double max_time_s=-5.*60.*60.;

  // max_time_s=10*60*60; //comment out to keep at default for run

  double check_bin_scale_factor=1000.; //Seconds


  double gif_time=100./(double)fps;

  int nbins=100;

  min_phd=1.5e4;
  max_phd=8.e4;

  //Yrange
  min_TBA=-0.9;
  max_TBA=0;
  //Zrange
  min_rate= 1.e-4;
  max_rate= 1.;

  TF1* corr_S1_TBA =  new TF1("corr_S1_TBA","pol4", min_TBA, max_TBA);
  corr_S1_TBA->SetParameter(0,45430.8);
  corr_S1_TBA->SetParameter(1, -48347.5);
  corr_S1_TBA->SetParameter(2,-152363);
  corr_S1_TBA->SetParameter(3,-257511);
  corr_S1_TBA->SetParameter(4,-123505);



  // This way to calculate livetime only applies to pulses collected with the S1 trigger and does not take into account efficiancy corrections
  TTree *tree= (TTree*) f_data->Get("Event");
  Double_t S1pulse_start_time_ns;
  Double_t first_evt_in_all_runs;
  Double_t runID;
  Double_t last_evt_in_all_runs;
  tree->SetBranchAddress("S1pulse_start_time_ns",&S1pulse_start_time_ns);
  tree->SetBranchAddress("run_id",&runID);
  int evts= tree->GetEntries();
  for (int i=0; i<evts; i++){
    tree->GetEntry(i);
    if (runID>7091) {
      S1pulse_start_time_ns=S1pulse_start_time_ns/1.e9;
      if (i==0) {
        first_evt_in_all_runs=S1pulse_start_time_ns;
        last_evt_in_all_runs=S1pulse_start_time_ns;
        cerr<< "first: "<<first_evt_in_all_runs<<endl;
      }
      if (S1pulse_start_time_ns<first_evt_in_all_runs)first_evt_in_all_runs=S1pulse_start_time_ns;
      if (S1pulse_start_time_ns>last_evt_in_all_runs)last_evt_in_all_runs=S1pulse_start_time_ns;
    }
  }
  if (max_time_s>0)  last_evt_in_all_runs=first_evt_in_all_runs+max_time_s;
  double livetime_s=last_evt_in_all_runs-first_evt_in_all_runs;

  cerr<< "first-last: "<<first_evt_in_all_runs<<"-"<<last_evt_in_all_runs<<endl;

  cerr<< "livetime="<<livetime_s<<" seconds"<<endl;
  const int n_check_bins=trunc(livetime_s/check_bin_scale_factor)+1;
  const int n_check_bins_sqrt=trunc(sqrt(livetime_s/check_bin_scale_factor))+1;
  bool checker_s1[n_check_bins];
  bool checker_rand[n_check_bins];
  bool checker_combo[n_check_bins];
  for (int i=0; i<n_check_bins;i++){
    checker_rand[i]=false;
    checker_s1[i]=false;
    checker_combo[i]=false;
  }

  cerr<< "made it?"<<endl;

  //Instantiate all frames and their properties
  const Int_t total_frames=fps*gif_length_s;
  const Int_t total_hours=ceil(livetime_s/time_window_size_s);
  TH1F* rate_per_h=new TH1F("","",total_hours, 0.,total_hours);
  TH1F* rate_per_h_Rn=new TH1F("","",total_hours, 0.,total_hours);
  TH1F* check_per_h_s1=new TH1F("","",total_hours, 0.,total_hours);
  double earliest_boundary[total_frames];
  double latest_boundary[total_frames];
  double time_bin_width[total_frames];
  TH2F* a_area_vs_tba[total_frames];
  // TH2F* a_larea_vs_tba[total_frames];
  for (int i=0; i<total_frames; i++){
    a_area_vs_tba[i]= new TH2F("", "", nbins, min_phd, max_phd, nbins, min_TBA, max_TBA);
    // a_larea_vs_tba[i]= new TH2F("", "", nbins, 0, 6, nbins, min_TBA, max_TBA);
    // BinLogX(a_larea_vs_tba[i]);
    earliest_boundary[i]=first_evt_in_all_runs;
    latest_boundary[i]=first_evt_in_all_runs+(i+1)*((livetime_s+time_window_size_s)/(double)(total_frames+1));
    if (latest_boundary[i]>first_evt_in_all_runs+time_window_size_s) earliest_boundary[i]=latest_boundary[i]-time_window_size_s;
    if (latest_boundary[i]>last_evt_in_all_runs) latest_boundary[i]= last_evt_in_all_runs;
    time_bin_width[i]=latest_boundary[i]-earliest_boundary[i];
    cerr<<Form("slice %i  %.9fs to  %.9fs (%.9f)",i, earliest_boundary[i]-first_evt_in_all_runs, latest_boundary[i]-first_evt_in_all_runs, time_bin_width[i])<<endl;
  }


  //Access data and fill frames with it
  double avg_rate_Rn_Hz=0;
  Double_t trigger_type;     tree->SetBranchAddress("trigger_type",&trigger_type);
  Double_t pulseTBA;     tree->SetBranchAddress("s1_TBA",&pulseTBA);
  Double_t evtID;     tree->SetBranchAddress("evt_id",&evtID);
  Double_t driftTime_ns;     tree->SetBranchAddress("driftTime_ns",&driftTime_ns);
  Double_t x_cm;     tree->SetBranchAddress("x_cm",&x_cm);
  Double_t y_cm;     tree->SetBranchAddress("y_cm",&y_cm);
  Double_t pulseArea;     tree->SetBranchAddress("s1Area_phd",&pulseArea);
  int entries= tree->GetEntries();
  double percentage_good_time_s1=0.;


  cerr<<endl;
  cerr<<entries<<" pulses"<<endl;
  cerr<<"Hours in rate plot: "<< total_hours<<endl;
  cerr<<"Absolute Pulse Rate: "<< (double)entries/livetime_s <<"Hz"<<endl;
  double prev_ID=0;
  double prev_S1=0;
  for (int i=0; i<entries; i++){
    tree->GetEntry(i);
    if (runID>7091) {
      S1pulse_start_time_ns=S1pulse_start_time_ns/1.e9;
      double hour_ii= (double)(S1pulse_start_time_ns-first_evt_in_all_runs)/60./60.;
      int second_ii=(int) trunc((S1pulse_start_time_ns-first_evt_in_all_runs)/check_bin_scale_factor) ;

      // cerr<< "huh"<< endl;

      // if (pulseArea>400.e3)  cerr<< "ding "<< " pulseArea "<<  pulseArea << " pulseTBA "<<pulseTBA<< " S1pulse_start_time_ns "<<S1pulse_start_time_ns-first_evt_in_all_runs<<endl;
      rate_per_h->Fill(hour_ii);
      if (true) {
        if (pulseArea>1.e1 &&pulseArea<80.e3){
          avg_rate_Rn_Hz++;
        }
        // cerr<<second_ii<<" "<<n_check_bins <<endl;

        if (!checker_s1[ second_ii]) {
          checker_s1[second_ii]=true;
          check_per_h_s1->Fill(hour_ii);
          percentage_good_time_s1+=1./n_check_bins;
        }
      }

      // cerr<< "huh2"<< endl;


      // if (trigger_type==128||trigger_type==129){
      if ((prev_S1!=pulseArea*(60000./corr_S1_TBA->Eval(pulseTBA)) )&& pulseTBA<max_TBA && pulseTBA> min_TBA&& pulseArea<max_phd && pulseArea>min_phd ){
        if (pulseArea>1000 ) pulseArea=pulseArea*(60000./corr_S1_TBA->Eval(pulseTBA));
        prev_S1=pulseArea;
        rate_per_h_Rn->Fill(hour_ii);

        for (int j=0; j<total_frames; j++){
          // if (i<1000) cerr<< S1pulse_start_time_ns << " "<<  earliest_boundary[j] << "-"<<latest_boundary[j]<<endl;
          if (S1pulse_start_time_ns<latest_boundary[j]&&S1pulse_start_time_ns>earliest_boundary[j]){
            // cerr<< "ding "<<j<< " "<<  pulseArea << " "<<pulseTBA<<endl;

            a_area_vs_tba[j]->Fill(pulseArea,pulseTBA);
            // a_larea_vs_tba[j]->Fill(pulseArea,pulseTBA);
          }
        }
      }
    }


    int percent = ceil(((double)i / (double)entries) * 100);
    if (percent % 5 == 0)
    {
        cout << "\r" << std::string(percent/5, '|') << percent << "%";
        cout.flush();
    }



  }
  avg_rate_Rn_Hz=avg_rate_Rn_Hz/livetime_s;

  check_per_h_s1->Scale( check_bin_scale_factor/60./60.);
  rate_per_h_Rn->Divide(check_per_h_s1);

  cerr<< "With a "<<check_bin_scale_factor <<" second window checker, the trigger efficiency is:"<<endl;
  cerr<< percentage_good_time_s1*100.<< "% efficient s1"<< endl;






  TCanvas *canv = new TCanvas("c", "c", 700, 700);
  canv->SetCanvasSize(800, 700);
  canv->SetWindowSize(1000, 800);
  canv->SetLeftMargin(0.15);
  canv->SetRightMargin(0.15);
  canv->SetBottomMargin(0.15);
  canv->SetLogz();
  canv->SetGrid();


  canv->Clear();
  rate_per_h_Rn->Scale(1.e3/60./60.);
  rate_per_h_Rn->Draw("E1");
  rate_per_h_Rn->SetLineColor(kRed);
  rate_per_h_Rn->SetLineWidth(2);
  rate_per_h_Rn->SetStats(0);
  rate_per_h_Rn->SetTitle("LZAP 5.4.2 (7223-7347) - (2h rolling time window)");
  // rate_per_h->SetTitle(Form("%.2e pulses 16-40 kphd (%.2e Hz)",  peak_area_trigger, peak_rate_trigger));
  rate_per_h_Rn->GetYaxis()->SetTitle("Rate [mHz]");
  rate_per_h_Rn->GetXaxis()->SetTitle("Time [h]");
  // rate_per_h_Rn->GetYaxis()->SetRangeUser(0,100);


  // TF1* g1 = new TF1("g1","pol0",0,total_hours-1);
  // rate_per_h_Rn->Fit(g1, "W E");


  TF1* expo2 = new TF1("expo2","[0]*exp(-[1]*x)",0,total_hours-1);
  rate_per_h_Rn->Fit(expo2, "R");
  // time_separation_of_pairs->Integral();
  expo2->Draw("same");
  // cerr<<"Halflife Pair = "<< log(2)/(expo2->GetParameter(1)) << " "<< log(2)/(expo2->GetParameter(1))*((expo2->GetParError(1))/expo2->GetParameter(1)) <<endl;


  // TText *t3 = new TText(.2,.2,Form("Best Fit %.1f (%.1f) mHz", log(2)/(expo2->GetParameter(1)),log(2)/(expo2->GetParameter(1))*((expo2->GetParError(1))/expo2->GetParameter(1))));
  // t3->SetNDC(true);
  // t3->SetTextAlign(11); //align at top left
  // t3->SetTextFont(42);
  // t3->Draw("same");

  canv->Print("output/rate.png");
  canv->Print("output/rate.pdf");
  canv->SetLogy(0);


  //Loop through frames. plot each fram with details on time and progeress bar
  for (int i=0; i<total_frames; i++){
    canv->Clear();

    a_area_vs_tba[i]->Scale(1/time_bin_width[i]);
    a_area_vs_tba[i]->Rebin2D();
    a_area_vs_tba[i]->Draw("colz");
    // gStyle->SetPalette(57);
    a_area_vs_tba[i]->GetYaxis()->SetRangeUser(min_TBA, max_TBA);
    a_area_vs_tba[i]->GetXaxis()->SetRangeUser(min_phd, max_phd);
    if (min_rate!=max_rate)
      a_area_vs_tba[i]->GetZaxis()->SetRangeUser(min_rate, max_rate);
    a_area_vs_tba[i]->SetStats(0);
    a_area_vs_tba[i]->SetTitle("LZAP 5.4.2 (7223-7347) - (2h rolling time window)");
    a_area_vs_tba[i]->GetZaxis()->SetTitle("Rate [Hz]");
    a_area_vs_tba[i]->GetYaxis()->SetTitle("TBA");
    a_area_vs_tba[i]->GetXaxis()->SetTitle("Pulse Area [phd]");
    a_area_vs_tba[i]->SetNdivisions(5,"X");


    canv->Update();

    // cerr<<canv->GetUxmax()<<"  "<<canv->GetUxmin()<<endl;

    //Progress bar
    // TLine *l=new TLine(canv->GetUxmin(),pow(10,canv->GetUymax()),(canv->GetUxmax())*((double)compress*(double)(i+1)/(double)gif_frames),pow(10,canv->GetUymax()));
    double loading_bar_begin=(canv->GetUxmax()-canv->GetUxmin())*(earliest_boundary[i]-first_evt_in_all_runs)/livetime_s;
    double loading_bar_end=(canv->GetUxmax()-canv->GetUxmin())*(latest_boundary[i]-first_evt_in_all_runs)/livetime_s;
    TLine *l=new TLine(loading_bar_begin, canv->GetUymax(),loading_bar_end,canv->GetUymax());
    l->SetLineColor(kRed);
    l->SetLineWidth(7);
    l->Draw("same");


    //Time indication
    double time_begin_m=(earliest_boundary[i]-first_evt_in_all_runs)/60.;
    double time_end_m=(latest_boundary[i]-first_evt_in_all_runs)/60.;
    double time_begin_h=floor(time_begin_m/60.);
    double time_end_h=floor(time_end_m/60.);
    time_begin_m=time_begin_m-time_begin_h*60;
    time_end_m=time_end_m-time_end_h*60;
    TText *t = new TText(.8,.85,Form("%02.0f:%02.0fh-%02.0f:%02.0fh ",time_begin_h,time_begin_m,time_end_h,time_end_m));
    t->SetNDC(true);
    t->SetTextAlign(31); //align at top left
    t->SetTextFont(42);
    t->SetTextColor(kRed);
    t->Draw("same");

    // Latex label
    // double rate_in_frame_Hz=a_area_vs_tba[i]->Integral(a_area_vs_tba[i]->GetXaxis()->FindBin(16.e3),
    //                                                     a_area_vs_tba[i]->GetXaxis()->FindBin(45.e3),
    //                                                     a_area_vs_tba[i]->GetYaxis()->FindBin(-1),
    //                                                     a_area_vs_tba[i]->GetYaxis()->FindBin(1),"");
    // TText *t2 = new TText(.8,.79,Form("%.1f mHz (16-45kphd)", 1000.*rate_in_frame_Hz));
    // t2->SetNDC(true);
    // t2->SetTextAlign(31); //align at top left
    // t2->SetTextFont(42);
    // t2->Draw("same");



    canv->Print(Form("output/data_areas.gif+%i",(int)gif_time));

    // canv->SetLogx();
    // a_larea_vs_tba[i]->Scale(1/time_bin_width[i]);
    // // a_larea_vs_tba[i]->Rebin2D();
    // a_larea_vs_tba[i]->Draw("colz");
    // // gStyle->SetPalette(57);
    // canv->SetLogx(0);
    // canv->Print(Form("output/data_lareas.gif+%i",(int)gif_time));

    if ((earliest_boundary[i]-first_evt_in_all_runs)>max_time_s && max_time_s>0) break;

  }



}
