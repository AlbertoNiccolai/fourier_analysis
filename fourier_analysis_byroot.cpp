#include "Evento.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TFile.h"
#include<iostream>
#include<format>
#include<math.h>
#include<vector>
#include <utility> // For std::pair
#include <TMath.h>
#include <cmath>  // For std::fmod

using namespace std;

void fourier_analysis_byroot(){

   // Importo i dati
   Evento e;
   Evento *p = &e;
   double waveform[16][1024], time[16][1024];

   TFile *input = new TFile("/Users/Alberto/Desktop/Magistrale/PSI/MuEDM/BeamTime_Analysis/run00001_analysis/data/run00001.root","READ");
   TTree *tree = (TTree*) input->Get("tree");
   tree->SetBranchAddress("e", &p);
   const int entries = 50000;

   vector<double> voltage_entrance(1024); //All the voltage values in entrance

   tree->GetEntry(48390);//6323, 1156
   for(int i = 0; i < 1024; i++){
      voltage_entrance[i] = e.channel[13].Volt[i];
   }

   //Numero di punti per waveform
   Int_t n=1024;
   Double_t x;

   //Definisco un array di lunghezza (N/2 + 1) * 2 dato che contiene sia la parte reale che immaginaria della trasformazione
   Double_t *in = new Double_t[2*((n+1)/2+1)];
   for (Int_t i=0; i<=n; i++){
      in[i] =  e.channel[6].Volt[i];
   }

   //Faccio la trasformata
   Int_t n_size = n+1;
   TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &n, "R2C ES K");
   if (!fft_own) return;
   fft_own->SetPoints(in);
   fft_own->Transform();

   //Copio l'output
   fft_own->GetPoints(in);

   //Preparo la canvas
   TCanvas *magc = new TCanvas("magc", "Fast Fourier Transform", 800, 600);
   
   TH1 *hr = nullptr;
   hr = TH1::TransformHisto(fft_own, hr, "MAG");
   hr->Scale(1.0 / sqrt(1024));
   double original_min = hr->GetXaxis()->GetXmin();  // Original minimum x value
   double original_max = hr->GetXaxis()->GetXmax();  // Original maximum x value

   // Set the new x-axis limits after scaling by `scale_factor`
   hr->GetXaxis()->SetLimits(original_min* 1.0 / 1.024e-6, original_max * 1. / 1.024e-6);
   hr->SetTitle("Magnitude of the array transform - Event 48390 - Yes signal");
   hr->GetXaxis()->SetLabelSize(0.05);
   hr->GetYaxis()->SetLabelSize(0.05);
   hr->Draw("hist");


   // Define the sampling frequency and calculate the index for 7 MHz
   double Fs = 1e9;  // Sampling frequency (in Hz)
   int N = 1024;     // Number of points in the signal

   const int n_frequencies = 2 * 50;
   double target_frequency[n_frequencies];
   for(int i = 0; i < n_frequencies/2; i++){
      target_frequency[i]= 2e6*i + 1e6;
   }
   for(int i = 0; i < n_frequencies/2; i++){
      target_frequency[i + n_frequencies/2]= 923e6 + 2e6*i + 1e6;
   }

   for(int elem : target_frequency){
      cout<<elem<<endl;
   }


   //const int n_frequencies = 12;
   //double target_frequency[n_frequencies] = {1e06, 3e6, 5e6, 7e6, 9e6, 11e6, 13e6, 15e6, 25e6, 35e6, 54e6, 58e6}; 
   double target_frequency1[n_frequencies];
   for (int i = 0; i<n_frequencies; i++){
      target_frequency1[i] = target_frequency[i] + 1e6; 
   }
   double target_frequency2[n_frequencies];
   for (int i = 0; i<n_frequencies; i++){
      target_frequency2[i] = target_frequency[i] - 1e6; 
   }
   

   // Calculate the index corresponding to the target frequency (7 MHz)
   int freq_index[n_frequencies];
   for (int i = 0; i<n_frequencies; i++){
      freq_index[i] = int(target_frequency[i] * N / Fs);
   }
   int freq_index1[n_frequencies];
   for (int i = 0; i<n_frequencies; i++){
      freq_index1[i] = int(target_frequency1[i] * N / Fs);
   }
   int freq_index2[n_frequencies];
   for (int i = 0; i<n_frequencies; i++){
      freq_index2[i] = int(target_frequency2[i] * N / Fs);
   }

   // Zero out the Fourier coefficients at freq_index and its symmetric component
   for (int i = 0; i < n_frequencies ; i++){
      if (freq_index[i] >= 0 && freq_index[i] < N/2) {
         in[2*freq_index[i]] = 0.002;     // Set real part to 0
         in[2*freq_index[i] + 1] = 0.002; // Set imaginary part to 0
         in[2*freq_index1[i]] = 0.002;     // Set real part to 0
         in[2*freq_index1[i] + 1] = 0.002; // Set imaginary part to 0
         in[2*freq_index2[i]] = 0.002;     // Set real part to 0
         in[2*freq_index2[i] + 1] = 0.002; // Set imaginary part to 0
      }
   }
   
   

   // Perform the inverse FFT to reconstruct the signal
   TVirtualFFT *ifft_own = TVirtualFFT::FFT(1, &N, "C2R M K");
   ifft_own->SetPoints(in);      // Set the modified spectrum
   ifft_own->Transform();        // Perform the inverse FFT

   
   
   // Retrieve the corrected signal (time domain)
   vector<double> re_out(n);
   vector<double> im_out(n);
   for(int i = 0; i< n; i++){
      ifft_own->GetPointComplex(i, re_out[i], im_out[i]);
   }
   
   
   for(int i = 0; i < re_out.size(); i++){
      re_out[i] = double(re_out[i]) / 1024.;
   }

   
   TCanvas *myc = new TCanvas("myc", "Fast Fourier Transform", 1600, 1200);
   myc->Divide(1,2);
   myc->cd(1);
   TGraph *g1 = new TGraph(1024, (double *) time[0], (double *) waveform[0]);
   g1->SetLineColor(kRed);
   g1->SetTitle("Data vs filtered array - Event 48390 ; Time [a.u.] ; Voltage [V]");
   g1->SetLineWidth(1);
   g1->SetLineStyle(2);
   g1->GetXaxis()->SetLimits(0, 1023);
   g1->GetXaxis()->SetLabelSize(0.05);
   g1->GetYaxis()->SetLabelSize(0.05);
   g1->GetYaxis()->SetTitleSize(0.05);
   g1->GetXaxis()->SetTitleSize(0.05);

   for(int j=0;  j<1024; j++){
      g1->SetPoint(j, j+1, e.channel[6].Volt[j+1]);
   }

   g1->Draw("");

   TGraph *g2 = new TGraph(1024, (double *) time[0], (double *) waveform[0]);
   g2->SetLineColor(kBlue+3);
   g2->SetLineWidth(2);
   g2->SetLineStyle(1);
   g2->GetXaxis()->SetLimits(0, 1023);

   for(int j=0;  j<1024; j++){
      g2->SetPoint(j, j, re_out[j]);
   }

   //hi->Draw("hist");
   g2->Draw("same");

   myc->Update();
   

   // Create a legend
   TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);  // Coordinates: (x1, y1, x2, y2)
   legend->AddEntry(g2, "Inverted spectrum", "l");
   legend->AddEntry(g1, "Data array", "l");
   legend->SetBorderSize(1);  // Add a thin border around the legend
   legend->SetFillColor(0);   // Make the legend background transparent
   legend->SetTextSize(0.025); // Set the text size
   legend->Draw();
   
   myc->Update();

   // ********** Third Pad: Plot Residuals **********
   myc->cd(2);  // Move to the third pad for residuals

   // Array for the residuals
   double residual[1024];
   double voltage_for_residual[1024];


   // Calculate the residuals (difference between original and corrected signals)
   for(int i = 0; i < 1023; i++){
      residual[i] = e.channel[6].Volt[i] - re_out[i];  // Subtract reconstructed signal from original
   }

   // Create the residual graph
   TGraph *residual_graph = new TGraph(1024);
   residual_graph->SetLineColor(kBlue);
   residual_graph->SetLineWidth(2);
   residual_graph->SetLineStyle(2);  // Dashed line for residual plot

   // Set residual points into the graph
   for(int i = 0; i < 1024; i++){
      residual_graph->SetPoint(i, i, residual[i]);  // X-axis: sample index, Y-axis: residual
   }

   // Draw the residual graph
   residual_graph->Draw("AL");  // "AL" draws with axes and line

   // Customize the residual plot
   residual_graph->SetTitle("Residuals: Original - Corrected Signal");
   residual_graph->GetXaxis()->SetTitle("Time [a.u.]");
   residual_graph->GetXaxis()->SetLimits(0, 1023);
   residual_graph->GetYaxis()->SetTitle("Residual Amplitude");
   residual_graph->GetXaxis()->SetLabelSize(0.05);
   residual_graph->GetYaxis()->SetLabelSize(0.05);
   residual_graph->GetYaxis()->SetTitleSize(0.05);
   residual_graph->GetXaxis()->SetTitleSize(0.05);

   // Update the canvas to show all pads
   myc->Update();

   
}

int main(){
   fourier_analysis_byroot();
}





