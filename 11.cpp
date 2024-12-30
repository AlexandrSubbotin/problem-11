#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TMinuit.h"
#include <cmath>
double f1;
static TH1F* h1 = new TH1F("h1", "Data 1;x;Counts", 100, 500, 600);
static TH1F* h2 = new TH1F("h2", "Data 2;x;Counts", 100, 500, 600);
static std::vector<double> data1,data2;

void log_min(int& npar, double* gin, double& f, double* par, int iflag) {
    double a = par[0];
    double m = par[1];
    double s = par[2]; 
    double const1 = par[3];
    
    double func_log = 0.0;
     for (int i = 1; i <= 100; ++i) {
         double model = const1 + a * TMath::Gaus(h1->GetBinCenter(i), m, s);
        if (model <= 0.0) {
             func_log += 2.0e30;
            continue;
        }
         double obs = h1->GetBinContent(i);
       if(obs > 0){
         func_log += -2.0 * (obs * TMath::Log(model / obs) - model + obs);
         } else {
         func_log +=2*model;
         }
        if(std::isnan(func_log)){
             std::cout << "NAN in bin " << i << " in h1" << std::endl;
           }
   }
     for (int i = 1; i <= 100; ++i) {
         double model2 = const1;
          if (model2 <= 0.0) {
             func_log += 2.0e30;
              continue;
           }
        double obs = h2->GetBinContent(i);
           if(obs > 0){
              func_log += -2.0 * (obs * TMath::Log(model2 / obs) - model2 + obs);
           } else {
            func_log+=2*model2;
           }
           if(std::isnan(func_log)){
             std::cout << "NAN in bin " << i << " in h2" << std::endl;
            }
    }
  f = func_log;
    f1=f;
}

void fitTwoHistograms() {
    std::ifstream file1("data_1.dat");
    std::ifstream file2("data_2.dat");
    double value;
     while (file1 >> value) data1.push_back(value);
     while (file2 >> value) data2.push_back(value);
     
    for (double val : data1) h1->Fill(val);
     for (double val : data2) h2->Fill(val);
    
     file1.close();
     file2.close();
    TMinuit minuit(4);
    minuit.SetFCN(log_min);

    minuit.DefineParameter(0, "Nsignal", 100, 10, 0, 1000);
    minuit.DefineParameter(1, "mu", 550, 1, 500, 600);
    minuit.DefineParameter(2, "sigma", 5, 0.5, 0.1, 20);
    minuit.DefineParameter(3, "const1", 10, 1, 0, 100);
   
    minuit.Migrad();

    double a, m, s, const1;
    double error[5];
    minuit.GetParameter(0, a, error[0]);
    minuit.GetParameter(1, m, error[1]);
    minuit.GetParameter(2, s, error[2]);
    minuit.GetParameter(3, const1, error[3]);
    

    TF1* fitFunc1 = new TF1("fitFunc1", "[3] + [0]*TMath::Gaus(x, [1], [2])", 500, 600);
    fitFunc1->SetParameters(a, m, s, const1);
    
    TF1* fitFunc2 = new TF1("fitFunc2", "[0]", 500, 600);
    fitFunc2->SetParameter(0, const1);
  
  
    TCanvas* c1 = new TCanvas("c1", "Histograms with Fit", 800, 600);
    c1->Divide(2,1);
     h1->SetLineColor(kBlue);
    h1->SetTitle("Data 1with Fit");
    h1->GetXaxis()->SetTitle("Value");
    h1->GetYaxis()->SetTitle("Counts");
    c1->cd(1);
    h1->Draw("E");
    h2->SetLineColor(kRed);
    h2->SetTitle("Data 2 with Fit");
    fitFunc1->SetLineColor(kRed);
    fitFunc1->Draw("SAME");
    c1->cd(2);
    h2->Draw("SAME E");
    fitFunc2->SetLineColor(kBlue);
    fitFunc2->Draw("SAME");

    

    TFile *outputFile = new TFile("p10.root", "UPDATE");
    h1->Write();
    h2->Write();
    fitFunc1->Write();
    fitFunc2->Write();
    double N=sqrt(2*M_PI)*s*a;
    double a_error1= sqrt(2*M_PI)*s*error[0];
    std::cout<< "count =  "<< N<< " +-"<< a_error1<<"\n";
    std::cout<< "chi2  "<<f1;
}

int main() {
    fitTwoHistograms();
    return 0;
}