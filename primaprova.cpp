// Compilazione:
// g++ primaprova.cpp `root-config --cflags --libs` -O2 -o primaprova

#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TStyle.h"

class Simulation {
 private:
  double k_;
  double phi_;
  double b_;
  double xmin_;
  double xmax_;
  int nBins_;

 public:
  Simulation(double k, double phi, double b, double xmin, double xmax,
             int nBins)
      : k_(k), phi_(phi), b_(b), xmin_(xmin), xmax_(xmax), nBins_(nBins) {}

  // Funzione teorica
  double f(double x) const {
    return std::cos(k_ * x + phi_) * std::cos(k_ * x + phi_) + b_;
  }

  // Funzione ROOT per disegno
  TF1* f_cos(double norm = 1.) const {
    TF1* f_cos = new TF1("Funzione coseno", "[3]*((cos([0]*x + [1]))^2 + [2])",
                         xmin_, xmax_);
    f_cos->SetParameters(k_, phi_, b_, norm);
    return f_cos;
  }

  // Prima prova 1. Disegna la funzione teorica f_cos
  void drawFunctionTheory(double norm = 1.) const {
    TF1* f = f_cos(norm);

    // Impostazioni grafiche
    f->SetLineColor(kBlue + 2);
    f->SetLineWidth(3);
    f->SetNpx(1000);
    f->SetTitle("Funzione teorica f(x) = cos^{2}(k x + #phi) + b; x; f(x)");

    TCanvas* c = new TCanvas("c_fun_only", "Funzione teorica", 900, 600);
    gStyle->SetOptStat(0);
    f->Draw();

    c->SaveAs("funzione_teorica.png");
  }

  // Prima prova 2. Genera eventi secondo f(x) (hit or miss)
  TGraph* generateEvents(int n, int b) {
    std::vector<double> vx;
    std::vector<double> vy;

    for (int i{0}; i < n; ++i) {
      double x = gRandom->Uniform(0., 0.6);
      double upper_bound = f_cos()->Eval(x);
      double y = gRandom->Uniform(0., 1.2);
      if (y <= upper_bound) {
        vx.push_back(x);
        vy.push_back(y);
      }
    }

    return new TGraph(vx.size(), &vx[0], &vy[0]);
  }

  void generateHistograms(int nEvents, int nBins) {
    // Creazione degli istogrammi
    TH1F* h_x =
        new TH1F(Form("h_x_%d", nEvents), "Distribuzione X", nBins, 0., 0.6);

    // Riempimento degli istogrammi
    for (int i = 0; i < nEvents; ++i) {
      double x = gRandom->Uniform(0., 0.6);
      double upper_bound = f_cos()->Eval(x);
      double y = gRandom->Uniform(0., 1.2);
      if (y <= upper_bound) {
        h_x->Fill(x);
      }
    }

    // Disegna gli istogrammi su un canvas separato
    TCanvas* c1 = new TCanvas(Form("c1_%d", nEvents), "Istogrammi", 800, 600);

    // Disegno dell'istogramma di x (colonne verticali)
    h_x->SetFillColor(kAzure + 2);  // Colore delle colonne
    h_x->Draw();                    // Disegna l'istogramma

    // Salvataggio delle immagini
    c1->SaveAs(Form("istogrammi_%d.png",
                    nEvents));  // Salva il canvas con gli istogrammi

    // Pulizia della memoria
    delete h_x;
    delete c1;
  }

  // Prima prova 5. Disegna insieme la funzione teorica e i punti generati
  // casualmente (hit or miss)
  void drawFunctionAndGenerated(int nEvents = 10000, int nBins = 100,
                                double norm = 1.) {
    // Ottieni la funzione teorica
    TF1* f = f_cos(norm);
    f->SetLineColor(kRed);
    f->SetLineWidth(3);
    f->SetNpx(1000);

    // Genera punti Monte Carlo
    TGraph* g = generateEvents(nEvents, nBins);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.5);
    g->SetMarkerColor(kAzure + 2);

    // Canvas
    TCanvas* c =
        new TCanvas("c_fun_vs_mc", "Funzione vs Generazione", 900, 600);
    gStyle->SetOptStat(0);

    // Disegno
    g->SetTitle("Confronto: Funzione teorica vs Generazione casuale; x; y");
    g->Draw("AP");    // A = assi, P = punti
    f->Draw("SAME");  // sovrappone la curva teorica

    // Legenda
    TLegend* leg = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg->AddEntry(f, "Funzione teorica f(x)", "l");
    leg->AddEntry(g, Form("Eventi generati (N = %d)", nEvents), "p");
    leg->Draw();

    c->SaveAs("funzione_vs_generazione.png");
  }

  /*

  void studyAgreementNB(const std::vector<int>& N_values,
                        const std::vector<int>& B_values) const {
    std::cout << "\n>>> [Studio dell’accordo al variare di N e B] <<<"
              << std::endl;

    TF1 f_theory(
        "f_theory", [this](double* x, double*) { return this->f(x[0]); }, xmin_,
        xmax_, 0);
    f_theory.SetNpx(1000);

    for (int N : N_values) {
      for (int B : B_values) {
        TH1D h(Form("h_%d_%d", N, B),
               Form("Accordo con N=%d, B=%d; x; Conteggi", N, B), B, xmin_,
               xmax_);

        // Generazione eventi
        TRandom3 rnd(0);
        for (int i = 0; i < N; ++i) h.Fill(f_theory.GetRandom());

        // Calcolo chi^2 rispetto alla funzione teorica
        double chi2 = 0;
        for (int i = 1; i <= B; ++i) {
          double x = h.GetBinCenter(i);
          double y_data = h.GetBinContent(i);
          double y_theo = f_theory.Eval(x) * h.GetBinWidth(1) * N;
          double err = std::sqrt(std::max(y_data, 1.0));
          chi2 += std::pow((y_data - y_theo) / err, 2);
        }

        double chi2ndf = chi2 / B;
        std::cout << "N=" << N << ", B=" << B << " → χ²/ndf = " << chi2ndf
                  << std::endl;
      }
    }
    std::cout << "----------------------------------------\n";
  } */
};

// === MAIN ===
int main() {
  double k = 5.2;
  double phi = 1.8;
  double b = 0.2;

  Simulation sim(k, phi, b, 0.0, 0.6, 100);

  sim.drawFunctionTheory();
  sim.drawFunctionAndGenerated();
  sim.generateEvents(10000, 100);
  sim.generateHistograms(10000, 100);

  return 0;
}
