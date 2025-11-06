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

  TH1F* generateHistograms(int nEvents, int nBins) {
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
    return h_x;
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

  void accordo(int nEvents = 1000, int nBins = 50) { // normalizza HIST e TF1
    TH1F* hist1 = (TH1F*)(generateHistograms(nEvents, nBins)->Clone("hist1"));
    hist1->Scale(1. / hist1->Integral()), "width";

    TF1* cos1     = (TF1*)(f_cos()->Clone("cos1"));
    double cosInt = cos1->Integral(0., 0.6);

    TF1* cosScaled = (TF1*)(f_cos(1 / cosInt)->Clone("cosScaled"));
    std::cout << "Hist integral: " << hist1->Integral() << "\n";
    std::cout << "Cos integral: " << cosScaled->Integral(0., 0.6) << "\n";

    std::vector<double> diff, sigma;
    auto binWidth = 0.6 / nBins;
    for (int i = 0; i < nBins; ++i) {
      double xlow        = hist1->GetBinLowEdge(i + 1);
      double xup         = hist1->GetBinLowEdge(i + 2);
            double cosIntegral = cosScaled->Integral(xlow, xup);
      diff.push_back(cosIntegral - hist1->GetBinContent(i + 1));
      sigma.push_back(cosIntegral);
    }

    double chiSquared;

    for (int i = 0; i < nBins; ++i) {
      chiSquared += std::pow(diff[i] / std::sqrt(sigma[i]), 2);
    }

    std::cout << "Chi quadro: " << chiSquared << "\n";

    TCanvas* c4 = new TCanvas("c4", "Hist scalato", 800, 600);
    hist1->Draw("HIST");
    c4->SaveAs("istrogramma_norm.png");
    TCanvas* c5 = new TCanvas("c5", "Coseno scalato", 800, 600);
    cosScaled->Draw();
    c5->SaveAs("cos_scalato.png");
  }
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

