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

void incertezzaRigenerazione(int nEvents = 10000, int nBins = 50, int nReps = 1000) {
  std::cout << "Calcolo dell’incertezza da rigenerazione (" << nReps
            << " ripetizioni)..." << std::endl;

  // Array per salvare la somma e la somma dei quadrati per ogni bin
  std::vector<double> sum(nBins, 0.0);
  std::vector<double> sum2(nBins, 0.0);

  // === Ciclo di rigenerazione ===
  for (int r = 0; r < nReps; ++r) {
    TH1F* h = generateHistograms(nEvents, nBins);
    for (int i = 1; i <= nBins; ++i) {
      double val = h->GetBinContent(i);
      sum[i - 1]  += val;
      sum2[i - 1] += val * val;
    }
    delete h;
  }

  // === Calcolo di media e deviazione standard per bin ===
  std::vector<double> mean(nBins), sigma(nBins), x(nBins), ex(nBins, 0.0);

  double binWidth = (xmax_ - xmin_) / nBins;
  for (int i = 0; i < nBins; ++i) {
    mean[i]  = sum[i] / nReps;
    double mean2 = sum2[i] / nReps;
    sigma[i] = std::sqrt(mean2 - mean[i] * mean[i]);
    x[i] = xmin_ + (i + 0.5) * binWidth;
  }

  // === Grafico con barre d’errore ===
  TGraphErrors* g_unc = new TGraphErrors(nBins, &x[0], &mean[0], &ex[0], &sigma[0]);
  g_unc->SetTitle("Incertezza statistica da rigenerazione; x; Conteggio medio");
  g_unc->SetMarkerStyle(20);
  g_unc->SetMarkerColor(kAzure + 2);
  g_unc->SetLineColor(kBlue + 1);

  // === Disegno ===
  TCanvas* c = new TCanvas("c_unc", "Incertezza da rigenerazione", 900, 600);
  g_unc->Draw("AP");
  c->SaveAs("incertezza_rigenerazione.png");

  // === Stampa risultati sintetici ===
  std::cout << "\nBin\t<x>\tConteggio medio\tSigma\n";
  for (int i = 0; i < nBins; ++i) {
    std::cout << i + 1 << "\t" << x[i]
              << "\t" << mean[i]
              << "\t" << sigma[i] << std::endl;
  }
}

void incertezzaBinSmeering(int nBins = 50, int nReps = 100) {
  std::cout << "=== Stima dell'incertezza da Bin-smeering ===" << std::endl;

  // === 1. Funzione teorica normalizzata ===
  TF1* f = f_cos();
  double norm = 1. / f->Integral(xmin_, xmax_);
  TF1* f_scaled = f_cos(norm);

  // === 2. Calcolo del valore teorico medio per bin ===
  std::vector<double> x(nBins), fVal(nBins), sigma(nBins, 0.), ex(nBins, 0.0);
  double binWidth = (xmax_ - xmin_) / nBins;

  for (int i = 0; i < nBins; ++i) {
    double xlow = xmin_ + i * binWidth;
    double xup  = xlow + binWidth;
    x[i] = (xlow + xup) / 2.0;
    fVal[i] = f_scaled->Integral(xlow, xup);  // valore teorico nel bin
  }

  // === 3. Preparazione per le fluttuazioni gaussiane ===
  std::vector<double> sum(nBins, 0.0);
  std::vector<double> sum2(nBins, 0.0);

  // === 4. Esegui nReps fluttuazioni gaussiane per ogni bin ===
  for (int rep = 0; rep < nReps; ++rep) {
    for (int i = 0; i < nBins; ++i) {
      // Definisci l’incertezza relativa, es. 10% oppure sqrt(f_i)
      double sigma_i = std::sqrt(std::max(fVal[i], 1e-8)); // statistica Poissoniana
      double fluctuated = gRandom->Gaus(fVal[i], sigma_i);

      sum[i]  += fluctuated;
      sum2[i] += fluctuated * fluctuated;
    }
  }

  // === 5. Calcolo media e deviazione standard per bin ===
  std::vector<double> mean(nBins);
  for (int i = 0; i < nBins; ++i) {
    mean[i] = sum[i] / nReps;
    double mean2 = sum2[i] / nReps;
    sigma[i] = std::sqrt(mean2 - mean[i] * mean[i]);
  }

  // === 6. Disegno grafico con barre d'errore ===
  TGraphErrors* g_smeering =
      new TGraphErrors(nBins, &x[0], &mean[0], &ex[0], &sigma[0]);
  g_smeering->SetTitle("Incertezza da Bin-smeering; x; Valore medio fluttuato");
  g_smeering->SetMarkerStyle(20);
  g_smeering->SetMarkerColor(kGreen + 2);
  g_smeering->SetLineColor(kGreen + 3);

  TCanvas* c = new TCanvas("c_smeering", "Bin-smeering", 900, 600);
  g_smeering->Draw("AP");
  c->SaveAs("incertezza_bin_smeering.png");

  // === 7. Stampa risultati numerici ===
  std::cout << "\nBin\t<x>\tValore medio\tSigma(fluttuazioni)" << std::endl;
  for (int i = 0; i < nBins; ++i) {
    std::cout << i + 1 << "\t" << x[i] << "\t" << mean[i] << "\t" << sigma[i] << std::endl;
  }

  std::cout << "\nFile salvato: incertezza_bin_smeering.png" << std::endl;
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
  sim.incertezzaRigenerazione(10000, 50, 1000);
  sim.incertezzaBinSmeering(50, 100);
  return 0;
}

