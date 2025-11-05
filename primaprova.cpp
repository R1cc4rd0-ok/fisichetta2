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
    TF1* f_cos = new TF1("Funzione coseno", "[3]*((cos([0]*x + [1]))^2 + [2])", xmin_, xmax_);
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
  TGraph* generateEvents(int n) {
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

  // Prima prova 3. Disegna insieme la funzione teorica e i punti generati casualmente
void drawFunctionAndGenerated(int nEvents = 10000, double norm = 1.) {
    // Ottieni la funzione teorica
    TF1* f = f_cos(norm);
    f->SetLineColor(kRed);
    f->SetLineWidth(3);
    f->SetNpx(1000);

    // Genera punti Monte Carlo
    TGraph* g = generateEvents(nEvents);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.5);
    g->SetMarkerColor(kAzure + 2);

    // Canvas
    TCanvas* c = new TCanvas("c_fun_vs_mc", "Funzione vs Generazione", 900, 600);
    gStyle->SetOptStat(0);

    // Disegno
    g->SetTitle("Confronto: Funzione teorica vs Generazione casuale; x; y");
    g->Draw("AP");  // A = assi, P = punti
    f->Draw("SAME"); // sovrappone la curva teorica

    // Legenda
    TLegend* leg = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg->AddEntry(f, "Funzione teorica f(x)", "l");
    leg->AddEntry(g, Form("Eventi generati (N = %d)", nEvents), "p");
    leg->Draw();

    c->SaveAs("funzione_vs_generazione.png");
}

void studyAgreementNB(const std::vector<int>& N_values,
                      const std::vector<int>& B_values) const {
    std::cout << "\n>>> [Studio dell’accordo al variare di N e B] <<<" << std::endl;

    TF1 f_theory("f_theory", [this](double *x, double *) { return this->f(x[0]); },
                 xmin_, xmax_, 0);
    f_theory.SetNpx(1000);

    for (int N : N_values) {
        for (int B : B_values) {
            TH1D h(Form("h_%d_%d", N, B),
                   Form("Accordo con N=%d, B=%d; x; Conteggi", N, B),
                   B, xmin_, xmax_);

            // Generazione eventi
            TRandom3 rnd(0);
            for (int i = 0; i < N; ++i)
                h.Fill(f_theory.GetRandom());

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
            std::cout << "N=" << N << ", B=" << B << " → χ²/ndf = " << chi2ndf << std::endl;
        }
    }
    std::cout << "----------------------------------------\n";
}

  // Studio incertezza da rigenerazione
  void studyRegenerationUncertainty(int N = 10000, int B = 50,
                                    int nRepeat = 100) const {
    std::cout << "\n>>> Studio dell'incertezza da rigenerazione <<<\n";

    TF1 f_gen(
        "f_gen", [this](double *x, double *) { return this->f_norm(x[0]); },
        xmin_, xmax_, 0);

    std::vector<double> sum(B, 0.0);
    std::vector<double> sum2(B, 0.0);

    for (int rep = 0; rep < nRepeat; ++rep) {
      TH1D h(Form("h_%d", rep), "Rigenerazioni", B, xmin_, xmax_);
      for (int i = 0; i < N; ++i) h.Fill(f_gen.GetRandom());

      for (int b = 1; b <= B; ++b) {
        double val = h.GetBinContent(b);
        sum[b - 1] += val;
        sum2[b - 1] += val * val;
      }
    }

    TH1D *h_sigma =
        new TH1D("h_sigma", "Incertezza da rigenerazione; x; #sigma(bin)", B,
                 xmin_, xmax_);
    for (int b = 1; b <= B; ++b) {
      double mean = sum[b - 1] / nRepeat;
      double var = (sum2[b - 1] / nRepeat) - mean * mean;
      if (var < 0) var = 0;
      h_sigma->SetBinContent(b, std::sqrt(var));
    }

    TCanvas *c = new TCanvas("c_reg", "Incertezza da rigenerazione", 900, 600);
    h_sigma->SetFillColorAlpha(kRed - 7, 0.4);
    h_sigma->Draw("HIST");
    c->SaveAs("uncertainty_regeneration.png");
  }

  // Studio incertezza da Bin-smearing
  void binSmearingUncertainty() const {
    std::cout << "\n>>> Studio dell'incertezza da Bin-smearing <<<\n";

    TH1D *h_func =
        new TH1D("h_func", "f(x) discreta; x; Conteggi", nBins_, xmin_, xmax_);
    for (int b = 1; b <= nBins_; ++b) {
      double x = h_func->GetBinCenter(b);
      h_func->SetBinContent(b, f(x));
    }

    TRandom3 rnd(0);
    TH1D *h_sigma = (TH1D *)h_func->Clone("h_sigma");
    h_sigma->Reset();
    for (int b = 1; b <= nBins_; ++b) {
      double val = h_func->GetBinContent(b);
      double flutt = rnd.Gaus(val, std::sqrt(val + 1e-6));
      h_sigma->SetBinContent(b, std::fabs(flutt - val));
    }

    TCanvas *c = new TCanvas("c_smear", "Incertezza da bin-smearing", 900, 600);
    h_sigma->SetLineColor(kGreen + 2);
    h_sigma->SetFillColorAlpha(kGreen + 1, 0.4);
    h_sigma->Draw("HIST");
    c->SaveAs("uncertainty_binsmearing.png");
  }

  // Propagazione incertezze sui parametri
  void parameterUncertainty(int nTrials = 500) const {
    TRandom3 rnd(0);
    int nPoints = 200;
    std::vector<double> xvals(nPoints);
    std::vector<double> mean(nPoints, 0.0);
    std::vector<double> sum2(nPoints, 0.0);
    double step = (xmax_ - xmin_) / (nPoints - 1);

    for (int i = 0; i < nPoints; ++i) xvals[i] = xmin_ + i * step;

    double sigma_k = 0.02 * k_;
    double sigma_phi = 0.05 * phi_;
    double sigma_b = 0.01 * b_;

    for (int t = 0; t < nTrials; ++t) {
      double k_rand = rnd.Gaus(k_, sigma_k);
      double phi_rand = rnd.Gaus(phi_, sigma_phi);
      double b_rand = rnd.Gaus(b_, sigma_b);
      for (int i = 0; i < nPoints; ++i) {
        double fx = std::cos(k_rand * xvals[i] + phi_rand);
        fx = fx * fx + b_rand;
        mean[i] += fx;
        sum2[i] += fx * fx;
      }
    }

    for (int i = 0; i < nPoints; ++i) mean[i] /= nTrials;

    TGraphErrors *g_unc = new TGraphErrors(nPoints);
    for (int i = 0; i < nPoints; ++i) {
      double sigma = std::sqrt(sum2[i] / nTrials - mean[i] * mean[i]);
      g_unc->SetPoint(i, xvals[i], mean[i]);
      g_unc->SetPointError(i, 0, sigma);
    }

    TF1 *f_central = new TF1(
        "f_central", [this](double *x, double *) { return this->f(x[0]); },
        xmin_, xmax_, 0);

    TCanvas *c =
        new TCanvas("c_param", "Propagazione incertezze parametri", 900, 600);
    f_central->SetLineColor(kBlue + 2);
    f_central->SetLineWidth(2);
    f_central->SetTitle("Propagazione incertezze parametri; x; f(x)");
    f_central->Draw();

    g_unc->SetFillColorAlpha(kRed, 0.3);
    g_unc->Draw("E3 SAME");

    TLegend *leg = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg->AddEntry(f_central, "Funzione nominale", "l");
    leg->AddEntry(g_unc, "Banda d'incertezza", "f");
    leg->Draw();

    c->SaveAs("param_uncertainty.png");
  }
};

// === MAIN ===
int main() {
  double k = 5.2;
  double phi = 1.8;
  double b = 0.2;

  Simulation sim(k, phi, b, 0.0, 0.6, 100);

  sim.drawFunction(1000);
  sim.generateEvents(10000);
  sim.studyRegenerationUncertainty(10000, 50, 200);
  sim.binSmearingUncertainty();
  sim.parameterUncertainty(500);

  return 0;
}
