// Per compilare: g++ primaprova.cpp `root-config --cflags --libs` -O2 -o primaprova

#include <iostream>
#include <cmath>
#include <vector>
#include "TCanvas.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraphErrors.h"

class Simulation {
private:
    double k_;
    double phi_;
    double b_;
    double xmin_;
    double xmax_;
    int nBins_;

public:
    Simulation(double k, double phi, double b, double xmin, double xmax, int nBins)
        : k_(k), phi_(phi), b_(b), xmin_(xmin), xmax_(xmax), nBins_(nBins) {}

    // Funzione teorica
    double f(double x) const {
        return std::cos(k_ * x + phi_) * std::cos(k_ * x + phi_) + b_;
    }

    // Funzione normalizzata (PDF)
    double f_norm(double x) const {
        TF1 f_int("f_int", [this](double *x, double *) { return this->f(x[0]); },
                  xmin_, xmax_, 0);
        return f(x) / f_int.Integral(xmin_, xmax_);
    }

    // Disegno della funzione normalizzata
    void drawFunctionNorm() const {
        TF1 *f = new TF1("f_norm", [this](double *x, double *) { return this->f_norm(x[0]); },
                         xmin_, xmax_, 0);
        f->SetTitle("f_{norm}(x) = cos^{2}(kx + #phi) + b (normalizzata); x; f(x)");
        f->SetLineColor(kBlue + 1);
        f->SetLineWidth(2);

        TCanvas *c1 = new TCanvas("c1", "Funzione Normalizzata", 800, 600);
        f->GetXaxis()->SetRangeUser(0.0, 5.0);
        f->GetYaxis()->SetRangeUser(0.0, 1.5);
        f->Draw();
        c1->SaveAs("funzionenorm.png");
    }

    // Disegno della funzione non normalizzata
    void drawFunction() const {
        TF1 *f = new TF1("f", [this](double *x, double *) { return this->f(x[0]); },
                         xmin_, xmax_, 0);
        f->SetTitle("f(x) = cos^{2}(kx + #phi) + b; x; f(x)");
        f->SetLineColor(kBlue + 1);
        f->SetLineWidth(2);

        TCanvas *c1 = new TCanvas("c1", "Funzione", 1200, 600);
        f->GetXaxis()->SetRangeUser(0.0, 5.0);
        f->GetYaxis()->SetRangeUser(0.0, 1.5);
        f->Draw();
        c1->SaveAs("funzione.png");
    }

    // Generazione di eventi secondo la distribuzione normalizzata
    void generateEvents(int N) const {
    // Funzione di densità di probabilità (normalizzata)
    TF1 *f_pdf = new TF1("f_pdf", [this](double *x, double *) { return this->f_norm(x[0]); },
                         xmin_, xmax_, 0);
    
    // Istogramma per gli eventi generati
    TH1D *h = new TH1D("h", "Eventi generati secondo f(x);x;Conteggi", nBins_, xmin_, xmax_);
    TRandom3 rand(0);

    for (int i = 0; i < N; ++i) {
        double x = f_pdf->GetRandom(xmin_, xmax_);
        h->Fill(x);
    }

    // Calcolo del fattore di scala per rendere comparabili le due curve
    double scale = h->Integral() * h->GetBinWidth(1);  // area istogramma ~ N_eventi

    // Funzione teorica scalata
    TF1 *f_scaled = new TF1("f_scaled", [this, scale](double *x, double *) {
        return this->f_norm(x[0]) * scale;
    }, xmin_, xmax_, 0);

    // Impostazioni grafiche
    TCanvas *c2 = new TCanvas("c2", "Distribuzione Generata", 1200, 600);
    gStyle->SetOptStat(0);

    h->SetLineColor(kBlue + 1);
    h->SetFillColorAlpha(kAzure + 7, 0.4);
    h->Draw("HIST");

    f_scaled->SetLineColor(kRed);
    f_scaled->SetLineWidth(3);
    f_scaled->SetNpx(1000);  // maggiore risoluzione
    f_scaled->Draw("SAME");

    // Legenda
    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(h, "Distribuzione generata (MC)", "f");
    leg->AddEntry(f_scaled, "Funzione teorica (scalata)", "l");
    leg->Draw();

    c2->SaveAs("distribuzione_generata.png");
}

/*
void studyRegenerationUncertainty(int N = 10000, int B = 50, int nRepeat = 100) const {
    std::cout << "\n>>> [3.2] Studio dell'incertezza da rigenerazione <<<" << std::endl;
    std::cout << "N eventi = " << N << ", bin = " << B << ", ripetizioni = " << nRepeat << std::endl;

    TF1 f_gen("f_gen", [this](double *x, double *) { return this->f_norm(x[0]); },
              xmin_, xmax_, 0);

    // Vettori per media e varianza per bin
    std::vector<double> sum(B, 0.0);
    std::vector<double> sum2(B, 0.0);

    for (int rep = 0; rep < nRepeat; ++rep) {
        TH1D h(Form("h_%d", rep), "Rigenerazioni", B, xmin_, xmax_);
        for (int i = 0; i < N; ++i) {
            double x = f_gen.GetRandom();
            h.Fill(x);
        }

        // Aggiorna somme e somme dei quadrati
        for (int b = 1; b <= B; ++b) {
            double val = h.GetBinContent(b);
            sum[b - 1]  += val;
            sum2[b - 1] += val * val;
        }
    }

    // Calcola media e deviazione standard per bin
    TH1D *h_sigma = new TH1D("h_sigma", "Incertezza da rigenerazione; x; #sigma(bin)", B, xmin_, xmax_);
    for (int b = 1; b <= B; ++b) {
        double mean = sum[b - 1] / nRepeat;
        double var  = (sum2[b - 1] / nRepeat) - mean * mean;
        if (var < 0) var = 0;
        double sigma = sqrt(var);
        h_sigma->SetBinContent(b, sigma);
    }

    // Disegno risultato
    TCanvas *c4 = new TCanvas("c4", "Incertezza da rigenerazione", 900, 600);
    h_sigma->SetLineColor(kRed + 1);
    h_sigma->SetFillColorAlpha(kRed - 7, 0.4);
    h_sigma->Draw("HIST");
    c4->SaveAs("uncertainty_regeneration.png");

    std::cout << "→ Salvato grafico: uncertainty_regeneration.png\n" << std::endl;
} */


// RICCARDO METTI QUA IL TUO CODICE

void binsUncertainty () const {
// Numero di rigenerazioni per stimare incertezza
    int nTrials = 500;
    TRandom3 rnd(0);

    // Salva le somme per calcolare media e sigma per ogni bin
    std::vector<double> sum(nBins, 0.0), sum2(nBins, 0.0);

    for (int t = 0; t < nTrials; ++t) {
        for (int i = 1; i <= nBins; ++i) {
            double y = h_theory.GetBinContent(i);
            // Fluttuazione gaussiana (±10% tipico)
            double y_fluct = rnd.Gaus(y, 0.1 * y);
            sum[i-1]  += y_fluct;
            sum2[i-1] += y_fluct * y_fluct;
        }
    }

    // Calcolo media e sigma per ogni bin
    TGraphErrors *g_unc = new TGraphErrors(nBins);
    for (int i = 1; i <= nBins_; ++i) {
        double mean = sum[i-1] / nTrials;
        double sigma = std::sqrt(sum2[i-1]/nTrials - mean*mean);
        double x = h_theory.GetBinCenter(i);
        g_unc->SetPoint(i-1, x, mean);
        g_unc->SetPointError(i-1, 0, sigma);
    }

    // Disegno
    TCanvas *c = new TCanvas("c", "Bin-smeering", 900, 600);
    h_theory.SetLineColor(kBlue + 1);
    h_theory.Draw("HIST");
    g_unc->SetFillColorAlpha(kRed, 0.3);
    g_unc->Draw("E3 SAME");
    c->SaveAs("bin_smeering.png");

    std::cout << "Salvato grafico: bin_smeering.png" << std::endl;
}

// FINE A QUA RICCARDO

    
};


int main() {
    double k = 5.2;
    double phi = 1.8;
    double b = 0.2;

    Simulation sim(k, phi, b, 0, 2 * M_PI, 100);

    sim.drawFunction();   // Punto 1
    sim.drawFunctionNorm();
    sim.generateEvents(10000); // Punto 2
    // sim.studyRegenerationUncertainty(10000, 50, 200);
    sim.binsUncertainty();

    return 0;
}
