#include <iostream>
#include <cmath>
#include "TCanvas.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"

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
        : k_(k), phi_(phi), b_(b), xmin_(xmin), xmax_(xmax), nBins_(nBins) 
    {
        // Calcolo dell’integrale una sola volta
        TF1 f_int("f_int", [this](double *x, double *) { return this->f(x[0]); },
                  xmin_, xmax_, 0);
        integral_ = f_int.Integral(xmin_, xmax_);
    }

    // Funzione teorica
    double f(double x) const {
        return std::cos(k_ * x + phi_) * std::cos(k_ * x + phi_) + b_;
    }

    // Funzione normalizzata (PDF)
    double f_norm(double x) const {
        return f(x) / integral_;
    }

    // Disegno della funzione normalizzata
    void drawFunction() const {
        TF1 *f = new TF1("f_norm", [this](double *x, double *) { return this->f_norm(x[0]); },
                         xmin_, xmax_, 0);
        f->SetTitle("f_{norm}(x) = cos^{2}(kx + #phi) + b (normalizzata); x; f(x)");
        f->SetLineColor(kBlue + 1);
        f->SetLineWidth(2);

        TCanvas *c1 = new TCanvas("c1", "Funzione Normalizzata", 800, 600);
        f->GetXaxis()->SetRangeUser(0.0, 4.5);
        f->GetYaxis()->SetRangeUser(0.0, 1.2);
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
    TCanvas *c2 = new TCanvas("c2", "Distribuzione Generata", 800, 600);
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
    
};


int main() {
    double k = 5.2;
    double phi = 1.8;
    double b = 0.2;

    Simulation sim(k, phi, b, 0, 2 * M_PI, 100);

    sim.drawFunction();   // Punto 1
    sim.generateEvents(10000); // Punto 2

    return 0;
}
