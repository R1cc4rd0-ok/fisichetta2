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
        : k_(k), phi_(phi), b_(b), xmin_(xmin), xmax_(xmax), nBins_(nBins) {}

     // Funzione teorica
    double f(double x) const {
        return std::cos(k_ * x + phi_) * std::cos(k_ * x + phi_) + b_;
    }

    // Disegno della funzione teorica
    void drawFunction() const {
        TF1 *f = new TF1("f", [this](double *x, double *) { return this->f(x[0]); },
                         xmin_, xmax_, 0);
        f->SetTitle("f(x) = cos^{2}(kx + #phi) + b; x; f(x)");
        f->SetLineColor(kBlue + 1);
        f->SetLineWidth(2);

        TCanvas *c1 = new TCanvas("c1", "Disegno della funzione", 800, 600);
        f->Draw();
        c1->SaveAs("funzione.png");
    }

    // Generazione casuale di eventi secondo la distribuzione
    void generateEvents(int N) const {
        TF1 *f = new TF1("f_gen", [this](double *x, double *) { return this->f(x[0]); },
                         xmin_, xmax_, 0);
        
        TH1D *h = new TH1D("h", "Eventi generati secondo f(x);x;Conteggi", nBins_, xmin_, xmax_);
        TRandom3 rand(0);

        for (int i = 0; i < N; ++i) {
            double x = f->GetRandom(xmin_, xmax_);
            h->Fill(x);
        }

        TCanvas *c2 = new TCanvas("c2", "Distribuzione Generata", 800, 600);
        gStyle->SetOptStat(0);

        f->SetLineColor(kRed);
        f->SetLineStyle(2);
        h->SetLineColor(kBlue + 1);
        h->SetFillColorAlpha(kAzure + 7, 0.4);

        h->Draw("HIST");
        f->Draw("SAME");

        TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
        leg->AddEntry(h, "Distribuzione generata (MC)", "f");
        leg->AddEntry(f, "Funzione teorica", "l");
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
