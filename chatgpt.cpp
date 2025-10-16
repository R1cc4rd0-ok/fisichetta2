// MC_generator_ROOT_didattico.cpp
// Versione didattica per la prova di laboratorio
// Linguaggio: C++ con librerie ROOT
// Obiettivo: generare eventi MC secondo f(x) = cos^2(k x + phi) + b
// Include i punti obbligatori e i due opzionali (hit-or-miss e Metropolis)

// Compilazione (da terminale):
// g++ MC_generator_ROOT_didattico.cpp `root-config --cflags --libs` -O2 -o mc_gen
// Esecuzione:
// ./mc_gen
// oppure in ROOT interactive: root -l MC_generator_ROOT_didattico.cpp

#include <TApplication.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TMath.h>
#include <TFile.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <random>

using namespace std;

// Funzione f(x)
Double_t f_formula(Double_t *x, Double_t *par){
    // par[0]=k  par[1]=phi  par[2]=b
    Double_t xx = x[0];
    Double_t k = par[0];
    Double_t phi = par[1];
    Double_t b = par[2];
    Double_t val = TMath::Cos(k*xx + phi);
    val *= val; // cos^2
    return val + b;
}

// Classe didattica per generatore e analisi
class DistributionGenerator {
public:
    // Parametri base
    double k, phi, b;
    double dk, dphi, db; // relative uncertainties (fractions)
    double xmin, xmax;
    int seed;
    TRandom3 rng;

    DistributionGenerator(double k_, double phi_, double b_,
                          double dk_frac, double dphi_frac, double db_frac,
                          double xmin_=0.0, double xmax_=TMath::Pi(), int seed_=12345)
        : k(k_), phi(phi_), b(b_), dk(dk_frac), dphi(dphi_frac), db(db_frac), xmin(xmin_), xmax(xmax_), seed(seed_), rng(seed_) {
        ;
    }

    // Disegna la funzione con TF1
    TF1* make_TF1(const char* name="f", int npoints=1000){
        TF1 *fun = new TF1(name, f_formula, xmin, xmax, 3);
        fun->SetParameters(k, phi, b);
        return fun;
    }

    // Generazione semplice usando TF1::GetRandom()
    TH1D* generate_TF1_random(int N, int NBins, TF1* fun, const char* hname="h_tf1"){
        TH1D *h = new TH1D(hname, Form("TF1::GetRandom generation N=%d; x; counts", N), NBins, xmin, xmax);
        for(int i=0;i<N;i++){
            double x = fun->GetRandom();
            h->Fill(x);
        }
        return h;
    }

    // Hit-or-miss (accettazione-rifiuto) generation
    TH1D* generate_hit_or_miss(int N, int NBins, TF1* fun, const char* hname="h_hom"){
        TH1D *h = new TH1D(hname, Form("Hit-or-Miss generation N=%d; x; counts", N), NBins, xmin, xmax);
        // assumiamo f_max = 1 + b (cos^2 max = 1)
        double fmax = 1.0 + b;
        int generated = 0;
        while(generated < N){
            double x = rng.Uniform(xmin, xmax);
            double y = rng.Uniform(0.0, fmax);
            double fx = fun->Eval(x);
            if(y <= fx){
                h->Fill(x);
                generated++;
            }
        }
        return h;
    }

    // Metropolis (Markov Chain) generation
    TH1D* generate_metropolis(int N, int NBins, TF1* fun, const char* hname="h_metro", double proposal_sigma=0.1){
        TH1D *h = new TH1D(hname, Form("Metropolis generation N=%d; x; counts", N), NBins, xmin, xmax);
        // inizializziamo a centro intervallo
        double x = 0.5*(xmin + xmax);
        double fx = fun->Eval(x);
        int accepted = 0;
        while(accepted < N){
            double x_prop = x + rng.Gaus(0, proposal_sigma);
            // riflettiamo ai bordi (simple boundary handling)
            if(x_prop < xmin) x_prop = xmin + (xmin - x_prop);
            if(x_prop > xmax) x_prop = xmax - (x_prop - xmax);
            double fx_prop = fun->Eval(x_prop);
            double alpha = 1.0;
            if(fx > 0) alpha = fx_prop / fx;
            if(alpha >= 1.0 || rng.Uniform() < alpha){
                x = x_prop;
                fx = fx_prop;
            }
            h->Fill(x);
            accepted++;
        }
        return h;
    }

    // Calcola incertezza da rigenerazione: ripete la generazione M volte e calcola stddev per bin
    vector<double> uncertainty_from_regeneration(int M, int N, int NBins, TF1* fun, const string &method="TF1"){
        vector<vector<double>> all_counts(NBins, vector<double>(M, 0.0));
        for(int m=0;m<M;m++){
            TH1D *h;
            if(method=="TF1") h = generate_TF1_random(N, NBins, fun, Form("h_tf1_%d", m));
            else if(method=="HOM") h = generate_hit_or_miss(N, NBins, fun, Form("h_hom_%d", m));
            else h = generate_metropolis(N, NBins, fun, Form("h_met_%d", m));
            for(int b=0;b<NBins;b++) all_counts[b][m] = h->GetBinContent(b+1);
            delete h;
        }
        vector<double> stds(NBins, 0.0);
        for(int b=0;b<NBins;b++){
            double mean = accumulate(all_counts[b].begin(), all_counts[b].end(), 0.0) / M;
            double s = 0.0;
            for(int m=0;m<M;m++) s += (all_counts[b][m] - mean)*(all_counts[b][m] - mean);
            s = (M>1) ? sqrt(s/(M-1)) : 0.0;
            stds[b] = s;
        }
        return stds;
    }

    // Calcola incertezza da bin-smearing: fluttua il valore teorico in ciascun bin
    vector<double> uncertainty_from_bin_smearing(int toys, int N, int NBins, TF1* fun){
        // Calcoliamo la normale di integrale per trasformare f(x) in counts
        double integral = fun->Integral(xmin, xmax);
        vector<double> counts_mean(NBins, 0.0);
        double binWidth = (xmax - xmin)/NBins;
        // valore teorico per bin: integrazione della funzione su ogni bin
        for(int i=0;i<NBins;i++){
            double x1 = xmin + i*binWidth;
            double x2 = x1 + binWidth;
            double ibin = fun->Integral(x1, x2);
            counts_mean[i] = N * (ibin / integral);
        }
        // toys: per ogni toy fluttuiamo ciascun bin con Gaussiana centrata su counts_mean[i]
        vector<vector<double>> toy_counts(NBins, vector<double>(toys, 0.0));
        for(int t=0;t<toys;t++){
            for(int i=0;i<NBins;i++){
                double mu = counts_mean[i];
                double sigma = sqrt(mu); // assumiamo poisson-like per sigma iniziale
                if(sigma<=0) sigma = 1.0; // sicurezza
                double val = rng.Gaus(mu, sigma);
                if(val < 0) val = 0;
                toy_counts[i][t] = val;
            }
        }
        vector<double> stds(NBins, 0.0);
        for(int i=0;i<NBins;i++){
            double mean = accumulate(toy_counts[i].begin(), toy_counts[i].end(), 0.0) / toys;
            double s=0.0;
            for(int t=0;t<toys;t++) s += (toy_counts[i][t]-mean)*(toy_counts[i][t]-mean);
            s = (toys>1) ? sqrt(s/(toys-1)) : 0.0;
            stds[i]=s;
        }
        return stds;
    }

    // Propagazione incertezze parametri: fluttua i parametri e ripeti i metodi 3.2 e 3.3
    // Restituisce le incertezze (stddev per bin) ottenute da ripetute generazioni con parametri fluttuati
    vector<double> uncertainty_from_param_fluctuations(int M, int N, int NBins, const string &method="TF1"){
        vector<vector<double>> all_counts(NBins, vector<double>(M, 0.0));
        for(int m=0;m<M;m++){
            // fluttua parametri
            double kp = rng.Gaus(k, fabs(dk*k));
            double phip = rng.Gaus(phi, fabs(dphi*phi));
            double bp = rng.Gaus(b, fabs(db*b));
            TF1 *funp = new TF1(Form("f_p_%d", m), f_formula, xmin, xmax, 3);
            funp->SetParameters(kp, phip, bp);
            TH1D *h;
            if(method=="TF1") h = generate_TF1_random(N, NBins, funp, Form("h_par_tf1_%d", m));
            else if(method=="HOM") h = generate_hit_or_miss(N, NBins, funp, Form("h_par_hom_%d", m));
            else h = generate_metropolis(N, NBins, funp, Form("h_par_met_%d", m));
            for(int bidx=0;bidx<NBins;bidx++) all_counts[bidx][m] = h->GetBinContent(bidx+1);
            delete h; delete funp;
        }
        vector<double> stds(NBins, 0.0);
        for(int bidx=0;bidx<NBins;bidx++){
            double mean = accumulate(all_counts[bidx].begin(), all_counts[bidx].end(), 0.0) / M;
            double s = 0.0;
            for(int m=0;m<M;m++) s += (all_counts[bidx][m]-mean)*(all_counts[bidx][m]-mean);
            s = (M>1) ? sqrt(s/(M-1)) : 0.0;
            stds[bidx] = s;
        }
        return stds;
    }

}; // end class


int main(int argc, char **argv){
    // ROOT application (serve per visualizzazione)
    TApplication app("app", &argc, argv);
    gStyle->SetOptStat(0);

    // Impostazioni iniziali (valori dal file istruzioni)
    double k = 5.2;
    double phi = 1.8;
    double b = 0.2;
    double dk_frac = 0.02; // 2%
    double dphi_frac = 0.05; // 5%
    double db_frac = 0.01; // 1%

    double xmin = 0.0;
    double xmax = TMath::Pi();

    DistributionGenerator gen(k, phi, b, dk_frac, dphi_frac, db_frac, xmin, xmax, 555);

    // Creiamo TF1 e disegniamo la funzione
    TF1 *f = gen.make_TF1("f_main");
    f->SetNpx(1000);
    f->SetLineWidth(2);

    TCanvas *c1 = new TCanvas("c1", "Funzione f(x)", 800, 600);
    f->Draw();

    // Salviamo immagine
    c1->SaveAs("f_function.png");

    // Parametri per le generazioni e gli studi
    int N = 10000; // numero eventi
    int NBins = 50;
    int M = 100; // ripetizioni per incertezza da rigenerazione
    int toys = 1000; // per bin-smearing

    // 2) Generazione semplice TF1::GetRandom
    TH1D *h_tf1 = gen.generate_TF1_random(N, NBins, f, "h_tf1_main");
    // Grafica sovrapposta confronto con funzione teorica (linea)
    TCanvas *c2 = new TCanvas("c2", "Confronto funzione e istogramma", 1000, 600);
    h_tf1->SetLineColor(kBlue);
    h_tf1->SetFillStyle(3004);
    h_tf1->SetFillColorAlpha(kBlue, 0.2);
    h_tf1->Draw("hist");
    // Disegniamo la funzione teorica normalizzata allo stesso numero di eventi
    double integral = f->Integral(xmin, xmax);
    double scale = N / integral;
    TF1 *f_scaled = new TF1("f_scaled", f_formula, xmin, xmax, 3);
    f_scaled->SetParameters(k, phi, b);
    f_scaled->SetNpx(1000);
    f_scaled->SetLineColor(kRed);
    // user function scaled to counts:
    TF1 *f_counts = new TF1("f_counts", [&](double *x, double *p)->Double_t{
        double par[3] = {k, phi, b};
        double val = f_formula(x, par);
        return val * scale;
    }, xmin, xmax, 0);
    f_counts->SetNpx(1000);
    f_counts->Draw("same");

    TLegend *leg = new TLegend(0.6,0.7,0.88,0.88);
    leg->AddEntry(h_tf1, "Istogramma generato", "f");
    leg->AddEntry(f_counts, "Funzione teorica (scalata)", "l");
    leg->Draw();
    c2->SaveAs("comparison_tf1_hist.png");

    // 3.1 Varie N e NBins: esempio rapido
    // (Il codice è didattico: l'utente può variare i parametri a piacere)

    // 3.2 Incertezza da rigenerazione
    cout << "Calcolo incertezza da rigenerazione (TF1) con M="<<M<<" ripetizioni..."<<endl;
    vector<double> stds_reg = gen.uncertainty_from_regeneration(M, N, NBins, f, "TF1");

    // 3.3 Incertezza da bin-smearing
    cout << "Calcolo incertezza da bin-smearing con toys="<<toys<<"..."<<endl;
    vector<double> stds_smear = gen.uncertainty_from_bin_smearing(toys, N, NBins, f);

    // 4. Propagazione incertezze sui parametri
    cout << "Propagazione incertezze sui parametri con fluttuazioni (M="<<M<<")..."<<endl;
    vector<double> stds_param = gen.uncertainty_from_param_fluctuations(M, N, NBins, "TF1");

    // Salviamo i risultati in un file di testo per eventuale analisi
    FILE *out = fopen("uncertainties_summary.txt", "w");
    fprintf(out, "BinCenter\tStdRegeneration\tStdBinSmear\tStdParamFluct\n");
    double binW = (xmax - xmin)/NBins;
    for(int i=0;i<NBins;i++){
        double center = xmin + (i+0.5)*binW;
        fprintf(out, "%g\t%g\t%g\t%g\n", center, stds_reg[i], stds_smear[i], stds_param[i]);
    }
    fclose(out);

    cout << "Risultati scritti in uncertainties_summary.txt e immagini PNG generate."<<endl;

    // 5. OPZIONALE 1: Hit-or-miss e confronto
    cout << "Eseguo generazione Hit-or-Miss (facoltativa)..."<<endl;
    TH1D *h_hom = gen.generate_hit_or_miss(N, NBins, f, "h_hom_main");
    TCanvas *c3 = new TCanvas("c3", "Hit-or-miss vs TF1", 1000, 600);
    h_hom->SetLineColor(kGreen+2);
    h_hom->Draw("hist");
    f_counts->Draw("same");
    c3->SaveAs("hit_or_miss_vs_tf1.png");

    // 6. OPZIONALE 2: Metropolis
    cout << "Eseguo generazione Metropolis (facoltativa)..."<<endl;
    TH1D *h_met = gen.generate_metropolis(N, NBins, f, "h_met_main", 0.05);
    TCanvas *c4 = new TCanvas("c4", "Metropolis vs TF1", 1000, 600);
    h_met->SetLineColor(kMagenta);
    h_met->Draw("hist");
    f_counts->Draw("same");
    c4->SaveAs("metropolis_vs_tf1.png");

    // 7. Confronti finali (semplice esempio di output numerico):
    // confronto tra std medi delle tre incertezze
    double mean_reg = accumulate(stds_reg.begin(), stds_reg.end(), 0.0)/NBins;
    double mean_sm = accumulate(stds_smear.begin(), stds_smear.end(), 0.0)/NBins;
    double mean_par = accumulate(stds_param.begin(), stds_param.end(), 0.0)/NBins;
    cout << "Media delle incertezze per bin:\n"
         << " - Rigenerazione: "<<mean_reg<<"\n"
         << " - Bin-smearing:  "<<mean_sm<<"\n"
         << " - Param-fluct:   "<<mean_par<<"\n";

    // Salviamo gli istogrammi in file ROOT
    TFile *fout = new TFile("mc_output.root","RECREATE");
    h_tf1->Write();
    h_hom->Write();
    h_met->Write();
    f_counts->Write("f_counts_as_histfunc");
    fout->Close();
    cout << "File ROOT mc_output.root scritto con gli istogrammi."<<endl;

    // Mostriamo l'app grafica fino alla chiusura manuale
    cout << "Tutte le operazioni completate. Le immagini PNG e il file ROOT sono nella cartella corrente."<<endl;
    cout << "Chiudere le finestre grafiche per terminare il programma."<<endl;

    app.Run();
    return 0;
}

