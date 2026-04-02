#ifndef DECAY_FIT_UTILS_H
#define DECAY_FIT_UTILS_H

#include <string>
#include <vector>
#include <utility>
#include <fstream>

#include "TF1.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TH2D.h"

// ==============================
// Config struct (GENERALIZED)
// ==============================
struct FitConfig {
    std::string filename;
    std::string histname;

    double xmin = 0, xmax = 0;
    int rebin = 1;

    // Optional gating
    int gate1_low = -1, gate1_high = -1;
    int gate2_low = -1, gate2_high = -1;

    // Multi-fit support
    std::vector<std::string> names;
    std::vector<double> init;
    std::vector<std::pair<double,double>> bounds;

    // Optional second fit
    std::vector<std::string> names2;
    std::vector<double> init2;
    std::vector<std::pair<double,double>> bounds2;

    // Channel toggles
    bool use_parent = true;
    bool use_daughter = true;
};

// ==============================
// Utility
// ==============================
std::string FormatNumber(double x);
double SafeExp(double x);
double LambdaToHalfLife(double lambda);

// ==============================
// Config
// ==============================
bool ReadConfig(const std::string& fname, FitConfig &cfg);

// ==============================
// ROOT helpers (NEW)
// ==============================
void SetupTF1(TF1* f,
              const std::vector<std::string>& names,
              const std::vector<double>& init,
              const std::vector<std::pair<double,double>>& bounds);

void CopyParameters(TF1* src, TF1* dest);

// ==============================
// Histogram utilities (NEW)
// ==============================
TH1D* MakeProjectionX(TH2D* h2, const std::string& name,
                      int ylow, int yhigh, int rebin);

TH1D* SubtractHistograms(TH1D* h1, TH1D* h2, const std::string& name);

// ==============================
// Fit helpers (NEW)
// ==============================
TFitResultPtr FitHistogram(TH1D* h, TF1* f,
                           double xmin, double xmax,
                           const std::string& options="S R M E");

std::pair<int,int> FindMaxCorrelation(TF1* fit, TFitResultPtr r);

// ==============================
// Results printing
// ==============================
void PrintFitResultsAppend(TF1* fit,
                          TFitResultPtr r,
                          std::ofstream &file,
                          const std::string &label,
                          const std::vector<double> &init,
                          const std::vector<std::pair<double,double>> &bounds);

// ==============================
// Visualization helpers
// ==============================
void DrawFitStatsBox(TF1* fit,
                     const std::vector<int>& lambda_indices,
                     double x1=0.35, double y1=0.65);

// ==============================
// Component builder (NEW)
// ==============================
TF1* MakeComponent(const std::string& name,
                   Double_t (*func)(double*, double*),
                   TF1* fit,
                   int npar,
                   double xmin, double xmax,
                   int color, int style);

#endif
