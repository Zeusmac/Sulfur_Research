#include "PeakTracker.h"
#include <cmath>

double PeakTracker::GetKey(double energy, double tol) {

    if (!std::isfinite(energy)) return -1;

    for (auto &kv : times) {
        if (std::fabs(kv.first - energy) < tol)
            return kv.first;
    }

    return energy; // DO NOT create anything here
}

void PeakTracker::Add(double energy, double time, double counts_val, double error, double sigma, double sigma_err) {
    if (!std::isfinite(energy) ||
        !std::isfinite(time)   ||
        !std::isfinite(counts_val))
        return;   // 🚨 kill bad data early
    double key = GetKey(energy);

    // If new key → create safely HERE
    if (times.find(key) == times.end()) {
        times[key] = {};
        counts[key] = {};
        errors[key] = {};
        sigmas[key] = {};
        sigma_errors[key] = {};
    }

    times[key].push_back(time);
    counts[key].push_back(counts_val);
    errors[key].push_back(error);
    sigmas[key].push_back(sigma);
    sigma_errors[key].push_back(sigma_err);
}

const std::map<double, std::vector<double>>& PeakTracker::GetTimes() const { return times; }
const std::map<double, std::vector<double>>& PeakTracker::GetCounts() const { return counts; }
const std::map<double, std::vector<double>>& PeakTracker::GetErrors() const { return errors; }
const std::map<double, std::vector<double>>& PeakTracker::GetSigmas() const { return sigmas; }
const std::map<double, std::vector<double>>& PeakTracker::GetSigmaErrors() const { return sigma_errors; }
