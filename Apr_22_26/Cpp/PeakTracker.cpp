#include "PeakTracker.h"
#include <cmath>

double PeakTracker::GetKey(double energy, double tol) {
    for (auto &kv : times) {
        if (fabs(kv.first - energy) < tol)
            return kv.first;
    }

    times[energy] = {};
    counts[energy] = {};
    errors[energy] = {};

    return energy;
}

void PeakTracker::Add(double energy, double time, double counts_val, double error) {
    double key = GetKey(energy);
    times[key].push_back(time);
    counts[key].push_back(counts_val);
    errors[key].push_back(error);
}

const std::map<double, std::vector<double>>& PeakTracker::GetTimes() const { return times; }
const std::map<double, std::vector<double>>& PeakTracker::GetCounts() const { return counts; }
const std::map<double, std::vector<double>>& PeakTracker::GetErrors() const { return errors; }
