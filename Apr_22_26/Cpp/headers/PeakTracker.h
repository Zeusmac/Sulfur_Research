#ifndef PEAKTRACKER_H
#define PEAKTRACKER_H

#include <map>
#include <vector>

class PeakTracker {
public:
    double GetKey(double energy, double tol = 3.0);
    void Add(double energy, double time, double counts, double error, double sigma, double sigma_err);

    const std::map<double, std::vector<double>>& GetTimes() const;
    const std::map<double, std::vector<double>>& GetCounts() const;
    const std::map<double, std::vector<double>>& GetErrors() const;
    const std::map<double, std::vector<double>>& GetSigmas() const;
    const std::map<double, std::vector<double>>& GetSigmaErrors() const;

private:
    std::map<double, std::vector<double>> times;
    std::map<double, std::vector<double>> counts;
    std::map<double, std::vector<double>> errors;
    std::map<double, std::vector<double>> sigmas;
    std::map<double, std::vector<double>> sigma_errors;
};

#endif
