#ifndef GAMMADB_H
#define GAMMADB_H

#include <vector>
#include <string>

struct GammaLine {
    std::string isotope;
    double energy;
};

struct MatchResult {
    std::string isotope;
    double energy;
    double deltaE;
};

class GammaDB {
public:
    bool Load(const std::string& filename);
    std::vector<MatchResult> Match(double energy, double tol = 3.0) const;

private:
    std::vector<GammaLine> db;
};

#endif
