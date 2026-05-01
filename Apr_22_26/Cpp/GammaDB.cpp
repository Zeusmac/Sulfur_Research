#include "GammaDB.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

bool GammaDB::Load(const string& filename) {
    ifstream in(filename);
    if (!in.is_open()) {
        cerr << "ERROR: Cannot open DB file\n";
        return false;
    }

    db.clear();

    string line;
    while (getline(in, line)) {

        if (line.empty()) continue;
        if (line[0] == '#') continue;

        string iso;
        double energy;

        stringstream ss(line);

        if (!(ss >> iso >> energy)) {
            cerr << "Skipping bad line: " << line << endl;
            continue;
        }

        db.push_back({iso, energy});
        cout << "Loaded: " << iso << " " << energy << endl;
    }
    
    cout << "Loaded DB size: " << db.size() << endl;
    return true;
}

vector<MatchResult> GammaDB::Match(double energy, double tol) const {

    vector<MatchResult> matches;

    for (auto &line : db) {
        double dE = fabs(line.energy - energy);
        if (dE <= tol) {
            matches.push_back({line.isotope, line.energy, dE});
        }
    }

    sort(matches.begin(), matches.end(),
         [](const MatchResult &a, const MatchResult &b) {
             return a.deltaE < b.deltaE;
         });

    return matches;
}
