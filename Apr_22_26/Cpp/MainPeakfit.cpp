#include "FitGammaAnalyzer.h"
#include <iostream>

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: ./gamma_fit input.root" << std::endl;
        return 1;
    }

    FitGammaAnalyzer analyzer;
    analyzer.Run(argv[1]);

    return 0;
}