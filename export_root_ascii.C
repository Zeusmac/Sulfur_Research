#include <TFile.h>
#include <TKey.h>
#include <TF1.h>
#include <TDirectory.h>
#include <iostream>
#include <fstream>

void export_root_ascii(const char* filename="myfile.root")
{
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file\n";
        return;
    }

    int npoints = 1000; // resolution of each function

    // Loop over all objects in the file
    TIter next(file->GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj();

        // Check if object is TF1
        if (obj->InheritsFrom("TF1")) {
            TF1 *f = (TF1*)obj;

            std::string fname = f->GetName();
            std::string outname = fname + ".dat";

            std::ofstream fout(outname);

            double xmin = f->GetXmin();
            double xmax = f->GetXmax();
            double step = (xmax - xmin) / npoints;

            for (int i = 0; i <= npoints; i++) {
                double x = xmin + i * step;
                double y = f->Eval(x);
                fout << x << " " << y << "\n";
            }

            fout.close();

            std::cout << "Wrote " << outname << "\n";
        }
    }

    file->Close();
}
