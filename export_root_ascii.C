#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TClass.h>

#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>

#include <fstream>
#include <iostream>
#include <string>

using std::string;

//////////////////////////////////////////////////////
// Utilities
//////////////////////////////////////////////////////

string sanitize(string name)
{
    for (auto &c : name) {
        if (c == '/' || c == ' ' || c == '.' || c == ':')
            c = '_';
    }
    return name;
}

string makeDir(const string& base, const string& sub)
{
    string path = base + "/" + sub;
    system(("mkdir -p " + path).c_str());
    return path;
}

//////////////////////////////////////////////////////
// TF1
//////////////////////////////////////////////////////

void exportTF1(TF1* f, const string& outdir)
{
    string name = sanitize(f->GetName());

    string dat = outdir + "/" + name + "_tf1.dat";
    string agr = outdir + "/" + name + "_tf1.agr";

    double xmin, xmax;
    f->GetRange(xmin, xmax);

    int npoints = 1000;

    std::ofstream out(dat);

    for (int i = 0; i <= npoints; i++) {
        double x = xmin + (xmax - xmin) * i / npoints;
        out << x << " " << f->Eval(x) << "\n";
    }
    out.close();

    std::ofstream a(agr);
    a << "@g0 on\n";
    a << "@type xy\n";
    a << "@title \"" << name << "\"\n";
    a << "@s0 line type 1\n";
    a << "@s0 symbol 0\n";
    a << "@s0 file \"" << dat << "\"\n";
    a.close();
}

//////////////////////////////////////////////////////
// TGraph
//////////////////////////////////////////////////////

void exportTGraph(TGraph* g, const string& outdir)
{
    string name = sanitize(g->GetName());

    string dat = outdir + "/" + name + "_graph.dat";
    string agr = outdir + "/" + name + "_graph.agr";

    std::ofstream out(dat);

    for (int i = 0; i < g->GetN(); i++) {
        double x, y;
        g->GetPoint(i, x, y);
        out << x << " " << y << "\n";
    }
    out.close();

    std::ofstream a(agr);
    a << "@g0 on\n";
    a << "@type xy\n";
    a << "@title \"" << name << "\"\n";
    a << "@s0 symbol 1\n";
    a << "@s0 line type 1\n";
    a << "@s0 file \"" << dat << "\"\n";
    a.close();
}

//////////////////////////////////////////////////////
// TGraphErrors
//////////////////////////////////////////////////////

void exportTGraphErrors(TGraphErrors* g, const string& outdir)
{
    string name = sanitize(g->GetName());

    string dat = outdir + "/" + name + "_gr.dat";
    string agr = outdir + "/" + name + "_gr.agr";

    std::ofstream out(dat);

    for (int i = 0; i < g->GetN(); i++) {
        double x, y;
        g->GetPoint(i, x, y);
        out << x << " " << y << " "
            << g->GetErrorX(i) << " "
            << g->GetErrorY(i) << "\n";
    }
    out.close();

    std::ofstream a(agr);
    a << "@g0 on\n";
    a << "@type xy\n";
    a << "@title \"" << name << "\"\n";
    a << "@s0 symbol 1\n";
    a << "@s0 line type 0\n";
    a << "@s0 errorbar on\n";
    a << "@s0 file \"" << dat << "\"\n";
    a.close();
}

//////////////////////////////////////////////////////
// TH1D (very useful for gamma spectra)
//////////////////////////////////////////////////////

void exportTH1(TH1* h, const string& outdir)
{
    string name = sanitize(h->GetName());

    string dat = outdir + "/" + name + "_h1.dat";
    string agr = outdir + "/" + name + "_h1.agr";

    std::ofstream out(dat);

    for (int i = 1; i <= h->GetNbinsX(); i++) {
        double x = h->GetBinCenter(i);
        double y = h->GetBinContent(i);
        out << x << " " << y << "\n";
    }
    out.close();

    std::ofstream a(agr);
    a << "@g0 on\n";
    a << "@type xy\n";
    a << "@title \"" << name << "\"\n";
    a << "@s0 symbol 1\n";
    a << "@s0 line type 1\n";
    a << "@s0 file \"" << dat << "\"\n";
    a.close();
}

//////////////////////////////////////////////////////
// Recursive directory walker
//////////////////////////////////////////////////////

void processDir(TDirectory* dir, const string& outdir)
{
    TIter next(dir->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {

        TObject* obj = key->ReadObj();
        string name = sanitize(obj->GetName());

        string subdir = makeDir(outdir, name);

        if (obj->InheritsFrom(TDirectory::Class())) {
            processDir((TDirectory*)obj, subdir);
        }
        else if (obj->InheritsFrom(TF1::Class())) {
            exportTF1((TF1*)obj, subdir);
        }
        else if (obj->InheritsFrom(TGraphErrors::Class())) {
            exportTGraphErrors((TGraphErrors*)obj, subdir);
        }
        else if (obj->InheritsFrom(TGraph::Class())) {
            exportTGraph((TGraph*)obj, subdir);
        }
        else if (obj->InheritsFrom(TH1::Class())) {
            exportTH1((TH1*)obj, subdir);
        }

        delete obj;
    }
}

//////////////////////////////////////////////////////
// MAIN
//////////////////////////////////////////////////////

void export_all_to_grace(const char* filename = "myfile.root")
{
    TFile* file = TFile::Open(filename);

    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file\n";
        return;
    }

    string outdir = "export";
    system("mkdir -p export");

    processDir(file, outdir);

    file->Close();

    std::cout << "Export complete → ./export/\n";
}
