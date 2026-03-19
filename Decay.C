#include "DecayFinder.C+"


void Decay(){

    TChain * ch = new TChain("tree");

    // ch -> Add("zzz__309.root");
    ch -> Add("1343LayeredPID2.root");
    //ch -> Add("zzz__246-255.root");

    DecayFinder * selector = new DecayFinder();

    ch -> Process(selector, "");
}