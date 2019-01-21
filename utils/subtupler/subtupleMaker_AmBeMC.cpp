// C++ Headers
#include <iostream>
#include <vector>
#include <map>
// ROOT Headers
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TFile.h>
#include <TString.h>
#include <Rtypes.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>
// RAT Headers
#include <RAT/SunUtil.hh>
//#include <RAT/DU/ReactorNuOsc.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/ReconCorrector.hh> 
#include <RAT/DU/ReconCalibrator.hh> 
// JSON Reader
#include <fstream>
//#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/json_parser.hpp>
using namespace std;

double shitZshift = -108;

// Order of operations:
// 1. Prune {fitValid, posr<6m}
// 2. Add to output tree: SunR & uDotR -- sun,dxdydz
// 3. Add header tree with: Accepted, rejected
class Parser
{
  public:
    string oname = "default.root";
    vector<string> flist;
    bool isData = false;

    Parser(int argc, char** argv)
    {
      vector<string> args(argv, argv+argc);
      for(int n=1; n<argc; n++)
      {
        if( args[n] == "--out" )
        {
          this->oname = args[n+1];
          n++;
        }
        else if( args[n] == "--data" )
          this->isData = true;
        else
          flist.push_back(args[n]);
      }
      if (flist.size() == 0)
      {
        cout << "Need at least one file" << endl;
        exit(-1);
      }
    }
};

int main(int argc, char** argv)
{
  Parser args(argc, argv); 
   
  // Check the files in argv
  TChain* chain = new TChain("output");
  cerr << "Saving to file " << args.oname << endl;
  cerr<<"loading tchains "<< args.flist.size() - 1 <<endl;
  for(auto& it : args.flist)
    chain->Add(it.c_str());
  //for(int i=1; i<argc; i++)
  //{
  //  chain->Add(argv[i]);
  //}
  // Correction to energy
  const RAT::DU::ReconCorrector &eCorr = RAT::DU::Utility::Get()->GetReconCorrector();
  const RAT::DU::ReconCalibrator &eCalib = RAT::DU::Utility::Get()->GetReconCalibrator();
  TFile* outfile = new TFile(args.oname.c_str(), "recreate");
  // Header tree
  TTree* header = new TTree("header", "MC Header Information");
  int accepted=0;
  double efficiency;
  int total;

  // Indicates how many events were cut with defined selection 
  header->Branch("efficiency", &efficiency);
  header->Branch("total", &total);

  TTree* outtree = chain->CloneTree(0);

  int evIndex;
  double energy, posr;
  double posx, posy, posz;
  double dirx, diry, dirz;
  int uTDays, uTSecs, uTNSecs;
  double itr, beta14;
  ULong64_t dcFlagged;
  double parentKE;
  int runID;

  bool fitValid;
  bool isCal;

  //TString* parentMeta = new TString(256); // Do i need parentMeta1 and 2?
  chain->SetBranchAddress("evIndex", &evIndex);
  chain->SetBranchAddress("runID", &runID);
  chain->SetBranchAddress("energy", &energy);
  chain->SetBranchAddress("posr", &posr);
  chain->SetBranchAddress("fitValid", &fitValid);
  chain->SetBranchAddress("isCal",&isCal); 
  chain->SetBranchAddress("posx", &posx);
  chain->SetBranchAddress("posy", &posy);
  chain->SetBranchAddress("posz", &posz);
  chain->SetBranchAddress("dirx", &dirx);
  chain->SetBranchAddress("diry", &diry);
  chain->SetBranchAddress("dirz", &dirz);
  chain->SetBranchAddress("uTDays", &uTDays);
  chain->SetBranchAddress("uTSecs", &uTSecs);
  chain->SetBranchAddress("uTNSecs", &uTNSecs);
  chain->SetBranchAddress("itr", &itr);
  chain->SetBranchAddress("beta14", &beta14);
  chain->SetBranchAddress("dcFlagged", &dcFlagged);
  // New Branches
  double udotr, sunct, sunx, suny, sunz, posr3;
  int nbc;
  double oscWeight;
  outtree->Branch("posr3", &posr3);
  outtree->Branch("udotr", &udotr);
  outtree->Branch("sunct", &sunct);
  outtree->Branch("sunx", &sunx);
  outtree->Branch("suny", &suny);
  outtree->Branch("sunz", &sunz);
  outtree->Branch("oscillate", &oscWeight);
  outtree->Branch("nbc", &nbc);
  unsigned int entries = chain->GetEntries();
  total = 0;

  cerr<<"Begin processing: "<<entries<<" events"<<endl;
  for(int i=0; i<entries; i++)
  {
    if(i%10000==0)
      cerr<<"Proc: "<<i/double(entries)*100<<"%\r"<<flush;
    chain->GetEvent(i);
    // Oscillate reactors
    total++;
    if( !fitValid )
      continue;
    energy = eCorr.CorrectEnergyRSP(energy);
    accepted++;
    double rho = sqrt( posx*posx + posy*posy);
    bool realData = false;
    energy = eCalib.CalibrateEnergyRSP(realData, energy, rho, posz);
    //Now correct z for AV shift; CalibrateEnergyRSP takes non-corr. Z
    posz = posz + shitZshift;
    posr = sqrt( posx*posx + posy*posy + posz*posz );
    udotr = (posx*dirx+posy*diry+posz*dirz)/sqrt(posx*posx+posy*posy+posz*posz)/sqrt(dirx*dirx+diry*diry+dirz*dirz);
    TVector3 Solar = RAT::SunDirection(uTDays, uTSecs, uTNSecs);
    sunct = (Solar.X()*dirx + Solar.Y()*diry + Solar.Z()*dirz);
    sunx = Solar.X();
    suny = Solar.Y();
    sunz = Solar.Z();
    posr3 = pow((posr/6000.0),3);

    outtree->Fill();

  }
  efficiency = double(accepted)/total;
  header->Fill();
  cerr<<endl<<"fin"<<endl;
  outfile->Write();
}
