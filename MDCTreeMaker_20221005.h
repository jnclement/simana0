#ifndef MDCTREEMAKER_H
#define MDCTREEMAKER_H

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <string>
#include <vector>

#include "TTree.h"
#include "TFile.h"

class PHCompositeNode;

class MDCTreeMaker : public SubsysReco
{
 public:

  MDCTreeMaker(const std::string &name = "MDCTreeMaker");

  virtual ~MDCTreeMaker();

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;


 private:

  TFile *_f;
  TTree *_tree;
  std::string _foutname;

  int pairp;
  int nlive;
  int nclus;
  float energy[100000];
  float eta[100000];
  float phi[100000];
  float p[3];
  float prob[1000];
  float chi2[1000]; 

};

#endif // MDCTREEMAKER
