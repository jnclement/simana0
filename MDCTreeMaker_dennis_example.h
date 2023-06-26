#ifndef MDCTREEMAKER_H
#define MDCTREEMAKER_H

#include <fun4all/SubsysReco.h>

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

  int _cluster_n;
  float _cluster_E[1000];
  float _cluster_phi[1000];
  float _cluster_eta[1000];
  int _cluster_ntower[1000];

};

#endif // MDCTREEMAKER
