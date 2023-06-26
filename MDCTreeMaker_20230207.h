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

  int truthjet_n;
  int sectorem;
  int sectorih;
  int sectoroh;
  float truthjet_pt[100];
  float truthjet_et[100];
  float truthjet_ph[100];
  int reco_jet_n;
  float reco_jet_pt[100];
  float reco_jet_et[100];
  float reco_jet_ph[100];
  int truthpar_n;
  float truthpar_px[10000];
  float truthpar_py[10000];
  float truthpar_pz[10000];
  int truthpar_n1;
  float truthpar_pt[10000];
  float truthpar_et[10000];
  float truthpar_ph[10000];
  int truthpar_id[10000];
  float truthpar_pt1[10000];
  float truthpar_et1[10000];
  float truthpar_ph1[10000];
  int truthpar_id1[10000];
  int truthpar_j[10000];
  float emfrac[100];
  int truthpar_em[10000];
  int truthpar_em1[10000];
  float emcalen[100000];
  float ihcalen[100000];
  float ohcalen[100000];
  float emcalet[100000];
  float ihcalet[100000];
  float ohcalet[100000];
  float emcalph[100000];
  float ihcalph[100000];
  float ohcalph[100000];
  float prob[10000];
  float chi2[10000];
  int nclus;
};

#endif // MDCTREEMAKER
