#ifndef MDCTREEMAKER_H
#define MDCTREEMAKER_H

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>
#include <string>
#include <vector>
#include <calowaveformsim/WaveformContainerv1.h>
#include "TTree.h"
#include "TFile.h"

class PHCompositeNode;
class MDCTreeMaker : public SubsysReco
{
 public:

  MDCTreeMaker(const std::string &name = "MDCTreeMaker", const int dataormc = 0, const int debug = 1, const int correct = 1);

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
  int _correct;
  int _debug;
  int _dataormc;
  TFile *_f;
  TTree *_tree;
  std::string _foutname;
  int sectorem;
  int sectorih;
  int sectoroh;
  int sectormb;
  float emcalen[24576];
  float ihcalen[1536];
  float ohcalen[1536];
  int emcalt[24576];
  int ihcalt[1536];
  int ohcalt[1536];
  int emcaladc[24576];
  int ihcaladc[1536];
  int ohcaladc[1536];
  float emcalpos[24576][3];
  float ihcalpos[1536][3];
  float ohcalpos[1536][3];
  int emcaletabin[24576];
  int ihcaletabin[1536];
  int ohcaletabin[1536];
  int emcalphibin[24576];
  int ihcalphibin[1536];
  int ohcalphibin[1536];
  float mbenrgy[256];
  int mbdtype[256];
  int mbdside[256];
  int mbdchan[256];
  int npart;
  int ncoll;
  float bimp;
  int truth_vertices;
  float track_vtx[3];
  float svtx_vtx[3];
  float mbd_vtx[3];
  float truth_vtx[3];
  int vtx_id;
  float emetacor[24576];
  float ihetacor[1536];
  float ohetacor[1536];
};

#endif // MDCTREEMAKER
