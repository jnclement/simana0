#include "MDCTreeMaker.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <phhepmc/PHHepMCGenEventMap.h>
#include <ffaobjects/EventHeaderv1.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4VtxPointv1.h>

#include <HepMC/GenEvent.h>

#include <jetbackground/TowerBackgroundv1.h>

#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>

#include <TLorentzVector.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <iostream>
#include <cmath>
#include <centrality/CentralityInfov1.h>

using namespace std;
int some_iterator = 0;
//____________________________________________________________________________..
MDCTreeMaker::MDCTreeMaker(const std::string &name):
  SubsysReco(name)
{
  _foutname = name;  
}

//____________________________________________________________________________..
MDCTreeMaker::~MDCTreeMaker()
{

}

//____________________________________________________________________________..
int MDCTreeMaker::Init(PHCompositeNode *topNode)
{
  
  _f = new TFile( _foutname.c_str(), "RECREATE");

  std::cout << " making a file = " <<  _foutname.c_str() << " , _f = " << _f << std::endl;
  
  _tree = new TTree("ttree","a persevering date tree");

  _tree->Branch("p",p,"p[3]/F");
  _tree->Branch("pairp",&pairp,"pairp/I");
  _tree->Branch("nlive",&nlive,"nlive/I");
  _tree->Branch("nclus",&nclus,"nclus/I");
  _tree->Branch("energy",energy,"energy[nlive]/F");
  _tree->Branch("eta",eta,"eta[nlive]/F");
  _tree->Branch("phi",phi,"phi[nlive]/F");
  _tree->Branch("prob",prob,"prob[nclus]/F");
  _tree->Branch("chi2",chi2,"chi2[nclus]/F");
  _tree->Branch("ohcalet",ohcalet,"ohcalet[sectoroh]/F");
  _tree->Branch("ohcalph",ohcalph,"ohcalph[sectoroh]/F");
  _tree->Branch("ohcalen",ohcalen,"ohcalen[sectoroh]/F");
  _tree->Branch("ihcalet",ihcalet,"ihcalet[sectorih]/F");
  _tree->Branch("ihcalph",ihcalph,"ihcalph[sectorih]/F");
  _tree->Branch("ihcalen",ihcalen,"ihcalen[sectorih]/F");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MDCTreeMaker::InitRun(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int MDCTreeMaker::process_event(PHCompositeNode *topNode)
{

  // here is where the code goes 
  p[0] = 0;
  p[1] = 0;
  p[2] = 0;
  nlive = 0;
  nclus = 0;
  pairp = 0;
  sectoroh = 0;
  sectorih = 0;
  for(int i=0; i<(int)(sizeof(energy)/sizeof(energy[0])); i++)
    {
      energy[i] = 0;
      ihcalen[i] = 0;
      ohcalen[i] = 0;
    }
  {

  RawTowerContainer *towersEM = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
  RawTowerContainer *towersIH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  RawTowerContainer *towersOH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");

  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  /*if(some_iterator == 787)
    {
      cout << "i = 787" << endl;
    }*/
  if(towersEM && geomEM)
    {
      RawTowerContainer::ConstRange begin_end = towersEM->getTowers();
      for(RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
	{
	  RawTower *tower = rtiter->second;

	  RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());
	  float en = tower->get_energy();
	  if(en > 0.04)
	    {
	      energy[nlive] = en;
	      float et = tower_geom->get_eta();
              //float theta = 2*atan(exp(eta));
              //float z = 1/tan(theta);
	      float ph = tower_geom->get_phi();
              //int sector = 0;
              //sector += 256*trunc((z+1.424)/(2.848/96));
              //sector += trunc((phi+M_PI)/(2*M_PI/256));
              eta[nlive] = et;
	      phi[nlive] = ph;
	      nlive++;
	    }
	}
    }
/*
  if(some_iterator == 787)
    {
      cout << "Made it past towers" << endl;
    }
*/
  if(towersIH && geomIH)
    {
      RawTowerContainer::ConstRange begin_end = towersIH->getTowers();
      for(RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
	{
	  RawTower *tower = rtiter->second;

	  RawTowerGeom *tower_geom = geomIH->get_tower_geometry(tower->get_key());
	  float energy = tower->get_energy();
	  if(energy > 0.04)
	    {
	      ihcalen[sectorih] = energy;
	      ihcalet[sectorih] = tower_geom->get_phi();
	      ihcalph[sectorih] = tower_geom->get_phi();
	      sectorih++;
	    }
	}
    }

  if(towersOH && geomOH)
    {
      RawTowerContainer::ConstRange begin_end = towersOH->getTowers();
      for(RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
	{
	  RawTower *tower = rtiter->second;

	  RawTowerGeom *tower_geom = geomOH->get_tower_geometry(tower->get_key());
	  float energy = tower->get_energy();
	  if(energy > 0.04)
	    {
	      ohcalen[sectoroh] = energy;
	      ohcalet[sectoroh] = tower_geom->get_eta();
	      ohcalph[sectoroh] = tower_geom->get_phi();
	      sectoroh++;
	    }
	}
    }
  }

  {
    RawClusterContainer *clusterEM = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");
    if(clusterEM)
      {
	RawClusterContainer::Range begin_end = clusterEM->getClusters();
	for(RawClusterContainer::Iterator iter = begin_end.first; iter != begin_end.second; ++iter)
	  {
	    RawCluster *clust = iter->second;
	    prob[nclus] = clust->get_prob();
	    chi2[nclus] = clust->get_chi2();
	    ++nclus;
	  }
      }
  }
/*
  if(some_iterator == 787)
    {
      cout << "Made it past clusters" << endl;
    }
*/
  {
    int decay_phots[2] = {0,0};
    int pairp_pi0[2] = {0,0};
    int count = 0;
    PHG4TruthInfoContainer* truthinfo = findNode::getClass <PHG4TruthInfoContainer> (topNode, "G4TruthInfo");
    PHG4TruthInfoContainer::Range range = truthinfo->GetParticleRange ();    
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
      {
	PHG4Particle* g4particle = iter->second;
        PHG4VtxPoint *vert = truthinfo->GetVtx(g4particle->get_vtx_id());
	if(vert)
	{
	if(g4particle->get_parent_id() == 1 && g4particle->get_pid() == 22)
	  {
	    decay_phots[count] = g4particle->get_primary_id();
	    if(count == 0) ++count;
	  }
	//if(some_iterator == 787) cout << "In loop" << endl;
	//cout << vert->get_x() << " " << vert->get_y() << endl;
        if(sqrt(vert->get_x()*vert->get_x() + vert->get_y()*vert->get_y()) < 80)
          {
	    if(abs(g4particle->get_pid()) == 11)
	      {
	        if(g4particle->get_parent_id() == 1)
	          {
		    pairp = 1;
		    //cout << "hit " << iter->first << " parent " << g4particle->get_parent_id() << " track " << g4particle->get_track_id() << " primary " << g4particle->get_primary_id() << endl;
	          }
	      }
	  }
	//if(some_iterator == 787) cout << "in loop 2" << endl;
	if(g4particle->get_parent_id() == 0)
	  {
	    p[0] = g4particle->get_px();
	    p[1] = g4particle->get_py();
	    p[2] = g4particle->get_pz();
	  }
	}
	else
	{
	p[0] = 0;
	p[1] = 0;
	p[2] = 0;
	pairp = -1;
	}
      }
/*
    if(some_iterator == 787)
      {
	cout << "Halfway through particles" << endl;
      }
*/
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
      {
        PHG4Particle* g4particle = iter->second;
        PHG4VtxPoint *vert = truthinfo->GetVtx(g4particle->get_vtx_id());
        int vtxid = g4particle->get_vtx_id();
	vert = truthinfo->GetVtx(vtxid);
	if(vert)
	{
	if(abs(g4particle->get_pid()) == 11)
	  {
	    if(g4particle->get_parent_id() == decay_phots[0])
	      {
		if(sqrt(vert->get_x()*vert->get_x() + vert->get_y()*vert->get_y()) < 80)
		  {
		    pairp_pi0[0] = 1;
		  }
	      }
	    if(g4particle->get_parent_id() == decay_phots[1])
              {
                if(sqrt(vert->get_x()*vert->get_x() + vert->get_y()*vert->get_y()) < 80)
                  {
                    pairp_pi0[1] = 1;
                  }
              }
	  }
	}
	else
	{
	pairp = -1;
	}
      }
    if(pairp_pi0[0] && pairp_pi0[1])
      {
	pairp = 1;
      }
  }
  _tree->Fill();
  //cout << "i: " << some_iterator << endl;
  ++some_iterator;
  return Fun4AllReturnCodes::EVENT_OK;
  
}

//int MDCTreeMaker::CreateNode(PHCompositeNode *topNode)
//{

//  return Fun4AllReturnCodes::EVENT_OK;
//}

//____________________________________________________________________________..
int MDCTreeMaker::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
  std::cout << "MDCTreeMaker::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MDCTreeMaker::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
  {
  std::cout << "MDCTreeMaker::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MDCTreeMaker::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
  std::cout << "MDCTreeMaker::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }

  _f->Write();
  _f->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MDCTreeMaker::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
 std::cout << "MDCTreeMaker::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void MDCTreeMaker::Print(const std::string &what) const
{
  std::cout << "MDCTreeMaker::Print(const std::string &what) const Printing info for " << what << std::endl;
}
