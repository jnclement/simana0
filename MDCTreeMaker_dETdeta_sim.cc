#include "MDCTreeMaker.h"
#include <ffaobjects/EventHeaderv1.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertex.h>
#include <calowaveformsim/WaveformContainerv1.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <phhepmc/PHHepMCGenEventMap.h>
#include <ffaobjects/EventHeaderv1.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <HepMC/GenEvent.h>

#include <jetbackground/TowerBackgroundv1.h>

#include <TLorentzVector.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <iostream>

#include <centrality/CentralityInfov1.h>

using namespace std;
//____________________________________________________________________________..
MDCTreeMaker::MDCTreeMaker(const std::string &name):
  SubsysReco((name+"test").c_str())
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
  
  _tree->Branch("track_vtx",track_vtx,"track_vtx[3]/F");
  _tree->Branch("sectorem",&sectorem,"sectorem/I"); //Number of hit sectors in the emcal                 
  _tree->Branch("sectorih",&sectorih,"sectorih/I"); // IHcal etc.                                        
  _tree->Branch("sectoroh",&sectoroh,"sectoroh/I");
  _tree->Branch("sectormb",&sectormb,"sectormb/I");
  _tree->Branch("mbenrgy",mbenrgy,"mbenrgy[sectormb]/F"); //MBD reported value (could be charge or time)
  _tree->Branch("emcalen",emcalen,"emcalen[sectorem]/F"); //energy per EMCal sector                      
  _tree->Branch("ihcalen",ihcalen,"ihcalen[sectorih]/F"); // per IHCal sector (etc.)                     
  _tree->Branch("ohcalen",ohcalen,"ohcalen[sectoroh]/F");
  _tree->Branch("emcalet",emcalet,"emcalet[sectorem]/I"); //eta of EMCal sector                          
  _tree->Branch("ihcalet",ihcalet,"ihcalet[sectorih]/I");
  _tree->Branch("ohcalet",ohcalet,"ohcalet[sectoroh]/I");
  _tree->Branch("emcalph",emcalph,"emcalph[sectorem]/I"); //phi of EMCal sector                          
  _tree->Branch("ihcalph",ihcalph,"ihcalph[sectorih]/I");
  _tree->Branch("ohcalph",ohcalph,"ihcalph[sectoroh]/I");
  _tree->Branch("mbdtype",mbdtype,"mbdtype[sectormb]/I"); //MBD type (charge or time)                    
  _tree->Branch("mbdside",mbdside,"mbdside[sectormb]/I"); //MBD side (N/S)                               
  _tree->Branch("mbdchan",mbdchan,"mbdchan[sectormb]/I"); //MBD channel number (0-63 on each side)       
  _tree->Branch("emcalt",emcalt,"emcalt[sectorem]/I"); //time value of EMCal sector                      
  _tree->Branch("ihcalt",ihcalt,"ihcalt[sectorih]/I");
  _tree->Branch("ohcalt",ohcalt,"ohcalt[sectoroh]/I");
  _tree->Branch("emcaladc",emcaladc,"emcaladc[sectorem]/I"); //time value of EMCal sector                
  _tree->Branch("ihcaladc",ihcaladc,"ihcaladc[sectorih]/I");
  _tree->Branch("ohcaladc",ohcaladc,"ohcaladc[sectoroh]/I");
  _tree->Branch("ihcalpos",ihcalpos,"ihcalpos[sectorih][3]/F"); //position (xyz) of EMCal sector center  
  _tree->Branch("emcalpos",emcalpos,"emcalpos[sectorem][3]/F");
  _tree->Branch("ohcalpos",ohcalpos,"ohcalpos[sectoroh][3]/F");
  _tree->Branch("emetacor",emetacor,"emetacor[sectorem]/F"); //corrected eta value **NOT A BIN INDEX**   
  _tree->Branch("ihetacor",ihetacor,"ihetacor[sectorih]/F");
  _tree->Branch("ohetacor",ohetacor,"ohetacor[sectoroh]/F");

  /*
  _tree->Branch("track_vtx",track_vtx,"track_vtx[3]/F");
  _tree->Branch("sectorem",&sectorem,"sectorem/I");
  _tree->Branch("sectorih",&sectorih,"sectorih/I");
  _tree->Branch("sectoroh",&sectoroh,"sectoroh/I");
  _tree->Branch("sectormb",&sectormb,"sectormb/I");
  _tree->Branch("mbenrgy",mbenrgy,"mbenrgy[sectormb]/F");
  _tree->Branch("emcalen",emcalen,"emcalen[sectorem]/F");
  _tree->Branch("ihcalen",ihcalen,"ihcalen[sectorih]/F");
  _tree->Branch("ohcalen",ohcalen,"ohcalen[sectoroh]/F");
  _tree->Branch("emcalet",emcalet,"emcalet[sectorem]/I");
  _tree->Branch("ihcalet",ihcalet,"ihcalet[sectorih]/I");
  _tree->Branch("ohcalet",ohcalet,"ohcalet[sectoroh]/I");
  _tree->Branch("emcalph",emcalph,"emcalph[sectorem]/I");
  _tree->Branch("ihcalph",ihcalph,"ihcalph[sectorih]/I");
  _tree->Branch("ohcalph",ohcalph,"ihcalph[sectoroh]/I");
  _tree->Branch("mbdtype",mbdtype,"mbdtype[sectormb]/I");
  _tree->Branch("mbdside",mbdside,"mbdside[sectormb]/I");
  _tree->Branch("mbdchan",mbdchan,"mbdchan[sectormb]/I");
  _tree->Branch("npart",&npart,"npart/I");
  _tree->Branch("ncoll",&ncoll,"ncoll/I");
  _tree->Branch("bimp",&bimp,"bimp/F");
  */
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
  sectorem = 0;
  sectorih = 0;
  sectoroh = 0;
  sectormb = 0;

  TowerInfoContainerv1 *towersEM = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_CEMC");
  TowerInfoContainerv1 *towersIH = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainerv1 *towersOH = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALOUT");

  //TowerInfoContainerv1 *towersEMuc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_CEMC");
  //TowerInfoContainerv1 *towersIHuc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_HCALIN");
  //TowerInfoContainerv1 *towersOHuc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_HCALOUT");

  //RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  //RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  //RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  //TowerInfoContainerv1 *towersMB = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_MBD");
  if(towersEM)
    {
      int nchannels = 24576;
      float emsume = 0;
      //TowerInfoContainerv1::TowerMap begin_end = towersEM->getTowers();
      //for(TowerInfoContainerv1::ConstIter rtiter = begin_end.begin; rtiter != begin_end.end; ++rtiter)
      for(int i=0; i<nchannels; ++i)
	{
	  
	  int skip = 0;
	  TowerInfov1 *tower = towersEM->get_tower_at_channel(i);
	  int key = towersEM->encode_key(i);
	  //int time = towersEMuc->get_tower_at_channel(i)->get_time();
	  //if(time > 10 || time < 7) continue;
	  int etabin = towersEM->getTowerEtaBin(key);
	  int phibin = towersEM->getTowerPhiBin(key);
	  /*
	  if(etabin == 37 && phibin == 0) continue;
	  if(etabin == 47 && phibin == 138) continue;
	  if(etabin == 80 && phibin == 228) continue;
	  */
	  /*
	  for(int j=0; j<hots; ++j)
	    {
	      if(hoteta[j] == etabin && hotpha[j] == phibin) 
		{
		  skip = 1;
		  break;
		}
	    }
	  
	  if(skip) continue;
	  */
	  //float eta = geomEM->get_etacenter(etabin);
	  //float phi = geomEM->get_phicenter(phibin);
	  //cout << "it: " << sector << " eta: " << eta << " eng: " << tower->get_energy() << endl;
	  emcalen[sectorem] = tower->get_energy();
	  //if(emcalen[sectorem] > 20) cout << etabin << " " << phibin << endl;
	  
	  //emsume += tower->get_energy();
	  emcalet[sectorem] = etabin;
	  
	  emcalph[sectorem] = phibin;
	  
	  sectorem++;
	  
	}
    }
  
  if(towersOH)
    {
      int nchannels = 1536;
      float ohsume = 0;
      //TowerInfoContainerv1::TowerMap begin_end = towersEM->getTowers();
      //for(TowerInfoContainerv1::ConstIter rtiter = begin_end.begin; rtiter != begin_end.end; ++rtiter)
      for(int i=0; i<nchannels; ++i)
	{
	  TowerInfov1 *tower = towersOH->get_tower_at_channel(i);
	  //int time = towersOHuc->get_tower_at_channel(i)->get_time();
	  //if(time > 7 || time < 4) continue;
	  int key = towersOH->encode_key(i);
	  int etabin = towersOH->getTowerEtaBin(key);
	  int phibin = towersOH->getTowerPhiBin(key);
	  //float eta = geomOH->get_etacenter(etabin);
	  //float phi = geomOH->get_phicenter(phibin);
	  //cout << "it: " << sector << " eta: " << eta << " eng: " << tower->get_energy() << endl;
	  ohsume += tower->get_energy();
	  ohcalen[sectoroh] = tower->get_energy();
	  ohcalet[sectoroh] = etabin;
	  ohcalph[sectoroh] = phibin;
	  sectoroh++;
	}
    }
  
  if(towersIH)
    {
      int nchannels = 1536;
      //TowerInfoContainerv1::TowerMap begin_end = towersIH->getTowers();
      //for(TowerInfoContainerv1::ConstIter rtiter = begin_end.begin; rtiter != begin_end.end; ++rtiter)
      for(int i=0; i<nchannels; ++i)
	{
	  TowerInfov1 *tower = towersIH->get_tower_at_channel(i);
	  //int time  = towersIHuc->get_tower_at_channel(i)->get_time();
	  //if(time > 7 || time < 4) continue;
	  int key = towersIH->encode_key(i);
	  int etabin = towersIH->getTowerEtaBin(key);
	  int phibin = towersIH->getTowerPhiBin(key);
	  //float eta = geomIH->get_etacenter(etabin);
	  //float phi = geomIH->get_phicenter(phibin);
	  //cout << "it: " << sector << " eta: " << eta << " eng: " << tower->get_energy() << endl;
	  ihcalen[sectorih] = tower->get_energy();
	  ihcalet[sectorih] = etabin;
	  ihcalph[sectorih] = phibin;
	  sectorih++;
	}
    }
  
  /*
  if(towersMB)
    {
      int nchannels = 256;
      float mbsume = 0;
      int printout = 0;
      int onenon = 0;
      //TowerInfoContainerv1::TowerMap begin_end = towersMB->getTowers();
      //for(TowerInfoContainerv1::ConstIter rtiter = begin_end.begin; rtiter != begin_end.end; ++rtiter)
      for(int i=0; i<nchannels; ++i)
	{
	  TowerInfov1 *tower = towersMB->get_tower_at_channel(i);
	  int key = towersMB->encode_key(i);
	  int type = (key >> 6) & 1;
	  int side = (key >> 7) & 1;
	  int chan = key & 0x3f;
	  //cout << "it: " << sector << " eta: " << eta << " eng: " << tower->get_energy() << endl;
	  mbenrgy[sectormb] = tower->get_energy();
	  if(tower->get_energy() == 0) printout = 1;
	  if(tower->get_energy() != 0) onenon = 1;
	  mbdtype[sectormb] = type;
	  mbdside[sectormb] = side;
	  mbdchan[sectormb] = chan;
	  sectormb++;
	  mbsume += tower->get_energy();
	}
      //if(mbsume == 0) cout << "zero mbsume" << endl;
      //if(printout) cout << printout << " " << onenon << endl;
    }
  else
    {
      cout << "no MBD towers!!!" << endl;
    }
  }
  */
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode,"GlobalVertexMap");
  auto iter = vertexmap->begin(); //z vertex getting                                                   
  if(iter != vertexmap->end()) {
    GlobalVertex *vtx = iter->second;
    track_vtx[0] = vtx->get_x();
    track_vtx[1] = vtx->get_y();
    track_vtx[2] = vtx->get_z();
    if (track_vtx[0] == 0 && track_vtx[1] == 0 && track_vtx[2] == 0) { //remove anything with identica\
      lly zero vertex (errors)                                                                                 
        return Fun4AllReturnCodes::EVENT_OK;
    }
    if (fabs(track_vtx[2]) > 30.0) { //cut on 30 - may widen in the future                             
      return Fun4AllReturnCodes::EVENT_OK;
    }
  } else {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  EventHeaderv1 *event_header = findNode::getClass<EventHeaderv1>(topNode, "EventHeader" );
  if ( event_header ) 
    {
      npart = event_header->get_intval("npart");
      ncoll = event_header->get_intval("ncoll");    
      bimp = event_header->get_floatval("bimp");    
      //std::cout << " npart / ncoll / bimp = " << npart << " / " << ncoll << " / " << bimp << std::endl;
    }
  //cout << "post header" << endl;
  _tree->Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;
  
}

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
