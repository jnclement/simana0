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
#include <g4main/PHG4VtxPoint.h>
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
#include <cmath>

#include <TLorentzVector.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <iostream>

#include <centrality/CentralityInfov1.h>

using namespace std;
//____________________________________________________________________________..
MDCTreeMaker::MDCTreeMaker(const std::string &name, const int dataormc, const int debug):
  SubsysReco((name+"test").c_str())
{
  _foutname = name;  
  _dataormc = dataormc;
  _debug = debug;
}

//____________________________________________________________________________..
MDCTreeMaker::~MDCTreeMaker()
{

}

//____________________________________________________________________________..
int MDCTreeMaker::Init(PHCompositeNode *topNode)
{


  _f = new TFile( _foutname.c_str(), "RECREATE");
  
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
  _tree->Branch("emcaletabin",emcaletabin,"emcaletabin[sectorem]/I"); //eta of EMCal sector
  _tree->Branch("ihcaletabin",ihcaletabin,"ihcaletabin[sectorih]/I");
  _tree->Branch("ohcaletabin",ohcaletabin,"ohcaletabin[sectoroh]/I");
  _tree->Branch("emcalphibin",emcalphibin,"emcalphibin[sectorem]/I"); //phi of EMCal sector
  _tree->Branch("ihcalphibin",ihcalphibin,"ihcalphibin[sectorih]/I");
  _tree->Branch("ohcalphibin",ohcalphibin,"ihcalphibin[sectoroh]/I");
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

  if(_dataormc) //get collision parameters for MC
    {
      _tree->Branch("truth_vertices",&truth_vertices,"truth_vertices/I");
      _tree->Branch("npart",&npart,"npart/I");
      _tree->Branch("ncoll",&ncoll,"ncoll/I");
      _tree->Branch("bimp",&bimp,"bimp/F");
      _tree->Branch("truth_vtx",truth_vtx,"truth_vtx[truth_vertices][3]/F");
      _tree->Branch("vtx_id",vtx_id,"vtx_id[truth_vertices]/I");
    }
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
  if(_debug) cout << "Beginning event processing" << endl;
  //reset lengths to all zero
  sectorem = 0;
  sectorih = 0;
  sectoroh = 0;
  sectormb = 0;
  /*
  for(int i=0; i<sizeof(emcalen)/sizeof(emcalen[0]); i++)
    {
      emcalen[i] = -1;
      emcalet[i] = 0;
      emcalph[i] = 0;
      ihcalen[i] = -1;
      ihcalet[i] = 0;
      ihcalph[i] = 0;
      ohcalen[i] = -1;
      ohcalet[i] = 0;
      ohcalph[i] = 0;
      mbenrgy[i] = -1;
      mbdtype[i] = -1;
      mbdside[i] = -1;
      mbdchan[i] = -1;
    }
  */
  {
     
    //Get towerinfocontainer objects from nodetree
    TowerInfoContainerv1 *towersEM = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_CEMC");
    TowerInfoContainerv1 *towersIH = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALIN");
    TowerInfoContainerv1 *towersOH = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALOUT");
    TowerInfoContainerv1 *towersEMuc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_CEMC");
    TowerInfoContainerv1 *towersIHuc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_HCALIN");
    TowerInfoContainerv1 *towersOHuc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_HCALOUT");
    //Geometry objects for determining eta and phi locations
    RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

    TowerInfoContainerv1 *towersMB = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_MBD");
    
    GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode,"GlobalVertexMap");
    if(_debug) cout << "Checking that all necessary objects exist" << endl;
    if(!towersEM || !towersIH || !towersOH || !towersEMuc || !towersIHuc || !towersOHuc || !geomEM || !geomIH || !geomOH || !vertexmap) return Fun4AllReturnCodes::EVENT_OK; //remove events which do not have all required information
    if(_debug) cout << "EM geomtry node: " << geomEM << endl;
    if(_debug) cout << "Getting vertex" << endl;
    auto iter = vertexmap->begin(); //z vertex getting
    if(iter != vertexmap->end()) {
      GlobalVertex *vtx = iter->second;
      track_vtx[0] = vtx->get_x();
      track_vtx[1] = vtx->get_y();
      track_vtx[2] = vtx->get_z();
      if (track_vtx[0] == 0 && track_vtx[1] == 0 && track_vtx[2] == 0) { //remove anything with identically zero vertex (errors)
	return Fun4AllReturnCodes::EVENT_OK;
      }
      if (fabs(track_vtx[2]) > 30.0) { //cut on 30 - may widen in the future
	return Fun4AllReturnCodes::EVENT_OK;
      }
    } else {
      return Fun4AllReturnCodes::EVENT_OK;
    }
    if(_debug) cout << "Getting EMCal info" << endl;
    if(towersEM) //get EMCal values
      {
	int nchannels = 24576; //channels in emcal
	for(int i=0; i<nchannels; ++i) //loop over channels
	  {
	    TowerInfov1 *tower = towersEM->get_tower_at_channel(i); //get EMCal tower
	    if(_debug) cout << "Tower " << i << ": " << tower << endl;
 	    int key = towersEM->encode_key(i);
	    int time = towersEM->get_tower_at_channel(i)->get_time(); //get uncalibrated tower
	    if(time > 7 || time < 5) continue; //timing cut
	    int etabin = towersEM->getTowerEtaBin(key); //get eta and phi indices
	    int phibin = towersEM->getTowerPhiBin(key);
	    emcalen[sectorem] = tower->get_energy(); //actual tower energy (calibrated)
	    const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, etabin, phibin);
	    RawTowerGeom *tower_geom = geomEM->get_tower_geometry(geomkey); //encode tower geometry
	    if(_debug) cout << "Tower geom: " << tower_geom << endl;
	    emcalt[sectorem] = time; //store time value
	    emcaladc[sectorem] = towersEMuc->get_tower_at_channel(i)->get_energy(); //emcal ADC value (uncalibrated "energy")
	    emcalpos[sectorem][0] = tower_geom->get_center_x(); //get positions of towers
	    emcalpos[sectorem][1] = tower_geom->get_center_y();
	    emcalpos[sectorem][2] = tower_geom->get_center_z();
	    double newz = emcalpos[sectorem][2] - track_vtx[2];
	    double newx = emcalpos[sectorem][0] - track_vtx[1];
	    double newy = emcalpos[sectorem][1] - track_vtx[0];
	    emetacor[sectorem] = asinh(newz/sqrt(newx*newx+newy*newy));
	    emcaletabin[sectorem] = etabin; //get eta and phi of towers
	    emcalphibin[sectorem] = phibin;
	    sectorem++;
	  }
      }
    if(_debug) cout << "Getting OHCal info" << endl;
    if(towersOH) //essentially the same as for EMCal, just with fewer sectors and named OH
      {
	int nchannels = 1536;
	for(int i=0; i<nchannels; ++i)
	  {
	    TowerInfov1 *tower = towersOH->get_tower_at_channel(i);
	    int time = towersOH->get_tower_at_channel(i)->get_time();
	    
	    if(time > 8 || time < 5) continue;
	    int key = towersOH->encode_key(i);
	    int etabin = towersOH->getTowerEtaBin(key);
	    int phibin = towersOH->getTowerPhiBin(key);
	    ohcalen[sectoroh] = tower->get_energy();
	    const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, etabin, phibin);
	    RawTowerGeom *tower_geom = geomOH->get_tower_geometry(geomkey);
	    ohcalt[sectoroh] = time;
	    ohcaladc[sectoroh] = towersOHuc->get_tower_at_channel(i)->get_energy();
	    ohcalpos[sectoroh][0] = tower_geom->get_center_x();
	    ohcalpos[sectoroh][1] = tower_geom->get_center_y();
	    ohcalpos[sectoroh][2] = tower_geom->get_center_z();
	    double newz = ohcalpos[sectoroh][2] - track_vtx[2];
	    double newx = ohcalpos[sectoroh][0] - track_vtx[1];
	    double newy = ohcalpos[sectoroh][1] - track_vtx[0];
	    ohetacor[sectoroh] = asinh(newz/sqrt(newx*newx+newy*newy));

	    ohcaletabin[sectoroh] = etabin;
	    ohcalphibin[sectoroh] = phibin;
	    sectoroh++;
	  }
      }
    if(_debug) cout << "Getting IHCal info" << endl;
    if(towersIH) //same
      {
	int nchannels = 1536;
	for(int i=0; i<nchannels; ++i)
	  {
	    TowerInfov1 *tower = towersIH->get_tower_at_channel(i);
	    int time = towersIH->get_tower_at_channel(i)->get_time();

	    if(time > 7 || time < 5) continue;
	    int key = towersIH->encode_key(i);
	    int etabin = towersIH->getTowerEtaBin(key);
	    int phibin = towersIH->getTowerPhiBin(key);
	    ihcalen[sectorih] = tower->get_energy();
	    const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, etabin, phibin);
	    RawTowerGeom *tower_geom = geomIH->get_tower_geometry(geomkey);
	    ihcalt[sectorih] = time;
	    ihcaladc[sectorih] = towersIHuc->get_tower_at_channel(i)->get_energy();
	    ihcalpos[sectorih][0] = tower_geom->get_center_x();
	    ihcalpos[sectorih][1] = tower_geom->get_center_y();
	    ihcalpos[sectorih][2] = tower_geom->get_center_z();
	    double newz = ihcalpos[sectorih][2] - track_vtx[2];
	    double newx = ihcalpos[sectorih][0] - track_vtx[1];
	    double newy = ihcalpos[sectorih][1] - track_vtx[0];
	    ihetacor[sectorih] = asinh(newz/sqrt(newx*newx+newy*newy));

	    ihcaletabin[sectorih] = etabin;
	    ihcalphibin[sectorih] = phibin;
	    sectorih++;
	  }
      }
    if(_debug) cout << "Getting MBD info" << endl;
    if(towersMB) //get MBD info
      {
	int nchannels = 256;
	for(int i=0; i<nchannels; ++i)
	  {
	    TowerInfov1 *tower = towersMB->get_tower_at_channel(i); //get mbd channel
	    int key = towersMB->encode_key(i);
	    int type = (key >> 6) & 1; //encodes type from key (only way to get it)
	    int side = (key >> 7) & 1; //encodes side from key (only way to get it)
	    int chan = key & 0x3f; //get MBD channel from key

	    mbenrgy[sectormb] = tower->get_energy(); //fill fields
	    mbdtype[sectormb] = type;
	    mbdside[sectormb] = side;
	    mbdchan[sectormb] = chan;
	    sectormb++;
	  }
      }
    else
      {
	cout << "no MBD towers!!!" << endl;
      }
  } 
  
  if(_dataormc) //get collision parameteres for MC
    {
      if(_debug) cout << "Getting event header info" << endl;
      EventHeaderv1 *event_header = findNode::getClass<EventHeaderv1>(topNode, "EventHeader" );
      if ( event_header )
	{
	  npart = event_header->get_intval("npart");
	  ncoll = event_header->get_intval("ncoll");
	  bimp = event_header->get_floatval("bimp");
	}
      PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
      if (!truthinfo && _debug) std::cout << PHWHERE << "PHG4TruthInfoContainer node is missing, can't collect additional truth particles"<< std::endl;
  
      PHG4TruthInfoContainer::VtxRange vtxrange = truthinfo->GetVtxRange();
      truth_vertices = 0;
      for (PHG4TruthInfoContainer::ConstVtxIterator iter = vtxrange.first; iter != vtxrange.second; ++iter)
	{
	  PHG4VtxPoint *vtx = iter->second;
	  truth_vtx[truth_vertices][0] = vtx->get_x();
	  truth_vtx[truth_vertices][1] = vtx->get_y();
	  truth_vtx[truth_vertices][2] = vtx->get_z();
	  vtx_id[truth_vertices] = vtx->get_id();
	  truth_vertices++;
	  if(_debug && truth_vertices > 998)
	    {
	      cout << "Too many truth vertices!" << endl;
	      break;
	    }
	}
    }
  
  if(_debug) cout << "Filling" << endl;

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
