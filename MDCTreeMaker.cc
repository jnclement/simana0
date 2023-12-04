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
#include <calobase/TowerInfoContainerv2.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertex.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <mbd/MbdPmtContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <ffaobjects/EventHeaderv1.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <HepMC/GenEvent.h>
#include <mbd/MbdPmtHit.h>
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
MDCTreeMaker::MDCTreeMaker(const std::string &name, const int dataormc, const int debug, const int correct):
  SubsysReco((name+"test").c_str())
{
  _evtct = 0;
  _correct = correct;
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
  
  _tree->Branch("track_vtx",track_vtx,"track_vtx[3]/F"); //svtx and mbd vtx
  _tree->Branch("mbd_vtx",mbd_vtx,"mbd_vtx[3]/F"); //mbd vtx
  _tree->Branch("sectorem",&sectorem,"sectorem/I"); //Number of hit sectors in the emcal
  _tree->Branch("sectorih",&sectorih,"sectorih/I"); // IHcal etc.
  _tree->Branch("sectoroh",&sectoroh,"sectoroh/I");
  _tree->Branch("emcalen",emcalen,"emcalen[sectorem]/F"); //energy per EMCal sector
  _tree->Branch("ihcalen",ihcalen,"ihcalen[sectorih]/F"); // per IHCal sector (etc.)
  _tree->Branch("ohcalen",ohcalen,"ohcalen[sectoroh]/F");
  _tree->Branch("emcaletabin",emcaletabin,"emcaletabin[sectorem]/I"); //eta of EMCal sector
  _tree->Branch("ihcaletabin",ihcaletabin,"ihcaletabin[sectorih]/I");
  _tree->Branch("ohcaletabin",ohcaletabin,"ohcaletabin[sectoroh]/I");
  _tree->Branch("emcalphibin",emcalphibin,"emcalphibin[sectorem]/I"); //phi of EMCal sector
  _tree->Branch("ihcalphibin",ihcalphibin,"ihcalphibin[sectorih]/I");
  _tree->Branch("ohcalphibin",ohcalphibin,"ohcalphibin[sectoroh]/I");
  _tree->Branch("sectormb",&sectormb,"sectormb/I");
  _tree->Branch("mbenrgy",mbenrgy,"mbenrgy[sectormb]/F"); //MBD reported value (could be charge or time)

  if(!_dataormc)
    {
      _tree->Branch("emcalt",emcalt,"emcalt[sectorem]/F"); //time value of EMCal sector
      _tree->Branch("ihcalt",ihcalt,"ihcalt[sectorih]/F");
      _tree->Branch("ohcalt",ohcalt,"ohcalt[sectoroh]/F");
      
      _tree->Branch("emhot",emhot,"emhot[sectorem]/I"); //time value of EMCal sector
      _tree->Branch("ihhot",ihhot,"ihhot[sectorih]/I");
      _tree->Branch("ohhot",ohhot,"ohhot[sectoroh]/I");

      _tree->Branch("emcaladc",emcaladc,"emcaladc[sectorem]/I"); //time value of EMCal sector
      _tree->Branch("ihcaladc",ihcaladc,"ihcaladc[sectorih]/I");
      _tree->Branch("ohcaladc",ohcaladc,"ohcaladc[sectoroh]/I");
      _tree->Branch("mbdtype",mbdtype,"mbdtype[sectormb]/I"); //MBD type (charge or time)
      _tree->Branch("mbdside",mbdside,"mbdside[sectormb]/I"); //MBD side (N/S)
      _tree->Branch("mbdchan",mbdchan,"mbdchan[sectormb]/I"); //MBD channel number (0-63 on each side)
      _tree->Branch("emchi2",emchi2,"emchi2[sectorem]/F");
      _tree->Branch("ihchi2",ihchi2,"ihchi2[sectorih]/F");
      _tree->Branch("ohchi2",ohchi2,"ohchi2[sectoroh]/F");
    }
  _tree->Branch("ihcalpos",ihcalpos,"ihcalpos[sectorih][3]/F"); //position (xyz) of EMCal sector center
  _tree->Branch("emcalpos",emcalpos,"emcalpos[sectorem][3]/F");
  _tree->Branch("ohcalpos",ohcalpos,"ohcalpos[sectoroh][3]/F");
  _tree->Branch("emetacor",emetacor,"emetacor[sectorem]/F"); //corrected eta value **NOT A BIN INDEX**
  _tree->Branch("ihetacor",ihetacor,"ihetacor[sectorih]/F");
  _tree->Branch("ohetacor",ohetacor,"ohetacor[sectoroh]/F");

  if(_dataormc) //get collision parameters for MC
    {
      _tree->Branch("truthpar_n",&truthpar_n,"truthpar_n/I");
      _tree->Branch("truthpar_pz",truthpar_pz,"truthpar_pz[truthpar_n]/F");
      _tree->Branch("truthpar_pt",truthpar_pt,"truthpar_pt[truthpar_n]/F");
      _tree->Branch("truthpar_e",truthpar_e,"truthpar_e[truthpar_n]/F");
      _tree->Branch("truthpar_eta",truthpar_eta,"truthpar_eta[truthpar_n]/F");
      _tree->Branch("truthpar_phi",truthpar_phi,"truthpar_phi[truthpar_n]/F");
      _tree->Branch("npart",&npart,"npart/I");
      _tree->Branch("ncoll",&ncoll,"ncoll/I");
      _tree->Branch("bimp",&bimp,"bimp/F");
      _tree->Branch("truth_vtx",truth_vtx,"truth_vtx[3]/F");
      _tree->Branch("truthpar_nh",&truthpar_nh,"truthpar_nh/I");
      _tree->Branch("truthparh_pz",truthparh_pz,"truthparh_pz[truthpar_nh]/F");
      _tree->Branch("truthparh_pt",truthparh_pt,"truthparh_pt[truthpar_nh]/F");
      _tree->Branch("truthparh_e",truthparh_e,"truthparh_e[truthpar_nh]/F");
      _tree->Branch("truthparh_eta",truthparh_eta,"truthparh_eta[truthpar_nh]/F");
      _tree->Branch("truthparh_phi",truthparh_phi,"truthparh_phi[truthpar_nh]/F");
      _tree->Branch("truthparh_id",truthparh_id,"truthparh_id[truthpar_nh]/I");
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
  _evtct++;
  if(_debug) cout << "Beginning event processing" << endl;
  if(_debug) cout << "Event " << _evtct << endl;
  float truthpar_sume = 0;
  float mbdq = 0;
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
    TowerInfoContainer *towersEM = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_CEMC");
    TowerInfoContainer *towersIH = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALIN");
    TowerInfoContainer *towersOH = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALOUT");
    TowerInfoContainer *towersEMuc = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_CEMC");
    TowerInfoContainer *towersIHuc = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_HCALIN");
    TowerInfoContainer *towersOHuc = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_HCALOUT");
    //Geometry objects for determining eta and phi locations
    RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");


    if(_dataormc)
      {
	towersEM = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_CEMC");
	towersIH = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALIN");
	towersOH = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALOUT");
	towersEMuc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_CEMC");
	towersIHuc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_HCALIN");
	towersOHuc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_HCALOUT");
      }
    //TowerInfoContainerv1 *towersMB = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_MBD");
    
    //GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode,"GlobalVertexMap");
    //GlobalVertexMap* danvtx = findNode::getClass<GlobalVertexMap>(topNode,"DansSpecialVertexMap");
    if(_debug) cout << "Checking that all necessary objects exist" << endl;
    if(!towersEM || !towersIH || !towersOH || !geomEM || !geomIH || !geomOH)// || !vertexmap)
      {
	cout << "em/ih/oh/gem/gih/goh/vtx: " << towersEM << " " << towersIH << " " << towersOH << " " << geomEM << " " << geomIH << " " << geomOH << endl; //" " << vertexmap << endl;
	return Fun4AllReturnCodes::EVENT_OK; //remove events which do not have all required information
      }
    if(!_dataormc && (!towersIHuc || !towersEMuc || !towersOHuc))// || !towersMB))// || !danvtx))
      {
	cout << "uce/uci/uco/mb/dvtx: " << " " << towersEMuc << " " << towersIHuc << " " << towersOHuc << endl; //" " << towersMB << endl;//" " << danvtx << endl;
	return Fun4AllReturnCodes::EVENT_OK;
      }
    if(_debug) cout << "EM geomtry node: " << geomEM << endl;
    if(_debug) cout << "Getting vertex" << endl;
    /*
    int mapnum = 0;
    //auto j=(_dataormc?vertexmap:danvtx)->begin();
    auto j = vertexmap->begin();
    //if(j!=(_dataormc?vertexmap:danvtx)->end())
    if(j!=vertexmap->end())
      {
	int vtxnum = 0;
	if(_debug) cout << "j iterator exists for vertexmap" << endl;
	GlobalVertex *vtx = j->second;
	if(_debug) cout << "Got j iterator.second for vertexmap" << endl;
	if(vtx)
	  {
	    track_vtx[0] = vtx->get_x();
	    track_vtx[1] = vtx->get_y();
	    track_vtx[2] = vtx->get_z();
	    if (track_vtx[0] == 0 && track_vtx[1] == 0 && track_vtx[2] == 0) 
	      { //remove anything with identically zero vertex (errors)
		if(_debug) cout << "zvtx: " << track_vtx[2] << " Vertex identically 0, skipping event" << endl;
		return Fun4AllReturnCodes::EVENT_OK;
	      }
	    if (fabs(track_vtx[2]) > 50.0) 
	      { //cut on 50 - may widen in the future
		if(_debug) cout << "Vertex outside of 50 cm, skipping event" << endl;
		return Fun4AllReturnCodes::EVENT_OK;
	      }
	  }
	else
	  {
	    if(_debug) cout << "NO VERTEX FOUND IN GLOBALVERTEXMAP AT ALL! CANCEL EVENT" << endl;
	    return Fun4AllReturnCodes::EVENT_OK;
	  }
	if(_debug) cout << "Mapnum/vtx/z: " << mapnum << " " << vtx << " " << vtx->get_z() << endl;
	for(auto i= vtx->begin_vtxids(); i!=vtx->end_vtxids(); ++i)
	  {
	    if(_debug) cout << "mapnum/vtxnum/type/id/x/y/z: "<< mapnum << " " <<vtxnum << " " << i->first << " " << i->second << endl;
	    vtxnum++;
	  }
	mapnum++;
      }
    else
      {
	if(_debug) cout << "no vertexmap, skipping event" << endl;
	return Fun4AllReturnCodes::EVENT_OK;
      }
    */
    /*
    auto iter = vertexmap->begin(); //z vertex getting
    if(iter != vertexmap->end()) 
      {
	//GlobalVertex *vtx = iter->second;
	GlobalVertex *vtx = vertexmap->get(GlobalVertex::MBD);
	if(!vtx)
	  {
	    if(_debug) cout << "no MBD VTX!" << endl;
	  }
	if(_debug) vertexmap->identify();
	if(_debug && vtx) cout << "zvtx: " << vtx->get_z() << endl;
	if(vtx)
	  {
	    mbd_vtx[0] = vtx->get_x();
	    mbd_vtx[1] = vtx->get_y();
	    mbd_vtx[2] = vtx->get_z();
	  } 
      }
    */
    //if(_dataormc && 1)
      {
	MbdVertexMap* mbdmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
	if(_debug) cout << "mbdmap: " << mbdmap << endl;
	if(!mbdmap || mbdmap->empty())
	  {
	    if(_debug) cout << "no MBD map!!" << endl;
	    return Fun4AllReturnCodes::EVENT_OK;
	  }
	auto it = mbdmap->begin();
	if(it == mbdmap->end())
	  {
	    if(_debug) cout << "Empty mbdmap!" << endl;
	    return Fun4AllReturnCodes::EVENT_OK;
	  }
	if(_debug) cout << "Made iterator for mbdmap" << endl;
	MbdVertex* mbdvtx = (*it).second;
	if(_debug) cout << "Got iterator.second from mbdmap" << endl;
	if(!mbdvtx)
	  {
	    if(_debug) cout << "no MBD vtx from MBDreco module!!" << endl;
	    return Fun4AllReturnCodes::EVENT_OK;
	  }
	if(_debug) cout << "about to get mbdvtx z value for mbdreco module with vertex pointer: " << mbdvtx << " type: " << typeid(mbdvtx).name() << endl;
	track_vtx[2] = mbdvtx->get_z() - (_dataormc?0:0);//18.3227);
	if(_debug) cout << "got mbdvtx z value for mbdreco module" << endl;
	if(track_vtx[2] == 0 || abs(track_vtx[2]) > 50)
	  {
	    if(_debug) cout << "Zero or very large MBD vtx in MBDreco module - skipping event." << endl;
	    return Fun4AllReturnCodes::EVENT_OK;
	  }
	track_vtx[0] = 0;
	track_vtx[1] = 0;
      }
      if(_debug)
	{
	  cout << "MBD zvtx: " << track_vtx[2] << endl;
	}
	/*
	  else
      {
	return Fun4AllReturnCodes::EVENT_OK;
      }
    */
      float EMCalEtot = 0;
    if(_debug) cout << "Getting EMCal info" << endl;
    if(towersEM) //get EMCal values
      {
	int nchannels = 24576; //channels in emcal
	for(int i=0; i<nchannels; ++i) //loop over channels
	  {
	    TowerInfo *tower = towersEM->get_tower_at_channel(i); //get EMCal tower
	    /*
	    if(_debug)
	      {
		TowerInfo* rawtower = towersEMuc->get_tower_at_channel(i);
		int rawkey = towersEMuc->encode_key(i);
		int raweta = towersEMuc->getTowerEtaBin(rawkey);
		int realeta = towersEM->getTowerEtaBin(towersEM->encode_key(i));
		if(raweta == 95 || raweta == 8)
		  {
		    cout << "raw/calib eta and raw/calib E: " << raweta << " " << realeta << " " << rawtower->get_energy() << " " << tower->get_energy() << endl;
		  }
	      }
	    */
	    if(!_dataormc)
	      {
		if(tower->get_chi2() > 9000*tower->get_energy()+43000) continue;
	      }
	    //if(_debug) cout << "Tower " << i << ": " << tower << endl;
 	    int key = towersEM->encode_key(i);
	    float time = towersEM->get_tower_at_channel(i)->get_time_float(); //get time
	    if(!_dataormc)
	      {
		if(towersEMuc->get_tower_at_channel(i)->get_energy() == 0)
		  {
		    if(_debug > 1)
		      {
			cout << "EMCal ADC 0 in tower " << i << endl;
		      }
		    continue;
		  }
	      }
	    if(!_dataormc && (time > 1.4 || time < -0.6))
	      {
		if(_debug > 2) cout << time << endl;
		continue; //timing cut
	      }
	    emhot[sectorem] = tower->get_isHot();
	    if(_debug > 2) cout << "i/time: " << i << " " << time << endl;
	    int etabin = towersEM->getTowerEtaBin(key); //get eta and phi indices
	    int phibin = towersEM->getTowerPhiBin(key);
	    if(tower->get_energy() < -0.1 && _debug) cout << "negative energy found!: i/E "<< i << " " << tower->get_energy() << endl;
	    //if(tower->get_energy() < -0.3 && _debug) continue;
	    emcalen[sectorem] = tower->get_energy(); //actual tower energy (calibrated)
	    emchi2[sectorem] = tower->get_chi2();
	    EMCalEtot += emcalen[sectorem];
	    const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, etabin, phibin);
	    RawTowerGeom *tower_geom = geomEM->get_tower_geometry(geomkey); //encode tower geometry
	    //if(_debug) cout << "Tower geom: " << tower_geom << endl;
	    emcalt[sectorem] = time; //store time value
	    if(!_dataormc) emcaladc[sectorem] = towersEMuc->get_tower_at_channel(i)->get_energy(); //emcal ADC value (uncalibrated "energy")
	    emcalpos[sectorem][0] = tower_geom->get_center_x(); //get positions of towers
	    emcalpos[sectorem][1] = tower_geom->get_center_y();
	    emcalpos[sectorem][2] = tower_geom->get_center_z();
	    double newz = emcalpos[sectorem][2] - track_vtx[2];
	    double newx = emcalpos[sectorem][0] - track_vtx[0];
	    double newy = emcalpos[sectorem][1] - track_vtx[1];
	    double oldx = tower_geom->get_center_x();
	    double oldy = tower_geom->get_center_y();
	    double oldz = tower_geom->get_center_z();
	    double neweta = (_correct?asinh(newz/sqrt(newx*newx+newy*newy)):asinh(oldz/sqrt(oldx*oldx+oldy*oldy)));
	    emetacor[sectorem] = neweta;
	    emcaletabin[sectorem] = etabin; //get eta and phi of towers
	    emcalphibin[sectorem] = phibin;
	    //if(_debug && abs(towersEMuc->get_tower_at_channel(i)->get_energy()) < 2) cout << etabin << " " << phibin << " " << towersEMuc->get_tower_at_channel(i)->get_energy() << " " << tower->get_isBadChi2() << " " << tower->get_isHot() << " " << tower->get_isBadTime() << " " << tower->get_isGood() << endl;
	    //if(_debug && etabin==0) cout << etabin << " " << phibin << endl;
	    sectorem++;
	  }
	if(_debug) cout << sectorem << endl;
      }
    if(_debug) cout << "total EMCal E: " << EMCalEtot << endl;
    if(_debug) cout << "Getting OHCal info" << endl;
    if(towersOH) //essentially the same as for EMCal, just with fewer sectors and named OH
      {
	int nchannels = 1536;
	for(int i=0; i<nchannels; ++i)
	  {
	    TowerInfo *tower = towersOH->get_tower_at_channel(i);
	    if(!_dataormc)
	      {
		if(tower->get_chi2() > 300*tower->get_energy()+250) continue;
	      }
	    float time = towersOH->get_tower_at_channel(i)->get_time_float();
	    if(!_dataormc){if(towersOHuc->get_tower_at_channel(i)->get_energy() == 0)
	      {
		if(_debug > 1)
                  {
                    cout << "OHCal ADC 0 in tower " << i << endl;
                  }
		continue;
	      }}
	    if(!_dataormc && (time > 2 || time < -2)) continue;
	    ohhot[sectoroh] = tower->get_isHot();
	    int key = towersOH->encode_key(i);
	    int etabin = towersOH->getTowerEtaBin(key);
	    int phibin = towersOH->getTowerPhiBin(key);
	    ohcalen[sectoroh] = tower->get_energy();
	    ohchi2[sectoroh] = tower->get_chi2();
	    const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, etabin, phibin);
	    RawTowerGeom *tower_geom = geomOH->get_tower_geometry(geomkey);
	    ohcalt[sectoroh] = time;
	    if(!_dataormc) ohcaladc[sectoroh] = towersOHuc->get_tower_at_channel(i)->get_energy();
	    ohcalpos[sectoroh][0] = tower_geom->get_center_x();
	    ohcalpos[sectoroh][1] = tower_geom->get_center_y();
	    ohcalpos[sectoroh][2] = tower_geom->get_center_z();
	    double newz = ohcalpos[sectoroh][2] - track_vtx[2];
	    double newx = ohcalpos[sectoroh][0] - track_vtx[0];
	    double newy = ohcalpos[sectoroh][1] - track_vtx[1];
	    double oldx = tower_geom->get_center_x();
	    double oldy = tower_geom->get_center_y();
	    double oldz = tower_geom->get_center_z();
	    double neweta = (_correct?asinh(newz/sqrt(newx*newx+newy*newy)):asinh(oldz/sqrt(oldx*oldx+oldy*oldy)));
	    ohetacor[sectoroh] = neweta;

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
	    TowerInfo *tower = towersIH->get_tower_at_channel(i);
	    if(!_dataormc)
	      {
		if(tower->get_chi2() > 3000) continue;
	      }
	    float time = towersIH->get_tower_at_channel(i)->get_time_float();
	    if(!_dataormc){if(towersIHuc->get_tower_at_channel(i)->get_energy() == 0)
	      {
		if(_debug > 1)
                  {
                    cout << "IHCal ADC 0 in tower " << i << endl;
                  }
		continue;
	      }}
	    if(!_dataormc && (time > 1 || time < -1)) continue;
	    ihhot[sectorih] = tower->get_isHot();
	    int key = towersIH->encode_key(i);
	    int etabin = towersIH->getTowerEtaBin(key);
	    int phibin = towersIH->getTowerPhiBin(key);
	    ihcalen[sectorih] = tower->get_energy();
	    ihchi2[sectorih] = tower->get_chi2();
	    const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, etabin, phibin);
	    RawTowerGeom *tower_geom = geomIH->get_tower_geometry(geomkey);
	    ihcalt[sectorih] = time;
	    if(!_dataormc) ihcaladc[sectorih] = towersIHuc->get_tower_at_channel(i)->get_energy();
	    ihcalpos[sectorih][0] = tower_geom->get_center_x();
	    ihcalpos[sectorih][1] = tower_geom->get_center_y();
	    ihcalpos[sectorih][2] = tower_geom->get_center_z();
	    double newz = ihcalpos[sectorih][2] - track_vtx[2];
	    double newx = ihcalpos[sectorih][0] - track_vtx[0];
	    double newy = ihcalpos[sectorih][1] - track_vtx[1];
	    double oldx = tower_geom->get_center_x();
	    double oldy = tower_geom->get_center_y();
	    double oldz = tower_geom->get_center_z();
	    double neweta = (_correct?asinh(newz/sqrt(newx*newx+newy*newy)):asinh(oldz/sqrt(oldx*oldx+oldy*oldy)));
	    ihetacor[sectorih] = neweta;

	    ihcaletabin[sectorih] = etabin;
	    ihcalphibin[sectorih] = phibin;
	    sectorih++;
	  }
      }
    if(_debug) cout << "Getting MBD info" << endl;
    /*
    if(towersMB && !_dataormc) //get MBD info
      {
	int nchannels = 256;
	for(int i=0; i<nchannels; ++i)
	  {
	    TowerInfo *tower = towersMB->get_tower_at_channel(i); //get mbd channel
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
    else if(!_dataormc)
      {
	cout << "no MBD towers!!!" << endl;
      }
    */
  }
      MbdPmtContainer *mbdtow = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
      if(!mbdtow)
	{
	  if(_debug) cout << "No MBD PMT Container found!" << endl;
	  return Fun4AllReturnCodes::EVENT_OK;
	}
      sectormb = 128;//mbdtow->get_npmt();
      //if(_debug) cout << "Got " << sectormb << " mbd sectors in sim." << endl;
      for(int i=0; i<sectormb; ++i)
	{
	  MbdPmtHit *mbdhit = mbdtow->get_pmt(i);
	  if(_debug > 1) cout << "PMT " << i << " address: " << mbdhit << " charge: " << mbdhit->get_q() << endl;
	  mbenrgy[i] = mbdhit->get_q();
	  mbdq += mbdhit->get_q();
	}

        
  if(_dataormc) //get collision parameters for MC
    {

      PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
      if (!hepmceventmap)
	{
	  if(_debug) std::cout << PHWHERE << "HEPMC event map node is missing, can't collected HEPMC truth particles"<< std::endl;
	  return Fun4AllReturnCodes::EVENT_OK;
	}
      for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin(); eventIter != hepmceventmap->end(); ++eventIter) {
     
	/// Get the event
	PHHepMCGenEvent *hepmcevent = eventIter->second;
      
	// To fill TTree, require that the event be the primary event (embedding_id > 0)
	if (hepmcevent && hepmcevent->get_embedding_id() == 0)
	  {
	    /// Get the event characteristics, inherited from HepMC classes
	    HepMC::GenEvent *truthevent = hepmcevent->getEvent();
	    if (!truthevent) {
	      if(_debug) std::cout << PHWHERE
			<< "no evt pointer under phhepmvgeneventmap found "
			<< std::endl;
	      return Fun4AllReturnCodes::EVENT_OK;
	    }

	    /// Loop over all the truth particles and get their information
	    truthpar_nh = 0;
	    for (HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin(); iter != truthevent->particles_end(); ++iter) {
	      if (!(*iter)->end_vertex() && (*iter)->status() == 1) {
		truthparh_e[truthpar_nh] = (*iter)->momentum().e();
		double px = (*iter)->momentum().px();
		double py = (*iter)->momentum().py();
		double pz = (*iter)->momentum().pz();
		truthparh_pt[truthpar_nh] = sqrt(px*px+py*py);
		truthparh_pz[truthpar_nh] = pz;
		truthparh_phi[truthpar_nh]=atan2(py,px);
		truthparh_eta[truthpar_nh]=atanh(pz/sqrt(px*px+py*py+pz*pz));
		truthparh_id[truthpar_nh]=(*iter)->pdg_id();
		/// Fill the truth tree
		truthpar_nh++;
	      }
	    }
	    break;
	  } 
      }


      PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
      if (!truthinfo && _debug) std::cout << PHWHERE << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"<< std::endl;
      if (!truthinfo) return Fun4AllReturnCodes::EVENT_OK;
      PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
      truthpar_n = 0;
      for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
	{
	  float pmass = 0.938272;
	  // Get truth particle
	  const PHG4Particle *truth = iter->second;
	  if (!truthinfo->is_primary(truth)) continue;
	  
	  /// Get this particles momentum, etc.
	  truthpar_pt[truthpar_n] = sqrt(truth->get_px() * truth->get_px()
					 + truth->get_py() * truth->get_py());
	  truthpar_pz[truthpar_n] = truth->get_pz();
	  float psquare = truthpar_pt[truthpar_n]*truthpar_pt[truthpar_n]+truthpar_pz[truthpar_n]*truthpar_pz[truthpar_n];
	  std::vector<int>::iterator parit = find(baryons.begin(), baryons.end(), truth->get_pid());
	  std::vector<int>::iterator apart = find(baryons.begin(), baryons.end(), -truth->get_pid());
	  if(parit != baryons.end())
	    {
	      truthpar_e[truthpar_n] = truth->get_e() - sqrt(truth->get_e()*truth->get_e()-psquare);
	    }
	  else if(apart != baryons.end())
	    {
	      truthpar_e[truthpar_n] = truth->get_e() - sqrt(truth->get_e()*truth->get_e()-psquare) + 2*pmass;
	    }
	  else
	    {
	      truthpar_e[truthpar_n] = truth->get_e();
	    }
	  truthpar_sume += truth->get_e();
	  truthpar_phi[truthpar_n] = atan2(truth->get_py(), truth->get_px());
	  truthpar_eta[truthpar_n] = atanh(truth->get_pz() / sqrt(truth->get_px()*truth->get_px()+truth->get_py()*truth->get_py()+truth->get_pz()*truth->get_pz()));
	  if (truthpar_eta[truthpar_n] != truthpar_eta[truthpar_n]) truthpar_eta[truthpar_n] = -999; // check for nans
	  truthpar_n++;
	  if(truthpar_n > 99999)
	    {
	      if(_debug) cout << "More than 100000 truth particles!" << endl;
	      return Fun4AllReturnCodes::EVENT_OK;
	    }
	}
      if(_debug) cout << "Getting event header info" << endl;
      EventHeaderv1 *event_header = findNode::getClass<EventHeaderv1>(topNode, "EventHeader" );
      if ( event_header )
	{
	  npart = event_header->get_intval("npart");
	  ncoll = event_header->get_intval("ncoll");
	  bimp = event_header->get_floatval("bimp");
	}
      else return Fun4AllReturnCodes::EVENT_OK;
      if (!truthinfo && _debug) std::cout << PHWHERE << "PHG4TruthInfoContainer node is missing, can't collect additional truth particles"<< std::endl;
      else if(!truthinfo) return Fun4AllReturnCodes::EVENT_OK;

      PHHepMCGenEventMap *phg = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
      if(phg)
	{
	  PHHepMCGenEvent *phe = (PHHepMCGenEvent*) phg->get(0);
	  if(phe)
	    {
	      truth_vtx[0] = phe->get_collision_vertex().x();
	      truth_vtx[1] = phe->get_collision_vertex().y();
	      truth_vtx[2] = phe->get_collision_vertex().z();
	    }
	  else
	    {
	       if(_debug) cout << "no PHHepMCGenEvent found for truth vertex getting" << endl;
	       return Fun4AllReturnCodes::EVENT_OK;
	    }
	}
      else
	{
	  if(_debug) cout << "no PHHepMCGenEventMap found for truth vertex getting" << endl;
	  return Fun4AllReturnCodes::EVENT_OK;
	}
    }
  
  if(_debug) cout << "Filling" << endl;

  _tree->Fill();
  if(_debug)
    {
      cout << "Total truth E |eta|<1: " << truthpar_sume << endl;
      cout << "Total mbd q: " << mbdq << endl;
      cout << "npart: " << npart << endl;
    }
  
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
