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

//#include <g4jets/JetMap.h>
//#include <g4jets/Jet.h>

#include <TLorentzVector.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <iostream>

#include <centrality/CentralityInfov1.h>

using namespace std;
int eventchecker = 0;
//sorry about the global - just for debug.
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


  hotfile = TFile::Open("/sphenix/user/jocl/projects/sandbox/run/hot_towers_21518_1.root");
  hottree = hotfile->Get<TTree>("T_hot_tower");
    
  _f = new TFile( _foutname.c_str(), "RECREATE");

  //std::cout << " making a file = " <<  _foutname.c_str() << " , _f = " << _f << std::endl;
  
  _tree = new TTree("ttree","a persevering date tree");
  //_tree->Branch("emfemevt",emfemevt,"emfemevt[48]/I");
  //_tree->Branch("empktfem",empktfem,"empktfem[128]/I");
  //_tree->Branch("emfemclk",emfemclk,"emfemclk[48]/I");
  //_tree->Branch("emxmtevt",emxmtevt,"emxmtevt[128]/I");
  //_tree->Branch("emxmtclk",emxmtclk,"emxmtclk[128]/I");
  //_tree->Branch("empktevt",empktevt,"empktevt[128]/I");
  //_tree->Branch("hofemevt",hofemevt,"hofemevt[48]/I");
  //_tree->Branch("hopktfem",hopktfem,"hopktfem[32]/I");
  //_tree->Branch("hofemclk",hofemclk,"hofemclk[12]/I");
  //_tree->Branch("hoxmtevt",hoxmtevt,"hoxmtevt[32]/I");
  //_tree->Branch("hoxmtclk",hoxmtclk,"hoxmtclk[32]/I");
  //_tree->Branch("hopktevt",hopktevt,"hopktevt[32]/I");
  //_tree->Branch("hifemevt",hifemevt,"hifemevt[48]/I");
  //_tree->Branch("hipktfem",hipktfem,"hipktfem[32]/I");
  //_tree->Branch("hifemclk",hifemclk,"hifemclk[12]/I");
  //_tree->Branch("hixmtevt",hixmtevt,"hixmtevt[32]/I");
  //_tree->Branch("hixmtclk",hixmtclk,"hixmtclk[32]/I");
  //_tree->Branch("hipktevt",hipktevt,"hipktevt[32]/I");
  //_tree->Branch("mbfemevt",mbfemevt,"mbfemevt[4]/I");
  //_tree->Branch("mbpktfem",mbpktfem,"mbpktfem[2]/I");
  //_tree->Branch("mbfemclk",mbfemclk,"mbfemclk[4]/I");
  //_tree->Branch("mbxmtevt",mbxmtevt,"mbxmtevt[2]/I");
  //_tree->Branch("mbxmtclk",mbxmtclk,"mbxmtclk[2]/I");
  //_tree->Branch("mbpktevt",mbpktevt,"mbpktevt[2]/I");
  /*
  _tree->Branch("emsize",&emsize,"emsize/I");
  _tree->Branch("hisize",&hisize,"hisize/I");
  _tree->Branch("hosize",&hosize,"hosize/I");
  _tree->Branch("emkey",&emkey);
  _tree->Branch("hikey",&hikey);
  _tree->Branch("hokey",&hokey);
  _tree->Branch("emwf",&emwf);
  _tree->Branch("hiwf",&hiwf);
  _tree->Branch("howf",&howf);
  */
  /*
  _tree->Branch("truthjet_n",&truthjet_n,"truthjet_n/I");
  _tree->Branch("truthjet_pt",truthjet_pt,"truthjet_pt[truthjet_n]/F");
  _tree->Branch("truthjet_et",truthjet_et,"truthjet_et[truthjet_n]/F");
  _tree->Branch("truthjet_ph",truthjet_ph,"truthjet_ph[truthjet_n]/F");
  _tree->Branch("reco_jet_n",&reco_jet_n,"reco_jet_n/I");
  _tree->Branch("reco_jet_pt",reco_jet_pt,"reco_jet_pt[reco_jet_n]/F");
  _tree->Branch("reco_jet_et",reco_jet_et,"reco_jet_et[reco_jet_n]/F");
  _tree->Branch("reco_jet_ph",reco_jet_ph,"reco_jet_ph[reco_jet_n]/F");
  _tree->Branch("truthpar_n",&truthpar_n,"truthpar_n/I");
  _tree->Branch("truthpar_pt",truthpar_pt,"truthpar_pt[truthpar_n]/F");
  _tree->Branch("truthpar_et",truthpar_et,"truthpar_et[truthpar_n]/F");
  _tree->Branch("truthpar_ph",truthpar_ph,"truthpar_ph[truthpar_n]/F");
  _tree->Branch("truthpar_id",truthpar_id,"truthpar_id[truthpar_n]/I");
  */
  // _tree->Branch("truthpar_n",&truthpar_n1,"truthpar_n1/I");
  // _tree->Branch("truthpar_pt",truthpar_pt1,"truthpar_pt1[truthpar_n1]/F");
  // _tree->Branch("truthpar_et",truthpar_et1,"truthpar_et1[truthpar_n1]/F");
  // _tree->Branch("truthpar_ph",truthpar_ph1,"truthpar_ph1[truthpar_n1]/F");
  // _tree->Branch("truthpar_id",truthpar_id1,"truthpar_id1[truthpar_n1]/I");
  /*
  _tree->Branch("truthpar_j",truthpar_j,"truthpar_j[truthpar_n]/I");
  _tree->Branch("emfrac",emfrac,"emfrac[truthjet_n]/F");
  _tree->Branch("truthpar_em",truthpar_em,"truthpar_em[truthpar_n]/I");
  _tree->Branch("truthpar_em1",truthpar_em1,"truthpar_em1[truthpar_n1]/I");
  */
   _tree->Branch("zvtx",&zvtx,"zvtx/F");
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
  // _tree->Branch("nclus",&nclus,"nclus/I");
  // _tree->Branch("prob",prob,"prob[nclus]/F");
  // _tree->Branch("chi2",chi2,"chi2[nclus]/F");
  // _tree->Branch("cluster_E",_cluster_E,"cluster_E[nclus]/F");
  // _tree->Branch("cluster_eta",_cluster_eta,"cluster_eta[nclus]/F");
  // _tree->Branch("cluster_phi",_cluster_phi,"cluster_phi[nclus]/F");
  // _tree->Branch("cluster_ntower",_cluster_ntower,"cluster_ntower[nclus]/I");
  // _tree->Branch("cluster_Ecore",cluster_Ecore,"cluster_Ecore[nclus]/F");
  //_tree->Branch("clustowE",clustowE,"clustowE[cluster_ntower[nclus]]/F");
  //_tree->Branch("clustowet",clustowet,"clustowet[cluster_ntower[nclus]]/F");
  //_tree->Branch("clustowph",clustowph,"clustowph[cluster_ntower[nclus]]/F");
  // _tree->Branch("npart",&npart,"npart/I");
  // _tree->Branch("ncoll",&ncoll,"ncoll/I");
  // _tree->Branch("bimp",&bimp,"bimp/F");
  // _tree->Branch("bestclus",&bestclus,"bestclus/I");
  // _tree->Branch("bcnt",&bcnt,"bcnt/I");
  // _tree->Branch("bcmt",&bcmt,"bcmt/I");
  // _tree->Branch("bctet",bctet,"bctet[bcnt]/F");
  // _tree->Branch("bctph",bctph,"bctph[bcnt]/F");
  // _tree->Branch("bcten",bcten,"bcten[bcnt]/F");
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MDCTreeMaker::InitRun(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
/*
void climbup(int depth, HepMC::GenVertex::particle_iterator part)
{
  ++depth;
  if(!(*part)->production_vertex())
    {
      for(int i=0; i<depth; i++)
	{
	  cout << "   ";
	}
      cout << "pid=" << (*part)->pdg_id() << "barcd=" << (*part)->barcode() << " top level particle" << endl;
    }
  else
    {
      for(HepMC::GenVertex::particle_iterator newpart = (*part)->production_vertex()->particles_begin(HepMC::parents);
	  newpart != (*part)->production_vertex()->particles_end(HepMC::parents); ++newpart)
	{
	  for(int i=0; i<depth; i++)
	    {
	      //cout << "   ";
	    }
	  //cout << "pid=" << (*newpart)->pdg_id() << " barcd=" << (*newpart)->barcode() <<" production position x y z = "
	  //     << (*part)->production_vertex()->position().x() << " " << (*part)->production_vertex()->position().y() <<
	  //  " " << (*part)->production_vertex()->position().z() << endl;
	  //climbup(depth, newpart);
	}
    }
}
*/
//____________________________________________________________________________..
int MDCTreeMaker::process_event(PHCompositeNode *topNode)
{
  int hotet;
  int hotph;
  int hotin;
  hottree->SetBranchAddress("hot_eta",&hotet);
  hottree->SetBranchAddress("hot_phi",&hotph);
  hottree->SetBranchAddress("index",&hotin);
  int hoteta[1000];
  int hotpha[1000];
  int hotina[1000];
  int hots = hottree->GetEntries();
  for(int i=0; i<hots; ++i)
    {
      hottree->GetEntry(i);
      hoteta[i] = hotet;
      hotpha[i] = hotph;
      hotina[i] = hotin;
    }
  /*
  emsize = 0;
  hisize = 0;
  hosize = 0;
  emkey.clear();
  hikey.clear();
  hokey.clear();
  emwf.clear();
  hiwf.clear();
  howf.clear();
  */
  /*
  for(int i=0; i<48; ++i)
     {
  //for(int j=0; j<3; ++j)
  //	{
       //emfemevt[i] = 0;
	  emfemclk[i] = -1;
	  if(i<12)
	    {
	  //hofemevt = 0;
	  hofemclk[i] = -1;
	  //hifemevt = 0;
	  hifemclk[i] = -1;
	    }
	  if(i<4)
	    {
	  //mbfemevt = 0;
	  mbfemclk[i] = -1;
	    }
	  //	}
	  //femit[i] = 0;
	  //evtit[i] = 0;
      }
  */
  // here is where the code goes 
//  truthjet_n = 0;
//  reco_jet_n = 0;
//  truthpar_n = 0;
//   truthpar_n1 = 0;
   sectorem = 0;
   sectorih = 0;
   sectoroh = 0;
   sectormb = 0;
//   nclus = 0;
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
//   //    truthpar_px[i] = 0;
//   //    truthpar_py[i] = 0;
//   //    truthpar_pz[i] = 0;
//   //    truthpar_pt[i] = 0;
//   //    truthpar_et[i] = 0;
//   //    truthpar_ph[i] = 0;
// //      truthpar_id[i] = 0;
//       truthpar_pt1[i] = 0;
//       truthpar_et1[i] = 0;
//       truthpar_ph1[i] = 0;
//       truthpar_id1[i] = 0;
//   //    truthpar_j[i] = 0;
//   //    truthpar_em[i] = 0;
//   //    truthpar_em1[i] = 0;
//       prob[i] = 0;
//       chi2[i] = 0;
     }
/*  
  for(int i=0; i<100; i++)
    {
      emfrac[i] = 0;
    }
  */
   /*
  WaveformContainerv1* emcwf = findNode::getClass<WaveformContainerv1>(topNode,"WAVEFORMS_CEMC");
  if(emcwf)
    {
      int i=0;
      // WaveformContainerv1::RangeInfo femevt = emcwf->get_fem_events();
      //WaveformContainerv1::IterInfo ievt = femevt.first;
      WaveformContainerv1::RangeInfo femclk = emcwf->get_fem_clocks();
      WaveformContainerv1::IterInfo iclk = femclk.first;
      /*
      for(; ievt != femevt.second; ++ievt)
	{
	  //int packet = ((ievt->first >> 16)) - 6001;
	  emfemevt[i][j] = ievt->second;
	  //evtit[]++;
	}
      *
      for(; iclk != femclk.second; ++iclk)
	{
	  //int packet = ((iclk->first >> 16)) -6001;
	  emfemclk[i] = iclk->second;
	  //femit[packet]++;
	  ++i;
	  if(i>48) cout << "i > 48 in emcal" << endl;
	}
    }
  else 
    {
      cout << "NO CEMC WF" << endl;
      exit(1);
    }
  /*
  for(int i=0; i<128; ++i)
    {
      femit[i] = 0;
      evtit[i] = 0;
    }
  */
  /*
    WaveformContainerv1* hicwf = findNode::getClass<WaveformContainerv1>(topNode,"WAVEFORMS_HCALIN");
  if(hicwf)
    {
      int i=0;
      //WaveformContainerv1::RangeInfo femevt = hicwf->get_fem_events();
      //WaveformContainerv1::IterInfo ievt = femevt.first;
      WaveformContainerv1::RangeInfo femclk = hicwf->get_fem_clocks();
      WaveformContainerv1::IterInfo iclk = femclk.first;
      //hifemevt = ievt->second;
      // hifemclk = iclk->second;
      /*
      for(; ievt != femevt.second; ++ievt)
	{
	  int packet = ((ievt->first >> 16) & 0xff) - 1;
	  
	  hifemevt[packet][evtit[packet]] = ievt->second;
	  evtit[packet]++;
	}
  *
      for(; iclk != femclk.second; ++iclk)
	{
	  //int packet = ((iclk->first >> 16) &0xff) -1;
	  
	  hifemclk[i] = iclk->second;
	  //femit[packet]++;
	  ++i;
	  if(i>12) cout << "i > 12 in hcali" << endl;
	}
    }
  else
    {
      cout << "NO HCALIN WF" << endl;
      exit(1);
    }
  /*
  for(int i=0; i<128; ++i)
    {
      femit[i] = 0;
      evtit[i] = 0;
    }
  *
    WaveformContainerv1* hocwf = findNode::getClass<WaveformContainerv1>(topNode,"WAVEFORMS_HCALOUT");
  if(hocwf)
    {
      int i=0;
      //WaveformContainerv1::RangeInfo femevt = hocwf->get_fem_events();
      //WaveformContainerv1::IterInfo ievt = femevt.first;
      WaveformContainerv1::RangeInfo femclk = hocwf->get_fem_clocks();
      WaveformContainerv1::IterInfo iclk = femclk.first;
      //hifemevt = ievt->second;
      //hifemclk = iclk->second;
      /*
      for(; ievt != femevt.second; ++ievt)
	{
	  int packet = ((ievt->first >> 16) & 0xff) - 1;
	  
	  hofemevt[packet][evtit[packet]] = ievt->second;
	  evtit[packet]++;
	}
      *
      for(; iclk != femclk.second; ++iclk)
	{
	  //int packet = ((iclk->first >> 16) &0xff) -1;
	  
	  hofemclk[i] = iclk->second;
	  ++i;
	  //femit[packet]++;
	  if(i>12) cout << "i > 12 in hcalo" <<endl;
	}
    }
  else
    {
      cout << "NO HCALOUT WF" << endl;
      exit(1);
    }
  /*
  for(int i=0; i<128; ++i)
    {
      femit[i] = 0;
      evtit[i] = 0;
    }
  *
  WaveformContainerv1* mbdwf = findNode::getClass<WaveformContainerv1>(topNode,"WAVEFORMS_MBD");
  if(mbdwf)
    {
      int i=0;
      //WaveformContainerv1::RangeInfo femevt = mbdwf->get_fem_events();
      //WaveformContainerv1::IterInfo ievt = femevt.first;
      WaveformContainerv1::RangeInfo femclk = mbdwf->get_fem_clocks();
      WaveformContainerv1::IterInfo iclk = femclk.first;
      //mbfemevt = ievt->second;
      //mbfemclk = iclk->second;
      /*
      for(; ievt != femevt.second; ++ievt)
	{
	  int packet = (((ievt->first) >> 16)) - 1001;
	  
	  mbfemevt[packet][evtit[packet]] = ievt->second;
	  evtit[packet]++;
	}
      *
      for(; iclk != femclk.second; ++iclk)
	{
	  //int packet = ((iclk->first >> 16)) -1001;
	  
	  mbfemclk[i] = iclk->second;
	  //femit[packet]++;
	  ++i;
	  if(i>4) cout << "i > 4 in mb" << endl;
	}
    }
  else
    {
      cout << "NO MBD WF" << endl;
      exit(1);
    }
  */

  {

  TowerInfoContainerv1 *towersEM = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_CEMC");
  TowerInfoContainerv1 *towersIH = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainerv1 *towersOH = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALOUT");
  TowerInfoContainerv1 *towersEMuc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_CEMC");
  TowerInfoContainerv1 *towersIHuc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_HCALIN");
  TowerInfoContainerv1 *towersOHuc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_HCALOUT");

  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  TowerInfoContainerv1 *towersMB = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_MBD");

  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode,"GlobalVertexMap");
  auto iter = vertexmap->begin();
  if(iter != vertexmap->end()) {
    GlobalVertex *vtx = iter->second;
    float x_vtx = vtx->get_x();
    float y_vtx = vtx->get_y();
    zvtx = vtx->get_z();
    if (x_vtx == 0 && y_vtx == 0 && zvtx == 0) {
      //event++;
      return Fun4AllReturnCodes::EVENT_OK;
    }
    //h_vertex->Fill(zvtx);
    if (fabs(zvtx) > 30.0) {
      //event++;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  } else {
    //event++;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  if(towersEM)
    {
      int nchannels = 24576;
      float emsume = 0;
      for(int i=0; i<nchannels; ++i)
	{
	  int skip = 0;
	  TowerInfov1 *tower = towersEM->get_tower_at_channel(i);
	  int key = towersEM->encode_key(i);
	  int time0 = towersEMuc->get_tower_at_channel(i)->get_time();
	  int time = towersEM->get_tower_at_channel(i)->get_time();
	  //if(time != time0) cout << i << " " << time << " emcal " << time0 << endl;
	  if(time > 7 || time < 5) continue;
	  int etabin = towersEM->getTowerEtaBin(key);
	  int phibin = towersEM->getTowerPhiBin(key);
	  /*
	  if(etabin == 37 && phibin == 0) continue;
	  if(etabin == 47 && phibin == 138) continue;
	  if(etabin == 80 && phibin == 228) continue;
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
	  emsume += tower->get_energy();
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
	  int time = towersOH->get_tower_at_channel(i)->get_time();
	  int time0 = towersOHuc->get_tower_at_channel(i)->get_time();
	  //if(time != time0) cout << i << " " << time << " hcal " << time0 << endl;
	  if(time > 8 || time < 5) continue;
	  int key = towersOH->encode_key(i);
	  int etabin = towersOH->getTowerEtaBin(key);
	  int phibin = towersOH->getTowerPhiBin(key);
	  /*
	  if (etabin == 23 && phibin == 2) continue;
	  if (etabin == 7 && phibin == 54) continue;
	  if (etabin == 22 && phibin == 26) continue;
	  */
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
	  int time0  = towersIHuc->get_tower_at_channel(i)->get_time();
	  int time = towersIH->get_tower_at_channel(i)->get_time();
	  //if(time != time0) cout << i << " " << time << " ihcal " << time0 << endl;
	  if(time > 7 || time < 5) continue;
	  int key = towersIH->encode_key(i);
	  int etabin = towersIH->getTowerEtaBin(key);
	  int phibin = towersIH->getTowerPhiBin(key);
	  /*
	  if (etabin == 3 && phibin == 15) continue;
	  if (etabin == 16 && phibin == 33) continue;
	  if (etabin == 19 && phibin == 1) continue;
	  if (etabin == 22 && phibin == 2) continue;
	  if (etabin == 22 && phibin == 3) continue;
	  */
	  //float eta = geomIH->get_etacenter(etabin);
	  //float phi = geomIH->get_phicenter(phibin);
	  //cout << "it: " << sector << " eta: " << eta << " eng: " << tower->get_energy() << endl;
	  ihcalen[sectorih] = tower->get_energy();
	  ihcalet[sectorih] = etabin;
	  ihcalph[sectorih] = phibin;
	  sectorih++;
	}
    }
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
  /*
  if(towersIH && geomIH)
    {
      TowerInfoContainerv1::ConstRange begin_end = towersIH->getTowers();
      for(TowerInfoContainerv1::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
	{
	  RawTower *tower = rtiter->second;

	  RawTowerGeom *tower_geom = geomIH->get_tower_geometry(tower->get_key());
	  float eta = tower_geom->get_eta();
	  ihcalen[sectorih] = tower->get_energy();
	  ihcalet[sectorih] = eta;
	  ihcalph[sectorih] = tower_geom->get_phi();
	  sectorih++;
	}
    }
  
  if(towersOH && geomOH)
    {
      TowerInfoContainerv1::ConstRange begin_end = towersOH->getTowers();
      for(TowerInfoContainerv1::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
	{
	  RawTower *tower = rtiter->second;

	  RawTowerGeom *tower_geom = geomOH->get_tower_geometry(tower->get_key());

	  ohcalen[sectoroh] = tower->get_energy();
	  ohcalet[sectoroh] = tower_geom->get_eta();
	  ohcalph[sectoroh] = tower_geom->get_phi();
	  sectoroh++;
	}
    }
  
  }
  //vector<int> pibarcodes;
  /*
  {
    
    JetMap* jets = findNode::getClass<JetMap>(topNode,"AntiKt_Truth_r04");

    PHG4TruthInfoContainer* truthinfo = findNode::getClass <PHG4TruthInfoContainer> (topNode, "G4TruthInfo");

    for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter) {

      Jet* this_jet = iter->second;
      
      float jet_pt = this_jet->get_pt();
      float jet_et = this_jet->get_eta();
      float jet_ph = this_jet->get_phi();

      truthjet_pt[ truthjet_n ] = jet_pt;
      truthjet_et[ truthjet_n ] = jet_et;
      truthjet_ph[ truthjet_n ] = jet_ph;
      //cout << "BEGIN NEW JET:" << endl;
      
      for (Jet::ConstIter comp = this_jet->begin_comp(); comp !=  this_jet->end_comp(); ++comp) {
        PHG4Particle* g4particle = truthinfo->GetParticle( (*comp).second  );
	
	if(g4particle->get_pid() == 111)
	  {
	    pibarcodes.push_back(g4particle->get_barcode());
	  }
	cout << "Jet particle: pid=" << g4particle->get_pid() << " Embedid= " << truthinfo->isEmbeded(g4particle->get_track_id())
	     << " pT=" << sqrt(g4particle->get_px()*g4particle->get_px() + g4particle->get_py()*g4particle->get_py())
	     << " trkid=" << g4particle->get_track_id() << " vtxid=" << g4particle->get_vtx_id() << " parid=" << g4particle->get_parent_id()
	     << " prmid=" << g4particle->get_primary_id() << " barcd=" << g4particle->get_barcode() << endl;
	
        TLorentzVector t;
        t.SetPxPyPzE (g4particle->get_px (), g4particle->get_py (), g4particle->get_pz (), g4particle->get_e ());
        float truth_pt = t.Pt ();
        float truth_eta = t.Eta ();
        float truth_phi = t.Phi ();
        int truth_pid = g4particle->get_pid ();

	if(truthinfo->isEmbeded(g4particle->get_track_id()) == 1)
	  {
	    truthpar_em[truthpar_n] = truthinfo->isEmbeded(g4particle->get_track_id());
	    truthpar_pt[truthpar_n] = truth_pt;
	    truthpar_et[truthpar_n] = truth_eta;
	    truthpar_ph[truthpar_n] = truth_phi;
	    truthpar_id[truthpar_n] = truth_pid;
	    truthpar_j[truthpar_n] = truthjet_n;
	    truthpar_n++;
	    if(truth_pid == 11 || truth_pid == 22 || truth_pid == 111)
	      {
		emfrac[truthjet_n] += truth_pt/jet_pt;
	      }
	  }
      }
      truthjet_n++;
    }
  }
  
  {
    JetMap* jets = findNode::getClass<JetMap>(topNode,"AntiKt_Tower_r04");
    
    for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter) {
      
      Jet* this_jet = iter->second;
      
      float jet_pt = this_jet->get_pt();
      float jet_et = this_jet->get_eta();
      float jet_ph = this_jet->get_phi();
      
      
      
      reco_jet_pt[ reco_jet_n ] = jet_pt;
      
      reco_jet_et[ reco_jet_n ] = jet_et;
      
      reco_jet_ph[ reco_jet_n ] = jet_ph;
      
      reco_jet_n++;
      
    }
  }
  */
  /*
  {    
    //PHHepMCGenEventMap *genEventMap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
    //PHHepMCGenEvent *genEvent = genEventMap->get(1);
    //HepMC::GenEvent *theEvent = genEvent->getEvent();
    PHG4TruthInfoContainer* truthinfo = findNode::getClass <PHG4TruthInfoContainer> (topNode, "G4TruthInfo");
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange ();    
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
      {
	PHG4Particle* g4particle = iter->second;
	
	if(truthinfo->isEmbeded(g4particle->get_track_id()) < 1) continue;
  */
	/*
	if(g4particle->get_pid() != 111) continue;
	cout << "Geant4 particle: pid=" << g4particle->get_pid() << " Embedid=" <<  truthinfo->isEmbeded(g4particle->get_track_id())
	     << " pT=" << sqrt(g4particle->get_px()*g4particle->get_px() + g4particle->get_py()*g4particle->get_py()) <<
	  " trkid=" << g4particle->get_track_id() << " vtxid=" << g4particle->get_vtx_id() << " parid=" << g4particle->get_parent_id()
	     << " prmid=" << g4particle->get_primary_id() << " barcd=" << g4particle->get_barcode() << " energy="
	     << g4particle->get_e() << endl;
	*/
  /*
	TLorentzVector t;
	t.SetPxPyPzE (g4particle->get_px (), g4particle->get_py (), g4particle->get_pz (), g4particle->get_e ());
	float truth_pt = t.Pt ();
	float truth_eta = t.Eta ();
	float truth_phi = t.Phi ();
	int truth_pid = g4particle->get_pid (); // particle species     
	// take only particles from primary Pythia event                                                            
//	truthpar_em1[truthpar_n1] = truthinfo->isEmbeded(g4particle->get_track_id());
	truthpar_pt1[truthpar_n1] = truth_pt;
	truthpar_et1[truthpar_n1] = truth_eta;
	truthpar_ph1[truthpar_n1] = truth_phi;
	truthpar_id1[truthpar_n1] = truth_pid;
	truthpar_n1++;
  */
	/*
	if(g4particle->get_pid() != 111) continue;
	if(truthinfo->isEmbeded(g4particle->get_track_id()) != 1) continue;
	counter++;
	for(HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin(); p!= theEvent->particles_end(); ++p)
	  {
	    assert(*p);
	    if((*p)->pdg_id() != 111) continue;
	    if(!(*p)->production_vertex()) continue;
	    if(((*p)->barcode() == g4particle->get_barcode()))
	      {
		if((*p)->production_vertex())
		  {
		    
		    //cout << "Production vertex position x: " << (*p)->production_vertex()->position().x() << endl;
		    //cout << "End vertex pointer: " << (*p)->end_vertex() << endl;
		    //cout << "Production vertex position y: " << (*p)->production_vertex()->position().y() << endl;
		    //cout << "Production vertex position z: " << (*p)->production_vertex()->position().z() << endl;
		    
		  }
		//cout << "Pythia particle: pid=" << (*p)->pdg_id() << " barcd=" << (*p)->barcode() << " pT=" << sqrt((*p)->momentum().px() * (*p)->momentum().px() + (*p)->momentum().py() * (*p)->momentum().py()) << " energy=" << (*p)->momentum().e() << endl;
		for(HepMC::GenVertex::particle_iterator mother = (*p)->production_vertex()->particles_begin(HepMC::parents); mother != (*p)->production_vertex()->particles_end(HepMC::parents); ++mother)
		  {
		    //cout << "   mother pid: " << (*mother)->pdg_id() << " barcd: " << (*mother)->barcode() << " production position x y z = " << (*mother)->production_vertex()->position().x() << " " << (*mother)->production_vertex()->position().y() << " " << (*mother)->production_vertex()->position().z() << endl;
		    climbup(1, mother);
		  }
	      }
	  }
	*/
  /*
      }
  }
  */
  /*
  {
    float maxE = 0;
    RawClusterContainer *clusterEM = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");
    TowerInfoContainerv1 *towersEM = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWER_CALIB_CEMC");
    RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if(clusterEM && geomEM && towersEM)
      {
	RawClusterContainer::Range begin_end = clusterEM->getClusters();
	for(RawClusterContainer::Iterator iter = begin_end.first; iter != begin_end.second; ++iter)
	  {
	    float theta = 3.14159 / 2.0 - atan2( iter->second->get_z() , iter->second->get_r() );
	    float eta = -1 * log( tan( theta / 2.0 ) );
	    RawCluster *clust = iter->second;
	    prob[nclus] = clust->get_prob();
	    chi2[nclus] = clust->get_chi2();
	    _cluster_eta[nclus] = eta;
	    _cluster_phi[nclus] = iter->second->get_phi();
	    _cluster_E[ nclus ] = iter->second->get_energy();
	    _cluster_ntower[ nclus ] =  iter->second->getNTowers();
	    cluster_Ecore[nclus] = iter->second->get_ecore();
	    if(_cluster_E[nclus] > maxE)
	      {
		bestclus = nclus;
		bcnt = _cluster_ntower[nclus];
		maxE = _cluster_E[nclus];
		int cit = 0;
		int tmaxE = 0;
		RawCluster::TowerConstRange ctowers = clust->get_towers();
		RawCluster::TowerConstIterator toweriter;
		for(toweriter = ctowers.first; toweriter != ctowers.second; toweriter++)
		  {
		    RawTower *tower = towersEM->getTower(toweriter->first);
		    RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());

		    bcten[cit] = tower->get_energy();
		    bctet[cit] = tower_geom->get_eta();
		    bctph[cit] = tower_geom->get_phi();
		    if(bcten[cit] > tmaxE) 
		      {
			bcmt = cit;
			tmaxE = bcten[cit];
		      }
		    ++cit;
		    if(cit >= 100)
		      {
			cout << "Warning: cluster found with > 100 towers, index " << nclus << endl;
			break;
		      }
		  }
	      }
	    ++nclus;
	  }
      }
  }
  EventHeaderv1 *event_header = findNode::getClass<EventHeaderv1>(topNode, "EventHeader" );
  if ( event_header ) 
    {
      npart = event_header->get_intval("npart");
      ncoll = event_header->get_intval("ncoll");    
      bimp = event_header->get_floatval("bimp");    
      std::cout << " npart / ncoll / bimp = " << npart << " / " << ncoll << " / " << bimp << std::endl;
    }

  */
  
  _tree->Fill();
  
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
