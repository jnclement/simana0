#include "MDCTreeMaker.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <iostream>


//____________________________________________________________________________..
MDCTreeMaker::MDCTreeMaker(const std::string &name):
  SubsysReco(name)
{
  //_foutname = name;  
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

  _tree->Branch("cluster_n",&_cluster_n,"cluster_n/I");
  _tree->Branch("cluster_eta",_cluster_eta,"cluster_eta[cluster_n]/F");
  _tree->Branch("cluster_phi",_cluster_phi,"cluster_phi[cluster_n]/F");
  _tree->Branch("cluster_E",_cluster_E,"cluster_E[cluster_n]/F");
  _tree->Branch("cluster_ntower",_cluster_ntower,"cluster_ntower[cluster_n]/I");

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

  std::cout << " I'm in MDCTreeMaker::process_event ! " << std::endl;

  // here is where the code goes 

  RawClusterContainer *clustersEM = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");

  {
    RawClusterContainer::ConstIterator hiter;
    RawClusterContainer::ConstRange begin_end = clustersEM->getClusters();
    for (hiter = begin_end.first; hiter != begin_end.second; ++hiter)
      {
	//float theta = 3.14159 / 2.0 - TMath::ATan2( hiter->second->get_z() , hiter->second->get_r() );
	float theta = 3.14159 / 2.0 - atan2( hiter->second->get_z() , hiter->second->get_r() );
	float eta = -1 * log( tan( theta / 2.0 ) );
	
	if ( hiter->second->get_energy() > 1 )
	  std::cout << " EM cluster eta / phi = " << eta << " / " << hiter->second->get_phi() << " , E = " << hiter->second->get_energy() << " , nTow = " << hiter->second->getNTowers()  << std::endl;
	else
	  continue;

	_cluster_eta[ _cluster_n ] = eta; 
	_cluster_phi[ _cluster_n ] = hiter->second->get_phi();
	_cluster_E[ _cluster_n ] = hiter->second->get_energy();
	_cluster_ntower[ _cluster_n ] =  hiter->second->getNTowers();
	
	_cluster_n ++;
      }
  }

  
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
