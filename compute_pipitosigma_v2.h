#ifndef _COMPUTE_PIPITOSIGMA_V2_H
#define _COMPUTE_PIPITOSIGMA_V2_H

#include<alg/a2a/mesonfield.h>
#include<alg/a2a/mesonfield_computemany.h>
#include<alg/a2a/inner_product.h>
#include<alg/a2a/mf_momcontainer.h>
#include<alg/a2a/compute_pion.h>
#include<alg/a2a/compute_sigma.h>
#include<alg/a2a/mf_productstore.h>
#include<alg/a2a/required_momenta.h>



#include<iostream>


CPS_START_NAMESPACE
#define NODE_LOCAL true

template<typename mf_Policies>
class computepipisigma_v2
{
public:
typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
template<typename PionMomentumPolicy, typename SigmaMomentumPolicy>
static void compute_v2(fMatrix<ScalarComplexType> &into, MesonFieldMomentumContainer<mf_Policies> &mf_ll_con, const PionMomentumPolicy &pion_mom, const SigmaMomentumPolicy &sigma_mom, MesonFieldProductStore<mf_Policies> &products, const int tsep, const int pidxpi, const int pidxsigma, const int traj, const std::string &work_dir, bool ground=true)
{
  bool test=false;
//  test=true; printf("WE ARE IN THE TEST MODE FOR PIPI TO SIGMA CALCULATION, ALL INFORMATION WILL BE SHOWN ON SCREEN!\n"); 
  int Lt = GJP.Tnodes()*GJP.TnodeSites();
//  std::string src_names[2] = {"1s","2s"};
  into.resize(Lt,Lt);
  // now we try to get pion-mesonfield
  ThreeMomentum p_pi_src = pion_mom.getMesonMomentum(pidxpi); //exp(-ipx)
  ThreeMomentum p_pi_snk = -p_pi_src; //exp(+ipx)
  //assert(mf_ll_con.contains(p_pi_src));
  //assert(mf_ll_con.contains(p_pi_snk));

  typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MesonFieldType;
  std::vector<MesonFieldType> mf_pi1(Lt);
  std::vector<MesonFieldType> mf_pi2(Lt);
//  std::vector<MesonFieldType> &mf_pi1 = mf_ll_con.get(p_pi_src);
//  std::vector<MesonFieldType> &mf_pi2 = mf_ll_con.get(p_pi_snk);

  std::ostringstream os;
  // we can rewrite the following line to do rad not equals to 2 or 2s wavefunction
if(ground)
{
  os << work_dir << "/traj_" << traj << "_pion_mf_mom" << p_pi_src.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
}
else
{
  os << work_dir << "/traj_" << traj << "_pion_mf_mom" << p_pi_src.file_str() << "_hyd" << "2s" << "_rad" << "2" << ".dat";
}
  MesonFieldType::read(os.str(),mf_pi1);
  os.str("");
if(ground)
{
  os << work_dir << "/traj_" << traj << "_pion_mf_mom" << p_pi_snk.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
}
else
{
  os << work_dir << "/traj_" << traj << "_pion_mf_mom" << p_pi_snk.file_str() << "_hyd" << "2s" << "_rad" << "2" << ".dat";
}
  MesonFieldType::read(os.str(),mf_pi2);
  os.str("");

  // now we try to get sigma-mesonfield
  //p_sigma_snk/src indicate the momentum for quark within sigma, minus sign comes from wdag(src)
  ThreeMomentum p_sigma_src = -sigma_mom.getWmom(pidxsigma); //exp(-ipx)
  ThreeMomentum p_sigma_snk = -p_sigma_src; //exp(+ipx)
//  assert(mf_ll_con.contains(p_sigma_src));
//  assert(mf_ll_con.contains(p_sigma_snk));

  std::vector<MesonFieldType> mf_sigma(Lt);
  //std::vector<MesonFieldType> &mf_sigma = mf_ll_con.get(p_sigma_src);
  // we can rewrite the following line to do rad not equals to 2 or 2s wavefunction
if(ground)
{
  os << work_dir << "/traj_" << traj << "_sigma_mfwv_mom" << p_sigma_src.file_str() << "_plus" << p_sigma_snk.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
}
else
{
  os << work_dir << "/traj_" << traj << "_sigma_mfwv_mom" << p_sigma_src.file_str() << "_plus" << p_sigma_snk.file_str() << "_hyd" << "2s" << "_rad" << "2" << ".dat";
}
  MesonFieldType::read(os.str(),mf_sigma);
  os.str("");

//#ifdef NODE_DISTRIBUTE_MESONFIELDS
//  if(!UniqueID()){ printf("Gathering meson fields\n");  fflush(stdout); }
//  nodeGetMany(3,&mf_pi1,&mf_pi2,&mf_sigma);
//  sync();
//#endif

  MesonFieldType prod;
//float time=0;
//float time1=0;
//t3-sigma; t2-pion1; t1-pion2; tdis=t3-t2; tsep=t2-t1
  int totwork=Lt * Lt;
  int start_pt, workpernode;
  bool do_work;
  int t3,t2,t1,tdis;
  getNodeWork(totwork,workpernode,start_pt,do_work);
  printf("workpernode=\t%d\tstart_pt=\t%d\tuid=\t%d\n",workpernode,start_pt,UniqueID());
  if(do_work)
  {
    for(int tonnode=start_pt; tonnode<start_pt + workpernode; tonnode++)
    {
      t3=tonnode / Lt;
      t2=tonnode % Lt;
      t1=(t2+Lt-tsep)%Lt;
      tdis=(t3+Lt-t2)%Lt;
      ScalarComplexType temp1(0,0), temp2(0,0);
std::cout<<"initial\t"<<t3<<'\t'<<t2<<'\t'<<"uid=\t"<<UniqueID()<<'\n';
      temp1 += trace(mf_sigma[t3]) * trace(mf_pi1[t2], mf_pi2[t1]);
//time=0; time-=dclock();
      temp1 *= ScalarComplexType(0.61237243569);
//time1=time+dclock();
std::cout<<"1st graph with temp1\t"<<temp1<<'\n';
      mult(prod, mf_pi2[t1], mf_sigma[t3],NODE_LOCAL);
//time1=time+dclock()-time1;
//std::cout<<"time for generate prod mesonfield is\t\t"<<time1<<'\n';
      temp2 += trace(prod, mf_pi1[t2]);
      temp2 *= ScalarComplexType(1.22474487139);
//time1=time+dclock()-time1;
std::cout<<"2nd graph with temp2\t"<<temp2<<'\n';
      into(t2, tdis) += temp1;
      into(t2, tdis) -= temp2;
//time1=0;
//time=0;
    }
  }
  sync();
  into.nodeSum();

  if(!UniqueID()){ printf("Finished trace\n");  fflush(stdout); }

//#ifdef NODE_DISTRIBUTE_MESONFIELDS
//    nodeDistributeMany(3,&mf_pi1,&mf_pi2,&mf_sigma);
//#endif  
}
};

CPS_END_NAMESPACE

#endif
