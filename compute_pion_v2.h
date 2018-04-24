#ifndef _COMPUTE_PION_V2_H
#define _COMPUTE_PION_V2_H


#include<alg/a2a/mesonfield_computemany.h>
#include<alg/a2a/inner_product.h>
#include<alg/a2a/mf_momcontainer.h>
#include<alg/a2a/compute_pion.h>

CPS_START_NAMESPACE

template<typename mf_Policies>
class computepion_v2
{
public:
typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
template<typename PionMomentumPolicy>
static void compute_v2(fMatrix<ScalarComplexType> &into, MesonFieldMomentumContainer<mf_Policies> &mf_ll_con, const PionMomentumPolicy &pion_mom, const int pidx, const int traj, const std::string &work_dir, bool ground=true)
{
  bool test =false;
//split_by_traj==false means we distribute work onto node according to time (split_by_time)
//split_by_traj==true means we did all trace calcualtion on one node, and is designed for split job by traj_num
  bool split_by_traj=false;
//test = true; std::cout<<"WE ARE IN 'pion cor' TEST MODE, ALL INFORMATION WILL BE SHOWN ON SCREEN!"<<'\n';

  int Lt = GJP.Tnodes()*GJP.TnodeSites();  
//std::string src_names[2] = {"1s","2s"};
  into.resize(Lt,Lt);
  
  ThreeMomentum p_pi_src = pion_mom.getMesonMomentum(pidx); //exp(-ipx)
  ThreeMomentum p_pi_snk = -p_pi_src; //exp(+ipx)
  //assert(mf_ll_con.contains(p_pi_src));
  //assert(mf_ll_con.contains(p_pi_snk));
  
  typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MesonFieldType;
  std::vector<MesonFieldType> mf_pi1(Lt);
  std::vector<MesonFieldType> mf_pi2(Lt);
  //std::vector<MesonFieldType> &mf_pi1 = mf_ll_con.get(p_pi_src);
  //std::vector<MesonFieldType> &mf_pi2 = mf_ll_con.get(p_pi_snk);
  std::ostringstream os; 
  // we can rewrite the following line to do rad not equals to 2 or 2s wavefunction
  if(!split_by_traj)
{
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
//#ifdef NODE_DISTRIBUTE_MESONFIELDS
//  if(!UniqueID()){ printf("Gathering meson fields\n");  fflush(stdout); }
//  nodeGetMany(2,&mf_pi1,&mf_pi2);
//  sync();
//#endif
  if(!UniqueID()){ printf("Starting trace\n");  fflush(stdout); }
  trace(into,mf_pi2,mf_pi1);

  into *= ScalarComplexType(0.5,0);
  sync();

//here we test mf_pi self loop, which should be 0 theoretically, and if two ways of trace calculations give the same results
if(test)
{
ScalarComplexType temp1=0, temp2=0, temp3=0;
for(int t1=0; t1<Lt; t1++)
{
  for(int t2=0; t2<Lt; t2++)
  {
    if(!UniqueID())
    {
      temp1= trace(mf_pi1[t1]);
      temp2= trace(mf_pi1[t2]);
      std::cout<<"t1= "<<t1<<"\tt2= "<<t2<<'\n'<<"temp1= "<<temp1<<'\n'<<"temp2= "<<temp2<<'\n';
      std::cout<<"temp1 * temp2=\t\t\t\t"<<temp1*temp2<<'\n';
      std::cout<<"the theoretical value is\t\t"<<into(t1,t2)*2.0<<'\n';
      temp3=trace(mf_pi2[t1],mf_pi1[t2]);
      std::cout<<"another way to get real value is\t"<<temp3<<'\n';
    }
  }
}
}
  rearrangeTsrcTsep(into); //rearrange temporal ordering
  if(!UniqueID()){ printf("Finished trace\n");  fflush(stdout); }

//#ifdef NODE_DISTRIBUTE_MESONFIELDS
//    nodeDistributeMany(2,&mf_pi1,&mf_pi2);
//#endif  
}


  else
{
  os << work_dir << "/traj_" << traj + UniqueID() * split_by_traj << "_pion_mf_mom" << p_pi_src.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
  MesonFieldType::read(os.str(),mf_pi1);
  os.str("");
  os << work_dir << "/traj_" << traj + UniqueID() * split_by_traj << "_pion_mf_mom" << p_pi_snk.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
  MesonFieldType::read(os.str(),mf_pi2);
  os.str("");
//#ifdef NODE_DISTRIBUTE_MESONFIELDS
//  if(!UniqueID()){ printf("Gathering meson fields\n");  fflush(stdout); }
//  nodeGetMany(2,&mf_pi1,&mf_pi2);
//  sync();
//#endif
  if(!UniqueID()){ printf("Starting trace\n");  fflush(stdout); }
  for(int t1=0; t1<Lt; t1++)
  {
    for(int t2=0; t2<Lt; t2++)
    {
      into(t1,t2) += trace(mf_pi2[t1] , mf_pi1[t2]);
    }
  }
  into *= ScalarComplexType(0.5,0);
  rearrangeTsrcTsep(into); //rearrange temporal ordering
  if(!UniqueID()){ printf("Finished trace\n");  fflush(stdout); }
}
}
};
CPS_END_NAMESPACE

#endif
