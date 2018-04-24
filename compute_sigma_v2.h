#ifndef _COMPUTE_SIGMA_V2_H
#define _COMPUTE_SIGMA_V2_H


#include<alg/a2a/mesonfield_computemany.h>
#include<alg/a2a/inner_product.h>
#include<alg/a2a/mf_momcontainer.h>
#include<alg/a2a/compute_sigma.h>

CPS_START_NAMESPACE

template<typename mf_Policies>
class computesigma_v2
{
  public:
    typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
    template<typename SigmaMomentumPolicy>
    static void compute_v2(fMatrix<ScalarComplexType> &into, fVector<ScalarComplexType> &intoself, MesonFieldMomentumContainer<mf_Policies> &mf_ll_con, const SigmaMomentumPolicy &sigma_mom, const int psrcidx, const int psnkidx, const int traj, const std::string &work_dir, bool ground=true, bool testsigma=true)
    {
//split_by_traj==false means we distribute work onto node according to time (split_by_time)
//split_by_traj==true means we did all trace calcualtion on one node, and is designed for split job by traj_num
//    bool split_by_traj=false;
      if(testsigma)
	std::cout<<"WE ARE IN TEST MODE FOR 'sigma cor', ALL INFORMATION WILL BE SHOWN ON SCREEN"<<'\n';
      int Lt = GJP.Tnodes()*GJP.TnodeSites();
//    std::string src_names[2] = {"1s","2s"};
      into.resize(Lt,Lt);
      intoself.resize(Lt);

//p_sigma_snk/src indicate the momentum for quark within sigma, minus sign comes from wdag
      ThreeMomentum p_sigma_src = -sigma_mom.getWmom(psrcidx); //exp(-ipx)
      ThreeMomentum p_sigma_snk = -sigma_mom.getWmom(psnkidx); //exp(+ipx)

      typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MesonFieldType;
      std::vector<MesonFieldType> mf_sigmasrc(Lt);
      std::vector<MesonFieldType> mf_sigmasnk(Lt);
      std::ostringstream os1;
      std::ostringstream os2;
//we can rewrite the following line to do rad not equals to 2 or 2s wavefunction
//      if(split_by_traj)
//      {
//        os << work_dir << "/traj_" << traj + UniqueID() * split_by_traj << "_sigma_mfwv_mom" << p_sigma_src.file_str() << "_plus" << p_sigma_snk.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
//      }
//      else
//      {
        if(ground)
        {
          os1<< work_dir << "/traj_" << traj << "_sigma_mfwv_mom" << p_sigma_src.file_str() << "_plus" << (-p_sigma_src).file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
	  os2<< work_dir << "/traj_" << traj << "_sigma_mfwv_mom" << p_sigma_snk.file_str() << "_plus" << (-p_sigma_snk).file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
        }
        else
        {
          os1<< work_dir << "/traj_" << traj << "_sigma_mfwv_mom" << p_sigma_src.file_str() << "_plus" << (-p_sigma_src).file_str() << "_hyd" << "2s" << "_rad" << "2" << ".dat";
          os2<< work_dir << "/traj_" << traj << "_sigma_mfwv_mom" << p_sigma_snk.file_str() << "_plus" << (-p_sigma_snk).file_str() << "_hyd" << "2s" << "_rad" << "2" << ".dat";
        }
//      }
//that seems redundant, maybe we don't need if statement ??????? 
      MesonFieldType::read(os1.str(),mf_sigmasrc);
      MesonFieldType::read(os2.str(),mf_sigmasnk);
      os1.str("");
      os2.str("");
//#ifdef NODE_DISTRIBUTE_MESONFIELDS
//  if(!UniqueID()){ printf("Gathering meson fields\n");  fflush(stdout); }
//  nodeGetMany(2,&mf_pi1,&mf_pi2);
//  sync();
//#endif
      if(!UniqueID()){ printf("Starting trace\n");  fflush(stdout); }
      int totwork=Lt * Lt;
      int start_pt, workpernode;
      bool do_work;
      getNodeWork(totwork,workpernode,start_pt,do_work);
      if(do_work)
      {
        for(int tonnode=start_pt; tonnode<start_pt + workpernode; tonnode++)
        {
          int t1=tonnode / Lt;
          int t2=tonnode % Lt;
	  int tdis = (t1-t2+Lt)%Lt;
          ScalarComplexType sigmatemp(0,0);
          sigmatemp += trace(mf_sigmasnk[t1]) * trace(mf_sigmasrc[t2]);
          if(testsigma)
          {
          ScalarComplexType temp(0,0);
          temp = trace(mf_sigmasnk[t1] , mf_sigmasrc[t2]);
          std::cout<<"t1= "<<t1<<"\tt2= "<<t2<<'\n';
          std::cout<<"disconnected=\t"<<sigmatemp<<'\n';
          std::cout<<"connected=\t"<<temp<<'\n';
          }
          sigmatemp -= trace(mf_sigmasnk[t1] , mf_sigmasrc[t2]);
          sigmatemp *= ScalarComplexType(0.5);
          into(t2,tdis) += sigmatemp;
        }
      }
      sync();
      into.nodeSum();
if(psnkidx==0)
{
      totwork=Lt;
      getNodeWork(totwork,workpernode,start_pt,do_work);
      if(do_work)
      {
        for(int tonnode=start_pt; tonnode<start_pt + workpernode; tonnode++)
        {
          intoself(tonnode) += trace(mf_sigmasrc[tonnode]);
          if(testsigma)
          {
          std::cout<<"self_matrix = "<<intoself(tonnode)<<'\n';
          }
        }
      }
      sync();
      intoself.nodeSum();
      if(!UniqueID()){ printf("Finished trace\n");  fflush(stdout); }
}
//#ifdef NODE_DISTRIBUTE_MESONFIELDS
//    nodeDistributeMany(2,&mf_pi1,&mf_pi2);
//#endif  
    }
};

CPS_END_NAMESPACE

#endif     
