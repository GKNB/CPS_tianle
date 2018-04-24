#ifndef _COMPUTE_SIGMATOPIPI_V2_H
#define _COMPUTE_SIGMATOPIPI_V2_H

#include<alg/a2a/mesonfield.h>
#include<alg/a2a/mesonfield_computemany.h>
#include<alg/a2a/inner_product.h>
#include<alg/a2a/mf_momcontainer.h>
#include<alg/a2a/compute_pion.h>
#include<alg/a2a/compute_sigma.h>
#include<alg/a2a/mf_productstore.h>
#include<alg/a2a/required_momenta.h>

CPS_START_NAMESPACE

template<typename mf_Policies>
class computesigmapipi_v2
{
public:
typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
template<typename PionMomentumPolicy, typename SigmaMomentumPolicy>
static void compute_v2(fMatrix<ScalarComplexType> &into, MesonFieldMomentumContainer<mf_Policies> &mf_ll_con, const PionMomentumPolicy &pion_mom, const SigmaMomentumPolicy &sigma_mom, MesonFieldProductStore<mf_Policies> &products, const int tsep, const int pidxpi, const int pidxsigma, const int traj, const std::string &work_dir)
{
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
  os << work_dir << "/traj_" << traj << "_pion_mf_mom" << p_pi_src.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
  MesonFieldType::read(os.str(),mf_pi1);
  os.str("");
  os << work_dir << "/traj_" << traj << "_pion_mf_mom" << p_pi_snk.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
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
  os << work_dir << "/traj_" << traj << "_sigma_mfwv_mom" << p_sigma_src.file_str() << "_plus" << p_sigma_snk.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
  MesonFieldType::read(os.str(),mf_sigma);
  os.str("");

#ifdef NODE_DISTRIBUTE_MESONFIELDS
  if(!UniqueID()){ printf("Gathering meson fields\n");  fflush(stdout); }
  nodeGetMany(3,&mf_pi1,&mf_pi2,&mf_sigma);
  sync();
#endif


//t3-sigma; t2-pion1; t1-pion2; tdis=t2-t3; tsep=t1-t2
  for(int t3=0;t3<Lt;t3++)
  {
    for(int t2=0;t2<Lt;t2++)
    {
      int t1=(t2+tsep)%Lt;
      int tdis=(t2+Lt-t3)%Lt;
      ScalarComplexType temp1(0,0), temp2(0,0);
      temp1 += trace(mf_sigma[t3]) * trace(mf_pi1[t2], mf_pi2[t1]);
      temp1 *= ScalarComplexType(0.61237243569);
      const MesonFieldType &prod = products.getProduct(mf_pi2[t1], mf_sigma[t3],NODE_LOCAL);
      temp2 += trace(prod, mf_pi1[t2]);
      temp2 *= ScalarComplexType(1.22474487139);
      into(t3, tdis) += temp1;
      into(t3, tdis) -= temp2;
    }
  }


  sync();
  if(!UniqueID()){ printf("Finished trace\n");  fflush(stdout); }

#ifdef NODE_DISTRIBUTE_MESONFIELDS
    nodeDistributeMany(3,&mf_pi1,&mf_pi2,&mf_sigma);
#endif  
}
};

CPS_END_NAMESPACE

#endif
