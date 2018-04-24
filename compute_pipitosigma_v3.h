#ifndef _COMPUTE_PIPITOSIGMA_V3_H
#define _COMPUTE_PIPITOSIGMA_V3_H


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
//define G5-pipisigma means do extra topology
#undef G5_PIPISIGMA	
template<typename mf_Policies>
class computepipisigma_v3
{
public:
  typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
  typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MesonFieldType;

  template<typename PionMomentumPolicy, typename SigmaMomentumPolicy>
  static void compute(fMatrix<ScalarComplexType> &into, MesonFieldMomentumContainer<mf_Policies> &mf_pion_con, MesonFieldMomentumPairContainer<mf_Policies> &mf_sigma_con,
                      const PionMomentumPolicy &pion_mom, const SigmaMomentumPolicy &sigma_mom, const int pidxpi, const int pidxsigma, const int tsep)
  {
    int Lt = GJP.Tnodes()*GJP.TnodeSites();
    into.resize(Lt,Lt);
  
    ThreeMomentum p_pi_src1 = pion_mom.getMesonMomentum(pidxpi); //define inner pion as pi1
    ThreeMomentum p_pi_src2 = -p_pi_src1; 
    ThreeMomentum p1_snk = -sigma_mom.getWmom(pidxsigma);    ThreeMomentum p2_snk = sigma_mom.getVmom(pidxsigma);
    assert(mf_pion_con.contains(p_pi_src1));
    assert(mf_pion_con.contains(p_pi_src2));
    assert(mf_sigma_con.contains(p1_snk,p2_snk));
  
    std::vector<MesonFieldType> &mf_pi1 = mf_pion_con.get(p_pi_src1);
    std::vector<MesonFieldType> &mf_pi2 = mf_pion_con.get(p_pi_src2);
    std::vector<MesonFieldType> &mf_sigma = mf_sigma_con.get(p1_snk,p2_snk);
#ifdef NODE_DISTRIBUTE_MESONFIELDS
    if(!UniqueID()){ printf("Gathering meson fields\n");  fflush(stdout); }
    nodeGetMany(3,&mf_pi1,&mf_pi2,&mf_sigma);
    cps::sync();
#endif

    if(!UniqueID()){ printf("Starting trace\n");  fflush(stdout); }
//t3-sigma; t2-pion1; t1-pion2; tdis=t3-t2; tsep=t2-t1
    int totwork=Lt * Lt;
    int start_pt, workpernode;
    bool do_work;
    getNodeWork(totwork,workpernode,start_pt,do_work);
    if(do_work)
    {
      for(int tonnode = start_pt; tonnode < start_pt + workpernode; tonnode++)
      {
        int t3 = tonnode % Lt;
        int t2 = tonnode / Lt;
        int t1 = (t2 + Lt - tsep) % Lt;
        int tdis =(t3 + Lt - t2) % Lt;

        ScalarComplexType temp1(0,0), temp2(0,0);
        temp1 += trace(mf_sigma[t3]) * trace(mf_pi1[t2], mf_pi2[t1]);
        temp1 *= (sqrt(6.0)/4.0);
#ifndef G5_PIPISIGMA
        MesonFieldType prod;
        mult(prod, mf_pi2[t1], mf_sigma[t3],NODE_LOCAL);
        temp2 += trace(prod, mf_pi1[t2]);
        temp2 *= (sqrt(6.0)/2.0);
#else
        MesonFieldType prod1, prod2;
        mult(prod1, mf_pi2[t1], mf_sigma[t3],NODE_LOCAL);
        mult(prod2, mf_sigma[t3], mf_pi2[t1],NODE_LOCAL);
        temp2 += trace(prod, mf_pi1[t2]);
        temp2 += trace(prod, mf_pi1[t2]);
        temp2 *= (sqrt(6.0)/4.0);
#endif
        into(t2, tdis) += temp1;
        into(t2, tdis) -= temp2;
      }
    }
    sync();
    into.nodeSum();
    if(!UniqueID()){ printf("Finished trace\n");  fflush(stdout); }

#ifdef NODE_DISTRIBUTE_MESONFIELDS
    nodeDistributeMany(3,&mf_pi1,&mf_pi2,&mf_sigma);
    cps::sync();
#endif  
  }
};

CPS_END_NAMESPACE

#endif
