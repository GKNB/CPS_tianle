#ifndef _COMPUTE_PION_V3_H
#define _COMPUTE_PION_V3_H


#include<alg/a2a/mesonfield.h>
#include<alg/a2a/mesonfield_computemany.h>
#include<alg/a2a/inner_product.h>
#include<alg/a2a/mf_momcontainer.h>
#include<alg/a2a/compute_pion.h>
#include<alg/a2a/required_momenta.h>

#include<iostream>
#include<string>

CPS_START_NAMESPACE
#define NODE_LOCAL true
template<typename mf_Policies>
struct computepion_v3
{
  typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
  typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MesonFieldType;

  std::array<std::string,2> pion_type{ {"1s","2s"} };
  template<typename PionMomentumPolicy>
  static void compute(fMatrix<ScalarComplexType> &into, MesonFieldMomentumContainer<mf_Policies> &mf_1s_con, MesonFieldMomentumContainer<mf_Policies> &mf_2s_con,
                      const PionMomentumPolicy &pion_mom, const int pidx, std::string src_type, std::string snk_type)
  {
    int Lt = GJP.Tnodes()*GJP.TnodeSites();
    into.resize(Lt,Lt);
    into.zero();

    if(src_type == pion_type[0])
      MesonFieldMomentumContainer<mf_Policies> &mf_src_con = mf_1s_con;
    else
    {
      assert(src_type == pion_type[1] && "Error: src_type can only be 1s or 2s\n");
      MesonFieldMomentumContainer<mf_Policies> &mf_src_con = mf_2s_con;
    }
    if(snk_type == pion_type[0])
      MesonFieldMomentumContainer<mf_Policies> &mf_snk_con = mf_1s_con;
    else
    {
      assert(snk_type == pion_type[1] && "Error: snk_type can only be 1s or 2s\n");
      MesonFieldMomentumContainer<mf_Policies> &mf_snk_con = mf_2s_con;
    }

    ThreeMomentum p_src = pion_mom.getMesonMomentum(pidx);
    ThreeMomentum p_snk = -p_src;
    assert(mf_src_con.contains(p_src));
    assert(mf_snk_con.contains(p_snk));

    std::vector<MesonFieldType> &mf_pi_src = mf_src_con.get(p_src);
    std::vector<MesonFieldType> &mf_pi_snk = mf_snk_con.get(p_snk);

    int totwork = Lt * Lt;
    int start_pt, workpernode;
    bool do_work;
    getNodeWork(totwork,workpernode,start_pt,do_work);
#ifdef NODE_DISTRIBUTE_MESONFIELDS
    nodeGetMany(2,&mf_pi_src,&mf_pi_snk);
    cps::sync();
#endif
    if(do_work)
    {
      for(int tonnode = start_pt; tonnode < start_pt + workpernode; tonnode++)
      {
	int tsnk = tonnode % Lt;
	int tsrc = tonnode / Lt;
	int tdis = (tsnk - tsrc + Lt ) % Lt;
	into(tsrc,tdis) += ScalarComplexType(0.5) * trace(mf_pi_src[tsrc], mf_pi_snk[tsnk]);
      }
    }
    sync();
    into.nodeSum();
    if(!UniqueID()){ printf("Finished trace\n");  fflush(stdout); }
#ifdef NODE_DISTRIBUTE_MESONFIELDS
    nodeDistributeMany(2,&mf_pi_src,&mf_pi_snk);
    cps::sync();
#endif
  }
};
CPS_END_NAMESPACE

#endif
