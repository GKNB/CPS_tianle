#ifndef _COMPUTE_SIGMA_V3_H
#define _COMPUTE_SIGMA_V3_H


#include<alg/a2a/mesonfield.h>
#include<alg/a2a/mesonfield_computemany.h>
#include<alg/a2a/inner_product.h>
#include<alg/a2a/mf_momcontainer.h>
#include<alg/a2a/compute_sigma.h>
#include<alg/a2a/mf_productstore.h>
#include<alg/a2a/required_momenta.h>

#include<iostream>

CPS_START_NAMESPACE
#define NODE_LOCAL true
template<typename mf_Policies>
class computesigma_v3
{
public:
  typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
  typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MesonFieldType;

  template<typename SigmaMomentumPolicy>
  static void compute(fMatrix<ScalarComplexType> &into, fVector<ScalarComplexType> &intoself, MesonFieldMomentumPairContainer<mf_Policies> &mf_sigma_con, 
                      const SigmaMomentumPolicy &sigma_mom, const int psrcidx, const int psnkidx)
  {
    int Lt = GJP.Tnodes()*GJP.TnodeSites();
    into.resize(Lt,Lt);
    intoself.resize(Lt);

//p_sigma_snk/src indicate the momentum for quark within sigma, minus sign comes from wdag
    ThreeMomentum p1_src = -sigma_mom.getWmom(psrcidx);    ThreeMomentum p2_src = sigma_mom.getVmom(psrcidx);
    ThreeMomentum p1_snk = -sigma_mom.getWmom(psnkidx);    ThreeMomentum p2_snk = sigma_mom.getVmom(psnkidx);
    assert(mf_sigma_con.contains(p1_src,p2_src));    assert(mf_sigma_con.contains(p1_snk,p2_snk));

    std::vector<MesonFieldType> &mf_sigmasrc = mf_sigma_con.get(p1_src,p2_src);
    std::vector<MesonFieldType> &mf_sigmasnk = mf_sigma_con.get(p1_snk,p2_snk);
#ifdef NODE_DISTRIBUTE_MESONFIELDS
    if(!UniqueID()){ printf("Gathering meson fields\n");  fflush(stdout); }
    nodeGetMany(2,&mf_sigmasrc,&mf_sigmasnk);
    cps::sync();
#endif

    if(!UniqueID()){ printf("Starting trace\n");  fflush(stdout); }
    int totwork = Lt * Lt;
    int start_pt, workpernode;
    bool do_work;
    getNodeWork(totwork,workpernode,start_pt,do_work);
    if(do_work)
    {
      for(int tonnode = start_pt; tonnode < start_pt + workpernode; tonnode++)
      {
        int t1 = tonnode % Lt;
        int t2 = tonnode / Lt;
        int tdis = (t1 - t2 + Lt) % Lt;
        ScalarComplexType sigmatemp(0,0);
        sigmatemp += trace(mf_sigmasnk[t1]) * trace(mf_sigmasrc[t2]);
        sigmatemp -= trace(mf_sigmasnk[t1] , mf_sigmasrc[t2]);
        sigmatemp *= ScalarComplexType(0.5);
        into(t2,tdis) += sigmatemp;
      }
    }
    sync();
    into.nodeSum();

    if(psnkidx == 0)
    {
      totwork = Lt;
      getNodeWork(totwork,workpernode,start_pt,do_work);
      if(do_work)
        for(int tonnode = start_pt; tonnode < start_pt + workpernode; tonnode++)
          intoself(tonnode) += trace(mf_sigmasrc[tonnode]);
      sync();
      intoself.nodeSum();
    }
    if(!UniqueID()){ printf("Finished trace\n");  fflush(stdout); }
#ifdef NODE_DISTRIBUTE_MESONFIELDS
    nodeDistributeMany(2,&mf_sigmasrc,&mf_sigmasnk);
    cps::sync();
#endif
  }
};

CPS_END_NAMESPACE

#endif     
