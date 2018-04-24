#ifndef COMPUTE_COMOVE_PIPI_V2_H
#define COMPUTE_COMOVE_PIPI_V2_H


#include<alg/a2a/mesonfield.h>
#include<alg/a2a/mesonfield_computemany.h>
#include<alg/a2a/inner_product.h>
#include<alg/a2a/mf_momcontainer.h>
#include<alg/a2a/compute_pion.h>
#include<alg/a2a/mf_productstore.h>
#include<alg/a2a/required_momenta.h>

#include<iostream>

CPS_START_NAMESPACE
#define NODE_LOCAL true
template<typename mf_Policies>
class compute_comove_pipi_v2
{
  public:
    typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
    typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MesonFieldType;

    //here into is a vector
    template<typename PionMomentumPolicy>
    static void compute_Vdis_v2(fVector<ScalarComplexType> &into, MesonFieldMomentumContainer<mf_Policies> &mf_pi_con,
	                        const PionMomentumPolicy &pion_mom, const int psrc, const int psnk, const int tsep)
    {
      int Lt = GJP.Tnodes()*GJP.TnodeSites();
      into.resize(Lt);
      into.zero();

      ThreeMomentum p_src = pion_mom.getMesonMomentum(psrc);
      ThreeMomentum p_snk = pion_mom.getMesonMomentum(psnk);
      assert(mf_pi_con.contains(p_src));
      assert(mf_pi_con.contains(p_snk));
      
      std::vector<MesonFieldType> &mf_pi_src = mf_pi_con.get(p_src);
      std::vector<MesonFieldType> &mf_pi_snk = mf_pi_con.get(p_snk);
#ifdef NODE_DISTRIBUTE_MESONFIELDS
      nodeGetMany(2,&mf_pi_src,&mf_pi_snk);
      cps::sync();
#endif
      int totwork = Lt;
      int start_pt, workpernode;
      bool do_work;
      getNodeWork(totwork,workpernode,start_pt,do_work);
      if(do_work)
      {
        for(int tsnk = start_pt; tsnk < start_pt + workpernode; tsnk++)
        {
          int tsrc = (tsnk - tsep + Lt) % Lt;
          into(tsnk) = ScalarComplexType(0.5) * trace(mf_pi_src[tsrc], mf_pi_snk[tsnk]);
        }
      }
      sync();
      into.nodeSum();
      if(!UniqueID()){ printf("Finished Vdis trace\n");  fflush(stdout); }
#ifdef NODE_DISTRIBUTE_MESONFIELDS
      nodeDistributeMany(2,&mf_pi_src,&mf_pi_snk);
#endif
    }

    //here into is an array of 3 matrix, they are C,D,R diagram
    //in type2,3, we need to do it generally
    template<typename PionMomentumPolicy>
    static void compute_type23_v2(std::array<fMatrix<ScalarComplexType>,3> &into, MesonFieldMomentumContainer<mf_Policies> &mf_pi_con,
	                         const PionMomentumPolicy &pion_mom, const int psrc1, const int psrc2, const int psnk1, const int psnk2,
				 const int tsep, const int tstep_src)
    {
      int Lt = GJP.Tnodes()*GJP.TnodeSites();
      for(int i=0; i<3; i++)
      {
        into[i].resize(Lt,Lt);
        into[i].zero();
      }
      if(Lt % tstep_src != 0) {assert( 0 && "Lt has to be an integer times tstep_src");}

      ThreeMomentum p_src1 = pion_mom.getMesonMomentum(psrc1);
      ThreeMomentum p_src2 = pion_mom.getMesonMomentum(psrc2);
      ThreeMomentum p_snk1 = pion_mom.getMesonMomentum(psnk1);
      ThreeMomentum p_snk2 = pion_mom.getMesonMomentum(psnk2);
      assert((p_src1 + p_src2) == (p_snk1 + p_snk2) || (p_src1 + p_src2) == -(p_snk1 + p_snk2));
      ThreeMomentum p_tmp = p_src1 + p_src2;
      assert(p_tmp.nzero() != 0 && p_tmp.nzero() < 3);
      assert(mf_pi_con.contains(p_src1));
      assert(mf_pi_con.contains(p_src2));
      assert(mf_pi_con.contains(p_snk1));
      assert(mf_pi_con.contains(p_snk2));

      int totwork = Lt * Lt / tstep_src;
      int start_pt, workpernode;
      bool do_work;
      getNodeWork(totwork,workpernode,start_pt,do_work);

      //t4=tsnk+tsep,t3=tsnk,t2=tsrc,t1=tsrc-tsep,tdis=t3-t2
      int t1=0,t2=0,t3=0,t4=0;
      int tdis=0;

      std::vector<MesonFieldType> &mf_pi_src1 = mf_pi_con.get(p_src1);
      std::vector<MesonFieldType> &mf_pi_src2 = mf_pi_con.get(p_src2);
      std::vector<MesonFieldType> &mf_pi_snk1 = mf_pi_con.get(p_snk1);
      std::vector<MesonFieldType> &mf_pi_snk2 = mf_pi_con.get(p_snk2);
#ifdef NODE_DISTRIBUTE_MESONFIELDS
      nodeGetMany(4,&mf_pi_src1,&mf_pi_src2,&mf_pi_snk1,&mf_pi_snk2);
      cps::sync();
#endif
      if(do_work)
      {
        for(int tonnode = start_pt; tonnode < start_pt + workpernode; tonnode++)
        {
          t3 = tonnode % Lt;
          t2 = (tonnode / Lt) * tstep_src;
          tdis = (t3 - t2 + Lt) % Lt;
          t1 = (t2 - tsep + Lt) % Lt;
          t4 = (t3 + tsep) % Lt;

          MesonFieldType prod1, prod2, prod7, prod8;
          mult(prod1, mf_pi_src2[t1], mf_pi_snk2[t4], NODE_LOCAL);
          mult(prod2, mf_pi_src1[t2], mf_pi_snk1[t3], NODE_LOCAL);
          into[0](t2, tdis) += ScalarComplexType(0.25) * trace(prod1, prod2);
          mult(prod7, mf_pi_snk2[t4], mf_pi_src2[t1], NODE_LOCAL);                  
          mult(prod8, mf_pi_snk1[t3], mf_pi_src1[t2], NODE_LOCAL);                  
          into[0](t2, tdis) += ScalarComplexType(0.25) * trace(prod7, prod8);

          ScalarComplexType temp(0,0);
          temp += trace(mf_pi_src2[t1], mf_pi_snk1[t3]) * trace(mf_pi_src1[t2], mf_pi_snk2[t4]);
          temp += trace(mf_pi_src2[t1], mf_pi_snk2[t4]) * trace(mf_pi_src1[t2], mf_pi_snk1[t3]);
          temp *= ScalarComplexType(0.125);
          into[1](t2, tdis) += temp;

          MesonFieldType prod3, prod4, prod5, prod6;
          ScalarComplexType temp2(0,0);
          mult(prod3, mf_pi_src1[t2], mf_pi_src2[t1], NODE_LOCAL);
          mult(prod4, mf_pi_snk1[t3], mf_pi_snk2[t4], NODE_LOCAL);
          temp2 += ScalarComplexType(0.125) * trace(prod3, prod4);
          mult(prod5, mf_pi_snk2[t4], mf_pi_snk1[t3], NODE_LOCAL);
          temp2 += ScalarComplexType(0.125) * trace(prod3, prod5);
          mult(prod6, mf_pi_src2[t1], mf_pi_src1[t2], NODE_LOCAL);
          temp2 += ScalarComplexType(0.125) * trace(prod6, prod4);
          temp2 += ScalarComplexType(0.125) * trace(prod6, prod5);
          into[2](t2, tdis) += temp2;
        }
      }
      sync();
      for(int i=0; i<3; i++)
        into[i].nodeSum();
      if(!UniqueID()){ printf("Finished trace\n");  fflush(stdout); }
#ifdef NODE_DISTRIBUTE_MESONFIELDS
      nodeDistributeMany(4,&mf_pi_src1,&mf_pi_src2,&mf_pi_snk1,&mf_pi_snk2);
#endif
    }

    //here into is an array of 3 matrix, they are C,D,R diagram
    //in type1, two src and two snk and the same, so we can simplify it 
    template<typename PionMomentumPolicy>
    static void compute_type1_v2(std::array<fMatrix<ScalarComplexType>,3> &into, MesonFieldMomentumContainer<mf_Policies> &mf_pi_con,
	                         const PionMomentumPolicy &pion_mom, const int psrc1, const int psrc2, const int psnk1, const int psnk2,
				 const int tsep, const int tstep_src)
    {
      int Lt = GJP.Tnodes()*GJP.TnodeSites();
      for(int i=0; i<3; i++)
      {
	into[i].resize(Lt,Lt);
	into[i].zero();
      }
      if(Lt % tstep_src != 0) {assert( 0 && "Lt has to be an integer times tstep_src");}
      
      assert(psrc1 == psrc2);
      assert(psnk1 == psnk2);

      ThreeMomentum p_src1 = pion_mom.getMesonMomentum(psrc1);
      ThreeMomentum p_snk1 = pion_mom.getMesonMomentum(psnk1);
      assert(p_src1 == p_snk1 || p_src1 == -p_snk1);
      assert(mf_pi_con.contains(p_src1));
      assert(mf_pi_con.contains(p_snk1));

      int totwork = Lt * Lt / tstep_src;
      int start_pt, workpernode;
      bool do_work;
      getNodeWork(totwork,workpernode,start_pt,do_work);

      //t4=tsnk+tsep,t3=tsnk,t2=tsrc,t1=tsrc-tsep,tdis=t3-t2
      int t1=0,t2=0,t3=0,t4=0;
      int tdis=0;

      std::vector<MesonFieldType> &mf_pi_src = mf_pi_con.get(p_src1);
      std::vector<MesonFieldType> &mf_pi_snk = mf_pi_con.get(p_snk1);
#ifdef NODE_DISTRIBUTE_MESONFIELDS
      nodeGetMany(2,&mf_pi_src,&mf_pi_snk);
      cps::sync();
#endif
      if(do_work)
      {
        for(int tonnode = start_pt; tonnode < start_pt + workpernode; tonnode++)
        {
          t3 = tonnode % Lt;
          t2 = (tonnode / Lt) * tstep_src;
          tdis = (t3 - t2 + Lt) % Lt;
          t1 = (t2 - tsep + Lt) % Lt;
          t4 = (t3 + tsep) % Lt;

          MesonFieldType prod1, prod2, prod7, prod8;
          mult(prod1, mf_pi_src[t1], mf_pi_snk[t4], NODE_LOCAL);
          mult(prod2, mf_pi_src[t2], mf_pi_snk[t3], NODE_LOCAL);
          into[0](t2, tdis) += ScalarComplexType(0.25) * trace(prod1, prod2);
          mult(prod7, mf_pi_snk[t4], mf_pi_src[t1], NODE_LOCAL);
          mult(prod8, mf_pi_snk[t3], mf_pi_src[t2], NODE_LOCAL);
          into[0](t2, tdis) += ScalarComplexType(0.25) * trace(prod7, prod8);

          ScalarComplexType temp(0,0);
          temp += trace(mf_pi_src[t1], mf_pi_snk[t3]) * trace(mf_pi_src[t2], mf_pi_snk[t4]);
          temp += trace(mf_pi_src[t1], mf_pi_snk[t4]) * trace(mf_pi_src[t2], mf_pi_snk[t3]);
          temp *= ScalarComplexType(0.125);
          into[1](t2, tdis) += temp;

          MesonFieldType prod3, prod4, prod5, prod6;
          ScalarComplexType temp2(0,0);
          mult(prod3, mf_pi_src[t2], mf_pi_src[t1], NODE_LOCAL);
          mult(prod4, mf_pi_snk[t3], mf_pi_snk[t4], NODE_LOCAL);
          temp2 += ScalarComplexType(0.125) * trace(prod3, prod4);
          mult(prod5, mf_pi_snk[t4], mf_pi_snk[t3], NODE_LOCAL);
          temp2 += ScalarComplexType(0.125) * trace(prod3, prod5);
          mult(prod6, mf_pi_src[t1], mf_pi_src[t2], NODE_LOCAL);
          temp2 += ScalarComplexType(0.125) * trace(prod6, prod4);
          temp2 += ScalarComplexType(0.125) * trace(prod6, prod5);
          into[2](t2, tdis) += temp2;
        }
      }
      sync();
      for(int i=0; i<3; i++)
        into[i].nodeSum();
      if(!UniqueID()){ printf("Finished trace\n");  fflush(stdout); }
#ifdef NODE_DISTRIBUTE_MESONFIELDS
      nodeDistributeMany(2,&mf_pi_src,&mf_pi_snk);
#endif
    }
};

CPS_END_NAMESPACE
#endif
