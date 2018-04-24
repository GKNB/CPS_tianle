#ifndef _COMPUTE_PIPI_V2_H
#define _COMPUTE_PIPI_V2_H


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
class computepipi_v2
{
  public:
    typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
    typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MesonFieldType;

    //here into is a vector
    template<typename PionMomentumPolicy>
    static void compute_Vdis_v2(fVector<ScalarComplexType> &into, MesonFieldMomentumContainer<mf_Policies> &mf_1s_con, MesonFieldMomentumContainer<mf_Policies> &mf_2s_con,
                                const PionMomentumPolicy &pion_mom, const int psnk, const int tsep, std::string mf_type)
    {
      const std::array<std::string,2> pion_type{ {"1s","2s"} };
      int Lt = GJP.Tnodes()*GJP.TnodeSites();
      into.resize(Lt);
      into.zero();

      MesonFieldMomentumContainer<mf_Policies> *mf_con_ptr;
      if(mf_type == pion_type[0])
        mf_con_ptr = &mf_1s_con;
      else
      {   
        assert(mf_type == pion_type[1] && "Error: mf_type can only be 1s or 2s\n");
        mf_con_ptr = &mf_2s_con;
      }

      ThreeMomentum p_snk = pion_mom.getMesonMomentum(psnk);
      ThreeMomentum p_src = -p_snk;
      assert(mf_con_ptr->contains(p_src));
      assert(mf_con_ptr->contains(p_snk));

      std::vector<MesonFieldType> &mf_pi_src = mf_con_ptr->get(p_src);
      std::vector<MesonFieldType> &mf_pi_snk = mf_con_ptr->get(p_snk);
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


    //tsep is b/w two pion, into is an array of 3 matrix, they are C,D,R diagram
    template<typename PionMomentumPolicy>
    static void compute_v2(std::array<fMatrix<ScalarComplexType>,3> &into, MesonFieldMomentumContainer<mf_Policies> &mf_1s_con, 
                           MesonFieldMomentumContainer<mf_Policies> &mf_2s_con, const PionMomentumPolicy &pion_mom, const int psrc, 
			   const int psnk, const int tsep, const int tstep_src, std::string src_type, std::string snk_type)
    {
      const std::array<std::string,2> pion_type{ {"1s","2s"} };
      int Lt = GJP.Tnodes()*GJP.TnodeSites();
      for(int i=0; i<3; i++)
      {
        into[i].resize(Lt,Lt);
        into[i].zero();
      }
      if(Lt % tstep_src != 0) {assert( 0 && "Lt has to be an integer times tstep_src");}

      MesonFieldMomentumContainer<mf_Policies> *mf_src_con_ptr, *mf_snk_con_ptr;
      if(src_type == pion_type[0])
        mf_src_con_ptr = &mf_1s_con;
      else
      {   
        assert(src_type == pion_type[1] && "Error: src_type can only be 1s or 2s\n");
        mf_src_con_ptr = &mf_2s_con;
      }   
      if(snk_type == pion_type[0])
        mf_snk_con_ptr = &mf_1s_con;
      else
      {   
        assert(snk_type == pion_type[1] && "Error: snk_type can only be 1s or 2s\n");
        mf_snk_con_ptr = &mf_2s_con;
      }

      ThreeMomentum p_src1 = pion_mom.getMesonMomentum(psrc);
      ThreeMomentum p_src2 = -p_src1;
      ThreeMomentum p_snk1 = pion_mom.getMesonMomentum(psnk);
      ThreeMomentum p_snk2 = -p_snk1;

      assert(mf_src_con_ptr->contains(p_src1));
      assert(mf_src_con_ptr->contains(p_src2));
      assert(mf_snk_con_ptr->contains(p_snk1));
      assert(mf_snk_con_ptr->contains(p_snk2));

      int totwork = Lt * Lt / tstep_src;
      int start_pt, workpernode;
      bool do_work;
      getNodeWork(totwork,workpernode,start_pt,do_work);

      //t4=tsnk+tsep,t3=tsnk,t2=tsrc,t1=tsrc-tsep,tdis=t3-t2
      int t1=0,t2=0,t3=0,t4=0;
      int tdis=0;

      std::vector<MesonFieldType> &mf_pi_src1 = mf_src_con_ptr->get(p_src1);
      std::vector<MesonFieldType> &mf_pi_src2 = mf_src_con_ptr->get(p_src2);
      std::vector<MesonFieldType> &mf_pi_snk1 = mf_snk_con_ptr->get(p_snk1);
      std::vector<MesonFieldType> &mf_pi_snk2 = mf_snk_con_ptr->get(p_snk2);
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

	  MesonFieldType prod1, prod2;
	  mult(prod1, mf_pi_src2[t1], mf_pi_snk2[t4], NODE_LOCAL);
	  mult(prod2, mf_pi_src1[t2], mf_pi_snk1[t3], NODE_LOCAL);
	  into[0](t2, tdis) += ScalarComplexType(0.5) * trace(prod1, prod2);
 
	  ScalarComplexType temp(0,0);
	  temp += trace(mf_pi_src2[t1], mf_pi_snk1[t3]) * trace(mf_pi_src1[t2], mf_pi_snk2[t4]);
	  temp += trace(mf_pi_src2[t1], mf_pi_snk2[t4]) * trace(mf_pi_src1[t2], mf_pi_snk1[t3]);
	  temp *= ScalarComplexType(0.125);
	  into[1](t2, tdis) += temp;

	  MesonFieldType prod3, prod4, prod5;
	  ScalarComplexType temp2(0,0);
	  mult(prod3, mf_pi_src1[t2], mf_pi_src2[t1], NODE_LOCAL);
	  mult(prod4, mf_pi_snk1[t3], mf_pi_snk2[t4], NODE_LOCAL);
	  temp2 += ScalarComplexType(0.25) * trace(prod3, prod4);
	  mult(prod5, mf_pi_snk2[t4], mf_pi_snk1[t3], NODE_LOCAL);
	  temp2 += ScalarComplexType(0.25) * trace(prod3, prod5);
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
};
CPS_END_NAMESPACE

#endif
