#ifndef _COMPUTE_COMOVE_PIPI_H
#define _COMPUTE_COMOVE_PIPI_H

#include<alg/a2a/mesonfield.h>
#include<alg/a2a/mesonfield_computemany.h>
#include<alg/a2a/inner_product.h>
#include<alg/a2a/mf_momcontainer.h>
#include<alg/a2a/compute_pion.h>
#include<alg/a2a/mf_productstore.h>
#include<alg/a2a/required_momenta.h>

#include<iostream>

using std::cout;
CPS_START_NAMESPACE
#define NODE_LOCAL true

template<typename mf_Policies>
class compute_comove_pipi
{
  public:
    typedef typename mf_Policies::ScalarComplexType ScalarComplexType;

    //tsep is b/w two pion
    //here psrc and psnk is defined so that tsnk > tsrc, into is evaluated at tsnk with tsep
    template<typename PionMomentumPolicy>
    static void compute_Vdis_v2(fVector<ScalarComplexType> &into, const PionMomentumPolicy &pion_mom, const int tsep, const int psrc, const int psnk, const int traj, const std::string &work_dir)
    {
      ThreeMomentum p_src = pion_mom.getMesonMomentum(psrc);
      ThreeMomentum p_snk = pion_mom.getMesonMomentum(psnk);
      int Lt = GJP.Tnodes()*GJP.TnodeSites();
      into.resize(Lt);
      into.zero();
      
      typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MesonFieldType;
      std::vector<MesonFieldType> mf_pi_src(Lt);
      std::vector<MesonFieldType> mf_pi_snk(Lt);

      std::ostringstream os;
      os << work_dir << "/traj_" << traj << "_pion_mf_mom" << p_src.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
      MesonFieldType::read(os.str(),mf_pi_src);
      os.str("");
      os << work_dir << "/traj_" << traj << "_pion_mf_mom" << p_snk.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
      MesonFieldType::read(os.str(),mf_pi_snk);
      os.str("");

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
    }


    template<typename PionMomentumPolicy>
    static void compute_v2(fMatrix<ScalarComplexType> &intoC, fMatrix<ScalarComplexType> &intoD, fMatrix<ScalarComplexType> &intoR, const PionMomentumPolicy &pion_mom, const int tsep, const int tstep_src, const int psrc1, const int psrc2, const int psnk1, const int psnk2, const int traj, const std::string &work_dir)
    {
      ThreeMomentum p_src1 = pion_mom.getMesonMomentum(psrc1);
      ThreeMomentum p_src2 = pion_mom.getMesonMomentum(psrc2);
      ThreeMomentum p_snk1 = pion_mom.getMesonMomentum(psnk1);
      ThreeMomentum p_snk2 = pion_mom.getMesonMomentum(psnk2);
#ifndef TIANLE_MOVING_PIPI
      assert(p_src1 + p_src2 == -(p_snk1 + p_snk2));
#endif
      int Lt = GJP.Tnodes()*GJP.TnodeSites();
      if(Lt % tstep_src != 0) {assert( 0 && "Lt has to be an integer times tstep_src");}

      intoC.resize(Lt,Lt);
      intoD.resize(Lt,Lt);
      intoR.resize(Lt,Lt);
      intoC.zero();
      intoD.zero();
      intoR.zero();

      typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MesonFieldType;
      std::vector<MesonFieldType> mf_pi_src1(Lt);
      std::vector<MesonFieldType> mf_pi_src2(Lt);
      std::vector<MesonFieldType> mf_pi_snk1(Lt);
      std::vector<MesonFieldType> mf_pi_snk2(Lt);

      std::ostringstream os;
      os << work_dir << "/traj_" << traj << "_pion_mf_mom" << p_src1.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
      MesonFieldType::read(os.str(),mf_pi_src1);
      os.str("");
      os << work_dir << "/traj_" << traj << "_pion_mf_mom" << p_snk1.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
      MesonFieldType::read(os.str(),mf_pi_snk1);
      os.str("");
      os << work_dir << "/traj_" << traj << "_pion_mf_mom" << p_src2.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
      MesonFieldType::read(os.str(),mf_pi_src2);
      os.str("");
      os << work_dir << "/traj_" << traj << "_pion_mf_mom" << p_snk2.file_str() << "_hyd" << "1s" << "_rad" << "2" << ".dat";
      MesonFieldType::read(os.str(),mf_pi_snk2);
      os.str("");

      int totwork = Lt * Lt / tstep_src;
      int start_pt, workpernode;
      bool do_work;
      //t4=tsnk+tsep,t3=tsnk,t2=tsrc,t1=tsrc-tsep,tdis=t3-t2
      int t1=0,t2=0,t3=0,t4=0;
      int tdis=0;
      getNodeWork(totwork,workpernode,start_pt,do_work);
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
	  intoC(t2, tdis) += ScalarComplexType(0.5) * trace(prod1, prod2); 
	    
          ScalarComplexType temp(0,0);
	  temp += trace(mf_pi_src2[t1], mf_pi_snk1[t3]) * trace(mf_pi_src1[t2], mf_pi_snk2[t4]);
	  temp += trace(mf_pi_src2[t1], mf_pi_snk2[t4]) * trace(mf_pi_src1[t2], mf_pi_snk1[t3]);
	  temp *= ScalarComplexType(0.125);
	  intoD(t2, tdis) += temp;
	    
          MesonFieldType prod3, prod4, prod5;
	  ScalarComplexType temp2(0,0);
	  mult(prod3, mf_pi_src1[t2], mf_pi_src2[t1], NODE_LOCAL);
	  mult(prod4, mf_pi_snk1[t3], mf_pi_snk2[t4], NODE_LOCAL);
	  temp2 += ScalarComplexType(0.25) * trace(prod3, prod4);
	  mult(prod5, mf_pi_snk2[t4], mf_pi_snk1[t3], NODE_LOCAL);
	  temp2 += ScalarComplexType(0.25) * trace(prod3, prod5);
	  intoR(t2, tdis) += temp2;
        }
      }
      sync();
      intoC.nodeSum();
      intoD.nodeSum();
      intoR.nodeSum();
      if(!UniqueID()){ printf("Finished trace\n");  fflush(stdout); }
    }
};
CPS_END_NAMESPACE

#endif
