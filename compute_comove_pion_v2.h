#ifndef _COMPUTE_COMOVE_PION_V2_H
#define _COMPUTE_COMOVE_PION_V2_H


#include<alg/a2a/mesonfield_computemany.h>
#include<alg/a2a/inner_product.h>
#include<alg/a2a/mf_momcontainer.h>
#include<alg/a2a/compute_pion.h>

CPS_START_NAMESPACE

template<typename mf_Policies>
class compute_comove_pion_v2
{
  public:
      typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
      template<typename PionMomentumPolicy>
      static void compute_v2(fMatrix<ScalarComplexType> &into, MesonFieldMomentumContainer<mf_Policies> &mf_ll_con, const PionMomentumPolicy &pion_mom, const int pidx1, const int pidx2, const int traj, const std::string &work_dir, bool ground=true)
      {
	int Lt = GJP.Tnodes()*GJP.TnodeSites();  
	into.resize(Lt,Lt);

	ThreeMomentum p_pi_src = pion_mom.getMesonMomentum(pidx1);
	ThreeMomentum p_pi_snk = pion_mom.getMesonMomentum(pidx2);

	typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MesonFieldType;
	std::vector<MesonFieldType> mf_pi1(Lt);
	std::vector<MesonFieldType> mf_pi2(Lt);

	std::ostringstream os; 
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

	if(!UniqueID()){ printf("Starting trace\n");  fflush(stdout); }
	trace(into,mf_pi2,mf_pi1);

	into *= ScalarComplexType(0.5,0);
	sync();

	rearrangeTsrcTsep(into); //rearrange temporal ordering
	if(!UniqueID()){ printf("Finished trace\n");  fflush(stdout); }

      }
};
CPS_END_NAMESPACE

#endif
