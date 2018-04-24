#ifndef _MAIN_H
#define _MAIN_H

//Header for main program

using namespace cps;

//Setup the A2A policy
#ifdef USE_DESTRUCTIVE_FFT

#ifdef A2A_PREC_DOUBLE
typedef A2ApoliciesDoubleManualAlloc A2Apolicies;
#elif defined(A2A_PREC_SINGLE)
typedef A2ApoliciesSingleManualAlloc A2Apolicies;
#elif defined(A2A_PREC_SIMD_DOUBLE)
typedef A2ApoliciesSIMDdoubleManualAlloc A2Apolicies;
#elif defined(A2A_PREC_SIMD_SINGLE)
typedef A2ApoliciesSIMDsingleManualAlloc A2Apolicies;
#else
#error "Must provide an A2A precision"
#endif

#else

#ifdef A2A_PREC_DOUBLE
typedef A2ApoliciesDoubleAutoAlloc A2Apolicies;
#elif defined(A2A_PREC_SINGLE)
typedef A2ApoliciesSingleAutoAlloc A2Apolicies;
#elif defined(A2A_PREC_SIMD_DOUBLE)
typedef A2ApoliciesSIMDdoubleAutoAlloc A2Apolicies;
#elif defined(A2A_PREC_SIMD_SINGLE)
typedef A2ApoliciesSIMDsingleAutoAlloc A2Apolicies;
#else
#error "Must provide an A2A precision"
#endif

#endif

//Defines for Grid/BFM wrapping
#ifdef USE_GRID_LANCZOS
  typedef GridLanczosWrapper<A2Apolicies> LanczosWrapper;
  typedef typename A2Apolicies::FgridGFclass LanczosLattice;
# define LANCZOS_LATARGS params.jp
# define LANCZOS_LATMARK isGridtype
# define LANCZOS_EXTRA_ARG *lanczos_lat
# define COMPUTE_EVECS_EXTRA_ARG_PASS NULL
# define COMPUTE_EVECS_EXTRA_ARG_GRAB void*
#else //USE_BFM_LANCZOS
  typedef BFMLanczosWrapper LanczosWrapper;
  typedef GwilsonFdwf LanczosLattice;
# define LANCZOS_LATARGS bfm_solvers
# define LANCZOS_LATMARK isBFMtype
# define LANCZOS_EXTRA_ARG bfm_solvers
# define COMPUTE_EVECS_EXTRA_ARG_PASS bfm_solvers
# define COMPUTE_EVECS_EXTRA_ARG_GRAB BFMsolvers &bfm_solvers
#endif

#ifdef USE_GRID_A2A
  typedef A2Apolicies::FgridGFclass A2ALattice;
# define A2A_LATARGS params.jp
# define A2A_LATMARK isGridtype
#else
  typedef GwilsonFdwf A2ALattice;
# define A2A_LATARGS bfm_solvers
# define A2A_LATMARK isBFMtype
#endif



//Command line argument store/parse
struct CommandLineArgs{
  int nthreads;
  bool randomize_vw; //rather than doing the Lanczos and inverting the propagators, etc, just use random vectors for V and W
  bool randomize_evecs; //skip Lanczos and just use random evecs for testing.
  bool tune_lanczos_light; //just run the light lanczos on first config then exit
  bool tune_lanczos_heavy; //just run the heavy lanczos on first config then exit
  bool skip_gauge_fix;
  bool double_latt; //most ancient 8^4 quenched lattices stored both U and U*. Enable this to read those configs
  bool mixed_solve; //do high mode inversions using mixed precision solves. Is disabled if we turn off the single-precision conversion of eigenvectors (because internal single-prec inversion needs singleprec eigenvectors)
  bool evecs_single_prec; //convert the eigenvectors to single precision to save memory
  bool do_kaon2pt;
  bool do_pion2pt;
  bool do_pion2ptv2;
  bool do_pion2ptv3;
  bool do_pion2ptcomovev2;
  bool do_comove_pipi;
  bool do_pipi;
  bool do_ktopipi;
  bool do_sigma;
  bool do_pipisigma;
  bool do_sigmapipi;
  bool do_fitting;
  bool do_pipiv2;
  bool do_1s;

  Float inner_cg_resid;
  Float *inner_cg_resid_p;

  bool do_split_job;
  int split_job_part;
  std::string checkpoint_dir; //directory the checkpointed data is stored in (could be scratch for example)
  
  CommandLineArgs(int argc, char **argv, int begin){
    nthreads = 1;
#if TARGET == BGQ
    nthreads = 64;
#endif
    randomize_vw = false;
    randomize_evecs = false;
    tune_lanczos_light = false; //just run the light lanczos on first config then exit
    tune_lanczos_heavy = false; //just run the heavy lanczos on first config then exit
    skip_gauge_fix = false;
    double_latt = false; //most ancient 8^4 quenched lattices stored both U and U*. Enable this to read those configs
    mixed_solve = true; //do high mode inversions using mixed precision solves. Is disabled if we turn off the single-precision conversion of eigenvectors (because internal single-prec inversion needs singleprec eigenvectors)
    evecs_single_prec = true; //convert the eigenvectors to single precision to save memory
    do_1s = true;
    do_kaon2pt = false;
    do_pion2pt = false;
    do_pion2ptv2 = false;
    do_pion2ptv3 = false;
    do_pion2ptcomovev2 = false;
    do_comove_pipi = false;
    do_pipi = false;
    do_pipiv2 = true;
    do_ktopipi = false;
    do_sigma = false;
    do_pipisigma = false;
    do_sigmapipi = false;
    do_fitting = false;
    inner_cg_resid;
    inner_cg_resid_p = NULL;

    do_split_job = false;
    
    parse(argc,argv,begin);
  }

  void parse(int argc, char **argv, int begin){
    if(!UniqueID()){ printf("Arguments:\n"); fflush(stdout); }
    for(int i=0;i<argc;i++){
      if(!UniqueID()){ printf("%d \"%s\"\n",i,argv[i]); fflush(stdout); }
    }
    
    const int ngrid_arg = 10;
    const std::string grid_args[ngrid_arg] = { "--debug-signals", "--dslash-generic", "--dslash-unroll", "--dslash-asm", "--shm", "--lebesgue", "--cacheblocking", "--comms-isend", "--comms-sendrecv", "--comms-overlap" };
    const int grid_args_skip[ngrid_arg] =    { 1                , 1                 , 1                , 1             , 2      , 1           , 2                , 1              , 1                 , 1 };

    int arg = begin;
    while(arg < argc){
      char* cmd = argv[arg];
      if( strncmp(cmd,"-nthread",8) == 0){
	if(arg == argc-1){ if(!UniqueID()){ printf("-nthread must be followed by a number!\n"); fflush(stdout); } exit(-1); }
	nthreads = strToAny<int>(argv[arg+1]);
	if(!UniqueID()){ printf("Setting number of threads to %d\n",nthreads); }
	arg+=2;
      }else if( strncmp(cmd,"-set_inner_resid",16) == 0){ //only for mixed CG
	if(arg == argc-1){ if(!UniqueID()){ printf("-set_inner_resid must be followed by a number!\n"); fflush(stdout); } exit(-1); }
	inner_cg_resid = strToAny<Float>(argv[arg+1]);
	inner_cg_resid_p = &inner_cg_resid;
	if(!UniqueID()){ printf("Setting inner CG initial residual to %g\n",inner_cg_resid); }
#ifdef USE_BFM_A2A
	ERR.General("","main","Changing initial inner CG residual not implemented for BFM version\n");
#endif      
	arg+=2;      
      }else if( strncmp(cmd,"-randomize_vw",15) == 0){
	randomize_vw = true;
	if(!UniqueID()){ printf("Using random vectors for V and W, skipping Lanczos and inversion stages\n"); fflush(stdout); }
	arg++;
      }else if( strncmp(cmd,"-randomize_evecs",15) == 0){
	randomize_evecs = true;
	if(!UniqueID()){ printf("Using random eigenvectors\n"); fflush(stdout); }
	arg++;      
      }else if( strncmp(cmd,"-tune_lanczos_light",15) == 0){
	tune_lanczos_light = true;
	if(!UniqueID()){ printf("Just tuning light lanczos on first config\n"); fflush(stdout); }
	arg++;
      }else if( strncmp(cmd,"-tune_lanczos_heavy",15) == 0){
	tune_lanczos_heavy = true;
	if(!UniqueID()){ printf("Just tuning heavy lanczos on first config\n"); fflush(stdout); }
	arg++;
      }else if( strncmp(cmd,"-double_latt",15) == 0){
	double_latt = true;
	if(!UniqueID()){ printf("Loading doubled lattices\n"); fflush(stdout); }
	arg++;
      }else if( strncmp(cmd,"-skip_gauge_fix",20) == 0){
	skip_gauge_fix = true;
	if(!UniqueID()){ printf("Skipping gauge fixing\n"); fflush(stdout); }
	arg++;
      }else if( strncmp(cmd,"-disable_evec_singleprec_convert",30) == 0){
	evecs_single_prec = false;
	mixed_solve = false;
	if(!UniqueID()){ printf("Disabling single precision conversion of evecs\n"); fflush(stdout); }
	arg++;
      }else if( strncmp(cmd,"-disable_mixed_prec_CG",30) == 0){
	mixed_solve = false;
	if(!UniqueID()){ printf("Disabling mixed-precision CG\n"); fflush(stdout); }
	arg++;
      }else if( strncmp(cmd,"-mf_outerblocking",15) == 0){
	int* b[3] = { &BlockedMesonFieldArgs::bi, &BlockedMesonFieldArgs::bj, &BlockedMesonFieldArgs::bp };
	for(int a=0;a<3;a++) *b[a] = strToAny<int>(argv[arg+1+a]);
	arg+=4;
      }else if( strncmp(cmd,"-mf_innerblocking",15) == 0){
	int* b[3] = { &BlockedMesonFieldArgs::bii, &BlockedMesonFieldArgs::bjj, &BlockedMesonFieldArgs::bpp };
	for(int a=0;a<3;a++) *b[a] = strToAny<int>(argv[arg+1+a]);
	arg+=4;
      }else if( strncmp(cmd,"-do_split_job",30) == 0){
	do_split_job = true;
	split_job_part = strToAny<int>(argv[arg+1]);
	checkpoint_dir = argv[arg+2];
	if(!UniqueID()) printf("Doing split job part %d with checkpoint directory %s\n",split_job_part,checkpoint_dir.c_str());
	arg+=3;       
      }else if( strncmp(cmd,"-skip_kaon2pt",30) == 0){
	do_kaon2pt = false;
	arg++;
      }else if( strncmp(cmd,"-skip_pion2pt",30) == 0){
	do_pion2pt = false;
	arg++;
      }else if( strncmp(cmd,"-skip_sigma",30) == 0){
	do_sigma = false;
	arg++;
      }else if( strncmp(cmd,"-skip_pipi",30) == 0){
	do_pipi = false;
	arg++;
      }else if( strncmp(cmd,"-skip_ktopipi",30) == 0){
	do_ktopipi = false;
	arg++;  
      }else{
	bool is_grid_arg = false;
	for(int i=0;i<ngrid_arg;i++){
	  if( std::string(cmd) == grid_args[i] ){
	    if(!UniqueID()){ printf("main.C: Ignoring Grid argument %s\n",cmd); fflush(stdout); }
	    arg += grid_args_skip[i];
	    is_grid_arg = true;
	    break;
	  }
	}
	if(!is_grid_arg){
	  if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
	  exit(-1);
	}
      }
    }
  }

};

//Store/read job parameters
struct Parameters{
  CommonArg common_arg;
  CommonArg common_arg2;
  DoArg do_arg;
  JobParams jp;
  MeasArg meas_arg;
  FixGaugeArg fix_gauge_arg;
  A2AArg a2a_arg;
  A2AArg a2a_arg_s;
  LancArg lanc_arg;
  LancArg lanc_arg_s;

  Parameters(const char* directory): common_arg("",""), common_arg2("",""){
    if(chdir(directory)!=0) ERR.General("Parameters","Parameters","Unable to switch to directory '%s'\n",directory);

    if(!do_arg.Decode("do_arg.vml","do_arg")){
      do_arg.Encode("do_arg.templ","do_arg");
      VRB.Result("Parameters","Parameters","Can't open do_arg.vml!\n");exit(1);
    }
    if(!jp.Decode("job_params.vml","job_params")){
      jp.Encode("job_params.templ","job_params");
      VRB.Result("Parameters","Parameters","Can't open job_params.vml!\n");exit(1);
    }
    if(!meas_arg.Decode("meas_arg.vml","meas_arg")){
      meas_arg.Encode("meas_arg.templ","meas_arg");
      std::cout<<"Can't open meas_arg!"<<std::endl;exit(1);
    }
    if(!a2a_arg.Decode("a2a_arg.vml","a2a_arg")){
      a2a_arg.Encode("a2a_arg.templ","a2a_arg");
      VRB.Result("Parameters","Parameters","Can't open a2a_arg.vml!\n");exit(1);
    }
    if(!a2a_arg_s.Decode("a2a_arg_s.vml","a2a_arg_s")){
      a2a_arg_s.Encode("a2a_arg_s.templ","a2a_arg_s");
      VRB.Result("Parameters","Parameters","Can't open a2a_arg_s.vml!\n");exit(1);
    }
    if(!lanc_arg.Decode("lanc_arg.vml","lanc_arg")){
      lanc_arg.Encode("lanc_arg.templ","lanc_arg");
      VRB.Result("Parameters","Parameters","Can't open lanc_arg.vml!\n");exit(1);
    }
    if(!lanc_arg_s.Decode("lanc_arg_s.vml","lanc_arg_s")){
      lanc_arg_s.Encode("lanc_arg_s.templ","lanc_arg_s");
      VRB.Result("Parameters","Parameters","Can't open lanc_arg_s.vml!\n");exit(1);
    }
    if(!fix_gauge_arg.Decode("fix_gauge_arg.vml","fix_gauge_arg")){
      fix_gauge_arg.Encode("fix_gauge_arg.templ","fix_gauge_arg");
      VRB.Result("Parameters","Parameters","Can't open fix_gauge_arg.vml!\n");exit(1);
    }

    common_arg.set_filename(meas_arg.WorkDirectory);
  }

};


void setupJob(int argc, char **argv, const Parameters &params, const CommandLineArgs &cmdline){
#ifdef NODE_DISTRIBUTE_MESONFIELDS
  if(!UniqueID()) printf("Using node distribution of meson fields\n");
#endif
#ifdef MEMTEST_MODE
  if(!UniqueID()) printf("Running in MEMTEST MODE (so don't expect useful results)\n");
#endif
  
#ifdef A2A_LANCZOS_SINGLE
  if(!cmdline.evecs_single_prec) ERR.General("",fname,"Must use single-prec eigenvectors when doing Lanczos in single precision\n");
#endif

  GJP.Initialize(params.do_arg);
  LRG.Initialize();

#if defined(USE_GRID_A2A) || defined(USE_GRID_LANCZOS)
  if(GJP.Gparity()){
#ifndef USE_GRID_GPARITY
    ERR.General("","","Must compile main program with flag USE_GRID_GPARITY to enable G-parity\n");
#endif
  }else{
#ifdef USE_GRID_GPARITY
    ERR.General("","","Must compile main program with flag USE_GRID_GPARITY off to disable G-parity\n");
#endif
  }      
#endif
  
  if(cmdline.double_latt) SerialIO::dbl_latt_storemode = true;

  if(!cmdline.tune_lanczos_light && !cmdline.tune_lanczos_heavy){ 
    assert(params.a2a_arg.nl <= params.lanc_arg.N_true_get);
    assert(params.a2a_arg_s.nl <= params.lanc_arg_s.N_true_get);
  }
#ifdef USE_BFM
  cps_qdp_init(&argc,&argv);
  //Chroma::initialize(&argc,&argv);
#endif
  omp_set_num_threads(cmdline.nthreads);

  if(!UniqueID()) printf("Initial memory post-initialize:\n");
  printMem();
}


void read_pion_mesonfield(MesonFieldMomentumContainer<A2Apolicies> &store, const StandardPionMomentaPolicy &pion_mom, 
                          const int traj, const Parameters &params, const std::string &work_dir)
{
  if(!UniqueID()) printf("Start doing read_pion_mesonfield\n");
  double time = -dclock();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  std::vector<A2AmesonField<A2Apolicies,A2AvectorWfftw,A2AvectorVfftw> > tmp(Lt);
  for(int p=0;p<pion_mom.nMom();p++)
  {   
    std::ostringstream os; 
    os << work_dir << "/traj_" << traj << "_pion_mf_mom" << pion_mom.getMesonMomentum(p).file_str() << "_hyd1s_rad" << params.jp.pion_rad << ".dat";
    A2AmesonField<A2Apolicies,A2AvectorWfftw,A2AvectorVfftw>::read(os.str(), tmp);    
    std::vector<A2AmesonField<A2Apolicies,A2AvectorWfftw,A2AvectorVfftw> > &stored = store.copyAdd(pion_mom.getMesonMomentum(p), tmp);
#ifdef NODE_DISTRIBUTE_MESONFIELDS
    nodeDistributeMany(1,&stored);
#endif   
  }
  time += dclock();
  print_time("main","read_pion_mesonfield",time);
  if(!UniqueID()) printf("Memory after read_pion_mesonfield:\n");
  printMem();
}

void read_2s_pion_mesonfield(MesonFieldMomentumContainer<A2Apolicies> &store, const StandardPionMomentaPolicy &pion_mom, 
                          const int traj, const Parameters &params, const std::string &work_dir)
{
  if(!UniqueID()) printf("Start doing read_2s_pion_mesonfield\n");
  double time = -dclock();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  std::vector<A2AmesonField<A2Apolicies,A2AvectorWfftw,A2AvectorVfftw> > tmp(Lt);
  for(int p=0;p<pion_mom.nMom();p++)
  {   
    std::ostringstream os; 
    os << work_dir << "/traj_" << traj << "_pion_mf_mom" << pion_mom.getMesonMomentum(p).file_str() << "_hyd2s_rad" << params.jp.pion_rad << ".dat";
    A2AmesonField<A2Apolicies,A2AvectorWfftw,A2AvectorVfftw>::read(os.str(), tmp);    
    std::vector<A2AmesonField<A2Apolicies,A2AvectorWfftw,A2AvectorVfftw> > &stored = store.copyAdd(pion_mom.getMesonMomentum(p), tmp);
#ifdef NODE_DISTRIBUTE_MESONFIELDS
    nodeDistributeMany(1,&stored);
#endif   
  }
  time += dclock();
  print_time("main","read_2s_pion_mesonfield",time);
  if(!UniqueID()) printf("Memory after read_2s_pion_mesonfield:\n");
  printMem();
}

void read_sigma_mesonfield(MesonFieldMomentumPairContainer<A2Apolicies> &store, const StationarySigmaMomentaPolicy &sigma_mom,
                           const int traj, const Parameters &params, const std::string &work_dir)
{
  if(!UniqueID()) printf("Start doing read_sigma_mesonfield\n");
  double time = -dclock();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  std::vector<A2AmesonField<A2Apolicies,A2AvectorWfftw,A2AvectorVfftw> > tmp(Lt);
  for(int p=0;p<sigma_mom.nMom();p++)
  {
    std::ostringstream os;
    os << work_dir << "/traj_" << traj << "_sigma_mfwv_mom" << sigma_mom.getWdagMom(p).file_str() << "_plus" << sigma_mom.getVmom(p).file_str() << "_hyd1s_rad" << params.jp.pion_rad << ".dat";
    A2AmesonField<A2Apolicies,A2AvectorWfftw,A2AvectorVfftw>::read(os.str(), tmp);
    std::vector<A2AmesonField<A2Apolicies,A2AvectorWfftw,A2AvectorVfftw> > &stored = store.copyAdd(sigma_mom.getWdagMom(p), sigma_mom.getVmom(p), tmp);
#ifdef NODE_DISTRIBUTE_MESONFIELDS
    nodeDistributeMany(1,&stored);
#endif
  }
  time += dclock();
  print_time("main","read_sigma_mesonfield",time);
  if(!UniqueID()) printf("Memory after read_sigma_mesonfield:\n");
  printMem();
}

void computepion2ptv2(MesonFieldMomentumContainer<A2Apolicies> &mf_ll_con, const StandardPionMomentaPolicy &pion_mom, const int conf, const Parameters &params, bool do_1s)
{
  const int nmom = pion_mom.nMom();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  std::cout<<"Lt = "<<Lt<<'\n';
  if(!UniqueID()) printf("Computing pion 2pt function\n");
  double time = -dclock();
  for(int p=0;p<nmom;p+=2) 
  {
    if(!UniqueID()) printf("Starting pidx %d\n",p);
    fMatrix<typename A2Apolicies::ScalarComplexType> pion(Lt,Lt);    
    computepion_v2<A2Apolicies>::compute_v2(pion, mf_ll_con, pion_mom, p, conf, params.meas_arg.WorkDirectory, do_1s);
#define DAIQIAN_PION_PHASE_CONVENTION

    std::ostringstream os; 
    os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_pioncorr_mom";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
    os << pion_mom.getMesonMomentum(p).file_str();  
#else
    os << (-pion_mom.getMesonMomentum(p)).file_str();
#endif
    if(!do_1s)
      os << "_2s";
    os << "_v2";
    pion.write(os.str());
#ifdef WRITE_HEX_OUTPUT
    os << ".hexfloat";
    pion.write(os.str(),true);
#endif
  }
  time += dclock();
  print_time("main","Pion 2pt function v2",time);

  if(!UniqueID()) printf("Memory after pion 2pt function v2 computation:\n");
  printMem();
}

void compute_comove_pion2ptv2(MesonFieldMomentumContainer<A2Apolicies> &mf_ll_con, const StandardPionMomentaPolicy &pion_mom, const int conf, const Parameters &params, bool do_1s)
{
  const int nmom = pion_mom.nMom();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  std::cout<<"Lt = "<<Lt<<'\n';
  if(!UniqueID()) printf("Computing pion 2pt function\n");
  double time = -dclock();
  for(int p1=0;p1<nmom;p1++)
  {
    for(int p2=0; p2<nmom; p2++)
    {
      if(!UniqueID()) printf("Starting pidx1 %d  and pidx2 %d\n",p1,p2);
      fMatrix<typename A2Apolicies::ScalarComplexType> pion(Lt,Lt);
      compute_comove_pion_v2<A2Apolicies>::compute_v2(pion, mf_ll_con, pion_mom, p1, p2, conf, params.meas_arg.WorkDirectory, do_1s);
#define DAIQIAN_PION_PHASE_CONVENTION

      std::ostringstream os;
      os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_comove_pioncorr_mom1";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
      os << pion_mom.getMesonMomentum(p1).file_str() ;
#else
      os << (-pion_mom.getMesonMomentum(p1)).file_str();
#endif
      os << "_comove_pioncorr_mom2";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
      os << pion_mom.getMesonMomentum(p2).file_str() ;
#else
      os << (-pion_mom.getMesonMomentum(p2)).file_str();
#endif
      if(!do_1s)
        os << "_2s";
      os << "_v2";
      pion.write(os.str());
#ifdef WRITE_HEX_OUTPUT
      os << ".hexfloat";
      pion.write(os.str(),true);
#endif
    }
  }

  time += dclock();
  print_time("main","comove_pion 2pt function v2",time);

  if(!UniqueID()) printf("Memory after comove_pion 2pt function v2 computation:\n");
  printMem();

}


void computesigma2ptv2(MesonFieldMomentumContainer<A2Apolicies> &mf_ll_con, const StationarySigmaMomentaPolicy &sigma_mom, const int conf, const Parameters &params, bool do_1s)
{
  const int nmom = sigma_mom.nMom();
//  const int nmom=8; 
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  
  if(!UniqueID()) printf("Computing sigma 2pt function v2\n");
  double time = -dclock();
  for(int psrc=0;psrc<nmom;psrc++) 
  {
    for(int psnk=0; psnk<nmom; psnk++)
    {
      if(!UniqueID()) printf("Starting psrc and psck with %d and %d\n",psrc,psnk);
      fMatrix<typename A2Apolicies::ScalarComplexType> sigma(Lt,Lt);    
      fVector<typename A2Apolicies::ScalarComplexType> sigmaself(Lt);
      computesigma_v2<A2Apolicies>::compute_v2(sigma, sigmaself, mf_ll_con, sigma_mom, psrc, psnk, conf, params.meas_arg.WorkDirectory, do_1s);
#define DAIQIAN_PION_PHASE_CONVENTION

      std::ostringstream os,os1; 
      os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_sigmacorr_mom";
      if(psnk==0)
        os1 << params.meas_arg.WorkDirectory << "/traj_" << conf << "_sigmaself_mom";

#ifndef DAIQIAN_PION_PHASE_CONVENTION
      os << "psrc" << (-sigma_mom.getWmom(psrc)).file_str() << "psnk" << (-sigma_mom.getWmom(psnk)).file_str();  
      if(psnk==0)
        os1 << (-sigma_mom.getWmom(psrc)).file_str();
#else
      os << "psrc" << sigma_mom.getWmom(psrc).file_str() << "psnk" << sigma_mom.getWmom(psnk).file_str();
      if(psnk==0)
        os1 << sigma_mom.getWmom(psrc).file_str();
#endif
      if(!do_1s)
        os << "_2s";
      os << "_v2";
      sigma.write(os.str());
      if(psnk==0)  
      {
        if(!do_1s)
          os1 << "_2s";
        os1 << "_v2";
        sigmaself.write(os1.str());
      }
#ifdef WRITE_HEX_OUTPUT
      os << ".hexfloat";
      sigma.write(os.str(),true);
      if(psnk==0)
      {
        os1 << ".hexfloat";
        sigmaself.write(os1.str(),true);
      }
#endif
    }
  }
  time += dclock();
  print_time("main","Sigma 2pt function v2",time);

  if(!UniqueID()) printf("Memory after sigma 2pt function v2 computation:\n");
  printMem();
}

void computePion2pt_v3(MesonFieldMomentumContainer<A2Apolicies> &mf_1s_con, MesonFieldMomentumContainer<A2Apolicies> &mf_2s_con, 
                       const StandardPionMomentaPolicy &pion_mom, const int conf, const Parameters &params, std::string src_type, std::string snk_type)
{
  const int nmom = pion_mom.nMom();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  if(!UniqueID()) printf("Doing computePion2pt_v3 with src_type: %s and snk_type: %s\n", src_type.c_str(), snk_type.c_str());
  double time = -dclock();
  for(int psrc = 0; psrc < nmom; psrc += 2)
  {
    if(!UniqueID()) printf("Starting pidx %d\n",psrc);
    fMatrix<typename A2Apolicies::ScalarComplexType> pion(Lt,Lt);
    computepion_v3<A2Apolicies>::compute(pion,mf_1s_con,mf_2s_con,pion_mom,psrc,src_type,snk_type);
#define DAIQIAN_PION_PHASE_CONVENTION
    std::ostringstream os; 
    os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_pioncorr_mom";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
    os << pion_mom.getMesonMomentum(psrc).file_str();  
#else
    os << (-pion_mom.getMesonMomentum(psrc)).file_str();
#endif
    os << "_src_" << src_type << "_snk_" << snk_type;
    os << "_v3";
    pion.write(os.str());
  }
  time += dclock();
  print_time("main","Pion 2pt function v3",time);
  if(!UniqueID()) printf("Memory after Pion 2pt function v3 computation:\n");
  printMem();
}

void computesigma2ptv3(MesonFieldMomentumPairContainer<A2Apolicies> &mf_sigma_con, const StationarySigmaMomentaPolicy &sigma_mom, const int conf, const Parameters &params)
{
  const int nmom = sigma_mom.nMom();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();

  if(!UniqueID()) printf("Computing sigma 2pt function v3\n");
  double time = -dclock();
  for(int psrc=0;psrc<nmom;psrc++)
  {
    for(int psnk=0; psnk<nmom; psnk++)
    {
      if(!UniqueID()) printf("Starting psrc and psck with %d and %d\n",psrc,psnk);
      fMatrix<typename A2Apolicies::ScalarComplexType> sigma(Lt,Lt);
      fVector<typename A2Apolicies::ScalarComplexType> sigmaself(Lt);
      computesigma_v3<A2Apolicies>::compute(sigma,sigmaself,mf_sigma_con,sigma_mom,psrc,psnk);
#define DAIQIAN_PION_PHASE_CONVENTION
      std::ostringstream os,os1;
      os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_sigmacorr_mom";
      if(psnk==0)
        os1 << params.meas_arg.WorkDirectory << "/traj_" << conf << "_sigmaself_mom";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
      os << "psrc" << (-sigma_mom.getWmom(psrc)).file_str() << "psnk" << (-sigma_mom.getWmom(psnk)).file_str();
      if(psnk==0)
        os1 << (-sigma_mom.getWmom(psrc)).file_str();
#else
      os << "psrc" << sigma_mom.getWmom(psrc).file_str() << "psnk" << sigma_mom.getWmom(psnk).file_str();
      if(psnk==0)
        os1 << sigma_mom.getWmom(psrc).file_str();
#endif
      os << "_v2";
      sigma.write(os.str());
      if(psnk==0)
      {
        os1 << "_v2";
        sigmaself.write(os1.str());
      }
    }
  }
  time += dclock();
  print_time("main","Sigma 2pt function v3",time);

  if(!UniqueID()) printf("Memory after sigma 2pt function v2 computation:\n");
  printMem();
}

void computepipisigmav2(MesonFieldMomentumContainer<A2Apolicies> &mf_ll_con, const StandardPionMomentaPolicy &pion_mom, const StationarySigmaMomentaPolicy &sigma_mom, const int conf, const Parameters &params, bool do_1s)
{
  const int nsigmamom = sigma_mom.nMom();
  const int npimom = pion_mom.nMom();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  if(!UniqueID()) printf("Computing pipitosigma v2\n");
  double time = -dclock();
  for(int psigma=0; psigma < nsigmamom; psigma++) 
  {
    for(int ppi=0; ppi < npimom; ppi++)
    {
      if(!UniqueID()){ printf("Doing pipitosigma v2, ppion=%d psigma=%d\n",ppi,psigma); fflush(stdout); }
      
      fMatrix<typename A2Apolicies::ScalarComplexType> pipisigma(Lt,Lt);
      MesonFieldProductStore<A2Apolicies> products;
      computepipisigma_v2<A2Apolicies>::compute_v2(pipisigma, mf_ll_con, pion_mom, sigma_mom, products, params.jp.pipi_separation, ppi, psigma, conf, params.meas_arg.WorkDirectory, do_1s);
#define DAIQIAN_PION_PHASE_CONVENTION
      std::ostringstream os; 
      os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_pipitosigma_sigmawdagmom";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
      os << (-sigma_mom.getWmom(psigma)).file_str() << "_pionmom" << pion_mom.getMesonMomentum(ppi).file_str();  
#else
      os << sigma_mom.getWmom(psigma).file_str() << "_pionmom" << (-pion_mom.getMesonMomentum(ppi)).file_str();
#endif
      if(!do_1s)
        os << "_2s";
      os << "_v2";
      pipisigma.write(os.str());
#ifdef WRITE_HEX_OUTPUT
      os << ".hexfloat";
      pipisigma.write(os.str(),true);
#endif
    }
  }
  time += dclock();
  print_time("main","pipitosigma v2",time);

  if(!UniqueID()) printf("Memory after pipitosigma v2 computation:\n");
  printMem();
}

void computepipitosigmav3(MesonFieldMomentumContainer<A2Apolicies> &mf_pion_con, MesonFieldMomentumPairContainer<A2Apolicies> &mf_sigma_con,
                          const StandardPionMomentaPolicy &pion_mom, const StationarySigmaMomentaPolicy &sigma_mom, const int conf, const Parameters &params)
{
  const int nsigmamom = sigma_mom.nMom();
  const int npimom = pion_mom.nMom();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  if(!UniqueID()) printf("Computing pipitosigma v3\n");
  double time = -dclock();
  for(int psigma=0; psigma < nsigmamom; psigma++)
  {
    for(int ppi=0; ppi < npimom; ppi++)
    {
      if(!UniqueID()){ printf("Doing pipitosigma v3, p_pion=%d p_sigma=%d\n",ppi,psigma); fflush(stdout); }
      fMatrix<typename A2Apolicies::ScalarComplexType> pipisigma(Lt,Lt);
      computepipisigma_v3<A2Apolicies>::compute(pipisigma,mf_pion_con,mf_sigma_con,pion_mom,sigma_mom,ppi,psigma,params.jp.pipi_separation);
#define DAIQIAN_PION_PHASE_CONVENTION
      std::ostringstream os;
      os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_pipitosigma_sigmawdagmom";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
      os << (-sigma_mom.getWmom(psigma)).file_str() << "_pionmom" << pion_mom.getMesonMomentum(ppi).file_str();
#else
      os << sigma_mom.getWmom(psigma).file_str() << "_pionmom" << (-pion_mom.getMesonMomentum(ppi)).file_str();
#endif
      os << "_v2";
      pipisigma.write(os.str());
    }
  }
  time += dclock();
  print_time("main","pipitosigma v3",time);

  if(!UniqueID()) printf("Memory after pipitosigma v3 computation:\n");
  printMem();
}

void computesigmapipiv2(MesonFieldMomentumContainer<A2Apolicies> &mf_ll_con, const StandardPionMomentaPolicy &pion_mom, const StationarySigmaMomentaPolicy &sigma_mom, const int conf, const Parameters &params)
{
  const int nsigmamom = sigma_mom.nMom();
  const int npimom = pion_mom.nMom();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  if(!UniqueID()) printf("Computing sigmatopipi v2\n");
  double time = -dclock();
  for(int psigma=0; psigma < nsigmamom; psigma++) 
  {
    for(int ppi=0; ppi < npimom; ppi++)
    {
      if(!UniqueID()){ printf("Doing sigmatopipi v2, ppion=%d psigma=%d\n",ppi,psigma); fflush(stdout); }
      
      fMatrix<typename A2Apolicies::ScalarComplexType> sigmapipi(Lt,Lt);
      MesonFieldProductStore<A2Apolicies> products;
      computesigmapipi_v2<A2Apolicies>::compute_v2(sigmapipi, mf_ll_con, pion_mom, sigma_mom, products, params.jp.pipi_separation, ppi, psigma, conf, params.meas_arg.WorkDirectory);
#define DAIQIAN_PION_PHASE_CONVENTION
      std::ostringstream os; 
      os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_sigmatopipi_sigmawdagmom";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
      os << (-sigma_mom.getWmom(psigma)).file_str() << "_pionmom" << pion_mom.getMesonMomentum(ppi).file_str();  
#else
      os << sigma_mom.getWmom(psigma).file_str() << "_pionmom" << (-pion_mom.getMesonMomentum(ppi)).file_str();
#endif
      os << "_v2";
      sigmapipi.write(os.str());
#ifdef WRITE_HEX_OUTPUT
      os << ".hexfloat";
      sigmapipi.write(os.str(),true);
#endif
    }
  }
  time += dclock();
  print_time("main","sigmatopipi v2",time);

  if(!UniqueID()) printf("Memory after sigmatopipi v2 computation:\n");
  printMem();
}

void computepipi2ptv2(MesonFieldMomentumContainer<A2Apolicies> &mf_1s_con, MesonFieldMomentumContainer<A2Apolicies> &mf_2s_con, 
                      const StandardPionMomentaPolicy &pion_mom, const int conf, const Parameters &params, std::string src_type, std::string snk_type)
{
  const int nmom = pion_mom.nMom();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  if(!UniqueID()) printf("Computing pipi2ptv2 function with src_type %s and snk_type %s\n", src_type.c_str(), snk_type.c_str());
  double time = -dclock();

  for(int psrc=0; psrc < nmom; psrc++)
  {
    ThreeMomentum p_pi1_src = pion_mom.getMesonMomentum(psrc);
    for(int psnk=0; psnk < nmom; psnk++)
    {
      ThreeMomentum p_pi1_snk = pion_mom.getMesonMomentum(psnk);
      std::array<fMatrix<typename A2Apolicies::ScalarComplexType>,3> pipi;
      if(!UniqueID()){ printf("Doing pipi2ptv2, psrc=%d psnk=%d\n",psrc,psnk); fflush(stdout); }
      computepipi_v2<A2Apolicies>::compute_v2(pipi, mf_1s_con, mf_2s_con, pion_mom, psrc, psnk, params.jp.pipi_separation, params.jp.tstep_pipi, src_type, snk_type);
      char diag[3] = {'C', 'D', 'R'};
      for(int i=0; i<3; i++)
      {
        std::ostringstream os;
        os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_Figure" << diag[i] << "_sep" << params.jp.pipi_separation;
#ifndef DAIQIAN_PION_PHASE_CONVENTION
        os << "_mom" << p_pi1_src.file_str(2) << "_mom" << p_pi1_snk.file_str(2);
#else
        os << "_mom" << (-p_pi1_src).file_str(2) << "_mom" << (-p_pi1_snk).file_str(2);
#endif
	os << "_src_" << src_type << "_snk_" << snk_type;
	os << "_v2";
        pipi[i].write(os.str());
        os.str("");
      }
    }
    if(src_type == snk_type)
    {
      if(!UniqueID()){ printf("Doing pipi2ptVdis, pidx=%d\n",psrc); fflush(stdout); }
      fVector<typename A2Apolicies::ScalarComplexType> pipi;
      computepipi_v2<A2Apolicies>::compute_Vdis_v2(pipi, mf_1s_con, mf_2s_con, pion_mom, psrc, params.jp.pipi_separation, src_type);
      std::ostringstream os;
      os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_FigureVdis_sep" << params.jp.pipi_separation;
#ifndef DAIQIAN_PION_PHASE_CONVENTION
      os << "_mom" << p_pi1_src.file_str(2);
#else
      os << "_mom" << (-p_pi1_src).file_str(2);
#endif
      os << "_src_" << src_type;
      os << "_v2";
      pipi.write(os.str());
      os.str("");
    }
  }
  time += dclock();
  print_time("main","pipi2pt_v2",time);

  if(!UniqueID()) printf("Memory after pipi2pt_v2 computation:\n");
  printMem();
}

void computecomovepipi_v2(MesonFieldMomentumContainer<A2Apolicies> &mf_pi_con, const StandardPionMomentaPolicy &pion_mom, const int conf, const Parameters &params)
{
  const int nmom = pion_mom.nMom();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  if(!UniqueID()) printf("Computing comove pipi2pt function\n");
  double time = -dclock();

  //here we first do the case when two src and two snk pion have the same momentum, which are 8 different cases
  //according to PHYSICAL REVIEW D 86, 094513 (2012), that number could be at least doubled, since p_cm_src can be opposite to p_cm_snk
  //TYPE1:
  {
    for(int psrc1=0; psrc1 < nmom; psrc1++)
    {
#ifdef P_WAVE_ELIMINATE
      int psrc2 = psrc1, psnk1 = psrc1, psnk2 = psrc1;
      if(!UniqueID()){ printf("Doing comove pipi2pt P wave eliminating part, psrc1=%d psrc2=%d psnk1=%d psnk2=%d\n",psrc1,psrc2,psnk1,psnk2); fflush(stdout); }
#else
      int psrc2 = psrc1, psnk1 = psrc1 + 1 - 2 * (psrc1 % 2), psnk2 = psrc1 + 1 - 2 * (psrc1 % 2);
      if(!UniqueID()){ printf("Doing comove pipi2pt, psrc1=%d psrc2=%d psnk1=%d psnk2=%d\n",psrc1,psrc2,psnk1,psnk2); fflush(stdout); }
#endif
      std::array<fMatrix<typename A2Apolicies::ScalarComplexType>,3> pipi;
      compute_comove_pipi_v2<A2Apolicies>::compute_type1_v2(pipi, mf_pi_con, pion_mom, psrc1, psrc2, psnk1, psnk2, params.jp.pipi_separation, params.jp.tstep_pipi);
      char diag[3] = {'C', 'D', 'R'};
      for(int i=0; i<3; i++)
      {
        std::ostringstream os;
        os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_Figure" << diag[i] <<"_sep" << params.jp.pipi_separation<<"_comove_pipi2pt";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
        os << "_srcmom1" << pion_mom.getMesonMomentum(psrc1).file_str() << "_srcmom2" << pion_mom.getMesonMomentum(psrc2).file_str();
        os << "_snkmom1" << pion_mom.getMesonMomentum(psnk1).file_str() << "_snkmom2" << pion_mom.getMesonMomentum(psnk2).file_str();
#else
        os << "_srcmom1" << (-pion_mom.getMesonMomentum(psrc1)).file_str() << "_srcmom2" << (-pion_mom.getMesonMomentum(psrc2)).file_str();
        os << "_snkmom1" << (-pion_mom.getMesonMomentum(psnk1)).file_str() << "_snkmom2" << (-pion_mom.getMesonMomentum(psnk2)).file_str();
#endif
        os << "_v2";
        pipi[i].write(os.str());
        os.str("");
      }
    }
  }

  //here we do the rest cases when p_CM has two/one component that are non zero
  //TYPE2:
  {
    for(int psrc1=0; psrc1 < nmom; psrc1++)
    {
      for(int psrc2=0; psrc2 < nmom; psrc2++)
      {
        ThreeMomentum p_src1 = pion_mom.getMesonMomentum(psrc1);
        ThreeMomentum p_src2 = pion_mom.getMesonMomentum(psrc2);
        ThreeMomentum p_tmp = p_src1 + p_src2;
        if(p_tmp.nzero()<1 || p_tmp.nzero()>2)
          continue;
        for(int psnk1=0; psnk1 < nmom; psnk1++)
        {
          for(int psnk2=0; psnk2 < nmom; psnk2++)
          {
            ThreeMomentum p_snk1 = pion_mom.getMesonMomentum(psnk1);
            ThreeMomentum p_snk2 = pion_mom.getMesonMomentum(psnk2);
#ifdef P_WAVE_ELIMINATE
            if(p_tmp != (p_snk1 + p_snk2))
              continue;
            if(!UniqueID()){ printf("Doing comove pipi2pt P wave eliminating part, psrc1=%d psrc2=%d psnk1=%d psnk2=%d\n",psrc1,psrc2,psnk1,psnk2); fflush(stdout);}
#else
            if(p_tmp != -(p_snk1 + p_snk2))
              continue; 
            if(!UniqueID()){ printf("Doing comove pipi2pt, psrc1=%d psrc2=%d psnk1=%d psnk2=%d\n",psrc1,psrc2,psnk1,psnk2); fflush(stdout); }
#endif
            std::array<fMatrix<typename A2Apolicies::ScalarComplexType>,3> pipi;
            compute_comove_pipi_v2<A2Apolicies>::compute_type23_v2(pipi, mf_pi_con, pion_mom, psrc1, psrc2, psnk1, psnk2, params.jp.pipi_separation, params.jp.tstep_pipi);
            char diag[3] = {'C', 'D', 'R'};
            for(int i=0; i<3; i++)
            {
              std::ostringstream os;
              os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_Figure" << diag[i] <<"_sep" << params.jp.pipi_separation<<"_comove_pipi2pt";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
              os << "_srcmom1" << pion_mom.getMesonMomentum(psrc1).file_str() << "_srcmom2" << pion_mom.getMesonMomentum(psrc2).file_str();
              os << "_snkmom1" << pion_mom.getMesonMomentum(psnk1).file_str() << "_snkmom2" << pion_mom.getMesonMomentum(psnk2).file_str();
#else
              os << "_srcmom1" << (-pion_mom.getMesonMomentum(psrc1)).file_str() << "_srcmom2" << (-pion_mom.getMesonMomentum(psrc2)).file_str();
              os << "_snkmom1" << (-pion_mom.getMesonMomentum(psnk1)).file_str() << "_snkmom2" << (-pion_mom.getMesonMomentum(psnk2)).file_str();
#endif
              os << "_v2";
              pipi[i].write(os.str());
              os.str("");
            }
          }
        }
      }
    }
  }
  for(int psrc=0; psrc<nmom; psrc++)
  {
    for(int psnk=0; psnk<nmom; psnk++)
    {
      if(!UniqueID()){ printf("Doing comove pipi2ptVdis, psrc=%d psnk=%d\n",psrc,psnk); fflush(stdout); }
      fVector<typename A2Apolicies::ScalarComplexType> pipi;
      compute_comove_pipi_v2<A2Apolicies>::compute_Vdis_v2(pipi, mf_pi_con, pion_mom, psrc, psnk, params.jp.pipi_separation);
      std::ostringstream os;
      os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_FigureVdis_sep" << params.jp.pipi_separation<<"_comove_pipi2pt";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
      os << "_srcmom" << pion_mom.getMesonMomentum(psrc).file_str() << "_snkmom" << pion_mom.getMesonMomentum(psnk).file_str();
#else
      os << "_srcmom" << (-pion_mom.getMesonMomentum(psrc)).file_str() << "_snkmom" << (-pion_mom.getMesonMomentum(psnk)).file_str();
#endif
      os << "_v2";
      pipi.write(os.str());
      os.str("");
    }
  }
  time += dclock();
  print_time("main","comove_pipi2pt",time);

  if(!UniqueID()) printf("Memory after comove_pipi2pt computation:\n");
  printMem();
}

void computecomovepipi(const StandardPionMomentaPolicy &pion_mom, const int conf, const Parameters &params)
{
  const int nmom = pion_mom.nMom();
  const int Lt = GJP.Tnodes() * GJP.TnodeSites();
  if(!UniqueID()) printf("Computing comove pipi2pt function\n");
  double time = -dclock();
  //here we first do the case when two src and two snk pion have the same momentum, which are 8 different cases
  for(int psrc1=0; psrc1 < nmom; psrc1++)
  {
#ifdef TIANLE_COMOVE_TEST
    int psrc2 = psrc1 + 1 - 2 * (psrc1 % 2), psnk1 = psrc1, psnk2 = psrc1 + 1 - 2 * (psrc1 % 2);
#else
    int psrc2 = psrc1, psnk1 = psrc1 + 1 - 2 * (psrc1 % 2), psnk2 = psrc1 + 1 - 2 * (psrc1 % 2);
#endif
    if(!UniqueID()){ printf("Doing comove pipi2pt, psrc1=%d psrc2=%d psnk1=%d psnk2=%d\n",psrc1,psrc2,psnk1,psnk2); fflush(stdout); }
    fMatrix<typename A2Apolicies::ScalarComplexType> pipiC(Lt,Lt);
    fMatrix<typename A2Apolicies::ScalarComplexType> pipiD(Lt,Lt);
    fMatrix<typename A2Apolicies::ScalarComplexType> pipiR(Lt,Lt);
    compute_comove_pipi<A2Apolicies>::compute_v2(pipiC, pipiD, pipiR, pion_mom, params.jp.pipi_separation, params.jp.tstep_pipi, psrc1, psrc2, psnk1, psnk2, conf, params.meas_arg.WorkDirectory);
    std::ostringstream os;
    os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_FigureC_sep" << params.jp.pipi_separation<<"_comove_pipi2pt";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
    os << "_srcmom1" << pion_mom.getMesonMomentum(psrc1).file_str() << "_srcmom2" << pion_mom.getMesonMomentum(psrc2).file_str();
    os << "_snkmom1" << pion_mom.getMesonMomentum(psnk1).file_str() << "_snkmom2" << pion_mom.getMesonMomentum(psnk2).file_str();
#else
    os << "_srcmom1" << (-pion_mom.getMesonMomentum(psrc1)).file_str() << "_srcmom2" << (-pion_mom.getMesonMomentum(psrc2)).file_str(); 
    os << "_snkmom1" << (-pion_mom.getMesonMomentum(psnk1)).file_str() << "_snkmom2" << (-pion_mom.getMesonMomentum(psnk2)).file_str();
#endif
    pipiC.write(os.str());
    os.str("");

    os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_FigureD_sep" << params.jp.pipi_separation<<"_comove_pipi2pt";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
    os << "_srcmom1" << pion_mom.getMesonMomentum(psrc1).file_str() << "_srcmom2" << pion_mom.getMesonMomentum(psrc2).file_str(); 
    os << "_snkmom1" << pion_mom.getMesonMomentum(psnk1).file_str() << "_snkmom2" << pion_mom.getMesonMomentum(psnk2).file_str();
#else
    os << "_srcmom1" << (-pion_mom.getMesonMomentum(psrc1)).file_str() << "_srcmom2" << (-pion_mom.getMesonMomentum(psrc2)).file_str();
    os << "_snkmom1" << (-pion_mom.getMesonMomentum(psnk1)).file_str() << "_snkmom2" << (-pion_mom.getMesonMomentum(psnk2)).file_str();
#endif
    pipiD.write(os.str());
    os.str("");

    os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_FigureR_sep" << params.jp.pipi_separation<<"_comove_pipi2pt";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
    os << "_srcmom1" << pion_mom.getMesonMomentum(psrc1).file_str() << "_srcmom2" << pion_mom.getMesonMomentum(psrc2).file_str(); 
    os << "_snkmom1" << pion_mom.getMesonMomentum(psnk1).file_str() << "_snkmom2" << pion_mom.getMesonMomentum(psnk2).file_str();
#else
    os << "_srcmom1" << (-pion_mom.getMesonMomentum(psrc1)).file_str() << "_srcmom2" << (-pion_mom.getMesonMomentum(psrc2)).file_str();
    os << "_snkmom1" << (-pion_mom.getMesonMomentum(psnk1)).file_str() << "_snkmom2" << (-pion_mom.getMesonMomentum(psnk2)).file_str();
#endif
    pipiR.write(os.str());
    os.str("");
  }

  for(int psrc=0; psrc<nmom; psrc++)
  {
    int psnk = psrc;
    if(!UniqueID()){ printf("Doing comove pipi2ptVdis, psrc=%d psnk=%d\n",psrc,psnk); fflush(stdout); }
    fVector<typename A2Apolicies::ScalarComplexType> pipi;
    compute_comove_pipi<A2Apolicies>::compute_Vdis_v2(pipi, pion_mom, params.jp.pipi_separation, psrc, psnk, conf, params.meas_arg.WorkDirectory);
    std::ostringstream os;
    os << params.meas_arg.WorkDirectory << "/traj_" << conf << "_FigureVdis_sep" << params.jp.pipi_separation<<"_comove_pipi2pt";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
    os << "_srcmom" << pion_mom.getMesonMomentum(psrc).file_str() << "_snkmom" << pion_mom.getMesonMomentum(psnk).file_str();
#else
    os << "_srcmom" << (-pion_mom.getMesonMomentum(psrc)).file_str() << "_snkmom" << (-pion_mom.getMesonMomentum(psnk)).file_str();
#endif
    pipi.write(os.str());
    os.str("");
  }
  time += dclock();
  print_time("main","comove_pipi2pt",time);

  if(!UniqueID()) printf("Memory after comove_pipi2pt computation:\n");
  printMem();
}


void doConfiguration(const int conf, Parameters &params, const CommandLineArgs &cmdline)
{
  params.meas_arg.TrajCur = conf;
  std::string dir(params.meas_arg.WorkDirectory);

  StandardPionMomentaPolicy pion_mom; 
  StationarySigmaMomentaPolicy sigma_mom;

  MesonFieldMomentumContainer<A2Apolicies> mf_ll_con; 
  MesonFieldMomentumContainer<A2Apolicies> mf_ll_con_2s; //Gparity only, might be a bad choice?
  MesonFieldMomentumPairContainer<A2Apolicies> mf_sigma_con;

  //------------------------------Read pion and sigma mesonfield------------------------------------------
  if(cmdline.do_comove_pipi || cmdline.do_pipisigma || cmdline.do_pipiv2 || cmdline.do_pion2ptv3)
    read_pion_mesonfield(mf_ll_con, pion_mom, conf, params, params.meas_arg.WorkDirectory);
  if(cmdline.do_sigma || cmdline.do_pipisigma)
    read_sigma_mesonfield(mf_sigma_con, sigma_mom, conf, params, params.meas_arg.WorkDirectory);
  if(cmdline.do_pion2ptv3 || cmdline.do_pipiv2)
    read_2s_pion_mesonfield(mf_ll_con_2s, pion_mom, conf, params, params.meas_arg.WorkDirectory);


  //----------------------------Compute the pion and moving pion two-point function---------------------------------
  if(cmdline.do_pion2ptv2) computepion2ptv2(mf_ll_con, pion_mom, conf, params, cmdline.do_1s);
  if(cmdline.do_pion2ptv3) computePion2pt_v3(mf_ll_con, mf_ll_con_2s, pion_mom, conf, params, std::string("1s"), std::string("1s"));
  if(cmdline.do_pion2ptv3) computePion2pt_v3(mf_ll_con, mf_ll_con_2s, pion_mom, conf, params, std::string("1s"), std::string("2s"));
  if(cmdline.do_pion2ptv3) computePion2pt_v3(mf_ll_con, mf_ll_con_2s, pion_mom, conf, params, std::string("2s"), std::string("2s"));
  if(cmdline.do_pion2ptcomovev2) compute_comove_pion2ptv2(mf_ll_con, pion_mom, conf, params, cmdline.do_1s);

  //----------------------------Compute the sigma two-point function and sigma/pipi matrix element---------------------------------
  if(cmdline.do_sigma) 
  {
    //computesigma2ptv2(mf_ll_con, sigma_mom, conf, params, cmdline.do_1s);
    computesigma2ptv3(mf_sigma_con, sigma_mom, conf, params);
  }
  if(cmdline.do_pipisigma) 
  {
    //computepipisigmav2(mf_ll_con, pion_mom, sigma_mom, conf, params, cmdline.do_1s);
    computepipitosigmav3(mf_ll_con, mf_sigma_con, pion_mom, sigma_mom, conf, params);
  }
//  if(cmdline.do_sigmapipi) computesigmapipiv2(mf_ll_con, pion_mom, sigma_mom, conf, params);
  
  //------------------------------I=0 and I=2 pipi and moving pipi two-point function---------------------------------
  if(cmdline.do_pipiv2) computepipi2ptv2(mf_ll_con, mf_ll_con_2s, pion_mom, conf, params, std::string("1s"), std::string("1s"));
  if(cmdline.do_pipiv2) computepipi2ptv2(mf_ll_con, mf_ll_con_2s, pion_mom, conf, params, std::string("1s"), std::string("2s"));
  if(cmdline.do_pipiv2) computepipi2ptv2(mf_ll_con, mf_ll_con_2s, pion_mom, conf, params, std::string("2s"), std::string("2s"));

  if(cmdline.do_comove_pipi)
  { 
    computecomovepipi_v2(mf_ll_con, pion_mom, conf, params);
//    computecomovepipi(pion_mom, conf, params);
  }
}


void checkWriteable(const std::string &dir,const int conf){
  std::string file;
  {
    std::ostringstream os; os << dir << "/writeTest.node" << UniqueID() << ".conf" << conf;
    file = os.str();
  }
  std::ofstream of(file);
  double fail = 0;
  if(!of.good()){ std::cout << "checkWriteable failed to open file for write: " << file << std::endl; std::cout.flush(); fail = 1; }

  of << "Test\n";
  if(!of.good()){ std::cout << "checkWriteable failed to write to file: " << file << std::endl; std::cout.flush(); fail = 1; }

  glb_sum_five(&fail);

  if(fail != 0.){
    if(!UniqueID()){ printf("Disk write check failed\n");  fflush(stdout); }
    exit(-1);
  }else{
    if(!UniqueID()){ printf("Disk write check passed\n"); fflush(stdout); }
  }
}


#endif
