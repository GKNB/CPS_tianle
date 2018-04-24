#define NODE_DISTRIBUTE_MESONFIELDS //Save memory by keeping meson fields only on single node until needed

#include <alg/alg_fix_gauge.h>
#include <alg/a2a/utils_main.h>
#include <alg/a2a/grid_wrappers.h>
#include <alg/a2a/bfm_wrappers.h>
#include <alg/a2a/compute_kaon.h>
#include <alg/a2a/compute_pion.h>
#include <alg/a2a/compute_sigma.h>
#include <alg/a2a/compute_pipi.h>
#include <alg/a2a/compute_ktopipi.h>
#include "compute_comove_pion_v2.h"
#include "compute_pion_v2.h"
#include "compute_sigma_v2.h"
#include "compute_sigma_v3.h"
#include "compute_pipitosigma_v2.h"
#include "compute_pipitosigma_v3.h"
#include "compute_sigmatopipi_v2.h"
#include "compute_pipi_v2.h"
#include "fit_v2.h"
#include "compute_comove_pipi.h"
#include "compute_comove_pipi_v2.h"
#include "main.h"

int main (int argc,char **argv )
{
  const char *fname="main(int,char**)";
  Start(&argc, &argv);

  const char *cname=argv[0];
  const int TrajStart = atoi(argv[2]);
  const int LessThanLimit = atoi(argv[3]);

  CommandLineArgs cmdline(argc,argv,4); //control the functioning of the program from the command line    
  Parameters params(argv[1]);

  setupJob(argc, argv, params, cmdline);

  const int Lt = GJP.Tnodes()*GJP.TnodeSites();

#if defined(USE_BFM_A2A) || defined(USE_BFM_LANCZOS)
  BFMsolvers bfm_solvers(cmdline.nthreads, 0.01, 1e-08, 20000, params.jp.solver, params.jp.mobius_scale); //for BFM holds a double and single precision bfm instance. Mass is not important as it is changed when necessary
#endif
  
  if(chdir(params.meas_arg.WorkDirectory)!=0) ERR.General("",fname,"Unable to switch to work directory '%s'\n",params.meas_arg.WorkDirectory);
  double time;

  if(!UniqueID()) printf("Memory prior to config loop:\n");
  printMem();

  //-------------------- Main Loop Begin! -------------------- //
  for(int conf = TrajStart; conf < LessThanLimit; conf += params.meas_arg.TrajIncrement) {
    double conf_time = -dclock();
    if(!UniqueID()) std::cout<<"Starting configuration "<<conf<< std::endl;
    doConfiguration(conf,params,cmdline);
    conf_time += dclock();
    print_time("main","Configuration total",conf_time);
  }//end of config loop

  if(!UniqueID()) printf("Done\n");
  End();
}

