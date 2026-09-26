#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

#include "src/facet.h"
#include "src/mesh/mesh.h"
#include "src/FE/basis_func.h"
#include "src/FE/local_calc.h"

// GABranch is a Voronkov correction status, not a kinetic branch.
// Solidification versus melting is decided only by the sign of dTkin.
static const int VOR_GA_NO_CORRECTION = 0;
static const int VOR_GA_FORCE_BALANCE_STOP = 1;
static const int VOR_GA_VELOCITY_CONTROLLED = 2;
static const int VOR_GA_INVALID_GEOMETRY = 3;

static const int PHIS_INVALID = 0;
static const int PHIS_ROUGH_ISOTROPIC = 1;
static const int PHIS_VICINAL_CH = 2;
static const int PHIS_EXACT_FACET_LIMIT_OR_SUPPRESSED = 3;

// -----------------------------------------------------------------------------
// Temperature sign convention
// -----------------------------------------------------------------------------
// The FE solver stores temperature relative to the melting point:
//     Tfe = T - Tmp.
// Thus:
//     Tfe < 0 : undercooled, solidification is thermodynamically possible.
//     Tfe > 0 : overheated, melting is thermodynamically possible.
//     Tfe = 0 : melting-point isotherm.
//
// Interface kinetics use the signed driving force:
//     dTkin = Tmp - Ti = -Tfe.
// Thus:
//     dTkin > 0 : solidification branch.
//     dTkin < 0 : melting branch.
//     dTkin = 0 : equilibrium.
//
// Never define the physical kinetic dT with max(0,...).  Use the signed
// dTkin to select the branch; one-sided names such as dTg and dTm are local
// aliases only after that branch is selected.
// -----------------------------------------------------------------------------

struct FEStepContext {
  LocalCalc local_calc;
  vector<LocalCalc> thread_calcs;
  vector<double> dTdt_tls;
  vector<double> cT_tls;
  GaussPointBFSet gp_set[NUM_OF_BF_SETS];
  GaussPointBFSet gp_set_corn;
  double dt;
  bool use_openmp;
};

struct AxisymmetricInitialState {
  int version;
  int crystalRadialCells;
  int outerRadialCells;
  int liquidVerticalCells;
  int solidVerticalCells;
  double bottomZ;
  double undisturbedMeltZ;
  double tplZ;
  double topZ;
  double crystalRadius;
  double crucibleRadius;
  double temperatureTop;
  double temperatureOut;
  double temperatureBottom;
  double targetTplTemperature;
  double bottomSensitivity;
  double diffusivityLiquid;
  double diffusivitySolid;
  double heatTransferCoefficient;
  double radiationCoefficient;
  double latentHeatCoefficient;
  double pullVelocity;
  double sigmaMG;
  double rhoG;
  double meltingTemperature;
  vector<double> liquidT;
  vector<double> solidT;
  vector<double> annulusT;
  vector<double> slR;
  vector<double> slZ;
  vector<double> meniscusR;
  vector<double> meniscusZ;
  vector<double> meniscusArc;
};

static int getEnvInt(const char *name, int fallback);
static double runFEtemperatureStep(Mesh *p_mesh, FEStepContext *p_fe_step);
static int ensureDirectoryExists(const char *path, ostream &outw);
static void initParam(Mesh *p_mesh, LocalCalc *p_local_calc,
                      const char *filename);
static void initFEStepStorage(FEStepContext *p_fe_step, Mesh *p_mesh,
                               ostream &outw);
static void writeMeshState(const char *filename, Mesh *p_mesh,
                           LocalCalc *p_local_calc,
                           int n, double d_crucible, double end_time,
                           double GAf, double GAr,
                           double Rc, double lc, double hSlab);
static int readMeshState(const char *filename, Mesh *p_mesh,
                         LocalCalc *p_local_calc,
                         int n, double d_crucible,
                         double GAf, double GAr,
                         ostream &outw);
static int readAxisymmetricInitialState(
    const char *filename, Mesh *p_mesh, LocalCalc *p_local_calc,
    double d_crucible, double Rc, ostream &outw,
    double *p_target_tpl_temperature, double *p_bottom_sensitivity);
static int initialize3dInitialState(
    const char *ssStateInputFile, const char *axisymmetricInitialStateFile,
    Mesh *p_mesh, LocalCalc *p_local_calc, int n, double d_crucible,
    double Rc, double GAf, double GAr, ostream &outw,
    int *p_restartFromStoredState,
    int *p_initializedFromAxisymmetricState,
    double *p_axisymmetricTplTarget,
    double *p_axisymmetricBottomSensitivity);
static void writeDriverOutput(Mesh *p_mesh, Mesh *p_metrics_mesh,
                              ostream &outw,
                              double time_value, double norm, double GAr,
                              double GAf, const char *outputDir,
                              int roughOnly, int writeDataFiles,
                              int writeHeader);
static void initCreateTplSubmesh(Mesh *p_local, Mesh *p_coarse,
                                  int supportLayers,
                                  int liquidRadialElements);
static void initCreateRefinedTplSubmesh(Mesh *p_local, Mesh *p_parent,
                                         int phiRefinement);
static void tplCopyPhysics(Mesh *p_local, const Mesh *p_coarse);
static void refreshTplSubmeshArtificialGeometry(Mesh *p_local,
                                                Mesh *p_parent);
static void tplMapLocalNodes(Mesh *p_local, const Mesh *p_parent);
static void transferTplSubmeshTemperature(Mesh *p_local,
                                           const Mesh *p_parent);
static double solveTplSubmeshTemperature(Mesh *p_local,
                                          const Mesh *p_parent,
                                         LocalCalc *p_local_calc,
                                         GaussPointBFSet *p_gp_brick,
                                         GaussPointBFSet *p_gp_quad,
                                         double physicalDt,
                                         const double *p_oldSolidGradient,
                                         const double *p_oldLiquidGradient,
                                         int *p_thermalSubsteps);
static void calcTplSubmeshGradients(Mesh *p_local,
                                    LocalCalc *p_local_calc,
                                    GaussPointBFSet *p_gp_corners);
static void deformTplSubmeshInterior(Mesh *p_local,
                                      const Mesh *p_parent);
static void moveTplSubmeshInterfaces(Mesh *p_local,
                                      const Mesh *p_parent,
                                     double physicalDt,
                                     double geometryDt,
                                     double GAr, double GAf);
static void restrictTplSubmeshToCoarse(Mesh *p_local, Mesh *p_parent);
static int writeTplSubmeshState(const char *filename, const Mesh *p_local,
                                int level);
static int writeTplCircularBenchmark(const char *filename, Mesh *p_local,
                                     int recursiveLevels,
                                     int halfColumns, ostream &outw);



int main() {
  const char *outputDir = "output/3dFEexp";
  if (!ensureDirectoryExists("output", cerr) ||
      !ensureDirectoryExists(outputDir, cerr))
    return 1;

  char outwFilename[100];
  sprintf(outwFilename, "%s/outw", outputDir);
  ofstream outw(outwFilename, ios::out);
  if (!outw) {
    cerr << "Failed to open " << outwFilename << "\n";
    return 1;
  }
  outw.setf(ios::unitbuf);
  streambuf *p_cout_buf = cout.rdbuf(outw.rdbuf());
  streambuf *p_cerr_buf = cerr.rdbuf(outw.rdbuf());

  const int n = getEnvInt("FEEXP_COARSE_N", 40);
  const double d_crucible = .1;
  const int max_iters = getEnvInt("FEEXP_STARTUP_ITERS", 20000);
  const int end_steps = getEnvInt("FEEXP_END_STEPS", -1);
  const int mesh_update_every = 5;
  const int relaxationOutputEvery = 500;
  // Resolve the narrow alpha interval betaRough/betaStep before the S/L
  // kinetics reaches the rough upper limit.
  // The downloaded n=35 run inverted the finest liquid cells near the facet
  // with 12 movement substeps. Keep the implicit thermal step unchanged and
  // reduce only the explicit kinetic/interface increment by a factor of four.
  const int nMoveSubsteps = 48;
  const double maxGeometrySubstep =
      .00025 / (double)nMoveSubsteps; // [s]
  // Source-level comparison switch.  Set to 0 for the original coarse-only
  // calculation without introducing a runtime mode layer.
  const int useTplSubmesh = 1;
  const int localSupportLayers = 2;
  const int localLiquidRadialElements = 8;
  const int recursiveSubmeshLevels = 3;
  const int recursivePhiRefinement = 2;
  const int ridgeBenchmarkHalfCoarseIntervals = 2;
  const int startupSurfaceRelaxIterations = 5;
  const int geometryRelaxIterations = 5000;
  const int submeshRelaxIterations = 2000;
  const double geometryLgMoveCoef = .01;
  // Pre-main L3 TPL+pull isotropic wall: accumulate geometryDt to this many
  // seconds (geometryLgMoveCoef per full-rate pass). Main-loop end_time below
  // still counts physicalDtMax thermal steps; this iso stage uses the existing
  // TPL/YL geometryDt clock, not a new mode flag.
  const double preMainIsotropicTime = 5.; // [s] geometryDt before transient
  const double output_interval = 1.;
  // Test horizon for Vpull=1.6 mm/min: 5 s main after 5 s pre-main iso.
  // The former 400 s horizon belonged to Vpull=2 mm/h (0.222 mm length).
  const double end_time = 5.; // [s]
  // The previous 0.005 s step produced persistent node-to-node kinetic
  // oscillations near the singular facet. Keep the physical step below the
  // local rough kinetic/thermal relaxation time.
  const double physicalDtMax = .001; // [s]
  const double ssCheckInterval = .1; // [s]
  const int ssRequiredChecks = 10;
  const double ssVlabTolerance = .1 / 1000. / 3600.; // [m/s], 0.00167 mm/min

  Mesh mesh = {};
  Mesh local_mesh = {};
  Mesh refined_mesh_1 = {};
  Mesh refined_mesh_2 = {};
  Mesh refined_mesh_3 = {};
  vector<double> coarseSolidGradientOld;
  vector<double> coarseLiquidGradientOld;
  vector<double> localSolidGradientOld;
  vector<double> localLiquidGradientOld;
  vector<double> refined1SolidGradientOld;
  vector<double> refined1LiquidGradientOld;
  vector<double> refined2SolidGradientOld;
  vector<double> refined2LiquidGradientOld;
  clock_t t;
  t = clock();
  FEStepContext fe_step = {};
  mesh.initCreateMesh3d1(n, d_crucible);
  coarseSolidGradientOld.resize(3 * mesh.number_of_nodes, 0.);
  coarseLiquidGradientOld.resize(3 * mesh.number_of_nodes, 0.);
  mesh.bc_set[mesh.BC_SOLID_LIQUID].relaxMoveCoef = .001;
  LocalCalc &local_calc = fe_step.local_calc;
  local_calc.LocalCalcInit();
  initParam(&mesh, &local_calc, "param.txt");
  fe_step.dt = .05 * mesh.dx;
  fe_step.use_openmp = true;

  mesh.anisotropicKinFlg = 0;
  mesh.rChangeFlg = 0;
  mesh.AniGaFlg = 0;

  const double GAf =
      180. * mesh.linearFacetAppAngle / PI; // [deg], calibrated linear tip GAapp
  const double GAr =
      180. * mesh.roughGrowthAngle / PI;    // [deg], calculated rough GAeq
  const double Rc = 0.25 * d_crucible;      // [m], crystal radius
  const double phi = mesh.roughGrowthAngle; // [rad]
  const double lc = sqrt(mesh.gam / mesh.rhog);
  const double sinPhi = sin(phi);
  const double hSlab =
      lc * sqrt(2.0 * max(0.0, 1.0 - sinPhi));
  // CZ second-order engineering approximation:
  // h(r)=h_slab*[1+a/(sqrt(2)*r)+a^2/(4*r^2)*(1+sin(alpha))/(1-sin(alpha))]^(-1/2).
  local_calc.h0 =
      hSlab /
      sqrt(1.0 + lc / (sqrt(2.0) * Rc) +
           lc * lc / (4.0 * Rc * Rc) *
           (1.0 + sinPhi) / (1.0 - sinPhi));

  char ssStateInputFile[100];
  const char *axisymmetricInitialStateFile =
      "input/axisymmetric_ic.dat";
  sprintf(ssStateInputFile, "input/ss_state_n%04d.dat", n);
  int restartFromStoredState = 0;
  int initializedFromAxisymmetricState = 0;
  double axisymmetricTplTarget = 0.;
  double axisymmetricBottomSensitivity = 0.;

  if (!initialize3dInitialState(
          ssStateInputFile, axisymmetricInitialStateFile,
          &mesh, &local_calc, n, d_crucible, Rc, GAf, GAr, outw,
          &restartFromStoredState, &initializedFromAxisymmetricState,
          &axisymmetricTplTarget, &axisymmetricBottomSensitivity)) {
    cout.rdbuf(p_cout_buf);
    cerr.rdbuf(p_cerr_buf);
    return 1;
  }
  if (restartFromStoredState)
    fe_step.dt = .05 * mesh.dx;

  initFEStepStorage(&fe_step, &mesh, outw);

  double dt = fe_step.dt;

  int omp_threads = 1;
#ifdef _OPENMP
  omp_threads = omp_get_max_threads();
#endif

  {
    const double r_crystal = 0.25 * d_crucible;
    const double r_relax_threshold = 0.5 * r_crystal;
    for (int ib(0); ib < mesh.bc_set[mesh.BC_SOLID_LIQUID].number_of_recordsNode; ib++) {
      const double x = mesh.bc_set[mesh.BC_SOLID_LIQUID].bc_at_node[ib].p_node->coord[0];
      const double y = mesh.bc_set[mesh.BC_SOLID_LIQUID].bc_at_node[ib].p_node->coord[1];
      const double r_node = sqrt(x * x + y * y);
      mesh.bc_set[mesh.BC_SOLID_LIQUID].bc_at_node[ib].relax_coef = 1. + 5. * (r_node - r_relax_threshold) /
                                                                            (r_crystal - r_relax_threshold);
      if (mesh.bc_set[mesh.BC_SOLID_LIQUID].bc_at_node[ib].relax_coef < 1.)
          mesh.bc_set[mesh.BC_SOLID_LIQUID].bc_at_node[ib].relax_coef = 1.;
      if (mesh.bc_set[mesh.BC_SOLID_LIQUID].bc_at_node[ib].relax_coef > 4.)
          mesh.bc_set[mesh.BC_SOLID_LIQUID].bc_at_node[ib].relax_coef = 4.;
    }
  }

  double relaxationNorm = 0.;

  // Only the explicit no-state path requires coarse thermal/isotropic startup.
  if (!restartFromStoredState) {
    if (initializedFromAxisymmetricState) {
      // Settle the interpolated 2D temperature and 3D interface geometry in
      // one coupled loop. Start interface motion far below its full value so
      // the first coarse FE gradients are established before appreciable
      // L/G or S/L displacement occurs. Growth-angle-driven TPL motion starts
      // only in the physical main loop.
      outw << "relaxation | coupled axisymmetric import"
           << " | move=0.001,0.01,0.1,1\n";
      double norm = 0.;
      for (int it(0); it < geometryRelaxIterations; it++) {
        norm = runFEtemperatureStep(&mesh, &fe_step);
        relaxationNorm = norm;

        double moveCoef = .001;
        if (it >= .25 * geometryRelaxIterations)
          moveCoef = .01;
        if (it >= .50 * geometryRelaxIterations)
          moveCoef = .1;
        if (it >= .75 * geometryRelaxIterations)
          moveCoef = 1.;

        mesh.bc_set[mesh.BC_LIQUID_GAS].MoveLiquidInterface(
            moveCoef * geometryLgMoveCoef, &mesh);
        mesh.bc_set[mesh.BC_SOLID_LIQUID].calcRelaxSL();
        mesh.bc_set[mesh.BC_SOLID_LIQUID].kinMoveInterface(
            moveCoef * dt, &mesh);
        mesh.bc_set[mesh.BC_SOLID_GAS].calcRelaxStr();
        mesh.bc_set[mesh.BC_LIQUID_GAS].calcRelaxStr();

        if (it % mesh_update_every == 0)
          mesh.calc3dmeshDef();
        if (it % relaxationOutputEvery == 0)
          writeDriverOutput(&mesh, &mesh, outw, (double)it, norm,
                            GAr, GAf, outputDir, 1, 0, it == 0);
      }
      outw << "startup | imported=axisymmetric"
           << " | coupled_steps=" << geometryRelaxIterations
           << " | norm=" << norm << "\n";
    } else {
  
      for (int i(0); i < mesh.number_of_nodes; i++)
      {
          mesh.node[i].dTdt = 0.;
          mesh.node[i].dTdt0 = 0.;
          mesh.node[i].cT = 0.;
          mesh.node[i].attribute[TEMPERATURE] = local_calc.Tbottom
          + (mesh.node[i].coord[2]/d_crucible) *
          (local_calc.Ttop - local_calc.Tbottom) ;
      }
  
      double norm ;
      fe_step.dt=2.5* dt;
      outw << "relaxation | initial meniscus\n";
      for (int it(0); it < max_iters+1; it++) {     
        norm = runFEtemperatureStep(&mesh, &fe_step);
        relaxationNorm = norm;
        if (it > .25*max_iters ) {
          for (int ib(0); ib < mesh.bc_set[mesh.BC_SOLID_LIQUID].number_of_recordsNode; ib++)
            if(mesh.bc_set[mesh.BC_SOLID_LIQUID].bc_at_node[ib].p_node->coord[2] >= mesh.ZI + local_calc.h0)
              mesh.bc_set[mesh.BC_SOLID_LIQUID].bc_at_node[ib].p_node->coord[2] -= .0005*mesh.dx ;
            else
              mesh.bc_set[mesh.BC_SOLID_LIQUID].bc_at_node[ib].p_node->coord[2] += .0005*mesh.dx ;
        }
  
        // The S/L level change moves the common TPL. Relax both attached
        // surfaces before applying the next L/G Young-Laplace update.
        for (int ir(0); ir < startupSurfaceRelaxIterations; ir++) {
          mesh.bc_set[mesh.BC_SOLID_GAS].calcRelaxStr();
          mesh.bc_set[mesh.BC_LIQUID_GAS].calcRelaxStr();
        }
        // The initial flat L/G mesh has a sharp corner at the lifted TPL. Use a
        // very small Young-Laplace relaxation rate only during the first quarter.
        double lgMoveCoef = .00001;
        if (it >= .25 * max_iters) lgMoveCoef = .001;
        if (it >= .5 * max_iters) lgMoveCoef = .1;
        mesh.bc_set[mesh.BC_LIQUID_GAS].MoveLiquidInterface(lgMoveCoef * dt,
                                                            &mesh);
  
        if (it % mesh_update_every == 0) mesh.calc3dmeshDef();

        if (it % relaxationOutputEvery == 0)
          writeDriverOutput(&mesh, &mesh, outw, (double)it, norm,
                            GAr, GAf, outputDir, 1, 0, it == 0);
  
          if (it > .5*max_iters && norm < 1.e-6)
          break;
      }
  
      fe_step.dt=dt;
  
      outw << "relaxation | isotropic interfaces\n";
      for (int it(0); it < max_iters+1; it++) {
  
        norm = runFEtemperatureStep(&mesh, &fe_step);
        relaxationNorm = norm;
        mesh.bc_set[mesh.BC_LIQUID_GAS].calcRelaxStr();
        mesh.bc_set[mesh.BC_LIQUID_GAS].MoveLiquidInterface(.01*dt, &mesh);
  
        double dtMoveCoef = .01;
        if (it > .5 * max_iters) dtMoveCoef = .1;
        if (it > .75 * max_iters) dtMoveCoef = 1.;
        mesh.bc_set[mesh.BC_SOLID_LIQUID].kinMoveInterface(dtMoveCoef * dt, &mesh);
        mesh.bc_set[mesh.BC_SOLID_LIQUID].calcRelaxSL();
        for (int ir(0); ir < startupSurfaceRelaxIterations; ir++) {
          mesh.bc_set[mesh.BC_SOLID_GAS].calcRelaxStr();
          mesh.bc_set[mesh.BC_LIQUID_GAS].calcRelaxStr();
        }
  
        if (it % mesh_update_every == 0) mesh.calc3dmeshDef();

        if (it % relaxationOutputEvery == 0)
          writeDriverOutput(&mesh, &mesh, outw, (double)it, norm,
                            GAr, GAf, outputDir, 1, 0, it == 0);
  
      }
      outw << "startup | imported=none"
           << " | isotropic_steps=" << max_iters + 1
           << " | norm=" << norm << "\n";
    }
  }

  if (useTplSubmesh) {
    // Prepare the coarse physical surfaces before cutting out L0. This stage
    // is applied after every startup source, including an axisymmetric import
    // and a 3D restart. The TPL remains fixed; there is no temperature step,
    // growth-angle motion, or crystal pulling. Normal S/L kinetics are applied
    // with a ramped moveCoef * dt schedule starting from a very small value.
    mesh.bc_set[mesh.BC_LIQUID_GAS].MoveLiquidInterface(0., &mesh);
    double initialLgResidual = 0.;
    for (int i(0); i < mesh.bc_set[mesh.BC_LIQUID_GAS].nS1; i++)
      for (int j(1); j < mesh.bc_set[mesh.BC_LIQUID_GAS].nS2 - 1; j++)
        if (mesh.bc_set[mesh.BC_LIQUID_GAS]
                .bc_at_nodeStr[i][j]->move_abs > initialLgResidual)
          initialLgResidual = mesh.bc_set[mesh.BC_LIQUID_GAS]
                                  .bc_at_nodeStr[i][j]->move_abs;
    outw << "relaxation | coarse interfaces"
         << " | move=0.001,0.01,0.1,1\n";
    for (int it(0); it < geometryRelaxIterations; it++) {
      double moveCoef = .001;
      if (it >= .25 * geometryRelaxIterations)
        moveCoef = .01;
      if (it >= .50 * geometryRelaxIterations)
        moveCoef = .1;
      if (it >= .75 * geometryRelaxIterations)
        moveCoef = 1.;

      mesh.bc_set[mesh.BC_LIQUID_GAS].MoveLiquidInterface(
          geometryLgMoveCoef, &mesh);
      mesh.bc_set[mesh.BC_SOLID_LIQUID].calcRelaxSL();
      mesh.bc_set[mesh.BC_SOLID_LIQUID].kinMoveInterface(
          moveCoef * dt, &mesh);
      mesh.bc_set[mesh.BC_SOLID_GAS].calcRelaxStr();
      mesh.bc_set[mesh.BC_LIQUID_GAS].calcRelaxStr();
      if (it % mesh_update_every == 0)
        mesh.calc3dmeshDef();
      if (it % relaxationOutputEvery == 0)
        writeDriverOutput(&mesh, &mesh, outw, (double)it,
                          relaxationNorm, GAr, GAf, outputDir, 1, 0,
                          it == 0);
    }
    mesh.calc3dmeshDef();
    mesh.bc_set[mesh.BC_LIQUID_GAS].MoveLiquidInterface(0., &mesh);
    mesh.bc_set[mesh.BC_TPL].MoveTPL(0., GAr, GAf, 0, &mesh);
    mesh.bc_set[mesh.BC_TPL].MoveTPLvoronkovGA(
        0., GAr, GAf, 0, &mesh);
    double maxLgResidual = 0.;
    for (int i(0); i < mesh.bc_set[mesh.BC_LIQUID_GAS].nS1; i++)
      for (int j(1); j < mesh.bc_set[mesh.BC_LIQUID_GAS].nS2 - 1; j++)
        if (mesh.bc_set[mesh.BC_LIQUID_GAS]
                .bc_at_nodeStr[i][j]->move_abs > maxLgResidual)
          maxLgResidual = mesh.bc_set[mesh.BC_LIQUID_GAS]
                              .bc_at_nodeStr[i][j]->move_abs;
    outw << "geometry | iterations=" << geometryRelaxIterations
         << " | LG_coef=" << geometryLgMoveCoef
         << " | GA="
         << 180. * mesh.bc_set[mesh.BC_TPL]
                        .bc_at_nodeStr[0][1]->GA / PI
         << " | SG="
         << 180. * mesh.bc_set[mesh.BC_TPL]
                        .bc_at_nodeStr[0][1]->SG / PI
         << " | SL="
         << 180. * mesh.bc_set[mesh.BC_TPL]
                        .bc_at_nodeStr[0][1]->theta_SL_TPL / PI
         << " | phiS="
         << 180. * mesh.bc_set[mesh.BC_TPL]
                        .bc_at_nodeStr[0][1]->phiS / PI
         << " | LG_initial=" << initialLgResidual
         << " | LG_final=" << maxLgResidual << "\n";
  }

  if (useTplSubmesh) {
    // Per-level create -> fixed-TPL isotropic -> TPL+pull isotropic -> child
    // (never batch-create then relax): initCreate* this level, T(+grads) from
    // parent, fixed-TPL isotropic on THIS level, then a second same-count stage
    // with MoveTPL/Voronkov + pullSGstr: TPL/YL geometryDt = moveCoef *
    // geometryLgMoveCoef (ramped 0.001->1 x 0.01); kin/pull on moveCoef*dt;
    // tip Tfe restored after MoveTPL (no skip-Taylor arg), then
    // initCreateRefined* for the next child.
    // Freeze = internal S/L bulk / artificial BC_SUB cuts (P->C refresh, incl.
    // Th). Child LG YL on: each pass refresh interpolates LG ends from parent
    // (TPL j=0 + outer j=nS2-1, same exact/sample/element patterns as bulk T),
    // then MoveLiquidInterface moves free-surface interior only. SL calcRelaxSL
    // + ramped kinMoveInterface and SG/LG calcRelaxStr on both stages. Fixed
    // stage has no MoveTPL; TPL+pull stage matches main-loop Voronkov selection.
    // Flag timing: AniGaFlg/anisotropicKinFlg stay 0 through fixed-TPL; enable
    // them on the level at the start of TPL+pull and leave on for main.
    Mesh *refined[3] = {
        &refined_mesh_1, &refined_mesh_2, &refined_mesh_3};
    Mesh *refinedParent[3] = {
        &local_mesh, &refined_mesh_1, &refined_mesh_2};
    Mesh *submesh[4] = {
        &local_mesh, &refined_mesh_1, &refined_mesh_2, &refined_mesh_3};
    Mesh *submeshParent[4] = {
        &mesh, &local_mesh, &refined_mesh_1, &refined_mesh_2};

    // Coarse geometry stage moved nodes without a temperature step; rebuild
    // phase gradients before the first child temperature transfer.
    for (int in(0); in < mesh.number_of_nodes; in++) {
      mesh.node[in].num_of_adj_elem_l = 0;
      mesh.node[in].num_of_adj_elem_s = 0;
      for (int co(0); co < 3; co++) {
        mesh.node[in].dTl[co] = 0.;
        mesh.node[in].dTs[co] = 0.;
      }
    }
    for (int ie(0); ie < mesh.number_of_elements; ie++)
      local_calc.LocalGetGrad3d(&mesh.element[ie], &fe_step.gp_set_corn);
    for (int in(0); in < mesh.number_of_nodes; in++) {
      if (mesh.node[in].num_of_adj_elem_l != 0)
        for (int co(0); co < 3; co++)
          mesh.node[in].dTl[co] /= mesh.node[in].num_of_adj_elem_l;
      if (mesh.node[in].num_of_adj_elem_s != 0)
        for (int co(0); co < 3; co++)
          mesh.node[in].dTs[co] /= mesh.node[in].num_of_adj_elem_s;
    }

    initCreateTplSubmesh(&local_mesh, &mesh, localSupportLayers,
                         localLiquidRadialElements);

    for (int level(0); level < recursiveSubmeshLevels + 1; level++) {
      Mesh *p_sub = submesh[level];
      Mesh *p_parent = submeshParent[level];

      // Parent-controlled bulk/artificial geometry (incl. Th) before T fill.
      refreshTplSubmeshArtificialGeometry(p_sub, p_parent);
      deformTplSubmeshInterior(p_sub, p_parent);

      // Force a full parent remap before the full-field T fill. Physical nodes
      // are normally mapped once at create and then skipped; after refresh /
      // deform (and later isotropic motion) that chart can leave phi-mid holes.
      for (int in(0); in < p_sub->number_of_nodes; in++)
        p_sub->node[in].parent_mapping_initialized = 0;
      tplMapLocalNodes(p_sub, p_parent);
      transferTplSubmeshTemperature(p_sub, p_parent);
      calcTplSubmeshGradients(p_sub, &local_calc, &fe_step.gp_set_corn);

      outw << "relaxation | submesh L" << level
           << " isotropic interfaces"
           << " | iterations=" << submeshRelaxIterations
           << " | freeze=internal bulk"
           << " | LG YL on; ends interp P->C"
           << " | SL kin+redistrib on"
           << " | SG/LG redistrib on"
           << " | TPL=fixed\n";
      for (int it(0); it < submeshRelaxIterations; it++) {
        // Internal bulk freeze + LG-end interpolate from parent (like bulk T),
        // then YL interior, SL kin+redistrib, SG/LG redistrib; no MoveTPL.
        refreshTplSubmeshArtificialGeometry(p_sub, p_parent);
        deformTplSubmeshInterior(p_sub, p_parent);
        double moveCoef = .001;
        if (it >= .25 * submeshRelaxIterations)
          moveCoef = .01;
        if (it >= .50 * submeshRelaxIterations)
          moveCoef = .1;
        if (it >= .75 * submeshRelaxIterations)
          moveCoef = 1.;
        p_sub->bc_set[p_sub->BC_SUB_LIQUID_GAS].MoveLiquidInterface(
            geometryLgMoveCoef, p_sub);
        p_sub->bc_set[p_sub->BC_SUB_SOLID_LIQUID].calcRelaxSL();
        p_sub->bc_set[p_sub->BC_SUB_SOLID_LIQUID].kinMoveInterface(
            moveCoef * dt, p_sub);
        p_sub->bc_set[p_sub->BC_SUB_SOLID_GAS].calcRelaxStr();
        p_sub->bc_set[p_sub->BC_SUB_LIQUID_GAS].calcRelaxStr();
        deformTplSubmeshInterior(p_sub, p_parent);
        // Parent already finished its interface+3d stage before this child
        // was created/relaxed. Restrict surface C->P only; 3d-relax this
        // child after its interface moves (not the parent mid-iter).
        restrictTplSubmeshToCoarse(p_sub, p_parent);
        deformTplSubmeshInterior(p_sub, p_parent);
        if (p_parent == &mesh && it % mesh_update_every == 0)
          mesh.calc3dmeshDef();
        if (level == recursiveSubmeshLevels &&
            it % relaxationOutputEvery == 0)
          writeDriverOutput(&mesh, p_sub, outw, (double)it,
                            relaxationNorm, GAr, GAf, outputDir, 1, 0,
                            it == 0);
      }
      refreshTplSubmeshArtificialGeometry(p_sub, p_parent);
      deformTplSubmeshInterior(p_sub, p_parent);

      // Enable main-loop Voronkov/anisotropic flags for this TPL+pull stage and
      // keep them for the physical main loop (later global enable is idempotent).
      p_sub->anisotropicKinFlg = 1;
      p_sub->AniGaFlg = 1;
      p_sub->bc_set[p_sub->BC_SUB_SOLID_LIQUID].calcFacetAngle();

      // L3 (finest) TPL+pull is the pre-main isotropic wall: accumulate
      // geometryDt to preMainIsotropicTime at full moveCoef. L0-L2 keep the
      // ramped iteration count used while building the hierarchy.
      const int timedPreMainIso = (level == recursiveSubmeshLevels);
      outw << "relaxation | submesh L" << level
           << " TPL+pull isotropic";
      if (timedPreMainIso)
        outw << " | target_s=" << preMainIsotropicTime
             << " | geom_dt=" << geometryLgMoveCoef
             << " | clock=geometryDt";
      else
        outw << " | iterations=" << submeshRelaxIterations
             << " | move=0.001,0.01,0.1,1";
      outw << " | freeze=internal bulk"
           << " | LG YL on; ends interp P->C"
           << " | SL kin+redistrib on"
           << " | SG/LG redistrib on"
           << " | TPL="
           << (p_sub->useVoronkovGA ? "MoveTPLvoronkovGA" : "MoveTPL")
           << "+pull"
           << " | geom=moveCoef*geometryLgMoveCoef"
           << " | tipT=restore after MoveTPL\n";
      double isoGeomTime = 0.;
      for (int it(0);; it++) {
        if (timedPreMainIso) {
          if (isoGeomTime >= preMainIsotropicTime)
            break;
        } else if (it >= submeshRelaxIterations) {
          break;
        }
        // Same bulk freeze + LG-end parent interp as fixed-TPL. TPL/YL use
        // moveCoef*geometryLgMoveCoef (ramped 0.001->1 x 0.01 on L0-L2), not
        // main-loop min(0.05*physicalDt, maxGeometrySubstep). L3 pre-main uses
        // moveCoef=1 so each pass advances geometryLgMoveCoef seconds until
        // preMainIsotropicTime. Kin/pull stay on moveCoef*dt. Mirror main-loop
        // MoveTPL+YL then helper(geometryDt=0) so tip Tfe can be restored after
        // MoveTPL Taylor (no skip-Taylor argument) before kin sees the tip.
        refreshTplSubmeshArtificialGeometry(p_sub, p_parent);
        deformTplSubmeshInterior(p_sub, p_parent);
        double moveCoef = 1.;
        if (!timedPreMainIso) {
          moveCoef = .001;
          if (it >= .25 * submeshRelaxIterations)
            moveCoef = .01;
          if (it >= .50 * submeshRelaxIterations)
            moveCoef = .1;
          if (it >= .75 * submeshRelaxIterations)
            moveCoef = 1.;
        }
        const double physicalSubstep = moveCoef * dt;
        const double geometryDt = moveCoef * geometryLgMoveCoef;
        tplCopyPhysics(p_sub, p_parent);
        BCSet *p_tpl = &p_sub->bc_set[p_sub->BC_SUB_TPL];
        vector<double> tipTfe(p_tpl->nS1);
        for (int i(0); i < p_tpl->nS1; i++)
          tipTfe[i] =
              p_tpl->bc_at_nodeStr[i][1]->p_node->attribute[TEMPERATURE];
        if (p_sub->useVoronkovGA)
          p_tpl->MoveTPLvoronkovGA(geometryDt, GAr, GAf,
                                   p_sub->AniGaFlg, p_sub);
        else
          p_tpl->MoveTPL(geometryDt, GAr, GAf, p_sub->AniGaFlg, p_sub);
        p_sub->bc_set[p_sub->BC_SUB_LIQUID_GAS].MoveLiquidInterface(
            geometryDt, p_sub);
        for (int i(0); i < p_tpl->nS1; i++) {
          Node *p_node = p_tpl->bc_at_nodeStr[i][1]->p_node;
          p_node->attribute[TEMPERATURE] = tipTfe[i];
          p_node->T = tipTfe[i];
        }
        moveTplSubmeshInterfaces(p_sub, p_parent, physicalSubstep,
                                 0., GAr, GAf);
        if (timedPreMainIso)
          isoGeomTime += geometryDt;
        restrictTplSubmeshToCoarse(p_sub, p_parent);
        deformTplSubmeshInterior(p_sub, p_parent);
        if (p_parent == &mesh && it % mesh_update_every == 0)
          mesh.calc3dmeshDef();
        if (level == recursiveSubmeshLevels &&
            it % relaxationOutputEvery == 0)
          writeDriverOutput(&mesh, p_sub, outw,
                            timedPreMainIso ? isoGeomTime : (double)it,
                            relaxationNorm, GAr, GAf, outputDir, 1, 0,
                            it == 0);
      }
      refreshTplSubmeshArtificialGeometry(p_sub, p_parent);
      deformTplSubmeshInterior(p_sub, p_parent);

      // After both isotropic stages, remap and refresh T (+ grads) so the
      // next child samples a finite parent field and L3 ends without nan holes.
      for (int in(0); in < p_sub->number_of_nodes; in++)
        p_sub->node[in].parent_mapping_initialized = 0;
      tplMapLocalNodes(p_sub, p_parent);
      transferTplSubmeshTemperature(p_sub, p_parent);
      calcTplSubmeshGradients(p_sub, &local_calc, &fe_step.gp_set_corn);

      if (level < recursiveSubmeshLevels)
        initCreateRefinedTplSubmesh(
            refined[level], refinedParent[level], recursivePhiRefinement);
    }
    mesh.calc3dmeshDef();

    // Reconnect every child to the final relaxed parent geometry, then refill
    // temperature with a forced remap so artificial cuts and phi-mids stay
    // coherent at the start of the first physical step.
    for (int level(0); level < recursiveSubmeshLevels + 1; level++) {
      refreshTplSubmeshArtificialGeometry(submesh[level],
                                          submeshParent[level]);
      deformTplSubmeshInterior(submesh[level], submeshParent[level]);
    }

    for (int in(0); in < mesh.number_of_nodes; in++) {
      mesh.node[in].num_of_adj_elem_l = 0;
      mesh.node[in].num_of_adj_elem_s = 0;
      for (int co(0); co < 3; co++) {
        mesh.node[in].dTl[co] = 0.;
        mesh.node[in].dTs[co] = 0.;
      }
    }
    for (int ie(0); ie < mesh.number_of_elements; ie++)
      local_calc.LocalGetGrad3d(&mesh.element[ie], &fe_step.gp_set_corn);
    for (int in(0); in < mesh.number_of_nodes; in++) {
      if (mesh.node[in].num_of_adj_elem_l != 0)
        for (int co(0); co < 3; co++)
          mesh.node[in].dTl[co] /= mesh.node[in].num_of_adj_elem_l;
      if (mesh.node[in].num_of_adj_elem_s != 0)
        for (int co(0); co < 3; co++)
          mesh.node[in].dTs[co] /= mesh.node[in].num_of_adj_elem_s;
    }
    for (int level(0); level < recursiveSubmeshLevels + 1; level++) {
      Mesh *p_sub = submesh[level];
      Mesh *p_parent = submeshParent[level];
      for (int in(0); in < p_sub->number_of_nodes; in++)
        p_sub->node[in].parent_mapping_initialized = 0;
      tplMapLocalNodes(p_sub, p_parent);
      transferTplSubmeshTemperature(p_sub, p_parent);
      calcTplSubmeshGradients(p_sub, &local_calc, &fe_step.gp_set_corn);
    }

    localSolidGradientOld.resize(3 * local_mesh.number_of_nodes, 0.);
    localLiquidGradientOld.resize(3 * local_mesh.number_of_nodes, 0.);
    refined1SolidGradientOld.resize(3 * refined_mesh_1.number_of_nodes, 0.);
    refined1LiquidGradientOld.resize(3 * refined_mesh_1.number_of_nodes, 0.);
    refined2SolidGradientOld.resize(3 * refined_mesh_2.number_of_nodes, 0.);
    refined2LiquidGradientOld.resize(3 * refined_mesh_2.number_of_nodes, 0.);
  }

  int steadyStateReached = 0;
  int finalIteration = 0;
  double finalTime = 0.;
  double finalNorm = 0.;
  int ssCheckCount = 0;
  int omp_max_threads = 1;
  int omp_num_procs = 1;
  int omp_probe_threads = 1;

#ifdef _OPENMP
  omp_max_threads = omp_get_max_threads();
  omp_num_procs = omp_get_num_procs();
  if (fe_step.use_openmp) {
#pragma omp parallel
    {
#pragma omp master
      omp_probe_threads = omp_get_num_threads();
    }
  }
#endif
  outw << "case | n=" << n
       << " | end_s=" << end_time
       << " | dt_s=" << physicalDtMax
       << " | move_substeps=" << nMoveSubsteps
       << " | move_dt_s=" << physicalDtMax / (double)nMoveSubsteps
       << " | dx_m=" << mesh.dx
       << " | GA=" << (mesh.useVoronkovGA ? "Voronkov" : "Linear")
       << "\n";
  outw << "mesh | submesh=" << useTplSubmesh
        << " | L0_layers=" << localSupportLayers
        << " | L0_radial=" << localLiquidRadialElements;
  if (useTplSubmesh)
    outw << " | L0_nodes=" << local_mesh.number_of_nodes
         << " | L0_elements=" << local_mesh.number_of_elements
         << " | levels=" << recursiveSubmeshLevels
         << " | phi_refine=" << recursivePhiRefinement
         << " | L0_SL=" << local_mesh.subSlRefinement
         << " | L1_layers=" << refined_mesh_1.tplSubmeshSupportLayers
         << " | L1_radial=" << refined_mesh_1.tplSubmeshOuterElements
         << " | L1_SL=" << refined_mesh_1.subSlRefinement
         << " | L2_layers=" << refined_mesh_2.tplSubmeshSupportLayers
         << " | L2_radial=" << refined_mesh_2.tplSubmeshOuterElements
         << " | L2_SL=" << refined_mesh_2.subSlRefinement
         << " | L3_layers=" << refined_mesh_3.tplSubmeshSupportLayers
         << " | L3_radial=" << refined_mesh_3.tplSubmeshOuterElements
         << " | L3_SL=" << refined_mesh_3.subSlRefinement
         << " | L1_nPhi=" << refined_mesh_1.nTPL
         << " | L2_nPhi=" << refined_mesh_2.nTPL
         << " | L3_nPhi=" << refined_mesh_3.nTPL;
  outw << "\n";
  outw << "state | input=" << ssStateInputFile
       << " | restart=" << restartFromStoredState
       << " | axisymmetric=" << initializedFromAxisymmetricState << "\n";
  outw << "thermal | top=" << fe_step.local_calc.Ttop
       << " | bottom=" << fe_step.local_calc.Tbottom
       << " | out=" << fe_step.local_calc.Tout << "\n";
  outw << "meniscus | h0=" << fe_step.local_calc.h0
       << " | lc=" << lc
       << " | hSlab=" << hSlab
       << " | Rc=" << Rc
       << " | GAr=" << GAr << "\n";
  outw << "kinetics | betaR=" << mesh.betaRough
       << " | betaS=" << mesh.betaStep
       << " | SLrelax="
       << mesh.bc_set[mesh.BC_SOLID_LIQUID].relaxMoveCoef
       << " | A2DN=" << mesh.kineticA2DN
       << " | B2DN=" << mesh.kineticB2DN << "\n";
  outw << "threads | active=" << omp_probe_threads
       << " | max=" << omp_max_threads
       << " | requested=" << omp_threads
       << " | procs=" << omp_num_procs
       << " | slots=" << fe_step.thread_calcs.size() << "\n";

  mesh.anisotropicKinFlg = 1;
  mesh.AniGaFlg = 1;
  mesh.bc_set[mesh.BC_SOLID_LIQUID].calcFacetAngle();
  if (useTplSubmesh) {
    local_mesh.anisotropicKinFlg = mesh.anisotropicKinFlg;
    local_mesh.AniGaFlg = mesh.AniGaFlg;
    local_mesh.bc_set[local_mesh.BC_SUB_SOLID_LIQUID].calcFacetAngle();
    refined_mesh_1.anisotropicKinFlg = mesh.anisotropicKinFlg;
    refined_mesh_1.AniGaFlg = mesh.AniGaFlg;
    refined_mesh_1.bc_set[refined_mesh_1.BC_SUB_SOLID_LIQUID]
        .calcFacetAngle();
    refined_mesh_2.anisotropicKinFlg = mesh.anisotropicKinFlg;
    refined_mesh_2.AniGaFlg = mesh.AniGaFlg;
    refined_mesh_2.bc_set[refined_mesh_2.BC_SUB_SOLID_LIQUID]
        .calcFacetAngle();
    refined_mesh_3.anisotropicKinFlg = mesh.anisotropicKinFlg;
    refined_mesh_3.AniGaFlg = mesh.AniGaFlg;
    refined_mesh_3.bc_set[refined_mesh_3.BC_SUB_SOLID_LIQUID]
        .calcFacetAngle();
  }
  if (mesh.useVoronkovGA)
    mesh.bc_set[mesh.BC_TPL].MoveTPLvoronkovGA(
        0., GAr, GAf, mesh.AniGaFlg, &mesh);
  else
    mesh.bc_set[mesh.BC_TPL].MoveTPL(
        0., GAr, GAf, mesh.AniGaFlg, &mesh);

  {
    char filename[100];
    const int ttOut = 0;
    sprintf(filename, "%s/dom%.6d.dat", outputDir, ttOut);
    mesh.writeDomBRICK_Tec(filename);
    sprintf(filename, "%s/BCcrysMelt%.6d.dat", outputDir, ttOut);
    mesh.writeBcQuad_Tec(&mesh.bc_set[mesh.BC_SOLID_LIQUID], filename);
    sprintf(filename, "%s/BCgasMelt%.6d.dat", outputDir, ttOut);
    mesh.bc_set[mesh.BC_LIQUID_GAS].writeIstr(filename);
    sprintf(filename, "%s/BCgasSolid%.6d.dat", outputDir, ttOut);
    mesh.bc_set[mesh.BC_SOLID_GAS].writeIstr(filename);
    sprintf(filename, "%s/BCtpl%.6d.dat", outputDir, ttOut);
    mesh.bc_set[mesh.BC_TPL].writeTPL(filename);
    if (useTplSubmesh) {
      Mesh *submesh[4] = {
          &local_mesh, &refined_mesh_1, &refined_mesh_2, &refined_mesh_3};
      for (int level(0); level < recursiveSubmeshLevels + 1; level++) {
        Mesh *p_sub = submesh[level];
        if (p_sub->useVoronkovGA)
          p_sub->bc_set[p_sub->BC_SUB_TPL].MoveTPLvoronkovGA(
              0., GAr, GAf, p_sub->AniGaFlg, p_sub);
        else
          p_sub->bc_set[p_sub->BC_SUB_TPL].MoveTPL(
              0., GAr, GAf, p_sub->AniGaFlg, p_sub);
        sprintf(filename, "%s/sub%d_before_Dom.dat", outputDir, level + 1);
        p_sub->writeDomBRICK_Tec(filename);
        sprintf(filename, "%s/sub%d_before_SL.dat", outputDir, level + 1);
        p_sub->writeBcQuad_Tec(
            &p_sub->bc_set[p_sub->BC_SUB_SOLID_LIQUID], filename);
        sprintf(filename, "%s/sub%d_before_LG.dat", outputDir, level + 1);
        p_sub->bc_set[p_sub->BC_SUB_LIQUID_GAS].writeIstr(filename);
        sprintf(filename, "%s/sub%d_before_SG.dat", outputDir, level + 1);
        p_sub->bc_set[p_sub->BC_SUB_SOLID_GAS].writeIstr(filename);
        sprintf(filename, "%s/sub%d_before_TPL.dat", outputDir, level + 1);
        p_sub->bc_set[p_sub->BC_SUB_TPL].writeTPL(filename);
      }
    }
  }

  outw << "transient\n";
  // The first transient row is written after the first complete thermal and
  // interface step, so its kinetic velocities are evaluated state values.
  double lastTableOutputTime = -1.;

  {
    double time_value = 0.;
    double nextSsCheck = ssCheckInterval;
    double nextOutput = output_interval;
    const double dtMoveCoef = 1.0 / (double)nMoveSubsteps;

    for (int it(1);
         time_value < end_time && (end_steps < 0 || it <= end_steps);
         it++) {
      dt = physicalDtMax;
      if (time_value + dt > end_time)
        dt = end_time - time_value;
      fe_step.dt = dt;
      time_value += dt;
      finalIteration = it;
      finalTime = time_value;

      if (!useTplSubmesh && it % mesh_update_every == 0)
        mesh.calc3dmeshDef();
      if (useTplSubmesh)
        for (int in(0); in < mesh.number_of_nodes; in++)
          for (int co(0); co < 3; co++) {
            coarseSolidGradientOld[3 * in + co] = mesh.node[in].dTs[co];
            coarseLiquidGradientOld[3 * in + co] = mesh.node[in].dTl[co];
          }
      if (useTplSubmesh)
        for (int in(0); in < local_mesh.number_of_nodes; in++)
          for (int co(0); co < 3; co++) {
            localSolidGradientOld[3 * in + co] =
                local_mesh.node[in].dTs[co];
            localLiquidGradientOld[3 * in + co] =
                local_mesh.node[in].dTl[co];
          }
      if (useTplSubmesh)
        for (int in(0); in < refined_mesh_1.number_of_nodes; in++)
          for (int co(0); co < 3; co++) {
            refined1SolidGradientOld[3 * in + co] =
                refined_mesh_1.node[in].dTs[co];
            refined1LiquidGradientOld[3 * in + co] =
                refined_mesh_1.node[in].dTl[co];
          }
      if (useTplSubmesh)
        for (int in(0); in < refined_mesh_2.number_of_nodes; in++)
          for (int co(0); co < 3; co++) {
            refined2SolidGradientOld[3 * in + co] =
                refined_mesh_2.node[in].dTs[co];
            refined2LiquidGradientOld[3 * in + co] =
                refined_mesh_2.node[in].dTl[co];
          }
      double norm = runFEtemperatureStep(&mesh, &fe_step);
      finalNorm = norm;

      if (useTplSubmesh) {
        int thermalSubsteps = 0;
        // Parent temperature and artificial-cut geometry flow downward.  The
        // physical interface geometry is retained until it is restricted
        // upward after movement.
        refreshTplSubmeshArtificialGeometry(&local_mesh, &mesh);
        deformTplSubmeshInterior(&local_mesh, &mesh);
        // Advance the local thermal field over the same physical interval as
        // the coarse FE step. The current one-step implicit solve uses the
        // new-time parent field on the artificial cut.
        solveTplSubmeshTemperature(
            &local_mesh, &mesh, &local_calc,
            &fe_step.gp_set[BRICK_8_NOD_LAGR],
            &fe_step.gp_set[QUADRATIC_4_NOD_LAGR],
            dt, coarseSolidGradientOld.data(),
            coarseLiquidGradientOld.data(),
            &thermalSubsteps);
        calcTplSubmeshGradients(&local_mesh, &local_calc,
                                &fe_step.gp_set_corn);

        refreshTplSubmeshArtificialGeometry(&refined_mesh_1, &local_mesh);
        deformTplSubmeshInterior(&refined_mesh_1, &local_mesh);
        solveTplSubmeshTemperature(
            &refined_mesh_1, &local_mesh, &local_calc,
            &fe_step.gp_set[BRICK_8_NOD_LAGR],
            &fe_step.gp_set[QUADRATIC_4_NOD_LAGR],
            dt, localSolidGradientOld.data(),
            localLiquidGradientOld.data(),
            &thermalSubsteps);
        calcTplSubmeshGradients(&refined_mesh_1, &local_calc,
                                &fe_step.gp_set_corn);

        refreshTplSubmeshArtificialGeometry(&refined_mesh_2,
                                            &refined_mesh_1);
        deformTplSubmeshInterior(&refined_mesh_2, &refined_mesh_1);
        solveTplSubmeshTemperature(
            &refined_mesh_2, &refined_mesh_1, &local_calc,
            &fe_step.gp_set[BRICK_8_NOD_LAGR],
            &fe_step.gp_set[QUADRATIC_4_NOD_LAGR],
            dt, refined1SolidGradientOld.data(),
            refined1LiquidGradientOld.data(),
            &thermalSubsteps);
        calcTplSubmeshGradients(&refined_mesh_2, &local_calc,
                                &fe_step.gp_set_corn);

        refreshTplSubmeshArtificialGeometry(&refined_mesh_3,
                                            &refined_mesh_2);
        deformTplSubmeshInterior(&refined_mesh_3, &refined_mesh_2);
        solveTplSubmeshTemperature(
            &refined_mesh_3, &refined_mesh_2, &local_calc,
            &fe_step.gp_set[BRICK_8_NOD_LAGR],
            &fe_step.gp_set[QUADRATIC_4_NOD_LAGR],
            dt, refined2SolidGradientOld.data(),
            refined2LiquidGradientOld.data(),
            &thermalSubsteps);
        calcTplSubmeshGradients(&refined_mesh_3, &local_calc,
                                &fe_step.gp_set_corn);

        // Node::T remains the old parent time level while each child is
        // solved. Advance it only after all recursive transfers finish.
        for (int in(0); in < local_mesh.number_of_nodes; in++)
          local_mesh.node[in].T =
              local_mesh.node[in].attribute[TEMPERATURE];
        for (int in(0); in < refined_mesh_1.number_of_nodes; in++)
          refined_mesh_1.node[in].T =
              refined_mesh_1.node[in].attribute[TEMPERATURE];
        for (int in(0); in < refined_mesh_2.number_of_nodes; in++)
          refined_mesh_2.node[in].T =
              refined_mesh_2.node[in].attribute[TEMPERATURE];
        for (int in(0); in < refined_mesh_3.number_of_nodes; in++)
          refined_mesh_3.node[in].T =
              refined_mesh_3.node[in].attribute[TEMPERATURE];
      }

      if (useTplSubmesh) {
        double physicalTimeMoved = 0.;
        while (physicalTimeMoved < dt) {
          double physicalSubstep =
              min(dtMoveCoef * dt, dt - physicalTimeMoved);
          // TPL and L/G motion is accelerated numerical relaxation. Scale it
          // with the physical substep so changing the thermal/kinetic step
          // does not change its relaxation rate per unit simulated time.
          const double geometrySubstep =
              min(.05 * physicalSubstep, maxGeometrySubstep);
          // Keep the TPL update explicit in the physical main loop. The
          // remaining helper receives geometryDt=0 below, so it cannot repeat
          // either this TPL motion or the adjacent L/G normal relaxation.
          // Parent (coarse) interfaces, then 3d, then update child BC and
          // child repeats: interfaces -> 3d -> next child. Do not move
          // children before the parent finishes this substep.
          if (mesh.useVoronkovGA)
            mesh.bc_set[mesh.BC_TPL].MoveTPLvoronkovGA(
                geometrySubstep, GAr, GAf, mesh.AniGaFlg, &mesh);
          else
            mesh.bc_set[mesh.BC_TPL].MoveTPL(
                geometrySubstep, GAr, GAf, mesh.AniGaFlg, &mesh);
          mesh.bc_set[mesh.BC_LIQUID_GAS].MoveLiquidInterface(
              geometrySubstep, &mesh);
          mesh.bc_set[mesh.BC_SOLID_LIQUID].calcRelaxSL();
          mesh.bc_set[mesh.BC_SOLID_LIQUID].kinMoveInterface(
              physicalSubstep, &mesh);
          mesh.bc_set[mesh.BC_SOLID_GAS].calcRelaxStr();
          mesh.bc_set[mesh.BC_LIQUID_GAS].calcRelaxStr();
          mesh.bc_set[mesh.BC_SOLID_GAS].pullSGstr(
              physicalSubstep * mesh.Vpull);
          if ((int)(physicalTimeMoved / physicalSubstep + .5) %
                  mesh_update_every ==
              0)
            mesh.calc3dmeshDef();

          Mesh *levelMesh[4] = {&local_mesh, &refined_mesh_1,
                                &refined_mesh_2, &refined_mesh_3};
          Mesh *levelParent[4] = {&mesh, &local_mesh, &refined_mesh_1,
                                  &refined_mesh_2};
          for (int level(0); level < recursiveSubmeshLevels + 1; level++) {
            Mesh *p_sub = levelMesh[level];
            Mesh *p_parent = levelParent[level];
            refreshTplSubmeshArtificialGeometry(p_sub, p_parent);
            tplCopyPhysics(p_sub, p_parent);
            if (p_sub->useVoronkovGA)
              p_sub->bc_set[p_sub->BC_SUB_TPL].MoveTPLvoronkovGA(
                  geometrySubstep, GAr, GAf, p_sub->AniGaFlg, p_sub);
            else
              p_sub->bc_set[p_sub->BC_SUB_TPL].MoveTPL(
                  geometrySubstep, GAr, GAf, p_sub->AniGaFlg, p_sub);
            p_sub->bc_set[p_sub->BC_SUB_LIQUID_GAS].MoveLiquidInterface(
                geometrySubstep, p_sub);
            moveTplSubmeshInterfaces(p_sub, p_parent, physicalSubstep, 0.,
                                     GAr, GAf);
            deformTplSubmeshInterior(p_sub, p_parent);
          }
          physicalTimeMoved += physicalSubstep;
        }
      } else {
        for (int s(0); s < nMoveSubsteps; s++) {
          if (mesh.useVoronkovGA)
            mesh.bc_set[mesh.BC_TPL].MoveTPLvoronkovGA(
                dtMoveCoef * dt, GAr, GAf, mesh.AniGaFlg, &mesh);
          else
            mesh.bc_set[mesh.BC_TPL].MoveTPL(
                dtMoveCoef * dt, GAr, GAf, mesh.AniGaFlg, &mesh);
          mesh.bc_set[mesh.BC_LIQUID_GAS].MoveLiquidInterface(
              dtMoveCoef * dt, &mesh);
          mesh.bc_set[mesh.BC_SOLID_LIQUID].calcRelaxSL();
          mesh.bc_set[mesh.BC_SOLID_LIQUID].kinMoveInterface(
              dtMoveCoef * dt, &mesh);
          mesh.bc_set[mesh.BC_SOLID_GAS].calcRelaxStr();
          mesh.bc_set[mesh.BC_LIQUID_GAS].calcRelaxStr();
          mesh.bc_set[mesh.BC_SOLID_GAS].pullSGstr(
              dtMoveCoef * dt * mesh.Vpull);
        }
      }

      if (useTplSubmesh) {
        // Geometry is restricted finest to coarsest.  Temperature remains on
        // the level where it was solved and is never passed upward.
        restrictTplSubmeshToCoarse(&refined_mesh_3, &refined_mesh_2);
        deformTplSubmeshInterior(&refined_mesh_2, &refined_mesh_1);
        restrictTplSubmeshToCoarse(&refined_mesh_2, &refined_mesh_1);
        deformTplSubmeshInterior(&refined_mesh_1, &local_mesh);
        restrictTplSubmeshToCoarse(&refined_mesh_1, &local_mesh);
        deformTplSubmeshInterior(&local_mesh, &mesh);
        restrictTplSubmeshToCoarse(&local_mesh, &mesh);
        if (it % mesh_update_every == 0)
          mesh.calc3dmeshDef();
        // Reconnect only the parent-controlled bulk cuts from parent to child.
        // Copying all coincident nodes here would erase the fine interface.
        refreshTplSubmeshArtificialGeometry(&local_mesh, &mesh);
        deformTplSubmeshInterior(&local_mesh, &mesh);
        refreshTplSubmeshArtificialGeometry(&refined_mesh_1, &local_mesh);
        deformTplSubmeshInterior(&refined_mesh_1, &local_mesh);
        refreshTplSubmeshArtificialGeometry(&refined_mesh_2,
                                            &refined_mesh_1);
        deformTplSubmeshInterior(&refined_mesh_2, &refined_mesh_1);
        refreshTplSubmeshArtificialGeometry(&refined_mesh_3,
                                            &refined_mesh_2);
        deformTplSubmeshInterior(&refined_mesh_3, &refined_mesh_2);
      }

      if (useTplSubmesh && time_value + .5 * dt >= nextSsCheck) {
        BCSet *p_ss_sl =
            &refined_mesh_3.bc_set[refined_mesh_3.BC_SUB_SOLID_LIQUID];
        BCAtNode *p_rough = p_ss_sl->bc_at_nodeStr[0][0];
        BCAtNode *p_facet =
            p_ss_sl->bc_at_nodeStr[p_ss_sl->nS1 / 4][0];
        const double roughVlab = fabs(p_rough->dcoord[2]);
        const double facetVlab = fabs(p_facet->dcoord[2]);
        if (roughVlab < ssVlabTolerance &&
            facetVlab < ssVlabTolerance)
          ssCheckCount++;
        else
          ssCheckCount = 0;
        if (ssCheckCount >= ssRequiredChecks) {
          steadyStateReached = 1;
          outw << "steady_state | t=" << time_value
               << " | rough_Vz=" << roughVlab * 1000. * 60.
               << " | facet_Vz=" << facetVlab * 1000. * 60. << "\n";
        }
        nextSsCheck += ssCheckInterval;
      }

      const int scheduledOutput =
          (time_value + .5 * dt >= nextOutput);
      const int finalRequestedStep =
          (time_value + .5 * dt >= end_time) ||
          (end_steps > 0 && it >= end_steps);
      const int writeFullOutput = scheduledOutput || finalRequestedStep;
      if (it == 1 || writeFullOutput) {
        writeDriverOutput(&mesh,
                          useTplSubmesh ? &refined_mesh_3 : &mesh,
                          outw, time_value, norm, GAr, GAf,
                          outputDir, 0, writeFullOutput, it == 1);
        lastTableOutputTime = time_value;
        if (writeFullOutput && useTplSubmesh) {
          char filename[100];
          const int ttOut = (int)floor(time_value + .555555);
          sprintf(filename, "%s/sub4_Dom%.6d.dat", outputDir, ttOut);
          refined_mesh_3.writeDomBRICK_Tec(filename);
          sprintf(filename, "%s/sub4_SL%.6d.dat", outputDir, ttOut);
          refined_mesh_3.writeBcQuad_Tec(
              &refined_mesh_3.bc_set[refined_mesh_3.BC_SUB_SOLID_LIQUID],
              filename);
          sprintf(filename, "%s/sub4_LG%.6d.dat", outputDir, ttOut);
          refined_mesh_3.bc_set[refined_mesh_3.BC_SUB_LIQUID_GAS]
              .writeIstr(filename);
          sprintf(filename, "%s/sub4_SG%.6d.dat", outputDir, ttOut);
          refined_mesh_3.bc_set[refined_mesh_3.BC_SUB_SOLID_GAS]
              .writeIstr(filename);
          sprintf(filename, "%s/sub4_TPL%.6d.dat", outputDir, ttOut);
          if (refined_mesh_3.useVoronkovGA)
            refined_mesh_3.bc_set[refined_mesh_3.BC_SUB_TPL]
                .MoveTPLvoronkovGA(
                    0., GAr, GAf, refined_mesh_3.AniGaFlg, &refined_mesh_3);
          else
            refined_mesh_3.bc_set[refined_mesh_3.BC_SUB_TPL].MoveTPL(
                0., GAr, GAf, refined_mesh_3.AniGaFlg, &refined_mesh_3);
          refined_mesh_3.bc_set[refined_mesh_3.BC_SUB_TPL]
              .writeTPL(filename);
        }
        if (scheduledOutput)
          nextOutput += output_interval;
      }
      if (steadyStateReached)
        break;
    }
  }

  if (finalIteration > 0 &&
      fabs(finalTime - lastTableOutputTime) > .5 * physicalDtMax) {
    writeDriverOutput(&mesh,
                      useTplSubmesh ? &refined_mesh_3 : &mesh,
                      outw, finalTime, finalNorm, GAr, GAf,
                      outputDir, 0, 1, 0);
    lastTableOutputTime = finalTime;
  }

  if (useTplSubmesh && finalIteration > 0) {
    if (refined_mesh_3.useVoronkovGA)
      refined_mesh_3.bc_set[refined_mesh_3.BC_SUB_TPL]
          .MoveTPLvoronkovGA(
              0., GAr, GAf, refined_mesh_3.AniGaFlg, &refined_mesh_3);
    else
      refined_mesh_3.bc_set[refined_mesh_3.BC_SUB_TPL].MoveTPL(
          0., GAr, GAf, refined_mesh_3.AniGaFlg, &refined_mesh_3);

    BCSet *p_tpl =
        &refined_mesh_3.bc_set[refined_mesh_3.BC_SUB_TPL];
    BCSet *p_sl =
        &refined_mesh_3.bc_set[refined_mesh_3.BC_SUB_SOLID_LIQUID];
    const int nPhi = p_tpl->nS1;
    const int tip = nPhi / 4;
    // The tip is the structured symmetry node nPhi/4.  Starting from this
    // known point, the two ridge bases are the first nodes where the local
    // undercooling reaches the rough steady kinetic value Vpull/betaRough.
    const double dTrough =
        fabs(refined_mesh_3.Vpull) / refined_mesh_3.betaRough;
    BCAtNode *p_tip = p_tpl->bc_at_nodeStr[tip][1];
    BCAtNode *p_tip_sl =
        &p_sl->bc_at_node[p_tip->parent_bc_node_id];
    const double dTtip =
        -p_tip_sl->p_node->attribute[TEMPERATURE];

    int leftBase = tip;
    int rightBase = tip;
    int leftSteps = 0;
    int rightSteps = 0;
    for (int step(1); step <= nPhi / 4; step++) {
      const int i = (tip - step + nPhi) % nPhi;
      BCAtNode *p_contact = p_tpl->bc_at_nodeStr[i][1];
      BCAtNode *p_node_sl =
          &p_sl->bc_at_node[p_contact->parent_bc_node_id];
      if (-p_node_sl->p_node->attribute[TEMPERATURE] <= dTrough) {
        leftBase = i;
        leftSteps = step;
        break;
      }
    }
    for (int step(1); step <= nPhi / 4; step++) {
      const int i = (tip + step) % nPhi;
      BCAtNode *p_contact = p_tpl->bc_at_nodeStr[i][1];
      BCAtNode *p_node_sl =
          &p_sl->bc_at_node[p_contact->parent_bc_node_id];
      if (-p_node_sl->p_node->attribute[TEMPERATURE] <= dTrough) {
        rightBase = i;
        rightSteps = step;
        break;
      }
    }

    const int ridgeResolved = leftSteps != 0 && rightSteps != 0 &&
                              dTtip > dTrough;
    double leftArc = 0.;
    double rightArc = 0.;
    if (ridgeResolved) {
      int previous = tip;
      for (int step(1); step <= leftSteps; step++) {
        const int i = (tip - step + nPhi) % nPhi;
        Node *p0 = p_tpl->bc_at_nodeStr[previous][1]->p_node;
        Node *p1 = p_tpl->bc_at_nodeStr[i][1]->p_node;
        const double dx = p1->coord[0] - p0->coord[0];
        const double dy = p1->coord[1] - p0->coord[1];
        const double dz = p1->coord[2] - p0->coord[2];
        leftArc += sqrt(dx * dx + dy * dy + dz * dz);
        previous = i;
      }
      previous = tip;
      for (int step(1); step <= rightSteps; step++) {
        const int i = (tip + step) % nPhi;
        Node *p0 = p_tpl->bc_at_nodeStr[previous][1]->p_node;
        Node *p1 = p_tpl->bc_at_nodeStr[i][1]->p_node;
        const double dx = p1->coord[0] - p0->coord[0];
        const double dy = p1->coord[1] - p0->coord[1];
        const double dz = p1->coord[2] - p0->coord[2];
        rightArc += sqrt(dx * dx + dy * dy + dz * dz);
        previous = i;
      }
    }

    ios::fmtflags ridgeFlags = outw.flags();
    streamsize ridgePrecision = outw.precision();
    if (ridgeResolved) {
      outw << "ridge final | K, mm, deg\n"
           << left
           << setw(7) << "kind" << "|"
           << setw(8) << "dTtip" << "|"
           << setw(8) << "dTbase" << "|"
           << setw(8) << "dr" << "|"
           << setw(8) << "dz" << "|"
           << setw(8) << "l" << "|"
           << setw(8) << "w" << "|"
           << setw(8) << "chi" << "\n";
      BCAtNode *p_left = p_tpl->bc_at_nodeStr[leftBase][1];
      BCAtNode *p_right = p_tpl->bc_at_nodeStr[rightBase][1];
      BCAtNode *p_left_sl =
          &p_sl->bc_at_node[p_left->parent_bc_node_id];
      BCAtNode *p_right_sl =
          &p_sl->bc_at_node[p_right->parent_bc_node_id];
      Node *p_tip_node = p_tip_sl->p_node;
      Node *p_left_node = p_left_sl->p_node;
      Node *p_right_node = p_right_sl->p_node;
      const double rTip = sqrt(p_tip_node->coord[0] * p_tip_node->coord[0] +
                               p_tip_node->coord[1] * p_tip_node->coord[1]);
      const double rLeft = sqrt(p_left_node->coord[0] * p_left_node->coord[0] +
                                p_left_node->coord[1] * p_left_node->coord[1]);
      const double rRight = sqrt(
          p_right_node->coord[0] * p_right_node->coord[0] +
          p_right_node->coord[1] * p_right_node->coord[1]);
      const double drMesh = rTip - .5 * (rLeft + rRight);
      const double dzMesh = p_tip_node->coord[2] -
                            .5 * (p_left_node->coord[2] +
                                  p_right_node->coord[2]);
      const double lMesh = sqrt(drMesh * drMesh + dzMesh * dzMesh);
      const double dTbase = -.5 *
          (p_left_node->attribute[TEMPERATURE] +
           p_right_node->attribute[TEMPERATURE]);
      const double gaBase = .5 * (p_left->GA + p_right->GA);
      const double gaAppBase = .5 * (p_left->GAv + p_right->GAv);
      const double chiMesh = fabs(p_tip->GA - gaBase);
      const double chiTheory = fabs(p_tip->GAv - gaAppBase);

      outw << left << setw(7) << "mesh" << "|" << right
           << defaultfloat << setprecision(4)
           << setw(8) << dTtip << "|"
           << setw(8) << dTbase << "|"
           << setw(8) << 1000. * drMesh << "|"
           << setw(8) << 1000. * dzMesh << "|"
           << setw(8) << 1000. * lMesh << "|"
           << setw(8) << 1000. * (leftArc + rightArc) << "|"
           << fixed << setprecision(2)
           << setw(8) << 180. * chiMesh / PI << "\n";

      const double facetSlopeFromVertical =
          fabs(refined_mesh_3.voronkovThetaFacet);
      const double rxy = sqrt(
          p_tip_node->coord[0] * p_tip_node->coord[0] +
          p_tip_node->coord[1] * p_tip_node->coord[1]);
      const double GrTip =
          (p_tip_node->dTs[0] * p_tip_node->coord[0] +
           p_tip_node->dTs[1] * p_tip_node->coord[1]) / rxy;
      const double Gfacet = fabs(
          GrTip * sin(facetSlopeFromVertical) +
          p_tip_node->dTs[2] * cos(facetSlopeFromVertical));
      const double dTridge = dTtip - dTrough;
      if (Gfacet > 0. && chiTheory > 0.) {
        const double lTheory = dTridge / Gfacet;
        const double drTheory = lTheory * sin(facetSlopeFromVertical);
        const double dzTheory = lTheory * cos(facetSlopeFromVertical);
        const double wTheory = 4. * fabs(drTheory) / chiTheory;
        outw << left << setw(7) << "theory" << "|" << right
             << defaultfloat << setprecision(4)
             << setw(8) << dTtip << "|"
             << setw(8) << dTrough << "|"
             << setw(8) << 1000. * drTheory << "|"
             << setw(8) << 1000. * dzTheory << "|"
             << setw(8) << 1000. * lTheory << "|"
             << setw(8) << 1000. * wTheory << "|"
             << fixed << setprecision(2)
             << setw(8) << 180. * chiTheory / PI << "\n";
      }
    }
    outw.flags(ridgeFlags);
    outw.precision(ridgePrecision);
  }

  char finalStateOutputFile[120];
  sprintf(finalStateOutputFile, "%s/final_state_n%04d.dat", outputDir, n);
  writeMeshState(finalStateOutputFile, &mesh, &local_calc,
                 n, d_crucible, finalTime, GAf, GAr, Rc, lc, hSlab);

  if (useTplSubmesh) {
    char filename[140];
    Mesh *submesh[4] = {
        &local_mesh, &refined_mesh_1, &refined_mesh_2, &refined_mesh_3};
    int finalStateComplete = 1;
    for (int level(0); level < recursiveSubmeshLevels + 1; level++) {
      Mesh *p_sub = submesh[level];
      if (p_sub->useVoronkovGA)
        p_sub->bc_set[p_sub->BC_SUB_TPL].MoveTPLvoronkovGA(
            0., GAr, GAf, p_sub->AniGaFlg, p_sub);
      else
        p_sub->bc_set[p_sub->BC_SUB_TPL].MoveTPL(
            0., GAr, GAf, p_sub->AniGaFlg, p_sub);

      sprintf(filename, "%s/sub%d_after_Dom.dat", outputDir, level + 1);
      p_sub->writeDomBRICK_Tec(filename);
      sprintf(filename, "%s/sub%d_after_SL.dat", outputDir, level + 1);
      p_sub->writeBcQuad_Tec(
          &p_sub->bc_set[p_sub->BC_SUB_SOLID_LIQUID], filename);
      sprintf(filename, "%s/sub%d_after_LG.dat", outputDir, level + 1);
      p_sub->bc_set[p_sub->BC_SUB_LIQUID_GAS].writeIstr(filename);
      sprintf(filename, "%s/sub%d_after_SG.dat", outputDir, level + 1);
      p_sub->bc_set[p_sub->BC_SUB_SOLID_GAS].writeIstr(filename);
      sprintf(filename, "%s/sub%d_after_TPL.dat", outputDir, level + 1);
      p_sub->bc_set[p_sub->BC_SUB_TPL].writeTPL(filename);
      sprintf(filename, "%s/sub%d_state.dat", outputDir, level + 1);
      finalStateComplete =
          writeTplSubmeshState(filename, p_sub, level + 1) &&
          finalStateComplete;
    }
    int benchmarkHalfColumns = ridgeBenchmarkHalfCoarseIntervals;
    for (int level(0); level < recursiveSubmeshLevels; level++)
      benchmarkHalfColumns *= recursivePhiRefinement;
    sprintf(filename, "%s/sub4_ridge_thermal_bridge.dat", outputDir);
    finalStateComplete =
        writeTplCircularBenchmark(filename, &refined_mesh_3,
                                  recursiveSubmeshLevels,
                                  benchmarkHalfColumns, outw) &&
        finalStateComplete;
    if (finalStateComplete) {
      const string manifestName =
          string(outputDir) + "/submesh_final_manifest.txt";
      const string temporaryName = manifestName + ".tmp";
      ofstream manifest(temporaryName.c_str());
      manifest << "FEEXP_TPL_SUBMESH_FINAL_V4\n"
               << "recursive_levels " << recursiveSubmeshLevels << "\n"
               << "steady_state " << steadyStateReached << "\n"
               << "final_time_s " << finalTime << "\n"
               << "sub1_nPhi " << local_mesh.nTPL << "\n"
               << "sub2_nPhi " << refined_mesh_1.nTPL << "\n"
               << "sub3_nPhi " << refined_mesh_2.nTPL << "\n"
               << "sub4_nPhi " << refined_mesh_3.nTPL << "\n"
               << "state sub1_state.dat\n"
               << "state sub2_state.dat\n"
               << "state sub3_state.dat\n"
               << "state sub4_state.dat\n"
               << "surfaces sub1_after_*.dat\n"
               << "surfaces sub2_after_*.dat\n"
               << "surfaces sub3_after_*.dat\n"
               << "surfaces sub4_after_*.dat\n"
               << "bridge sub4_ridge_thermal_bridge.dat\n"
               << "COMPLETE\n";
      manifest.close();
      if (!manifest || rename(temporaryName.c_str(), manifestName.c_str()) != 0) {
        remove(manifestName.c_str());
        if (rename(temporaryName.c_str(), manifestName.c_str()) != 0)
          outw << "Failed to publish recursive submesh final manifest\n";
      }
    } else {
      outw << "Recursive submesh final state is incomplete\n";
    }
  }

  outw << "run_end | t=" << finalTime
       << " | steady=" << steadyStateReached
       << " | state=" << finalStateOutputFile << "\n";

  t = clock() - t;
  outw << "end runtime_min=" << ((float)t / 60.) / CLOCKS_PER_SEC << "\n";
  cout.rdbuf(p_cout_buf);
  cerr.rdbuf(p_cerr_buf);
  return 0;
}

static int getEnvInt(const char *name, int fallback) {
  const char *value = getenv(name);
  if (value == NULL || value[0] == '\0')
    return fallback;

  char *end_ptr = NULL;
  long parsed = strtol(value, &end_ptr, 10);
  if (end_ptr == value)
    return fallback;

  return (int)parsed;
}

static void writeMeshState(const char *filename, Mesh *p_mesh,
                           LocalCalc *p_local_calc,
                           int n, double d_crucible, double end_time,
                           double GAf, double GAr,
                           double Rc, double lc, double hSlab) {
  FILE *ifwri = fopen(filename, "w");
  if (ifwri == NULL) {
    cerr << "Failed to open mesh state file: " << filename << "\n";
    return;
  }

  fprintf(ifwri, "FEEXP_SS_MESH_STATE_V2\n");
  fprintf(ifwri, "N %d\n", n);
  fprintf(ifwri, "COUNTS %d %d\n",
          p_mesh->number_of_nodes,
          p_mesh->number_of_node_attributes);
  fprintf(ifwri,
          "GEOM %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
          d_crucible, end_time, GAf, GAr, Rc, lc, hSlab, p_local_calc->h0);
  fprintf(ifwri,
          "LOCAL %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
          p_local_calc->Ttop, p_local_calc->Tout, p_local_calc->Tmp,
          p_local_calc->Tbottom, p_local_calc->hh, p_local_calc->Rad,
          p_local_calc->PeS, p_local_calc->cond[Mesh::SOLID],
          p_local_calc->cond[Mesh::LIQUID], p_local_calc->Vp);
  fprintf(ifwri,
          "MESH %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
          p_mesh->dx, p_mesh->Vpull, p_mesh->betaRough, p_mesh->betaStep,
          p_mesh->kineticA2DN, p_mesh->kineticB2DN,
          0., 0., 0., 0., p_mesh->rhog, p_mesh->gam);
  fprintf(ifwri,
          "VOR %.17e %.17e %.17e %.17e %d %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
          p_mesh->voronkovTm, p_mesh->voronkovQ,
          p_mesh->voronkovLambdaSG, p_mesh->voronkovAtomicDensity,
          p_mesh->voronkovTemperatureIsAbsolute,
          p_mesh->voronkovAlphaTransition,
          p_mesh->voronkovSigmaSLFacet, p_mesh->voronkovSigmaSLRough,
          p_mesh->voronkovSigmaSLPrimeFacet, p_mesh->voronkovSigmaSG,
          p_mesh->voronkovSigmaMG, p_mesh->voronkovThetaFacet,
          p_mesh->voronkovThetaStepCorrection,
          p_mesh->voronkovVelocityToCm, p_mesh->voronkovKV);
  fprintf(ifwri, "NODES\n");

  for (int i(0); i < p_mesh->number_of_nodes; i++) {
    fprintf(ifwri, "%d", i);
    for (int co(0); co < 3; co++)
      fprintf(ifwri, " %.17e", p_mesh->node[i].coord[co]);
    for (int a(0); a < p_mesh->number_of_node_attributes; a++)
      fprintf(ifwri, " %.17e", p_mesh->node[i].attribute[a]);
    fprintf(ifwri, "\n");
  }

  fclose(ifwri);
}

static int readMeshState(const char *filename, Mesh *p_mesh,
                         LocalCalc *p_local_calc,
                         int n, double d_crucible,
                         double GAf, double GAr,
                         ostream &outw) {
  FILE *ifre = fopen(filename, "r");
  if (ifre == NULL) {
    outw << "Failed to open mesh state file: " << filename << "\n";
    return 0;
  }

  char tag[80];
  char key[80];
  int file_n = -1;
  int n_nodes = 0;
  int n_attr = 0;
  if (fscanf(ifre, "%79s", tag) != 1) {
    outw << "Bad mesh state header: " << filename << "\n";
    fclose(ifre);
    return 0;
  }
  if (strcmp(tag, "FEEXP_SS_MESH_STATE_V2") != 0) {
    outw << "Unsupported mesh state version " << tag
         << "; regenerate with the corrected V2 mesh: " << filename << "\n";
    fclose(ifre);
    return 0;
  }
  if (fscanf(ifre, "%79s %d", key, &file_n) != 2 ||
      strcmp(key, "N") != 0 ||
      fscanf(ifre, "%79s %d %d", key, &n_nodes, &n_attr) != 3 ||
      strcmp(key, "COUNTS") != 0) {
    outw << "Bad mesh state header: " << filename << "\n";
    fclose(ifre);
    return 0;
  }

  if (file_n != n) {
    outw << "Mesh state n mismatch: file n=" << file_n
         << " current n=" << n << "\n";
    fclose(ifre);
    return 0;
  }

  if (n_nodes != p_mesh->number_of_nodes ||
      n_attr != p_mesh->number_of_node_attributes) {
    outw << "Mesh state size mismatch: file nodes=" << n_nodes
         << " attrs=" << n_attr
         << " current nodes=" << p_mesh->number_of_nodes
         << " attrs=" << p_mesh->number_of_node_attributes << "\n";
    fclose(ifre);
    return 0;
  }

  double file_d_crucible = 0.;
  double file_end_time = 0.;
  double file_GAf = 0.;
  double file_GAr = 0.;
  double file_Rc = 0.;
  double file_lc = 0.;
  double file_hSlab = 0.;
  double file_h0 = 0.;
  if (fscanf(ifre, "%79s %lf %lf %lf %lf %lf %lf %lf %lf",
             key, &file_d_crucible, &file_end_time, &file_GAf,
             &file_GAr, &file_Rc, &file_lc, &file_hSlab, &file_h0) != 9 ||
      strcmp(key, "GEOM") != 0) {
    outw << "Bad mesh state GEOM block: " << filename << "\n";
    fclose(ifre);
    return 0;
  }

  if (fabs(file_d_crucible - d_crucible) > 1.e-12 ||
      fabs(file_GAf - GAf) > 1.e-12 ||
      fabs(file_GAr - GAr) > 1.e-12) {
    outw << "Mesh state parameter mismatch in " << filename << "\n";
    fclose(ifre);
    return 0;
  }
  (void)file_end_time;
  (void)file_Rc;
  (void)file_lc;
  (void)file_hSlab;

  double file_Ttop = 0.;
  double file_Tout = 0.;
  double file_Tmp = 0.;
  double file_Tbottom = 0.;
  double file_hh = 0.;
  double file_Rad = 0.;
  double file_PeS = 0.;
  double file_cond_solid = 0.;
  double file_cond_liquid = 0.;
  double file_Vp = 0.;
  if (fscanf(ifre, "%79s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
             key, &file_Ttop, &file_Tout, &file_Tmp, &file_Tbottom,
             &file_hh, &file_Rad, &file_PeS, &file_cond_solid,
             &file_cond_liquid, &file_Vp) != 11 ||
      strcmp(key, "LOCAL") != 0) {
    outw << "Bad mesh state LOCAL block: " << filename << "\n";
    fclose(ifre);
    return 0;
  }
  // Four retained movement-control fields remain in the V2 state layout.
  // They are no longer authoritative and are read only to preserve the
  // layout of states written by this driver.
  double file_dx = 0.;
  double file_Vpull = 0.;
  double file_betaRough = 0.;
  double file_betaStep = 0.;
  double file_kineticA2DN = 0.;
  double file_kineticB2DN = 0.;
  double legacyMaxNodeMoveDx = 0.;
  double legacyKinMaxNodeMoveDx = 0.;
  double legacyLgMoveCoef = 0.;
  double legacyLgMaxNodeMoveDx = 0.;
  double file_rhog = 0.;
  double file_gam = 0.;
  if (fscanf(ifre,
             "%79s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
             key, &file_dx, &file_Vpull, &file_betaRough,
             &file_betaStep, &file_kineticA2DN,
             &file_kineticB2DN, &legacyMaxNodeMoveDx,
             &legacyKinMaxNodeMoveDx, &legacyLgMoveCoef,
             &legacyLgMaxNodeMoveDx, &file_rhog, &file_gam) != 13 ||
      strcmp(key, "MESH") != 0) {
    outw << "Bad mesh state MESH block: " << filename << "\n";
    fclose(ifre);
    return 0;
  }
  (void)legacyMaxNodeMoveDx;
  (void)legacyKinMaxNodeMoveDx;
  (void)legacyLgMoveCoef;
  (void)legacyLgMaxNodeMoveDx;

  double file_voronkovTm = 0.;
  double file_voronkovQ = 0.;
  double file_voronkovLambdaSG = 0.;
  double file_voronkovAtomicDensity = 0.;
  int file_voronkovTemperatureIsAbsolute = 0;
  double file_voronkovAlphaTransition = 0.;
  double file_voronkovSigmaSLFacet = 0.;
  double file_voronkovSigmaSLRough = 0.;
  double file_voronkovSigmaSLPrimeFacet = 0.;
  double file_voronkovSigmaSG = 0.;
  double file_voronkovSigmaMG = 0.;
  double file_voronkovThetaFacet = 0.;
  double file_voronkovThetaStepCorrection = 0.;
  double file_voronkovVelocityToCm = 0.;
  double file_voronkovKV = 0.;
  if (fscanf(ifre,
             "%79s %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
             key, &file_voronkovTm, &file_voronkovQ,
             &file_voronkovLambdaSG, &file_voronkovAtomicDensity,
             &file_voronkovTemperatureIsAbsolute,
             &file_voronkovAlphaTransition, &file_voronkovSigmaSLFacet,
             &file_voronkovSigmaSLRough, &file_voronkovSigmaSLPrimeFacet,
             &file_voronkovSigmaSG, &file_voronkovSigmaMG,
             &file_voronkovThetaFacet, &file_voronkovThetaStepCorrection,
             &file_voronkovVelocityToCm, &file_voronkovKV) != 16 ||
      strcmp(key, "VOR") != 0 ||
      fscanf(ifre, "%79s", key) != 1 ||
      strcmp(key, "NODES") != 0) {
    outw << "Bad mesh state VOR/NODES block: " << filename << "\n";
    fclose(ifre);
    return 0;
  }

  // Restart files carry a snapshot of the physical model for provenance.
  // param.txt and the current driver geometry remain authoritative: accept a
  // state only when every stored physical value matches, and never overwrite
  // the active parameters with values from the file.
  const double restart_tol = 1.e-12;
  const double expected_Rc = 0.25 * d_crucible;
  const double expected_lc = sqrt(p_mesh->gam / p_mesh->rhog);
  const double sin_phi = sin(PI * GAr / 180.);
  const double expected_hSlab =
      expected_lc * sqrt(2.0 * max(0.0, 1.0 - sin_phi));
  if (fabs(file_Rc - expected_Rc) >
          restart_tol * max(1.e-30, max(fabs(file_Rc), fabs(expected_Rc))) ||
      fabs(file_lc - expected_lc) >
          restart_tol * max(1.e-30, max(fabs(file_lc), fabs(expected_lc))) ||
      fabs(file_hSlab - expected_hSlab) >
          restart_tol * max(1.e-30, max(fabs(file_hSlab), fabs(expected_hSlab))) ||
      fabs(file_h0 - p_local_calc->h0) >
          restart_tol * max(1.e-30, max(fabs(file_h0), fabs(p_local_calc->h0))) ||
      fabs(file_Ttop - p_local_calc->Ttop) >
          restart_tol * max(1.e-30, max(fabs(file_Ttop), fabs(p_local_calc->Ttop))) ||
      fabs(file_Tout - p_local_calc->Tout) >
          restart_tol * max(1.e-30, max(fabs(file_Tout), fabs(p_local_calc->Tout))) ||
      fabs(file_Tmp - p_local_calc->Tmp) >
          restart_tol * max(1.e-30, max(fabs(file_Tmp), fabs(p_local_calc->Tmp))) ||
      fabs(file_Tbottom - p_local_calc->Tbottom) >
          restart_tol * max(1.e-30, max(fabs(file_Tbottom), fabs(p_local_calc->Tbottom))) ||
      fabs(file_hh - p_local_calc->hh) >
          restart_tol * max(1.e-30, max(fabs(file_hh), fabs(p_local_calc->hh))) ||
      fabs(file_Rad - p_local_calc->Rad) >
          restart_tol * max(1.e-30, max(fabs(file_Rad), fabs(p_local_calc->Rad))) ||
      fabs(file_PeS - p_local_calc->PeS) >
          restart_tol * max(1.e-30, max(fabs(file_PeS), fabs(p_local_calc->PeS))) ||
      fabs(file_cond_solid - p_local_calc->cond[Mesh::SOLID]) >
          restart_tol * max(1.e-30, max(fabs(file_cond_solid),
                                        fabs(p_local_calc->cond[Mesh::SOLID]))) ||
      fabs(file_cond_liquid - p_local_calc->cond[Mesh::LIQUID]) >
          restart_tol * max(1.e-30, max(fabs(file_cond_liquid),
                                        fabs(p_local_calc->cond[Mesh::LIQUID]))) ||
      fabs(file_Vp - p_local_calc->Vp) >
          restart_tol * max(1.e-30, max(fabs(file_Vp), fabs(p_local_calc->Vp))) ||
      fabs(file_dx - p_mesh->dx) >
          restart_tol * max(1.e-30, max(fabs(file_dx), fabs(p_mesh->dx))) ||
      fabs(file_Vpull - p_mesh->Vpull) >
          restart_tol * max(1.e-30, max(fabs(file_Vpull), fabs(p_mesh->Vpull))) ||
      fabs(file_betaRough - p_mesh->betaRough) >
          restart_tol * max(1.e-30, max(fabs(file_betaRough), fabs(p_mesh->betaRough))) ||
      fabs(file_betaStep - p_mesh->betaStep) >
          restart_tol * max(1.e-30, max(fabs(file_betaStep), fabs(p_mesh->betaStep))) ||
      fabs(file_kineticA2DN - p_mesh->kineticA2DN) >
          restart_tol * max(1.e-30, max(fabs(file_kineticA2DN), fabs(p_mesh->kineticA2DN))) ||
      fabs(file_kineticB2DN - p_mesh->kineticB2DN) >
          restart_tol * max(1.e-30, max(fabs(file_kineticB2DN), fabs(p_mesh->kineticB2DN))) ||
      fabs(file_rhog - p_mesh->rhog) >
          restart_tol * max(1.e-30, max(fabs(file_rhog), fabs(p_mesh->rhog))) ||
      fabs(file_gam - p_mesh->gam) >
          restart_tol * max(1.e-30, max(fabs(file_gam), fabs(p_mesh->gam))) ||
      fabs(file_voronkovTm - p_mesh->voronkovTm) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovTm), fabs(p_mesh->voronkovTm))) ||
      fabs(file_voronkovQ - p_mesh->voronkovQ) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovQ), fabs(p_mesh->voronkovQ))) ||
      fabs(file_voronkovLambdaSG - p_mesh->voronkovLambdaSG) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovLambdaSG), fabs(p_mesh->voronkovLambdaSG))) ||
      fabs(file_voronkovAtomicDensity - p_mesh->voronkovAtomicDensity) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovAtomicDensity), fabs(p_mesh->voronkovAtomicDensity))) ||
      file_voronkovTemperatureIsAbsolute !=
          p_mesh->voronkovTemperatureIsAbsolute ||
      fabs(file_voronkovAlphaTransition - p_mesh->voronkovAlphaTransition) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovAlphaTransition), fabs(p_mesh->voronkovAlphaTransition))) ||
      fabs(file_voronkovSigmaSLFacet - p_mesh->voronkovSigmaSLFacet) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovSigmaSLFacet), fabs(p_mesh->voronkovSigmaSLFacet))) ||
      fabs(file_voronkovSigmaSLRough - p_mesh->voronkovSigmaSLRough) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovSigmaSLRough), fabs(p_mesh->voronkovSigmaSLRough))) ||
      fabs(file_voronkovSigmaSLPrimeFacet - p_mesh->voronkovSigmaSLPrimeFacet) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovSigmaSLPrimeFacet), fabs(p_mesh->voronkovSigmaSLPrimeFacet))) ||
      fabs(file_voronkovSigmaSG - p_mesh->voronkovSigmaSG) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovSigmaSG), fabs(p_mesh->voronkovSigmaSG))) ||
      fabs(file_voronkovSigmaMG - p_mesh->voronkovSigmaMG) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovSigmaMG), fabs(p_mesh->voronkovSigmaMG))) ||
      fabs(file_voronkovThetaFacet - p_mesh->voronkovThetaFacet) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovThetaFacet), fabs(p_mesh->voronkovThetaFacet))) ||
      fabs(file_voronkovThetaStepCorrection -
           p_mesh->voronkovThetaStepCorrection) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovThetaStepCorrection),
                                        fabs(p_mesh->voronkovThetaStepCorrection))) ||
      fabs(file_voronkovVelocityToCm - p_mesh->voronkovVelocityToCm) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovVelocityToCm), fabs(p_mesh->voronkovVelocityToCm))) ||
      fabs(file_voronkovKV - p_mesh->voronkovKV) >
          restart_tol * max(1.e-30, max(fabs(file_voronkovKV), fabs(p_mesh->voronkovKV)))) {
    outw << "Mesh state physical-parameter mismatch; current param.txt and "
            "driver geometry are authoritative: " << filename << "\n";
    fclose(ifre);
    return 0;
  }

  for (int i(0); i < p_mesh->number_of_nodes; i++) {
    int node_id = -1;
    if (fscanf(ifre, "%d", &node_id) != 1 ||
        node_id < 0 || node_id >= p_mesh->number_of_nodes) {
      outw << "Bad mesh state node id in " << filename << "\n";
      fclose(ifre);
      return 0;
    }

    for (int co(0); co < 3; co++) {
      if (fscanf(ifre, "%lf", &p_mesh->node[node_id].coord[co]) != 1) {
        outw << "Bad mesh state coordinate in " << filename << "\n";
        fclose(ifre);
        return 0;
      }
    }

    for (int a(0); a < p_mesh->number_of_node_attributes; a++) {
      if (fscanf(ifre, "%lf", &p_mesh->node[node_id].attribute[a]) != 1) {
        outw << "Bad mesh state attribute in " << filename << "\n";
        fclose(ifre);
        return 0;
      }
    }

    p_mesh->node[node_id].T = p_mesh->node[node_id].attribute[TEMPERATURE];
    p_mesh->node[node_id].dTdt = 0.;
    p_mesh->node[node_id].dTdt0 = 0.;
    p_mesh->node[node_id].cT = 0.;
  }

  fclose(ifre);
  outw << "read mesh state " << filename << "\n";
  return 1;
}

static int readAxisymmetricInitialState(const char *filename, Mesh *p_mesh,
                                        LocalCalc *p_local_calc,
                                        double d_crucible, double Rc,
                                        ostream &outw,
                                        double *p_target_tpl_temperature,
                                        double *p_bottom_sensitivity) {
  ifstream in(filename);
  if (!in) {
    outw << "Failed to open axisymmetric IC: " << filename << "\n";
    return 0;
  }

  AxisymmetricInitialState state = {};
  string key;
  if (!(in >> key) ||
      (key != "FEEXP_AXISYMMETRIC_IC_V1" &&
       key != "FEEXP_AXISYMMETRIC_IC_V2")) {
    outw << "Bad axisymmetric IC version: " << filename << "\n";
    return 0;
  }
  state.version = (key == "FEEXP_AXISYMMETRIC_IC_V2") ? 2 : 1;
  if (!(in >> key) || key != "DIMS" ||
      !(in >> state.crystalRadialCells >> state.outerRadialCells >>
            state.liquidVerticalCells >> state.solidVerticalCells) ||
      !(in >> key) || key != "GEOM" ||
      !(in >> state.bottomZ >> state.undisturbedMeltZ >> state.tplZ >>
            state.topZ >> state.crystalRadius >> state.crucibleRadius) ||
      !(in >> key) || key != "THERMAL" ||
      !(in >> state.temperatureTop >> state.temperatureOut >>
            state.temperatureBottom)) {
    outw << "Bad axisymmetric IC header: " << filename << "\n";
    return 0;
  }
  if (state.version >= 2) {
    if (!(in >> key) || key != "CALIBRATION" ||
        !(in >> state.targetTplTemperature >> state.bottomSensitivity)) {
      outw << "Bad axisymmetric IC calibration: " << filename << "\n";
      return 0;
    }
  }
  if (!(in >> key) || key != "MODEL" ||
      !(in >> state.diffusivityLiquid >> state.diffusivitySolid >>
            state.heatTransferCoefficient >> state.radiationCoefficient >>
            state.latentHeatCoefficient >> state.pullVelocity >>
            state.sigmaMG >> state.rhoG >> state.meltingTemperature)) {
    outw << "Bad axisymmetric IC header: " << filename << "\n";
    return 0;
  }

  const int nc = state.crystalRadialCells;
  const int no = state.outerRadialCells;
  const int nl = state.liquidVerticalCells;
  const int ns = state.solidVerticalCells;
  if (nc < 2 || no < 2 || nl < 2 || ns < 2 ||
      nc > 10000 || no > 10000 || nl > 10000 || ns > 10000) {
    outw << "Invalid axisymmetric IC dimensions: " << filename << "\n";
    return 0;
  }

  state.liquidT.resize((nl + 1) * (nc + 1));
  state.solidT.resize((ns + 1) * (nc + 1));
  state.annulusT.resize((nl + 1) * (no + 1));

  if (!(in >> key) || key != "LIQUID") {
    outw << "Missing LIQUID block in " << filename << "\n";
    return 0;
  }
  for (size_t i = 0; i < state.liquidT.size(); i++)
    if (!(in >> state.liquidT[i])) {
      outw << "Bad LIQUID block in " << filename << "\n";
      return 0;
    }

  if (!(in >> key) || key != "SOLID") {
    outw << "Missing SOLID block in " << filename << "\n";
    return 0;
  }
  for (size_t i = 0; i < state.solidT.size(); i++)
    if (!(in >> state.solidT[i])) {
      outw << "Bad SOLID block in " << filename << "\n";
      return 0;
    }

  if (!(in >> key) || key != "ANNULUS") {
    outw << "Missing ANNULUS block in " << filename << "\n";
    return 0;
  }
  for (size_t i = 0; i < state.annulusT.size(); i++)
    if (!(in >> state.annulusT[i])) {
      outw << "Bad ANNULUS block in " << filename << "\n";
      return 0;
    }

  if (state.version >= 2) {
    int slPoints = 0;
    if (!(in >> key >> slPoints) || key != "SL_PROFILE" ||
        slPoints != nc + 1) {
      outw << "Missing SL_PROFILE block in " << filename << "\n";
      return 0;
    }
    state.slR.resize(slPoints);
    state.slZ.resize(slPoints);
    for (int i = 0; i < slPoints; i++)
      if (!(in >> state.slR[i] >> state.slZ[i])) {
        outw << "Bad SL_PROFILE block in " << filename << "\n";
        return 0;
      }
  } else {
    state.slR.resize(nc + 1);
    state.slZ.resize(nc + 1, state.tplZ);
    for (int i = 0; i <= nc; i++)
      state.slR[i] = state.crystalRadius * i / nc;
  }

  int meniscusPoints = 0;
  if (!(in >> key >> meniscusPoints) || key != "MENISCUS" ||
      meniscusPoints < no + 1 || meniscusPoints > 1000000) {
    outw << "Missing MENISCUS block in " << filename << "\n";
    return 0;
  }
  state.meniscusR.resize(meniscusPoints);
  state.meniscusZ.resize(meniscusPoints);
  state.meniscusArc.resize(meniscusPoints, 0.);
  for (int i = 0; i < meniscusPoints; i++)
    if (!(in >> state.meniscusR[i] >> state.meniscusZ[i])) {
      outw << "Bad MENISCUS block in " << filename << "\n";
      return 0;
    }
  if (!(in >> key) || key != "END") {
    outw << "Missing END marker in " << filename << "\n";
    return 0;
  }

  const double geometryTol = 1.e-10;
  if (fabs(state.bottomZ) > geometryTol ||
      fabs(state.topZ - d_crucible) > geometryTol ||
      fabs(state.crystalRadius - Rc) > geometryTol ||
      fabs(state.crucibleRadius - .5 * d_crucible) > geometryTol ||
      fabs(state.undisturbedMeltZ - p_mesh->ZI) > geometryTol ||
      fabs(state.temperatureTop - p_local_calc->Ttop) > 1.e-10 ||
      fabs(state.temperatureOut - p_local_calc->Tout) > 1.e-10 ||
      fabs(state.diffusivityLiquid -
           p_local_calc->cond[Mesh::LIQUID]) > 1.e-14 ||
      fabs(state.diffusivitySolid -
           p_local_calc->cond[Mesh::SOLID]) > 1.e-14 ||
      fabs(state.heatTransferCoefficient - p_local_calc->hh) > 1.e-14 ||
      fabs(state.radiationCoefficient - p_local_calc->Rad) > 1.e-18 ||
      fabs(state.latentHeatCoefficient - p_local_calc->PeS) > 1.e-10 ||
      fabs(state.pullVelocity - p_mesh->Vpull) > 1.e-14 ||
      fabs(state.sigmaMG - p_mesh->gam) > 1.e-12 ||
      fabs(state.rhoG - p_mesh->rhog) > 1.e-8 ||
      fabs(state.meltingTemperature - p_local_calc->Tmp) > 1.e-10 ||
      (state.version >= 2 &&
       (state.bottomSensitivity <= 0.0 ||
        fabs(state.targetTplTemperature +
             fabs(state.pullVelocity) / p_mesh->betaRough) > 1.e-6)) ||
      !(state.bottomZ < state.undisturbedMeltZ &&
        state.undisturbedMeltZ < state.tplZ && state.tplZ < state.topZ)) {
    outw << "Axisymmetric IC parameters do not match the 3D case: "
         << filename << "\n";
    return 0;
  }

  for (int i = 0; i <= nc; i++) {
    if (!isfinite(state.slR[i]) || !isfinite(state.slZ[i]) ||
        state.slZ[i] <= state.bottomZ || state.slZ[i] >= state.topZ ||
        (i > 0 && state.slR[i] <= state.slR[i - 1])) {
      outw << "Invalid axisymmetric S/L profile: " << filename << "\n";
      return 0;
    }
  }
  if (fabs(state.slR.front()) > geometryTol ||
      fabs(state.slR.back() - Rc) > geometryTol ||
      fabs(state.slZ.back() - state.tplZ) > geometryTol) {
    outw << "Axisymmetric S/L endpoints do not match the 3D case: "
         << filename << "\n";
    return 0;
  }

  for (int i = 0; i < meniscusPoints; i++) {
    if (!isfinite(state.meniscusR[i]) || !isfinite(state.meniscusZ[i]) ||
        (i > 0 && state.meniscusR[i] <= state.meniscusR[i - 1])) {
      outw << "Invalid axisymmetric meniscus profile: " << filename << "\n";
      return 0;
    }
    if (i > 0) {
      const double dr = state.meniscusR[i] - state.meniscusR[i - 1];
      const double dz = state.meniscusZ[i] - state.meniscusZ[i - 1];
      state.meniscusArc[i] =
          state.meniscusArc[i - 1] + sqrt(dr * dr + dz * dz);
    }
  }
  if (fabs(state.meniscusR.front() - Rc) > geometryTol ||
      fabs(state.meniscusR.back() - state.crucibleRadius) > geometryTol ||
      fabs(state.meniscusZ.front() - state.tplZ) > geometryTol ||
      fabs(state.meniscusZ.back() - state.undisturbedMeltZ) > geometryTol) {
    outw << "Axisymmetric meniscus endpoints do not match the 3D case: "
         << filename << "\n";
    return 0;
  }

  const auto meniscusZ = [&](double r) {
    if (r <= state.meniscusR.front()) return state.meniscusZ.front();
    if (r >= state.meniscusR.back()) return state.meniscusZ.back();
    const int i = (int)(upper_bound(state.meniscusR.begin(),
                                    state.meniscusR.end(), r) -
                        state.meniscusR.begin()) - 1;
    const double f = (r - state.meniscusR[i]) /
                     (state.meniscusR[i + 1] - state.meniscusR[i]);
    return state.meniscusZ[i] +
           f * (state.meniscusZ[i + 1] - state.meniscusZ[i]);
  };

  const auto solidLiquidZ = [&](double r) {
    if (r <= state.slR.front()) return state.slZ.front();
    if (r >= state.slR.back()) return state.slZ.back();
    const int i = (int)(upper_bound(state.slR.begin(), state.slR.end(), r) -
                        state.slR.begin()) - 1;
    const double f =
        (r - state.slR[i]) / (state.slR[i + 1] - state.slR[i]);
    return state.slZ[i] + f * (state.slZ[i + 1] - state.slZ[i]);
  };

  const double initialSLZ =
      p_mesh->bc_set[p_mesh->BC_SOLID_LIQUID].bc_at_node[0].p_node->coord[2];
  const double initialZI = p_mesh->ZI;
  for (int i = 0; i < p_mesh->number_of_nodes; i++) {
    Node *p_node = &p_mesh->node[i];
    const double r = sqrt(p_node->coord[0] * p_node->coord[0] +
                          p_node->coord[1] * p_node->coord[1]);
    const double z = p_node->coord[2];
    if (r <= Rc + geometryTol) {
      const double localSLZ = solidLiquidZ(min(r, Rc));
      if (z <= initialSLZ)
        p_node->coord[2] = state.bottomZ +
            (z - state.bottomZ) * (localSLZ - state.bottomZ) /
                (initialSLZ - state.bottomZ);
      else
        p_node->coord[2] = localSLZ +
            (z - initialSLZ) * (state.topZ - localSLZ) /
                (state.topZ - initialSLZ);
    } else {
      const double surfaceZ = meniscusZ(r);
      p_node->coord[2] = state.bottomZ +
          (z - state.bottomZ) * (surfaceZ - state.bottomZ) /
              (initialZI - state.bottomZ);
    }
  }

  // Place every generated L/G row at equal arc length on the imported
  // Young-Laplace profile. This is also the spacing sought by calcRelaxStr(),
  // so mesh redistribution does not immediately distort the initial profile.
  BCSet *p_lg = &p_mesh->bc_set[p_mesh->BC_LIQUID_GAS];
  for (int j = 0; j < p_lg->nS2; j++)
    for (int i = 0; i < p_lg->nS1; i++) {
      Node *p_node = p_lg->bc_at_nodeStr[i][j]->p_node;
      double r = state.meniscusR.front();
      if (j == p_lg->nS2 - 1) {
        r = state.meniscusR.back();
      } else if (j > 0) {
        const double targetArc =
            state.meniscusArc.back() * j / (p_lg->nS2 - 1);
        const int profile =
            (int)(upper_bound(state.meniscusArc.begin(),
                              state.meniscusArc.end(), targetArc) -
                  state.meniscusArc.begin()) - 1;
        const double f =
            (targetArc - state.meniscusArc[profile]) /
            (state.meniscusArc[profile + 1] -
             state.meniscusArc[profile]);
        r = state.meniscusR[profile] +
            f * (state.meniscusR[profile + 1] -
                 state.meniscusR[profile]);
      }
      const double oldRadius =
          sqrt(p_node->coord[0] * p_node->coord[0] +
               p_node->coord[1] * p_node->coord[1]);
      p_node->coord[0] *= r / oldRadius;
      p_node->coord[1] *= r / oldRadius;
      p_node->coord[2] = meniscusZ(r);
    }

  BCSet *p_sl = &p_mesh->bc_set[p_mesh->BC_SOLID_LIQUID];
  for (int i = 0; i < p_sl->number_of_recordsNode; i++) {
    Node *p_node = p_sl->bc_at_node[i].p_node;
    const double r = sqrt(p_node->coord[0] * p_node->coord[0] +
                          p_node->coord[1] * p_node->coord[1]);
    p_node->coord[2] = solidLiquidZ(min(r, Rc));
  }

  BCSet *p_tpl = &p_mesh->bc_set[p_mesh->BC_TPL];
  for (int i = 0; i < p_tpl->nS1; i++) {
    Node *p_node = p_tpl->bc_at_nodeStr[i][1]->p_node;
    const double r = sqrt(p_node->coord[0] * p_node->coord[0] +
                          p_node->coord[1] * p_node->coord[1]);
    p_node->coord[0] *= Rc / r;
    p_node->coord[1] *= Rc / r;
    p_node->coord[2] = state.tplZ;
  }

  const auto bilinear = [](const vector<double> &field, int nx,
                           int i, int j, double fr, double fz) {
    const double t00 = field[j * nx + i];
    const double t10 = field[j * nx + i + 1];
    const double t01 = field[(j + 1) * nx + i];
    const double t11 = field[(j + 1) * nx + i + 1];
    return (1. - fz) * ((1. - fr) * t00 + fr * t10) +
           fz * ((1. - fr) * t01 + fr * t11);
  };

  for (int n = 0; n < p_mesh->number_of_nodes; n++) {
    Node *p_node = &p_mesh->node[n];
    double r = sqrt(p_node->coord[0] * p_node->coord[0] +
                    p_node->coord[1] * p_node->coord[1]);
    double z = p_node->coord[2];
    if (r > state.crucibleRadius + geometryTol ||
        z < state.bottomZ - geometryTol || z > state.topZ + geometryTol) {
      outw << "3D node outside the axisymmetric IC domain: node=" << n
           << " r=" << r << " z=" << z << "\n";
      return 0;
    }

    double temperature = 0.;
    if (r <= Rc + geometryTol) {
      if (r > Rc) r = Rc;
      double ur = nc * r / Rc;
      int i = (int)floor(ur);
      if (i >= nc) i = nc - 1;
      const double fr = ur - i;
      const double localSLZ = solidLiquidZ(r);

      if (z <= localSLZ) {
        if (z < state.bottomZ) z = state.bottomZ;
        double uz = nl * (z - state.bottomZ) /
                    (localSLZ - state.bottomZ);
        int j = (int)floor(uz);
        if (j >= nl) j = nl - 1;
        temperature = bilinear(state.liquidT, nc + 1, i, j,
                               fr, uz - j);
      } else {
        double uz = ns * (z - localSLZ) /
                    (state.topZ - localSLZ);
        int j = (int)floor(uz);
        if (j >= ns) j = ns - 1;
        temperature = bilinear(state.solidT, nc + 1, i, j,
                               fr, uz - j);
      }
    } else {
      if (r > state.crucibleRadius) r = state.crucibleRadius;
      double ur = no * (r - Rc) / (state.crucibleRadius - Rc);
      int i = (int)floor(ur);
      if (i >= no) i = no - 1;
      const double fr = ur - i;
      const double surfaceZ = meniscusZ(r);
      if (z > surfaceZ + geometryTol) {
        outw << "3D liquid node above imported meniscus: node=" << n
             << " r=" << r << " z=" << z
             << " surface=" << surfaceZ << "\n";
        return 0;
      }
      if (z > surfaceZ) z = surfaceZ;
      double uz = nl * (z - state.bottomZ) /
                  (surfaceZ - state.bottomZ);
      int j = (int)floor(uz);
      if (j >= nl) j = nl - 1;
      temperature = bilinear(state.annulusT, no + 1, i, j,
                             fr, uz - j);
    }

    if (!isfinite(temperature)) {
      outw << "Non-finite imported temperature at 3D node " << n << "\n";
      return 0;
    }
    p_node->attribute[TEMPERATURE] = temperature;
    p_node->T = temperature;
    p_node->dTdt = 0.;
    p_node->dTdt0 = 0.;
    p_node->cT = 0.;
  }

  p_mesh->ZI = state.undisturbedMeltZ;
  p_local_calc->h0 = state.tplZ - state.undisturbedMeltZ;
  p_local_calc->Tbottom = state.temperatureBottom;
  *p_target_tpl_temperature = state.targetTplTemperature;
  *p_bottom_sensitivity = state.bottomSensitivity;
  outw << "axisymmetric IC=" << filename
       << " version=" << state.version
       << " h_YL_mm=" << 1000. * p_local_calc->h0
       << " Tbottom_fit=" << p_local_calc->Tbottom
       << " dZSL_um="
       << 1.e6 * (*max_element(state.slZ.begin(), state.slZ.end()) -
                  *min_element(state.slZ.begin(), state.slZ.end()))
       << " grids=" << nc << "x" << nl << "/"
       << nc << "x" << ns << "/" << no << "x" << nl << "\n";
  return state.version;
}

static int initialize3dInitialState(
    const char *ssStateInputFile, const char *axisymmetricInitialStateFile,
    Mesh *p_mesh, LocalCalc *p_local_calc, int n, double d_crucible,
    double Rc, double GAf, double GAr, ostream &outw,
    int *p_restartFromStoredState,
    int *p_initializedFromAxisymmetricState,
    double *p_axisymmetricTplTarget,
    double *p_axisymmetricBottomSensitivity) {
  *p_restartFromStoredState = 0;
  *p_initializedFromAxisymmetricState = 0;
  *p_axisymmetricTplTarget = 0.;
  *p_axisymmetricBottomSensitivity = 0.;

  struct stat stateInfo;
  if (stat(ssStateInputFile, &stateInfo) == 0) {
    if (!readMeshState(ssStateInputFile, p_mesh, p_local_calc,
                       n, d_crucible, GAf, GAr, outw))
      return 0;

    *p_restartFromStoredState = 1;
    // Restart snapshots store geometry and fields but not derived S/L
    // velocities. Reconstruct those diagnostics without moving the mesh.
    p_mesh->bc_set[p_mesh->BC_SOLID_LIQUID]
        .kinMoveInterface(0., p_mesh);
    return 1;
  }

  if (stat(axisymmetricInitialStateFile, &stateInfo) == 0) {
    *p_initializedFromAxisymmetricState =
        readAxisymmetricInitialState(
            axisymmetricInitialStateFile, p_mesh, p_local_calc,
            d_crucible, Rc, outw, p_axisymmetricTplTarget,
            p_axisymmetricBottomSensitivity);
    return (*p_initializedFromAxisymmetricState != 0);
  }

  outw << "axisymmetric IC not found: " << axisymmetricInitialStateFile
       << "; using linear axial temperature and engineering meniscus lift\n";
  return 1;
}

static void initFEStepStorage(FEStepContext *p_fe_step, Mesh *p_mesh,
                              ostream &outw) {
  LocalCalc &local_calc = p_fe_step->local_calc;

  p_fe_step->gp_set[QUADRATIC_4_NOD_LAGR].gp_flag = GP25_2D_QUAD;
  p_fe_step->gp_set[QUADRATIC_4_NOD_LAGR].bf_flag = QUADRATIC_4_NOD_LAGR;
  p_fe_step->gp_set[QUADRATIC_4_NOD_LAGR].BasisFuncAtGaussPointsInit();

  p_fe_step->gp_set[BRICK_8_NOD_LAGR].gp_flag = GP125_3D_BRICK;
  p_fe_step->gp_set[BRICK_8_NOD_LAGR].bf_flag = BRICK_8_NOD_LAGR;
  p_fe_step->gp_set[BRICK_8_NOD_LAGR].BasisFuncAtGaussPointsInit();

  outw << "gp brick=" << p_fe_step->gp_set[BRICK_8_NOD_LAGR].num_of_bf
       << "/" << p_fe_step->gp_set[BRICK_8_NOD_LAGR].num_of_gp
       << " quad=" << p_fe_step->gp_set[QUADRATIC_4_NOD_LAGR].num_of_bf
       << "/" << p_fe_step->gp_set[QUADRATIC_4_NOD_LAGR].num_of_gp << "\n";

  p_fe_step->gp_set_corn.gp_flag = GP8_3D_CORNERS;
  p_fe_step->gp_set_corn.bf_flag = BRICK_8_NOD_LAGR;
  p_fe_step->gp_set_corn.BasisFuncAtGaussPointsInit();

  int num_threads = 1;
#ifdef _OPENMP
  num_threads = omp_get_max_threads();
#endif

  p_fe_step->thread_calcs.resize(num_threads);
  for (int i(0); i < num_threads; i++) {
    memset(&p_fe_step->thread_calcs[i], 0, sizeof(p_fe_step->thread_calcs[i]));
    p_fe_step->thread_calcs[i].assembly = (int *)malloc(MAX_RESIDUAL * sizeof(int));
    p_fe_step->thread_calcs[i].local_res =
        (double *)malloc(MAX_RESIDUAL * sizeof(double));
    p_fe_step->thread_calcs[i].local_res_V1 =
        (double *)malloc(MAX_RESIDUAL * sizeof(double));
    p_fe_step->thread_calcs[i].local_res_V2 =
        (double *)malloc(MAX_RESIDUAL * sizeof(double));
    p_fe_step->thread_calcs[i].local_res_V3 =
        (double *)malloc(MAX_RESIDUAL * sizeof(double));
    p_fe_step->thread_calcs[i].local_res_P =
        (double *)malloc(MAX_RESIDUAL * sizeof(double));

    p_fe_step->thread_calcs[i].cond[Mesh::SOLID] =
        local_calc.cond[Mesh::SOLID];
    p_fe_step->thread_calcs[i].cond[Mesh::LIQUID] =
        local_calc.cond[Mesh::LIQUID];
    p_fe_step->thread_calcs[i].Rad = local_calc.Rad;
    p_fe_step->thread_calcs[i].hh = local_calc.hh;
    p_fe_step->thread_calcs[i].Tout = local_calc.Tout;
    p_fe_step->thread_calcs[i].Ttop = local_calc.Ttop;
    p_fe_step->thread_calcs[i].Tbottom = local_calc.Tbottom;
    p_fe_step->thread_calcs[i].Tmp = local_calc.Tmp;
    p_fe_step->thread_calcs[i].h0 = local_calc.h0;
    p_fe_step->thread_calcs[i].PeS = local_calc.PeS;
    p_fe_step->thread_calcs[i].Vp = local_calc.Vp;
  }

  p_fe_step->dTdt_tls.resize(num_threads * p_mesh->number_of_nodes, 0.);
  p_fe_step->cT_tls.resize(num_threads * p_mesh->number_of_nodes, 0.);
}

static int ensureDirectoryExists(const char *path, ostream &outw) {
#ifdef _WIN32
  if (_mkdir(path) == 0)
    return 1;
#else
  if (mkdir(path, 0777) == 0)
    return 1;
#endif

  struct stat info;
  if (stat(path, &info) == 0 && S_ISDIR(info.st_mode))
    return 1;

  outw << "Failed to create output directory: " << path << "\n";
  return 0;
}

static double runFEtemperatureStep(Mesh *p_mesh, FEStepContext *p_fe_step) {
  LocalCalc *p_local_calc = &p_fe_step->local_calc;
  // Quasi-steady latent heat uses the imposed pulling velocity, not the
  // instantaneous local S/L kinetic velocity.
  p_local_calc->Vp = p_mesh->Vpull;
  const double step_dt = p_fe_step->dt;
  const bool use_openmp = p_fe_step->use_openmp;
  const int num_nodes = p_mesh->number_of_nodes;
  const int cn_iters = 3;
  const double temp_min = (p_local_calc->Ttop < p_local_calc->Tout)
                              ? p_local_calc->Ttop
                              : p_local_calc->Tout;
  const double temp_max = p_local_calc->Tbottom;
  const int num_threads = (int)p_fe_step->thread_calcs.size();
  const int tls_size = num_threads * num_nodes;
  GaussPointBFSet *p_gp_brick = &p_fe_step->gp_set[BRICK_8_NOD_LAGR];
  GaussPointBFSet *p_gp_quad = &p_fe_step->gp_set[QUADRATIC_4_NOD_LAGR];
  vector<LocalCalc> &thread_calcs = p_fe_step->thread_calcs;
  vector<double> &dTdt_tls = p_fe_step->dTdt_tls;
  vector<double> &cT_tls = p_fe_step->cT_tls;

  for (int i(0); i < p_mesh->bc_set[p_mesh->BC_HOT].number_of_recordsNode; i++)
    p_mesh->bc_set[p_mesh->BC_HOT]
        .bc_at_node[i]
        .p_node->attribute[TEMPERATURE] = p_local_calc->Tbottom;
  for (int i(0); i < p_mesh->bc_set[p_mesh->BC_COLD].number_of_recordsNode; i++)
    p_mesh->bc_set[p_mesh->BC_COLD]
        .bc_at_node[i]
        .p_node->attribute[TEMPERATURE] = p_local_calc->Ttop;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_openmp)
#endif
  for (int i = 0; i < num_nodes; i++)
    p_mesh->node[i].T = p_mesh->node[i].attribute[TEMPERATURE];

  for (int it(0); it < num_threads; it++) {
    thread_calcs[it].cond[Mesh::SOLID] = p_local_calc->cond[Mesh::SOLID];
    thread_calcs[it].cond[Mesh::LIQUID] = p_local_calc->cond[Mesh::LIQUID];
    thread_calcs[it].Rad = p_local_calc->Rad;
    thread_calcs[it].hh = p_local_calc->hh;
    thread_calcs[it].Tout = p_local_calc->Tout;
    thread_calcs[it].Ttop = p_local_calc->Ttop;
    thread_calcs[it].Tbottom = p_local_calc->Tbottom;
    thread_calcs[it].Tmp = p_local_calc->Tmp;
    thread_calcs[it].h0 = p_local_calc->h0;
    thread_calcs[it].PeS = p_local_calc->PeS;
    thread_calcs[it].Vp = p_mesh->Vpull;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_openmp)
#endif
  for (int i = 0; i < tls_size; i++) {
    dTdt_tls[i] = 0.;
    cT_tls[i] = 0.;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_openmp)
#endif
  for (int elem_num = 0; elem_num < p_mesh->number_of_elements; elem_num++) {
    int tid(0);
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    LocalCalc *p_calc = &thread_calcs[tid];
    double *p_dTdt_local = &dTdt_tls[tid * num_nodes];
    double *p_cT_local = &cT_tls[tid * num_nodes];

    p_calc->LocalGetResJac3dDin(&p_mesh->element[elem_num], p_gp_brick);
    for (int i(0); i < 8; i++) {
      const int node_num = p_mesh->element[elem_num].p_node[i]->node_num;
      p_dTdt_local[node_num] += p_calc->local_res[i];
      p_cT_local[node_num] += p_calc->local_res_V1[i];
    }
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_openmp)
#endif
  for (int i = 0; i < num_nodes; i++) {
    p_mesh->node[i].dTdt = 0.;
    p_mesh->node[i].cT = 0.;
    for (int tid(0); tid < num_threads; tid++) {
      p_mesh->node[i].dTdt += dTdt_tls[tid * num_nodes + i];
      p_mesh->node[i].cT += cT_tls[tid * num_nodes + i];
    }
  }

  for (int ib(0); ib < p_mesh->bc_set[p_mesh->BC_LIQUID_GAS].number_of_recordsFace;
       ib++) {
    p_local_calc->LocalGetResJac2dOutwall(
        &p_mesh->bc_set[p_mesh->BC_LIQUID_GAS].bc_at_face[ib], p_gp_quad);
    for (int i(0); i < 4; i++)
      p_mesh->bc_set[p_mesh->BC_LIQUID_GAS]
          .bc_at_face[ib]
          .p_node[i]
          ->dTdt += p_local_calc->local_res[i];
  }
  for (int ib(0); ib < p_mesh->bc_set[p_mesh->BC_SOLID_GAS].number_of_recordsFace;
       ib++) {
    p_local_calc->LocalGetResJac2dOutwall(
        &p_mesh->bc_set[p_mesh->BC_SOLID_GAS].bc_at_face[ib], p_gp_quad);
    for (int i(0); i < 4; i++)
      p_mesh->bc_set[p_mesh->BC_SOLID_GAS]
          .bc_at_face[ib]
          .p_node[i]
          ->dTdt += p_local_calc->local_res[i];
  }
  for (int ib(0);
       ib < p_mesh->bc_set[p_mesh->BC_SOLID_LIQUID].number_of_recordsFace;
       ib++) {
    p_local_calc->LocalGetResidual2dInterfEntalp(
        &p_mesh->bc_set[p_mesh->BC_SOLID_LIQUID].bc_at_face[ib], p_gp_quad);
    for (int i(0); i < 4; i++)
      p_mesh->bc_set[p_mesh->BC_SOLID_LIQUID]
          .bc_at_face[ib]
          .p_node[i]
          ->dTdt += p_local_calc->local_res[i];
  }

  for (int i(0); i < p_mesh->bc_set[p_mesh->BC_COLD].number_of_recordsNode; i++)
    p_mesh->bc_set[p_mesh->BC_COLD].bc_at_node[i].p_node->dTdt = 0.;
  for (int i(0); i < p_mesh->bc_set[p_mesh->BC_HOT].number_of_recordsNode; i++)
    p_mesh->bc_set[p_mesh->BC_HOT].bc_at_node[i].p_node->dTdt = 0.;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_openmp)
#endif
  for (int i = 0; i < num_nodes; i++) {
    Node *p_node = &p_mesh->node[i];
    p_node->dTdt0 = p_node->dTdt;
    p_node->IntBf = p_node->cT;
    const double denom = p_node->IntBf + 1.e-12 * step_dt;
    double temperature;
    if (!isfinite(p_node->dTdt0) || !isfinite(p_node->IntBf) ||
        !isfinite(denom) || fabs(denom) <= 1.e-18)
      temperature = p_node->T;
    else
      temperature = p_node->T + step_dt * p_node->dTdt0 / denom;
    if (temperature < temp_min)
      temperature = temp_min;
    if (temperature > temp_max)
      temperature = temp_max;
    p_node->attribute[TEMPERATURE] = temperature;
  }

  for (int iter(0); iter < cn_iters; iter++) {
    for (int i(0); i < p_mesh->bc_set[p_mesh->BC_HOT].number_of_recordsNode;
         i++)
      p_mesh->bc_set[p_mesh->BC_HOT]
          .bc_at_node[i]
          .p_node->attribute[TEMPERATURE] = p_local_calc->Tbottom;
    for (int i(0); i < p_mesh->bc_set[p_mesh->BC_COLD].number_of_recordsNode;
         i++)
      p_mesh->bc_set[p_mesh->BC_COLD]
          .bc_at_node[i]
          .p_node->attribute[TEMPERATURE] = p_local_calc->Ttop;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_openmp)
#endif
    for (int i = 0; i < tls_size; i++) {
      dTdt_tls[i] = 0.;
      cT_tls[i] = 0.;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_openmp)
#endif
    for (int elem_num = 0; elem_num < p_mesh->number_of_elements; elem_num++) {
      int tid(0);
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      LocalCalc *p_calc = &thread_calcs[tid];
      double *p_dTdt_local = &dTdt_tls[tid * num_nodes];
      double *p_cT_local = &cT_tls[tid * num_nodes];

      p_calc->LocalGetResJac3dDin(&p_mesh->element[elem_num], p_gp_brick);
      for (int i(0); i < 8; i++) {
        const int node_num = p_mesh->element[elem_num].p_node[i]->node_num;
        p_dTdt_local[node_num] += p_calc->local_res[i];
        p_cT_local[node_num] += p_calc->local_res_V1[i];
      }
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_openmp)
#endif
    for (int i = 0; i < num_nodes; i++) {
      p_mesh->node[i].dTdt = 0.;
      p_mesh->node[i].cT = 0.;
      for (int tid(0); tid < num_threads; tid++) {
        p_mesh->node[i].dTdt += dTdt_tls[tid * num_nodes + i];
        p_mesh->node[i].cT += cT_tls[tid * num_nodes + i];
      }
    }

    for (int ib(0);
         ib < p_mesh->bc_set[p_mesh->BC_LIQUID_GAS].number_of_recordsFace;
         ib++) {
      p_local_calc->LocalGetResJac2dOutwall(
          &p_mesh->bc_set[p_mesh->BC_LIQUID_GAS].bc_at_face[ib], p_gp_quad);
      for (int i(0); i < 4; i++)
        p_mesh->bc_set[p_mesh->BC_LIQUID_GAS]
            .bc_at_face[ib]
            .p_node[i]
            ->dTdt += p_local_calc->local_res[i];
    }
    for (int ib(0);
         ib < p_mesh->bc_set[p_mesh->BC_SOLID_GAS].number_of_recordsFace;
         ib++) {
      p_local_calc->LocalGetResJac2dOutwall(
          &p_mesh->bc_set[p_mesh->BC_SOLID_GAS].bc_at_face[ib], p_gp_quad);
      for (int i(0); i < 4; i++)
        p_mesh->bc_set[p_mesh->BC_SOLID_GAS]
            .bc_at_face[ib]
            .p_node[i]
            ->dTdt += p_local_calc->local_res[i];
    }
    for (int ib(0);
         ib < p_mesh->bc_set[p_mesh->BC_SOLID_LIQUID].number_of_recordsFace;
         ib++) {
      p_local_calc->LocalGetResidual2dInterfEntalp(
          &p_mesh->bc_set[p_mesh->BC_SOLID_LIQUID].bc_at_face[ib], p_gp_quad);
      for (int i(0); i < 4; i++)
        p_mesh->bc_set[p_mesh->BC_SOLID_LIQUID]
            .bc_at_face[ib]
            .p_node[i]
            ->dTdt += p_local_calc->local_res[i];
    }

    for (int i(0); i < p_mesh->bc_set[p_mesh->BC_COLD].number_of_recordsNode;
         i++)
      p_mesh->bc_set[p_mesh->BC_COLD].bc_at_node[i].p_node->dTdt = 0.;
    for (int i(0); i < p_mesh->bc_set[p_mesh->BC_HOT].number_of_recordsNode;
         i++)
      p_mesh->bc_set[p_mesh->BC_HOT].bc_at_node[i].p_node->dTdt = 0.;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_openmp)
#endif
    for (int i = 0; i < num_nodes; i++) {
      Node *p_node = &p_mesh->node[i];
      const double residual_new = p_node->dTdt;
      const double capacity_new = p_node->cT;
      const double capacity_mid = 0.5 * (p_node->IntBf + capacity_new);
      const double denom = capacity_mid + 1.e-12 * step_dt;
      double temperature;
      if (!isfinite(p_node->dTdt0) || !isfinite(residual_new) ||
          !isfinite(p_node->IntBf) || !isfinite(capacity_new) ||
          !isfinite(capacity_mid) || !isfinite(denom) || fabs(denom) <= 1.e-18)
        temperature = p_node->T;
      else
        temperature = p_node->T +
                      0.5 * step_dt * (p_node->dTdt0 + residual_new) / denom;
      if (temperature < temp_min)
        temperature = temp_min;
      if (temperature > temp_max)
        temperature = temp_max;
      p_node->attribute[TEMPERATURE] = temperature;
    }
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_openmp)
#endif
  for (int i = 0; i < num_nodes; i++) {
    double temperature = p_mesh->node[i].attribute[TEMPERATURE];
    if (!isfinite(temperature))
      temperature = p_mesh->node[i].T;
    if (temperature < temp_min)
      temperature = temp_min;
    if (temperature > temp_max)
      temperature = temp_max;
    p_mesh->node[i].attribute[TEMPERATURE] = temperature;
  }
  for (int i(0); i < p_mesh->bc_set[p_mesh->BC_HOT].number_of_recordsNode; i++)
    p_mesh->bc_set[p_mesh->BC_HOT]
        .bc_at_node[i]
        .p_node->attribute[TEMPERATURE] = p_local_calc->Tbottom;
  for (int i(0); i < p_mesh->bc_set[p_mesh->BC_COLD].number_of_recordsNode; i++)
    p_mesh->bc_set[p_mesh->BC_COLD]
        .bc_at_node[i]
        .p_node->attribute[TEMPERATURE] = p_local_calc->Ttop;

  double norm(0.);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : norm) if (use_openmp)
#endif
  for (int i = 0; i < num_nodes; i++) {
    const double dTstep =
        p_mesh->node[i].attribute[TEMPERATURE] - p_mesh->node[i].T;
    norm += dTstep * dTstep;
    p_mesh->node[i].dTdt0 = p_mesh->node[i].dTdt;
    p_mesh->node[i].dTdt = 0.;
    p_mesh->node[i].cT = 0.;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_openmp)
#endif
  for (int i = 0; i < num_nodes; i++) {
    p_mesh->node[i].num_of_adj_elem_l = 0;
    p_mesh->node[i].num_of_adj_elem_s = 0;
    for (int iglob = 0; iglob < 3; iglob++)
      p_mesh->node[i].dTl[iglob] = 0.;
    for (int iglob = 0; iglob < 3; iglob++)
      p_mesh->node[i].dTs[iglob] = 0.;
  }
  for (int elem_num(0); elem_num < p_mesh->number_of_elements; elem_num++)
    p_local_calc->LocalGetGrad3d(&p_mesh->element[elem_num],
                                 &p_fe_step->gp_set_corn);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_openmp)
#endif
  for (int i = 0; i < p_mesh->number_of_nodes; i++) {
    if (p_mesh->node[i].num_of_adj_elem_l != 0)
      for (int iglob = 0; iglob < 3; iglob++)
        p_mesh->node[i].dTl[iglob] /= p_mesh->node[i].num_of_adj_elem_l;
    if (p_mesh->node[i].num_of_adj_elem_s != 0)
      for (int iglob = 0; iglob < 3; iglob++)
        p_mesh->node[i].dTs[iglob] /= p_mesh->node[i].num_of_adj_elem_s;
    for (int iglob = 0; iglob < 3; iglob++) {
      if (!isfinite(p_mesh->node[i].dTl[iglob]))
        p_mesh->node[i].dTl[iglob] = 0.;
      if (!isfinite(p_mesh->node[i].dTs[iglob]))
        p_mesh->node[i].dTs[iglob] = 0.;
    }
  }

  return norm;
}


#include "src/mesh/3d_driver_common.cpp"
#include "src/mesh/tpl_submesh.cpp"
