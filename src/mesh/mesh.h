/*****************************************************************************/
/*  (mesh.h)                                                                 */
/*****************************************************************************/
#include <cstdint>
#include <unordered_map>
#include <vector>

struct mesh;

#define NUM_OF_BF_SETS 20
#define NUM_OF_NODE_MARKERS 2
#define NUM_OF_NODE_ATTRIBUTES 8
#define MAX_RESIDUAL 20
#define PI 3.1415926536

const int n_columns(10);   // x y z xx yy xy xz yz 1.
const int n_columns7(7);   // x y z xx yy 1.
const int n_rows_max(280); //
int indxc[n_columns], indxr[n_columns], ipiv[n_columns];

double **E;
double **H;
double **Y;
double **A;
double **At;
double *T;
double *K;

typedef enum {
  EDGE,
  QUADRATERIAL,
  TRIANGLE,
  BRICK,
  PRIZM,
  TETRAGEDRON,
  PYRAMID
} ElementGeomTypes;

typedef enum {
  EDGE_3_NOD_LAGR_TRI_BOUNDARY_0,
  EDGE_3_NOD_LAGR_TRI_BOUNDARY_1,
  EDGE_3_NOD_LAGR_TRI_BOUNDARY_2,
  EDGE_2_NOD_LAGR,
  EDGE_3_NOD_LAGR,
  EDGE_4_NOD_LAGR,
  EDGE_4_NOD_HERM,
  BRICK_8_NOD_LAGR,
  QUADRATIC_9_NOD_LAGR,
  QUADRATIC_4_NOD_LAGR,
  TRIANGLE_3_NOD_LAGR,
  TRIANGLE_6_NOD_LAGR,
  TRIANGLE_9_NOD_LAGR,
  TETRAHEDRON_10_NOD_LAGR,
  TETRAHEDRON_4_NOD_LAGR
} BF_SetTypes;

typedef enum {
  GP8_1D_TRI_BOUNDARY_0,
  GP8_1D_TRI_BOUNDARY_1,
  GP8_1D_TRI_BOUNDARY_2,
  GP8_1D,
  GP25_2D_QUAD,
  GP12_2D_TRIANGLE,
  GP4_2D_TRIANGLE,
  GP125_3D_TETRAHEDRON_GAUSS_RADAU,
  GP29_3D_TETRAHEDRON,
  GP125_3D_BRICK,
  GP8_3D_CORNERS
} GP_SetTypes;

typedef enum { BC_NODE, BC_FACE } BC_SetTypes;

typedef enum {
  TEMPERATURE,
  CONCENTR,
  VELOCITY1,
  VELOCITY2,
  VELOCITY3,
  PRESSURE
} NodeAttributes;

// typedef enum {
//   INSIDE,
//   TOP,
//   BOTTOM,
//   SIDE,
//   INTERFACE
//   } BoundaryMarker;
//
// typedef enum {
//   SOLID,
//   MELT
//   } DomainMarker;
//
typedef struct coordinate {
  double coord[3];
} Coordinate;

void Error(char *s) {
  printf("%s\n", s);
  exit(99);
}

typedef struct node {
  double coord[3];
  double delta[3];
  double dZl[2];
  double dZs[2];
  //  double     Kt[n_rows_max];
  // int        ijk[n_rows_max][3];
  // double     H[n_columns][280];
  int bound_mark;
  int bulk_mark;
  int n_rows;
  double *attribute;
  double K;
  double Ka;
  double Kref;
  double Vz;
  double Vr;
  double Vn;

  double theta, theta0;
  double Ns[3];
  double Nl[3];
  double N[3];
  double dTl[3];
  double dTs[3];

  double dTdt, dTdt0, T, cT, zz0;
  double IntBf;
  double Q;
  double Z0;
  int InodeNum;
  // TPL submesh coupling. parent_nat uses the brick coordinates
  // (eta,nu,ksi) of BasisFuncBF3d8nodLagr(). coarse_node_id is retained as a
  // compatibility name; it identifies an exact node of the immediate parent,
  // which can itself be another submesh.
  int parent_element_id;
  double parent_nat[3];
  int coarse_node_id;
  int local_i;
  int local_j;
  int local_k;
  int local_block;
  int artificial_boundary;
  int parent_mapping_initialized;
  int n_facet;
  int *marker;
  int node_num;
  int tmp_local_node_num;
  int num_of_adj_elem;
  int num_of_adj_tri;
  int num_of_adj_elem_l;
  int num_of_adj_elem_s;
  int num_of_adj_nodes;
  //  node      *p_adj_node[8];
  int n_adj_node[8];
  double l_adj_node[8];
  int n_adj_tri[8];
  int LBMbcNum;
  int LBMrectNum;
} Node;

typedef struct element {
  int number_of_nodes;
  int element_num;
  int group_num;
  int element_type_flag;
  BF_SetTypes bf_set_flag;
  int *assembly;
  Node **p_node;
} Element;

typedef struct Triangle {
  Node *p_node[3];
  int group_num;
  double n[3];
  double centerCoord[3];
  double cx;
  double cy;
} Tri;

typedef struct {
  int assembly[2];
  Node *p_node[2];
  int n_adj_elem;
  int n_adj_elem_ids[2];
  double length;
} Edge;

typedef struct {
  Node *p_node;
  Node *p_adj_sg_node;
  Node *p_adj_lg_node;
  double value;
  double K;
  int node_num;
  double node_sg[3];
  double t_sl[3];
  double t_lg[3];
  double t_sg[3];
  double theta_sg;
  double N[3];
  double Zf, Th;
  double Vn, GA, GAs, GAv, SG, LG;
  double V_SL_n;
  double V_SG_n_local;
  double V_SG_t_far;
  double V_SG_t_local;
  double V_vo;
  // Voronkov force-balance threshold in cm/s, matching the theory closure.
  double Vmin_SG_far;
  double V_pull;
  double trace_dot_t_SG_far;
  // Signed local thermal driving force for TPL diagnostics:
  // DeltaT = dTkin = Tm - Tlocal. Positive is solidification, negative melting.
  double DeltaT;
  double chi0_abs;
  double chi_geom;
  double chi_use;
  double theta_g_0;
  double theta_g_app;
  double phiS;
  int phiS_mode;
  double alpha_TPL;
  double beta_step;
  double dTloc;
  double dTsolid;
  double sigmaSG_diag;
  double sigmaSL_diag;
  double sigmaLG_diag;
  double sigmaSL_f;
  double sigmaSL_alpha;
  double SigmaSL_CH;
  double psiSL_CH;
  double deltaS_CH;
  double theta_SG_far;
  double theta_SL_TPL;
  double theta_SG_stop;
  double chiV;
  double chiMax_FB;
  double theta_SG_TPL;
  // Voronkov correction status: 0 no correction, 1 force-balance stop,
  // 2 velocity-controlled Voronkov, 3 invalid geometry. Not a kinetic branch.
  int GABranch;
  int parent_face_id;
  int parent_bc_node_id;
  int coarse_bc_node_id;
  double parent_xi;
  double parent_eta;
  double parent_zeta;
  double interp_T_image;
  double interp_T_corr;
  double refine_coef;
  double relax_coef;
  double dcoord[3];
  double relax_sum[3];
  double move_abs;
  double curvatureRate; // [1/m^2], absolute cotangent row sum
  int num_of_adj_nodes;
  //  node      *p_adj_node[8];
  int n_adj_node[16];
  int num_of_adj_elem;
  int n_adj_elem[6];
  int i_adj_n;
  //  \  |
  //   \t|
  //    \|
} BCAtNode;

typedef struct {
  int i;
  int j;
  int k;
  int nBC;
  float Asurfel;
  float Vel;
} StrElem;

typedef struct {
  Element *p_elem;
  float n[3];
  BF_SetTypes bf_set_flag;
  int face_number;
  int *assembly;
  int *local_assembly;
  int number_of_nodes;
  Node **p_node;
  StrElem *strElem;
  Node *p_internal_node;
  int nStrElem;
} BCAtFace;

typedef struct {
  int assembly[2];
  Node *p_node[2];
  int n_adj_face;
  int n_adj_face_ids[2];
  double length;
} BCEdge;

typedef struct {
  BC_SetTypes bc_flag;
  int mesh_bc_id;
  int number_of_recordsNode;
  int number_of_recordsFace;
  int number_of_recordsEdge;
  BCAtNode *bc_at_node;
  BCAtFace *bc_at_face;
  BCEdge *bc_at_edge;
  BCAtNode ***bc_at_nodeStr;
  int nS1, nS2;
  int structured_wrap_columns;
  int structured_z_ghost_columns;
  // Optional nonperiodic structured patch policy. The active recursive TPL
  // meshes are circular and keep this zero.
  int structured_fixed_side_columns;
  double relaxMoveCoef;
  void writeIstr(char *);
  void writeTPL(char *);
  void calcNandK();
  void calcNandKstr();
  void MoveLiquidInterface(double, mesh *);
  void calcRelaxStr();
  void calcRelaxSL();
  void calcFacetAngle();
  void kinMoveInterface(double, mesh *);
  void MoveTPL(double , double , double , int, mesh *);
  void MoveTPLvoronkovGA(double , double , double , int, mesh *);
  void pullSGstr(double);
} BCSet;

typedef struct mesh {
public:
  // Node markers
  static const int SOLID_LIQUID = 0;
  static const int LIQUID_GAS = 1;
  static const int CRUSIBLE_LIQUID_GAS = 2;
  static const int SOLID_LIQUID_GAS = 3;

  static const int BC_COLD = 0;
  static const int BC_HOT = 1;
  static const int BC_LIQUID_GAS = 2;
  static const int BC_SOLID_GAS = 3;
  static const int BC_SOLID_LIQUID = 4;
  static const int BC_BOTTOM = 5;
  static const int BC_TPL = 6;
  static const int BC_SUB_SOLID_LIQUID = 7;
  static const int BC_SUB_LIQUID_GAS = 8;
  static const int BC_SUB_TPL = 9;
  static const int BC_SUB_SOLID_GAS = 10;
  static const int BC_SUB = 11;
  int number_of_bc_sets;
  int number_of_edges;

  static const int LIQUID = 0;
  static const int SOLID = 1;

  int anisotropicKinFlg;
  int rChangeFlg;
  int AniGaFlg;
  int useVoronkovGA;

  int nTPL;
  int number_of_nodes;
  int number_of_elements;
  int number_of_coord_dir;

  int number_of_node_attributes;
  int number_of_node_markers;
  int number_of_groups;
  int number_of_nods_in_group[2];
  int number_of_Pnods;
  int number_of_variables;
  double sdx[3], sdT;
  int sdnode;

  double x_max[3];
  double x_min[3];

  Node *node;
  Node *sub_node;
  Element *element;
  Edge *edge;
  BCSet *bc_set;
  BCSet *bc_setV;

  facet *fac;
  int NumOfFacets;
  int num_of_tri;
  Tri *tri;

  double ZI;
  float theta;
  float theta0;
  double gam;
  double rhog;
  double lgRelaxCfl; // dimensionless explicit Young-Laplace relaxation bound
  double betaRough; // [m/(s K)], rough S/L normal kinetic coefficient
  double betaStep;  // [m/(s K)], lateral wide-terrace step coefficient
  double kineticA2DN;
  double kineticB2DN;
  // Legacy fields retained only so the vector-comparison snapshot still
  // compiles. Active drivers and common movement routines do not use them.
  double maxNodeMoveDx;
  double kinMaxNodeMoveDx;
  double lgMoveCoef;
  double lgMaxNodeMoveDx;
  double voronkovKV;
  double voronkovVelocityToCm;
  double voronkovTm;
  double voronkovQ;
  double voronkovLambdaSG;
  double voronkovAtomicDensity;
  double voronkovAlphaTransition;
  int voronkovTemperatureIsAbsolute;
  double voronkovSigmaSLFacet;
  double voronkovSigmaSLRough;
  double voronkovSigmaSLPrimeFacet;
  double voronkovSigmaSG;
  double voronkovSigmaMG;
  double roughGrowthAngle; // [rad], rough GAeq from free-energy balance
  double facetGrowthAngle; // [rad], stepped GAeq from free-energy balance
  double linearRoughUndercooling; // [K], Vpull/betaRough
  double linearTipUndercooling; // [K], steady 2DN inversion at the tip
  double linearTipChi; // [rad], force-limited tip Voronkov rotation
  double linearFacetAppAngle; // [rad], reduced linear-model tip GAapp
  double voronkovThetaFacet;
  double voronkovThetaStepCorrection;
  int subSlRefinement; // number of radial S/L elements in the local TPL belt
  int subLgRefinement;
  int tplSubmeshSupportLayers;
  int tplSubmeshOuterElements;
  double cDG, dx;
  double Vpull;
  void calc3dmeshDef();
  void calcNandK();
  void calcNtri();
  void calcNbcSet(BCSet *);
  void calcNnode();
  void calcNandK_Meyer();
  void calcNandKlin();
  void CalcGaussJo(double **, int /*, float **b, int m */);
  void calcFacetAngle();
  void calcRelax();
  void calcRelax1();
  void calcRelax2();
  void calcArtificialT(double);

  void kinMoveInterface(double, double, double, double, double);
  void kinMoveInterfaceAni(double, double, double, double, double);

  void LocalGetK(Element *);
  void LocalGetN(Element *);
  void BasisFuncDBF2d6nodLagr(double, double, Coordinate *);

  void initCreateMesh3d(int, double);
  void initCreateMesh3d1(int, double d_crucible);
  void initAdjNodesQuad(BCSet *);
  int initInterfaceFacet111f6();
  int initInterfaceFacet111f3();
  void initXminmax();
  void init(int);
  void init3d(int);
  void initBoundCondAssembly();
  void initAdjNodes1();
  void initAdjNodes();
  void initEdges();
  void initBCEdges(BCSet *);
  void initTri();
  void initTri1();

  void writeDomBRICK_Tec(char *);
  void writeBcQuad_Tec(BCSet *, char *);
  void writeTecplotTet(char *);
  void writeTecplotTri(char *);
  void writeVTKtri(char *);
  void MeshWriteTecplot(char *);
  void MeshWrite2dBoundTecplot(BCSet *, char *);
  void writeTecplotTriCenter(char *);
  void writeReadGAMBITfile(char *);
  void writeReadGAMBITfile1(char *);
  void write(char *);
  void writeRead(char *);
  void writeCGStri(char *);
} Mesh;

#include "calculate.c"
#include "init.c"
#include "kin.c"
#include "local_calc.c"
#include "write_read.c"
