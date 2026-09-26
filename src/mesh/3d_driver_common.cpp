#include "../research/params.cpp"

static void initParam(Mesh *p_mesh, LocalCalc *p_local_calc,
                      const char *filename) {
  const Params p = readParams(filename);

  // LocalCalc contains coefficients used by the FE heat equation and its
  // thermal boundary terms. Material group 0 is liquid and group 1 is solid.
  p_local_calc->Ttop = p.temperatureTop;
  p_local_calc->Tout = p.temperatureOut;
  p_local_calc->Tmp = p.meltingTemperature;
  p_local_calc->Tbottom = p.temperatureBottom;

  // LocalCalc stores the heat equation after division by rho*cp. The
  // external heat-transfer, radiation, and interface-enthalpy terms retain
  // the rhoSi*cpLiquid denominator (2520*946) of the original FE formulation.
  p_local_calc->hh =
      p.heatTransferCoefficient / (p.rhoSi * p.cpLiquid);
  p_local_calc->Rad =
      p.stefanBoltzmann * p.surfaceEmissivity /
      (p.rhoSi * p.cpLiquid);
  p_local_calc->PeS =
      p.latentHeatVolume / (p.rhoSi * p.cpLiquid);

  // Mesh material groups are 0=liquid and 1=solid.
  p_local_calc->cond[Mesh::LIQUID] =
      p.conductivityLiquid / (p.rhoSi * p.cpLiquid);
  p_local_calc->cond[Mesh::SOLID] =
      p.conductivitySolid / (p.rhoSi * p.cpSolid);

  // Mesh contains capillary, S/L kinetic, pulling, and Voronkov TPL
  // parameters used by the interface-motion routines.
  p_local_calc->Vp = p.pullVelocity;
  p_mesh->Vpull = p.pullVelocity;
  p_mesh->useVoronkovGA = (p.useVoronkovGA != 0);
  p_mesh->gam = p.sigmaMG;
  p_mesh->rhog = p.rhoSi * p.g;
  p_mesh->lgRelaxCfl = .5;
  p_mesh->betaRough = p.betaRough;
  p_mesh->betaStep = p.betaStep;
  p_mesh->kineticA2DN = p.kineticA2DN;
  p_mesh->kineticB2DN = p.kineticB2DN;

  // Voronkov's published coefficient is evaluated in cgs. Surface energies
  // read in J/m^2 are converted to erg/cm^2; tau/h has the same conversion.
  const double surfaceEnergyToCgs = 1000.0;
  p_mesh->voronkovTm = p.voronkovTm;
  p_mesh->voronkovQ = p.voronkovQ;
  p_mesh->voronkovLambdaSG = p.voronkovLambdaSG;
  p_mesh->voronkovAtomicDensity = p.voronkovAtomicDensity;
  p_mesh->voronkovTemperatureIsAbsolute =
      p.voronkovTemperatureIsAbsolute;
  p_mesh->voronkovAlphaTransition = p.voronkovAlphaTransition;
  p_mesh->voronkovSigmaSLFacet =
      surfaceEnergyToCgs * p.sigmaSLFacet;
  p_mesh->voronkovSigmaSLRough =
      surfaceEnergyToCgs * p.sigmaSLRough;
  p_mesh->voronkovSigmaSG = surfaceEnergyToCgs * p.sigmaSG;
  p_mesh->voronkovSigmaMG = surfaceEnergyToCgs * p.sigmaMGSurface;
  const double sigmaDiff =
      p_mesh->voronkovSigmaSG - p_mesh->voronkovSigmaMG;
  p_mesh->roughGrowthAngle =
      2.0 * asin(0.5 *
                 sqrt(p_mesh->voronkovSigmaSLRough *
                          p_mesh->voronkovSigmaSLRough -
                      sigmaDiff * sigmaDiff) /
                 sqrt(p_mesh->voronkovSigmaSG *
                      p_mesh->voronkovSigmaMG));
  p_mesh->voronkovSigmaSLPrimeFacet =
      surfaceEnergyToCgs * p.tauStep / p.hStep;
  p_mesh->facetGrowthAngle =
      2.0 * asin(0.5 *
                 sqrt(p_mesh->voronkovSigmaSLFacet *
                          p_mesh->voronkovSigmaSLFacet +
                      p_mesh->voronkovSigmaSLPrimeFacet *
                          p_mesh->voronkovSigmaSLPrimeFacet -
                      sigmaDiff * sigmaDiff) /
                 sqrt(p_mesh->voronkovSigmaSG *
                      p_mesh->voronkovSigmaMG));
  p_mesh->voronkovKV =
      (p.voronkovQ / p.voronkovTm) *
      pow(p.voronkovLambdaSG, 1.0 / 3.0) *
      pow(p_mesh->voronkovSigmaSG * p.voronkovAtomicDensity,
          -2.0 / 3.0);
  p_mesh->voronkovThetaFacet = p.facetSlope;
  p_mesh->voronkovThetaStepCorrection =
      p.voronkovThetaStepCorrection;
  p_mesh->voronkovVelocityToCm = 100.0;

  p_mesh->linearRoughUndercooling =
      (p_mesh->Vpull != 0.0)
          ? fabs(p_mesh->Vpull) / p_mesh->betaRough
          : 0.0;
  p_mesh->linearTipUndercooling = 0.0;
  p_mesh->linearTipChi = 0.0;
  p_mesh->linearFacetAppAngle = p_mesh->facetGrowthAngle;
  const double VnSLTip =
      fabs(p_mesh->Vpull) * sin(fabs(p_mesh->voronkovThetaFacet));
  if (p_mesh->Vpull != 0.0 && VnSLTip > 0.0 &&
      p_mesh->kineticB2DN > 0.0) {
    double dTlo = 1.e-12;
    double dThi = 1.0;
    while (p_mesh->kineticB2DN * pow(dThi, 5. / 6.) *
               exp(-p_mesh->kineticA2DN / dThi) < VnSLTip)
      dThi *= 2.0;
    for (int it(0); it < 100; it++) {
      const double dTmid = .5 * (dTlo + dThi);
      const double V2dn =
          p_mesh->kineticB2DN * pow(dTmid, 5. / 6.) *
          exp(-p_mesh->kineticA2DN / dTmid);
      if (V2dn < VnSLTip)
        dTlo = dTmid;
      else
        dThi = dTmid;
    }
    p_mesh->linearTipUndercooling = .5 * (dTlo + dThi);
    const double alphaTip =
        VnSLTip /
        (p_mesh->betaStep * p_mesh->linearTipUndercooling);
    const double sigmaSLPrime =
        -fabs(p_mesh->voronkovSigmaSLPrimeFacet);
    const double sigmaSLTip =
        p_mesh->voronkovSigmaSLFacet * cos(alphaTip) +
        sigmaSLPrime * sin(alphaTip);
    const double sigmaSLAlphaTip =
        -p_mesh->voronkovSigmaSLFacet * sin(alphaTip) +
        sigmaSLPrime * cos(alphaTip);
    const double SigmaSLTip =
        sqrt(sigmaSLTip * sigmaSLTip +
             sigmaSLAlphaTip * sigmaSLAlphaTip);
    const double psiSLTip = atan2(sigmaSLAlphaTip, sigmaSLTip);
    double cosDeltaS =
        (p_mesh->voronkovSigmaSG * p_mesh->voronkovSigmaSG +
         SigmaSLTip * SigmaSLTip -
         p_mesh->voronkovSigmaMG * p_mesh->voronkovSigmaMG) /
        (2.0 * p_mesh->voronkovSigmaSG * SigmaSLTip);
    cosDeltaS = max(-1.0, min(1.0, cosDeltaS));
    const double phiSTip = psiSLTip + acos(cosDeltaS);
    const double chiMaxTip =
        fabs(p_mesh->voronkovThetaFacet) + alphaTip - phiSTip;
    const double VtSGTip =
        fabs(p_mesh->Vpull) * p_mesh->voronkovVelocityToCm;
    const double chiTip =
        p_mesh->voronkovKV * p_mesh->linearTipUndercooling *
        pow(VtSGTip, -1.0 / 3.0);
    if (chiMaxTip > 0.0)
      p_mesh->linearTipChi = min(chiTip, chiMaxTip);
    p_mesh->linearFacetAppAngle =
        p_mesh->facetGrowthAngle - p_mesh->linearTipChi;
  }
}

static void writeDriverOutput(Mesh *p_mesh, Mesh *p_metrics_mesh,
                              ostream &outw,
                              double time_value, double norm, double GAr,
                              double GAf, const char *outputDir,
                              int roughOnly, int writeDataFiles,
                              int writeHeader) {
  Mesh &mesh = *p_mesh;
  Mesh &metrics = *p_metrics_mesh;
  char filename[100];

  const int ttOut((int)floor(time_value + .555555));
  const int metricsTplBc =
      (p_metrics_mesh == p_mesh) ? metrics.BC_TPL : metrics.BC_SUB_TPL;
  const int metricsSlBc =
      (p_metrics_mesh == p_mesh)
          ? metrics.BC_SOLID_LIQUID
          : metrics.BC_SUB_SOLID_LIQUID;

  if (metrics.useVoronkovGA)
    metrics.bc_set[metricsTplBc].MoveTPLvoronkovGA(
        0., GAr, GAf, metrics.AniGaFlg, &metrics);
  else
    metrics.bc_set[metricsTplBc].MoveTPL(
        0., GAr, GAf, metrics.AniGaFlg, &metrics);

  BCSet *p_tpl_metrics = &metrics.bc_set[metricsTplBc];
  const int tpl_ring_size = p_tpl_metrics->nS1;
  // Fixed symmetry directions: i=0 is the rough base and i=nS1/4 is the
  // selected facet tip. Selecting by nearest GAs can jump between several
  // symmetry-equivalent exact-facet nodes and create false ridge motion.
  BCAtNode *p_rough = p_tpl_metrics->bc_at_nodeStr[0][1];
  BCAtNode *p_facet = p_tpl_metrics->bc_at_nodeStr[tpl_ring_size / 4][1];

  BCSet *p_sl_metrics = &metrics.bc_set[metricsSlBc];
  BCAtNode *p_sl_rough = p_rough;
  BCAtNode *p_sl_facet = p_facet;
  if (p_rough->parent_bc_node_id >= 0 &&
      p_rough->parent_bc_node_id < p_sl_metrics->number_of_recordsNode)
    p_sl_rough =
        &p_sl_metrics->bc_at_node[p_rough->parent_bc_node_id];
  if (p_facet->parent_bc_node_id >= 0 &&
      p_facet->parent_bc_node_id < p_sl_metrics->number_of_recordsNode)
    p_sl_facet =
        &p_sl_metrics->bc_at_node[p_facet->parent_bc_node_id];

  Node *p_node_rough = p_sl_rough->p_node;
  Node *p_node_facet = p_sl_facet->p_node;
  const double Tfe_rough_K = p_node_rough->attribute[TEMPERATURE];
  const double Tfe_facet_K = p_node_facet->attribute[TEMPERATURE];

  const double velocity_to_mm_min = 60000.;
  const double grad_to_K_cm = .01;
  const double coord_to_mm = 1000.;
  const int timeWidth = roughOnly ? 6 : 7;

  const double radius_rough_m =
      sqrt(p_node_rough->coord[0] * p_node_rough->coord[0] +
           p_node_rough->coord[1] * p_node_rough->coord[1]);
  const double radius_facet_m =
      sqrt(p_node_facet->coord[0] * p_node_facet->coord[0] +
           p_node_facet->coord[1] * p_node_facet->coord[1]);
  const double height_rough_m = p_node_rough->coord[2] - metrics.ZI;
  const double height_facet_m = p_node_facet->coord[2] - metrics.ZI;
  // Signed vertical translation speed implied by the local S/L kinetic normal
  // velocity. The calculated normal contributes only through |Nz|: the stored
  // kinetic law has Vn < 0 for solidification and Vn > 0 for melting, hence
  // VtKin > 0 for solidification and VtKin < 0 for melting.
  const double VtKin_rough_m_s =
      -p_sl_rough->Vn / (fabs(p_sl_rough->N[2]) + 1.e-16);
  const double VtKin_facet_m_s =
      -p_sl_facet->Vn / (fabs(p_sl_facet->N[2]) + 1.e-16);

  const double GzS_rough_K_cm = p_node_rough->dTs[2] * grad_to_K_cm;
  const double GrS_rough_K_cm =
      (p_node_rough->dTs[0] * p_node_rough->coord[0] +
       p_node_rough->dTs[1] * p_node_rough->coord[1]) /
      radius_rough_m * grad_to_K_cm;
  const double GzL_rough_K_cm = p_node_rough->dTl[2] * grad_to_K_cm;
  const double GrL_rough_K_cm =
      (p_node_rough->dTl[0] * p_node_rough->coord[0] +
       p_node_rough->dTl[1] * p_node_rough->coord[1]) /
      radius_rough_m * grad_to_K_cm;

  const double GzS_facet_K_cm = p_node_facet->dTs[2] * grad_to_K_cm;
  const double GrS_facet_K_cm =
      (p_node_facet->dTs[0] * p_node_facet->coord[0] +
       p_node_facet->dTs[1] * p_node_facet->coord[1]) /
      radius_facet_m * grad_to_K_cm;
  const double GzL_facet_K_cm = p_node_facet->dTl[2] * grad_to_K_cm;
  const double GrL_facet_K_cm =
      (p_node_facet->dTl[0] * p_node_facet->coord[0] +
       p_node_facet->dTl[1] * p_node_facet->coord[1]) /
      radius_facet_m * grad_to_K_cm;

  ios::fmtflags outw_flags = outw.flags();
  streamsize outw_precision = outw.precision();

  if (writeHeader)
    outw << right
         << setw(1) << "k" << "|"
         << setw(timeWidth) << (roughOnly ? "it" : "t") << "|"
         << setw(8) << "norm" << "|"
         << setw(6) << "GA" << "|"
         << setw(6) << "GAeq" << "|"
         << setw(6) << "GAapp" << "|"
         << setw(8) << "T-Tm" << "|"
         << setw(8) << "VtKin" << "|"
         << setw(8) << "VtSG" << "|"
         << setw(6) << "chi" << "|"
         << setw(7) << "GzS" << "|"
         << setw(7) << "GrS" << "|"
         << setw(7) << "GzL" << "|"
         << setw(7) << "GrL" << "|"
         << setw(8) << "Z" << "|"
         << setw(8) << "R" << "\n";

  outw << right
       << setw(1) << "r" << "|"
       << fixed << setprecision(roughOnly ? 1 : 3)
       << setw(timeWidth) << time_value << "|"
       << defaultfloat << setprecision(3) << setw(8) << sqrt(norm) << "|"
       << fixed << setprecision(2)
       << setw(6) << 180. * p_rough->GA / PI << "|"
       << setw(6) << 180. * p_rough->GAs / PI << "|"
       << setw(6) << 180. * p_rough->GAv / PI << "|"
       << defaultfloat << setprecision(3)
       << setw(8) << Tfe_rough_K << "|"
       << setw(8) << VtKin_rough_m_s * velocity_to_mm_min << "|"
       << setw(8) << fabs(p_rough->V_vo) * velocity_to_mm_min << "|"
       << fixed << setprecision(2)
       << setw(6) << 180. * p_rough->chi_use / PI << "|"
       << defaultfloat << setprecision(3)
       << setw(7) << GzS_rough_K_cm << "|"
       << setw(7) << GrS_rough_K_cm << "|"
       << setw(7) << GzL_rough_K_cm << "|"
       << setw(7) << GrL_rough_K_cm << "|"
       << fixed << setprecision(4)
       << setw(8) << height_rough_m * coord_to_mm << "|"
       << setw(8) << radius_rough_m * coord_to_mm << "\n";

  if (roughOnly) {
    outw.flags(outw_flags);
    outw.precision(outw_precision);
    return;
  }

  outw << setw(1) << "f" << "|"
       << fixed << setprecision(3) << setw(timeWidth) << time_value << "|"
       << defaultfloat << setprecision(3) << setw(8) << sqrt(norm) << "|"
       << fixed << setprecision(2)
       << setw(6) << 180. * p_facet->GA / PI << "|"
       << setw(6) << 180. * p_facet->GAs / PI << "|"
       << setw(6) << 180. * p_facet->GAv / PI << "|"
       << defaultfloat << setprecision(3)
       << setw(8) << Tfe_facet_K << "|"
       << setw(8) << VtKin_facet_m_s * velocity_to_mm_min << "|"
       << setw(8) << fabs(p_facet->V_vo) * velocity_to_mm_min << "|"
       << fixed << setprecision(2)
       << setw(6) << 180. * p_facet->chi_use / PI << "|"
       << defaultfloat << setprecision(3)
       << setw(7) << GzS_facet_K_cm << "|"
       << setw(7) << GrS_facet_K_cm << "|"
       << setw(7) << GzL_facet_K_cm << "|"
       << setw(7) << GrL_facet_K_cm << "|"
       << fixed << setprecision(4)
       << setw(8) << height_facet_m * coord_to_mm << "|"
       << setw(8) << radius_facet_m * coord_to_mm << "\n";

  outw.flags(outw_flags);
  outw.precision(outw_precision);

  if (!writeDataFiles)
    return;

  sprintf(filename, "%s/dom%.6d.dat", outputDir, ttOut);
  mesh.writeDomBRICK_Tec(filename);
  sprintf(filename, "%s/BCcrysMelt%.6d.dat", outputDir, ttOut);
  mesh.writeBcQuad_Tec(&mesh.bc_set[mesh.BC_SOLID_LIQUID], filename);
  sprintf(filename, "%s/BCgasMelt%.6d.dat", outputDir, ttOut);
  mesh.bc_set[mesh.BC_LIQUID_GAS].writeIstr(filename);
  sprintf(filename, "%s/BCgasSolid%.6d.dat", outputDir, ttOut);
  mesh.bc_set[mesh.BC_SOLID_GAS].writeIstr(filename);
  sprintf(filename, "%s/BCtpl%.6d.dat", outputDir, ttOut);
  if (mesh.useVoronkovGA)
    mesh.bc_set[mesh.BC_TPL].MoveTPLvoronkovGA(
        0., GAr, GAf, mesh.AniGaFlg, &mesh);
  else
    mesh.bc_set[mesh.BC_TPL].MoveTPL(
        0., GAr, GAf, mesh.AniGaFlg, &mesh);
  mesh.bc_set[mesh.BC_TPL].writeTPL(filename);
}

/*
  The detailed per-node kinetic, force-balance, and Voronkov diagnostics remain
  in BCtpl*.dat and sub4_TPL*.dat. The console table above is intentionally
  limited to the variables needed to compare the FE phases and the selected
  GA closure.
*/

void BCSet::writeIstr(char *filename) {
  FILE *ifwri;
  ifwri = fopen(filename, "w");
  if (ifwri == NULL) {
    cerr << "Failed to open structured output file: " << filename << "\n";
    exit(1);
  }
  // fprintf(ifwri, "VARIABLES = \"X\", \"Y\", \"Z\", \"f\", \"T\" " );
  fprintf(ifwri, "VARIABLES=\"X\" \"Y\" \"Z\"    \"Nx\" \"Ny\" \"Nz\"     "
                 "\"K\" \"Th\"    \"T\"  \n");
  const bool wrap_export =
      structured_fixed_side_columns == 0 &&
      structured_z_ghost_columns == 0 &&
      ((structured_wrap_columns > 0 && structured_wrap_columns < nS1) ||
       ((mesh_bc_id == Mesh::BC_LIQUID_GAS ||
         mesh_bc_id == Mesh::BC_SOLID_GAS ||
         mesh_bc_id == Mesh::BC_SUB_LIQUID_GAS ||
         mesh_bc_id == Mesh::BC_SUB_SOLID_GAS) &&
        nS1 > 1));
  const int output_cols = structured_z_ghost_columns != 0
                              ? nS1 - 2
                              : (wrap_export ? nS1 + 1 : nS1);
  fprintf(ifwri, "ZONE I=%d J=%d  F=POINT  \n", nS2, output_cols);

  for (int i(0); i < output_cols; i++)
    for (int j(0); j < nS2; j++) {
      const int src_i = structured_z_ghost_columns != 0
                            ? i + 1
                            : ((nS1 > 0) ? (i % nS1) : i);
      fprintf(ifwri, "%e  %e  %e      %e  %e  %e      %e  %e  %e   \n",
              bc_at_nodeStr[src_i][j]->p_node->coord[0],
              bc_at_nodeStr[src_i][j]->p_node->coord[1],
              bc_at_nodeStr[src_i][j]->p_node->coord[2],
              bc_at_nodeStr[src_i][j]->N[0],
              bc_at_nodeStr[src_i][j]->N[1], // bc_at_nodeStr[i][j]->Zf,
              bc_at_nodeStr[src_i][j]->N[2], bc_at_nodeStr[src_i][j]->K,
              bc_at_nodeStr[src_i][j]->Th,
              bc_at_nodeStr[src_i][j]->p_node->attribute[TEMPERATURE]);
    }
  fclose(ifwri);
};

void BCSet::writeTPL(char *filename) {
  FILE *ifwri = fopen(filename, "w");
  if (ifwri == NULL) {
    cerr << "Failed to open TPL output file: " << filename << "\n";
    exit(1);
  }

  if (bc_at_nodeStr == NULL || nS1 < 1 || nS2 < 3) {
    fprintf(ifwri,
            "VARIABLES=\"X_m\" \"Y_m\" \"Z_m\" \"T_minus_Tm_K\" \"Vn_m_s\" \"GA_deg\" \"GAs_deg\" \"GAv_deg\" \"SG_deg\" \"LG_deg\" \"GABranch\" \"dR_m\" \"L_m\" \"V_SL_n_m_s\" \"V_SG_n_local_m_s\" \"V_SG_t_far_m_s\" \"V_SG_t_local_m_s\" \"V_vo_m_s\" \"Vmin_SG_far_cm_s\" \"V_pull_m_s\" \"trace_dot_t_SG_far\" \"DeltaT_K\" \"chi0_abs_rad\" \"chi_geom_rad\" \"chi_use_rad\" \"theta_g_0_rad\" \"theta_g_app_rad\" \"phiS_deg\" \"phiS_mode\" \"alpha_TPL_deg\" \"beta_step_m_s_K\" \"dTloc_K\" \"dTsolid_K\" \"sigmaSG_erg_cm2\" \"sigmaSL_erg_cm2\" \"sigmaLG_erg_cm2\" \"sigmaSL_f_erg_cm2\" \"sigmaSL_alpha_erg_cm2\" \"SigmaSL_CH_erg_cm2\" \"psiSL_CH_deg\" \"deltaS_CH_deg\" \"theta_SG_far_deg\" \"theta_SL_TPL_deg\" \"theta_SG_stop_deg\" \"chiV_deg\" \"chiMax_FB_deg\" \"theta_SG_TPL_deg\"\n");
    fprintf(ifwri, "ZONE I=0 J=1 F=POINT\n");
    fclose(ifwri);
    return;
  }

  double mean_radius = 0.;
  for (int i(0); i < nS1; i++) {
    BCAtNode *p_contact = bc_at_nodeStr[i][1];
    const double x = p_contact->p_node->coord[0];
    const double y = p_contact->p_node->coord[1];
    mean_radius += sqrt(x * x + y * y);
  }
  mean_radius /= (double)nS1;

  fprintf(ifwri,
          "VARIABLES=\"X_m\" \"Y_m\" \"Z_m\" \"T_minus_Tm_K\" \"Vn_m_s\" \"GA_deg\" \"GAs_deg\" \"GAv_deg\" \"SG_deg\" \"LG_deg\" \"GABranch\" \"dR_m\" \"L_m\" \"V_SL_n_m_s\" \"V_SG_n_local_m_s\" \"V_SG_t_far_m_s\" \"V_SG_t_local_m_s\" \"V_vo_m_s\" \"Vmin_SG_far_cm_s\" \"V_pull_m_s\" \"trace_dot_t_SG_far\" \"DeltaT_K\" \"chi0_abs_rad\" \"chi_geom_rad\" \"chi_use_rad\" \"theta_g_0_rad\" \"theta_g_app_rad\" \"phiS_deg\" \"phiS_mode\" \"alpha_TPL_deg\" \"beta_step_m_s_K\" \"dTloc_K\" \"dTsolid_K\" \"sigmaSG_erg_cm2\" \"sigmaSL_erg_cm2\" \"sigmaLG_erg_cm2\" \"sigmaSL_f_erg_cm2\" \"sigmaSL_alpha_erg_cm2\" \"SigmaSL_CH_erg_cm2\" \"psiSL_CH_deg\" \"deltaS_CH_deg\" \"theta_SG_far_deg\" \"theta_SL_TPL_deg\" \"theta_SG_stop_deg\" \"chiV_deg\" \"chiMax_FB_deg\" \"theta_SG_TPL_deg\"\n");
  fprintf(ifwri, "ZONE I=%d J=1 F=POINT\n", nS1);

  double length = 0.;
  for (int i(0); i < nS1; i++) {
    BCAtNode *p_contact = bc_at_nodeStr[i][1];
    if (i > 0) {
      Node *p0 = bc_at_nodeStr[i - 1][1]->p_node;
      Node *p1 = p_contact->p_node;
      const double dx = p1->coord[0] - p0->coord[0];
      const double dy = p1->coord[1] - p0->coord[1];
      const double dz = p1->coord[2] - p0->coord[2];
      length += sqrt(dx * dx + dy * dy + dz * dz);
    }
    const double radius = sqrt(p_contact->p_node->coord[0] *
                                   p_contact->p_node->coord[0] +
                               p_contact->p_node->coord[1] *
                                   p_contact->p_node->coord[1]);
    const double lg_deg = 180. * p_contact->LG / PI;
    const double ga_deg = 180. * p_contact->GA / PI;
    const double gas_deg = 180. * p_contact->GAs / PI;
    const double gav_deg = 180. * p_contact->GAv / PI;
    const double sg_deg = 180. * p_contact->SG / PI;
    const double dr = radius - mean_radius;
    fprintf(ifwri, "%e %e %e %e %e %e %e %e %e %e %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            p_contact->p_node->coord[0], p_contact->p_node->coord[1],
            p_contact->p_node->coord[2],
            p_contact->p_node->attribute[TEMPERATURE], p_contact->Vn,
            ga_deg, gas_deg, gav_deg, sg_deg, lg_deg, p_contact->GABranch,
            dr, length, p_contact->V_SL_n, p_contact->V_SG_n_local,
            p_contact->V_SG_t_far, p_contact->V_SG_t_local,
            p_contact->V_vo, p_contact->Vmin_SG_far, p_contact->V_pull,
            p_contact->trace_dot_t_SG_far, p_contact->DeltaT,
            p_contact->chi0_abs, p_contact->chi_geom, p_contact->chi_use,
            p_contact->theta_g_0, p_contact->theta_g_app,
            180. * p_contact->phiS / PI, p_contact->phiS_mode,
            180. * p_contact->alpha_TPL / PI, p_contact->beta_step,
            p_contact->dTloc, p_contact->dTsolid,
            p_contact->sigmaSG_diag, p_contact->sigmaSL_diag,
            p_contact->sigmaLG_diag, p_contact->sigmaSL_f,
            p_contact->sigmaSL_alpha, p_contact->SigmaSL_CH,
            180. * p_contact->psiSL_CH / PI,
            180. * p_contact->deltaS_CH / PI,
            180. * p_contact->theta_SG_far / PI,
            180. * p_contact->theta_SL_TPL / PI,
            180. * p_contact->theta_SG_stop / PI,
            180. * p_contact->chiV / PI,
            180. * p_contact->chiMax_FB / PI,
            180. * p_contact->theta_SG_TPL / PI);
  }

  fclose(ifwri);
}

void BCSet::calcNandKstr() {
  const int v[8][2] = {{0, 1},  {1, 1},   {1, 0},  {1, -1},
                       {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}};
  const int active_cols = nS1 - structured_wrap_columns;
  const int z_ghost_columns =
      structured_z_ghost_columns != 0 && active_cols > 2;
  const int fixed_side_columns =
      structured_fixed_side_columns != 0 && active_cols > 2;
  const int i_begin = (z_ghost_columns || fixed_side_columns) ? 1 : 0;
  const int i_end = (z_ghost_columns || fixed_side_columns)
                        ? active_cols - 1
                        : active_cols;

  if (active_cols < 1 || nS2 < 2)
    return;

  if (z_ghost_columns) {
#pragma omp parallel for schedule(static)
    for (int j = 0; j < nS2; j++) {
      // Local structured patch: columns 1...n are physical. Only Z is
      // continued through the two ghost columns.
      bc_at_nodeStr[0][j]->p_node->coord[2] =
          bc_at_nodeStr[active_cols - 2][j]->p_node->coord[2];
      bc_at_nodeStr[active_cols - 1][j]->p_node->coord[2] =
          bc_at_nodeStr[1][j]->p_node->coord[2];
    }
  }

#pragma omp parallel for collapse(2) schedule(static)
  for (int i = i_begin; i < i_end; i++)
    for (int j = 1; j < nS2 - 1; j++) {
      double a[3], b[3], c[3], n[3], nK[3];
      double area(0.);
      double cotangentSum = 0.;
      int i0(i);
      bc_at_nodeStr[i0][j]->K = 0.;
      for (int co = 0; co < 3; co++) {
        bc_at_nodeStr[i0][j]->N[co] = 0.;
        nK[co] = 0.;
      }

      for (int d(0); d < 8; d++) {
        int dn(d + 1);
        if (dn > 7)
          dn = 0;
        int ia = (i0 + v[d][0] + active_cols) % active_cols;
        int ib = (i0 + v[dn][0] + active_cols) % active_cols;

        a[0] = -bc_at_nodeStr[i0][j]->p_node->coord[0] +
               bc_at_nodeStr[ia][j + v[d][1]]->p_node->coord[0];
        a[1] = -bc_at_nodeStr[i0][j]->p_node->coord[1] +
               bc_at_nodeStr[ia][j + v[d][1]]->p_node->coord[1];
        a[2] = -bc_at_nodeStr[i0][j]->p_node->coord[2] +
               bc_at_nodeStr[ia][j + v[d][1]]->p_node->coord[2];

        b[0] = -bc_at_nodeStr[i0][j]->p_node->coord[0] +
               bc_at_nodeStr[ib][j + v[dn][1]]->p_node->coord[0];
        b[1] = -bc_at_nodeStr[i0][j]->p_node->coord[1] +
               bc_at_nodeStr[ib][j + v[dn][1]]->p_node->coord[1];
        b[2] = -bc_at_nodeStr[i0][j]->p_node->coord[2] +
               bc_at_nodeStr[ib][j + v[dn][1]]->p_node->coord[2];

        c[0] = bc_at_nodeStr[ia][j + v[d][1]]->p_node->coord[0] -
               bc_at_nodeStr[ib][j + v[dn][1]]->p_node->coord[0];
        c[1] = bc_at_nodeStr[ia][j + v[d][1]]->p_node->coord[1] -
               bc_at_nodeStr[ib][j + v[dn][1]]->p_node->coord[1];
        c[2] = bc_at_nodeStr[ia][j + v[d][1]]->p_node->coord[2] -
               bc_at_nodeStr[ib][j + v[dn][1]]->p_node->coord[2];

        n[0] = a[1] * c[2] - a[2] * c[1];
        n[1] = a[2] * c[0] - a[0] * c[2];
        n[2] = a[0] * c[1] - a[1] * c[0];

        const double tri_norm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (tri_norm <= 1.e-16)
          continue;
        double cota = (a[0] * c[0] + a[1] * c[1] + a[2] * c[2]) / tri_norm;
        double cotb = (b[0] * c[0] + b[1] * c[1] + b[2] * c[2]) / tri_norm;
        cotangentSum += fabs(cota) + fabs(cotb);

        for (int co = 0; co < 3; co++) {
          nK[co] += a[co] * cotb - b[co] * cota;
          bc_at_nodeStr[i0][j]->N[co] += n[co];
        }
        area += .5 * tri_norm;
      }

      double Ll;

      Ll = sqrt(nK[0] * nK[0] + nK[1] * nK[1] + nK[2] * nK[2]);
      bc_at_nodeStr[i0][j]->K = (area > 1.e-16) ? (-2. * Ll / area) : 0.;
      bc_at_nodeStr[i0][j]->curvatureRate =
          (area > 1.e-16) ? 2. * cotangentSum / area : 0.;

      Ll = (sqrt(bc_at_nodeStr[i0][j]->N[0] * bc_at_nodeStr[i0][j]->N[0] +
                 bc_at_nodeStr[i0][j]->N[1] * bc_at_nodeStr[i0][j]->N[1] +
                 bc_at_nodeStr[i0][j]->N[2] * bc_at_nodeStr[i0][j]->N[2]));
      if (Ll > 1.e-16)
        for (int co = 0; co < 3; co++)
          bc_at_nodeStr[i0][j]->N[co] /= -Ll;
      else {
        bc_at_nodeStr[i0][j]->N[0] = 0.;
        bc_at_nodeStr[i0][j]->N[1] = 0.;
        bc_at_nodeStr[i0][j]->N[2] = 0.;
      }

      if (nK[0] * bc_at_nodeStr[i0][j]->N[0] +
              nK[1] * bc_at_nodeStr[i0][j]->N[1] +
              nK[2] * bc_at_nodeStr[i0][j]->N[2] <
          0)
        bc_at_nodeStr[i0][j]->K *= -1;
    }

#pragma omp parallel for schedule(static)
  for (int i = 0; i < active_cols; i++)
    bc_at_nodeStr[i][0]->K = bc_at_nodeStr[i][1]->K;

#pragma omp parallel for schedule(static)
  for (int i = active_cols; i < nS1; i++) {
    const int src_i = i % active_cols;
    for (int j(0); j < nS2; j++) {
      bc_at_nodeStr[i][j]->K = bc_at_nodeStr[src_i][j]->K;
      for (int co = 0; co < 3; co++)
        bc_at_nodeStr[i][j]->N[co] = bc_at_nodeStr[src_i][j]->N[co];
    }
  }
};

void BCSet::calcNandK() {
#pragma omp parallel for schedule(static)
  for (int ib = 0; ib < number_of_recordsNode; ib++) {
    double a[3], b[3], c[3], n[3], nK[3];
    double area(0.);
    bc_at_node[ib].K = 0.;
    for (int co = 0; co < 3; co++) {
      bc_at_node[ib].N[co] = 0.;
      nK[co] = 0.;
    }
    for (int d(0); d < bc_at_node[ib].num_of_adj_elem; d++) {
      int i0(0);
      int ne(bc_at_node[ib].n_adj_elem[d]);
      for (int i(0); i < 4; i++)
        if (bc_at_face[ne].assembly[i] == ib)
          i0 = i;

      for (int i(0); i < 2; i++) {
        int i1(i0 + 1);
        if (i1 > 3)
          i1 = 0;
        int i2(i1 + 1);
        if (i2 > 3)
          i2 = 0;
        int ib1(bc_at_face[ne].assembly[i1]);
        int ib2(bc_at_face[ne].assembly[i2]);

        for (int co = 0; co < 3; co++) {
          a[co] = -bc_at_node[ib].p_node->coord[co] +
                  bc_at_node[ib1].p_node->coord[co];
          b[co] = -bc_at_node[ib].p_node->coord[co] +
                  bc_at_node[ib2].p_node->coord[co];
          c[co] = -bc_at_node[ib1].p_node->coord[co] +
                  bc_at_node[ib2].p_node->coord[co];
        }
        n[0] = a[1] * c[2] - a[2] * c[1];
        n[1] = a[2] * c[0] - a[0] * c[2];
        n[2] = a[0] * c[1] - a[1] * c[0];

        const double tri_norm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (tri_norm <= 1.e-16)
          continue;
        double cota = (a[0] * c[0] + a[1] * c[1] + a[2] * c[2]) / tri_norm;
        double cotb = (b[0] * c[0] + b[1] * c[1] + b[2] * c[2]) / tri_norm;

        for (int co = 0; co < 3; co++) {
          nK[co] += a[co] * cotb - b[co] * cota;
          bc_at_node[ib].N[co] += n[co];
        }
        area += .5 * tri_norm;
      }
    }
    double Ll(sqrt(nK[0] * nK[0] + nK[1] * nK[1] + nK[2] * nK[2]));
    if ((bc_at_node[ib].p_node->bound_mark == Mesh::BC_SOLID_LIQUID ||
         bc_at_node[ib].p_node->bound_mark == Mesh::BC_SUB_SOLID_LIQUID) &&
        area > 1.e-16)
      bc_at_node[ib].K = -Ll / (area);
    Ll = (sqrt(bc_at_node[ib].N[0] * bc_at_node[ib].N[0] +
               bc_at_node[ib].N[1] * bc_at_node[ib].N[1] +
               bc_at_node[ib].N[2] * bc_at_node[ib].N[2]));
    if (Ll > 1.e-16)
      for (int co = 0; co < 3; co++)
        bc_at_node[ib].N[co] /= -Ll;
    else {
      bc_at_node[ib].N[0] = 0.;
      bc_at_node[ib].N[1] = 0.;
      bc_at_node[ib].N[2] = 0.;
    }
    if (nK[0] * bc_at_node[ib].N[0] + nK[1] * bc_at_node[ib].N[1] +
            nK[2] * bc_at_node[ib].N[2] <
        0)
      bc_at_node[ib].K *= -1;
    if (!isfinite(bc_at_node[ib].K))
      bc_at_node[ib].K = 0.;
  }
};

void BCSet::calcFacetAngle() {
  double Fac[] = {0., -sqrt(2.) / sqrt(3.), 1. / sqrt(3.)};
#pragma omp parallel for schedule(static)
  for (int n_node = 0; n_node < number_of_recordsNode; n_node++) {
    bc_at_node[n_node].Zf = bc_at_node[n_node].p_node->coord[0] * Fac[0] +
                            bc_at_node[n_node].p_node->coord[1] * Fac[1] +
                            bc_at_node[n_node].p_node->coord[2] * Fac[2];
  }

#pragma omp parallel for schedule(static)
  for (int n_node = 0; n_node < number_of_recordsNode; n_node++) {
    double a[3];
    int dd;
    double hh;
    bc_at_node[n_node].Th = 0.;
    dd = (-1);
    hh = (0.);
    for (int n_adj(0); n_adj < bc_at_node[n_node].num_of_adj_nodes; n_adj++) {
      int adj_num = bc_at_node[n_node].n_adj_node[n_adj];

      a[0] = bc_at_node[n_node].p_node->coord[0] -
             bc_at_node[adj_num].p_node->coord[0];
      a[1] = bc_at_node[n_node].p_node->coord[1] -
             bc_at_node[adj_num].p_node->coord[1];
      a[2] = bc_at_node[n_node].p_node->coord[2] -
             bc_at_node[adj_num].p_node->coord[2];
      double ll = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

      double Th(-(bc_at_node[adj_num].Zf - bc_at_node[n_node].Zf) / ll);
      if (hh < Th) {
        hh = Th;
        dd = n_adj;
      }
    }
    
    bc_at_node[n_node].i_adj_n=dd;
    if (dd != -1) {
        bc_at_node[n_node].Th = asin(hh);
    }
  }
};

void BCSet::MoveLiquidInterface(double dt, Mesh *p_mesh) {

  calcNandKstr();

  const int active_cols = nS1 - structured_wrap_columns;
  const int z_ghost_columns =
      structured_z_ghost_columns != 0 && active_cols > 2;
  const int fixed_side_columns =
      structured_fixed_side_columns != 0 && active_cols > 2;
  const int i_begin = (z_ghost_columns || fixed_side_columns) ? 1 : 0;
  const int i_end = (z_ghost_columns || fixed_side_columns)
                        ? active_cols - 1
                        : active_cols;
  if (active_cols < 1 || nS2 < 3)
    return;

  double maxRelaxRate = 0.;
#pragma omp parallel for collapse(2) schedule(static) reduction(max : maxRelaxRate)
  for (int i = i_begin; i < i_end; i++)
    // j=0 is the TPL and j=nS2-1 is the parent-controlled outer cut.
    // Structured relaxation redistributes only the intervening nodes.
    for (int j = 1; j < nS2 - 1; j++) {
      BCAtNode *p_bc_node = bc_at_nodeStr[i][j];
      // Each structured i-column is a fixed azimuthal ray.  Its TPL endpoint
      // defines the radial direction; TPL motion itself is radial.
      const Node *p_path_node = bc_at_nodeStr[i][0]->p_node;
      const double path_radius =
          sqrt(p_path_node->coord[0] * p_path_node->coord[0] +
               p_path_node->coord[1] * p_path_node->coord[1]);
      const double erx = p_path_node->coord[0] / path_radius;
      const double ery = p_path_node->coord[1] / path_radius;
      // The L/G surface is relaxed to Young-Laplace equilibrium, not advanced
      // with a physical interface kinetic law.  Dividing the pressure residual
      // by rho*g gives its geometric height residual [m].  The supplied dt is
      // therefore the explicit relaxation step for this quasi-static surface.
      // calcNandKstr() gives K>0 for the raised Cz meniscus.  The hydrostatic
      // and Laplace pressures balance as rho*g*(z-ZI) - sigma*K = 0.
      const double V =
          -(p_bc_node->p_node->coord[2] - p_mesh->ZI) +
          p_mesh->gam * p_bc_node->K / p_mesh->rhog;
      const double normal_radial =
          p_bc_node->N[0] * erx + p_bc_node->N[1] * ery;
      const double normal_in_path_sq =
          normal_radial * normal_radial +
          p_bc_node->N[2] * p_bc_node->N[2];
      p_bc_node->Vn = V;
      p_bc_node->p_node->Vn = V;
      // Scale the permitted in-plane displacement so that dX dot N = V.
      p_bc_node->dcoord[0] =
          V * normal_radial * erx / normal_in_path_sq;
      p_bc_node->dcoord[1] =
          V * normal_radial * ery / normal_in_path_sq;
      p_bc_node->dcoord[2] = V * p_bc_node->N[2] / normal_in_path_sq;
      p_bc_node->move_abs = fabs(V) / sqrt(normal_in_path_sq);
      const double relaxRate =
          (1. + p_mesh->gam / p_mesh->rhog * p_bc_node->curvatureRate) /
          sqrt(normal_in_path_sq);
      maxRelaxRate = max(maxRelaxRate, relaxRate);
    }

  // Explicit curvature relaxation scales with the square of local spacing.
  // Use one global step on this surface, including both radial and azimuthal
  // cotangent weights; keep all free nodes and the Young-Laplace residual.
  dt = min(dt, p_mesh->lgRelaxCfl / maxRelaxRate);

#pragma omp parallel for collapse(2) schedule(static)
  for (int i = i_begin; i < i_end; i++)
    for (int j = 1; j < nS2 - 1; j++) {
      BCAtNode *p_bc_node = bc_at_nodeStr[i][j];
      p_bc_node->p_node->coord[0] += dt * p_bc_node->dcoord[0];
      p_bc_node->p_node->coord[1] += dt * p_bc_node->dcoord[1];
      p_bc_node->p_node->coord[2] += dt * p_bc_node->dcoord[2];
      const Node *p_path_node = bc_at_nodeStr[i][0]->p_node;
      const double path_radius =
          sqrt(p_path_node->coord[0] * p_path_node->coord[0] +
               p_path_node->coord[1] * p_path_node->coord[1]);
      const double node_radius =
          sqrt(p_bc_node->p_node->coord[0] * p_bc_node->p_node->coord[0] +
               p_bc_node->p_node->coord[1] * p_bc_node->p_node->coord[1]);
      p_bc_node->p_node->coord[0] =
          node_radius * p_path_node->coord[0] / path_radius;
      p_bc_node->p_node->coord[1] =
          node_radius * p_path_node->coord[1] / path_radius;
    }

  if (z_ghost_columns) {
#pragma omp parallel for schedule(static)
    for (int j = 0; j < nS2; j++) {
      bc_at_nodeStr[0][j]->p_node->coord[2] =
          bc_at_nodeStr[active_cols - 2][j]->p_node->coord[2];
      bc_at_nodeStr[active_cols - 1][j]->p_node->coord[2] =
          bc_at_nodeStr[1][j]->p_node->coord[2];
    }
  }
}

void BCSet::calcRelaxStr() {

  const int active_cols = nS1 - structured_wrap_columns;
  const int z_ghost_columns =
      structured_z_ghost_columns != 0 && active_cols > 2;
  const int fixed_side_columns =
      structured_fixed_side_columns != 0 && active_cols > 2;
  const int i_begin = (z_ghost_columns || fixed_side_columns) ? 1 : 0;
  const int i_end = (z_ghost_columns || fixed_side_columns)
                        ? active_cols - 1
                        : active_cols;
  if (active_cols < 1 || nS2 < 3)
    return;

#pragma omp parallel for collapse(2) schedule(static)
  for (int i = i_begin; i < i_end; i++)
    for (int j = 1; j < nS2 - 1; j++) {
      BCAtNode *p_bc_node = bc_at_nodeStr[i][j];
      // Redistribute a structured column in its fixed radial-z plane.  Using
      // radial distances here prevents mesh relaxation from changing azimuth.
      const Node *p_path_node = bc_at_nodeStr[i][0]->p_node;
      const double path_radius =
          sqrt(p_path_node->coord[0] * p_path_node->coord[0] +
               p_path_node->coord[1] * p_path_node->coord[1]);
      const double erx = p_path_node->coord[0] / path_radius;
      const double ery = p_path_node->coord[1] / path_radius;
      const double radius =
          sqrt(bc_at_nodeStr[i][j]->p_node->coord[0] *
                   bc_at_nodeStr[i][j]->p_node->coord[0] +
               bc_at_nodeStr[i][j]->p_node->coord[1] *
                   bc_at_nodeStr[i][j]->p_node->coord[1]);
      const double radius_forward =
          sqrt(bc_at_nodeStr[i][j + 1]->p_node->coord[0] *
                   bc_at_nodeStr[i][j + 1]->p_node->coord[0] +
               bc_at_nodeStr[i][j + 1]->p_node->coord[1] *
                   bc_at_nodeStr[i][j + 1]->p_node->coord[1]);
      const double drf = radius - radius_forward;
      const double dzf = bc_at_nodeStr[i][j]->p_node->coord[2] -
                         bc_at_nodeStr[i][j + 1]->p_node->coord[2];
      const double llf = sqrt(drf * drf + dzf * dzf);

      const double radius_backward =
          sqrt(bc_at_nodeStr[i][j - 1]->p_node->coord[0] *
                   bc_at_nodeStr[i][j - 1]->p_node->coord[0] +
               bc_at_nodeStr[i][j - 1]->p_node->coord[1] *
                   bc_at_nodeStr[i][j - 1]->p_node->coord[1]);
      const double drb = radius - radius_backward;
      const double dzb = bc_at_nodeStr[i][j]->p_node->coord[2] -
                         bc_at_nodeStr[i][j - 1]->p_node->coord[2];
      const double llb = sqrt(drb * drb + dzb * dzb);

      const double relax_len = llf + llb;

      double relax_coef = 1.;
      //if (bc_at_nodeStr[i][j]->relax_coef > 0.) relax_coef = bc_at_nodeStr[i][j]->relax_coef;
      const double k = .2 * relax_coef * fabs(llf - llb) / relax_len;

      if (llf > llb) {
        p_bc_node->dcoord[0] = -drf * k * erx;
        p_bc_node->dcoord[1] = -drf * k * ery;
        p_bc_node->dcoord[2] = -dzf * k;
      } else {
        p_bc_node->dcoord[0] = -drb * k * erx;
        p_bc_node->dcoord[1] = -drb * k * ery;
        p_bc_node->dcoord[2] = -dzb * k;
      }
    }

#pragma omp parallel for collapse(2) schedule(static)
  for (int i = i_begin; i < i_end; i++)
    for (int j = 1; j < nS2 - 1; j++) {
      BCAtNode *p_bc_node = bc_at_nodeStr[i][j];
      p_bc_node->p_node->coord[0] += p_bc_node->dcoord[0];
      p_bc_node->p_node->coord[1] += p_bc_node->dcoord[1];
      p_bc_node->p_node->coord[2] += p_bc_node->dcoord[2];
      const Node *p_path_node = bc_at_nodeStr[i][0]->p_node;
      const double path_radius =
          sqrt(p_path_node->coord[0] * p_path_node->coord[0] +
               p_path_node->coord[1] * p_path_node->coord[1]);
      const double node_radius =
          sqrt(p_bc_node->p_node->coord[0] * p_bc_node->p_node->coord[0] +
               p_bc_node->p_node->coord[1] * p_bc_node->p_node->coord[1]);
      p_bc_node->p_node->coord[0] =
          node_radius * p_path_node->coord[0] / path_radius;
      p_bc_node->p_node->coord[1] =
          node_radius * p_path_node->coord[1] / path_radius;
    }

  if (z_ghost_columns) {
#pragma omp parallel for schedule(static)
    for (int j = 0; j < nS2; j++) {
      bc_at_nodeStr[0][j]->p_node->coord[2] =
          bc_at_nodeStr[active_cols - 2][j]->p_node->coord[2];
      bc_at_nodeStr[active_cols - 1][j]->p_node->coord[2] =
          bc_at_nodeStr[1][j]->p_node->coord[2];
    }
  }
}

void BCSet::pullSGstr(double dZ) {

  const int i_begin =
      (structured_fixed_side_columns != 0 && nS1 > 2) ? 1 : 0;
  const int i_end =
      (structured_fixed_side_columns != 0 && nS1 > 2) ? nS1 - 1 : nS1;
  for (int i(i_begin); i < i_end ; i++)
  {
      Node *p_path_node = bc_at_nodeStr[i][0]->p_node;
      const double path_radius =
          sqrt(p_path_node->coord[0] * p_path_node->coord[0] +
               p_path_node->coord[1] * p_path_node->coord[1]);
      const double erx = p_path_node->coord[0] / path_radius;
      const double ery = p_path_node->coord[1] / path_radius;
      bc_at_nodeStr[i][0]->p_node->K = bc_at_nodeStr[i][0]->p_node->attribute[TEMPERATURE];
      for (int j(nS2-1); j > 0  ; j--)
      {
        if (bc_at_nodeStr[i][j]->p_node->bound_mark == Mesh::BC_SUB)
          continue;
       const double radius =
           sqrt(bc_at_nodeStr[i][j]->p_node->coord[0] *
                    bc_at_nodeStr[i][j]->p_node->coord[0] +
                bc_at_nodeStr[i][j]->p_node->coord[1] *
                    bc_at_nodeStr[i][j]->p_node->coord[1]);
       const double radius_back =
           sqrt(bc_at_nodeStr[i][j - 1]->p_node->coord[0] *
                    bc_at_nodeStr[i][j - 1]->p_node->coord[0] +
                bc_at_nodeStr[i][j - 1]->p_node->coord[1] *
                    bc_at_nodeStr[i][j - 1]->p_node->coord[1]);
       const double drb = -(radius - radius_back);
       double dzb = -(bc_at_nodeStr[i][j]->p_node->coord[2] -
                           bc_at_nodeStr[i][j - 1]->p_node->coord[2]);
     // if (dzb < .000001) (dzb) = .000001;
      
      
      bc_at_nodeStr[i][j]->p_node->K -= (bc_at_nodeStr[i][j-1]->p_node->K 
                                        - bc_at_nodeStr[i][j]->p_node->K) * dZ / dzb;
      const double radius_new = radius - drb * dZ / dzb;
      bc_at_nodeStr[i][j]->p_node->coord[0] = radius_new * erx;
      bc_at_nodeStr[i][j]->p_node->coord[1] = radius_new * ery;
    }
  }
}

void BCSet::calcRelaxSL() {
  calcNandK();

#pragma omp parallel for schedule(static)
  for (int n_node = 0; n_node < number_of_recordsNode; n_node++) {
    for (int co = 0; co < 3; co++) {
      bc_at_node[n_node].relax_sum[co] = 0.;
      bc_at_node[n_node].dcoord[co] = 0.;
    }
    for (int n_adj(0); n_adj < bc_at_node[n_node].num_of_adj_nodes; n_adj++) {
      const int adj_num = bc_at_node[n_node].n_adj_node[n_adj];
      double coef = bc_at_node[adj_num].relax_coef;
      for (int co(0); co < 3; co++)
        bc_at_node[n_node].relax_sum[co] +=
            coef*(bc_at_node[n_node].p_node->coord[co] -
             bc_at_node[adj_num].p_node->coord[co]);
    }
  }

#pragma omp parallel for schedule(static)
  for (int n_node = 0; n_node < number_of_recordsNode; n_node++) {
    const double Nd = bc_at_node[n_node].relax_sum[0] *
                          bc_at_node[n_node].N[0] +
                      bc_at_node[n_node].relax_sum[1] *
                          bc_at_node[n_node].N[1] +
                      bc_at_node[n_node].relax_sum[2] *
                          bc_at_node[n_node].N[2];
    bc_at_node[n_node].dcoord[0] =
        -relaxMoveCoef * (bc_at_node[n_node].relax_sum[0] -
                 bc_at_node[n_node].N[0] * Nd);
    bc_at_node[n_node].dcoord[1] =
        -relaxMoveCoef * (bc_at_node[n_node].relax_sum[1] -
                 bc_at_node[n_node].N[1] * Nd);
    bc_at_node[n_node].dcoord[2] =
        -relaxMoveCoef * (bc_at_node[n_node].relax_sum[2] -
                 bc_at_node[n_node].N[2] * Nd);
  }

#pragma omp parallel for schedule(static)
  for (int n_node = 0; n_node < number_of_recordsNode; n_node++) {
    // In the local S/L belt, logical row j=0 is the physical TPL. Its
    // position is controlled by MoveTPL and S/L kinetics, not mesh-quality
    // relaxation.  The opposite BC_SUB endpoint is also excluded by the
    // bound_mark test below and remains fixed by the immediate parent.
    if (mesh_bc_id == Mesh::BC_SUB_SOLID_LIQUID &&
        bc_at_node[n_node].p_node->local_j == 0)
      continue;
    if (bc_at_node[n_node].p_node->bound_mark == Mesh::BC_SOLID_LIQUID ||
        bc_at_node[n_node].p_node->bound_mark == Mesh::BC_SUB_SOLID_LIQUID) {
      bc_at_node[n_node].p_node->coord[0] += bc_at_node[n_node].dcoord[0];
      bc_at_node[n_node].p_node->coord[1] += bc_at_node[n_node].dcoord[1];
      bc_at_node[n_node].p_node->coord[2] += bc_at_node[n_node].dcoord[2];
    }
  }
}

void BCSet::kinMoveInterface(double dt, Mesh *p_mesh) {

  calcNandK();
  if (p_mesh->anisotropicKinFlg != 0) calcFacetAngle();

#pragma omp parallel for schedule(static)
  for (int n_node = 0; n_node < number_of_recordsNode; n_node++) {
    if (bc_at_node[n_node].p_node->bound_mark == Mesh::BC_SUB) {
      bc_at_node[n_node].Vn = 0.;
      bc_at_node[n_node].dcoord[0] = 0.;
      bc_at_node[n_node].dcoord[1] = 0.;
      bc_at_node[n_node].dcoord[2] = 0.;
      bc_at_node[n_node].move_abs = 0.;
      continue;
    }
    // FE convention: Tfe = T - Tmp.
    // Kinetic convention: dTkin = Tmp - Ti = -Tfe.
    // Do not clamp dTkin here. Its sign selects growth or melting.
    const double Tfe = bc_at_node[n_node].p_node->attribute[TEMPERATURE];
    const double dTkin = -Tfe;
    double Vn = 0.0;
    double dVn_dTfe = 0.0;

    if (dTkin > 0.0) {
      // Vn is the S/L velocity relative to the translating crystal.  In the
      // current Cz mesh convention solidification is Vn < 0, toward lower Z.
      const double dTg = dTkin;
      Vn = -p_mesh->betaRough * dTg;
      dVn_dTfe = p_mesh->betaRough;

      if (p_mesh->anisotropicKinFlg != 0) {
        const double betaKin = p_mesh->betaStep * bc_at_node[n_node].Th;
        const double V2dn =
            p_mesh->kineticB2DN * pow(dTg, 5. / 6.) *
            exp(-p_mesh->kineticA2DN / dTg);
        double Vkin = betaKin * dTg;
        double betaActive = betaKin;

        // Old kinetic ordering: 2DN supplies the lower velocity and the rough
        // interface mobility supplies the upper velocity.
        if (Vkin < V2dn) {
          Vkin = V2dn;
          betaActive =
              V2dn * (5. / (6. * dTg) +
                       p_mesh->kineticA2DN / (dTg * dTg));
        }
        if (Vkin > p_mesh->betaRough * dTg) {
          Vkin = p_mesh->betaRough * dTg;
          betaActive = p_mesh->betaRough;
        }
        Vn = -Vkin;
        dVn_dTfe = betaActive;
      }
    } else if (dTkin < 0.0) {
      // Melting/overheating. Use reversible rough-interface melting only;
      // growth-only faceted SD and 2DN laws are not evaluated.
      const double dTm = -dTkin;
      Vn = p_mesh->betaRough * dTm;
      dVn_dTfe = p_mesh->betaRough;
    }

    bc_at_node[n_node].dcoord[0] = 0.;
    bc_at_node[n_node].dcoord[1] = 0.;
    // Lab-frame node velocity = motion relative to the crystal + positive-Z
    // crystal pull.  A growing node can therefore still move upward when the
    // pull is faster than the projected negative crystallization velocity.
    const double NzAbs = fabs(bc_at_node[n_node].N[2]) + 1.e-16;
    // In this Cz configuration the crystal is above the melt and its resolved
    // axial thermal direction is negative. Recovered nodal gradients can flip
    // sign locally while a refined element is being deformed; using that
    // transient sign in the Taylor update creates nonphysical positive
    // feedback. Retain the FE magnitude but impose the known physical axial
    // direction for the local implicit kinetic step.
    const double Gz = -fabs(bc_at_node[n_node].p_node->dTs[2]);
    const double denominator =
        1. - dt * dVn_dTfe * Gz / NzAbs;
    bc_at_node[n_node].dcoord[2] =
        (Vn / NzAbs + p_mesh->Vpull) / denominator;
    bc_at_node[n_node].Vn =
        (bc_at_node[n_node].dcoord[2] - p_mesh->Vpull) * NzAbs;
    bc_at_node[n_node].move_abs = fabs(bc_at_node[n_node].dcoord[2]);
  }

#pragma omp parallel for schedule(static)
  for (int n_node = 0; n_node < number_of_recordsNode; n_node++) {
    if (bc_at_node[n_node].p_node->bound_mark != Mesh::BC_SUB) {
      const double dz = dt * bc_at_node[n_node].dcoord[2];
      const double Gz = -fabs(bc_at_node[n_node].p_node->dTs[2]);
      bc_at_node[n_node].p_node->coord[2] += dz;
      bc_at_node[n_node].p_node->attribute[TEMPERATURE] += Gz * dz;
      bc_at_node[n_node].p_node->T =
          bc_at_node[n_node].p_node->attribute[TEMPERATURE];
    }
  }
};


void BCSet::MoveTPL(double coef, double GAr, double GAf, int AniGaFlg,
                    Mesh *p_mesh) {
  // Simplified static GA routine.  This path uses a linear estimate between
  // rough and faceted static angles; MoveTPLvoronkovGA() remains available for
  // manual comparison with the full Voronkov correction.
  const double GArRad = PI * GAr / 180.;
  const double GAfRad = PI * GAf / 180.;
  double delt, xTPL[3], xLG[3], xLG2[3], xSG[3];
  BCSet *p_sl_set = &p_mesh->bc_set[p_mesh->BC_SOLID_LIQUID];
  BCSet *p_lg_set = &p_mesh->bc_set[p_mesh->BC_LIQUID_GAS];
  if (mesh_bc_id == p_mesh->BC_SUB_TPL)
    p_sl_set = &p_mesh->bc_set[p_mesh->BC_SUB_SOLID_LIQUID];
  if (mesh_bc_id == p_mesh->BC_SUB_TPL)
    p_lg_set = &p_mesh->bc_set[p_mesh->BC_SUB_LIQUID_GAS];

  // The TPL tangent uses the TPL node and the first two L/G nodes. One
  // current-state evaluation is therefore sufficient here; the physical main
  // loop supplies the repeated TPL motion through its movement substeps.
  for (int i(0); i < nS1; i++) {
      BCAtNode *p_contact = bc_at_nodeStr[i][1];
      BCAtNode *p_sl =
          &p_sl_set->bc_at_node[p_contact->parent_bc_node_id];
      // p_contact and p_sl can be different mesh records. The parent id is set
      // from initialized S/L geometry, so TPL uses S/L kinetics at the same
      // physical contact-line position instead of assuming equal array indexes.
      p_contact->K = p_sl->K;
      p_contact->Vn = p_sl->Vn;
      p_contact->Th = 0.;
      if (AniGaFlg == 1)
        p_contact->Th = p_sl->Th;

      for (int co = 0; co < 3; co++) {
        xTPL[co] = p_contact->p_node->coord[co];
        xLG[co] = p_contact->p_adj_lg_node->coord[co];
        xLG2[co] =
            p_lg_set
                ->bc_at_nodeStr[i +
                                (p_lg_set->structured_z_ghost_columns != 0)]
                               [2]
                ->p_node->coord[co];
        xSG[co] = p_contact->p_adj_sg_node->coord[co];
      }
      double dz(xTPL[2] - xLG[2]);
      double rTPL(sqrt(xTPL[0] * xTPL[0] + xTPL[1] * xTPL[1]));
      double rLG(sqrt(xLG[0] * xLG[0] + xLG[1] * xLG[1]));
      double rLG2(sqrt(xLG2[0] * xLG2[0] + xLG2[1] * xLG2[1]));
      double rSG(sqrt(xSG[0] * xSG[0] + xSG[1] * xSG[1]));
      double dr(rLG - rTPL);
      const double dsLG1 = sqrt(dr * dr + dz * dz);
      const double lgChord1 = asin(dr / dsLG1);
      const double drLG2 = rLG2 - rLG;
      const double dzLG2 = xLG[2] - xLG2[2];
      const double dsLG2 = sqrt(drLG2 * drLG2 + dzLG2 * dzLG2);
      const double lgChord2 = asin(drLG2 / dsLG2);
      // A chord gives the inclination averaged over one L/G edge.  Extrapolate
      // the first two chord inclinations linearly in arc length to obtain the
      // endpoint tangent used by the TPL growth-angle condition.
      p_contact->LG =
          lgChord1 - (lgChord2 - lgChord1) * dsLG1 / (dsLG1 + dsLG2);

      dr = (rSG - rTPL);
      dz = (xTPL[2] - xSG[2]);
      p_contact->SG = -asin(dr / sqrt(dr * dr + dz * dz));

      p_contact->GA = p_contact->LG - p_contact->SG;
      p_contact->theta_sg = p_contact->SG;

      double gaWeight = 0.0;
      if (AniGaFlg == 1) {
        const double thetaFacet = fabs(p_mesh->voronkovThetaFacet);
        // calcFacetAngle() stores misorientation from the selected facet:
        // Th=0 is the exact facet and Th=thetaFacet is the rough orientation.
        if (thetaFacet > 1.e-12)
          gaWeight = 1.0 - fabs(p_contact->Th) / thetaFacet;
        if (gaWeight < 0.0)
          gaWeight = 0.0;
        if (gaWeight > 1.0)
          gaWeight = 1.0;
      }

      // gaWeight is the faceted fraction: GAf on the facet, GAr when rough.
      const double gaTarget = GArRad + gaWeight * (GAfRad - GArRad);
      p_contact->GAs = gaTarget;
      p_contact->GAv = gaTarget;
      p_contact->GABranch = VOR_GA_NO_CORRECTION;
      p_contact->V_SL_n = p_contact->Vn;
      p_contact->V_SG_n_local = 0.0;
      p_contact->V_SG_t_far = 0.0;
      p_contact->V_SG_t_local = 0.0;
      p_contact->V_vo = 0.0;
      p_contact->Vmin_SG_far = 0.0;
      p_contact->V_pull = p_mesh->Vpull;
      p_contact->trace_dot_t_SG_far = 0.0;
      const double Tstored = p_contact->p_node->attribute[TEMPERATURE];
      const double dTkinStored =
          (p_mesh->voronkovTemperatureIsAbsolute != 0)
              ? (p_mesh->voronkovTm - Tstored)
              : (-Tstored);
      // DeltaT stores the signed local thermal driving force dTkin=Tm-Tlocal.
      // Its sign, not a clamped auxiliary variable, selects the physical branch.
      p_contact->DeltaT = dTkinStored;
      p_contact->chi0_abs = 0.0;
      p_contact->chi_geom = 0.0;
      p_contact->chi_use = 0.0;
      p_contact->theta_g_0 = gaTarget;
      p_contact->theta_g_app = gaTarget;
      p_contact->phiS = 0.0;
      p_contact->phiS_mode = PHIS_INVALID;
      p_contact->alpha_TPL = fabs(p_contact->Th);
      p_contact->beta_step = p_mesh->betaStep;
      p_contact->dTloc = -dTkinStored;
      p_contact->dTsolid = (dTkinStored > 0.0) ? dTkinStored : 0.0;
      p_contact->sigmaSG_diag = p_mesh->voronkovSigmaSG;
      p_contact->sigmaSL_diag = p_mesh->voronkovSigmaSLRough;
      p_contact->sigmaLG_diag = p_mesh->voronkovSigmaMG;
      p_contact->sigmaSL_f = p_mesh->voronkovSigmaSLRough;
      p_contact->sigmaSL_alpha = 0.0;
      p_contact->SigmaSL_CH = p_mesh->voronkovSigmaSLRough;
      p_contact->psiSL_CH = 0.0;
      p_contact->deltaS_CH = 0.0;
      p_contact->theta_SG_far = p_contact->SG;
      p_contact->theta_SL_TPL = 0.0;
      p_contact->theta_SG_stop = 0.0;
      p_contact->chiV = 0.0;
      p_contact->chiMax_FB = 0.0;
      p_contact->theta_SG_TPL = p_contact->SG;

      // One update with 0.01 preserves the first-order movement rate of the
      // former 100 updates with coefficient 0.0001.
      delt = coef * .01 * (p_contact->GAv - p_contact->GA);
      if (!p_contact->p_node->artificial_boundary && rTPL > 1.e-16) {
        const double dx = -delt * xTPL[0] / rTPL;
        const double dy = -delt * xTPL[1] / rTPL;
        p_contact->p_node->coord[0] += dx;
        p_contact->p_node->coord[1] += dy;
        // The next local FE solve follows after all interface substeps. Keep
        // the S/L kinetic temperature at the displaced TPL node consistent
        // during those substeps using the resolved solid-side gradient.
        p_contact->p_node->attribute[TEMPERATURE] +=
            p_contact->p_node->dTs[0] * dx +
            p_contact->p_node->dTs[1] * dy;
        p_contact->p_node->T =
            p_contact->p_node->attribute[TEMPERATURE];
      }
  }
}

void BCSet::MoveTPLvoronkovGA(double coef, double GAr, double GAf,
                              int AniGaFlg, Mesh *p_mesh) {
  double delt, xTPL[3], xLG[3], xLG2[3], xSG[3];
  BCSet *p_sl_set = &p_mesh->bc_set[p_mesh->BC_SOLID_LIQUID];
  BCSet *p_lg_set = &p_mesh->bc_set[p_mesh->BC_LIQUID_GAS];
  if (mesh_bc_id == p_mesh->BC_SUB_TPL)
    p_sl_set = &p_mesh->bc_set[p_mesh->BC_SUB_SOLID_LIQUID];
  if (mesh_bc_id == p_mesh->BC_SUB_TPL)
    p_lg_set = &p_mesh->bc_set[p_mesh->BC_SUB_LIQUID_GAS];

  const double angleTol = 1.0e-12;
  const double velocityTol = 1.0e-20;
  const double surfaceTol = 1.0e-30;
  const double vectorTol = 1.0e-14;
  const double domainTol = 1.0e-12;
  // GAr and GAf remain the signed targets of the adjacent linear MoveTPL()
  // alternative. The Voronkov closure calculates GAeq from the stored
  // physical surface and step energies.
  (void)GAr;
  (void)GAf;

  // The TPL tangent uses the TPL node and the first two L/G nodes. Evaluate
  // the Voronkov closure once for the current state; the main-loop movement
  // substeps provide the time relaxation and retain physical TPL motion.
  for (int i(0); i < nS1; i++) {
      BCAtNode *p_contact = bc_at_nodeStr[i][1];
      BCAtNode *p_sl =
          &p_sl_set->bc_at_node[p_contact->parent_bc_node_id];
      // p_contact and p_sl can be different mesh records. The parent id is set
      // from initialized S/L geometry, so TPL uses S/L kinetics at the same
      // physical contact-line position instead of assuming equal array indexes.
      p_contact->K = p_sl->K;
      p_contact->Vn = p_sl->Vn;
      const double alphaTPL = (AniGaFlg == 1) ? fabs(p_sl->Th) : 0.;
      p_contact->Th = alphaTPL;

      for (int co = 0; co < 3; co++) {
        xTPL[co] = p_contact->p_node->coord[co];
        xLG[co] = p_contact->p_adj_lg_node->coord[co];
        xLG2[co] =
            p_lg_set
                ->bc_at_nodeStr[i +
                                (p_lg_set->structured_z_ghost_columns != 0)]
                               [2]
                ->p_node->coord[co];
        xSG[co] = p_contact->p_adj_sg_node->coord[co];
      }
      double dz(xTPL[2] - xLG[2]);
      double rTPL(sqrt(xTPL[0] * xTPL[0] + xTPL[1] * xTPL[1]));
      double rLG(sqrt(xLG[0] * xLG[0] + xLG[1] * xLG[1]));
      double rLG2(sqrt(xLG2[0] * xLG2[0] + xLG2[1] * xLG2[1]));
      double rSG(sqrt(xSG[0] * xSG[0] + xSG[1] * xSG[1]));
      double dr(rLG - rTPL);
      const double dsLG1 = sqrt(max(1.e-30, dr * dr + dz * dz));
      const double lgChord1 =
          asin(max(-1.0, min(1.0, dr / dsLG1)));
      const double drLG2 = rLG2 - rLG;
      const double dzLG2 = xLG[2] - xLG2[2];
      const double dsLG2 =
          sqrt(max(1.e-30, drLG2 * drLG2 + dzLG2 * dzLG2));
      const double lgChord2 =
          asin(max(-1.0, min(1.0, drLG2 / dsLG2)));
      // A chord gives the inclination averaged over one L/G edge.  Extrapolate
      // the first two chord inclinations linearly in arc length to obtain the
      // endpoint tangent used by the TPL growth-angle condition.
      p_contact->LG =
          lgChord1 - (lgChord2 - lgChord1) * dsLG1 / (dsLG1 + dsLG2);

      dr = (rSG - rTPL);
      dz = (xTPL[2] - xSG[2]);
      p_contact->SG =
          -asin(max(-1.0, min(1.0, dr / sqrt(max(1.e-30, dr * dr + dz * dz)))));

      p_contact->GA = p_contact->LG - p_contact->SG;
      p_contact->theta_sg = p_contact->SG;

      const double Tloc = p_contact->p_node->attribute[TEMPERATURE];
      const double dTkin =
          (p_mesh->voronkovTemperatureIsAbsolute != 0)
              ? (p_mesh->voronkovTm - Tloc)
              : (-Tloc);

      // Legacy fields:
      //   GAs  = GAeq, the S/L capillary equilibrium angle.
      //   GAv  = GAapp, the apparent angle after Voronkov rotation.
      //   V_vo = Vt_SG_far, the tangential sweep speed of the already formed
      //          far S/G surface.  It is kept only as a compatibility alias.
      //
      // GA is the L/G--S/G growth angle, GA=theta_LG_TPL-theta_SG_TPL.
      // phiS is the solid-side S/G--S/L dihedral angle.  It is not a growth
      // angle and is not used as GAeq.  It only defines theta_SG_stop and the
      // force-balance allowed rotation chiMax_FB.
      //
      // The S/L capillary selection of GAeq, the S/L curvature stiffness, and
      // the S/G Voronkov dynamic rotation are separate mechanisms.  Do not use
      // small Vn_SL_TPL to force GAeq=GA_f(0).  If the local stepped source is
      // active, the TPL can lie on a growth hill with alpha_TPL>0 even when
      // the far facet has alpha_far=0 and Vn_SL_TPL is nearly zero.
      int stepSourceActive = 0;
      // The stepped growth source is a solidification mechanism. At
      // equilibrium or during melting it is suppressed, GAeq is evaluated
      // from the rough static branch, and the Voronkov rotation remains zero.
      if (AniGaFlg == 1 && dTkin > 0.0) {
        // Use the same local 2DN/step/rough ordering as kinMoveInterface().
        // Reaching the rough upper kinetic limit means a rough S/L state;
        // otherwise the local stepped/vicinal source remains active.  This
        // state selection is separate from the Voronkov correction branch.
        const double V2dn =
            p_mesh->kineticB2DN * pow(dTkin, 5. / 6.) *
            exp(-p_mesh->kineticA2DN / dTkin);
        double Vkin = p_mesh->betaStep * alphaTPL * dTkin;
        if (Vkin < V2dn)
          Vkin = V2dn;
        if (Vkin < p_mesh->betaRough * dTkin)
          stepSourceActive = 1;
      }
      const double alphaLocal = stepSourceActive ? alphaTPL : 0.0;
      double GAeq = 0.0;
      const double dTloc = -dTkin;
      const double dTsolid = (dTkin > 0.0) ? dTkin : 0.0;
      p_contact->GAs = 0.0;
      p_contact->GAv = 0.0;
      p_contact->GABranch = VOR_GA_NO_CORRECTION;
      p_contact->V_SL_n = p_contact->Vn; // Vn_SL_TPL, diagnostic only here.
      p_contact->V_SG_n_local = 0.0;
      p_contact->V_SG_t_far = 0.0;
      p_contact->V_SG_t_local = 0.0;
      p_contact->V_vo = 0.0;
      p_contact->Vmin_SG_far = 0.0;
      p_contact->V_pull = p_mesh->Vpull;
      p_contact->trace_dot_t_SG_far = 0.0;
      p_contact->chi0_abs = 0.0;
      p_contact->chi_geom = 0.0;
      p_contact->chi_use = 0.0;
      p_contact->theta_g_0 = GAeq;
      p_contact->theta_g_app = GAeq;
      // DeltaT stores the signed local thermal driving force dTkin=Tm-Tlocal.
      p_contact->DeltaT = dTkin;
      p_contact->phiS = 0.0;
      p_contact->phiS_mode = PHIS_INVALID;
      p_contact->alpha_TPL = alphaLocal;
      p_contact->beta_step = p_mesh->betaStep;
      p_contact->dTloc = dTloc;
      p_contact->dTsolid = dTsolid;
      p_contact->sigmaSG_diag = p_mesh->voronkovSigmaSG;
      p_contact->sigmaSL_diag = p_mesh->voronkovSigmaSLRough;
      p_contact->sigmaLG_diag = p_mesh->voronkovSigmaMG;
      p_contact->sigmaSL_f = p_mesh->voronkovSigmaSLRough;
      p_contact->sigmaSL_alpha = 0.0;
      p_contact->SigmaSL_CH = p_mesh->voronkovSigmaSLRough;
      p_contact->psiSL_CH = 0.0;
      p_contact->deltaS_CH = 0.0;
      p_contact->theta_SG_far = 0.0;
      p_contact->theta_SL_TPL = 0.0;
      p_contact->theta_SG_stop = 0.0;
      p_contact->chiV = 0.0;
      p_contact->chiMax_FB = 0.0;
      p_contact->theta_SG_TPL = 0.0;

      if (stepSourceActive) {
        // voronkovSigmaSLPrimeFacet stores |tauStep/hStep|.  The selected
        // local facet direction is toward negative e2, so its oriented
        // derivative is negative in the same basis used for theta_SL_TPL.
        const double sigmaSLPrimeFacet =
            -fabs(p_mesh->voronkovSigmaSLPrimeFacet);
        p_contact->sigmaSL_f =
            p_mesh->voronkovSigmaSLFacet * cos(alphaLocal) +
            sigmaSLPrimeFacet * sin(alphaLocal);
        p_contact->sigmaSL_alpha =
            -p_mesh->voronkovSigmaSLFacet * sin(alphaLocal) +
            sigmaSLPrimeFacet * cos(alphaLocal);
        p_contact->SigmaSL_CH =
            sqrt(p_contact->sigmaSL_f * p_contact->sigmaSL_f +
                 p_contact->sigmaSL_alpha * p_contact->sigmaSL_alpha);
        p_contact->psiSL_CH =
            atan2(p_contact->sigmaSL_alpha, p_contact->sigmaSL_f);
        p_contact->sigmaSL_diag = p_contact->sigmaSL_f;

        // Faceted/vicinal solid-side S/G--S/L dihedral angle.  The S/L
        // side is represented by its Cahn-Hoffman vector:
        //   SigmaSL_CH = sqrt(sigmaSL^2 + sigmaSL_alpha^2)
        //   psiSL_CH   = atan2(sigmaSL_alpha, sigmaSL)
        // deltaS_CH is the isotropic triangle angle between sigmaSG and
        // SigmaSL_CH, opposite sigmaMG.  In this local sign convention the
        // actual solid-side angle from the S/L tangent to S/G is
        // phiS = psiSL_CH + deltaS_CH.
        const double denomDeltaS =
            2.0 * p_mesh->voronkovSigmaSG * p_contact->SigmaSL_CH;
        if (denomDeltaS > surfaceTol) {
          double cosDeltaS =
              (p_mesh->voronkovSigmaSG * p_mesh->voronkovSigmaSG +
               p_contact->SigmaSL_CH * p_contact->SigmaSL_CH -
               p_mesh->voronkovSigmaMG * p_mesh->voronkovSigmaMG) /
              denomDeltaS;
          if (cosDeltaS > 1.0 + domainTol ||
              cosDeltaS < -1.0 - domainTol) {
            p_contact->GABranch = VOR_GA_INVALID_GEOMETRY;
            p_contact->phiS_mode = PHIS_INVALID;
          } else {
            if (cosDeltaS > 1.0) cosDeltaS = 1.0;
            if (cosDeltaS < -1.0) cosDeltaS = -1.0;
            p_contact->deltaS_CH = acos(cosDeltaS);
            p_contact->phiS = p_contact->psiSL_CH + p_contact->deltaS_CH;
            p_contact->phiS_mode = PHIS_VICINAL_CH;
          }
        } else {
          p_contact->GABranch = VOR_GA_INVALID_GEOMETRY;
          p_contact->phiS_mode = PHIS_INVALID;
        }
      } else {
        // Rough isotropic solid-side S/G--S/L dihedral angle:
        // cos(phiS) =
        //   (sigmaSG^2 + sigmaSL^2 - sigmaMG^2)/(2*sigmaSG*sigmaSL).
        // This is not GAeq.  It only sets the S/G stop tangent.
        const double denomPhiS =
            2.0 * p_mesh->voronkovSigmaSG * p_mesh->voronkovSigmaSLRough;
        if (denomPhiS > surfaceTol) {
          double cosPhiS =
              (p_mesh->voronkovSigmaSG * p_mesh->voronkovSigmaSG +
               p_mesh->voronkovSigmaSLRough *
                   p_mesh->voronkovSigmaSLRough -
               p_mesh->voronkovSigmaMG * p_mesh->voronkovSigmaMG) /
              denomPhiS;
          if (cosPhiS > 1.0 + domainTol || cosPhiS < -1.0 - domainTol) {
            p_contact->GABranch = VOR_GA_INVALID_GEOMETRY;
            p_contact->phiS_mode = PHIS_INVALID;
          } else {
            if (cosPhiS > 1.0) cosPhiS = 1.0;
            if (cosPhiS < -1.0) cosPhiS = -1.0;
            p_contact->phiS = acos(cosPhiS);
            p_contact->phiS_mode = PHIS_ROUGH_ISOTROPIC;
          }
        } else {
          p_contact->GABranch = VOR_GA_INVALID_GEOMETRY;
          p_contact->phiS_mode = PHIS_INVALID;
        }
      }

      // Voronkov free-energy-vector balance for the static equilibrium
      // growth angle. All surface quantities here use the same cgs units.
      // The stepped branch includes the physical step energy tauStep/hStep.
      if (p_contact->GABranch != VOR_GA_INVALID_GEOMETRY) {
        const double sigmaDiff =
            p_mesh->voronkovSigmaSG - p_mesh->voronkovSigmaMG;
        const double radicand =
            p_contact->sigmaSL_f * p_contact->sigmaSL_f +
            p_contact->sigmaSL_alpha * p_contact->sigmaSL_alpha -
            sigmaDiff * sigmaDiff;
        const double denominator =
            sqrt(p_mesh->voronkovSigmaSG * p_mesh->voronkovSigmaMG);

        if (radicand < -domainTol || denominator <= surfaceTol) {
          p_contact->GABranch = VOR_GA_INVALID_GEOMETRY;
        } else {
          double asinArgument =
              0.5 * sqrt(max(0.0, radicand)) / denominator;
          if (asinArgument > 1.0 + domainTol) {
            p_contact->GABranch = VOR_GA_INVALID_GEOMETRY;
          } else {
            if (asinArgument > 1.0) asinArgument = 1.0;
            GAeq = 2.0 * asin(asinArgument);
          }
        }
      }

      p_contact->GAs = GAeq;
      p_contact->GAv = GAeq;
      p_contact->theta_g_0 = GAeq;
      p_contact->theta_g_app = GAeq;

      {
        int vt_ok = 0;
        double Vt_SG_far = 0.0;
        double traceDotSGFar = 0.0;
        int sgFormationReversed = 0;

        int localBasisOk = 1;
        double eTPL[3] = {0.0, 0.0, 0.0};
        if (nS1 > 2) {
          int i_prev = (i + nS1 - 1) % nS1;
          int i_next = (i + 1) % nS1;
          if (p_lg_set->structured_z_ghost_columns != 0 ||
              p_lg_set->structured_fixed_side_columns != 0) {
            i_prev = (i > 0) ? i - 1 : i;
            i_next = (i + 1 < nS1) ? i + 1 : i;
          }
          Node *p_prev = bc_at_nodeStr[i_prev][1]->p_node;
          Node *p_next = bc_at_nodeStr[i_next][1]->p_node;
          for (int co = 0; co < 3; co++)
            eTPL[co] = p_next->coord[co] - p_prev->coord[co];
        } else if (rTPL > vectorTol) {
          eTPL[0] = -xTPL[1] / rTPL;
          eTPL[1] = xTPL[0] / rTPL;
          eTPL[2] = 0.0;
        }
        double nTPLv =
            sqrt(eTPL[0] * eTPL[0] + eTPL[1] * eTPL[1] +
                 eTPL[2] * eTPL[2]);
        if (nTPLv <= vectorTol)
          localBasisOk = 0;
        else
          for (int co = 0; co < 3; co++)
            eTPL[co] /= nTPLv;

        double tSG_far[3] = {0.0, 0.0, 0.0};
        for (int co = 0; co < 3; co++)
          tSG_far[co] = xSG[co] - xTPL[co];
        double nSGfar =
            sqrt(tSG_far[0] * tSG_far[0] + tSG_far[1] * tSG_far[1] +
                 tSG_far[2] * tSG_far[2]);
        if (nSGfar <= vectorTol)
          localBasisOk = 0;
        else
          for (int co = 0; co < 3; co++)
            tSG_far[co] /= nSGfar;

        double e1[3] = {0.0, 0.0, 0.0};
        double e2[3] = {0.0, 0.0, 0.0};
        if (localBasisOk) {
          const double dotSG =
              tSG_far[0] * eTPL[0] + tSG_far[1] * eTPL[1] +
              tSG_far[2] * eTPL[2];
          for (int co = 0; co < 3; co++)
            e1[co] = tSG_far[co] - dotSG * eTPL[co];
          double n1 =
              sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
          if (n1 <= vectorTol)
            localBasisOk = 0;
          else {
            for (int co = 0; co < 3; co++)
              e1[co] /= n1;
            e2[0] = eTPL[1] * e1[2] - eTPL[2] * e1[1];
            e2[1] = eTPL[2] * e1[0] - eTPL[0] * e1[2];
            e2[2] = eTPL[0] * e1[1] - eTPL[1] * e1[0];
          }
        }

        if (localBasisOk) {
          // Present quasi-steady Cz trace construction:
          //   P_l = I - t_TPL (x) t_TPL,
          //   tbar_SG_far = P_l t_SG_far / |P_l t_SG_far| = e1,
          //   v_trace = Vpull e_z,
          //   Vt_SG_far = |(P_l v_trace) . tbar_SG_far|.
          // Because e1 is normal to t_TPL, the final dot product is also
          // v_trace.e1.  This gives |Vpull| for a vertical S/G surface at a
          // symmetric ridge tip and includes the local TPL inclination on a
          // ridge shoulder.  Vn_SL_TPL is not used.  The inverse relation
          // V=Vn/sin(chi) applies only if Vn is the newly formed S/G normal
          // rate, Vn_SG_TPL, and is diagnostic rather than an input here.
          const double vTrace[3] = {0.0, 0.0, p_mesh->Vpull};
          const double dotTraceTPL =
              vTrace[0] * eTPL[0] + vTrace[1] * eTPL[1] +
              vTrace[2] * eTPL[2];
          double vTracePlane[3];
          for (int co = 0; co < 3; co++)
            vTracePlane[co] = vTrace[co] - dotTraceTPL * eTPL[co];
          traceDotSGFar =
              vTracePlane[0] * e1[0] + vTracePlane[1] * e1[1] +
              vTracePlane[2] * e1[2];
          Vt_SG_far = fabs(traceDotSGFar);
          if (traceDotSGFar < -velocityTol)
            sgFormationReversed = 1;
          else if (traceDotSGFar > velocityTol)
            vt_ok = 1;
        }

        p_contact->V_vo = Vt_SG_far;
        p_contact->V_SG_t_far = Vt_SG_far;
        p_contact->trace_dot_t_SG_far = traceDotSGFar;
        p_contact->V_SG_t_local = Vt_SG_far;

        double tSL_TPL[3] = {0.0, 0.0, 0.0};
        double bestInward = -1.e100;
        double bestLength = 0.0;
        const double erx = (rTPL > vectorTol) ? xTPL[0] / rTPL : 1.0;
        const double ery = (rTPL > vectorTol) ? xTPL[1] / rTPL : 0.0;
        for (int n_adj(0); n_adj < p_sl->num_of_adj_nodes; n_adj++) {
          BCAtNode *p_adj = &p_sl_set->bc_at_node[p_sl->n_adj_node[n_adj]];
          double v[3];
          for (int co = 0; co < 3; co++)
            v[co] = p_adj->p_node->coord[co] - p_sl->p_node->coord[co];
          const double outward = v[0] * erx + v[1] * ery;
          const double inward = -outward;
          const double len =
              sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
          if (len > vectorTol && inward > bestInward) {
            bestInward = inward;
            bestLength = len;
            for (int co = 0; co < 3; co++)
              tSL_TPL[co] = v[co] / len;
          }
        }
        if (bestLength <= vectorTol)
          localBasisOk = 0;

        double thetaSGFar = 0.0;
        double thetaSLTPL = 0.0;
        if (localBasisOk) {
          const double dotSL =
              tSL_TPL[0] * eTPL[0] + tSL_TPL[1] * eTPL[1] +
              tSL_TPL[2] * eTPL[2];
          double tSLp[3];
          for (int co = 0; co < 3; co++)
            tSLp[co] = tSL_TPL[co] - dotSL * eTPL[co];
          double nSLp =
              sqrt(tSLp[0] * tSLp[0] + tSLp[1] * tSLp[1] +
                   tSLp[2] * tSLp[2]);
          if (nSLp <= vectorTol)
            localBasisOk = 0;
          else {
            for (int co = 0; co < 3; co++)
              tSLp[co] /= nSLp;
            const double sl2 =
                tSLp[0] * e2[0] + tSLp[1] * e2[1] + tSLp[2] * e2[2];
            const double sl1 =
                tSLp[0] * e1[0] + tSLp[1] * e1[1] + tSLp[2] * e1[2];
            thetaSLTPL = atan2(sl2, sl1);
            while (thetaSLTPL > 0.5 * PI) thetaSLTPL -= PI;
            while (thetaSLTPL < -0.5 * PI) thetaSLTPL += PI;
            if (thetaSLTPL > 0.0)
              thetaSLTPL = -thetaSLTPL;
          }
        }

        if (localBasisOk && stepSourceActive) {
          // e1 is the far S/G direction and negative e2 points toward the
          // S/L side.  The two sides of a stepped ridge have opposite local
          // constructions about the exact facet:
          //   symmetric-tip side: theta_SL_TPL = -theta_f + alpha_TPL,
          //   shoulder side:      theta_SL_TPL = -theta_f - alpha_TPL.
          // Select the candidate closest to the S/L tangent resolved from the
          // actual mesh vectors above.  An exact tie selects the tip branch;
          // this is also the branch required when a future one-sided subgrid
          // alpha is supplied at an exactly symmetric facet point.
          const double thetaFacet =
              fabs(p_mesh->voronkovThetaFacet);
          const double thetaSLResolved = thetaSLTPL;
          const double thetaSLTip = -thetaFacet + alphaLocal;
          const double thetaSLShoulder = -thetaFacet - alphaLocal;
          if (fabs(thetaSLResolved - thetaSLTip) <=
              fabs(thetaSLResolved - thetaSLShoulder))
            thetaSLTPL = thetaSLTip;
          else
            thetaSLTPL = thetaSLShoulder;
        }

        p_contact->theta_SG_far = thetaSGFar;
        p_contact->theta_SL_TPL = thetaSLTPL;
        p_contact->theta_SG_stop = thetaSLTPL + p_contact->phiS;
        p_contact->chiMax_FB =
            thetaSGFar - p_contact->theta_SG_stop;
        p_contact->chi_geom = p_contact->chiMax_FB;

        double chiUse = 0.0;
        if (!localBasisOk || p_contact->phiS_mode == PHIS_INVALID ||
            p_contact->GABranch == VOR_GA_INVALID_GEOMETRY) {
          p_contact->GABranch = VOR_GA_INVALID_GEOMETRY;
          chiUse = 0.0;
        } else if (dTkin <= 0.0) {
          // The classical Voronkov equation is a solidification closure.
          // At equilibrium or during melting VtKin may be zero or negative,
          // but no dynamic S/G rotation is applied: chi_use=0 and GAapp=GAeq.
          chiUse = 0.0;
          p_contact->GABranch = VOR_GA_NO_CORRECTION;
        } else if (sgFormationReversed) {
          chiUse = 0.0;
          p_contact->GABranch = VOR_GA_NO_CORRECTION;
        } else if (p_contact->chiMax_FB <= angleTol) {
          chiUse = 0.0;
          p_contact->GABranch = VOR_GA_NO_CORRECTION;
        } else {
          const double A =
              (p_mesh->voronkovQ / p_mesh->voronkovTm) *
              pow(p_mesh->voronkovLambdaSG, 1.0 / 3.0) *
              pow(p_mesh->voronkovSigmaSG *
                  p_mesh->voronkovAtomicDensity, -2.0 / 3.0) *
              dTkin;
          const double Vmin_SG_far = pow(A / p_contact->chiMax_FB, 3.0);
          p_contact->Vmin_SG_far = Vmin_SG_far;

          if (!vt_ok) {
            if (stepSourceActive) {
              chiUse = p_contact->chiMax_FB;
              p_contact->GABranch = VOR_GA_FORCE_BALANCE_STOP;
            } else
              p_contact->GABranch = VOR_GA_INVALID_GEOMETRY;
          } else {
            const double Vt_SG_far_voronkov =
                Vt_SG_far * p_mesh->voronkovVelocityToCm;
            const double chiRaw =
                A * pow(Vt_SG_far_voronkov, -1.0 / 3.0);
            p_contact->chi0_abs = chiRaw;
            p_contact->chiV = chiRaw;

            // Vmin_SG_far and Vt_SG_far_voronkov are both in cm/s here.
            // The threshold branch is written explicitly so the stored status
            // records the physical force-balance stop.
            if (Vt_SG_far_voronkov < Vmin_SG_far) {
              chiUse = p_contact->chiMax_FB;
              p_contact->GABranch = VOR_GA_FORCE_BALANCE_STOP;
            } else if (chiRaw >= 0.0) {
              chiUse = chiRaw;
              p_contact->GABranch = VOR_GA_VELOCITY_CONTROLLED;
            } else {
              chiUse = 0.0;
              p_contact->GABranch = VOR_GA_INVALID_GEOMETRY;
            }
          }
        }

        if (localBasisOk) {
          const double thetaSGTPLLocal = thetaSGFar - chiUse;
          double tSG_TPL[3];
          for (int co = 0; co < 3; co++)
            tSG_TPL[co] = cos(thetaSGTPLLocal) * e1[co] +
                           sin(thetaSGTPLLocal) * e2[co];
          const double dotTPL =
              tSG_TPL[0] * eTPL[0] + tSG_TPL[1] * eTPL[1] +
              tSG_TPL[2] * eTPL[2];
          const double normTPL =
              sqrt(tSG_TPL[0] * tSG_TPL[0] +
                   tSG_TPL[1] * tSG_TPL[1] +
                   tSG_TPL[2] * tSG_TPL[2]);
          if (fabs(dotTPL) > vectorTol || fabs(normTPL - 1.0) > vectorTol) {
            chiUse = 0.0;
            p_contact->GABranch = VOR_GA_INVALID_GEOMETRY;
          }

          if (p_contact->GABranch == VOR_GA_FORCE_BALANCE_STOP) {
            double tSG_stop[3];
            for (int co = 0; co < 3; co++)
              tSG_stop[co] = cos(p_contact->theta_SG_stop) * e1[co] +
                              sin(p_contact->theta_SG_stop) * e2[co];
            double dotStop =
                tSG_TPL[0] * tSG_stop[0] +
                tSG_TPL[1] * tSG_stop[1] +
                tSG_TPL[2] * tSG_stop[2];
            if (dotStop > 1.0) dotStop = 1.0;
            if (dotStop < -1.0) dotStop = -1.0;
            const double stopError = acos(dotStop);
            if (stopError > angleTol) {
              chiUse = 0.0;
              p_contact->GABranch = VOR_GA_INVALID_GEOMETRY;
            }
          }
        }

        p_contact->chi_use = chiUse;
        p_contact->theta_SG_TPL = thetaSGFar - chiUse;
        p_contact->GAv =
            (p_contact->GABranch == VOR_GA_INVALID_GEOMETRY)
                ? p_contact->GA
                : GAeq - chiUse;
        p_contact->theta_g_app = p_contact->GAv;
        p_contact->V_SG_n_local = Vt_SG_far * sin(-chiUse);
        p_contact->V_SG_t_local = Vt_SG_far * cos(-chiUse);
      }

      // One update with 0.01 preserves the first-order movement rate of the
      // former 100 updates with coefficient 0.0001.
      delt = coef * .01 * (p_contact->GAv - p_contact->GA);
      if (!p_contact->p_node->artificial_boundary && rTPL > 1.e-16) {
        const double dx = -delt * xTPL[0] / rTPL;
        const double dy = -delt * xTPL[1] / rTPL;
        p_contact->p_node->coord[0] += dx;
        p_contact->p_node->coord[1] += dy;
        // The next local FE solve follows after all interface substeps. Keep
        // the S/L kinetic temperature at the displaced TPL node consistent
        // during those substeps using the resolved solid-side gradient.
        p_contact->p_node->attribute[TEMPERATURE] +=
            p_contact->p_node->dTs[0] * dx +
            p_contact->p_node->dTs[1] * dy;
        p_contact->p_node->T =
            p_contact->p_node->attribute[TEMPERATURE];
      }
  }
}

void Mesh::calc3dmeshDef() {
  for (int ib(0); ib < number_of_nodes; ib++)
    for (int co = 0; co < 3; co++)
      node[ib].delta[co] = 0.;

  for (int ib(0); ib < number_of_elements; ib++) {
    for (int i(0); i < 8; i++) {
      for (int j(0); j < 8; j++) {
        float coef(1.);
        if (element[ib].p_node[j]->bound_mark != -1 &&
            element[ib].p_node[j]->bound_mark != Mesh::BC_COLD)
          coef = 2.;
        for (int co = 0; co < 3; co++)
          element[ib].p_node[i]->delta[co] +=
              .01 * coef *
              (element[ib].p_node[j]->coord[co] -
               element[ib].p_node[i]->coord[co]);
      }
    }
  }
  for (int ib(0); ib < number_of_nodes; ib++)
    if (node[ib].bound_mark == -1) {

      node[ib].coord[2] += node[ib].delta[2];
      node[ib].coord[1] += node[ib].delta[1];
      node[ib].coord[0] += node[ib].delta[0];
    }
}
