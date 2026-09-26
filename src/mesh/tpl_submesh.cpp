// Nested structured TPL submesh coupling.
//
// This file is included by the two active 3-D drivers after mesh.h,
// basis_func.h, and local_calc.h.  The local mesh overlaps the coarse mesh but
// is never assembled into the coarse FE system.

static int tplInnerNodeId(const Mesh *p_local, int i, int j, int q) {
  const int nPhi = p_local->nTPL;
  const int ni = p_local->subSlRefinement;
  return i + nPhi * (j + (ni + 1) * q);
}

static int tplOuterNodeId(const Mesh *p_local, int i, int j, int k) {
  const int nPhi = p_local->nTPL;
  const int ns = p_local->tplSubmeshSupportLayers;
  if (j == 0)
    return tplInnerNodeId(p_local, i, 0, ns - k);
  const int innerCount =
      nPhi * (p_local->subSlRefinement + 1) * (2 * ns + 1);
  const int nr = p_local->tplSubmeshOuterElements;
  return innerCount + i + nPhi * ((j - 1) + nr * k);
}

static void tplBrickShape(double eta, double nu, double ksi,
                          double *N, double dN[8][3]) {
  const double se[8] = {-1., -1., -1., -1., 1., 1., 1., 1.};
  const double sn[8] = {-1., 1., 1., -1., -1., 1., 1., -1.};
  const double sk[8] = {-1., -1., 1., 1., -1., -1., 1., 1.};
  for (int i(0); i < 8; i++) {
    N[i] = .125 * (1. + se[i] * eta) *
           (1. + sn[i] * nu) * (1. + sk[i] * ksi);
    dN[i][0] = .125 * se[i] * (1. + sn[i] * nu) *
               (1. + sk[i] * ksi);
    dN[i][1] = .125 * sn[i] * (1. + se[i] * eta) *
               (1. + sk[i] * ksi);
    dN[i][2] = .125 * sk[i] * (1. + se[i] * eta) *
               (1. + sn[i] * nu);
  }
}

static void tplInterpolateElementPoint(const Element *p_element,
                                       const double nat[3],
                                       double point[3]) {
  double N[8];
  double dN[8][3];
  tplBrickShape(nat[0], nat[1], nat[2], N, dN);
  point[0] = 0.;
  point[1] = 0.;
  point[2] = 0.;
  for (int i(0); i < 8; i++)
    for (int co(0); co < 3; co++)
      point[co] += N[i] * p_element->p_node[i]->coord[co];
}

// Second-order reconstruction from the structured parent-brick nodes. The
// one-half gradient correction preserves nodal temperatures, is exact for a
// linear field, and reproduces a quadratic field when the recovered nodal
// gradients are exact.
static double tplInterpolateElementTemperatureSecondOrder(
    const Element *p_element, const double nat[3], double timeFraction,
    const double *p_oldSolidGradient,
    const double *p_oldLiquidGradient) {
  double N[8];
  double dN[8][3];
  tplBrickShape(nat[0], nat[1], nat[2], N, dN);
  double point[3] = {0., 0., 0.};
  for (int i(0); i < 8; i++)
    for (int co(0); co < 3; co++)
      point[co] += N[i] * p_element->p_node[i]->coord[co];

  const int solid = (p_element->group_num == Mesh::SOLID);
  double temperature = 0.;
  for (int i(0); i < 8; i++) {
    Node *p_node = p_element->p_node[i];
    const double nodalTemperature =
        (1. - timeFraction) * p_node->T +
        timeFraction * p_node->attribute[TEMPERATURE];
    double gradientCorrection = 0.;
    for (int co(0); co < 3; co++) {
      const double newGradient = solid ? p_node->dTs[co] : p_node->dTl[co];
      const double *p_oldGradient =
          solid ? p_oldSolidGradient : p_oldLiquidGradient;
      const double oldGradient =
          (p_oldGradient == NULL)
              ? newGradient
              : p_oldGradient[3 * p_node->node_num + co];
      const double gradient =
          (1. - timeFraction) * oldGradient +
          timeFraction * newGradient;
      gradientCorrection +=
          gradient * (point[co] - p_node->coord[co]);
    }
    temperature +=
        N[i] * (nodalTemperature + .5 * gradientCorrection);
  }
  return temperature;
}

static int tplInverseBrick(const Element *p_element, const double point[3],
                           double nat[3], const double *p_start) {
  static const double start[9][3] = {
      {0., 0., 0.},
      {-.9, -.9, -.9}, {-.9, .9, -.9},
      {-.9, .9, .9}, {-.9, -.9, .9},
      {.9, -.9, -.9}, {.9, .9, -.9},
      {.9, .9, .9}, {.9, -.9, .9}};
  double scale = 0.;
  for (int i(0); i < 8; i++)
    for (int co(0); co < 3; co++)
      scale = max(scale, fabs(p_element->p_node[i]->coord[co]));
  const double errorTol = 1.e-9 * max(1., scale);
  for (int is(-1); is < 9; is++) {
    double trial[3];
    if (is < 0) {
      if (p_start == NULL)
        continue;
      trial[0] = p_start[0];
      trial[1] = p_start[1];
      trial[2] = p_start[2];
    } else {
      trial[0] = start[is][0];
      trial[1] = start[is][1];
      trial[2] = start[is][2];
    }
    for (int it(0); it < 30; it++) {
      double N[8];
      double dN[8][3];
      tplBrickShape(trial[0], trial[1], trial[2], N, dN);
      double image[3] = {0., 0., 0.};
      double jac[3][3] = {};
      for (int i(0); i < 8; i++)
        for (int co(0); co < 3; co++) {
          image[co] += N[i] * p_element->p_node[i]->coord[co];
          for (int iloc(0); iloc < 3; iloc++)
            jac[co][iloc] +=
                dN[i][iloc] * p_element->p_node[i]->coord[co];
        }
      double residual[3];
      for (int co(0); co < 3; co++)
        residual[co] = image[co] - point[co];
      const double det =
          jac[0][0] * (jac[1][1] * jac[2][2] -
                       jac[1][2] * jac[2][1]) -
          jac[0][1] * (jac[1][0] * jac[2][2] -
                       jac[1][2] * jac[2][0]) +
          jac[0][2] * (jac[1][0] * jac[2][1] -
                       jac[1][1] * jac[2][0]);
      if (!isfinite(det) || fabs(det) < 1.e-30)
        break;
      double delta[3];
      delta[0] =
          ((jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) *
               residual[0] +
           (jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]) *
               residual[1] +
           (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]) *
               residual[2]) /
          det;
      delta[1] =
          ((jac[1][2] * jac[2][0] - jac[1][0] * jac[2][2]) *
               residual[0] +
           (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]) *
               residual[1] +
           (jac[0][2] * jac[1][0] - jac[0][0] * jac[1][2]) *
               residual[2]) /
          det;
      delta[2] =
          ((jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]) *
               residual[0] +
           (jac[0][1] * jac[2][0] - jac[0][0] * jac[2][1]) *
               residual[1] +
           (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]) *
               residual[2]) /
          det;
      double deltaNorm = sqrt(delta[0] * delta[0] +
                              delta[1] * delta[1] +
                              delta[2] * delta[2]);
      const double damping = (deltaNorm > 2.) ? 2. / deltaNorm : 1.;
      for (int co(0); co < 3; co++)
        trial[co] -= damping * delta[co];
      if (deltaNorm * damping < 1.e-11)
        break;
      if (fabs(trial[0]) > 10. || fabs(trial[1]) > 10. ||
          fabs(trial[2]) > 10.)
        break;
    }
    double image[3];
    tplInterpolateElementPoint(p_element, trial, image);
    const double error =
        sqrt((image[0] - point[0]) * (image[0] - point[0]) +
             (image[1] - point[1]) * (image[1] - point[1]) +
             (image[2] - point[2]) * (image[2] - point[2]));
    const double excess =
        max(fabs(trial[0]), max(fabs(trial[1]), fabs(trial[2])));
    if (error <= errorTol && excess <= 4.) {
      for (int co(0); co < 3; co++)
        nat[co] = trial[co];
      return 1;
    }
  }
  return 0;
}

static int tplFindCoarseParent(const Mesh *p_coarse, const double point[3],
                               int preferredGroup, int previousElement,
                               const double *p_start, int *p_element_id,
                               double nat[3]) {
  // Liquid nodes stay in liquid; solid nodes stay in solid. Do not inverse-
  // map a constructed melt-side cut into gas through the TPL wedge.
  if (previousElement >= 0 && previousElement < p_coarse->number_of_elements &&
      (preferredGroup < 0 ||
       p_coarse->element[previousElement].group_num == preferredGroup) &&
      tplInverseBrick(&p_coarse->element[previousElement], point, nat,
                      p_start)) {
    *p_element_id = previousElement;
    return 1;
  }

  int nearest = -1;
  double nearestDist = 1.e300;
  for (int ie(0); ie < p_coarse->number_of_elements; ie++) {
    const Element *p_element = &p_coarse->element[ie];
    if (p_element->number_of_nodes != 8)
      continue;
    if (preferredGroup >= 0 && p_element->group_num != preferredGroup)
      continue;

    double xmin[3] = {1.e300, 1.e300, 1.e300};
    double xmax[3] = {-1.e300, -1.e300, -1.e300};
    for (int i(0); i < 8; i++)
      for (int co(0); co < 3; co++) {
        xmin[co] = min(xmin[co], p_element->p_node[i]->coord[co]);
        xmax[co] = max(xmax[co], p_element->p_node[i]->coord[co]);
      }
    double dist = 0.;
    int inBox = 1;
    for (int co(0); co < 3; co++) {
      const double tol = 1.e-9 + 1.e-7 * (xmax[co] - xmin[co]);
      if (point[co] < xmin[co] - tol || point[co] > xmax[co] + tol)
        inBox = 0;
      if (point[co] < xmin[co])
        dist += (xmin[co] - point[co]) * (xmin[co] - point[co]);
      else if (point[co] > xmax[co])
        dist += (point[co] - xmax[co]) * (point[co] - xmax[co]);
    }
    if (dist < nearestDist) {
      nearestDist = dist;
      nearest = ie;
    }
    if (!inBox)
      continue;
    if (tplInverseBrick(p_element, point, nat, NULL)) {
      *p_element_id = ie;
      return 1;
    }
  }
  if (nearest < 0)
    return 0;
  if (tplInverseBrick(&p_coarse->element[nearest], point, nat, NULL)) {
    *p_element_id = nearest;
    return 1;
  }
  return 0;
}

static int tplNodePreferredGroup(const Node *p_node) {
  if (p_node->local_block == 1)
    return Mesh::SOLID;
  if (p_node->local_block == 2 || p_node->local_block == 3)
    return Mesh::LIQUID;
  return -1;
}

static int tplParentIsTplSubmesh(const Mesh *p_parent) {
  return p_parent->bc_set != NULL &&
         p_parent->number_of_bc_sets > p_parent->BC_SUB_TPL &&
         p_parent->bc_set[p_parent->BC_TPL].bc_at_nodeStr == NULL &&
         p_parent->bc_set[p_parent->BC_SUB_TPL].bc_at_nodeStr != NULL &&
         p_parent->bc_set[p_parent->BC_SUB_TPL].nS1 == p_parent->nTPL &&
         p_parent->bc_set[p_parent->BC_SUB_TPL].nS2 == 3 &&
         p_parent->nTPL >= 3 &&
         p_parent->subSlRefinement >= 2 &&
         p_parent->tplSubmeshSupportLayers >= 2;
}

// Recursive children are defined on the parent grid. Liquid lives only
// below L/G and below S/L; solid lives only on the solid side of S/G.
// Assign the parent brick from structured indices so cut nodes are not
// inverse-mapped into the empty gas.
static int tplAssignStructuredParent(Mesh *p_local, const Mesh *p_parent,
                                     Node *p_node) {
  const int nPhi = p_parent->nTPL;
  const int ni = p_parent->subSlRefinement;
  const int ns = p_parent->tplSubmeshSupportLayers;
  const int nr = p_parent->tplSubmeshOuterElements;
  const int childNs = p_local->tplSubmeshSupportLayers;
  const int phiRefinement = p_local->nTPL / nPhi;
  const double phi =
      (double)p_node->local_i / (double)phiRefinement;
  int i = (int)floor(phi);
  if (i >= nPhi)
    i = nPhi - 1;
  const double fracI = phi - (double)i;
  const double ksi = -1. + 2. * fracI;
  const int nSolid = ns * ni * nPhi;
  int element = -1;
  double eta = 0.;
  double nu = 0.;
  if (p_node->local_block != 3) {
    const double path = .5 * (double)p_node->local_j;
    int j = (int)floor(path);
    if (j >= ni)
      j = ni - 1;
    const double fracJ = path - (double)j;
    nu = 1. - 2. * fracJ;
    const int q = p_node->local_k + childNs;
    const int childDistance = abs(q - childNs);
    const double parentDistance = .5 * (double)childDistance;
    if (p_node->local_block == 1) {
      const double solidPath = parentDistance;
      int k = (int)floor(solidPath);
      if (k >= ns)
        k = ns - 1;
      const double fracK = solidPath - (double)k;
      eta = -1. + 2. * fracK;
      element = k * (ni * nPhi) + j * nPhi + i;
    } else {
      const double liquidQ = (double)ns - parentDistance;
      int q0 = (int)floor(liquidQ);
      if (q0 >= ns)
        q0 = ns - 1;
      const double fracK = liquidQ - (double)q0;
      const int k = ns - 1 - q0;
      eta = -1. + 2. * fracK;
      element = nSolid + k * (ni * nPhi) + j * nPhi + i;
    }
  } else {
    const int jNode = p_node->local_j;
    const int kNode = p_node->local_k;
    const double path = .5 * (double)jNode;
    const double depth = .5 * (double)kNode;
    int j = (int)floor(path);
    if (j >= nr)
      j = nr - 1;
    const double fracJ = path - (double)j;
    nu = -1. + 2. * fracJ;
    int k = (int)floor(depth);
    if (k >= ns)
      k = ns - 1;
    const double fracK = depth - (double)k;
    eta = 1. - 2. * fracK;
    element = nSolid + nSolid + k * (nr * nPhi) + j * nPhi + i;
  }
  if (element < 0 || element >= p_parent->number_of_elements)
    return 0;
  p_node->parent_element_id = element;
  p_node->parent_nat[0] = eta;
  p_node->parent_nat[1] = nu;
  p_node->parent_nat[2] = ksi;
  p_node->parent_mapping_initialized = 1;
  return 1;
}

static void tplMapLocalNodes(Mesh *p_local, const Mesh *p_coarse) {
  int failed = 0;
  for (int in(0); in < p_local->number_of_nodes; in++) {
    Node *p_node = &p_local->node[in];
    // A physical child node is mapped once when the submesh is created.  Its
    // temperature then evolves from the local FE equation; only artificial
    // Dirichlet-cut nodes need a new parent interpolation after deformation.
    if (!p_node->artificial_boundary &&
        p_node->parent_mapping_initialized)
      continue;
    if (tplParentIsTplSubmesh(p_coarse)) {
      if (!tplAssignStructuredParent(p_local, p_coarse, p_node)) {
        if (failed == 0)
          cerr << "TPL parent miss node=" << in
               << " i=" << p_node->local_i
               << " j=" << p_node->local_j
               << " k=" << p_node->local_k
               << " block=" << p_node->local_block
               << " xyz=" << p_node->coord[0] << ","
               << p_node->coord[1] << "," << p_node->coord[2] << "\n";
        failed++;
      }
      continue;
    }
    if (p_node->parent_element_id >= 0 &&
        p_node->parent_element_id < p_coarse->number_of_elements) {
      double image[3];
      tplInterpolateElementPoint(
          &p_coarse->element[p_node->parent_element_id],
          p_node->parent_nat, image);
      if (isfinite(image[0]) && isfinite(image[1]) && isfinite(image[2])) {
        const double error =
            sqrt((image[0] - p_node->coord[0]) *
                     (image[0] - p_node->coord[0]) +
                 (image[1] - p_node->coord[1]) *
                     (image[1] - p_node->coord[1]) +
                 (image[2] - p_node->coord[2]) *
                     (image[2] - p_node->coord[2]));
        if (error <= 1.e-10) {
          p_node->parent_mapping_initialized = 1;
          continue;
        }
      }
    }
    int parent = -1;
    double nat[3];
    const double *p_start =
        (p_node->parent_element_id >= 0) ? p_node->parent_nat : NULL;
    if (!tplFindCoarseParent(p_coarse, p_node->coord,
                             tplNodePreferredGroup(p_node),
                             p_node->parent_element_id, p_start, &parent,
                             nat)) {
      if (p_node->parent_mapping_initialized &&
          p_node->parent_element_id >= 0 &&
          p_node->parent_element_id < p_coarse->number_of_elements) {
        double image[3];
        tplInterpolateElementPoint(
            &p_coarse->element[p_node->parent_element_id],
            p_node->parent_nat, image);
        if (isfinite(image[0]) && isfinite(image[1]) && isfinite(image[2]))
          continue;
      }
      if (failed == 0)
        cerr << "TPL parent miss node=" << in
             << " i=" << p_node->local_i
             << " j=" << p_node->local_j
             << " k=" << p_node->local_k
             << " block=" << p_node->local_block
             << " xyz=" << p_node->coord[0] << ","
             << p_node->coord[1] << "," << p_node->coord[2] << "\n";
      failed++;
      continue;
    }
    p_node->parent_element_id = parent;
    for (int co(0); co < 3; co++)
      p_node->parent_nat[co] = nat[co];
    p_node->parent_mapping_initialized = 1;
  }
  if (failed != 0) {
    cerr << "TPL submesh parent search failed for " << failed
         << " local nodes\n";
    exit(EXIT_FAILURE);
  }
}

static void tplCopyPhysics(Mesh *p_local, const Mesh *p_coarse) {
  p_local->ZI = p_coarse->ZI;
  p_local->dx = p_coarse->dx;
  p_local->gam = p_coarse->gam;
  p_local->rhog = p_coarse->rhog;
  p_local->lgRelaxCfl = p_coarse->lgRelaxCfl;
  p_local->Vpull = p_coarse->Vpull;
  p_local->betaRough = p_coarse->betaRough;
  p_local->betaStep = p_coarse->betaStep;
  p_local->kineticA2DN = p_coarse->kineticA2DN;
  p_local->kineticB2DN = p_coarse->kineticB2DN;
  p_local->anisotropicKinFlg = p_coarse->anisotropicKinFlg;
  p_local->AniGaFlg = p_coarse->AniGaFlg;
  p_local->rChangeFlg = p_coarse->rChangeFlg;
  p_local->useVoronkovGA = p_coarse->useVoronkovGA;
  p_local->voronkovKV = p_coarse->voronkovKV;
  p_local->voronkovVelocityToCm = p_coarse->voronkovVelocityToCm;
  p_local->voronkovTm = p_coarse->voronkovTm;
  p_local->voronkovQ = p_coarse->voronkovQ;
  p_local->voronkovLambdaSG = p_coarse->voronkovLambdaSG;
  p_local->voronkovAtomicDensity = p_coarse->voronkovAtomicDensity;
  p_local->voronkovAlphaTransition = p_coarse->voronkovAlphaTransition;
  p_local->voronkovTemperatureIsAbsolute =
      p_coarse->voronkovTemperatureIsAbsolute;
  p_local->voronkovSigmaSLFacet = p_coarse->voronkovSigmaSLFacet;
  p_local->voronkovSigmaSLRough = p_coarse->voronkovSigmaSLRough;
  p_local->voronkovSigmaSLPrimeFacet =
      p_coarse->voronkovSigmaSLPrimeFacet;
  p_local->voronkovSigmaSG = p_coarse->voronkovSigmaSG;
  p_local->voronkovSigmaMG = p_coarse->voronkovSigmaMG;
  p_local->roughGrowthAngle = p_coarse->roughGrowthAngle;
  p_local->facetGrowthAngle = p_coarse->facetGrowthAngle;
  p_local->linearRoughUndercooling =
      p_coarse->linearRoughUndercooling;
  p_local->linearTipUndercooling = p_coarse->linearTipUndercooling;
  p_local->linearTipChi = p_coarse->linearTipChi;
  p_local->linearFacetAppAngle = p_coarse->linearFacetAppAngle;
  p_local->voronkovThetaFacet = p_coarse->voronkovThetaFacet;
  p_local->voronkovThetaStepCorrection =
      p_coarse->voronkovThetaStepCorrection;
}

static int tplMostInwardSLRecord(BCSet *p_sl, int record,
                                 int excludedRecord) {
  BCAtNode *p_record = &p_sl->bc_at_node[record];
  int best = -1;
  double bestRadius = 1.e300;
  double bestDirection = -2.;
  const double sourceRadius =
      sqrt(p_record->p_node->coord[0] * p_record->p_node->coord[0] +
           p_record->p_node->coord[1] * p_record->p_node->coord[1]);
  const double sourceX = p_record->p_node->coord[0] / sourceRadius;
  const double sourceY = p_record->p_node->coord[1] / sourceRadius;
  for (int ia(0); ia < p_record->num_of_adj_nodes; ia++) {
    const int candidate = p_record->n_adj_node[ia];
    if (candidate == record || candidate == excludedRecord ||
        candidate < 0 || candidate >= p_sl->number_of_recordsNode)
      continue;
    Node *p_node = p_sl->bc_at_node[candidate].p_node;
    const double radius = sqrt(p_node->coord[0] * p_node->coord[0] +
                               p_node->coord[1] * p_node->coord[1]);
    if (radius >= sourceRadius - 1.e-12)
      continue;
    const double direction =
        (p_node->coord[0] * sourceX + p_node->coord[1] * sourceY) / radius;
    if (direction > bestDirection + 1.e-12 ||
        (fabs(direction - bestDirection) <= 1.e-12 &&
         radius < bestRadius)) {
      bestDirection = direction;
      bestRadius = radius;
      best = candidate;
    }
  }
  if (best < 0) {
    cerr << "Failed to follow the coarse S/L belt inward from record "
         << record << "\n";
    exit(EXIT_FAILURE);
  }
  return best;
}

static void tplSetNode(Node *p_node, int nodeId, int i, int j, int k,
                       int block, const double point[3], int coarseNodeId) {
  p_node->node_num = nodeId;
  p_node->InodeNum = nodeId;
  p_node->coord[0] = point[0];
  p_node->coord[1] = point[1];
  p_node->coord[2] = point[2];
  p_node->Z0 = point[2];
  p_node->bound_mark = -1;
  p_node->bulk_mark = (block == 1) ? Mesh::SOLID : Mesh::LIQUID;
  p_node->parent_element_id = -1;
  p_node->coarse_node_id = coarseNodeId;
  p_node->local_i = i;
  p_node->local_j = j;
  p_node->local_k = k;
  p_node->local_block = block;
  p_node->artificial_boundary = 0;
  p_node->parent_mapping_initialized = 0;
}

static void tplUnitRadial(const double point[3], double *erx, double *ery) {
  const double radius = sqrt(point[0] * point[0] + point[1] * point[1]);
  if (radius <= 1.e-20) {
    *erx = 1.;
    *ery = 0.;
  } else {
    *erx = point[0] / radius;
    *ery = point[1] / radius;
  }
}

static void tplCurveNormal(const double a[3], const double b[3],
                           int orientation, double normal[3]) {
  double erx, ery;
  tplUnitRadial(a, &erx, &ery);
  const double ra = sqrt(a[0] * a[0] + a[1] * a[1]);
  const double rb = sqrt(b[0] * b[0] + b[1] * b[1]);
  const double dr = rb - ra;
  const double dz = b[2] - a[2];
  const double length = sqrt(dr * dr + dz * dz);
  if (length <= 1.e-20) {
    normal[0] = 0.;
    normal[1] = 0.;
    normal[2] = (orientation >= 0) ? 1. : -1.;
    return;
  }
  const double nr = orientation * dz / length;
  const double nz = -orientation * dr / length;
  normal[0] = nr * erx;
  normal[1] = nr * ery;
  normal[2] = nz;
}

static void tplNormalize3(double v[3]) {
  const double length = sqrt(v[0] * v[0] + v[1] * v[1] +
                             v[2] * v[2]);
  if (length <= 1.e-20)
    return;
  for (int co(0); co < 3; co++)
    v[co] /= length;
}

static void tplBuildBCNodes(BCSet *p_bc, int count) {
  p_bc->number_of_recordsNode = count;
  p_bc->bc_at_node = (BCAtNode *)calloc(count, sizeof(BCAtNode));
  for (int ib(0); ib < count; ib++) {
    p_bc->bc_at_node[ib].parent_face_id = -1;
    p_bc->bc_at_node[ib].parent_bc_node_id = -1;
    p_bc->bc_at_node[ib].coarse_bc_node_id = -1;
  }
}

static void tplBuildStructuredPointers(BCSet *p_bc, int nS1, int nS2) {
  p_bc->nS1 = nS1;
  p_bc->nS2 = nS2;
  p_bc->structured_wrap_columns = 0;
  p_bc->structured_z_ghost_columns = 0;
  p_bc->structured_fixed_side_columns = 0;
  p_bc->bc_at_nodeStr =
      (BCAtNode ***)malloc(nS1 * sizeof(BCAtNode **));
  for (int i(0); i < nS1; i++)
    p_bc->bc_at_nodeStr[i] =
        (BCAtNode **)malloc(nS2 * sizeof(BCAtNode *));
}

static void tplBuildQuadFaces(BCSet *p_bc, int nPhi, int nPath,
                              int periodic) {
  const int nPhiElements = periodic ? nPhi : nPhi - 1;
  p_bc->number_of_recordsFace = nPhiElements * (nPath - 1);
  p_bc->bc_at_face =
      (BCAtFace *)calloc(p_bc->number_of_recordsFace, sizeof(BCAtFace));
  int face = 0;
  for (int j(0); j < nPath - 1; j++)
    for (int i(0); i < nPhiElements; i++) {
      const int i1 = periodic ? (i + 1) % nPhi : i + 1;
      BCAtFace *p_face = &p_bc->bc_at_face[face];
      p_face->number_of_nodes = 4;
      p_face->assembly = (int *)malloc(4 * sizeof(int));
      p_face->p_node = (Node **)malloc(4 * sizeof(Node *));
      p_face->assembly[0] = i + j * nPhi;
      p_face->assembly[1] = i + (j + 1) * nPhi;
      p_face->assembly[2] = i1 + (j + 1) * nPhi;
      p_face->assembly[3] = i1 + j * nPhi;
      for (int n(0); n < 4; n++)
        p_face->p_node[n] =
            p_bc->bc_at_node[p_face->assembly[n]].p_node;
      face++;
    }
}

static void initCreateTplSubmesh(Mesh *p_local, Mesh *p_coarse,
                                  int supportLayers,
                                  int liquidRadialElements) {
  if (supportLayers < 2 || liquidRadialElements < 2) {
    cerr << "TPL submesh requires at least two supporting layers and "
            "two radial liquid elements\n";
    exit(EXIT_FAILURE);
  }
  memset(p_local, 0, sizeof(Mesh));
  tplCopyPhysics(p_local, p_coarse);
  p_local->nTPL = p_coarse->bc_set[p_coarse->BC_TPL].nS1;
  // Two coarse S/L belt rows are split so L0 has four S/L elements from the
  // TPL. The inner cut remains the coarse inner S/L ring.
  p_local->subSlRefinement = 4;
  p_local->tplSubmeshSupportLayers = supportLayers;
  p_local->tplSubmeshOuterElements = liquidRadialElements;
  p_local->number_of_coord_dir = 3;
  p_local->number_of_bc_sets = 12;
  p_local->number_of_node_attributes = p_coarse->number_of_node_attributes;
  p_local->number_of_node_markers = p_coarse->number_of_node_markers;
  p_local->number_of_groups = 2;
  const int nPhi = p_local->nTPL;
  const int ni = p_local->subSlRefinement;
  const int ns = supportLayers;
  const int nr = liquidRadialElements;
  p_local->number_of_nodes =
      nPhi * ((ni + 1) * (2 * ns + 1) + nr * (ns + 1));
  p_local->number_of_elements = (2 * ni + nr) * nPhi * ns;
  p_local->node =
      (Node *)calloc(p_local->number_of_nodes, sizeof(Node));
  p_local->element =
      (Element *)calloc(p_local->number_of_elements, sizeof(Element));
  p_local->bc_set =
      (BCSet *)calloc(p_local->number_of_bc_sets, sizeof(BCSet));
  for (int ibc(0); ibc < p_local->number_of_bc_sets; ibc++) {
    p_local->bc_set[ibc].mesh_bc_id = ibc;
    p_local->bc_set[ibc].bc_flag = BC_FACE;
  }
  for (int in(0); in < p_local->number_of_nodes; in++)
    p_local->node[in].attribute =
        (double *)calloc(p_local->number_of_node_attributes, sizeof(double));

  BCSet *p_coarse_tpl = &p_coarse->bc_set[p_coarse->BC_TPL];
  BCSet *p_coarse_sl = &p_coarse->bc_set[p_coarse->BC_SOLID_LIQUID];
  BCSet *p_coarse_sg = &p_coarse->bc_set[p_coarse->BC_SOLID_GAS];
  BCSet *p_coarse_lg = &p_coarse->bc_set[p_coarse->BC_LIQUID_GAS];
  vector<int> slRecord((ni + 1) * nPhi, -1);
  vector<int> lgRecord((nr + 1) * nPhi, -1);
  for (int i(0); i < nPhi; i++) {
    const int outer =
        p_coarse_tpl->bc_at_nodeStr[i][1]->parent_bc_node_id;
    const int middle = tplMostInwardSLRecord(p_coarse_sl, outer, -1);
    const int inner = tplMostInwardSLRecord(p_coarse_sl, middle, outer);
    slRecord[i] = outer;
    slRecord[i + nPhi] = -1;
    slRecord[i + 2 * nPhi] = middle;
    slRecord[i + 3 * nPhi] = -1;
    slRecord[i + 4 * nPhi] = inner;

    Node *p_tpl_node = p_coarse_lg->bc_at_nodeStr[i][0]->p_node;
    const double tplRadius = sqrt(
        p_tpl_node->coord[0] * p_tpl_node->coord[0] +
        p_tpl_node->coord[1] * p_tpl_node->coord[1]);
    int firstOutside = -1;
    for (int j(1); j < p_coarse_lg->nS2; j++) {
      Node *p_candidate = p_coarse_lg->bc_at_nodeStr[i][j]->p_node;
      const double radius = sqrt(
          p_candidate->coord[0] * p_candidate->coord[0] +
          p_candidate->coord[1] * p_candidate->coord[1]);
      if (radius > tplRadius + 1.e-10) {
        firstOutside = j;
        break;
      }
    }
    if (firstOutside < 0 || firstOutside + nr - 1 >= p_coarse_lg->nS2) {
      cerr << "TPL submesh needs " << nr
           << " coarse L/G rings outside the TPL\n";
      exit(EXIT_FAILURE);
    }
    lgRecord[i] = i;
    for (int j(1); j <= nr; j++)
      lgRecord[i + j * nPhi] =
          i + (firstOutside + j - 1) * nPhi;
  }

  // L0 nests artificial solid-top and liquid-bottom cuts on the coarse volume
  // node layer next to S/L when that neighbor is bulk, like L1+ nesting.
  // Vertical neighbors use brick-column stride
  // BC_SOLID_LIQUID.number_of_recordsNode. At the TPL outer ring the first
  // neighbor is on S/G or L/G; those free-surface marks are not used as
  // exact-copy cuts (they collapse under TPL/meniscus motion). Instead use
  // the proven geometric supportDepth=0.2*dx offset into the bulk, with
  // parent-element refresh (coarse_node_id=-1). TPL j=0 solid top still
  // attaches to S/G row 1 (exact parent node on S/G).
  const int slStride = p_coarse_sl->number_of_recordsNode;
  const int slTopLayer =
      p_coarse->bc_set[p_coarse->BC_COLD].bc_at_node[0].p_node->node_num /
      slStride;
  const int slStackNodes = slStride * (slTopLayer + 1);
  const int firstOutside0 = lgRecord[nPhi] / nPhi;
  const int outerRadialNodes = p_coarse_lg->nS2 - firstOutside0;
  const int outerVertStride = nPhi * outerRadialNodes;
  const double supportDepth = .2 * p_coarse->dx;
  vector<double> slPoint((ni + 1) * nPhi * 3, 0.);
  vector<double> solidTop((ni + 1) * nPhi * 3, 0.);
  vector<double> liquidBottom((ni + 1) * nPhi * 3, 0.);
  vector<double> lgPoint((nr + 1) * nPhi * 3, 0.);
  vector<double> outerBottom((nr + 1) * nPhi * 3, 0.);
  vector<int> solidTopNode((ni + 1) * nPhi, -1);
  vector<int> liquidBottomNode((ni + 1) * nPhi, -1);
  vector<int> outerBottomNode((nr + 1) * nPhi, -1);

  for (int i(0); i < nPhi; i++) {
    for (int j(0); j <= ni; j++) {
      const int rec = slRecord[i + j * nPhi];
      if (rec < 0)
        continue;
      const int ij = i + j * nPhi;
      Node *p_surface = p_coarse_sl->bc_at_node[rec].p_node;
      for (int co(0); co < 3; co++)
        slPoint[3 * ij + co] = p_surface->coord[co];
      const int solidCand = p_surface->node_num + slStride;
      const int liquidCand = p_surface->node_num - slStride;
      const int solidOk =
          solidCand >= 0 && solidCand < p_coarse->number_of_nodes &&
          p_coarse->node[solidCand].bound_mark != Mesh::BC_SOLID_GAS &&
          p_coarse->node[solidCand].bound_mark != Mesh::BC_LIQUID_GAS;
      const int liquidOk =
          liquidCand >= 0 && liquidCand < p_coarse->number_of_nodes &&
          p_coarse->node[liquidCand].bound_mark != Mesh::BC_SOLID_GAS &&
          p_coarse->node[liquidCand].bound_mark != Mesh::BC_LIQUID_GAS;
      if (solidOk) {
        solidTopNode[ij] = solidCand;
        for (int co(0); co < 3; co++)
          solidTop[3 * ij + co] = p_coarse->node[solidCand].coord[co];
      }
      if (liquidOk) {
        liquidBottomNode[ij] = liquidCand;
        for (int co(0); co < 3; co++)
          liquidBottom[3 * ij + co] =
              p_coarse->node[liquidCand].coord[co];
      }
    }
    for (int j(0); j <= ni; j++) {
      if (slRecord[i + j * nPhi] >= 0)
        continue;
      const int ij = i + j * nPhi;
      const int ij0 = i + (j - 1) * nPhi;
      const int ij1 = i + (j + 1) * nPhi;
      for (int co(0); co < 3; co++) {
        slPoint[3 * ij + co] =
            .5 * (slPoint[3 * ij0 + co] + slPoint[3 * ij1 + co]);
      }
    }

    for (int j(0); j <= ni; j++) {
      if (slRecord[i + j * nPhi] < 0)
        continue;
      const int ij = i + j * nPhi;
      if (solidTopNode[ij] >= 0 && liquidBottomNode[ij] >= 0)
        continue;
      const int j0 = (j == 0) ? 0 : j - 1;
      const int j1 = (j == ni) ? ni : j + 1;
      double a[3];
      double b[3];
      for (int co(0); co < 3; co++) {
        a[co] = slPoint[3 * (i + j0 * nPhi) + co];
        b[co] = slPoint[3 * (i + j1 * nPhi) + co];
      }
      double liquidNormal[3];
      tplCurveNormal(a, b, -1, liquidNormal);
      if (j == 0) {
        double lg0[3];
        double lg1[3];
        for (int co(0); co < 3; co++) {
          lg0[co] = p_coarse_lg->bc_at_node[lgRecord[i]].p_node->coord[co];
          lg1[co] =
              p_coarse_lg->bc_at_node[lgRecord[i + nPhi]].p_node->coord[co];
        }
        double lgLiquidNormal[3];
        tplCurveNormal(lg0, lg1, 1, lgLiquidNormal);
        for (int co(0); co < 3; co++)
          liquidNormal[co] += lgLiquidNormal[co];
        tplNormalize3(liquidNormal);
      }
      for (int co(0); co < 3; co++) {
        const double surface = slPoint[3 * ij + co];
        if (liquidBottomNode[ij] < 0)
          liquidBottom[3 * ij + co] =
              surface + supportDepth * liquidNormal[co];
        if (solidTopNode[ij] < 0)
          solidTop[3 * ij + co] =
              surface - supportDepth * liquidNormal[co];
      }
    }

    // TPL column solid top: exact next S/G parent node (keep S/G attachment).
    solidTopNode[i] = p_coarse_sg->bc_at_nodeStr[i][1]->p_node->node_num;
    for (int co(0); co < 3; co++)
      solidTop[3 * i + co] =
          p_coarse_sg->bc_at_nodeStr[i][1]->p_node->coord[co];

    // Odd S/L rows split parent belt elements: average nested neighbor cuts.
    for (int j(1); j < ni; j += 2) {
      const int ij = i + j * nPhi;
      const int ij0 = i + (j - 1) * nPhi;
      const int ij1 = i + (j + 1) * nPhi;
      solidTopNode[ij] = -1;
      liquidBottomNode[ij] = -1;
      for (int co(0); co < 3; co++) {
        liquidBottom[3 * ij + co] =
            .5 * (liquidBottom[3 * ij0 + co] +
                  liquidBottom[3 * ij1 + co]);
        solidTop[3 * ij + co] =
            .5 * (solidTop[3 * ij0 + co] +
                  solidTop[3 * ij1 + co]);
      }
    }

    for (int j(0); j <= nr; j++)
      for (int co(0); co < 3; co++)
        lgPoint[3 * (i + j * nPhi) + co] =
            p_coarse_lg->bc_at_node[lgRecord[i + j * nPhi]]
                .p_node->coord[co];

    outerBottomNode[i] = liquidBottomNode[i];
    for (int co(0); co < 3; co++)
      outerBottom[3 * i + co] = liquidBottom[3 * i + co];
    double outerA[3];
    double outerB[3];
    for (int co(0); co < 3; co++) {
      outerA[co] = lgPoint[3 * (i + (nr - 1) * nPhi) + co];
      outerB[co] = lgPoint[3 * (i + nr * nPhi) + co];
    }
    double outerLiquidNormal[3];
    tplCurveNormal(outerA, outerB, 1, outerLiquidNormal);
    double tplLiquidOffset[3];
    for (int co(0); co < 3; co++)
      tplLiquidOffset[co] =
          liquidBottom[3 * i + co] - lgPoint[3 * i + co];
    for (int j(1); j <= nr; j++) {
      const int ij = i + j * nPhi;
      Node *p_lg =
          p_coarse_lg->bc_at_node[lgRecord[ij]].p_node;
      const int belowStride =
          (p_lg->node_num < slStackNodes) ? slStride : outerVertStride;
      const int belowCand = p_lg->node_num - belowStride;
      const int belowOk =
          belowCand >= 0 && belowCand < p_coarse->number_of_nodes &&
          p_coarse->node[belowCand].bound_mark != Mesh::BC_SOLID_GAS &&
          p_coarse->node[belowCand].bound_mark != Mesh::BC_LIQUID_GAS;
      if (belowOk) {
        outerBottomNode[ij] = belowCand;
        for (int co(0); co < 3; co++)
          outerBottom[3 * ij + co] =
              p_coarse->node[belowCand].coord[co];
      } else {
        outerBottomNode[ij] = -1;
        const double u = (double)j / (double)nr;
        for (int co(0); co < 3; co++)
          outerBottom[3 * ij + co] =
              lgPoint[3 * ij + co] +
              (1. - u) * tplLiquidOffset[co] +
              u * supportDepth * outerLiquidNormal[co];
      }
    }
  }

  // Inner solid/liquid block. q=ns is S/L, q<ns is liquid, q>ns solid.
  for (int q(0); q < 2 * ns + 1; q++)
    for (int j(0); j <= ni; j++)
      for (int i(0); i < nPhi; i++) {
        const int id = tplInnerNodeId(p_local, i, j, q);
        double point[3];
        int block = 0;
        if (q == ns) {
          for (int co(0); co < 3; co++)
            point[co] = slPoint[3 * (i + j * nPhi) + co];
        } else if (q > ns) {
          const double s = (double)(q - ns) / (double)ns;
          block = 1;
          for (int co(0); co < 3; co++)
            point[co] = (1. - s) * slPoint[3 * (i + j * nPhi) + co] +
                        s * solidTop[3 * (i + j * nPhi) + co];
        } else {
          const double s = (double)(ns - q) / (double)ns;
          block = 2;
          for (int co(0); co < 3; co++)
            point[co] = (1. - s) * slPoint[3 * (i + j * nPhi) + co] +
                        s * liquidBottom[3 * (i + j * nPhi) + co];
        }
        const int ij = i + j * nPhi;
        int coarseNode = -1;
        if (q == ns && slRecord[ij] >= 0)
          coarseNode =
              p_coarse_sl->bc_at_node[slRecord[ij]].p_node->node_num;
        else if (q == 2 * ns)
          coarseNode = solidTopNode[ij];
        else if (q == 0)
          coarseNode = liquidBottomNode[ij];
        tplSetNode(&p_local->node[id], id, i, j, q - ns, block,
                   point, coarseNode);
      }

  // Outer liquid block. j=0 reuses the inner block's TPL-depth line.
  for (int k(0); k < ns + 1; k++)
    for (int j(1); j <= nr; j++)
      for (int i(0); i < nPhi; i++) {
        const int id = tplOuterNodeId(p_local, i, j, k);
        const double s = (double)k / (double)ns;
        double point[3];
        for (int co(0); co < 3; co++)
          point[co] = (1. - s) * lgPoint[3 * (i + j * nPhi) + co] +
                      s * outerBottom[3 * (i + j * nPhi) + co];
        const int ij = i + j * nPhi;
        int coarseNode = -1;
        if (k == 0)
          coarseNode =
              p_coarse_lg->bc_at_node[lgRecord[ij]].p_node->node_num;
        else if (k == ns)
          coarseNode = outerBottomNode[ij];
        tplSetNode(&p_local->node[id], id, i, j, k, 3, point,
                   coarseNode);
      }

  for (int ie(0); ie < p_local->number_of_elements; ie++) {
    Element *p_element = &p_local->element[ie];
    p_element->element_num = ie;
    p_element->number_of_nodes = 8;
    p_element->bf_set_flag = BRICK_8_NOD_LAGR;
    p_element->assembly = (int *)malloc(8 * sizeof(int));
    p_element->p_node = (Node **)malloc(8 * sizeof(Node *));
  }

  int element = 0;
  // Solid: lower face is S/L, upper face is the solid cut.
  for (int k(0); k < ns; k++)
    for (int j(0); j < ni; j++)
      for (int i(0); i < nPhi; i++) {
        const int i1 = (i + 1) % nPhi;
        Element *p_element = &p_local->element[element++];
        p_element->group_num = Mesh::SOLID;
        p_element->assembly[0] = tplInnerNodeId(p_local, i, j + 1, ns + k);
        p_element->assembly[1] = tplInnerNodeId(p_local, i, j, ns + k);
        p_element->assembly[2] = tplInnerNodeId(p_local, i1, j, ns + k);
        p_element->assembly[3] = tplInnerNodeId(p_local, i1, j + 1, ns + k);
        p_element->assembly[4] = tplInnerNodeId(p_local, i, j + 1, ns + k + 1);
        p_element->assembly[5] = tplInnerNodeId(p_local, i, j, ns + k + 1);
        p_element->assembly[6] = tplInnerNodeId(p_local, i1, j, ns + k + 1);
        p_element->assembly[7] = tplInnerNodeId(p_local, i1, j + 1, ns + k + 1);
      }
  // Inner liquid: lower face is farther from S/L.
  for (int k(0); k < ns; k++)
    for (int j(0); j < ni; j++)
      for (int i(0); i < nPhi; i++) {
        const int i1 = (i + 1) % nPhi;
        Element *p_element = &p_local->element[element++];
        p_element->group_num = Mesh::LIQUID;
        p_element->assembly[0] = tplInnerNodeId(p_local, i, j + 1, ns - k - 1);
        p_element->assembly[1] = tplInnerNodeId(p_local, i, j, ns - k - 1);
        p_element->assembly[2] = tplInnerNodeId(p_local, i1, j, ns - k - 1);
        p_element->assembly[3] = tplInnerNodeId(p_local, i1, j + 1, ns - k - 1);
        p_element->assembly[4] = tplInnerNodeId(p_local, i, j + 1, ns - k);
        p_element->assembly[5] = tplInnerNodeId(p_local, i, j, ns - k);
        p_element->assembly[6] = tplInnerNodeId(p_local, i1, j, ns - k);
        p_element->assembly[7] = tplInnerNodeId(p_local, i1, j + 1, ns - k);
      }
  // Outer liquid: top is L/G and k increases into liquid.
  for (int k(0); k < ns; k++)
    for (int j(0); j < nr; j++)
      for (int i(0); i < nPhi; i++) {
        const int i1 = (i + 1) % nPhi;
        Element *p_element = &p_local->element[element++];
        p_element->group_num = Mesh::LIQUID;
        p_element->assembly[0] = tplOuterNodeId(p_local, i, j, k + 1);
        p_element->assembly[1] = tplOuterNodeId(p_local, i, j + 1, k + 1);
        p_element->assembly[2] = tplOuterNodeId(p_local, i1, j + 1, k + 1);
        p_element->assembly[3] = tplOuterNodeId(p_local, i1, j, k + 1);
        p_element->assembly[4] = tplOuterNodeId(p_local, i, j, k);
        p_element->assembly[5] = tplOuterNodeId(p_local, i, j + 1, k);
        p_element->assembly[6] = tplOuterNodeId(p_local, i1, j + 1, k);
        p_element->assembly[7] = tplOuterNodeId(p_local, i1, j, k);
      }
  if (element != p_local->number_of_elements) {
    cerr << "TPL submesh element count mismatch\n";
    exit(EXIT_FAILURE);
  }
  for (int ie(0); ie < p_local->number_of_elements; ie++)
    for (int n(0); n < 8; n++)
      p_local->element[ie].p_node[n] =
          &p_local->node[p_local->element[ie].assembly[n]];

  BCSet *p_sl = &p_local->bc_set[p_local->BC_SUB_SOLID_LIQUID];
  tplBuildBCNodes(p_sl, (ni + 1) * nPhi);
  tplBuildStructuredPointers(p_sl, nPhi, ni + 1);
  for (int j(0); j <= ni; j++)
    for (int i(0); i < nPhi; i++) {
      const int record = i + j * nPhi;
      p_sl->bc_at_node[record].p_node =
          &p_local->node[tplInnerNodeId(p_local, i, j, ns)];
      p_sl->bc_at_node[record].node_num = record;
      p_sl->bc_at_node[record].coarse_bc_node_id =
          slRecord[i + j * nPhi];
      p_sl->bc_at_node[record].relax_coef = 1.;
      p_sl->bc_at_node[record].refine_coef = 1.;
      p_sl->bc_at_nodeStr[i][j] = &p_sl->bc_at_node[record];
      p_sl->bc_at_node[record].p_node->bound_mark =
          (j == ni) ? Mesh::BC_SUB : Mesh::BC_SUB_SOLID_LIQUID;
    }
  tplBuildQuadFaces(p_sl, nPhi, ni + 1, 1);
  p_local->initAdjNodesQuad(p_sl);
  p_sl->relaxMoveCoef =
      p_coarse->bc_set[p_coarse->BC_SOLID_LIQUID].relaxMoveCoef;

  BCSet *p_sg = &p_local->bc_set[p_local->BC_SUB_SOLID_GAS];
  tplBuildBCNodes(p_sg, nPhi * (ns + 1));
  tplBuildStructuredPointers(p_sg, nPhi, ns + 1);
  for (int k(0); k < ns + 1; k++)
    for (int i(0); i < nPhi; i++) {
      const int record = i + k * nPhi;
      p_sg->bc_at_node[record].p_node =
          &p_local->node[tplInnerNodeId(p_local, i, 0, ns + k)];
      p_sg->bc_at_node[record].node_num = record;
      p_sg->bc_at_node[record].relax_coef = 1.;
      p_sg->bc_at_node[record].refine_coef = 1.;
      p_sg->bc_at_nodeStr[i][k] = &p_sg->bc_at_node[record];
      if (k > 0 && k < ns)
        p_sg->bc_at_node[record].p_node->bound_mark =
            Mesh::BC_SUB_SOLID_GAS;
      if (k == ns) {
        p_sg->bc_at_node[record].p_node->bound_mark = Mesh::BC_SUB;
        p_sg->bc_at_node[record].p_node->artificial_boundary = 1;
      }
    }
  tplBuildQuadFaces(p_sg, nPhi, ns + 1, 1);

  BCSet *p_lg = &p_local->bc_set[p_local->BC_SUB_LIQUID_GAS];
  tplBuildBCNodes(p_lg, (nr + 1) * nPhi);
  tplBuildStructuredPointers(p_lg, nPhi, nr + 1);
  for (int j(0); j <= nr; j++)
    for (int i(0); i < nPhi; i++) {
      const int record = i + j * nPhi;
      p_lg->bc_at_node[record].p_node =
          &p_local->node[tplOuterNodeId(p_local, i, j, 0)];
      p_lg->bc_at_node[record].node_num = record;
      p_lg->bc_at_node[record].coarse_bc_node_id =
          lgRecord[i + j * nPhi];
      p_lg->bc_at_node[record].relax_coef = 1.;
      p_lg->bc_at_node[record].refine_coef = 1.;
      p_lg->bc_at_nodeStr[i][j] = &p_lg->bc_at_node[record];
      if (j > 0 && j < nr)
        p_lg->bc_at_node[record].p_node->bound_mark =
            Mesh::BC_SUB_LIQUID_GAS;
      if (j == nr) {
        p_lg->bc_at_node[record].p_node->bound_mark = Mesh::BC_SUB;
        p_lg->bc_at_node[record].p_node->artificial_boundary = 1;
      }
    }
  tplBuildQuadFaces(p_lg, nPhi, nr + 1, 1);

  BCSet *p_tpl = &p_local->bc_set[p_local->BC_SUB_TPL];
  tplBuildBCNodes(p_tpl, 3 * nPhi);
  tplBuildStructuredPointers(p_tpl, nPhi, 3);
  for (int i(0); i < nPhi; i++) {
    p_tpl->bc_at_node[i].p_node = p_sg->bc_at_nodeStr[i][1]->p_node;
    p_tpl->bc_at_node[i + nPhi].p_node = p_sl->bc_at_nodeStr[i][0]->p_node;
    p_tpl->bc_at_node[i + 2 * nPhi].p_node = p_lg->bc_at_nodeStr[i][1]->p_node;
  }
  for (int j(0); j < 3; j++)
    for (int i(0); i < nPhi; i++) {
      const int record = i + j * nPhi;
      p_tpl->bc_at_node[record].node_num = record;
      p_tpl->bc_at_node[record].parent_face_id = -1;
      p_tpl->bc_at_node[record].parent_bc_node_id = -1;
      p_tpl->bc_at_node[record].coarse_bc_node_id = -1;
      p_tpl->bc_at_node[record].relax_coef = 1.;
      p_tpl->bc_at_node[record].refine_coef = 1.;
      p_tpl->bc_at_nodeStr[i][j] = &p_tpl->bc_at_node[record];
    }
  for (int i(0); i < nPhi; i++) {
    BCAtNode *p_contact = p_tpl->bc_at_nodeStr[i][1];
    p_contact->p_adj_sg_node = p_sg->bc_at_nodeStr[i][1]->p_node;
    p_contact->p_adj_lg_node = p_lg->bc_at_nodeStr[i][1]->p_node;
    p_contact->parent_bc_node_id = i;
    p_contact->coarse_bc_node_id = i + nPhi;
  }

  // Mark all remaining artificial-cut nodes.  The complete ring has no
  // azimuthal side cuts.
  vector<int> artificial;
  vector<int> seen(p_local->number_of_nodes, 0);
  for (int i(0); i < nPhi; i++) {
    for (int j(0); j <= ni; j++) {
      const int solidTopId = tplInnerNodeId(p_local, i, j, 2 * ns);
      const int liquidBottomId = tplInnerNodeId(p_local, i, j, 0);
      if (!seen[solidTopId]) {
        seen[solidTopId] = 1;
        artificial.push_back(solidTopId);
      }
      if (!seen[liquidBottomId]) {
        seen[liquidBottomId] = 1;
        artificial.push_back(liquidBottomId);
      }
    }
    for (int q(0); q < 2 * ns + 1; q++) {
      const int innerCutId = tplInnerNodeId(p_local, i, ni, q);
      if (!seen[innerCutId]) {
        seen[innerCutId] = 1;
        artificial.push_back(innerCutId);
      }
    }
    for (int j(0); j <= nr; j++) {
      const int lowerId = tplOuterNodeId(p_local, i, j, ns);
      if (!seen[lowerId]) {
        seen[lowerId] = 1;
        artificial.push_back(lowerId);
      }
    }
    for (int k(0); k < ns + 1; k++) {
      const int outerId = tplOuterNodeId(p_local, i, nr, k);
      if (!seen[outerId]) {
        seen[outerId] = 1;
        artificial.push_back(outerId);
      }
    }
  }
  BCSet *p_artificial = &p_local->bc_set[p_local->BC_SUB];
  tplBuildBCNodes(p_artificial, (int)artificial.size());
  for (int ib(0); ib < (int)artificial.size(); ib++) {
    Node *p_node = &p_local->node[artificial[ib]];
    p_node->artificial_boundary = 1;
    p_node->bound_mark = Mesh::BC_SUB;
    p_artificial->bc_at_node[ib].p_node = p_node;
    p_artificial->bc_at_node[ib].node_num = ib;
  }

  tplMapLocalNodes(p_local, p_coarse);
  for (int ib(0); ib < p_artificial->number_of_recordsNode; ib++) {
    Node *p_node = p_artificial->bc_at_node[ib].p_node;
    p_artificial->bc_at_node[ib].parent_face_id =
        p_node->parent_element_id;
    p_artificial->bc_at_node[ib].parent_xi = p_node->parent_nat[0];
    p_artificial->bc_at_node[ib].parent_eta = p_node->parent_nat[1];
    p_artificial->bc_at_node[ib].parent_zeta = p_node->parent_nat[2];
  }
  p_local->initBCEdges(p_sl);
  p_local->initBCEdges(p_sg);
  p_local->initBCEdges(p_lg);
  p_local->initBCEdges(p_tpl);
  p_local->initEdges();
}

static int tplExactGridIndex(double value) {
  const int index = (int)floor(value + .5);
  return (fabs(value - index) < 1.e-12) ? index : -1;
}

static void tplLinearGridStencil(double coordinate, int extent,
                                 int periodic, int index[2],
                                 double weight[2]) {
  // At a nonperiodic cut, coordinate 0 uses the first available interval
  // (forward form) and coordinate extent uses the last interval (backward
  // form). Interior points use their enclosing parent interval.
  int lower = (int)floor(coordinate);
  double x = coordinate - lower;
  if (!periodic) {
    if (lower < 0) {
      lower = 0;
      x = 0.;
    } else if (lower >= extent) {
      lower = extent - 1;
      x = 1.;
    }
  }
  index[0] = lower;
  index[1] = lower + 1;
  weight[0] = 1. - x;
  weight[1] = x;
  if (periodic)
    for (int a(0); a < 2; a++)
      index[a] = (index[a] % extent + extent) % extent;
  else
    for (int a(0); a < 2; a++)
      index[a] = max(0, min(extent, index[a]));
}

static void tplSampleInnerGrid(const Mesh *p_parent, double phi,
                               double path, double depth,
                               double point[3], int *p_parentNode) {
  const int nPhi = p_parent->nTPL;
  const int ni = p_parent->subSlRefinement;
  const int ns = p_parent->tplSubmeshSupportLayers;
  int ii[2], jj[2], qq[2];
  double wi[2], wj[2], wq[2];
  tplLinearGridStencil(phi, nPhi, 1, ii, wi);
  tplLinearGridStencil(path, ni, 0, jj, wj);
  tplLinearGridStencil(depth, 2 * ns, 0, qq, wq);
  point[0] = point[1] = point[2] = 0.;
  for (int ai(0); ai < 2; ai++)
    for (int aj(0); aj < 2; aj++)
      for (int aq(0); aq < 2; aq++) {
        const double weight = wi[ai] * wj[aj] * wq[aq];
        const Node *p_node =
            &p_parent->node[tplInnerNodeId(
                p_parent, ii[ai], jj[aj], qq[aq])];
        for (int co(0); co < 3; co++)
          point[co] += weight * p_node->coord[co];
      }
  const int ie = tplExactGridIndex(phi);
  const int je = tplExactGridIndex(path);
  const int qe = tplExactGridIndex(depth);
  *p_parentNode = (ie >= 0 && ie < nPhi &&
                   je >= 0 && je <= ni && qe >= 0 && qe <= 2 * ns)
                      ? tplInnerNodeId(p_parent, ie, je, qe)
                      : -1;
}

static void tplSampleOuterGrid(const Mesh *p_parent, double phi,
                               double path, double depth,
                               double point[3], int *p_parentNode) {
  const int nPhi = p_parent->nTPL;
  const int ns = p_parent->tplSubmeshSupportLayers;
  const int nr = p_parent->tplSubmeshOuterElements;
  int ii[2], jj[2], kk[2];
  double wi[2], wj[2], wk[2];
  tplLinearGridStencil(phi, nPhi, 1, ii, wi);
  tplLinearGridStencil(path, nr, 0, jj, wj);
  tplLinearGridStencil(depth, ns, 0, kk, wk);
  point[0] = point[1] = point[2] = 0.;
  for (int ai(0); ai < 2; ai++)
    for (int aj(0); aj < 2; aj++)
      for (int ak(0); ak < 2; ak++) {
        const double weight = wi[ai] * wj[aj] * wk[ak];
        const Node *p_node =
            &p_parent->node[tplOuterNodeId(
                p_parent, ii[ai], jj[aj], kk[ak])];
        for (int co(0); co < 3; co++)
          point[co] += weight * p_node->coord[co];
      }
  const int ie = tplExactGridIndex(phi);
  const int je = tplExactGridIndex(path);
  const int ke = tplExactGridIndex(depth);
  *p_parentNode = (ie >= 0 && ie < nPhi && je >= 0 && ke >= 0)
                      ? tplOuterNodeId(p_parent, ie, je, ke)
                      : -1;
}

// Refine the complete periodic ring from an existing structured TPL submesh.
// The child splits selected immediate-parent S/L elements next to the TPL.
// Volume blocks stop at the parent bulk node next to each parent cut, so
// each child block sits inside its parent. Physical S/L, S/G, and L/G stay
// on the parent surfaces. Exact child cut nodes coincide with parent volume
// nodes; intermediate coordinates use shape-preserving linear interpolation.
static void initCreateRefinedTplSubmesh(
    Mesh *p_local, Mesh *p_parent, int phiRefinement) {
  if (p_parent->bc_set == NULL ||
      p_parent->bc_set[p_parent->BC_SUB_TPL].bc_at_nodeStr == NULL ||
      p_parent->nTPL < 3 || phiRefinement < 2) {
    cerr << "Invalid recursive circular TPL-submesh request\n";
    exit(EXIT_FAILURE);
  }

  memset(p_local, 0, sizeof(Mesh));
  tplCopyPhysics(p_local, p_parent);
  const int nPhi = p_parent->nTPL * phiRefinement;
  const int parentNs = p_parent->tplSubmeshSupportLayers;
  const int parentNr = p_parent->tplSubmeshOuterElements;
  const int parentNi = p_parent->subSlRefinement;
  // Nested volume: stop at the parent bulk node next to each parent-controlled
  // cut. Surfaces stay on the parent S/L, S/G, and L/G. L0 has four S/L
  // elements; L1 splits three into six; L2 splits five into ten; L3 splits
  // nine into eighteen. The same parentN-1 then 2*(parentN-1) rule is used
  // for the normal and outer-liquid blocks.
  const int parentInnerElements = parentNi - 1;
  const int parentInnerLayers = parentNs - 1;
  const int parentInnerRadial = parentNr - 1;
  if (parentInnerElements < 1 || parentInnerLayers < 1 ||
      parentInnerRadial < 1) {
    cerr << "Recursive TPL submesh needs a parent bulk node next to "
            "each artificial cut\n";
    exit(EXIT_FAILURE);
  }
  const int ni = 2 * parentInnerElements;
  const int ns = 2 * parentInnerLayers;
  const int nr = 2 * parentInnerRadial;
  p_local->nTPL = nPhi;
  p_local->subSlRefinement = ni;
  p_local->tplSubmeshSupportLayers = ns;
  p_local->tplSubmeshOuterElements = nr;
  p_local->number_of_coord_dir = 3;
  p_local->number_of_bc_sets = 12;
  p_local->number_of_node_attributes = p_parent->number_of_node_attributes;
  p_local->number_of_node_markers = p_parent->number_of_node_markers;
  p_local->number_of_groups = 2;
  p_local->number_of_nodes =
      nPhi * ((ni + 1) * (2 * ns + 1) + nr * (ns + 1));
  p_local->number_of_elements = (2 * ni + nr) * nPhi * ns;
  p_local->node =
      (Node *)calloc(p_local->number_of_nodes, sizeof(Node));
  p_local->element =
      (Element *)calloc(p_local->number_of_elements, sizeof(Element));
  p_local->bc_set =
      (BCSet *)calloc(p_local->number_of_bc_sets, sizeof(BCSet));
  for (int ibc(0); ibc < p_local->number_of_bc_sets; ibc++) {
    p_local->bc_set[ibc].mesh_bc_id = ibc;
    p_local->bc_set[ibc].bc_flag = BC_FACE;
  }
  for (int in(0); in < p_local->number_of_nodes; in++)
    p_local->node[in].attribute =
        (double *)calloc(p_local->number_of_node_attributes, sizeof(double));

  vector<double> parentPhi(nPhi, 0.);
  for (int i(0); i < nPhi; i++)
    parentPhi[i] = (double)i / (double)phiRefinement;

  for (int q(0); q < 2 * ns + 1; q++)
    for (int j(0); j <= ni; j++)
      for (int i(0); i < nPhi; i++) {
        // Split each selected immediate-parent S/L element in two. The child
        // inner cut is the last relaxed parent node, next to the parent
        // BC_SUB cut.
        const double parentPath = .5 * (double)j;
        const int childDistance = abs(q - ns);
        const double parentDistance = .5 * (double)childDistance;
        const double parentDepth =
            parentNs + ((q >= ns) ? parentDistance : -parentDistance);
        double point[3];
        int parentNode = -1;
        tplSampleInnerGrid(p_parent, parentPhi[i], parentPath,
                           parentDepth, point, &parentNode);
        const int block = (q > ns) ? 1 : ((q < ns) ? 2 : 0);
        const int id = tplInnerNodeId(p_local, i, j, q);
        tplSetNode(&p_local->node[id], id, i, j, q - ns, block,
                   point, parentNode);
      }

  for (int k(0); k < ns + 1; k++)
    for (int j(1); j <= nr; j++)
      for (int i(0); i < nPhi; i++) {
        // Split the parent bulk, stopping at the parent node next to the
        // outer artificial cut.
        const double parentPath = .5 * (double)j;
        const double parentDepth = .5 * (double)k;
        double point[3];
        int parentNode = -1;
        tplSampleOuterGrid(p_parent, parentPhi[i], parentPath,
                           parentDepth, point, &parentNode);
        const int id = tplOuterNodeId(p_local, i, j, k);
        tplSetNode(&p_local->node[id], id, i, j, k, 3,
                   point, parentNode);
      }

  for (int ie(0); ie < p_local->number_of_elements; ie++) {
    Element *p_element = &p_local->element[ie];
    p_element->element_num = ie;
    p_element->number_of_nodes = 8;
    p_element->bf_set_flag = BRICK_8_NOD_LAGR;
    p_element->assembly = (int *)malloc(8 * sizeof(int));
    p_element->p_node = (Node **)malloc(8 * sizeof(Node *));
  }

  int element = 0;
  for (int k(0); k < ns; k++)
    for (int j(0); j < ni; j++)
      for (int i(0); i < nPhi; i++) {
        const int i1 = (i + 1) % nPhi;
        Element *p_element = &p_local->element[element++];
        p_element->group_num = Mesh::SOLID;
        p_element->assembly[0] = tplInnerNodeId(p_local, i, j + 1, ns + k);
        p_element->assembly[1] = tplInnerNodeId(p_local, i, j, ns + k);
        p_element->assembly[2] = tplInnerNodeId(p_local, i1, j, ns + k);
        p_element->assembly[3] = tplInnerNodeId(p_local, i1, j + 1, ns + k);
        p_element->assembly[4] = tplInnerNodeId(p_local, i, j + 1, ns + k + 1);
        p_element->assembly[5] = tplInnerNodeId(p_local, i, j, ns + k + 1);
        p_element->assembly[6] = tplInnerNodeId(p_local, i1, j, ns + k + 1);
        p_element->assembly[7] = tplInnerNodeId(p_local, i1, j + 1, ns + k + 1);
      }
  for (int k(0); k < ns; k++)
    for (int j(0); j < ni; j++)
      for (int i(0); i < nPhi; i++) {
        const int i1 = (i + 1) % nPhi;
        Element *p_element = &p_local->element[element++];
        p_element->group_num = Mesh::LIQUID;
        p_element->assembly[0] = tplInnerNodeId(p_local, i, j + 1, ns - k - 1);
        p_element->assembly[1] = tplInnerNodeId(p_local, i, j, ns - k - 1);
        p_element->assembly[2] = tplInnerNodeId(p_local, i1, j, ns - k - 1);
        p_element->assembly[3] = tplInnerNodeId(p_local, i1, j + 1, ns - k - 1);
        p_element->assembly[4] = tplInnerNodeId(p_local, i, j + 1, ns - k);
        p_element->assembly[5] = tplInnerNodeId(p_local, i, j, ns - k);
        p_element->assembly[6] = tplInnerNodeId(p_local, i1, j, ns - k);
        p_element->assembly[7] = tplInnerNodeId(p_local, i1, j + 1, ns - k);
      }
  for (int k(0); k < ns; k++)
    for (int j(0); j < nr; j++)
      for (int i(0); i < nPhi; i++) {
        const int i1 = (i + 1) % nPhi;
        Element *p_element = &p_local->element[element++];
        p_element->group_num = Mesh::LIQUID;
        p_element->assembly[0] = tplOuterNodeId(p_local, i, j, k + 1);
        p_element->assembly[1] = tplOuterNodeId(p_local, i, j + 1, k + 1);
        p_element->assembly[2] = tplOuterNodeId(p_local, i1, j + 1, k + 1);
        p_element->assembly[3] = tplOuterNodeId(p_local, i1, j, k + 1);
        p_element->assembly[4] = tplOuterNodeId(p_local, i, j, k);
        p_element->assembly[5] = tplOuterNodeId(p_local, i, j + 1, k);
        p_element->assembly[6] = tplOuterNodeId(p_local, i1, j + 1, k);
        p_element->assembly[7] = tplOuterNodeId(p_local, i1, j, k);
      }
  if (element != p_local->number_of_elements) {
    cerr << "Recursive circular TPL-submesh element count mismatch\n";
    exit(EXIT_FAILURE);
  }
  for (int ie(0); ie < p_local->number_of_elements; ie++)
    for (int n(0); n < 8; n++)
      p_local->element[ie].p_node[n] =
          &p_local->node[p_local->element[ie].assembly[n]];

  BCSet *p_sl = &p_local->bc_set[p_local->BC_SUB_SOLID_LIQUID];
  tplBuildBCNodes(p_sl, (ni + 1) * nPhi);
  tplBuildStructuredPointers(p_sl, nPhi, ni + 1);
  for (int j(0); j <= ni; j++)
    for (int i(0); i < nPhi; i++) {
      const int record = i + j * nPhi;
      p_sl->bc_at_node[record].p_node =
          &p_local->node[tplInnerNodeId(p_local, i, j, ns)];
      p_sl->bc_at_node[record].node_num = record;
      const int pi = tplExactGridIndex(parentPhi[i]);
      const double parentPath = .5 * (double)j;
      const int pj = tplExactGridIndex(parentPath);
      p_sl->bc_at_node[record].coarse_bc_node_id =
          (pi >= 0 && pj >= 0 && pj <= parentNi)
              ? pi + pj * p_parent->nTPL
              : -1;
      p_sl->bc_at_node[record].relax_coef = 1.;
      p_sl->bc_at_node[record].refine_coef = 1.;
      p_sl->bc_at_nodeStr[i][j] = &p_sl->bc_at_node[record];
      p_sl->bc_at_node[record].p_node->bound_mark =
          (j == ni) ? Mesh::BC_SUB : Mesh::BC_SUB_SOLID_LIQUID;
    }
  tplBuildQuadFaces(p_sl, nPhi, ni + 1, 1);
  p_local->initAdjNodesQuad(p_sl);
  p_sl->relaxMoveCoef =
      p_parent->bc_set[p_parent->BC_SUB_SOLID_LIQUID].relaxMoveCoef;

  BCSet *p_sg = &p_local->bc_set[p_local->BC_SUB_SOLID_GAS];
  tplBuildBCNodes(p_sg, nPhi * (ns + 1));
  tplBuildStructuredPointers(p_sg, nPhi, ns + 1);
  for (int k(0); k < ns + 1; k++)
    for (int i(0); i < nPhi; i++) {
      const int record = i + k * nPhi;
      p_sg->bc_at_node[record].p_node =
          &p_local->node[tplInnerNodeId(p_local, i, 0, ns + k)];
      p_sg->bc_at_node[record].node_num = record;
      p_sg->bc_at_node[record].relax_coef = 1.;
      p_sg->bc_at_node[record].refine_coef = 1.;
      p_sg->bc_at_nodeStr[i][k] = &p_sg->bc_at_node[record];
      if (k > 0 && k < ns)
        p_sg->bc_at_node[record].p_node->bound_mark =
            Mesh::BC_SUB_SOLID_GAS;
      if (k == ns) {
        p_sg->bc_at_node[record].p_node->bound_mark = Mesh::BC_SUB;
        p_sg->bc_at_node[record].p_node->artificial_boundary = 1;
      }
    }
  tplBuildQuadFaces(p_sg, nPhi, ns + 1, 1);

  BCSet *p_lg = &p_local->bc_set[p_local->BC_SUB_LIQUID_GAS];
  tplBuildBCNodes(p_lg, (nr + 1) * nPhi);
  tplBuildStructuredPointers(p_lg, nPhi, nr + 1);
  for (int j(0); j <= nr; j++)
    for (int i(0); i < nPhi; i++) {
      const int record = i + j * nPhi;
      p_lg->bc_at_node[record].p_node =
          &p_local->node[tplOuterNodeId(p_local, i, j, 0)];
      p_lg->bc_at_node[record].node_num = record;
      const int pi = tplExactGridIndex(parentPhi[i]);
      const double parentPath = .5 * (double)j;
      const int pj = tplExactGridIndex(parentPath);
      p_lg->bc_at_node[record].coarse_bc_node_id =
          (pi >= 0 && pj >= 0) ? pi + pj * p_parent->nTPL : -1;
      p_lg->bc_at_node[record].relax_coef = 1.;
      p_lg->bc_at_node[record].refine_coef = 1.;
      p_lg->bc_at_nodeStr[i][j] = &p_lg->bc_at_node[record];
      if (j > 0 && j < nr)
        p_lg->bc_at_node[record].p_node->bound_mark =
            Mesh::BC_SUB_LIQUID_GAS;
      if (j == nr) {
        p_lg->bc_at_node[record].p_node->bound_mark = Mesh::BC_SUB;
        p_lg->bc_at_node[record].p_node->artificial_boundary = 1;
      }
    }
  tplBuildQuadFaces(p_lg, nPhi, nr + 1, 1);

  BCSet *p_tpl = &p_local->bc_set[p_local->BC_SUB_TPL];
  tplBuildBCNodes(p_tpl, 3 * nPhi);
  tplBuildStructuredPointers(p_tpl, nPhi, 3);
  for (int i(0); i < nPhi; i++) {
    p_tpl->bc_at_node[i].p_node = p_sg->bc_at_nodeStr[i][1]->p_node;
    p_tpl->bc_at_node[i + nPhi].p_node = p_sl->bc_at_nodeStr[i][0]->p_node;
    p_tpl->bc_at_node[i + 2 * nPhi].p_node = p_lg->bc_at_nodeStr[i][1]->p_node;
  }
  for (int j(0); j < 3; j++)
    for (int i(0); i < nPhi; i++) {
      const int record = i + j * nPhi;
      p_tpl->bc_at_node[record].node_num = record;
      p_tpl->bc_at_node[record].parent_face_id = -1;
      p_tpl->bc_at_node[record].parent_bc_node_id = -1;
      p_tpl->bc_at_node[record].coarse_bc_node_id = -1;
      p_tpl->bc_at_node[record].relax_coef = 1.;
      p_tpl->bc_at_node[record].refine_coef = 1.;
      p_tpl->bc_at_nodeStr[i][j] = &p_tpl->bc_at_node[record];
    }
  for (int i(0); i < nPhi; i++) {
    BCAtNode *p_contact = p_tpl->bc_at_nodeStr[i][1];
    p_contact->p_adj_sg_node = p_sg->bc_at_nodeStr[i][1]->p_node;
    p_contact->p_adj_lg_node = p_lg->bc_at_nodeStr[i][1]->p_node;
    p_contact->parent_bc_node_id = i;
    const int pi = tplExactGridIndex(parentPhi[i]);
    p_contact->coarse_bc_node_id =
        (pi >= 0) ? pi + p_parent->nTPL : -1;
  }

  vector<int> artificial;
  vector<int> seen(p_local->number_of_nodes, 0);
  for (int i(0); i < nPhi; i++) {
    for (int j(0); j <= ni; j++) {
      const int top = tplInnerNodeId(p_local, i, j, 2 * ns);
      const int bottom = tplInnerNodeId(p_local, i, j, 0);
      if (!seen[top]) { seen[top] = 1; artificial.push_back(top); }
      if (!seen[bottom]) { seen[bottom] = 1; artificial.push_back(bottom); }
    }
    for (int q(0); q < 2 * ns + 1; q++) {
      const int inner = tplInnerNodeId(p_local, i, ni, q);
      if (!seen[inner]) { seen[inner] = 1; artificial.push_back(inner); }
    }
    for (int j(0); j <= nr; j++) {
      const int lower = tplOuterNodeId(p_local, i, j, ns);
      if (!seen[lower]) { seen[lower] = 1; artificial.push_back(lower); }
    }
    for (int k(0); k < ns + 1; k++) {
      const int outer = tplOuterNodeId(p_local, i, nr, k);
      if (!seen[outer]) { seen[outer] = 1; artificial.push_back(outer); }
    }
  }
  BCSet *p_artificial = &p_local->bc_set[p_local->BC_SUB];
  tplBuildBCNodes(p_artificial, (int)artificial.size());
  for (int ib(0); ib < (int)artificial.size(); ib++) {
    Node *p_node = &p_local->node[artificial[ib]];
    p_node->artificial_boundary = 1;
    p_node->bound_mark = Mesh::BC_SUB;
    p_artificial->bc_at_node[ib].p_node = p_node;
    p_artificial->bc_at_node[ib].node_num = ib;
  }

  tplMapLocalNodes(p_local, p_parent);
  for (int ib(0); ib < p_artificial->number_of_recordsNode; ib++) {
    Node *p_node = p_artificial->bc_at_node[ib].p_node;
    p_artificial->bc_at_node[ib].parent_face_id =
        p_node->parent_element_id;
    p_artificial->bc_at_node[ib].parent_xi = p_node->parent_nat[0];
    p_artificial->bc_at_node[ib].parent_eta = p_node->parent_nat[1];
    p_artificial->bc_at_node[ib].parent_zeta = p_node->parent_nat[2];
  }
  p_local->initBCEdges(p_sl);
  p_local->initBCEdges(p_sg);
  p_local->initBCEdges(p_lg);
  p_local->initBCEdges(p_tpl);
  p_local->initEdges();
}

[[maybe_unused]] static void refreshTplSubmeshArtificialGeometry(
    Mesh *p_local, Mesh *p_parent) {
  // Boundary geometry is parent controlled. Exact coincident nodes are
  // copied; noncoincident recursive nodes use shape-preserving trilinear
  // interpolation of the immediate parent. Temperature interpolation remains
  // second order and is evaluated separately below.
  BCSet *p_artificial = &p_local->bc_set[p_local->BC_SUB];
  const int parentIsTplSubmesh = tplParentIsTplSubmesh(p_parent);
  const int phiRefinement =
      parentIsTplSubmesh ? p_local->nTPL / p_parent->nTPL : 1;
  const int parentNs = p_parent->tplSubmeshSupportLayers;
  const int childNs = p_local->tplSubmeshSupportLayers;
  for (int ib(0); ib < p_artificial->number_of_recordsNode; ib++) {
    Node *p_node = p_artificial->bc_at_node[ib].p_node;
    if (p_node->coarse_node_id >= 0 &&
        p_node->coarse_node_id < p_parent->number_of_nodes) {
      Node *p_src = &p_parent->node[p_node->coarse_node_id];
      if (isfinite(p_src->coord[0]) && isfinite(p_src->coord[1]) &&
          isfinite(p_src->coord[2]))
        for (int co(0); co < 3; co++)
          p_node->coord[co] = p_src->coord[co];
    }
    else if (parentIsTplSubmesh) {
      const double parentPhi =
          (double)p_node->local_i / (double)phiRefinement;
      double point[3];
      int parentNode = -1;
      if (p_node->local_block != 3) {
        const int q = p_node->local_k + childNs;
        const int childDistance = abs(q - childNs);
        const double parentDistance = .5 * (double)childDistance;
        const double parentDepth =
            parentNs + ((q >= childNs) ? parentDistance
                                       : -parentDistance);
        tplSampleInnerGrid(
            p_parent, parentPhi, .5 * p_node->local_j,
            parentDepth, point, &parentNode);
      } else {
        const double parentPath = .5 * (double)p_node->local_j;
        const double parentDepth = .5 * (double)p_node->local_k;
        tplSampleOuterGrid(
            p_parent, parentPhi, parentPath, parentDepth,
            point, &parentNode);
      }
      if (isfinite(point[0]) && isfinite(point[1]) && isfinite(point[2]))
        for (int co(0); co < 3; co++)
          p_node->coord[co] = point[co];
    } else if (p_node->parent_element_id >= 0 &&
               p_node->parent_element_id < p_parent->number_of_elements) {
      double point[3];
      tplInterpolateElementPoint(
          &p_parent->element[p_node->parent_element_id],
          p_node->parent_nat, point);
      if (isfinite(point[0]) && isfinite(point[1]) && isfinite(point[2]))
        for (int co(0); co < 3; co++)
          p_node->coord[co] = point[co];
    }
  }

  // Bulk/artificial S/L cut remains parent controlled (same class as bulk
  // liquid BC). Beside coordinates (via BC_SUB above), refresh T and surface
  // slope Th (plus Vn/K) from the immediate parent. Exact coincident BC
  // records are copied; noncoincident recursive cut nodes use bilinear
  // interpolation on the parent S/L belt.
  BCSet *p_sl = &p_local->bc_set[p_local->BC_SUB_SOLID_LIQUID];
  const int parentSlBc =
      (p_parent->number_of_bc_sets > p_parent->BC_SUB_SOLID_LIQUID &&
       p_parent->bc_set[p_parent->BC_SUB_SOLID_LIQUID].bc_at_node != NULL)
          ? p_parent->BC_SUB_SOLID_LIQUID
          : p_parent->BC_SOLID_LIQUID;
  BCSet *p_parent_sl = &p_parent->bc_set[parentSlBc];
  for (int ib(0); ib < p_sl->number_of_recordsNode; ib++) {
    BCAtNode *p_bc = &p_sl->bc_at_node[ib];
    if (!p_bc->p_node->artificial_boundary)
      continue;
    const int parentRecord = p_bc->coarse_bc_node_id;
    if (parentRecord >= 0 &&
        parentRecord < p_parent_sl->number_of_recordsNode) {
      BCAtNode *p_src = &p_parent_sl->bc_at_node[parentRecord];
      p_bc->Vn = p_src->Vn;
      p_bc->K = p_src->K;
      p_bc->Th = p_src->Th;
      p_bc->p_node->Vn = p_src->Vn;
      const double Tsrc = p_src->p_node->attribute[TEMPERATURE];
      p_bc->p_node->attribute[TEMPERATURE] = Tsrc;
      p_bc->p_node->T = Tsrc;
      p_bc->interp_T_image = Tsrc;
      continue;
    }
    if (!parentIsTplSubmesh || p_parent_sl->bc_at_nodeStr == NULL)
      continue;
    const double parentPhi =
        (double)p_bc->p_node->local_i / (double)phiRefinement;
    const double parentPath = .5 * (double)p_bc->p_node->local_j;
    int ii[2], jj[2];
    double wi[2], wj[2];
    tplLinearGridStencil(parentPhi, p_parent->nTPL, 1, ii, wi);
    tplLinearGridStencil(parentPath, p_parent->subSlRefinement, 0, jj,
                         wj);
    double Th = 0., Vn = 0., K = 0., Tsum = 0.;
    for (int ai(0); ai < 2; ai++)
      for (int aj(0); aj < 2; aj++) {
        const double weight = wi[ai] * wj[aj];
        const BCAtNode *p_src =
            p_parent_sl->bc_at_nodeStr[ii[ai]][jj[aj]];
        Th += weight * p_src->Th;
        Vn += weight * p_src->Vn;
        K += weight * p_src->K;
        Tsum += weight * p_src->p_node->attribute[TEMPERATURE];
      }
    p_bc->Th = Th;
    p_bc->Vn = Vn;
    p_bc->K = K;
    p_bc->p_node->Vn = Vn;
    p_bc->p_node->attribute[TEMPERATURE] = Tsum;
    p_bc->p_node->T = Tsum;
    p_bc->interp_T_image = Tsum;
  }

  // LG outer ring (j=nS2-1) is bulk liquid BC; TPL end j=0 is shared. Both
  // stay parent-driven each pass (exact-copy / linear-sample / parent-element
  // like bulk T and BC_SUB geometry). Free-surface interior j=1..nS2-2 is
  // left for MoveLiquidInterface on every submesh level.
  BCSet *p_lg = &p_local->bc_set[p_local->BC_SUB_LIQUID_GAS];
  if (p_lg->bc_at_nodeStr != NULL && p_lg->nS2 >= 3) {
    const int parentLgBc =
        (p_parent->number_of_bc_sets > p_parent->BC_SUB_LIQUID_GAS &&
         p_parent->bc_set[p_parent->BC_SUB_LIQUID_GAS].bc_at_node != NULL &&
         p_parent->bc_set[p_parent->BC_SUB_LIQUID_GAS].bc_at_nodeStr != NULL)
            ? p_parent->BC_SUB_LIQUID_GAS
            : p_parent->BC_LIQUID_GAS;
    BCSet *p_parent_lg = &p_parent->bc_set[parentLgBc];
    const int active_cols = p_lg->nS1 - p_lg->structured_wrap_columns;
    for (int i(0); i < active_cols; i++) {
      for (int end(0); end < 2; end++) {
        const int j = end ? (p_lg->nS2 - 1) : 0;
        Node *p_node = p_lg->bc_at_nodeStr[i][j]->p_node;
        double point[3] = {0., 0., 0.};
        int have = 0;
        if (p_node->coarse_node_id >= 0 &&
            p_node->coarse_node_id < p_parent->number_of_nodes) {
          Node *p_src = &p_parent->node[p_node->coarse_node_id];
          for (int co(0); co < 3; co++)
            point[co] = p_src->coord[co];
          have = 1;
        } else if (parentIsTplSubmesh &&
                   p_parent_lg->bc_at_nodeStr != NULL &&
                   p_parent_lg->nS2 >= 3) {
          const double parentPhi =
              (double)i / (double)phiRefinement;
          const double parentPath =
              .5 * (double)j;
          int ii[2], jj[2];
          double wi[2], wj[2];
          tplLinearGridStencil(parentPhi, p_parent->nTPL, 1, ii, wi);
          tplLinearGridStencil(parentPath, p_parent_lg->nS2 - 1, 0, jj,
                               wj);
          for (int ai(0); ai < 2; ai++)
            for (int aj(0); aj < 2; aj++) {
              const double weight = wi[ai] * wj[aj];
              const Node *p_src =
                  p_parent_lg->bc_at_nodeStr[ii[ai]][jj[aj]]->p_node;
              for (int co(0); co < 3; co++)
                point[co] += weight * p_src->coord[co];
            }
          have = 1;
        } else if (p_node->parent_element_id >= 0 &&
                   p_node->parent_element_id <
                       p_parent->number_of_elements) {
          tplInterpolateElementPoint(
              &p_parent->element[p_node->parent_element_id],
              p_node->parent_nat, point);
          have = 1;
        }
        if (!have)
          continue;
        if (!isfinite(point[0]) || !isfinite(point[1]) ||
            !isfinite(point[2]))
          continue;
        for (int co(0); co < 3; co++)
          p_node->coord[co] = point[co];
      }
    }
  }
}

[[maybe_unused]] static void transferTplSubmeshTemperature(
    Mesh *p_local, const Mesh *p_parent) {
  // The complete initial child field comes from its immediate parent.  Later
  // thermal coupling uses the same parent on the artificial Dirichlet cut.
  for (int in(0); in < p_local->number_of_nodes; in++) {
    Node *p_node = &p_local->node[in];
    double temperature = 0.;
    if (p_node->coarse_node_id >= 0)
      temperature =
          p_parent->node[p_node->coarse_node_id].attribute[TEMPERATURE];
    else {
      const int parent = p_node->parent_element_id;
      temperature = tplInterpolateElementTemperatureSecondOrder(
          &p_parent->element[parent], p_node->parent_nat, 1., NULL, NULL);
    }
    p_node->T = temperature;
    p_node->attribute[TEMPERATURE] = temperature;
  }
  BCSet *p_artificial = &p_local->bc_set[p_local->BC_SUB];
  for (int ib(0); ib < p_artificial->number_of_recordsNode; ib++)
    p_artificial->bc_at_node[ib].interp_T_image =
        p_artificial->bc_at_node[ib].p_node->attribute[TEMPERATURE];
}

static void imposeTplArtificialTemperature(Mesh *p_local,
                                            const Mesh *p_parent,
                                            double timeFraction,
                                            const double *p_oldSolidGradient,
                                            const double *p_oldLiquidGradient) {
  BCSet *p_artificial = &p_local->bc_set[p_local->BC_SUB];
  for (int ib(0); ib < p_artificial->number_of_recordsNode; ib++) {
    BCAtNode *p_bc = &p_artificial->bc_at_node[ib];
    Node *p_node = p_bc->p_node;
    const int parent = p_node->parent_element_id;
    const double temperature = tplInterpolateElementTemperatureSecondOrder(
        &p_parent->element[parent], p_node->parent_nat, timeFraction,
        p_oldSolidGradient, p_oldLiquidGradient);
    p_bc->interp_T_image = temperature;
    p_node->attribute[TEMPERATURE] = temperature;
  }
}

[[maybe_unused]] static double solveTplSubmeshTemperature(
    Mesh *p_local, const Mesh *p_parent,
    LocalCalc *p_local_calc,
    GaussPointBFSet *p_gp_brick,
    GaussPointBFSet *p_gp_quad,
    double physicalDt,
    const double *p_oldSolidGradient,
    const double *p_oldLiquidGradient,
    int *p_thermalSubsteps) {
  const int nNodes = p_local->number_of_nodes;
  vector<double> capacity(nNodes, 0.);
  if (physicalDt <= 0.) {
    imposeTplArtificialTemperature(
        p_local, p_parent, 1.,
        p_oldSolidGradient, p_oldLiquidGradient);
    *p_thermalSubsteps = 0;
    return 0.;
  }

  vector<double> oldTemperature(nNodes, 0.);
  for (int in(0); in < nNodes; in++)
    oldTemperature[in] = p_local->node[in].attribute[TEMPERATURE];

  // Geometry and material coefficients are fixed during this local interval.
  // Assemble the same brick diffusion operator and lumped capacity used by
  // the explicit FE residual, then solve one backward-Euler step. This avoids
  // a CFL time step dictated by the 12-um recursive cells without changing
  // the heat equation or the physical duration advanced.
  vector<vector<int> > rowColumns(nNodes);
  for (int ie(0); ie < p_local->number_of_elements; ie++) {
    Element *p_element = &p_local->element[ie];
    for (int row(0); row < 8; row++) {
      const int rowNode = p_element->p_node[row]->node_num;
      for (int column(0); column < 8; column++)
        rowColumns[rowNode].push_back(
            p_element->p_node[column]->node_num);
    }
  }
  vector<int> operatorRow(nNodes + 1, 0);
  vector<int> operatorColumn;
  for (int row(0); row < nNodes; row++) {
    sort(rowColumns[row].begin(), rowColumns[row].end());
    rowColumns[row].erase(
        unique(rowColumns[row].begin(), rowColumns[row].end()),
        rowColumns[row].end());
    operatorRow[row] = (int)operatorColumn.size();
    operatorColumn.insert(operatorColumn.end(),
                          rowColumns[row].begin(), rowColumns[row].end());
  }
  operatorRow[nNodes] = (int)operatorColumn.size();
  vector<double> operatorValue(operatorColumn.size(), 0.);

  // Same linear diffusion operator as LocalGetResJac3dDin: one Gauss sweep
  // builds K_ij = cond * ∫ ∇Ni·∇Nj and the lumped capacity, instead of
  // eight unit-source residual probes per brick.
  for (int ie(0); ie < p_local->number_of_elements; ie++) {
    Element *p_element = &p_local->element[ie];
    const double conductivity = p_local_calc->cond[p_element->group_num];
    int localColumnIndex[8][8];
    for (int row(0); row < 8; row++) {
      const int rowNode = p_element->p_node[row]->node_num;
      for (int column(0); column < 8; column++) {
        const int columnNode = p_element->p_node[column]->node_num;
        localColumnIndex[row][column] =
            (int)(lower_bound(rowColumns[rowNode].begin(),
                              rowColumns[rowNode].end(), columnNode) -
                  rowColumns[rowNode].begin());
      }
    }
    for (int ig(0); ig < p_gp_brick->num_of_gp; ig++) {
      double jac[3][3] = {};
      for (int iglob(0); iglob < 3; iglob++)
        for (int iloc(0); iloc < 3; iloc++)
          for (int ibf(0); ibf < 8; ibf++)
            jac[iloc][iglob] +=
                p_gp_brick->gp[ig].dbf[ibf].coord[iloc] *
                p_element->p_node[ibf]->coord[iglob];
      const double determjac =
          -jac[0][2] * jac[1][1] * jac[2][0] +
           jac[0][1] * jac[1][2] * jac[2][0] +
           jac[0][2] * jac[1][0] * jac[2][1] -
           jac[0][0] * jac[1][2] * jac[2][1] -
           jac[0][1] * jac[1][0] * jac[2][2] +
           jac[0][0] * jac[1][1] * jac[2][2];
      double invjac[3][3];
      invjac[0][0] = (-jac[1][2] * jac[2][1] + jac[1][1] * jac[2][2]) / determjac;
      invjac[0][1] = ( jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]) / determjac;
      invjac[0][2] = (-jac[0][2] * jac[1][1] + jac[0][1] * jac[1][2]) / determjac;
      invjac[1][0] = ( jac[1][2] * jac[2][0] - jac[1][0] * jac[2][2]) / determjac;
      invjac[1][1] = (-jac[0][2] * jac[2][0] + jac[0][0] * jac[2][2]) / determjac;
      invjac[1][2] = ( jac[0][2] * jac[1][0] - jac[0][0] * jac[1][2]) / determjac;
      invjac[2][0] = (-jac[1][1] * jac[2][0] + jac[1][0] * jac[2][1]) / determjac;
      invjac[2][1] = ( jac[0][1] * jac[2][0] - jac[0][0] * jac[2][1]) / determjac;
      invjac[2][2] = (-jac[0][1] * jac[1][0] + jac[0][0] * jac[1][1]) / determjac;
      double dbfg[8][3];
      for (int i(0); i < 8; i++)
        for (int iglob(0); iglob < 3; iglob++)
          dbfg[i][iglob] =
              p_gp_brick->gp[ig].dbf[i].coord[0] * invjac[iglob][0] +
              p_gp_brick->gp[ig].dbf[i].coord[1] * invjac[iglob][1] +
              p_gp_brick->gp[ig].dbf[i].coord[2] * invjac[iglob][2];
      const double weight = fabs(determjac) * p_gp_brick->gp[ig].wo;
      for (int row(0); row < 8; row++) {
        const int rowNode = p_element->p_node[row]->node_num;
        capacity[rowNode] +=
            weight * p_gp_brick->gp[ig].bf[row] * p_gp_brick->gp[ig].bf[row];
        for (int column(0); column < 8; column++) {
          const double gradDot =
              dbfg[row][0] * dbfg[column][0] +
              dbfg[row][1] * dbfg[column][1] +
              dbfg[row][2] * dbfg[column][2];
          operatorValue[operatorRow[rowNode] + localColumnIndex[row][column]] +=
              -weight * conductivity * gradDot;
        }
      }
    }
  }
  const vector<double> intervalStartTemperature = oldTemperature;
  const int implicitSubsteps = 1;
  const double implicitDt = physicalDt / (double)implicitSubsteps;
  vector<double> surfaceResidual(nNodes, 0.);
  vector<double> rhs(nNodes, 0.);
  vector<double> diagonal(nNodes, 1.);
  vector<double> temperature(nNodes, 0.);
  vector<double> cgResidual(nNodes, 0.);
  vector<double> preconditioned(nNodes, 0.);
  vector<double> direction(nNodes, 0.);
  vector<double> matrixDirection(nNodes, 0.);
  int totalIterations = 0;

  for (int substep(0); substep < implicitSubsteps; substep++) {
    for (int in(0); in < nNodes; in++)
      oldTemperature[in] = p_local->node[in].attribute[TEMPERATURE];
    imposeTplArtificialTemperature(
        p_local, p_parent, (double)(substep + 1) / implicitSubsteps,
        p_oldSolidGradient, p_oldLiquidGradient);

    fill(surfaceResidual.begin(), surfaceResidual.end(), 0.);
    BCSet *physical[2] = {
        &p_local->bc_set[p_local->BC_SUB_SOLID_GAS],
        &p_local->bc_set[p_local->BC_SUB_LIQUID_GAS]};
    for (int is(0); is < 2; is++)
      for (int ib(0); ib < physical[is]->number_of_recordsFace; ib++) {
        BCAtFace *p_face = &physical[is]->bc_at_face[ib];
        p_local_calc->LocalGetResJac2dOutwall(p_face, p_gp_quad);
        for (int n(0); n < 4; n++)
          surfaceResidual[p_face->p_node[n]->node_num] +=
              p_local_calc->local_res[n];
      }
    BCSet *p_sl = &p_local->bc_set[p_local->BC_SUB_SOLID_LIQUID];
    for (int ib(0); ib < p_sl->number_of_recordsFace; ib++) {
      BCAtFace *p_face = &p_sl->bc_at_face[ib];
      p_local_calc->LocalGetResidual2dInterfEntalp(p_face, p_gp_quad);
      for (int n(0); n < 4; n++)
        surfaceResidual[p_face->p_node[n]->node_num] +=
            p_local_calc->local_res[n];
    }

    temperature = oldTemperature;
    fill(rhs.begin(), rhs.end(), 0.);
    fill(diagonal.begin(), diagonal.end(), 1.);
    fill(cgResidual.begin(), cgResidual.end(), 0.);
    fill(preconditioned.begin(), preconditioned.end(), 0.);
    fill(direction.begin(), direction.end(), 0.);
    for (int row(0); row < nNodes; row++) {
      Node *p_row = &p_local->node[row];
      if (p_row->artificial_boundary) {
        temperature[row] = p_row->attribute[TEMPERATURE];
        continue;
      }
      rhs[row] = capacity[row] * oldTemperature[row] +
                 implicitDt * surfaceResidual[row];
      diagonal[row] = capacity[row];
      for (int entry(operatorRow[row]); entry < operatorRow[row + 1]; entry++) {
        const int column = operatorColumn[entry];
        const double matrixValue = -implicitDt * operatorValue[entry];
        if (column == row)
          diagonal[row] += matrixValue;
        if (p_local->node[column].artificial_boundary)
          rhs[row] -= matrixValue *
                      p_local->node[column].attribute[TEMPERATURE];
      }
    }

    double rhsNorm = 0.;
    double rz = 0.;
    double residualNorm = 0.;
    for (int row(0); row < nNodes; row++) {
      if (p_local->node[row].artificial_boundary)
        continue;
      double matrixTemperature = capacity[row] * temperature[row];
      for (int entry(operatorRow[row]); entry < operatorRow[row + 1]; entry++) {
        const int column = operatorColumn[entry];
        if (!p_local->node[column].artificial_boundary)
          matrixTemperature -= implicitDt * operatorValue[entry] *
                               temperature[column];
      }
      cgResidual[row] = rhs[row] - matrixTemperature;
      preconditioned[row] = cgResidual[row] / diagonal[row];
      direction[row] = preconditioned[row];
      rz += cgResidual[row] * preconditioned[row];
      residualNorm += cgResidual[row] * cgResidual[row];
      rhsNorm += rhs[row] * rhs[row];
    }

    const double tolerance2 = 1.e-20 * rhsNorm + 1.e-60;
    const int maxIterations = 8 * nNodes;
    int iterations = 0;
    for (; iterations < maxIterations && residualNorm > tolerance2;
         iterations++) {
      double directionMatrixDirection = 0.;
      for (int row(0); row < nNodes; row++) {
        if (p_local->node[row].artificial_boundary)
          continue;
        matrixDirection[row] = capacity[row] * direction[row];
        for (int entry(operatorRow[row]); entry < operatorRow[row + 1]; entry++) {
          const int column = operatorColumn[entry];
          if (!p_local->node[column].artificial_boundary)
            matrixDirection[row] -= implicitDt * operatorValue[entry] *
                                    direction[column];
        }
        directionMatrixDirection +=
            direction[row] * matrixDirection[row];
      }
      if (!isfinite(directionMatrixDirection) ||
          directionMatrixDirection <= 0.) {
        *p_thermalSubsteps = -1;
        return NAN;
      }
      const double alpha = rz / directionMatrixDirection;
      double residual2 = 0.;
      for (int row(0); row < nNodes; row++) {
        if (p_local->node[row].artificial_boundary)
          continue;
        temperature[row] += alpha * direction[row];
        cgResidual[row] -= alpha * matrixDirection[row];
        residual2 += cgResidual[row] * cgResidual[row];
      }
      if (residual2 <= tolerance2) {
        iterations++;
        residualNorm = residual2;
        break;
      }
      double rzNew = 0.;
      for (int row(0); row < nNodes; row++) {
        if (p_local->node[row].artificial_boundary)
          continue;
        preconditioned[row] = cgResidual[row] / diagonal[row];
        rzNew += cgResidual[row] * preconditioned[row];
      }
      const double beta = rzNew / rz;
      for (int row(0); row < nNodes; row++)
        if (!p_local->node[row].artificial_boundary)
          direction[row] = preconditioned[row] + beta * direction[row];
      rz = rzNew;
      residualNorm = residual2;
    }
    totalIterations += iterations;
    for (int in(0); in < nNodes; in++)
      if (!p_local->node[in].artificial_boundary)
        p_local->node[in].attribute[TEMPERATURE] = temperature[in];
  }
  *p_thermalSubsteps = totalIterations;

  double norm = 0.;
  for (int in(0); in < nNodes; in++) {
    if (p_local->node[in].artificial_boundary)
      continue;
    const double change =
        p_local->node[in].attribute[TEMPERATURE] -
        intervalStartTemperature[in];
    norm += change * change;
  }
  return norm;
}

[[maybe_unused]] static void calcTplSubmeshGradients(
    Mesh *p_local, LocalCalc *p_local_calc,
    GaussPointBFSet *p_gp_corners) {
  for (int in(0); in < p_local->number_of_nodes; in++) {
    p_local->node[in].num_of_adj_elem_l = 0;
    p_local->node[in].num_of_adj_elem_s = 0;
    for (int co(0); co < 3; co++) {
      p_local->node[in].dTl[co] = 0.;
      p_local->node[in].dTs[co] = 0.;
    }
  }
  for (int ie(0); ie < p_local->number_of_elements; ie++)
    p_local_calc->LocalGetGrad3d(&p_local->element[ie], p_gp_corners);
  for (int in(0); in < p_local->number_of_nodes; in++) {
    if (p_local->node[in].num_of_adj_elem_l > 0)
      for (int co(0); co < 3; co++)
        p_local->node[in].dTl[co] /=
            p_local->node[in].num_of_adj_elem_l;
    if (p_local->node[in].num_of_adj_elem_s > 0)
      for (int co(0); co < 3; co++)
        p_local->node[in].dTs[co] /=
            p_local->node[in].num_of_adj_elem_s;
  }
}

static void tplCoonsPoint(const double bottom[3], const double top[3],
                          const double outer[3], const double inner[3],
                          const double corner00[3], const double corner10[3],
                          const double corner01[3], const double corner11[3],
                          double u, double v, double point[3]) {
  for (int co(0); co < 3; co++) {
    const double bilinear =
        (1. - u) * (1. - v) * corner00[co] +
        u * (1. - v) * corner10[co] +
        (1. - u) * v * corner01[co] + u * v * corner11[co];
    point[co] = (1. - v) * bottom[co] + v * top[co] +
                (1. - u) * outer[co] + u * inner[co] - bilinear;
  }
}

static void deformTplSubmeshInterior(Mesh *p_local,
                                     const Mesh *p_coarse) {
  (void)p_coarse;
  const int nPhi = p_local->nTPL;
  const int ni = p_local->subSlRefinement;
  const int ns = p_local->tplSubmeshSupportLayers;
  const int nr = p_local->tplSubmeshOuterElements;
  const int fixedSides =
      p_local->bc_set[p_local->BC_SUB_LIQUID_GAS]
          .structured_fixed_side_columns;
  const int iBegin = fixedSides ? 1 : 0;
  const int iEnd = fixedSides ? nPhi - 1 : nPhi;
  // The liquid block interface shared at the TPL is an explicit straight
  // column between its physical TPL endpoint and parent-driven lower cut.
  for (int i(iBegin); i < iEnd; i++) {
    Node *p_top = &p_local->node[tplInnerNodeId(p_local, i, 0, ns)];
    Node *p_bottom = &p_local->node[tplInnerNodeId(p_local, i, 0, 0)];
    for (int k(1); k < ns; k++) {
      const double v = (double)k / (double)ns;
      Node *p_node = &p_local->node[tplOuterNodeId(p_local, i, 0, k)];
      for (int co(0); co < 3; co++)
        p_node->coord[co] =
            (1. - v) * p_top->coord[co] + v * p_bottom->coord[co];
    }
  }

  for (int i(iBegin); i < iEnd; i++) {
    // Solid Coons patch between the physical TPL line and the immediate-
    // parent inner cut.
    for (int j(1); j < ni; j++)
      for (int k(1); k < ns; k++) {
        const double u = (double)j / (double)ni;
        const double v = (double)k / (double)ns;
        double point[3];
        tplCoonsPoint(
            p_local->node[tplInnerNodeId(p_local, i, j, ns)].coord,
            p_local->node[tplInnerNodeId(p_local, i, j, 2 * ns)].coord,
            p_local->node[tplInnerNodeId(p_local, i, 0, ns + k)].coord,
            p_local->node[tplInnerNodeId(p_local, i, ni, ns + k)].coord,
            p_local->node[tplInnerNodeId(p_local, i, 0, ns)].coord,
            p_local->node[tplInnerNodeId(p_local, i, ni, ns)].coord,
            p_local->node[tplInnerNodeId(p_local, i, 0, 2 * ns)].coord,
            p_local->node[tplInnerNodeId(p_local, i, ni, 2 * ns)].coord,
            u, v, point);
        Node *p_node =
            &p_local->node[tplInnerNodeId(p_local, i, j, ns + k)];
        for (int co(0); co < 3; co++)
          p_node->coord[co] = point[co];
      }

    // Inner liquid Coons patch.
    for (int j(1); j < ni; j++)
      for (int k(1); k < ns; k++) {
        const double u = (double)j / (double)ni;
        const double v = (double)k / (double)ns;
        double point[3];
        tplCoonsPoint(
            p_local->node[tplInnerNodeId(p_local, i, j, ns)].coord,
            p_local->node[tplInnerNodeId(p_local, i, j, 0)].coord,
            p_local->node[tplOuterNodeId(p_local, i, 0, k)].coord,
            p_local->node[tplInnerNodeId(p_local, i, ni, ns - k)].coord,
            p_local->node[tplInnerNodeId(p_local, i, 0, ns)].coord,
            p_local->node[tplInnerNodeId(p_local, i, ni, ns)].coord,
            p_local->node[tplInnerNodeId(p_local, i, 0, 0)].coord,
            p_local->node[tplInnerNodeId(p_local, i, ni, 0)].coord,
            u, v, point);
        Node *p_node =
            &p_local->node[tplInnerNodeId(p_local, i, j, ns - k)];
        for (int co(0); co < 3; co++)
          p_node->coord[co] = point[co];
      }

    // Outer liquid Coons patch.
    for (int j(1); j < nr; j++)
      for (int k(1); k < ns; k++) {
        const double u = (double)j / (double)nr;
        const double v = (double)k / (double)ns;
        double point[3];
        tplCoonsPoint(
            p_local->node[tplOuterNodeId(p_local, i, j, 0)].coord,
            p_local->node[tplOuterNodeId(p_local, i, j, ns)].coord,
            p_local->node[tplOuterNodeId(p_local, i, 0, k)].coord,
            p_local->node[tplOuterNodeId(p_local, i, nr, k)].coord,
            p_local->node[tplOuterNodeId(p_local, i, 0, 0)].coord,
            p_local->node[tplOuterNodeId(p_local, i, nr, 0)].coord,
            p_local->node[tplOuterNodeId(p_local, i, 0, ns)].coord,
            p_local->node[tplOuterNodeId(p_local, i, nr, ns)].coord,
            u, v, point);
        Node *p_node = &p_local->node[tplOuterNodeId(p_local, i, j, k)];
        for (int co(0); co < 3; co++)
          p_node->coord[co] = point[co];
      }
  }
}

static void moveTplSubmeshInterfaces(Mesh *p_local, const Mesh *p_coarse,
                                     double physicalDt,
                                     double geometryDt,
                                     double GAr, double GAf) {
  tplCopyPhysics(p_local, p_coarse);
  BCSet *p_tpl = &p_local->bc_set[p_local->BC_SUB_TPL];
  BCSet *p_lg = &p_local->bc_set[p_local->BC_SUB_LIQUID_GAS];
  BCSet *p_sl = &p_local->bc_set[p_local->BC_SUB_SOLID_LIQUID];
  BCSet *p_sg = &p_local->bc_set[p_local->BC_SUB_SOLID_GAS];
  if (geometryDt > 0.) {
    if (p_local->useVoronkovGA)
      p_tpl->MoveTPLvoronkovGA(geometryDt, GAr, GAf,
                               p_local->AniGaFlg, p_local);
    else
      p_tpl->MoveTPL(geometryDt, GAr, GAf,
                     p_local->AniGaFlg, p_local);
    p_lg->MoveLiquidInterface(geometryDt, p_local);
  }
  // Same order at the full-ring and every recursive level: one physical TPL
  // update from the three-node L/G endpoint tangent, L/G motion, tangential
  // S/L redistribution, physical S/L kinetic relaxation, structured S/G and
  // L/G redistribution, crystal pull, then interior Coons deformation. The
  // BC_SUB endpoints are parent controlled and are not moved by these
  // relaxation or kinetic routines.
  p_sl->calcRelaxSL();
  p_sl->kinMoveInterface(physicalDt, p_local);
  p_sg->calcRelaxStr();
  p_lg->calcRelaxStr();
  p_sg->pullSGstr(physicalDt * p_local->Vpull);
  deformTplSubmeshInterior(p_local, p_coarse);
}

static void restrictTplSubmeshToCoarse(Mesh *p_local, Mesh *p_parent) {
  // Restriction is geometry-only and surface-only. Bulk/artificial coincident
  // nodes remain parent controlled (P->C refresh). Temperature and its
  // gradients remain on the level where they were solved.
  for (int in(0); in < p_local->number_of_nodes; in++) {
    Node *p_local_node = &p_local->node[in];
    if (p_local_node->coarse_node_id < 0)
      continue;
    if (p_local_node->artificial_boundary)
      continue;
    Node *p_parent_node = &p_parent->node[p_local_node->coarse_node_id];
    for (int co(0); co < 3; co++)
      p_parent_node->coord[co] = p_local_node->coord[co];
    p_parent_node->Vn = p_local_node->Vn;
  }

  BCSet *p_local_sl =
      &p_local->bc_set[p_local->BC_SUB_SOLID_LIQUID];
  const int parentSlBc =
      (p_parent->number_of_bc_sets > p_parent->BC_SUB_SOLID_LIQUID &&
       p_parent->bc_set[p_parent->BC_SUB_SOLID_LIQUID].bc_at_node != NULL)
          ? p_parent->BC_SUB_SOLID_LIQUID
          : p_parent->BC_SOLID_LIQUID;
  BCSet *p_parent_sl = &p_parent->bc_set[parentSlBc];
  for (int ib(0); ib < p_local_sl->number_of_recordsNode; ib++) {
    if (p_local_sl->bc_at_node[ib].p_node->artificial_boundary)
      continue;
    const int parentRecord = p_local_sl->bc_at_node[ib].coarse_bc_node_id;
    if (parentRecord < 0 ||
        parentRecord >= p_parent_sl->number_of_recordsNode)
      continue;
    p_parent_sl->bc_at_node[parentRecord].Vn =
        p_local_sl->bc_at_node[ib].Vn;
    p_parent_sl->bc_at_node[parentRecord].K =
        p_local_sl->bc_at_node[ib].K;
    p_parent_sl->bc_at_node[parentRecord].Th =
        p_local_sl->bc_at_node[ib].Th;
  }
}

[[maybe_unused]] static void setTplSubmeshAnalyticalField(
    Mesh *p_local, double T0, double Gz,
    double Gr, double r0, double z0) {
  for (int in(0); in < p_local->number_of_nodes; in++) {
    Node *p_node = &p_local->node[in];
    const double radius = sqrt(p_node->coord[0] * p_node->coord[0] +
                               p_node->coord[1] * p_node->coord[1]);
    const double temperature =
        T0 + Gz * (p_node->coord[2] - z0) + Gr * (radius - r0);
    p_node->T = temperature;
    p_node->attribute[TEMPERATURE] = temperature;
    const double erx = (radius > 1.e-20) ? p_node->coord[0] / radius : 1.;
    const double ery = (radius > 1.e-20) ? p_node->coord[1] / radius : 0.;
    p_node->dTs[0] = Gr * erx;
    p_node->dTs[1] = Gr * ery;
    p_node->dTs[2] = Gz;
    p_node->dTl[0] = Gr * erx;
    p_node->dTl[1] = Gr * ery;
    p_node->dTl[2] = Gz;
  }
}

[[maybe_unused]] static int writeTplSubmeshState(
    const char *filename, const Mesh *p_local, int level) {
  const string temporary = string(filename) + ".tmp";
  ofstream file(temporary.c_str());
  if (!file) {
    cerr << "Failed to open TPL submesh state file: " << filename << "\n";
    return 0;
  }
  file << "FEEXP_TPL_SUBMESH_STATE_V2\n";
  file << "LEVEL " << level << "\n";
  file << "COUNTS " << p_local->number_of_nodes << " "
       << p_local->number_of_elements << " "
       << p_local->number_of_node_attributes << "\n";
  file << "LAYOUT " << p_local->nTPL << " "
       << p_local->subSlRefinement << " "
       << p_local->tplSubmeshSupportLayers << " "
       << p_local->tplSubmeshOuterElements << " "
       << p_local->bc_set[p_local->BC_SUB_TPL]
              .structured_fixed_side_columns
       << "\n";
  file << scientific << setprecision(17);
  file << "PHYSICS " << p_local->ZI << " " << p_local->dx << " "
       << p_local->Vpull << " " << p_local->gam << " "
       << p_local->rhog << " " << p_local->betaRough << " "
       << p_local->betaStep << " " << p_local->kineticA2DN << " "
       << p_local->kineticB2DN << " " << p_local->anisotropicKinFlg << " "
       << p_local->AniGaFlg << " " << p_local->useVoronkovGA << "\n";
  file << "VORONKOV " << p_local->voronkovTm << " "
       << p_local->voronkovQ << " " << p_local->voronkovLambdaSG << " "
       << p_local->voronkovAtomicDensity << " "
       << p_local->voronkovAlphaTransition << " "
       << p_local->voronkovSigmaSLFacet << " "
       << p_local->voronkovSigmaSLRough << " "
       << p_local->voronkovSigmaSLPrimeFacet << " "
       << p_local->voronkovSigmaSG << " "
       << p_local->voronkovSigmaMG << " "
       << p_local->voronkovThetaFacet << " "
       << p_local->voronkovThetaStepCorrection << " "
       << p_local->voronkovVelocityToCm << " "
       << p_local->voronkovKV << "\n";
  file << "NODES\n";
  for (int in(0); in < p_local->number_of_nodes; in++) {
    const Node *p_node = &p_local->node[in];
    file << in;
    for (int co(0); co < 3; co++) file << " " << p_node->coord[co];
    for (int a(0); a < p_local->number_of_node_attributes; a++)
      file << " " << p_node->attribute[a];
    file << " " << p_node->T;
    for (int co(0); co < 3; co++) file << " " << p_node->dTs[co];
    for (int co(0); co < 3; co++) file << " " << p_node->dTl[co];
    file << " " << p_node->Vn
         << " " << p_node->bulk_mark
         << " " << p_node->bound_mark
         << " " << p_node->artificial_boundary
         << " " << p_node->local_i
         << " " << p_node->local_j
         << " " << p_node->local_k
         << " " << p_node->local_block
         << " " << p_node->parent_element_id;
    for (int co(0); co < 3; co++) file << " " << p_node->parent_nat[co];
    file << " " << p_node->coarse_node_id << "\n";
  }
  file << "ELEMENTS\n";
  for (int ie(0); ie < p_local->number_of_elements; ie++) {
    const Element *p_element = &p_local->element[ie];
    file << ie << " " << p_element->group_num;
    for (int n(0); n < 8; n++) file << " " << p_element->assembly[n];
    file << "\n";
  }
  file << "END\n";
  file.close();
  if (!file) {
    cerr << "Incomplete TPL submesh state file: " << filename << "\n";
    return 0;
  }
  if (rename(temporary.c_str(), filename) != 0) {
    remove(filename);
    if (rename(temporary.c_str(), filename) != 0) {
      cerr << "Failed to publish TPL submesh state file: "
           << filename << "\n";
      return 0;
    }
  }
  return 1;
}

// ArtRidge reads this V1 trace directly. Every row is a solved node of the
// finest circular TPL mesh. Only a wrapped window around one facet direction
// is exported because ArtRidge models one ridge and its shoulders.
static int writeTplCircularBenchmark(
    const char *filename, Mesh *p_local, int recursiveLevels,
    int halfColumns, ostream &outw) {
  BCSet *p_tpl = &p_local->bc_set[p_local->BC_SUB_TPL];
  BCSet *p_sl = &p_local->bc_set[p_local->BC_SUB_SOLID_LIQUID];
  const int nPhi = p_tpl->nS1;
  const int count = 2 * halfColumns + 1;
  if (halfColumns < 1 || count > nPhi ||
      p_tpl->structured_fixed_side_columns != 0) {
    outw << "ridge benchmark: finest mesh is not a valid circular ring\n";
    return 0;
  }
  const int centerColumn = nPhi / 4;
  const int centerRow = halfColumns;
  vector<int> column(count, 0);
  for (int row(0); row < count; row++)
    column[row] = (centerColumn - halfColumns + row + nPhi) % nPhi;
  vector<double> arc(count, 0.);
  for (int row(centerRow + 1); row < count; row++) {
    Node *p0 = p_tpl->bc_at_nodeStr[column[row - 1]][1]->p_node;
    Node *p1 = p_tpl->bc_at_nodeStr[column[row]][1]->p_node;
    const double dx = p1->coord[0] - p0->coord[0];
    const double dy = p1->coord[1] - p0->coord[1];
    const double dz = p1->coord[2] - p0->coord[2];
    arc[row] = arc[row - 1] + sqrt(dx * dx + dy * dy + dz * dz);
  }
  for (int row(centerRow - 1); row >= 0; row--) {
    Node *p0 = p_tpl->bc_at_nodeStr[column[row]][1]->p_node;
    Node *p1 = p_tpl->bc_at_nodeStr[column[row + 1]][1]->p_node;
    const double dx = p1->coord[0] - p0->coord[0];
    const double dy = p1->coord[1] - p0->coord[1];
    const double dz = p1->coord[2] - p0->coord[2];
    arc[row] = arc[row + 1] - sqrt(dx * dx + dy * dy + dz * dz);
  }

  ofstream file(filename);
  if (!file) {
    outw << "ridge benchmark: failed to open " << filename << "\n";
    return 0;
  }
  file << "TITLE=\"FEexp circular-submesh to ArtRidge thermal bridge\"\n";
  file << "VARIABLES=\"X_m\" \"Y_m\" \"Z_minus_ZI_m\" \"R_m\" "
          "\"S_m\" \"T_minus_Tm_K\" \"GsR_K_m\" \"GsZ_K_m\" "
          "\"GlR_K_m\" \"GlZ_K_m\"\n";
  file << "DATASETAUXDATA FEEXP_BRIDGE=\"V1 " << recursiveLevels
       << " " << count << " " << centerRow << "\"\n";
  file << "ZONE T=\"finest solved circular TPL window\", I=" << count
       << ", F=POINT\n";
  file << scientific << setprecision(16);
  for (int row(0); row < count; row++) {
    const int i = column[row];
    BCAtNode *p_contact = p_tpl->bc_at_nodeStr[i][1];
    BCAtNode *p_sl_node =
        &p_sl->bc_at_node[p_contact->parent_bc_node_id];
    Node *p_node = p_sl_node->p_node;
    const double radius =
        sqrt(p_node->coord[0] * p_node->coord[0] +
             p_node->coord[1] * p_node->coord[1]);
    const double erx = p_node->coord[0] / radius;
    const double ery = p_node->coord[1] / radius;
    const double gsR = p_node->dTs[0] * erx + p_node->dTs[1] * ery;
    const double glR = p_node->dTl[0] * erx + p_node->dTl[1] * ery;
    file << p_node->coord[0] << " " << p_node->coord[1] << " "
         << p_node->coord[2] - p_local->ZI << " " << radius << " "
         << arc[i] << " " << p_node->attribute[TEMPERATURE] << " "
         << gsR << " " << p_node->dTs[2] << " "
         << glR << " " << p_node->dTl[2] << "\n";
  }
  file.close();
  if (!file) {
    outw << "ridge benchmark: incomplete output " << filename << "\n";
    return 0;
  }
  outw << "ridge benchmark | file=" << filename
       << " | recursive_levels=" << recursiveLevels
       << " | samples=" << count
       << " | range_mm=" << 1000. * arc.front()
       << ":" << 1000. * arc.back() << "\n";
  return 1;
}
