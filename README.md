# Ridge / FEexp3dSubR1

Research code for Czochralski crystal growth, ridge formation, and recursive TPL submeshes.

## Snapshot status

**This is a partial source handoff, not a standalone buildable release.**

Imported unchanged from `FEexp3dSubR1_Pro_handoff.zip` (14 files). The original archive is preserved under `archive/`.

Archive SHA256: `aaca63aa6a80bbfec9acd8eff2d7e0374d4196490b3c8593fd9b45e6f5bec7a0`.

The build was attempted on 26 September 2026 using GCC 13.3.0 on Linux. `sh build.sh` exited with status 1: missing `src/facet.h`. Additional sources and inputs are absent. No numerical simulation ran, and no physical validation is claimed.

See [the build report](reports/Ridge_build_check_2026-09-26.md) and [compiler log](reports/build.log) for evidence and the directly identified missing files.

## Contents

- `3dFEexp.cpp`, `src/mesh/`: supplied C++ sources.
- `AGENTS.md`, project/context documents, and `todo`: original handoff documentation; some documents describe different historical states.
- `build.sh`: original build commands, requiring missing dependencies.
- `archive/FEexp3dSubR1_Pro_handoff.zip`: byte-for-byte original archive.
- `reports/`: actual build outcome and integrity manifest.

## Next step

Add the complete matching FEexp3dSubR1 source/input tree, including `param.txt`, the axisymmetric initializer for the documented case, and the missing 2D driver. Then rebuild and run the documented bounded cases in separate working directories. Do not substitute unrelated historical dependencies or fabricate physical inputs.
