# Ridge / FEexp3dSubR1

Research code for Czochralski crystal growth, ridge formation, and recursive TPL submeshes.

## Status

**Partial source snapshot. Not a standalone buildable release.**

This repository contains the five C++ source/header/build files available in the latest located handoff, copied byte-for-byte. The handoff itself also lacks required dependencies and physical inputs.

The original build was attempted on 26 September 2026 with GCC 13.3.0 on Linux. `sh build.sh` exited 1 at missing `src/facet.h`. No numerical simulation ran; no physical validation is claimed. See [build report](reports/Ridge_build_check_2026-09-26.md) and [compiler log](reports/build.log).

## Contents

- `3dFEexp.cpp` and `src/mesh/`: unchanged supplied C++ source/header files.
- `build.sh`: original build script; still needs the missing dependencies.
- `archive/Ridge_partial_sources_2026-09-26.zip`: ZIP of those five files only; NOT the original handoff archive.
- `reports/SHA256SUMS`: checksums for the published source/build files, ZIP, report and log.

Original project/agent documents and the original ZIP are withheld after automatic approval review rejected publication of detailed agent instructions and a local filesystem path in `todo`. This does not affect the source/build file contents.

## Next step

Supply the complete matching FEexp3dSubR1 source/input tree, including `param.txt`, the required headers/implementations, the axisymmetric initializer for the documented case, and the missing 2D driver. Then rebuild and run the bounded cases in separate case directories. Do not substitute unrelated historical dependencies or fabricate physical inputs.
