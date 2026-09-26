# Ridge archive: build check, 26 September 2026

## Source and target

- Latest matching archive found in Gmail: `FEexp3dSubR1_Pro_handoff.zip`, 119807 bytes, forwarded 24 September 2026; original message dated 20 September.
- SHA256: `aaca63aa6a80bbfec9acd8eff2d7e0374d4196490b3c8593fd9b45e6f5bec7a0`.
- ZIP integrity check passed; 14 entries extracted without changing source files.
- Candidate repository: https://github.com/OlegWeinstein/ridges . It is PUBLIC. Public upload was authorized by the user after this build check. This report accompanies the authorized source-only import.
- The accompanying PRO_DEBUG_PROMPT.md explicitly identifies this as a partial handoff, not a standalone package.

## Actual build result

Command: `sh build.sh` in the extracted source directory.

Exit status: **1**.

```text
3dFEexp.cpp:25:10: fatal error: src/facet.h: No such file or directory
   25 | #include "src/facet.h"
      |          ^~~~~~~~~~~~~
compilation terminated.
```

No executable was produced. No numerical simulation or physical validation ran. The second build target was not reached because the script stops on errors.

## Missing files directly identified from the available sources

- src/facet.h
- src/FE/basis_func.h
- src/FE/local_calc.h
- src/mesh/calculate.c
- src/mesh/init.c
- src/mesh/kin.c
- src/mesh/local_calc.c
- src/mesh/write_read.c
- src/research/params.cpp

This is the direct include list, not a complete transitive dependency inventory.

Also absent: `2d_axisymmetric_fe.cpp` referenced by build.sh, `param.txt` read by the main driver, and `input/axisymmetric_ic.dat` needed for the documented imported-initial-condition test. The driver has other startup paths, but no substitute inputs were fabricated.

## Remaining work

Provide the complete current FEexp3dSubR1 source/input tree. A Drive search located only an older FEexpDef10 folder from 2023, which was not substituted for current dependencies.

After the complete package is available, build it unchanged and run the documented bounded cases, subject to its current instructions:

```sh
OMP_NUM_THREADS=6 FEEXP_COARSE_N=9 FEEXP_END_STEPS=2 ./a.out
OMP_NUM_THREADS=6 FEEXP_COARSE_N=11 FEEXP_END_STEPS=2 ./a.out
```

Use separate case directories with the appropriate input files. These commands are proposed follow-up checks and HAVE NOT RUN.

The GitHub import preserves five source/build files byte-for-byte, a ZIP of those five files, and this build report. Original project/agent documents and the original handoff archive are withheld: automatic approval review rejected the todo file as containing detailed agent instructions and a local filesystem path outside its interpretation of the user authorization. The source-only ZIP is not the original archive. No code fix has been made or claimed.
