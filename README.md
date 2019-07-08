# LPDM-postprocessing
Utilities for postprocessing the output of the Uliasz (1994) LPDM into
forms more useful for inversions.

`carsurf_loop.py` replaces `carsurf_loop_v2.pro`, and
`run_carbounds.pro` copies additional information from the LPDM run
configuration file before delegating to `carsurf_loop.py` for the
heavy lifting.
