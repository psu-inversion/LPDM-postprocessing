# LPDM-postprocessing
Utilities for postprocessing the output of the Uliasz (1994) LPDM into
forms more useful for inversions.

`carsurf_loop.py` replaces `carsurf_loop_v2.pro`, and
`run_carbounds.pro` copies additional information from the LPDM run
configuration file before delegating to `carsurf_loop.py` for the
heavy lifting.

## Initial setup

There is some setup to get everything working.

The most straightforward version is set up access to anaconda python, then execute:
```bash
conda create -n lpdm_post dateutil numpy scipy netcdf4 pyproj -c conda-forge
```
The next step is to change the lines in `run_carbounds.pro` saying
```bash
module switch python/3.4.0
. ~/python34
```
or
```bash
. ~/anaconda36
```
with directions to set up that access to anaconda, followed by
```bash
source activate lpdm_post
```

The other option is to set up your favorite python, netcdf, and compiler modules, then execute
```bash
python -m pip install --user dateutil numpy scipy netCDF4 pyproj
```
and changing the same lines in `run_carbounds.pro` as above to activate those same python, netcdf, and compiler modules.
