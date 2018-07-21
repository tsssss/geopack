# geopack and Tsyganenko models in Python
Author: Sheng Tian, Univ. of Minnesota, tianx138@umn.edu

This python `geopack` has integrated two modules originally written in Fortran: the `geopack` and the Tsyganenko models (T89, T96, T01, and T04). The Fortran `geopack` is available at https://ccmc.gsfc.nasa.gov/modelweb/magnetos/data-based/Geopack_2005.html and its DLM in IDL is available at http://ampere.jhuapl.edu/code/idl_geopack.html. As a crucial complement to `geopack`, the Tsyganenko models are available in Fortran at https://ccmc.gsfc.nasa.gov/models/modelinfo.php?model=Tsyganenko%20Magnetic%20Field.

A test is attached in `./python_geopack.md` to demonstrate that the Python `geopack` returns the same outputs as the Fortran and IDL counterparts. However, invisible to the user, several improvements have been implemented:
1. The latest IGRF coefficients are used, which cover the time range from 1900 to 2020. Years beyond this range are valid inputs and the corresponding IGRF coefficients will be extrapolated.
2. The IGRF coefficients in the Python `geopack` are smooth functions of the time (at milli-second accuray), whereas in the Fortran `geopack`, the coefficients are of daily resolution.
3. `igrf_gsm` is changed to a wrapper of `igrf_geo` plus the proper coordinate transforms. This is an obvious place for code-reusing, but there are many similar places in `goepack` and Tsyganenko models where there are repeated codes.
4. All `goto` statements in the Fortran `geopack` and Tsyganenko models are eliminated.

## Installation
The Python `geopack` is currently available at https://github.com/tsssss/geopack. Users are welcomed to download the source codes. I also plan to publish on PIP so it can be installed through command line: `$ pip install geopack`.

I've only tested the Python `geopack` on Mac OS in Python 3.6. Performance on other platform and other versions of Python is unclear.

## Notes on `geopack_08` and `T07d`
Strictly speaking, the Fortran `geopack` implemented here in Python is the `geopack_05`. A new version of `geopack_08` has been released, where the main change is to replace the widely used `GSM` coordinate with a newly defined `GSW` coordinate. Similarly, a new Tsyganenko `T07d` model has been released with a new algorithm. However, I decide to skip them for now. If people want the updates, please email me and let me know (tianx138@umn.edu).

## Package Interface
The Python `geopack` follows the Python way: function parameters are all input parameters and the outputs are returned. (This is *different* from the Fortran and IDL.)

* `ps = recalc(ut)`
* `bxgsm,bygsm,bzgsm = igrf_gsm(xgsm,ygsm,zgsm)`
Similarly `br,btheta,bphi = igrf_geo(r,theta,phi)`
* `bxgsm,bygsm,bzgsm = dip(xgsm,ygsm,zgsm)`
* `b1,b2,b3 = gsmgse(p1,p2,p3, flag)`
`flag` controls the direction of the transform: `flag>0` for gsm to gse; `flag<0` for gse to gsm.
The 6 available functions are `geomag`, `geigeo`, `magsm`, `gsmgse`, `smgsm`, `geogsm`. For definitions of the coordinates, please check Hapgood (1992).
* `b1,b2,b3 = sphcar(p1,p2,p3, flag)`
* `bx,by,bz = bspcar(theta,phi, br,btheta,bphi)`
* `br,theta,bphi = bcarsp(x,y,z, bx,by,bz)`
* `x1gsm,y1gsm,z1gsm = trace(x0gsm,y0gsm,z0gsm, dir, rlim, r0, par, exname, inname)`
`trace` will trace a given point `x0gsm,y0gsm,z0gsm` parallel (`dir=-1`) or anti-parallel (`dir=1`), until reaching a boundary.
`rlim` sets the radius of the outer boundary in Re, and `r0` sets the radius of the inner boundary in Re.
`inname` is a string sets either the dipole or IGRF model is used to provide the "internal" magnetic field of the earth.
`exname` is a string sets which Tsyganenko model is used to provide the "external" magnetic field, affected mainly by the solar wind.
`par` provides the solar wind info to the Tsyganenko models. For T89, `par` is an integer maps to the Kp index:
| `par` | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|-------|---|---|---|---|---|---|---|
| Kp    |0,0+|1-,1,1+|2-,2,2+|3-,3,3+|4-,4,4+|5-,5,5+|>6-|
* `bxgsm,bygsm,bzgsm = t96_mgnp(xn_pd, vel, xgsm,ygsm,zgsm)`
* `bxgms,bygsm,bzgsm = shuetal_mgnp(xn_pd, vel, bzimf, xgsm,ygsm,zgsm)`

## References
Hapgood, M. A. (1992). Space physics coordinate transformations: A user guide. Planetary and Space Science, 40(5), 711–717. http://doi.org/10.1016/0032-0633(92)90012-D


