# nlodisfit
BK evolved dipole amplitude fitted to HERA data at NLO

This library provides easy access to the dipole-target scattering amplitude obtained by fitting the HERA reduced cross section data at NLO accuracy.

## Reference
When using this code, please cite
[1] G. Beuf, H. Hänninen, T. Lappi, H. Mäntysaari, Color Glass Condensate at next-to-leading order meets HERA data. e-Print: 2007.01645 [hep-ph]

Questions and comments to heikki.mantysaari@jyu.fi and henri.j.hanninen@jyu.fi

## Building
Requires
* Cmake
* GSL

How to compile:
```
 mkdir build
 cd build
 cmake ..
 make
```

This generates a library build/lib/libamplitude.a that you can link in your own program. The necessary header files can be found from the src directory. For example, see how ```src/dipole_amplitude.cpp``` is complied (see ```src/CMakeLists.txt```).

## Usage
See example ```src/dipole_amplitude.cpp```

* Add ```#include "amplitudelib.hpp"
* Read datafile: ```AmplitudeLib N(datafile);```
* Initialize interpolation, BK evolution over ```Y``` units of rapidity (optional, but makes evaluation of the dipole amplitude faster): ```N.InitializeInterpolation(Y);```
* Evaluate dipole: ```N.DipoleAmplitude(r,Y);```

Note that r has units [1/GeV].

## Datafiles
The different datafiles can be found as ```data/dipole-BKTYPE-FITDATA-COUPLING-Y0.dip```

The different BK equations (```BKTYPE```) are
* resumbk
* kcbk
* tbk

The coupling constants (```COUPLING```) are
* parent: Parent dipole
* Bal+SD: Balitsky+smallest

The initial rapidity where the BK evolution is started is ```Y0``` (in case of TBK, the number refers to the target rapidity eta_0).

These abbreviations are the same as used in the article referred above.

## Evolution rapidity
The initial condition is parametrized at ```Y=Y0``` in case of resumbk and kcbk, and at ```eta=eta0``` in case of tbk. Consequently, the dipole in the rapidity region ```[0,Y0]``` (or ```[0,eta0]```) is constant. 

When using tbk solution, the rapidity argument in the code refers to the evolution rapidity eta. In case of resumbk and kcbk, it refers to evolution rapidity Y. They are defined in Ref. [1] in Eqs. (16) and (26).

The datafiles contain the dipole amplitudes up to evolution rapidity 20.
