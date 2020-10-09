# nlodisfit
BK evolved dipole amplitude fitted to HERA data at NLO

This library provides easy access to the dipole-target scattering amplitude obtained by fitting the HERA reduced cross section data at NLO accuracy.

## Reference
When using this code, please cite
G. Beuf, H. Hänninen, T. Lappi, H. Mäntysaari, Color Glass Condensate at next-to-leading order meets HERA data. e-Print: 2007.01645 [hep-ph]

Questions and comments to heikki.mantysaari@jyu.fi and henri.j.hanninen@student.jyu.fi

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

This generates a library build/lib/libamplitude.a that you can link in your own program. The necessary header files can be found from the src directory. For example, see how ```src/dipole_amplitude``` is complied (see ```src/CMakeLists.txt```)

## Usage
See example ```src/dipole_amplitude.cpp```

* Read datafile: ```AmplitudeLib N(datafile);```
* Initialize interpolation, BK evolution over ```Y``` units of rapidity (evolution rapidity is ```Y+Y_{0,BK}```): ```N.InitializeInterpolation(Y);
* Evaluate dipole: ```N.DipoleAmplitude(r,Y);```

Note that r has units [1/GeV]

## Datafiles
The different datafils can be found from dirs ```data/dipole-BKTYPE-FITDATA-COUPLING-Y0.dip```
The different BK equations (```BKTYPE```) are
* resumbk
* kcbk
* tbk

The coupling constants are
* parent: Parent dipole
* Bal+SD: Balitsky+smallest

The initial rapidity where the BK evolution is started is ```Y0``` (in case of TBK, the number refers to the target rapidity eta_0)

These abbreviations are the same as used in the article referred above.
