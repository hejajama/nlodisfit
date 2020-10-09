/*
 * Example code: how to read BK solution and extract dipole amplitude
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2019
 */

#include "amplitudelib.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <ctime>
#include <unistd.h>

using namespace std;



int main(int argc, char* argv[])
{
    // Case 1: Use fit where Y0=0
    
    string datafile="data/dipole-resumbk-hera-parent-0.00.dip";
    // Read data
    AmplitudeLib N(datafile);
    cout << "File 1: " << N.GetString() << endl;
    
    const double Ns =  1.0-std::exp(-0.5);
    cout << "# Q_s^2(Y=0) = " << N.SaturationScale(0, Ns) << " GeV^2" <<endl;
    cout << "# Q_s^2(Y=2) = " << N.SaturationScale(2, Ns) << " GeV^2" << endl;

    
    // Note: If I were to evaluate N many times at the same evolution
    // rapidity, I should do N.InitializeInterpolation(Y);
    
    
    cout << "# r [1/GeV]   N(r,Y=0)  N(r,Y=2)" << endl;
    double minr = N.MinR()*1.01; double maxr=N.MaxR()*0.99;
    for (double r=minr; r<maxr; r*=3)
    {
        cout << std::scientific << std::setprecision(9) << r << " " << N.DipoleAmplitude(r, 0) << " " << N.DipoleAmplitude(r, 2)  << endl;
    }
   
    
    cout << "====" << endl;
    // Case 2:
    // Here Y0=4.61
    // So dipole is frozen in the region 0<Y<Y0
    datafile="data/dipole-resumbk-hera-parent-4.61.dip";
    AmplitudeLib N2(datafile);
     cout << "File 2: " << N2.GetString() << endl;
    cout << "# Q_s^2(Y=0) = " << N2.SaturationScale(0, Ns) << " GeV^2" <<endl;
    cout << "# Q_s^2(Y=2) = " << N2.SaturationScale(2, Ns) << " GeV^2" << endl;
    cout << "# Q_s^2(Y=6) = " << N2.SaturationScale(6, Ns) << " GeV^2" << endl;
    
    cout << "# r [1/GeV]   N(r,Y=0)  N(r,Y=2)  N(r,Y=6)" << endl;
    for (double r=minr; r<maxr; r*=3)
    {
       cout << std::scientific << std::setprecision(9) << r << " " << N2.DipoleAmplitude(r, 0) << " " << N2.DipoleAmplitude(r, 2) << " " << N2.DipoleAmplitude(r,6) << endl;
    }
     
    return 0;
}
