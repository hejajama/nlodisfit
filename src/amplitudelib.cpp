/*
 * BK evolved dipole amplitude
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>,2020
 */

#include "amplitudelib.hpp"
#include "datafile.hpp"
#include "tools.hpp"
#include <string>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>
#include <string>
#include <sstream>


#include <algorithm>

using std::isinf;
using std::isnan;

inline double SQR(double x) { return x*x; }
#define LINEINFO __FILE__ << ":" << __LINE__

/*
 * Load data from a given file
 * Format is specified in file bk/README
 */
AmplitudeLib::AmplitudeLib(std::string datafile)
{
    // Read BK solution
    datafilename = datafile;
    DataFile data(datafile);
    data.GetData(n, yvals);

    // Initialize
    out_of_range_errors=true;
    interpolation_method = SPLINE_LINEAR;
    
    minr = data.MinR();
    rmultiplier = data.RMultiplier();
    rpoints = data.RPoints();

    tmprarray = new double[data.RPoints()];
    tmpnarray = new double[data.RPoints()];
    for (int i=0; i<rpoints; i++)
    {
        double tmpr = std::log(minr*std::pow(rmultiplier, i));
        lnrvals.push_back(tmpr);
        rvals.push_back(std::exp(lnrvals[i]));
        tmprarray[i] = std::exp(lnrvals[i]);
    }


    interpolator_Y=-1;  // if >=0, interpolator is initialized, must free
    // memory (delete tmprarray and tmpnarray at the end)


	std::stringstream ss;
    ss << "Data read from file " << datafile << ", minr: " << minr
        << " maxr: " << MaxR() << " rpoints: " << rpoints << " Y0=" << MinY() << ", max Y="
        << yvals[yvals.size()-1]
        << " Q_{s,0}^2(Y=Y0) = " << SaturationScale(0, 0.393469) << " GeV^2 [ N(r^2=2/Q_s^2, Y=0) = 0.3934]"
        << " (" << AMPLITUDELIB_VERSION << ")" ;
    info_string = ss.str();
    
    
}

/*
 * Constructor, dipole data is given in an array
 *
 */
AmplitudeLib::AmplitudeLib(std::vector< std::vector< double > > data, std::vector<double> yvals_, std::vector<double> rvals_)
{
    interpolation_method = SPLINE_LINEAR;
    out_of_range_errors = true;
    minr = rvals_[0];
    rmultiplier = rvals_[1] / rvals_[0];
    rpoints = rvals_.size();
    
    rvals = rvals_;
    n = data;
    yvals = yvals_;
    
    tmprarray = new double[rpoints];
    tmpnarray = new double[rpoints];
    for (int i=0; i<rpoints; i++)
    {
        double tmpr = std::log(minr*std::pow(rmultiplier, i));
        double tmpr2 = rvals[i];
        if ( abs(std::exp(tmpr)/tmpr2-1.0) > 0.01)
        {
            cerr << "WARNING: Dipole sizes does not form a logarithmic grid! tmpr " << tmpr << " r2 " << rvals[i] << " AmplitudeLib::AmplitudeLib()" << endl;
        }
        lnrvals.push_back(tmpr);
        tmprarray[i] = std::exp(lnrvals[i]);
    }
    
   
    interpolator_Y=-1;
    
    
    std::stringstream ss;
    ss << "#AmplitudeLib initialized , minr: " << minr
    << " maxr: " << MaxR() << " rpoints: " << rpoints << " Y0=" << MinY() << ", max Y=" << MaxY()
    << " Q_{s,0}^2(Y=0) = " << SaturationScale(0, 0.393469) << " GeV^2 [ N(r^2=2/Q_s^2, Y=0) = 0.3934]"
    << " (v. " << AMPLITUDELIB_VERSION << ")" ;
    info_string = ss.str();
    
    cout << ss.str() << endl;
    
    
}


/*
 * Release reserved memory
 */

AmplitudeLib::~AmplitudeLib()
{
    if (interpolator_Y>=0)
    {
        delete interpolator;
    }
    delete[] tmpnarray;
    delete[] tmprarray;
}

/*
 * Calculate amplitude interpolating data
 * Interpolate rapidity linearly and r using spline
 */
double AmplitudeLib::DipoleAmplitude(double r, double y)
{
    if (isnan(r) or isinf(r))
    {
        cerr << "r=" << r << " at AmplitudeLib::N(): " << LINEINFO << endl;
        exit(1);
    }
    
    if (r < MinR() or r > MaxR() )
    {
        if (out_of_range_errors)
            cerr << "r must be between limits [" << MinR() << ", " << MaxR() << "]"
                << " asked r=" << r << " " << LINEINFO
                << endl;
        if (r<MinR()) r=MinR()*1.000001; else if (r>MaxR()) r=MaxR()*0.999999;
    }
    
    // If Y0>0, the dipole is frozen in [0,Y0]
    if (y>=0 and y<MinY())
        y=MinY();
    
    
    if (y<MinY() or y>MaxY() )
    {
        //if (out_of_range_errors)
            cerr << "y must be between limits [" << 0 << ", "
                << yvals[yvals.size()-1] << "], asked y=" << y << " datafile " << datafilename << " "
                << LINEINFO << endl;
            exit(1);
    }

    
    // Use already initialized interpolator
    if (InterpolatorInitialized(y))
    {
        double result=0;

        // Can't interpolate (too large dipole), return 1.0
        if (r >= maxr_interpolate and maxr_interpolate>0)
        {
            return 1.0;
        }

        result = interpolator->Evaluate(r);
        
        if (result>1) return 1;              // limit N(r)<=1
        if (result<0) return 0;              // limit N >= 0
        return result;

    }
    
    int rind = FindIndex(r, rvals);
    int yind = FindIndex(y, yvals);
    
    // Check interpolation method, use bspline if needed
    if (interpolation_method == LINEAR_LINEAR)
    {
        int rind2 = rind+1;
        int yind2 = yind+1;
        if (rind2 >= rvals.size()) return 1.0;
        if (yind2 >= yvals.size()) return n[yind][rind]; // this should never happen
        double r1_ = rvals[rind];
        double r2_ = rvals[rind2];
        double y1_ = yvals[yind];
        double y2_ = yvals[yind2];
        double bilin =  (1.0/( (r2_ - r1_)*(y2_ - y1_) ))*( (n[yind][rind])*(r2_ - r)*(y2_ - y) + (n[yind][rind2])*(r - r1_)*(y2_ - y) + (n[yind2][rind])*(r2_ - r)*(y - y1_) + (n[yind2][rind2])*(r - r1_)*(y - y1_) );
        
        return bilin;
        
    }

    // Initialize new interpolator and use it
    // Note: This is relatively slow. User should
    // usually InitializeInterpolation and then
    // multiple evaluations are fast at the same Y
    int interpolation_points = INTERPOLATION_POINTS;

    int interpolation_start, interpolation_end;
    if (rind - interpolation_points/2 < 0)
    {
        interpolation_start=0;
        interpolation_end=interpolation_points;
    }
    else if (rind + interpolation_points/2 > RPoints()-1 )
    {
        interpolation_end = RPoints()-1;
        interpolation_start = RPoints()-interpolation_points-3;
    }
    else
    {
        interpolation_start = rind - interpolation_points/2;
        interpolation_end = rind + interpolation_points/2;
    }

    int interpo_points = interpolation_end - interpolation_start+1;
    double *tmparray = new double[interpo_points];
    double *tmpxarray = new double[interpo_points];
    for (int i=interpolation_start; i<= interpolation_end; i++)
    {
        tmpxarray[i-interpolation_start]=rvals[i];

        tmparray[i-interpolation_start] = n[yind][i];

        // Interpolate in y if possible
        if (yind < static_cast<int>(yvals.size()-1) )
        {
            tmparray[i-interpolation_start]=n[yind][i] 
                + (y - yvals[yind]) * (n[yind+1][i] - n[yind][i])
                / (yvals[yind+1]-yvals[yind]);
        }
    }
    
    Interpolator interp(tmpxarray, tmparray, interpo_points);
    interp.Initialize();
    double result=0;

    result = interp.Evaluate(r);

    
    delete[] tmparray;
    delete[] tmpxarray;

    return result;

}

bool AmplitudeLib::InterpolatorInitialized(double Y)
{
    // Check Y=0 separately
    if (std::abs(Y) < 0.001 and std::abs(interpolator_Y) < 0.001)
        return true;
    
    if (interpolator_Y >= 0 and std::abs(Y - interpolator_Y)/std::min(Y, interpolator_Y) < 0.001)
        return true;
    else
        return false;
}


/*
 * Calculate saturation scale, definition is
 * N(r_s, y) = Ns = (for example) 0.5
 */
struct SatscaleSolverHelper
{
    double Y;
    double Ns;
    AmplitudeLib* N;
};
double SatscaleSolverHelperf(double r, void* p)
{
    SatscaleSolverHelper* par = (SatscaleSolverHelper*)p;
    return par->N->DipoleAmplitude(r, par->Y) - par->Ns;
}


double AmplitudeLib::SaturationScale(double Y, double Ns)
{
    
    SatscaleSolverHelper helper;
    helper.Y=Y; helper.Ns=Ns; helper.N=this;
    const int MAX_ITER = 1000;
    const double ROOTFINDACCURACY = 0.00001;
    gsl_function f;
    f.params = &helper;
        
    f.function = &SatscaleSolverHelperf;

    const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
        
    gsl_root_fsolver_set(s, &f, MinR()*1.0001, MaxR()*0.999);
    int iter=0; int status; double min,max;
    do
    {
        iter++;
        gsl_root_fsolver_iterate(s);
        min = gsl_root_fsolver_x_lower(s);
        max = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(min, max, 0, ROOTFINDACCURACY);    
    } while (status == GSL_CONTINUE and iter < MAX_ITER);

    if (iter>=MAX_ITER)
        cerr << "Solving failed at Y=" << Y << " " << LINEINFO << endl;


    double sat = gsl_root_fsolver_root(s);

    gsl_root_fsolver_free(s);

    return 2.0/(sat*sat);
    
}




/*
 * Initializes dipole interpolation at given rapidity
 */
void AmplitudeLib::InitializeInterpolation(double Y)
{
    // Between 0<Y<MinY() use MinY
    if (Y>=0 and Y<=MinY()) Y=MinY();
    
	if (Y < 0)
	{
		cerr << "Asked to initialize interpolator at negative Y=" << Y
             << ", I don't know what to do, panicking... " << LINEINFO << endl;
		exit(1);
	}

	if (Y>MaxY())
	{
		cerr << "Asked to initialize interpolator with too large y=" << Y <<", maxy=" << MaxY() << endl;
		exit(1);
	}
	
    if (InterpolatorInitialized(Y)) return;    // Already done it
    
    // Remove old interpolator if constructed
    if (interpolator_Y>=0)
    {
        interpolator_Y=-1;
        delete interpolator;
    }
    
    for (int i=0; i<rpoints; i++)
    {
        double tmpr = tmprarray[i];
        if (i==0) tmpr*=1.0001; if (i==rpoints-1) tmpr*=0.9999;
        tmpnarray[i] = DipoleAmplitude(tmpr, Y);
    }
    interpolator = new Interpolator(tmprarray, tmpnarray, rpoints);

    interpolator->Initialize();
    interpolator_Y = Y;
    
    maxr_interpolate=-1;
    int iter=0;
    
    // Find rmax s.t. N(r) \approx 1 at r>rmax
    double step=2; double prevr=0.1;
    const int MAXITER=40;
    for (double r=0.01; r<=MaxR(); r+=step)
    {
        if (DipoleAmplitude(r,Y)>=0.99999)		//
        {
            if (step<1e-2)
            {
                maxr_interpolate = r;
                break;
            }
            step /= 1.5;
            r=prevr;
        }
        prevr=r;
        iter++;
        if (iter > MAXITER)
        {
            maxr_interpolate=-1;
            break;  // Didn't find, dont force any upper limit
        }
    }
}



/*
 * Initializes interpolator and returns it
 */
Interpolator* AmplitudeLib::ConstructInterpolator(double Y)
{
    std::vector<double> tmpnvals;
    for (int i=0; i<rpoints; i++)
    {
        tmpnvals.push_back(DipoleAmplitude(rvals[i], Y));
    }
    Interpolator* inter = new Interpolator(rvals, tmpnvals);
    inter->Initialize();

    return inter;
        
}


int AmplitudeLib::RPoints()
{
    return rpoints;
}
double AmplitudeLib::MinR()
{
    return minr;
}

double AmplitudeLib::MaxR()
{
    return minr*std::pow(rmultiplier, rpoints-1);
}

int AmplitudeLib::YPoints()
{
    return yvals.size();
}

double AmplitudeLib::MaxY()
{
    return yvals[yvals.size()-1];
}

double AmplitudeLib::MinY()
{
    return yvals[0];
}


bool AmplitudeLib::SetOutOfRangeErrors(bool er)
{
    bool outofrange = out_of_range_errors;
    out_of_range_errors=er;
    return outofrange;    
}


std::string AmplitudeLib::GetString()
{
	return info_string;
}

std::string AmplitudeLib::Version()
{
	std::stringstream s;
	s << AMPLITUDELIB_VERSION << " (build " <<  __DATE__ << " " << __TIME__ << ")";
	return s.str();
}


/*
 * Return ith rapidity value
 */
double AmplitudeLib::YValue(int yind)
{
    if (yind >= YPoints() or yind<0)
    {
        cerr << "Asked rapidity of index " << yind <<", but the maximum index is " << YPoints()-1 << " " << LINEINFO << endl;
        exit(1);
    }
    return yvals[yind];
}
