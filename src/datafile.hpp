/*
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2020
 */

#ifndef _DATAFILE_H
#define _DATAFILE_H


#include <sstream>
#include <vector>

/**
 * Read BK equation solution from a given file.
 */
class DataFile
{
    public:
        /**
         * Read BK equation solution from the given file
         */
        DataFile(std::string fname);
        double MinR();
        double RMultiplier();
        int RPoints();
        
        double MaxY();
        double Y0();

        /**
         * Return data in a given vector and evolution y values
         *
         * @param n n[yind][rind]: dipole amplitude values
         * @param rapidities: evolution rapidity values
         */
        void GetData(std::vector< std::vector<double> > &n,
            std::vector<double> &rapidities);

    private:
        std::string filename;
        std::vector<std::vector <double> > data;
        std::vector<double> yvals;
        double minr;
        double r_multiplier;
        int rpoints;
        double y0;

};

#endif
