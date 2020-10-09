/*
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2020
 */

#ifndef _TOOLS_H
#define _TOOLS_H

#include <string>
#include <vector>

double StrToReal(std::string str);
int StrToInt(std::string str);

int FindIndex(double val, std::vector<double> &array);

#endif

