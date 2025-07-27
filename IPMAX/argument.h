//
// Created by DBL-XQ on 2024/10/15.
//

#ifndef IPMAX_ARGUMENT_H
#define IPMAX_ARGUMENT_H

#include <limits.h>

class Argument
{
public:
    double _k = 0.1; // seed set size k. k<α*n
    int _R = 1000; // number of MC/TRIS.
    int _RSCC = 200; // number of SCC_DAG.
    int _theta = 5; // persistent time span.
    double _alpha = 0.2; // influence threshold percentage.
    int _W = 10; // size of temporal window.
//    int _Ts = 0; // start time interval
//    int _Te = INT32_MAX; // end time interval
    int _p = 100; // percentage of nodes in graph
    std::string _dataset = "../dataset/emailcore/";
    std::string _alg = "BA";
    // Algorithm. Default is BasicGreedy (BA,BA+MC),
    // ImproveGreedy (IM, BA+RIS), FillGreedy (FI, FI+RIS)
    // Combination (CO, CO+MC)

    Argument(int argc, char * argv[])
    {
        std::string param, value;
        for (int ind = 1; ind < argc; ind++)
        {
            if (argv[ind][0] != '-') {
                break;
            }
            std::stringstream sstr(argv[ind]);
            getline(sstr, param, '=');
            getline(sstr, value, '=');
            if (!param.compare("-k")) {
                _k = stod(value);
            }
            else if (!param.compare("-theta")) {
                _theta = stoi(value);
            }
            else if (!param.compare("-alpha")) {
                _alpha = stod(value);
            }
            else if (!param.compare("-R")) {
                _R = stoi(value);
            }
            else if (!param.compare("-RSCC")) {
                _RSCC = stoi(value);
            }
            else if (!param.compare("-W")) {
                _W = stoi(value);
            }
//            else if (!param.compare("-Ts")) {
//                _Ts = stoi(value);
//            }
//            else if (!param.compare("-Te")) {
//                _Te = stoi(value);
//            }
            else if (!param.compare("-p")) {
                _p = stoi(value);
            }
            else if (!param.compare("-dataset")) {
                _dataset = value;
            }
            else if (!param.compare("-alg")) {
                _alg = value;
            }
        }
    }
};

using TArgument = Argument;
using PArgument = std::shared_ptr<TArgument>;



#endif //IPMAX_ARGUMENT_H
