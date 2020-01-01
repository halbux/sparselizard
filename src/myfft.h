// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

#ifndef MYFFT_H
#define MYFFT_H

#include <iostream>
#include <cmath>
#include "densematrix.h"
#include <vector>
#include "harmonic.h"

namespace myfft
{
    // The fft is computed on every column of the input matrix.
    // Every input matrix row corresponds to a time evaluation 
    // of a flattened mym x myn matrix (rows concatenated).
    std::vector<std::vector<densematrix>> fft(densematrix input, int mym, int myn);
    
    // Remove the harmonics that are 'threshold' times smaller than the max(abs()) harmonic.
    void removeroundoffnoise(std::vector<std::vector<densematrix>>& input, double threshold = 1e-8);
    
    // The inverse fft takes a vector of harmonic values as input and
    // outputs a matrix where each row corresponds to a time value and 
    // each column to a data point of the densematrix inputs. 
    // input[harm][0] gives the values of harmonic 'harm'. 
    // If input[harm] is empty its value is assumed to be 0.
    densematrix inversefft(std::vector<std::vector<densematrix>>& input, int numtimevals, int mym, int myn);
        
    // Take a densematrix in the format given by 'inversefft', i.e.
    // with row i corresponding to time evaluation i. Output a 
    // densematrix in which row i corresponds to element i.
    densematrix toelementrowformat(densematrix timestepsinrows, int numberofelements);
    
    // Make sure notsame[i] has the same harmonic content for every i (fill with zeros if not).
    void sameharmonics(std::vector<std::vector<std::vector<densematrix>>>& notsame);
};

#endif
