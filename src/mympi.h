// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

#ifndef MYMPI_H
#define MYMPI_H

#include <iostream>
#include <vector>

namespace mympi
{
    bool isavailable(void);
    
    void errornompi(void);

    void initialize(void);
    void finalize(void);

    int getrank(void);
    int count(void);
    
    void barrier(void);

    void send(int destination, int tag, std::vector<int>& data);
    void send(int destination, int tag, std::vector<double>& data);
    // Data vector must be preallocated to the correct size:
    void receive(int source, int tag, std::vector<int>& data);
    void receive(int source, int tag, std::vector<double>& data);
    
    // Broadcast from the broadcaster rank and receive in 'data' on all ranks:
    void broadcast(int broadcaster, std::vector<int>& data);
    void broadcast(int broadcaster, std::vector<double>& data);

    // Gather fixed size fragments in the gatherer:
    void gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered);
    void gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered);
    // Gather variable size messages to the gatherer (fragment size only needed on gatherer):
    void gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered, std::vector<int>& fragsizes);
    void gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered, std::vector<int>& fragsizes);
    
    // Scatter fixed size messages:
    void scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment);
    void scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment);
    // Scatter variable size messages from the scatterer (fragment size only needed on scatterer).
    // Fragment vectors must be preallocated to the correct size.
    void scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment, std::vector<int>& fragsizes);
    void scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment, std::vector<int>& fragsizes);
    
    // Send + receive time for 'messagesize' doubles. Timings in ns are returned on rank 0:
    std::vector<double> ping(int messagesize = 1e6, int verbosity = 1);
};

#endif

