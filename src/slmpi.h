// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

#ifndef SLMPI_H
#define SLMPI_H

#include <iostream>
#include <vector>

namespace slmpi
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
    
    // Sum values from all ranks and distribute the result back to all ranks:
    void sum(std::vector<int>& data);
    void sum(std::vector<double>& data);
    
    // Broadcast from the broadcaster rank and receive in 'data' on all ranks.
    // Data vector must be preallocated to the correct size:
    void broadcast(int broadcaster, std::vector<int>& data);
    void broadcast(int broadcaster, std::vector<double>& data);

    // Returns on the gatherer a vector with an int value from all ranks:
    std::vector<int> gather(int gatherer, int value);

    // Gather fixed size fragments in the gatherer:
    void gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered);
    void gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered);
    // Gather variable size messages to the gatherer (fragment size only needed on gatherer):
    void gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered, std::vector<int>& fragsizes);
    void gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered, std::vector<int>& fragsizes);
    
    // Scatter fixed size messages.
    // Fragment vectors must be preallocated to the correct size:
    void scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment);
    void scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment);
    // Scatter variable size messages from the scatterer (fragment size only needed on scatterer).
    // Fragment vectors must be preallocated to the correct size:
    void scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment, std::vector<int>& fragsizes);
    void scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment, std::vector<int>& fragsizes);
    
    // Send + receive time for 'messagesize' doubles. Timings in ns are returned on rank 0:
    std::vector<double> ping(int messagesize, int verbosity = 1);
};

#endif

