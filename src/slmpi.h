// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

#ifndef SLMPI_H
#define SLMPI_H

#include <iostream>
#include <vector>
#include "logs.h"

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
    void sum(int len, int* data);
    void sum(int len, long long int* data);
    void sum(int len, double* data);
    void sum(std::vector<int>& data);
    void sum(std::vector<double>& data);
    
    // Take the max of the values from all ranks and distribute the result back to all ranks:
    void max(int len, int* data);
    void max(int len, double* data);
    void max(std::vector<int>& data);
    void max(std::vector<double>& data);
    
    // Take the min of the values from all ranks and distribute the result back to all ranks:
    void min(int len, int* data);
    void min(int len, double* data);
    void min(std::vector<int>& data);
    void min(std::vector<double>& data);
    
    // Broadcast from the broadcaster rank and receive in 'data' on all ranks.
    // Data vector must be preallocated to the correct size:
    void broadcast(int broadcaster, std::vector<int>& data);
    void broadcast(int broadcaster, std::vector<double>& data);

    // Gather fixed size fragments in the gatherer:
    void gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered);
    void gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered);
    // Gather variable size messages to the gatherer (fragment sizes only needed on gatherer):
    void gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered, std::vector<int>& fragsizes);
    void gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered, std::vector<int>& fragsizes);
    
    // Gather fixed size fragments from all ranks in all ranks:
    void allgather(std::vector<int>& fragment, std::vector<int>& gathered);
    void allgather(std::vector<double>& fragment, std::vector<double>& gathered);
    // Gather variable size messages from all ranks in all ranks (fragment sizes needed on all ranks):
    void allgather(std::vector<int>& fragment, std::vector<int>& gathered, std::vector<int>& fragsizes);
    void allgather(std::vector<double>& fragment, std::vector<double>& gathered, std::vector<int>& fragsizes);
    
    // Scatter fixed size messages.
    // Fragment vectors must be preallocated to the correct size:
    void scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment);
    void scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment);
    // Scatter variable size messages from the scatterer (fragment sizes only needed on scatterer).
    // Fragment vectors must be preallocated to the correct size:
    void scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment, std::vector<int>& fragsizes);
    void scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment, std::vector<int>& fragsizes);
    
    // Exchange a fixed number of values with each unique target:
    void exchange(std::vector<int> targetranks, std::vector<int>& sendvalues, std::vector<int>& receivevalues);
    void exchange(std::vector<int> targetranks, std::vector<double>& sendvalues, std::vector<double>& receivevalues);
    // Exchange a message with each unique target (send and receive data must be preallocated to the correct size):
    void exchange(std::vector<int> targetranks, std::vector<std::vector<int>>& sends, std::vector<std::vector<int>>& receives);
    void exchange(std::vector<int> targetranks, std::vector<std::vector<double>>& sends, std::vector<std::vector<double>>& receives);
    void exchange(std::vector<int> targetranks, std::vector<int> sendlens, std::vector<int*> sendbuffers, std::vector<int> receivelens, std::vector<int*> receivebuffers);
    void exchange(std::vector<int> targetranks, std::vector<int> sendlens, std::vector<double*> sendbuffers, std::vector<int> receivelens, std::vector<double*> receivebuffers);
    
    // Send + receive time for 'messagesize' doubles. Timings in ns are returned on rank 0:
    std::vector<double> ping(int messagesize, int verbosity = 1);
};

#endif

