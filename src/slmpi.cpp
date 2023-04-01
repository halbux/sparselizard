#include "slmpi.h"
#include "wallclock.h"


void slmpi::errornompi(void)
{
    logs log;
    log.msg() << "Error in 'slmpi' namespace: MPI is not available" << std::endl;
    log.error();
}


#ifndef HAVE_MPI
bool slmpi::isavailable(void) { return false; }
void slmpi::initialize(void) {}
void slmpi::finalize(void) {}
int slmpi::getrank(void) { return 0; }
int slmpi::count(void) { return 1; }
void slmpi::barrier(void) {}
void slmpi::send(int destination, int tag, std::vector<int>& data) { errornompi(); }
void slmpi::send(int destination, int tag, std::vector<double>& data) { errornompi(); }
void slmpi::receive(int source, int tag, std::vector<int>& data) { errornompi(); }
void slmpi::receive(int source, int tag, std::vector<double>& data) { errornompi(); }
void slmpi::sum(int len, int* data) {}
void slmpi::sum(int len, long long int* data) {}
void slmpi::sum(int len, double* data) {}
void slmpi::sum(std::vector<int>& data) {}
void slmpi::sum(std::vector<double>& data) {}
void slmpi::max(int len, int* data) {}
void slmpi::max(int len, double* data) {}
void slmpi::max(std::vector<int>& data) {}
void slmpi::max(std::vector<double>& data) {}
void slmpi::min(int len, int* data) {}
void slmpi::min(int len, double* data) {}
void slmpi::min(std::vector<int>& data) {}
void slmpi::min(std::vector<double>& data) {}
void slmpi::broadcast(int broadcaster, std::vector<int>& data) { errornompi(); }
void slmpi::broadcast(int broadcaster, std::vector<double>& data) { errornompi(); }
void slmpi::gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered) { errornompi(); }
void slmpi::gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered) { errornompi(); }
void slmpi::gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered, std::vector<int>& fragsizes) { errornompi(); }
void slmpi::gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered, std::vector<int>& fragsizes) { errornompi(); }
void slmpi::allgather(std::vector<int>& fragment, std::vector<int>& gathered) { errornompi(); }
void slmpi::allgather(std::vector<double>& fragment, std::vector<double>& gathered) { errornompi(); }
void slmpi::allgather(std::vector<int>& fragment, std::vector<int>& gathered, std::vector<int>& fragsizes) { errornompi(); }
void slmpi::allgather(std::vector<double>& fragment, std::vector<double>& gathered, std::vector<int>& fragsizes) { errornompi(); }
void slmpi::scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment) { errornompi(); }
void slmpi::scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment) { errornompi(); }
void slmpi::scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment, std::vector<int>& fragsizes) { errornompi(); }
void slmpi::scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment, std::vector<int>& fragsizes) { errornompi(); }
void slmpi::exchange(std::vector<int> targetranks, std::vector<int>& sendvalues, std::vector<int>& receivevalues) { errornompi(); }
void slmpi::exchange(std::vector<int> targetranks, std::vector<double>& sendvalues, std::vector<double>& receivevalues) { errornompi(); }
void slmpi::exchange(std::vector<int> targetranks, std::vector<std::vector<int>>& sends, std::vector<std::vector<int>>& receives) { errornompi(); }
void slmpi::exchange(std::vector<int> targetranks, std::vector<std::vector<double>>& sends, std::vector<std::vector<double>>& receives) { errornompi(); }
void slmpi::exchange(std::vector<int> targetranks, std::vector<int> sendlens, std::vector<int*> sendbuffers, std::vector<int> receivelens, std::vector<int*> receivebuffers) { errornompi(); }
void slmpi::exchange(std::vector<int> targetranks, std::vector<int> sendlens, std::vector<double*> sendbuffers, std::vector<int> receivelens, std::vector<double*> receivebuffers) { errornompi(); }
std::vector<double> slmpi::ping(int messagesize, int verbosity) { errornompi(); throw std::runtime_error(""); }
#endif


#ifdef HAVE_MPI

#include "mpi.h"

bool slmpi::isavailable(void) { return true; }

void slmpi::initialize(void)
{
    MPI_Init(NULL, NULL);
}

void slmpi::finalize(void)
{
    MPI_Finalize();
}

int slmpi::getrank(void)
{
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    return world_rank;
}

int slmpi::count(void)
{
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    return world_size;
}

void slmpi::barrier(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
}


void slmpi::send(int destination, int tag, std::vector<int>& data)
{
    MPI_Send(data.data(), data.size(), MPI_INT, destination, tag, MPI_COMM_WORLD);
}

void slmpi::send(int destination, int tag, std::vector<double>& data)
{
    MPI_Send(data.data(), data.size(), MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
}


void slmpi::receive(int source, int tag, std::vector<int>& data)
{
    MPI_Recv(data.data(), data.size(), MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void slmpi::receive(int source, int tag, std::vector<double>& data)
{
    MPI_Recv(data.data(), data.size(), MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


void slmpi::sum(int len, int* data)
{
    MPI_Allreduce(MPI_IN_PLACE, data, len, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void slmpi::sum(int len, long long int* data)
{
    MPI_Allreduce(MPI_IN_PLACE, data, len, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
}

void slmpi::sum(int len, double* data)
{
    MPI_Allreduce(MPI_IN_PLACE, data, len, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void slmpi::sum(std::vector<int>& data)
{
    MPI_Allreduce(MPI_IN_PLACE, data.data(), data.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void slmpi::sum(std::vector<double>& data)
{
    MPI_Allreduce(MPI_IN_PLACE, data.data(), data.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


void slmpi::max(int len, int* data)
{
    MPI_Allreduce(MPI_IN_PLACE, data, len, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
}

void slmpi::max(int len, double* data)
{
    MPI_Allreduce(MPI_IN_PLACE, data, len, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

void slmpi::max(std::vector<int>& data)
{
    MPI_Allreduce(MPI_IN_PLACE, data.data(), data.size(), MPI_INT, MPI_MAX, MPI_COMM_WORLD);
}

void slmpi::max(std::vector<double>& data)
{
    MPI_Allreduce(MPI_IN_PLACE, data.data(), data.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}


void slmpi::min(int len, int* data)
{
    MPI_Allreduce(MPI_IN_PLACE, data, len, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
}

void slmpi::min(int len, double* data)
{
    MPI_Allreduce(MPI_IN_PLACE, data, len, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}

void slmpi::min(std::vector<int>& data)
{
    MPI_Allreduce(MPI_IN_PLACE, data.data(), data.size(), MPI_INT, MPI_MIN, MPI_COMM_WORLD);
}

void slmpi::min(std::vector<double>& data)
{
    MPI_Allreduce(MPI_IN_PLACE, data.data(), data.size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}


void slmpi::broadcast(int broadcaster, std::vector<int>& data)
{
    MPI_Bcast(data.data(), data.size(), MPI_INT, broadcaster, MPI_COMM_WORLD);
}

void slmpi::broadcast(int broadcaster, std::vector<double>& data)
{
    MPI_Bcast(data.data(), data.size(), MPI_DOUBLE, broadcaster, MPI_COMM_WORLD);
}


void slmpi::gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered)
{
    gathered = {};
    
    if (getrank() == gatherer)
        gathered.resize(count()*fragment.size());

    MPI_Gather(fragment.data(), fragment.size(), MPI_INT, gathered.data(), fragment.size(), MPI_INT, gatherer, MPI_COMM_WORLD); 
}

void slmpi::gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered)
{
    gathered = {};
    
    if (getrank() == gatherer)
        gathered.resize(count()*fragment.size());

    MPI_Gather(fragment.data(), fragment.size(), MPI_DOUBLE, gathered.data(), fragment.size(), MPI_DOUBLE, gatherer, MPI_COMM_WORLD); 
}


void slmpi::gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered, std::vector<int>& fragsizes)
{
    gathered = {};
    
    int totlen = 0;
    std::vector<int> shifts(fragsizes.size());
    for (int i = 0; i < fragsizes.size(); i++)
    {
        shifts[i] = totlen;
        totlen += fragsizes[i];
    }
    if (getrank() == gatherer)
        gathered.resize(totlen);

    MPI_Gatherv(fragment.data(), fragment.size(), MPI_INT, gathered.data(), fragsizes.data(), shifts.data(), MPI_INT, gatherer, MPI_COMM_WORLD); 
}

void slmpi::gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered, std::vector<int>& fragsizes)
{
    gathered = {};
    
    int totlen = 0;
    std::vector<int> shifts(fragsizes.size());
    for (int i = 0; i < fragsizes.size(); i++)
    {
        shifts[i] = totlen;
        totlen += fragsizes[i];
    }
    if (getrank() == gatherer)
        gathered.resize(totlen);

    MPI_Gatherv(fragment.data(), fragment.size(), MPI_DOUBLE, gathered.data(), fragsizes.data(), shifts.data(), MPI_DOUBLE, gatherer, MPI_COMM_WORLD); 
}


void slmpi::allgather(std::vector<int>& fragment, std::vector<int>& gathered)
{
    gathered.resize(count()*fragment.size());

    MPI_Allgather(fragment.data(), fragment.size(), MPI_INT, gathered.data(), fragment.size(), MPI_INT, MPI_COMM_WORLD);
}

void slmpi::allgather(std::vector<double>& fragment, std::vector<double>& gathered)
{
    gathered.resize(count()*fragment.size());

    MPI_Allgather(fragment.data(), fragment.size(), MPI_DOUBLE, gathered.data(), fragment.size(), MPI_DOUBLE, MPI_COMM_WORLD);
}


void slmpi::allgather(std::vector<int>& fragment, std::vector<int>& gathered, std::vector<int>& fragsizes)
{
    std::vector<int> shifts(fragsizes.size(), 0);
    for (int i = 1; i < fragsizes.size(); i++)
        shifts[i] = shifts[i-1]+fragsizes[i-1];
        
    gathered.resize(shifts[fragsizes.size()-1]+fragsizes[fragsizes.size()-1]);

    MPI_Allgatherv(fragment.data(), fragment.size(), MPI_INT, gathered.data(), fragsizes.data(), shifts.data(), MPI_INT, MPI_COMM_WORLD);
}

void slmpi::allgather(std::vector<double>& fragment, std::vector<double>& gathered, std::vector<int>& fragsizes)
{
    std::vector<int> shifts(fragsizes.size(), 0);
    for (int i = 1; i < fragsizes.size(); i++)
        shifts[i] = shifts[i-1]+fragsizes[i-1];
        
    gathered.resize(shifts[fragsizes.size()-1]+fragsizes[fragsizes.size()-1]);

    MPI_Allgatherv(fragment.data(), fragment.size(), MPI_DOUBLE, gathered.data(), fragsizes.data(), shifts.data(), MPI_DOUBLE, MPI_COMM_WORLD);
}
    

void slmpi::scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment)
{
    MPI_Scatter(toscatter.data(), fragment.size(), MPI_INT, fragment.data(), fragment.size(), MPI_INT, scatterer, MPI_COMM_WORLD); 
}

void slmpi::scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment)
{
    MPI_Scatter(toscatter.data(), fragment.size(), MPI_DOUBLE, fragment.data(), fragment.size(), MPI_DOUBLE, scatterer, MPI_COMM_WORLD); 
}


void slmpi::scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment, std::vector<int>& fragsizes)
{
    std::vector<int> shifts(fragsizes.size(), 0);
    for (int i = 1; i < fragsizes.size(); i++)
        shifts[i] = shifts[i-1]+fragsizes[i-1];

    MPI_Scatterv(toscatter.data(), fragsizes.data(), shifts.data(), MPI_INT, fragment.data(), fragment.size(), MPI_INT, scatterer, MPI_COMM_WORLD); 
}

void slmpi::scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment, std::vector<int>& fragsizes)
{
    std::vector<int> shifts(fragsizes.size(), 0);
    for (int i = 1; i < fragsizes.size(); i++)
        shifts[i] = shifts[i-1]+fragsizes[i-1];

    MPI_Scatterv(toscatter.data(), fragsizes.data(), shifts.data(), MPI_DOUBLE, fragment.data(), fragment.size(), MPI_DOUBLE, scatterer, MPI_COMM_WORLD); 
}


void slmpi::exchange(std::vector<int> targetranks, std::vector<int>& sendvalues, std::vector<int>& receivevalues)
{
    int numtargets = targetranks.size();
 
    receivevalues = {};
    if (numtargets == 0 || sendvalues.size() == 0)
        return;
    
    int numvalues = sendvalues.size()/numtargets;
    receivevalues = std::vector<int>(numtargets*numvalues);

    std::vector<MPI_Request> sendrequests(numtargets);
    std::vector<MPI_Request> receiverequests(numtargets);

    for (int i = 0; i < numtargets; i++)
        MPI_Isend(&sendvalues[numvalues*i], numvalues, MPI_INT, targetranks[i], 0, MPI_COMM_WORLD, &sendrequests[i]);

    for (int i = 0; i < numtargets; i++)
        MPI_Irecv(&receivevalues[numvalues*i], numvalues, MPI_INT, targetranks[i], 0, MPI_COMM_WORLD, &receiverequests[i]);

    MPI_Waitall(numtargets, &sendrequests[0], MPI_STATUSES_IGNORE);
    MPI_Waitall(numtargets, &receiverequests[0], MPI_STATUSES_IGNORE);
}

void slmpi::exchange(std::vector<int> targetranks, std::vector<double>& sendvalues, std::vector<double>& receivevalues)
{
    int numtargets = targetranks.size();
 
    receivevalues = {};
    if (numtargets == 0 || sendvalues.size() == 0)
        return;
    
    int numvalues = sendvalues.size()/numtargets;
    receivevalues = std::vector<double>(numtargets*numvalues);

    std::vector<MPI_Request> sendrequests(numtargets);
    std::vector<MPI_Request> receiverequests(numtargets);

    for (int i = 0; i < numtargets; i++)
        MPI_Isend(&sendvalues[numvalues*i], numvalues, MPI_DOUBLE, targetranks[i], 0, MPI_COMM_WORLD, &sendrequests[i]);

    for (int i = 0; i < numtargets; i++)
        MPI_Irecv(&receivevalues[numvalues*i], numvalues, MPI_DOUBLE, targetranks[i], 0, MPI_COMM_WORLD, &receiverequests[i]);

    MPI_Waitall(numtargets, &sendrequests[0], MPI_STATUSES_IGNORE);
    MPI_Waitall(numtargets, &receiverequests[0], MPI_STATUSES_IGNORE);
}

void slmpi::exchange(std::vector<int> targetranks, std::vector<std::vector<int>>& sends, std::vector<std::vector<int>>& receives)
{
    int numtargets = targetranks.size();

    if (numtargets == 0)
        return;

    std::vector<int> sendlens(numtargets), reclens(numtargets);
    std::vector<int*> sendbuffers(numtargets), recbuffers(numtargets);
    
    for (int i = 0; i < numtargets; i++)
    {
        sendlens[i] = sends[i].size();
        reclens[i] = receives[i].size();
        
        sendbuffers[i] = sends[i].data();
        recbuffers[i] = receives[i].data();
    }

    exchange(targetranks, sendlens, sendbuffers, reclens, recbuffers);
}

void slmpi::exchange(std::vector<int> targetranks, std::vector<std::vector<double>>& sends, std::vector<std::vector<double>>& receives)
{
    int numtargets = targetranks.size();

    if (numtargets == 0)
        return;

    std::vector<int> sendlens(numtargets), reclens(numtargets);
    std::vector<double*> sendbuffers(numtargets), recbuffers(numtargets);
    
    for (int i = 0; i < numtargets; i++)
    {
        sendlens[i] = sends[i].size();
        reclens[i] = receives[i].size();
        
        sendbuffers[i] = sends[i].data();
        recbuffers[i] = receives[i].data();
    }

    exchange(targetranks, sendlens, sendbuffers, reclens, recbuffers);
}
    
void slmpi::exchange(std::vector<int> targetranks, std::vector<int> sendlens, std::vector<int*> sendbuffers, std::vector<int> receivelens, std::vector<int*> receivebuffers)
{
    int numtargets = targetranks.size();

    if (numtargets == 0)
        return;

    std::vector<MPI_Request> sendrequests(numtargets);
    std::vector<MPI_Request> receiverequests(numtargets);

    for (int i = 0; i < numtargets; i++)
        MPI_Isend(sendbuffers[i], sendlens[i], MPI_INT, targetranks[i], 0, MPI_COMM_WORLD, &sendrequests[i]);

    for (int i = 0; i < numtargets; i++)
        MPI_Irecv(receivebuffers[i], receivelens[i], MPI_INT, targetranks[i], 0, MPI_COMM_WORLD, &receiverequests[i]);

    MPI_Waitall(numtargets, &sendrequests[0], MPI_STATUSES_IGNORE);
    MPI_Waitall(numtargets, &receiverequests[0], MPI_STATUSES_IGNORE);
}

void slmpi::exchange(std::vector<int> targetranks, std::vector<int> sendlens, std::vector<double*> sendbuffers, std::vector<int> receivelens, std::vector<double*> receivebuffers)
{
    int numtargets = targetranks.size();

    if (numtargets == 0)
        return;

    std::vector<MPI_Request> sendrequests(numtargets);
    std::vector<MPI_Request> receiverequests(numtargets);

    for (int i = 0; i < numtargets; i++)
        MPI_Isend(sendbuffers[i], sendlens[i], MPI_DOUBLE, targetranks[i], 0, MPI_COMM_WORLD, &sendrequests[i]);

    for (int i = 0; i < numtargets; i++)
        MPI_Irecv(receivebuffers[i], receivelens[i], MPI_DOUBLE, targetranks[i], 0, MPI_COMM_WORLD, &receiverequests[i]);

    MPI_Waitall(numtargets, &sendrequests[0], MPI_STATUSES_IGNORE);
    MPI_Waitall(numtargets, &receiverequests[0], MPI_STATUSES_IGNORE);
}
    

std::vector<double> slmpi::ping(int messagesize, int verbosity)
{
    std::vector<double> output = {};
    std::vector<double> datavec(messagesize);
    for (int i = 0; i < messagesize; i++)
        datavec[i] = i;
    
    barrier();

    if (getrank() == 0)
    {
        output = std::vector<double>(count(), 0);
        
        if (verbosity > 0)
            std::cout << "Data round trip time (" << messagesize << " doubles):" << std::endl;
        for (int i = 1; i < count(); i++)
        {
            wallclock clk;
            send(i, 0, datavec);
            receive(i, 0, datavec);
            output[i] = clk.toc();
            if (verbosity > 0)
                clk.print("0 <> "+std::to_string(i)+":");
        }
    }
    else
    {
        receive(0, 0, datavec);
        send(0, 0, datavec);
    }
    
    // Timings should not be influenced by following computations:
    barrier();
    
    return output;
}

#endif

