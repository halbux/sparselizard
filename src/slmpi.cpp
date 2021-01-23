#include "slmpi.h"
#include "wallclock.h"


void slmpi::errornompi(void)
{
    std::cout << "Error in 'slmpi' namespace: MPI is not available" << std::endl;
    abort();
}


#ifndef HAVE_MPI
bool slmpi::isavailable(void) { return false; }
void slmpi::initialize(void) {}
void slmpi::finalize(void) {}
int slmpi::getrank(void) { return 0; }
int slmpi::count(void) { return 1; }
void slmpi::barrier(void) {}
void slmpi::send(int destination, int tag, int len, int* data) { errornompi(); }
void slmpi::send(int destination, int tag, int len, double* data) { errornompi(); }
void slmpi::send(int destination, int tag, std::vector<int>& data) { errornompi(); }
void slmpi::send(int destination, int tag, std::vector<double>& data) { errornompi(); }
void slmpi::receive(int source, int tag, int len, int* data) { errornompi(); }
void slmpi::receive(int source, int tag, int len, double* data) { errornompi(); }
void slmpi::receive(int source, int tag, std::vector<int>& data) { errornompi(); }
void slmpi::receive(int source, int tag, std::vector<double>& data) { errornompi(); }
void slmpi::exchange(std::vector<int> targetranks, std::vector<int> sendlens, std::vector<int*> sendbuffers, std::vector<int> receivelens, std::vector<int*> receivebuffers) { errornompi(); }
void slmpi::exchange(std::vector<int> targetranks, std::vector<int> sendlens, std::vector<double*> sendbuffers, std::vector<int> receivelens, std::vector<double*> receivebuffers) { errornompi(); }
void slmpi::sum(int len, int* data) {}
void slmpi::sum(int len, double* data) {}
void slmpi::sum(std::vector<int>& data) {}
void slmpi::sum(std::vector<double>& data) {}
void slmpi::broadcast(int broadcaster, std::vector<int>& data) { errornompi(); }
void slmpi::broadcast(int broadcaster, std::vector<double>& data) { errornompi(); }
void slmpi::gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered) { errornompi(); }
void slmpi::gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered) { errornompi(); }
void slmpi::gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered, std::vector<int>& fragsizes) { errornompi(); }
void slmpi::gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered, std::vector<int>& fragsizes) { errornompi(); }
void slmpi::scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment) { errornompi(); }
void slmpi::scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment) { errornompi(); }
void slmpi::scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment, std::vector<int>& fragsizes) { errornompi(); }
void slmpi::scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment, std::vector<int>& fragsizes) { errornompi(); }
std::vector<double> slmpi::ping(int messagesize, int verbosity) { errornompi(); abort(); }
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


void slmpi::send(int destination, int tag, int len, int* data)
{
    MPI_Send(data, len, MPI_INT, destination, tag, MPI_COMM_WORLD);
}

void slmpi::send(int destination, int tag, int len, double* data)
{
    MPI_Send(data, len, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
}

void slmpi::send(int destination, int tag, std::vector<int>& data)
{
    MPI_Send(&data[0], data.size(), MPI_INT, destination, tag, MPI_COMM_WORLD);
}

void slmpi::send(int destination, int tag, std::vector<double>& data)
{
    MPI_Send(&data[0], data.size(), MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
}


void slmpi::receive(int source, int tag, int len, int* data)
{
    MPI_Recv(data, len, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void slmpi::receive(int source, int tag, int len, double* data)
{
    MPI_Recv(data, len, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void slmpi::receive(int source, int tag, std::vector<int>& data)
{
    MPI_Recv(&data[0], data.size(), MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void slmpi::receive(int source, int tag, std::vector<double>& data)
{
    MPI_Recv(&data[0], data.size(), MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


void slmpi::exchange(std::vector<int> targetranks, std::vector<int> sendlens, std::vector<int*> sendbuffers, std::vector<int> receivelens, std::vector<int*> receivebuffers)
{
    int numtargets = targetranks.size();
    
    if (numtargets == 0)
        return;

    int totbytelen = 0;
    for (int i = 0; i < numtargets; i++)
        totbytelen += sendlens[i]*sizeof(int) + MPI_BSEND_OVERHEAD;

    std::vector<char> sendbuffer(totbytelen); // a char is one byte long
    MPI_Buffer_attach(&sendbuffer[0], totbytelen);
     
    for (int i = 0; i < numtargets; i++)
        MPI_Bsend(sendbuffers[i], sendlens[i], MPI_INT, targetranks[i], 0, MPI_COMM_WORLD);
     
    for (int i = 0; i < numtargets; i++)
        MPI_Recv(receivebuffers[i], receivelens[i], MPI_INT, targetranks[i], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     
    MPI_Buffer_detach(&sendbuffer[0], &totbytelen);
}

void slmpi::exchange(std::vector<int> targetranks, std::vector<int> sendlens, std::vector<double*> sendbuffers, std::vector<int> receivelens, std::vector<double*> receivebuffers)
{
    int numtargets = targetranks.size();

    if (numtargets == 0)
        return;
        
    int totbytelen = 0;
    for (int i = 0; i < numtargets; i++)
        totbytelen += sendlens[i]*sizeof(double) + MPI_BSEND_OVERHEAD;

    std::vector<char> sendbuffer(totbytelen); // a char is one byte long
    MPI_Buffer_attach(&sendbuffer[0], totbytelen);
     
    for (int i = 0; i < numtargets; i++)
        MPI_Bsend(sendbuffers[i], sendlens[i], MPI_DOUBLE, targetranks[i], 0, MPI_COMM_WORLD);
     
    for (int i = 0; i < numtargets; i++)
        MPI_Recv(receivebuffers[i], receivelens[i], MPI_DOUBLE, targetranks[i], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     
    MPI_Buffer_detach(&sendbuffer[0], &totbytelen);
}


void slmpi::sum(int len, int* data)
{
    MPI_Allreduce(MPI_IN_PLACE, data, len, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void slmpi::sum(int len, double* data)
{
    MPI_Allreduce(MPI_IN_PLACE, data, len, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void slmpi::sum(std::vector<int>& data)
{
    MPI_Allreduce(MPI_IN_PLACE, &data[0], data.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void slmpi::sum(std::vector<double>& data)
{
    MPI_Allreduce(MPI_IN_PLACE, &data[0], data.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


void slmpi::broadcast(int broadcaster, std::vector<int>& data)
{
    MPI_Bcast(&data[0], data.size(), MPI_INT, broadcaster, MPI_COMM_WORLD);
}

void slmpi::broadcast(int broadcaster, std::vector<double>& data)
{
    MPI_Bcast(&data[0], data.size(), MPI_DOUBLE, broadcaster, MPI_COMM_WORLD);
}


void slmpi::gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered)
{
    gathered = {};
    
    if (getrank() == gatherer)
        gathered.resize(count()*fragment.size());

    MPI_Gather(&fragment[0], fragment.size(), MPI_INT, &gathered[0], fragment.size(), MPI_INT, gatherer, MPI_COMM_WORLD); 
}

void slmpi::gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered)
{
    gathered = {};
    
    if (getrank() == gatherer)
        gathered.resize(count()*fragment.size());

    MPI_Gather(&fragment[0], fragment.size(), MPI_DOUBLE, &gathered[0], fragment.size(), MPI_DOUBLE, gatherer, MPI_COMM_WORLD); 
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

    MPI_Gatherv(&fragment[0], fragment.size(), MPI_INT, &gathered[0], &fragsizes[0], &shifts[0], MPI_INT, gatherer, MPI_COMM_WORLD); 
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

    MPI_Gatherv(&fragment[0], fragment.size(), MPI_DOUBLE, &gathered[0], &fragsizes[0], &shifts[0], MPI_DOUBLE, gatherer, MPI_COMM_WORLD); 
}


void slmpi::scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment)
{
    MPI_Scatter(&toscatter[0], fragment.size(), MPI_INT, &fragment[0], fragment.size(), MPI_INT, scatterer, MPI_COMM_WORLD); 
}

void slmpi::scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment)
{
    MPI_Scatter(&toscatter[0], fragment.size(), MPI_DOUBLE, &fragment[0], fragment.size(), MPI_DOUBLE, scatterer, MPI_COMM_WORLD); 
}


void slmpi::scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment, std::vector<int>& fragsizes)
{
    std::vector<int> shifts(fragsizes.size(), 0);
    for (int i = 1; i < fragsizes.size(); i++)
        shifts[i] = shifts[i-1]+fragsizes[i-1];

    MPI_Scatterv(&toscatter[0], &fragsizes[0], &shifts[0], MPI_INT, &fragment[0], fragment.size(), MPI_INT, scatterer, MPI_COMM_WORLD); 
}

void slmpi::scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment, std::vector<int>& fragsizes)
{
    std::vector<int> shifts(fragsizes.size(), 0);
    for (int i = 1; i < fragsizes.size(); i++)
        shifts[i] = shifts[i-1]+fragsizes[i-1];

    MPI_Scatterv(&toscatter[0], &fragsizes[0], &shifts[0], MPI_DOUBLE, &fragment[0], fragment.size(), MPI_DOUBLE, scatterer, MPI_COMM_WORLD); 
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



// The functions below do not call MPI themselves and therefore do not need to be in the mpi available check:

std::vector<int> slmpi::broadcastgathered(std::vector<int> valuestoinclude, std::vector<int>& fragment, std::vector<int>& data)
{
    int numranks = count();

    int il = valuestoinclude.size();
    int fl = fragment.size();
    int tl = il+fl;

    std::vector<int> tosend(tl);
    for (int i = 0; i < il; i++)
        tosend[i] = valuestoinclude[i];
    for (int i = 0; i < fl; i++)
        tosend[il+i] = fragment[i];

    gather(0, tosend, data);
    data.resize(numranks*tl); // only rank 0 has correct data size
    broadcast(0, data);

    std::vector<int> output(numranks*il);
    
    for (int r = 0; r < numranks; r++)
    {
        for (int i = 0; i < il; i++)
            output[r*il+i] = data[r*tl + i];
        for (int i = 0; i < fl; i++)
            data[r*fl+i] = data[r*tl + il+i];
    }
    data.resize(numranks*fl);

    return output;
}

std::vector<double> slmpi::broadcastgathered(std::vector<double> valuestoinclude, std::vector<double>& fragment, std::vector<double>& data)
{
    int numranks = count();

    int il = valuestoinclude.size();
    int fl = fragment.size();
    int tl = il+fl;

    std::vector<double> tosend(tl);
    for (int i = 0; i < il; i++)
        tosend[i] = valuestoinclude[i];
    for (int i = 0; i < fl; i++)
        tosend[il+i] = fragment[i];

    gather(0, tosend, data);
    data.resize(numranks*tl); // only rank 0 has correct data size
    broadcast(0, data);

    std::vector<double> output(numranks*il);
    
    for (int r = 0; r < numranks; r++)
    {
        for (int i = 0; i < il; i++)
            output[r*il+i] = data[r*tl + i];
        for (int i = 0; i < fl; i++)
            data[r*fl+i] = data[r*tl + il+i];
    }
    data.resize(numranks*fl);

    return output;
}

