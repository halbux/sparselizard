#include "slmpi.h"
#include "myalgorithm.h"
#include "wallclock.h"


void slmpi::errornompi(void)
{
    std::cout << "Error in 'slmpi' namespace: MPI is not available" << std::endl;
    abort();
}


#ifndef HAVE_MPI
bool slmpi::isavailable(void) { return false; }
void slmpi::initialize(void) { errornompi(); }
void slmpi::finalize(void) { errornompi(); }
int slmpi::getrank(void) { errornompi(); abort(); }
int slmpi::count(void) { errornompi(); abort(); }
void slmpi::barrier(void) { errornompi(); }
void slmpi::send(int destination, int tag, std::vector<int>& data) { errornompi(); }
void slmpi::send(int destination, int tag, std::vector<double>& data) { errornompi(); }
void slmpi::receive(int source, int tag, std::vector<int>& data) { errornompi(); }
void slmpi::receive(int source, int tag, std::vector<double>& data) { errornompi(); }
void slmpi::broadcast(int broadcaster, std::vector<int>& data) { errornompi(); }
void slmpi::broadcast(int broadcaster, std::vector<double>& data) { errornompi(); }
std::vector<int> slmpi::gather(int gatherer, int value) { errornompi(); abort(); }
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


void slmpi::send(int destination, int tag, std::vector<int>& data)
{
    MPI_Send(&data[0], data.size(), MPI_INT, destination, tag, MPI_COMM_WORLD);
}

void slmpi::send(int destination, int tag, std::vector<double>& data)
{
    MPI_Send(&data[0], data.size(), MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
}


void slmpi::receive(int source, int tag, std::vector<int>& data)
{
    MPI_Recv(&data[0], data.size(), MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void slmpi::receive(int source, int tag, std::vector<double>& data)
{
    MPI_Recv(&data[0], data.size(), MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


void slmpi::broadcast(int broadcaster, std::vector<int>& data)
{
    MPI_Bcast(&data[0], data.size(), MPI_INT, broadcaster, MPI_COMM_WORLD);
}

void slmpi::broadcast(int broadcaster, std::vector<double>& data)
{
    MPI_Bcast(&data[0], data.size(), MPI_DOUBLE, broadcaster, MPI_COMM_WORLD);
}

std::vector<int> slmpi::gather(int gatherer, int value)
{
    std::vector<int> gathered = {};
    
    if (getrank() == gatherer)
        gathered.resize(count());

    MPI_Gather(&value, 1, MPI_INT, &gathered[0], 1, MPI_INT, gatherer, MPI_COMM_WORLD); 
    
    return gathered;
}


void slmpi::gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered)
{
    if (getrank() == gatherer)
        gathered.resize(count()*fragment.size());

    MPI_Gather(&fragment[0], fragment.size(), MPI_INT, &gathered[0], fragment.size(), MPI_INT, gatherer, MPI_COMM_WORLD); 
}

void slmpi::gather(int gatherer, std::vector<double>& fragment, std::vector<double>& gathered)
{
    if (getrank() == gatherer)
        gathered.resize(count()*fragment.size());

    MPI_Gather(&fragment[0], fragment.size(), MPI_DOUBLE, &gathered[0], fragment.size(), MPI_INT, gatherer, MPI_COMM_WORLD); 
}


void slmpi::gather(int gatherer, std::vector<int>& fragment, std::vector<int>& gathered, std::vector<int>& fragsizes)
{
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
    fragment.resize(toscatter.size()/count());
        
    MPI_Scatter(&toscatter[0], fragment.size(), MPI_INT, &fragment[0], fragment.size(), MPI_INT, scatterer, MPI_COMM_WORLD); 
}

void slmpi::scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment)
{
    fragment.resize(toscatter.size()/count());

    MPI_Scatter(&toscatter[0], fragment.size(), MPI_DOUBLE, &fragment[0], fragment.size(), MPI_DOUBLE, scatterer, MPI_COMM_WORLD); 
}


void slmpi::scatter(int scatterer, std::vector<int>& toscatter, std::vector<int>& fragment, std::vector<int>& fragsizes)
{
    std::vector<int> shifts(fragsizes.size(), 0);
    for (int i = 1; i < fragsizes.size(); i++)
        shifts[i] = shifts[i-1]+fragsizes[i-1];

    MPI_Scatterv(&toscatter[0], &fragsizes[0], &shifts[0], MPI_INT, &fragment[0], fragsizes.size(), MPI_INT, scatterer, MPI_COMM_WORLD); 
}

void slmpi::scatter(int scatterer, std::vector<double>& toscatter, std::vector<double>& fragment, std::vector<int>& fragsizes)
{
    std::vector<int> shifts(fragsizes.size(), 0);
    for (int i = 1; i < fragsizes.size(); i++)
        shifts[i] = shifts[i-1]+fragsizes[i-1];

    MPI_Scatterv(&toscatter[0], &fragsizes[0], &shifts[0], MPI_DOUBLE, &fragment[0], fragsizes.size(), MPI_DOUBLE, scatterer, MPI_COMM_WORLD); 
}
    

std::vector<double> slmpi::ping(int messagesize, int verbosity)
{
    std::vector<double> output(count(), 0);
    std::vector<double> datavec(messagesize);
    for (int i = 0; i < messagesize; i++)
        datavec[i] = i;
    
    barrier();

    if (getrank() == 0)
    {
        std::cout << "Send + receive time for " << messagesize << " doubles:" << std::endl;
        for (int i = 1; i < count(); i++)
        {
            wallclock clk;
            send(i, 0, datavec);
            receive(i, 0, datavec);
            output[i] = clk.toc();
            clk.print("0->"+std::to_string(i)+"->0:");
        }
        
        return output;
    }
    else
    {
        receive(0, 0, datavec);
        send(0, 0, datavec);
    }
    
    return {};
}

#endif

