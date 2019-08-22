#include "legendre.h"

std::vector<polynomial> legendre::l(int maxn, polynomial x)
{
    std::vector<polynomial> l(maxn+1);
    // l[0] is 1, l[1] is x:
    l[0].set({{{1.0}}});
    if (maxn > 0)
        l[1] = x;
    // Define the remaining orders recursively:
    for (int n = 1; n < maxn; n++)
        l[n+1] = 1.0/(n+1)*( (2*n+1)*l[n]*x-n*l[n-1] );
        
    return l;
}

std::vector<polynomial> legendre::L(int maxn, polynomial x)
{
    if (maxn == 0)
        return std::vector<polynomial>(0);

    std::vector<polynomial> L(maxn+1);
    // L[1] is x, L[2] is 1/2(x^2-1):
    L[1] = x;
    if (maxn > 1)
        L[2] = 0.5*(x*x-1);
    // Define the remaining orders recursively:
    for (int n = 2; n < maxn; n++)
        L[n+1] = 1.0/(n+1)*( (2*n-1)*x*L[n]-(n-2)*L[n-1] );    
        
    return L;
}

std::vector<polynomial> legendre::ls(int maxn, polynomial x, polynomial t)
{
    std::vector<polynomial> ls(maxn+1);
    // ls[0] is 1, ls[1] is x:
    ls[0].set({{{1.0}}});
    if (maxn > 0)
        ls[1] = x;
    // Define the remaining orders recursively:
    for (int n = 1; n < maxn; n++)
        ls[n+1] = 1.0/(n+1) * ( (2*n+1)*x*ls[n]-n*t*t*ls[n-1] );
        
    return ls;
}

std::vector<polynomial> legendre::Ls(int maxn, polynomial x, polynomial t)
{
    if (maxn == 0)
        return std::vector<polynomial>(0);

    std::vector<polynomial> Ls(maxn+1);
    // Ls[1] is x, Ls[2] is 1/2(x^2-t^2):
    Ls[1] = x;
    if (maxn > 1)
        Ls[2] = 0.5*(x*x-t*t);
    // Define the remaining orders recursively:
    for (int n = 2; n < maxn; n++)
        Ls[n+1] = 1.0/(n+1) * ( (2*n-1)*x*Ls[n]-(n-2)*t*t*Ls[n-1] );    
        
    return Ls;
}
