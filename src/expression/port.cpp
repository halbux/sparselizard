#include "port.h"

port::port(void)
{
    rawportptr = std::shared_ptr<rawport>(new rawport({1}, false));
}

port::port(std::vector<int> harmonicnumbers)
{
    rawportptr = std::shared_ptr<rawport>(new rawport(harmonicnumbers, true));
}

port::port(std::shared_ptr<rawport> rp)
{
    rawportptr = rp;
}

void port::setname(std::string name) { rawportptr->setname(name); }

std::string port::getname(void) { return rawportptr->getname(); }

std::vector<int> port::getharmonics(void)
{
    return rawportptr->getharmonics();
}

port port::harmonic(int harmonicnumber)
{
    return port(rawportptr->harmonic({harmonicnumber}));
}

port port::harmonic(std::vector<int> harmonicnumbers)
{
    return port(rawportptr->harmonic(harmonicnumbers));
}

port port::sin(int freqindex) { return harmonic(2*freqindex); }
port port::cos(int freqindex) { return harmonic(2*freqindex+1); }

std::shared_ptr<rawport> port::getpointer(void) { return rawportptr; }


expression port::operator+(void) { return (expression)*this; }
expression port::operator-(void) { return -(expression)*this; }

expression port::operator+(port inputport) { return (expression)*this + inputport; }
expression port::operator-(port inputport) { return (expression)*this - inputport; }
expression port::operator*(port inputport) { return (expression)*this * inputport; }
expression port::operator/(port inputport) { return (expression)*this / inputport; }

expression port::operator+(double val) { return (expression)*this + val; }
expression port::operator-(double val) { return (expression)*this - val; }
expression port::operator*(double val) { return (expression)*this * val; }
expression port::operator/(double val) { return (expression)*this / val; } 


expression operator+(double val, port inputport) { return (expression)val + inputport; }
expression operator-(double val, port inputport) { return (expression)val - inputport; }
expression operator*(double val, port inputport) { return (expression)val * inputport; }
expression operator/(double val, port inputport) { return (expression)val / inputport; }

