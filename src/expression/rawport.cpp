#include "rawport.h"

rawport::rawport(std::vector<int> harmonicnumbers, bool ismultiharm)
{
    amimultiharmonic = ismultiharm;

    if (harmonicnumbers.size() != 0)
    {
        myharmonics.resize(*max_element(harmonicnumbers.begin(), harmonicnumbers.end())+1);
        for (int h = 0; h < harmonicnumbers.size(); h++)
            myharmonics[harmonicnumbers[h]] = { std::shared_ptr<rawport>(new rawport({}, ismultiharm)) };
    }    
    // Nothing to preallocate otherwise
}

void rawport::setvalue(double portval)
{
    if (myharmonics.size() == 0)
        myvalue = portval;
    else
    {
        if (myharmonics.size() == 2 && myharmonics[1].size() > 0)
            myharmonics[1][0]->setvalue(portval);
        else
        {
            logs log;
            log.msg() << "Error in 'rawport' object: cannot set the value of a multiharmonic port (only constant harmonic 1)" << std::endl;
            log.error();
        }
    }
}

double rawport::getvalue(void)
{
    if (myharmonics.size() == 0)
        return myvalue;
    else
    {
        if (myharmonics.size() == 2 && myharmonics[1].size() > 0)
            return myharmonics[1][0]->getvalue();
        else
        {
            logs log;
            log.msg() << "Error in 'rawport' object: cannot get the value of a multiharmonic port (only constant harmonic 1)" << std::endl;
            log.error();
        }
    }
    
    throw std::runtime_error(""); // fix return warning
}

void rawport::setname(std::string name)
{
    myname = name;
    
    for (int i = 0; i < myharmonics.size(); i++)
    {
        if (myharmonics[i].size() > 0)
            myharmonics[i][0]->setname(name);
    }
}

std::string rawport::getname(void) { return myname; }

bool rawport::isharmonicone(void)
{
    if (myharmonics.size() == 0 || myharmonics.size() == 2 && myharmonics[1].size() == 1)
        return true;
    else
        return false;
}

std::vector<int> rawport::getharmonics(void)
{
    std::vector<int> harms = {};

    for (int h = 0; h < myharmonics.size(); h++)
    {
        if (myharmonics[h].size() > 0)
            harms.push_back(h);
    }
    if (myharmonics.size() == 0)
        harms = {1};

    return harms;
}

std::shared_ptr<rawport> rawport::harmonic(int harmonicnumber)
{
    return harmonic(std::vector<int>{harmonicnumber});
}

std::shared_ptr<rawport> rawport::harmonic(std::vector<int> harmonicnumbers)
{
    if (myharmonics.size() != 0)
    {
        if (harmonicnumbers.size() == 1)
        {
            // If there is a single harmonic we can not put it into a new rawport!
            if (harmonicnumbers[0] < myharmonics.size() && myharmonics[harmonicnumbers[0]].size() > 0)
                return myharmonics[harmonicnumbers[0]][0];
            else
            {
                logs log;
                log.msg() << "Error in 'rawport' object: in .harmonic cannot get harmonic " << harmonicnumbers[0] << " (does not exist)" << std::endl; 
                log.error();
            }
        }
        
        // In case we want several harmonics a new raw port container is created.
        int maxharmnum = *max_element(harmonicnumbers.begin(), harmonicnumbers.end());

        std::shared_ptr<rawport> harmsrawport(new rawport());
        *harmsrawport = *this;

        // Set a brand new harmonic vector:
        harmsrawport->myharmonics = std::vector<std::vector<std::shared_ptr<rawport>>>(maxharmnum+1, std::vector<std::shared_ptr<rawport>>(0));
        
        // Add only the requested harmonics:
        for (int i = 0; i < harmonicnumbers.size(); i++)
        {
            if (harmonicnumbers[i] < myharmonics.size() && myharmonics[harmonicnumbers[i]].size() > 0)
                harmsrawport->myharmonics[harmonicnumbers[i]] = {myharmonics[harmonicnumbers[i]][0]};
            else
            {
                logs log;
                log.msg() << "Error in 'rawport' object: in .harmonic cannot get harmonic " << harmonicnumbers[i] << " (does not exist)" << std::endl; 
                log.error();
            }
        }
        return harmsrawport;
    }
        
    // In case there is no harmonic it is a constant (harmonic 1).
    if (harmonicnumbers.size() == 1 && harmonicnumbers[0] == 1)
        return shared_from_this();
    else
    {
        logs log;
        log.msg() << "Error in 'rawport' object: in .harmonic cannot get harmonic in constant port (does not exist)" << std::endl; 
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}
        
bool rawport::isprimal(void)
{
    return myisprimal;
}

int rawport::getphysicalregion(void)
{
    return myphysreg;
}

std::shared_ptr<rawfield> rawport::getrawfield(void)
{
    if (myrawfield.expired())
    {
        logs log;
        log.msg() << "Error in 'rawport' object: the associated rawfield is needed but it was destroyed" << std::endl;
        log.error();
    }
    
    return myrawfield.lock();
}

std::shared_ptr<rawport> rawport::getprimal(void)
{
    if (mybrother.expired())
    {
        logs log;
        log.msg() << "Error in 'rawport' object: the associated rawport is needed but it was destroyed" << std::endl;
        log.error();
    }
    
    if (myisprimal)
        return shared_from_this();
    else
        return mybrother.lock();
}

std::shared_ptr<rawport> rawport::getdual(void)
{
    if (mybrother.expired())
    {
        logs log;
        log.msg() << "Error in 'rawport' object: the associated rawport is needed but it was destroyed" << std::endl;
        log.error();
    }
    
    if (myisprimal)
        return mybrother.lock();
    else
        return shared_from_this();
}

bool rawport::isassociated(void)
{
    return (myphysreg != -1);
}

void rawport::associate(bool isprim, std::shared_ptr<rawport> bro, int physreg, std::shared_ptr<rawfield> rf)
{
    myisprimal = isprim;
    mybrother = bro;
    myphysreg = physreg;
    myrawfield = rf;
}

