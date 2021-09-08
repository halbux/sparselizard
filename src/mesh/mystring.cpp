#include "mystring.h"


mystring::mystring(std::string inputstring)
{
    stringtoprocess = inputstring;
    stringindex = 0;
}

void mystring::setstringindex(int inputindex)
{
    stringindex = inputindex;
}

std::string mystring::getstring(void)
{
    return stringtoprocess;
}

int mystring::getstringindex(void)
{
    return stringindex;
}

std::string mystring::getstringtonextwhitespace(void)
{
    int indexbackup = stringindex;
    while (stringindex < stringtoprocess.length() && stringtoprocess[stringindex] != ' ')
        stringindex = stringindex + 1;
    
    stringindex = stringindex + 1;
    
    return stringtoprocess.substr(indexbackup, stringindex-indexbackup-1);
}

void mystring::jumptonextwhitespace(void)
{
    while (stringindex < stringtoprocess.length() && stringtoprocess[stringindex] != ' ')
        stringindex = stringindex + 1;
    stringindex = stringindex + 1;    
}

std::string mystring::getstringtonextcomma(void)
{
    int indexbackup = stringindex;
    while (stringindex < stringtoprocess.length() && stringtoprocess[stringindex] != ',')
        stringindex = stringindex + 1;
    
    stringindex = stringindex + 1;
    
    return stringtoprocess.substr(indexbackup, stringindex-indexbackup-1);
}

std::string mystring::getstringwhileletter(void)
{
    int indexbackup = stringindex;
    while (stringindex < stringtoprocess.length() && std::isalpha(stringtoprocess[stringindex]))
        stringindex = stringindex + 1;
    
    return stringtoprocess.substr(indexbackup, stringindex-indexbackup);
}
