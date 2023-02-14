// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <iostream>
#include <string.h>

class slexception: public std::exception
{
protected:
    std::string msg_;

public:
    explicit slexception(const std::string& message)
        : msg_(message) {}

    virtual ~slexception() noexcept {}

    virtual const char* what() const noexcept {
       return msg_.c_str();
    }
};

#endif