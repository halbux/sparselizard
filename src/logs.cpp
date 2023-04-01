#include "logs.h"


std::ostream& logs::msg(void)
{
    return message;
}

void logs::error(void)
{
    throw std::runtime_error(message.str());
}
