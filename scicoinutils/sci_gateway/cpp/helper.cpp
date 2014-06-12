#include <string.h>

#include <CoinMessageHandler.hpp>

extern "C" {
#include <stack-c.h>
#include <sciprint.h>
}

#include <helper.hpp>

/////////////////////
// Message handler //
/////////////////////

int DerivedHandler::print()
{
  sciprint("%s\n", messageBuffer());
  return 0;
}

