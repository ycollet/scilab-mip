// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpInterfacesRegOp.hpp 1327 2008-09-18 19:01:17Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#ifndef __IPINTERFACESREGOP_HPP__
#define __IPINTERFACESREGOP_HPP__

#include "IpSmartPtr.hpp"

namespace Ipopt
{
  class RegisteredOptions;

  void RegisterOptions_Interfaces(const SmartPtr<RegisteredOptions>& roptions);

} // namespace Ipopt

#endif
