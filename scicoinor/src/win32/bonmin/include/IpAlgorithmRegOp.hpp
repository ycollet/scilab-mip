// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpAlgorithmRegOp.hpp 1327 2008-09-18 19:01:17Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#ifndef __IPALGORITHMREGOP_HPP__
#define __IPALGORITHMREGOP_HPP__

#include "IpSmartPtr.hpp"

namespace Ipopt
{
  class RegisteredOptions;

  void RegisterOptions_Algorithm(const SmartPtr<RegisteredOptions>& roptions);

} // namespace Ipopt

#endif
