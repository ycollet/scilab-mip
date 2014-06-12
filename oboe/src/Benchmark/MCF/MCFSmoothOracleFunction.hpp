// Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
//
// OBOE is published under the Common Public License.
//
// Authors :
// Nidhi Sawhney <nsawhney@yahoo.com>
// The OBOE team
//

#ifndef MCF_SMOOTH_ORACLE_FUNCTION_HPP
#define MCF_SMOOTH_ORACLE_FUNCTION_HPP

#include "MCFOracleFunction.hpp"

class MCFSmoothOracleFunction : public OracleFunction {
 
 private:
  const MCFData *_data;
 
 public:
  MCFSmoothOracleFunction(const MCFData *data);
  virtual ~MCFSmoothOracleFunction();
  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, 
		   AccpmGenMatrix &subGradients, AccpmGenMatrix *info);
};

#endif
