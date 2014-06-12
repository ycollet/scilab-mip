#ifndef MANAGE_PARAMS_HPP
#define MANAGE_PARAMS_HPP
#include <IpIpoptApplication.hpp>

int manage_bonmin_params(Ipopt::SmartPtr<Ipopt::OptionsList> options, int * param_in_addr, int Log);
int manage_ipopt_params(Ipopt::SmartPtr<Ipopt::OptionsList> options, int * param_in_addr, int Log);
#endif
