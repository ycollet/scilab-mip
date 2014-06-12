#ifndef CGLAPL_H
#define CGLAPL_H

void diskhelmholtz();
void helmholtzsolve(Grid& g, Vector& sol, Vector& dif, Vector& dir, Vector& neu, 
		    Vector& rhs,Vector& rob,Vector& vis);
#endif
