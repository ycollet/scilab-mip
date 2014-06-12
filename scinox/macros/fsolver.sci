function [x,fval,exitflag,status] = fsolver(userfun,x0,options)
    [nargout,nargin] = argn();
    //printf('Got nargout %d nargin %d\n',nargout,nargin);
    if nargin < 2
        printf('Usage: fsolver(userfun,x0)\n');
    end 

    if nargin < 3
      //Define default parameters
      options = init_param();
      options = add_param(options,'Tol',1e-6);
      options = add_param(options,'Jacobian','off');
      options = add_param(options,'DerivativeCheck','off');
    end 

    // Get the size of the system 
    n = size(x0,1); 
    
    // Define the number of function evals that occur inside this 
    // wrapper routine. These need to be added to the number of 
    // function evaluations that occur inside of the solver. 
    nfeval = 0;
    
    [jactype  , err] = get_param(options,'Jacobian','off');
    [dervcheck, err] = get_param(options,'DerivativeCheck','off');
    
    
    // Test if we need to check the derivatives 
    if dervcheck == 'on'
      // We need to compute a finite-difference 
      // approximation to the Jacobian to compare 
      // with the user-supplied analytical Jacobian
      // How we compute this approximation depends on 
      // whether the user's Jacobian is full or sparse
      if jactype == 'off'
	error('To check derivatives you must supply a function to compute the Jacobian\n' + ...
	      'and the option ''Jacobian'' must be set to ''Full'' or ' + ...
	    '''Sparse''.\n');
      elseif jactype == 'on'
	// compute the user-supplied J at x0
	[f,user_deriv] = userfun(x0);
	nfeval = nfeval + 1; 

	if issparse(user_deriv) 
	  // the user supplied a sparse Jacobian perform sparse 
	  // finite-differences 
	  
	  // compute a coloring of the column graph. 
	  // use the sparsity structure from the user-supplied
	  // Jacobian. Note this means we won't be able to 
	  // detect errors in the user-supplied sparsity pattern.
	  // That is if the user misssed out a nonzero we cannot
	  // catch it. 
	  [coloring,numcolors,mincolors] = columncolor(user_deriv);
	  
	  printf('Using %d function evaluations to estimate the ' ...
	      + 'Jacobian, at least %d required\n',numcolors,mincolors);

	  // compute the finitie-difference approximation to J at x0
	  fd_deriv = sfdjac(userfun,x0,user_deriv,coloring);
	  nfeval = nfeval + numcolors + 1;
	else 
	  // the user supplied a dense Jacobian perform dense 
	  // finite-differences
	  	  
	  // compute the finite-difference approximation to J at x0 
	  fd_deriv   = fdjac(userfun,x0);
	  nfeval = nfeval + n + 1; 

	end
	
	// compare the users Jacobian with that generated via finite-differences
	checkderivative(fd_deriv,user_deriv);
      else 
	error('Unknown option for ''Jacobian''.\n');
      end
    end
    
    // Call the solver 
    [x,fval,exitflag,status] = fsolver_simple(userfun,x0,options);

    //we need to add nfeval to status('nfev') here
    nfeval = nfeval + get_param(status,'nfev',0);
    status = set_param(status,'nfev',nfeval);

endfunction
