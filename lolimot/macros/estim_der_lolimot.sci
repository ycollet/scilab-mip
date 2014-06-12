//////////////////////////////////////////////////////////////////////////
// estim_der_lolimot                                                    //
// For a learnt Lolimot model, computes the estimation of a given point //
//////////////////////////////////////////////////////////////////////////

function f = estim_der_lolimot(x,lolModel)
if lolModel('type')=='exp' then
  //////////////////////////////////////////////////////////////////////////////////
  // Expression of the partial derivative for the exponential membership function //
  //////////////////////////////////////////////////////////////////////////////////
  
  u  = 0;
  v  = 0;
  for i=1:size(lolModel('listofcutplus'),1)
    // u(x) = sum(i=1,Np) (Li(x) . Ei(x))
    xc = (lolModel('listofcutplus')(i,:) + lolModel('listofcutinf')(i,:))/2;
    dx = (lolModel('listofcutplus')(i,:) - lolModel('listofcutinf')(i,:))*lolModel('sigma');
    Li = lolModel('listofmod')(i,1) + sum(x(:)' .* lolModel('listofmod')(i,2:$));
    
    u = u + Li * exp(- sum((x(:)' - xc(:)') ./ dx(:)').^2);
    // v(x) = sum(i=1,Np) (Ei(x))
    v = v + exp(- sum((x(:)' - xc(:)') ./ dx(:)').^2);
  end

  for i=1:length(x)
    dv = 0;
    du = 0;
    for j=1:size(lolModel('listofcutplus'),1)
      xc = (lolModel('listofcutplus')(j,:) + lolModel('listofcutinf')(j,:))/2;
      dx = ((lolModel('listofcutplus')(j,:) - lolModel('listofcutinf')(j,:))*lolModel('sigma'));
      Lj = (lolModel('listofmod')(j,1) + sum(x(:)' .* lolModel('listofmod')(j,2:$)));
      
      // dvk(x) = -2 * sum(i=1,Np) (((xk-xk_centre)/sigmak^2) . Ei(x))
      dv = dv - 2 * ((x(i) - xc(i)) / dx(i)^2) * exp(- sum((x(:)' - xc(:)') ./ dx(:)').^2);

      // duk(x) = sum(i=1,Np)((ak - 2 * ((xk-xk_centre)/sigmak^2).Li(x)).Ei(x))
      du = du + (lolModel('listofmod')(j,i+1) - 2*((x(i) - xc(i)) / dx(i)^2) * Lj) * exp(- sum((x(:)' - xc(:)') ./ dx(:)').^2);
    end

    f(i,1) = (du * v - dv * u) / max(v^2,%eps);
  end
else
  ///////////////////////////////////////////////////////////////////////////////////////
  // Expression of the partial derivative for the piecewise linear membership function //
  ///////////////////////////////////////////////////////////////////////////////////////
  
  k_overlap = 0.5/lolModel('sigma');
  L   = [];
  Mf  = [];
  dMf = [];
  U = 0; dU = [];
  V = 0; dV = [];

  for i=1:size(lolModel('listofcutinf'),1)
    xc = (lolModel('listofcutinf')(i,:) + lolModel('listofcutplus')(i,:))/2;
    dx = (lolModel('listofcutplus')(i,:) - lolModel('listofcutinf')(i,:))*lolModel('sigma');
    
    // Compute the value of membership functions
    // Mf = [nb_part x nb_dim]  Mf(i,j) = 1;
    
    for j=1:length(xc)
      if ((xc(j) - dx(j)/2 < x(j))&(x(j) < xc(j) + dx(j)/2)) then
        phi_aux = 1.0;
      elseif ((lolModel('listofcutinf')(i,j) - k_overlap * dx(j)/2 <= x(j))&(x(j) <= xc(j) - dx(j)/2)) then
        phi_aux = (x(j) - (lolModel('listofcutinf')(i,j) - k_overlap * dx(j)/2)) / ...
                  ((xc(j) - dx(j)/2) - (lolModel('listofcutinf')(i,j) - k_overlap * dx(j)/2));
      elseif ((xc(j) + dx(j)/2 <= x(j))&(x(j) <= lolModel('listofcutplus')(i,j) + k_overlap * dx(j)/2)) then
        phi_aux = (x(j) - (lolModel('listofcutplus')(i,j) + k_overlap * dx(j)/2)) / ...
                  ((xc(j) + dx(j)/2) - (lolModel('listofcutplus')(i,j) + k_overlap * dx(j)/2));
      else
        phi_aux = 0;
      end
      Mf(i,j) = phi_aux;
    end

    // Compute the value of the linear models
    // L = [nb_part x 1]
    
    L(i,1) = lolModel('listofmod')(i,1) + sum(x(:)' .* lolModel('listofmod')(i,2:$));
    
    // Compute the derivative of the membership function
    // dMf = [nb_part x nb_dim]

    // Which dimension of the current partition has the smallest value of membership function ?
    [_tmp,Index] = min(Mf(i,:)); Index = Index(1);

    dMf_tmp = 0;
    if ((lolModel('listofcutinf')(i,Index) - k_overlap * dx(Index)/2 <= x(Index))&(x(Index) <= xc(Index) - dx(Index)/2)) then
      dMf_tmp = 1 / ((xc(Index) - dx(Index)/2) - (lolModel('listofcutinf')(i,Index) - k_overlap * dx(Index)/2));
    end
    if ((xc(Index) + dx(Index)/2 <= x(Index))&(x(Index) <= lolModel('listofcutplus')(i,Index) + k_overlap * dx(Index)/2)) then
      dMf_tmp = 1 / ((xc(Index) + dx(Index)/2) - (lolModel('listofcutplus')(i,Index) + k_overlap * dx(Index)/2));
    end
    dMf(i,:)     = zeros(1,length(x));
    dMf(i,Index) = dMf_tmp;

    // Compute U
    U = U + Mf(i,Index)*L(i);
    // Compute V
    V = V + Mf(i,Index);
    // Compute dU/dxk
    dU = dU + lolModel('listofmod')(i,2:$) * Mf(i,Index) + dMf(i,:) * L(i);
    // Compute dV/dxk
    dV = dV + dMf(i,:);
  end
  // Compute df/dxk
  f = (dU*V-U*dV)/max(V^2,%eps);
end
endfunction

