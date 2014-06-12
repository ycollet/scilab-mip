function x = doe_sobol(n)
// x = doe_sobol(n)
//   quasirandom number generator
//
//   initialize with
//     sobol(-1)
//   then repeatedly use
//     x = sobol(n)
//   to generate new random vector x of n components
//   where n = 1,2,3,4,5,6
global ip mdeg ix iv in fac
maxdim=6; maxbit=30;
x = [];
if n<0  // initialization
  ip=[0 1 1 2 1 4]; mdeg=[1 2 3 3 4 4]; ix=zeros(1,6);
  iv=zeros(maxdim,maxbit);
  iv(:)=[ones(1,6),3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,zeros(1,156)];

  for k=1:maxdim
    for j=1:mdeg(k)
      iv(k,j) = iv(k,j)*2^(maxbit-j);
    end
    for j=(mdeg(k)+1):maxbit
      ipp = ip(k);
      i = iv(k,j-mdeg(k));
      i = bitxor(i,floor(i/2^mdeg(k)));
      for l=(mdeg(k)-1):-1:1
        if bitand(ipp,1) ~= 0
          i = bitxor(i,iv(k,j-l));
        end
        ipp = floor(ipp/2);
      end
      iv(k,j) = i;
    end
  end
  fac = 1/2^maxbit;
  in = 0;
else
  im = in;
  flag = 1;
  for j=1:maxbit
    if bitand(im,1) == 0
      flag = 0;
      break
    end
    im = floor(im/2);
  end
  if flag
    error('maxbit too small in sobol')
  end
  im = (j-1)*maxdim;
  for k=1:min(n,maxdim)
    ix(k) = bitxor(ix(k),iv(im+k));
    x(k) = ix(k)*fac;
  end
  in = in + 1;
end
endfunction
