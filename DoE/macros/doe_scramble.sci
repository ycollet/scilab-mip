function H = doe_scramble(H1, N)
if (~isdef('N','local')) then N = size(H1,1); end

for i=1:N
  Index1 = ceil(rand(1,1)*size(H1,1));
  Index2 = ceil(rand(1,1)*size(H1,1));
  Aux          = H1(Index1,:);
  H1(Index1,:) = H1(Index2,:);
  H1(Index2,:) = Aux;
end

H = H1;
endfunction
