function TSPVar_out = scramble_tsp(TSPVar_in, Iteration)
if (size(TSPVar_in,1)==1) then
  TSPVar_in = TSPVar_in';
end
for i=1:Iteration
  Index1 = ceil(rand(1,1)*size(TSPVar_in,1));
  Index2 = ceil(rand(1,1)*size(TSPVar_in,1));
  Aux               = TSPVar_in(Index1);
  TSPVar_in(Index1) = TSPVar_in(Index2);
  TSPVar_in(Index2) = Aux;
end
TSPVar_out = TSPVar_in;
endfunction
