// Function which computes the inverse of a matrix 
// Given the inverse of a matrix and the differences between the desired matrix and the given matrix
// Work done by Ms. Ibrahimcha, Nkamani, Hammad.

function [Kinv_new]= sherwoodmor(K_inv,Kdiff)
NonNullLines = [];
for i=1:size(Kdiff,1)
  Index = find(Kdiff(i,:)~=0);
  if ~isempty(Index) then
    NonNullLines = [NonNullLines, i];
  end
end
Vs = Kdiff(NonNullLines,:);
Us = zeros(size(Kdiff,1),length(NonNullLines));
Us(NonNullLines,:) = eye(length(NonNullLines),length(NonNullLines));

tmp = K_inv*Us;
Kinv_new = K_inv - tmp*inv(eye(length(NonNullLines),length(NonNullLines)) + Vs*tmp)*Vs*K_inv;
endfunction 

