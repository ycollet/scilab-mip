function L = localise2d(t)
//
// L = Localise(T)
// T : table de connectivites de l'element  
// L : table de localisation 
//
// A. Seghir, 08/08/04 

//Preferable d'utiliser
//L = kron(2*t-1,[1,0])+kron(2*t,[0,1]);
L = [];
nne = length(t);
for i= 1:nne
  L([2*i-1 2*i]) = [2*t(i)-1 2*t(i)];
end
endfunction


 
