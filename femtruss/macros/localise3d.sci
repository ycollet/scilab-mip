function L = localise3d(t)
//
// L = Localise(T)
// T : table de connectivites de l'element  
// L : table de localisation 
//
// A. Seghir, 08/08/04 
L = [];
nne = length(t);
for i= 1:nne
  L([3*i-2 3*i-1 3*i]) = [3*t(i)-2 3*t(i)-1 3*t(i)];
end
endfunction


 
