function plot_truss(p,t)
//
// plot_truss(U,p,t,s) : plot the shape of the truss structure
// p : table des coordonees nodales
// t : table des connectivites des elements

drawlater;
plotmesh(t,p,0,0,'green');
drawnow;
endfunction
