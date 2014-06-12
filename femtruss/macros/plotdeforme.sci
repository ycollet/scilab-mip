function plotdeforme(u,p,t,s)
//
// plotdef(U,p,t,s) : plot deformed shap  
// U : FE solution k . u = f    ( u = k \ f )
// p : table des coordonees nodales
// t : table des connectivites des elements
// s : scal factor

drawlater;
plotmesh(t,p,0,0,'green');

_3D_problem = (size(p,2)==3);

if _3D_problem then
  p = p + s * matrix([u(1:3:$),u(2:3:$),u(3:3:$)],size(p,1),size(p,2));
else
  p = p + s * matrix([u(1:2:$),u(2:2:$)],size(p,1),size(p,2));
end

plotmesh(t,p,0,0,'red'); 
legends(['before loading','after loading'], [color('green') color('red')], 'ur');
drawnow;
// axis equal;
endfunction
