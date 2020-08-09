function [] = ff_adaptmesh(str,hmin,hmax,err,errg,nbvx,nbsmooth,nbjacoby,ratio,omega,iso,abserror,cutoff,verbosity,inquire,splitpbedge,maxsubdiv,anisomax,rescaling)
  if exists('hmin') then
    str = str + ' hmin=' + string(hmin) + ',';
  end
  if exists('hmax') then
    str = str + ' hmax=' + string(hmax) + ',';
  end
  if exists('err') then
    str = str + ' err=' + string(err) + ',';
  end
  if exists('errg') then
    str = str + ' errg=' + string(errg) + ',';
  end
  if exists('nbvx') then
    str = str + ' nbvx=' + string(nbvx) + ',';
  end
  if exists('nbsmooth') then
    str = str + ' nbsmooth=' + string(nbsmooth) + ',';
  end
  if exists('nbjacoby') then
    str = str + ' nbjacoby=' + string(nbjacoby) + ',';
  end
  if exists('ratio') then
    str = str + ' ratio=' + string(ratio) + ',';
  end
  if exists('omega') then
    str = str + ' omega=' + string(omega) + ',';
  end
  if exists('iso') then
    str = str + ' iso=' + string(iso) + ',';
  end
  if exists('abserror') then
    str = str + ' abserror=' + string(abserror) + ',';
  end
  if exists('cutoff') then
    str = str + ' cutoff=' + string(cutoff) + ',';
  end
  if exists('verbosity') then
    str = str + ' verbosity=' + string(verbosity) + ',';
  end
  if exists('inquire') then
    str = str + ' inquire=' + string(inquire) + ',';
  end
  if exists('splitpbedge') then
    str = str + ' splitpbedge=' + string(splitpbedge) + ',';
  end
  if exists('maxsubdiv') then
    str = str + ' maxsubdiv=' + string(maxsubdiv) + ',';
  end
  if exists('anisomax') then
    str = str + ' anisomax=' + string(anisomax) + ',';
  end
  if exists('rescaling') then
    str = str + ' rescaling=' + string(rescaling) + ',';
  end
  
  str = str + '!';
  s   = strsubst(str,',!',';');
  ff_exec(s)
endfunction
