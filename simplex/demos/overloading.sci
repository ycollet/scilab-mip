function x = ask(this)
  execstr('x = %' + typeof(this) + '_ask(this)')
endfunction

function this = tell(this, x, y)
  execstr('this = %' + typeof(this) + '_tell(this, x, y)')
endfunction

function [yopt, xopt] = best(this)
  execstr('[yopt, xopt] = %' + typeof(this) + '_best(this)')
endfunction
