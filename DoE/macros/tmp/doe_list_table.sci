function List = doe_list_table(type_doe, nb_exp, nb_var, nb_level)
// type_doe : can be 
// - ssd (super saturated design)
// - oa (orthogonal array)
// - Had (hadamard array)
// - ma (mixed array)
// nb_exp : the number of experiments
// nb_var : the number of variables in the doe
// nb_level : the number of levels per variable

if (~isdef('type_doe','local')) then
  type_doe = '*';
end
if (~isdef('nb_exp','local')) then
  nb_exp = '*';
end
if (~isdef('nb_var','local')) then
  nb_var = '*';
end
if (~isdef('nb_level','local')) then
  nb_level = '*';
end

File_WildCard = type_doe + '_' + nb_exp + 'r_' + nb_var + 'v_' + nb_level + 'l*';

List = list_files(File_WildCard);
endfunction
