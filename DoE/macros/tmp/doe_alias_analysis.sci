function permut_list = gen_permute(nb_var, resolution)
// This function generate a list of resolution symbols in a list of nb_var symbols
// We count in binary and if there are resolution bits at one, we select the corresponding symbols as a permutation
permut_list = [];
Count = zeros(1,nb_var);
AuxVar = 1:nb_var;
permut_list = list();
for i=1:2^nb_var -1
  Count(1) = Count(1) + 1;
  for j=2:nb_var
    if (Count(j-1)>1) then
      Count(j)   = Count(j) + 1;
      Count(j-1) = 0;
    end
  end
  if (sum(Count)==resolution) then
    I = find(Count==1);
    permut_list($+1) = AuxVar(I);
  end
end
endfunction

function alias_list = doe_alias_analysis(model, nb_var, resolution)
if (resolution>=nb_var) then
  error('doe_alias_analysis: resolution must be strictly inferior to nb_var');
end
// This alias analysis is only dedicated to 2-levels DoE with interactions
// With this alias generator, we have to find all the alias to a given interaction
// I = A.B.C.D
// First, locate all the interactions in the model
// Second, build the alias to a given interaction
// Third, find if this aliased interaction is present

// Generation of the list of alias generator
AliasGenerator = gen_permute(nb_var, resolution);

// Localisation of the interaction monomials
Index = 0;
ListMonom = tokens(model,' ');
for i=1:size(ListMonom,1)
  row_p = grep(ListMonom(i),'^'); // If there is a power sign in the monom, it's not an interaction
  row_m = grep(ListMonom(i),'*'); // If there is a multiply sign in the monom, it's an interaction
  if (isempty(row_p)&~isempty(row_m)) then // It's an interaction
    Index = Index + 1;
    ListInteraction(Index) = ListMonom(i);
  end
end

// Transformation of the interaction monomials into an integer representation
Alias = 1:resolution;
for i=1:size(ListInteraction,1)
  ListInteraction(i) = strsubst(ListInteraction(i), '*', ' ');
  ListInteraction(i) = strsubst(ListInteraction(i), 'x', '');
end

// Construction of the alias to a given interaction
ListAlias = list();
for i=1:size(ListInteraction,1)
  // If the size of the alias generator is smaller than the size of the interaction, we skip this comparison
  AuxListInter = eval(tokens(ListInteraction(i),' '));
  if (size(AliasGenerator)<=size(AuxListInter,2)) then continue; end
  
  for j=1:size(AliasGenerator)
    // We have something like [1 4], we transform this into [1 0 0 4] to get easily the complementary
    AuxLI = zeros(1, nb_var);
    AuxLI(AuxListInter) = AuxListInter';
    AuxAG = zeros(1, nb_var);
    AuxAG(AliasGenerator(j)) = AliasGenerator(j);
    Res = find((~(AuxAG&AuxLI)).*(1:nb_var)>0);
    if (~isempty(Res)) then
      ListAlias($+1) = Res;
    end
  end
end

// Reconstruction of the integer form into a string representation
ListAlias2 = [];
for i=1:size(ListAlias)
  AuxString = '';
  for j=1:size(ListAlias(i),2)
    AuxString = AuxString + 'x' + string(ListAlias(i)(j))+'*';
  end
  AuxString = part(AuxString, 1:length(AuxString)-1); // We remove the last '*'
  ListAlias2(i) = AuxString;
end
clear ListAlias;

// Initialisation of the output structure
// This structure will be a list. The first element is the alias generator and the following elements are
// the aliases identified in the model
alias_list = list();
for i=1:size(ListAlias2,1)
  alias_list($+1) = list(ListAlias2(i));
end

// Locating the alias
for i=1:size(ListAlias2,1)
  row = grep(ListMonom, ListAlias2(i));
  if (~isempty(row)) then
    for j=1:size(row,2)
      alias_list(i)($+1) = ListMonom(row(j));
    end
  end
end
endfunction

