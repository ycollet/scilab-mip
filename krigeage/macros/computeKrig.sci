function [Estim, Var] = computeKrig(modelList, vector_value)
  // Definition of the kriging structure
  // modelList('model')    : name of the model
  // modelList('tendency') : list of the tendency parameters
  // modelList('vars')     : list of the variables
  // modelList('data')     : learning data set
  // modelList('K')        : K matrix
  // modelList('param')    : parameters of the covariance

  DataMeas = modelList('data');
  K = modelList('K');
  X = DataMeas(:,1:$-1);
  Y = DataMeas(:,$);
  NbPoints = size(DataMeas,1);
  
  // We split the tendancy into a list of terms which will be used to build the matrixes P and p_x
  ListOfTerms = tokens(modelList('tendency'),' ');
  ListOfVar   = tokens(modelList('var'),' ');
  
  Model    = modelList('model');
  parm_sol = modelList('param');
  
  ieee_flag = ieee();
  ieee(2);
      
  // Building matrix P
  // We replace all the le xj by X(i,j)
  for i=1:size(ListOfTerms,1)
    for j=1:size(ListOfVar,1)
      ListOfTerms(i) = strsubst(ListOfTerms(i),ListOfVar(j),'X(i,' + string(j) + ')');
    end
  end

  P = zeros(NbPoints,size(ListOfTerms,1));
  for i=1:NbPoints
    for j=1:size(ListOfTerms,1)
      P(i,j) = eval(ListOfTerms(j));
    end
  end

  x = vector_value;

  // Building tendency for vector x
  // We replace all the le xj by x(j)
  for i=1:size(ListOfTerms,1)
    for j=1:size(ListOfVar,1)
      x_Tendency(i) = eval(strsubst(ListOfTerms(i),ListOfVar(j),'x(' + string(j) + ')'));
    end
  end
    
  p_x = eval(ListOfTerms);

  Index_xi = 0;
  for i=1:NbPoints
    Index_xi = Index_xi + 1;
    dist_xXi = (x - X(i,:))*(x - X(i,:))';
    dist_xXi = sqrt(dist_xXi);
      
    execstr('k_x(Index_xi) = correl_' + Model + '(dist_xXi,parm_sol)');
  end
  
  null_mat=zeros(size(ListOfTerms,1), size(ListOfTerms,1));
  
  big_mat=[K  , P;
           P' , null_mat];

  small_mat=[k_x;p_x];
  weight = inv(big_mat)*small_mat;
  
  execstr('Var = correl_' + Model + '(0,parm_sol)');
  
  Estim = 0.0;
  
  for i=1:NbPoints
    Estim = Estim + weight(i)*Y(i);
    Var   = Var - weight(i)*k_x(i);
  end

  // Ajout de la tendance a  l'estimation.
  for i=1:size(ListOfTerms,1)
    Var = Var - weight(NbPoints+i)*x_Tendency(i);
  end
  ieee(ieee_flag);
  Var = abs(Var); //YC: a revoir. Il doit y avoir un bug dans la formule de la variance
endfunction

