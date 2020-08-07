function kModelOut = learnKrig(kModelIn, Plot, StepDistance, Horizon, NbRestart)
// Definition of the kriging structure
// kModel('model')    : name of the model
// kModel('tendency') : list of the tendency parameters
// kModel('var')      : list of the variables
// kModel('data')     : learning data set
// kModel('K')        : K matrix
// kModel('param')    : parameters of the covariance

if ~isdef('Horizon','local') then
  Horizon = 1;
end
if ~isdef('NbRestart','local') then
  NbRestart = 10;
end
if ~isdef('StepDistance','local') then
  StepDistance = 1;
end

NbPoints = size(kModelIn('data'),1);

//
// Computation of the variogram
//

if (Plot) then
  WId = waitbar(0, 'Learning the model');
end

// Data to be estimated
X = kModelIn('data')(:,1:$-1);
Y = kModelIn('data')(:,$);

// Variogram : semi-variance

Z       = [];
V       = []; // Contains the values of the variogram
C       = []; // Contains the values of the co-variance
delta   = []; // Contains the distances
MaxDist = 0;

ieee_flag = ieee();
ieee(2);

// Computation of the matrix of distances
// Building of a list of distances between points Distance(i,j) - F(i) - F(j)
List_Dist = zeros(NbPoints, 3);
Index = 0;
for i=1:NbPoints
  if (Plot) & ~modulo(NbPoints, floor(NbPoints/10)) then
    waitbar(min([1 floor(1000*i/(NbPoints-1))/1000]),WId);
  end

  for j=1:NbPoints
    Index = Index + 1;
    List_Dist(Index,1) = norm(X(j,:) - X(i,:));
    List_Dist(Index,2) = Y(j);
    List_Dist(Index,3) = Y(i);
    if (MaxDist<List_Dist(Index,1)) then
      MaxDist = List_Dist(Index,1);
    end
  end
end

List_Dist = lex_sort(List_Dist,1,'unique');
//List_Dist = lex_sort(List_Dist,1);

// Computation of the semi-variance and the co-variance

Index      = 1;
Index_List = 0;
while(Index~=size(List_Dist,1))
  if (Plot) then
    waitbar(min([1 floor(1000*Index/size(List_Dist,1))/1000]),WId);
  end

  Count  = 0;
  Count2 = 0;
  SemVa  = 0;
  CoVa   = 0;
  while((abs(List_Dist(Index,1)-List_Dist(Index+Count,1))<MaxDist*StepDistance) & ...
        ((Index+Count)~=size(List_Dist,1)) & ...
        (List_Dist(Index+Count,1)<MaxDist*Horizon))
    SemVa = SemVa + (List_Dist(Index+Count,2) - List_Dist(Index+Count,3))^2;
    CoVa  = CoVa  + (List_Dist(Index+Count,2) * List_Dist(Index+Count,3));
    Count = Count + 1;
  end
  // We skip the remaining points
  while((abs(List_Dist(Index,1)-List_Dist(Index+Count2,1))<MaxDist*StepDistance) & ...
        ((Index+Count+Count2)~=size(List_Dist,1)))
    Count2 = Count2 + 1;
  end
  if (Count~=0) then
    Index_List = Index_List + 1;
    V(Index_List) = SemVa / (2*(Count+1));
    C(Index_List) = CoVa / ((Count+1));
    delta(Index_List) = List_Dist(Index,1);
    Index = Index + Count + Count2;
  else
    Index = Index + Count + Count2;
  end
end

if (Plot) then
  hg = scf();
  set(hg,'figure_name','Krigeage - Learning - Plot %d');
  drawlater;  
  subplot(2,1,1);
  plot(delta,C,'r');
  plot(delta,V,'g');
  plot(delta,C+V,'k'); // YC: a revoir
  xtitle('Variogramme','x','Amplitude');
  legend(['Co-Variance','Semi-variance','Sum Cov+Semi-Var']);
end

Z=[C,delta];

kModelOut = kModelIn;

// Learning of covariance parameters via restarting
MinError    = %inf;
Min_par_sol = [];
for i=1:NbRestart
  execstr('[parm_sol,err] = learn_correl_' + kModelIn('model') + '(Z'')');
  if (err<MinError) then
    MinError     = err;
    Min_parm_sol = parm_sol;
    Parm_Index   = i;
  end
  if Plot then
    printf('learnKrig: Restart number %d / %d - learning error = %f - Min learning error = %f\n', i, NbRestart, err, MinError);
  end
end

if Plot then
  printf('learnKrig: parameters set %d chosen\n',Parm_Index);
end

parm_sol = Min_parm_sol;

kModelOut('param') = parm_sol;

//
// Computation of the K matrix
//

kModelOut('K') = zeros(NbPoints, NbPoints);

for i=1:NbPoints
  if (Plot) & ~modulo(i,floor(NbPoints/10)) then
    waitbar(min([1 floor(1000*i/NbPoints)/1000]),WId);
  end
  for j=i:NbPoints       
    dist_XiXj = norm(X(i,:) - X(j,:));
    dist_XiXj = max([%eps dist_XiXj]);
      
    execstr('kModelOut(''K'')(i, j) = correl_' + kModelOut('model') + '(dist_XiXj,parm_sol)');
    kModelOut('K')(j,i) = kModelOut('K')(i,j);
  end
end

if (Plot) then
  winclose(WId);
  subplot(2,1,2);
  execstr('Y_correl = correl_' + kModelIn('model') + '(delta,parm_sol)');
  
  plot(delta, Y_correl, 'r', delta, Z(:,1), 'g', delta, Y_correl-Z(:,1), 'k');

  xtitle('Plot of the covariance function','x','Covariance');
  legend(['Model' 'Covariance' 'Error']);
  drawnow;
  printf('Error mean     = %f\n', mean(Y_correl-Z(:,1)));
  printf('Error variance = %f\n', stdev(Y_correl-Z(:,1)));
end

ieee(ieee_flag);
endfunction

