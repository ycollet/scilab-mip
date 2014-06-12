NbParam = 10;

//selmeth = "aic";
selmeth = "datavalid";

ListOfFiles = dir('data_mult_lasso_??.dat');

for i=1:size(ListOfFiles(2),1)
  Data = [];
  Data = read(ListOfFiles(2)(i),-1,2);
  ListOfMean_Yhat(i) = mean(Data(:,2));
  ListOfVar_Yhat(i)  = stdev(Data(:,2));
  ListOfR2(i)        = 1 - sum((Data(:,1) - Data(:,2)).^2) / sum((Data(:,1)-mean(Data(:,1))).^2);
end

X = 1:size(ListOfFiles(2),1);

Param = read('data_mult_lasso_listofparam.dat',-1,NbParam);
ListOfColors = int(rand(3,NbParam)*255);

SelFunc = read('data_mult_lasso_selparam.dat',-1,3);

Bound = SelFunc(:,1);

drawlater;
subplot(2,2,1);
plot(X, ListOfMean_Yhat, 'k', X, ListOfMean_Yhat + ListOfVar_Yhat, 'r.-', X, ListOfMean_Yhat - ListOfVar_Yhat, 'r.-');
xtitle('Estimation','N° pb', 'Mean / Var');
legend(['Mean', 'Mean+Std','Mean-Std']);

subplot(2,2,2);
plot(X, max(min(ListOfR2,1),0), 'k');
xtitle('R2','N° pb', 'R2');

subplot(2,2,3);
for i=1:NbParam
  plot2d(X, Param(:,i),style=color(ListOfColors(1,i), ListOfColors(2,i), ListOfColors(3,i)));
end
xtitle('Evolution of the parameters','N° pb','Bound');

subplot(2,2,4);
if (selmeth=="aic") then
  plot(SelFunc(:,1), SelFunc(:,2), 'k');
end
if (selmeth=="datavalid") then
  plot(SelFunc(:,1), SelFunc(:,3), 'k');
end
xtitle('Value of the selection parameter','Bound', 'FObj');
drawnow;

