my_handle = scf(100001);
clf(my_handle,'reset');
demo_viewCode('lars_test.sce');

lines(0);
old_funcprot = funcprot();
funcprot(0);

// Example script for performing LASSO regression on the diabetes data set.

loadmatfile('diabetes.mat');

X = diabetes.x; size(X)
X = normalize(X);
y = diabetes.y;
y = center(y);
[n p] = size(X);

b1 = lars(X, y, 'lasso', 0, 0, [], 1);
s1 = sum(abs(b1),2)/sum(abs(b1(size(b1,1),:)));

plot(s1, b1, '-');
plot(s1, b1, 'ro');

xtitle('LARS','lambda','param');

[s_opt, b_opt, res_mean, res_std] = crossvalidate(lars, 10, 1000, X, y, 'lasso', 0, 0, [], 0);

cvplot(s_opt, res_mean, res_std);

funcprot(old_funcprot);

