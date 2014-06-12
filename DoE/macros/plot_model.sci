function plot_model(meas_learn, estim_learn, meas_valid, estim_valid)
scf();
learn_resid = meas_learn - estim_learn;
max_learn_resid = max(learn_resid);
min_learn_resid = min(learn_resid);

valid_resid = meas_valid - estim_valid;
max_valid_resid = max(valid_resid);
min_valid_resid = min(valid_resid);

drawlater;

xsetech([0   0   2/3 2/3]); // Regression line
plot(meas_learn, estim_learn, 'k+');
xtitle('Regression line','Measure','Estimation');

xsetech([2/3 0   1/3 1/3]); // Residual distribution
learn_resid = meas_learn - estim_learn;
max_learn_resid = max(learn_resid);
min_learn_resid = min(learn_resid);
x_resid = min_learn_resid:(max_learn_resid-min_learn_resid)/10:max_learn_resid;
histplot(x_resid, learn_resid);
xtitle('Residual distribution','X','Res. amplitude');

xsetech([2/3 1/3 1/3 1/3]); // Textual information
xstring(0.01, 0.01, sprintf('Fisher test : %f',0.0));
xstring(0.01, 1/3+0.01, sprintf('T test : %f', 0.0));
xstring(0.01, 2/3+0.01, sprintf('T test : %f', 0.0));
xtitle('Statistical tests');

xsetech([2/3 2/3-0.01 1/3 1/3]); // L2 for learning / validation
histo_list = sum(learn_resid.^2)
for i=1:size(meas_valid,2)
  histo_list = [histo_list sum((valid_resid(:,i)).^2)];
end
histplot(1+size(meas_valid,2), histo_list);
xtitle('L2 residuals','Data','L2');

xsetech([0   2/3-0.01 1/3 1/3]); // L1 for learning / validation
histo_list = sum(abs(learn_resid))
for i=1:size(meas_valid,2)
  histo_list = [histo_list sum(abs(valid_resid(:,i)))];
end
histplot(1+size(meas_valid,2), histo_list);
xtitle('L1 residuals','Data','L1');

xsetech([1/3 2/3-0.01 1/3 1/3]); // R2 for learning / validation
histplot(1+size(meas_valid,2), [sum((meas_learn - estim_learn).^2) sum((meas_valid - estim_valid).^2)]);
xtitle('R2','Data','R2');
drawnow;
endfunction
