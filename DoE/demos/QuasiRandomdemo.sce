my_handle = scf(100001);
clf(my_handle,'reset');
demo_viewCode('QuasiRandomdemo.sce');

lines(0);
old_funcprot = funcprot();
funcprot(0);

Dimension = 2;
NbPts     = 100;

List_Halton     = doe_halton(Dimension, NbPts);
List_Hammersley = doe_hammersley(Dimension, NbPts);
List_Faure      = doe_faure(Dimension, NbPts);

subplot(2,2,1);
plot(List_Halton(:,1), List_Halton(:,2), 'k.');
xtitle('Halton random sequence','x1','x2');
subplot(2,2,2);
plot(List_Hammersley(:,1), List_Hammersley(:,2), 'k.');
xtitle('Hammersley random sequence','x1','x2');
subplot(2,2,3);
plot(List_Faure(:,1), List_Faure(:,2), 'k.');
xtitle('Faure random sequence','x1','x2');

my_handle2 = scf(100002);
clf(my_handle2,'reset');

// Now plot a 3D cloud of points

Dimension = 3;
NbPts     = 100;

List_Halton     = doe_halton(Dimension, NbPts);
List_Hammersley = doe_hammersley(Dimension, NbPts);
List_Faure      = doe_faure(Dimension, NbPts);

subplot(2,2,1);
param3d1(List_Halton(:,1), List_Halton(:,2), list(List_Halton(:,3),-1));
e = gce();
e.mark_size = 2;
xtitle('Halton random sequence','x1','x2','x3');

subplot(2,2,2);
param3d1(List_Hammersley(:,1), List_Hammersley(:,2), list(List_Hammersley(:,3),-1));
e = gce();
e.mark_size = 2;
xtitle('Hammersley random sequence','x1','x2','x3');

subplot(2,2,3);
param3d1(List_Faure(:,1), List_Faure(:,2), list(List_Faure(:,3),-1));
e = gce();
e.mark_size = 2;
xtitle('Faure random sequence','x1','x2','x3');

funcprot(old_funcprot);

