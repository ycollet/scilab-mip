%%testgamma
%fonction de lancement de tests pour le test gamma
%écrite le 29/11/04 par A. Schmied

clear all
close all
rand('state',sum(100*clock))

%variables de test
nb_essais=1000
P=nb_essais/10
dim=20

%essais+bruit
x=rand(nb_essais,dim);

bruit(:,1)=0.75*(2*rand(length(x),1)-1);
s(1)=var(bruit(:,1))

bruit(:,2)=0.05*(2*rand(length(x),1)-1);
s(2)=var(bruit(:,2))

f(:,1)=sum(x,2)+bruit(:,1);
f(:,2)=sum(x,2).^2+0.01.*bruit(:,2);

%Boucle sur le nombre de prototypes
h=waitbar(0,'je calcule...');
for(k=1:P)
    
    [A(k,:),B(k,:)]=gammaTest(x,f,k);
    waitbar(k/P);
    
end;

close(h);
r=1:P
plot(r,A(:,1))
hold on
C=s(1)*ones(P,1)
plot(r,C,'*g')

figure
r=1:P
plot(r,A(:,2))
hold on
C=s(2)*ones(P,1)
plot(r,C,'*g')