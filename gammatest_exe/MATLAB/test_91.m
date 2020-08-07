%%testgamma
%fonction de lancement de tests pour le test gamma
%écrite le 29/11/04 par A. Schmied

clear all
close all
rand('state',sum(100*clock))



%essais
x=load('C:\Documents and Settings\a066399.CORP\Desktop\GammaTest\données\91 collage\FACTEURS.txt');
f=load('C:\Documents and Settings\a066399.CORP\Desktop\GammaTest\données\91 collage\REPONSES.txt');

temp1=size(f)




temp=size(x)
%variables de test
nb_essais=temp(1,1)
P=40
dim=temp(1,2)

%Boucle sur le nombre de prototypes
h=waitbar(0,'je calcule...');
for(k=5:P)
    
    [A(k-4,:),B(k-4,:)]=gammaTest(x,f,k);
    waitbar(k/P);
    
end;

close(h);
for i=1:temp1(1,2)
    figure
    r=1:P-4
    plot(r,B(:,i))
end;
A_mean=mean(B,1)