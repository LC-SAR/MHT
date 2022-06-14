function [] = K_alpha(datafile)
% compute alpha and critical value and store them
% [] = K_alpha('limburg_tsx_project.mat')
% 2014-06-09
% to compute q = 1, 2, 3, 4
% updated on 17-Feb-2015
load(datafile);

alpha = 1/2/length(Btemp);
Kritical = chi2inv(1-alpha,length(Btemp)-1);

gamma = 0.5; % power, 0.8 or 0.5
[lambda1] = lambda0est(alpha,1,gamma); % this is the first initial value when alpha is alpha_0
% search new alpha


q = 1;
crivalue_new1=ncx2inv(1-gamma,q,lambda1);
q= 2;
crivalue_new2=ncx2inv(1-gamma,q,lambda1);
q=3;
crivalue_new3=ncx2inv(1-gamma,q,lambda1);
q=4;
crivalue_new4=ncx2inv(1-gamma,q,lambda1);
q = length(Btemp)-1;
crivalue_new5=ncx2inv(1-gamma,q,lambda1);

K_a = [alpha crivalue_new1 crivalue_new2 crivalue_new3 crivalue_new4 crivalue_new5];
save('K_a.mat','K_a');
