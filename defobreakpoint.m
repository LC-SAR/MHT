function [x e] = defobreakpoint(defo,vuw,doplot,Bperp,Btemp,dates,i_Tr_q5,Q,lambda)
% 2015-09-30
if nargin==0,help defobreakpoint ;return;end
if nargin==1, doplot = 'n'          ;end

[lin colum] = size(defo);
if lin == 1
    defo = defo;
else if colum == 1
        defo =defo';
    end
end


%D = [Btemp(1:i_Tr_q5+1)' zeros(1,length(Btemp)-1-i_Tr_q5); zeros(1,i_Tr_q5+1) Btemp(i_Tr_q5+2:end)']';
dp = [zeros(1, i_Tr_q5+1) (Btemp(i_Tr_q5+2:end)- Btemp(i_Tr_q5+1))']';
D = [Btemp dp]
x = inv(D'*inv(Q)*D)*D'*inv(Q)*defo';
% error
e = defo' - D*x;
% posterior variance
sigma_2 = e'*inv(Q/vuw)*e/(length(Btemp)-2);

