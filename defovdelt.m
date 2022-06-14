function [x e] = defovdelt(defo,vuw,doplot,Bperp,Btemp,dates,i_Tr_q3,Q,lambda)
% defo model base 2
% y = b_t * v + Delta * heavyside + e
% e.g. [x e Tr_q3] = defovdelt(defo,3,'y')
% 2013-04-15
if nargin==0,help defovdelt ;return;end
if nargin==1, doplot = 'n'          ;end

[lin colum] = size(defo);
if lin == 1
    defo = defo;
else if colum == 1
        defo =defo';
    end
end
    
 
% heavyside, one offset delta location test
M_7_temp = ones([length(Btemp)-1 length(Btemp)-1]);
M_7_temp = tril(M_7_temp);
M_7 = [zeros(1,length(Btemp)-1); M_7_temp];

% 2014-02-22
M_7 = M_7(:,i_Tr_q3);
D = [Btemp M_7 ones(colum, 1)];
x = inv(D'*inv(Q)*D)*D'*inv(Q)*defo';
% error
e = defo' - D*x;
% posterior variance
sigma_2 = e'*inv(Q/vuw)*e/(length(Btemp)-2);
% 2014-02-22




% data snooping 2014-05-08
 Q_x = inv(D'*inv(Q)*D);
 Q_e = Q - D*Q_x*D';
%  for i = 1:length(Btemp)
%      w(i) = e(i)/Q_e(i,i);
%  end

 if strcmp(doplot,'y'),
     f = figure;
    %set(f,'name','defo base model 4: v + Delta','numbertitle','off');
    %errorbar(Btemp,defo,e,'r');
    plot(Btemp,defo,'r+-');
    hold on
   % exam = find(sigma_2 == min(sigma_2)); %show given delta number
    exam = i_Tr_q3;
    %plot(Btemp,D{exam}*x(:,exam),'-');
    plot(Btemp,D*x,'k-');
    % 2014-05-08 plot the ambiguity
    hold on
    plot(Btemp,D*x+1000*lambda/4,'g-')
    hold on
    plot(Btemp,D*x-1000*lambda/4,'b-')
    legend('obs','est',['+' '\pi'],['-' '\pi']);
    % 2014-05-08
    %legend('obs','est');
    DATES=datenum(dates(1:end-1,:),'dd-mmm-yyyy');
    BBB2=[Btemp(1):1:Btemp(end)];
    BBB=[DATES(1):365:DATES(end)];  
    set(gca, 'XTick',BBB2);
    set(gca, 'XTickLabel', datestr(BBB,'yyyy/mm'));
    grid on
    %xlabel('Time');
    xlabel(['Time   v = ' num2str(x(1)) ', \Delta = '  num2str(x(2)) ', at epoch ' num2str(i_Tr_q3), ', \sigma^2 = ' num2str(sigma_2)]);
    ylabel('deformation [mm]');
    %title('defo base function 4: B_t * v + t*Delta');
    axis tight
    

 end
  
