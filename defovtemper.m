function [x e] = defovtemper(defo,vuw,doplot,temp,Bperp,Btemp,dates,Q,lambda)
% defo model base 15
% y = v*t + temperature*eta + e
% e.g. [x e Tr_q1] = defovtemper(defo,1,'y','temperature_test.mat')
% 2013-05-20
% no need nonlinear LSE
% temp_file = temperaturefile;
% load(temp_file);

if nargin==0,help  defovtemper;return;end
if nargin==1, doplot = 'n'          ;end

[lin colum] = size(defo);
if lin == 1
    defo = defo;
else if colum == 1
        defo =defo';
    end
end
    

D = [Btemp temp ones(colum,1)]; %  

x = inv(D'*inv(Q)*D)*D'*inv(Q)*defo';
e = defo' - D*x;
% posterior variance
sigma_2 = e'*inv(Q/vuw)*e/(length(Btemp)-2);
% data snooping 2014-05-08
 Q_x = inv(D'*inv(Q)*D);
 Q_e = Q - D*Q_x*D';
%  for i = 1:length(Btemp)
%      w(i) = e(i)/Q_e(i,i);
%  end

if strcmp(doplot,'y'),
    f = figure;
    %set(f,'name','defo base model 2: v + eta','numbertitle','off');
    %figure(1);
    %errorbar(Btemp,defo,e,'r');
    plot(Btemp,defo,'r+');
    hold on
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
    xlabel(['Time   v = ' num2str(x(1)) ', \eta = ' num2str(x(2)) ', \sigma^2 = ' num2str(sigma_2)]);
    ylabel('deformation [mm]');
   % title('defo base function 15: v + temperature');
    %title('defo base function 2: B_t * v +temperature * \eta');
    axis tight
    

    
end

