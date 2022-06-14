function [x e] = defovtemperdelt(defo,vuw,doplot,temp,Bperp,Btemp,dates,i_Tr_q2,Q,lambda)
% defo model base 17
% y = v*t + temperature*eta + Delta * heavyside + e
% e.g. [x e Tr_q2] = defovtemperdelt(defo,3,'y','temperature_test.mat')
% 2013-05-20
% no need nonlinear LSE
% temp_file = temperaturefile;
% load(temp_file);

if nargin==0,help  defovtemperdelt;return;end
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
M_7 = M_7(:,i_Tr_q2);
D = [Btemp M_7 temp ones(colum,1)];
x = inv(D'*inv(Q)*D)*D'*inv(Q)*defo';
e = defo' - D*x;
sigma_2 = e'*inv(Q/vuw)*e/(length(Btemp)-3);


% data snooping 2014-05-08
 Q_x = inv(D'*inv(Q)*D);
 Q_e = Q - D*Q_x*D';
%  for i = 1:length(Btemp)
%      w(i) = e(i)/Q_e(i,i);
%  end

 if strcmp(doplot,'y'),
    f = figure;
    %set(f,'name','defo base model 3: v + eta + Delta','numbertitle','off');
    %errorbar(Btemp,defo,e,'r');
    plot(Btemp,defo,'r+');
    hold on
    %exam = find(sigma_2 == min(sigma_2)); %show given delta number
    exam = i_Tr_q2;
    plot(Btemp,D*x,'k-');
    % 2014-05-08 plot the ambiguity
    hold on
    plot(Btemp,D*x+1000*lambda/4,'g-')
    hold on
    plot(Btemp,D*x-1000*lambda/4,'b-')
    legend('obs','est',['+' '\pi'],['-' '\pi']);
    hold on
    plot([Btemp(i_Tr_q2) Btemp(i_Tr_q2)],[min(defo) max(defo)],'k--');
    % 2014-05-08
    %legend('obs','est');
    DATES=datenum(dates(1:end-1,:),'dd-mmm-yyyy');
    BBB2=[Btemp(1):1:Btemp(end)];
    BBB=[DATES(1):365:DATES(end)];  
    set(gca, 'XTick',BBB2);
    set(gca, 'XTickLabel', datestr(BBB,'yyyy/mm'));
    grid on
    xlabel(['Time   v = ' num2str(x(1)) ', \eta = ' num2str(x(3)), ', \Delta = '  num2str(x(2)) ', at epoch ' num2str(i_Tr_q2) ', \sigma^2 = ' num2str(sigma_2)]);
    ylabel('deformation [mm]');
    %title('defo base function 3: B_t * v +temperature * eta + t*Delta');
   % title('defo base function 17: v + temperature + offset_{heavyside}');
    axis tight
    

 end
  
 
 
