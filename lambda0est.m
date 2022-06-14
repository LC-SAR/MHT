%2014-07-25
%estimate the lambda_0 based on gamma alpha and q
function [lambda] = lambda0est(a,q,g)
% a alpha, g gamma
% e.g.  lambda0est(0.001,1,0.8);lambda0est(0.001,1,[0.05:0.1:0.95]) %
% 17.075 lambda0est(0.001,1,0.8)
% a = 0.001;
% q = 1;
% g = 0.8;
crtvlue=chi2inv(1-a,q)
beta0= 1 - g;
lambda = nan(1,length(g));
Err=5E-5;
for i = 1:length(g)
    
    
    low=0; high=100;
    mid=(low+high)/2;
    while high-low>Err 
        if ncx2cdf(crtvlue,q,mid)>beta0(i)
           low=mid;
           mid=(high+low)/2;
        elseif ncx2cdf(crtvlue,q,mid)<beta0(i)
           high=mid;
           mid=(high+low)/2;
        else
           high=mid;
           low=mid; 
        end
 
    end
        lambda(i) = mid;
end

%figure;plot(g,lambda,'k');ylabel('\lambda','fontsize',12);xlabel('\gamma','fontsize',12);
%grid on
%axis tight;
%title(['\alpha = ' num2str(a) ', q =' num2str(q)],'fontsize',12);
%set(gca,'fontsize',12)
