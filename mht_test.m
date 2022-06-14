function [] = mht_test(defo,vuw,temperaturefile,datafile,K_a)
%mht_test('ez_limburg_tsx_fp_ps_standard_nl_high_first_10_lines.csv',9,'limburg_tsx_temperature.mat','limburg_tsx_project.mat','K_a.mat')

tic

if nargin==0,help mht_speed_test ;return;end

ts = csvread(defo,1,1);
defo = ts(:,10:end); % ignore the reference epoch, i.e.the first acquisition

if exist(datafile, 'file') == 2;
    load(datafile);
else
    fprintf('please define the datafile');
    fprintf('\n');
    
end

[lin column] = size(defo);


if column == length(Btemp)
    defo = defo*1e3;  % in case the deformation time series is stored in meters, we convert it to millimeter 
else
    defo =defo(:,2:end)*1e3;
end
[lin column] = size(defo);



if exist([project_id '_defo_model_param.csv'],'file') ~= 2;
    fid = fopen([project_id '_defo_model_param.csv'],'w');
    fclose(fid);
elseif exist([project_id '_defo_model_param.csv'],'file') == 2;
    fid = fopen([project_id '_defo_model_param.csv'],'a+');
    fclose(fid);
end

temp_file = temperaturefile;
load(temp_file);

load(K_a);

%Q= diag(Bperp); % cov-var defined by Bperp, note that ratio is better to be applied based on aprior knowledge, option A
Q= vuw*diag(ones(1,length(Bperp))); % uniform weight, option B

% default model
for i = 1:lin
      v(i,:) = inv([Btemp ones(column,1)]'*inv(Q)*[Btemp ones(column,1)])*[Btemp ones(column, 1)]'*inv(Q)*defo(i,:)';
      e_0(:,i) = defo(i,:)' - [Btemp ones(column,1)]*v(i,:)';  
end
 
   %%%%%%correct the error of the reference point%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for i = 1:length(Btemp)
       mean_e(i) = mean(e_0(i,:));
   end
   for i = 1:lin
       defo(i,:) = defo(i,:) - mean_e;
   end
   % recompute the residual of default model
   for i = 1:lin
       v(i,:) = inv([Btemp ones(column,1)]'*inv(Q)*[Btemp ones(column,1)])*[Btemp ones(column, 1)]'*inv(Q)*defo(i,:)';
       e_0(:,i) = defo(i,:)' - [Btemp ones(column,1)]*v(i,:)'; %  
       Q_x = inv([Btemp ones(column,1)]'*inv(Q)*[Btemp ones(column,1)]); % same for all the PS points, as Q is identical, but Q can be customized
       Q_e = Q - [Btemp ones(column,1)]*Q_x*[Btemp ones(column,1)]'; % same for all the PS points
   end
  
for i = 1:lin
    sigma_0p(i) = e_0(:,i)'*inv(Q)*e_0(:,i)/(length(Btemp)-1); % there is one parameter v, if you consider intercept as the second unknown, it should be 2 instead of 1
    T_1(i) = e_0(:,i)'*inv(Q)*e_0(:,i);
    %   alpha = 1/2/length(Btemp);
    %   Kritical = chi2inv(1-alpha,length(Btemp)-1);
end

% Ling's matrix, same for all PS points for the same Ha
% v + eta
% Cy for Btemp + temperature
Ling2 = inv(Q)*temp*inv(temp'*inv(Q)*Q_e*inv(Q)*temp)*temp'*inv(Q);

M_7_temp = ones([length(Btemp)-1 length(Btemp)-1]);
M_7_temp = tril(M_7_temp);
M_7 = [zeros(1,length(Btemp)-1); M_7_temp];
for i = 1:length(Btemp)-1
    % v + delta
    % Cy for Btemp + delta
    % heavyside, one offset delta location test
    Ling4(i,:,:) = inv(Q)*M_7(:,i)*inv(M_7(:,i)'*inv(Q)*Q_e*inv(Q)*M_7(:,i))*M_7(:,i)'*inv(Q);
    % v + eta + delta
    % Cy for Btemp + temperature + delta
    Ling3(i,:,:) = inv(Q)*[temp M_7(:,i)]*inv([temp M_7(:,i)]'*inv(Q)*Q_e*inv(Q)*[temp M_7(:,i)])*[temp M_7(:,i)]'*inv(Q);
end

Tr_q1 = zeros(lin,1);
alpha = K_a(1);

%%%%%%add (one) breakpoint model  %%%%%
% assume the velocity difference is the additional parameter (v_new - v)
% H_0: y = [t1 t2]v + e; H_j: y = [t1 t2]v + [0 t2](v_new - v) + e
for i = 1:length(Btemp)-3
    M_bp(i,:,:) = [zeros(1,i+1) (Btemp(i+2:end)-Btemp(i+1))'];
Ling5(i,:,:) = inv(Q)* reshape(M_bp(i,:,:),length(Btemp),1)*inv(reshape(M_bp(i,:,:),length(Btemp),1)'*inv(Q)*Q_e*inv(Q)*reshape(M_bp(i,:,:),length(Btemp),1))*reshape(M_bp(i,:,:),length(Btemp),1)'*inv(Q);
end
%%%%%%add (one) breakpoint model %%%%%

for i = 1:lin
    %if T_1(i) < K_a(2)%Kritical
    if    T_1(i) < K_a(end)
        fid = fopen([project_id '_defo_model_param.csv'],'a+');
        para = [i; 1; v(i); 0; 0; 0; 0; 0; 0];
        fprintf(fid,'%12.0f, %12.0f, %12.8f, %12.8f, %12.8f, %12.8f, %12.8f,%12.8f, %12.0f\n',para);
        fclose(fid);
        %      clc;
        %      fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        %      fprintf(['T_0 = ' num2str(T_1(i)) ', K = ' num2str(Kritical(i))]);
        %      fprintf('\n');
        %      fprintf('%%%%%%%%%%%%%%%%%%%% go to alternatives %%%%%%%%%%%%%%%%%%%%');
        %
    elseif     T_1(i) > K_a(end) %T_1(i) > K_a(2)%Kritical
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % v + eta
        % Cy for Btemp + temperature
        crivalue_new = K_a(2);
        T_q1(i) = e_0(:,i)'*Ling2*e_0(:,i);
        % ratio test
        Tr_q1(i) = T_q1(i)/crivalue_new;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % v + delta
        for p = 1:length(Btemp)-1
            T_q3(p) = e_0(:,i)'*reshape(Ling4(p,:,:),length(Btemp),length(Btemp))*e_0(:,i);
            % ratio test
            Tr_q3(p) = T_q3(p)/crivalue_new;
        end
        i_Tr_q3 = find(Tr_q3 == max(Tr_q3));
         %%%%%%2015-09-30%%%%%%%%%%%%%%%%%%% 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % v + delata_v breakpoint
        for p = 1:length(Btemp)-3
            T_q5(p) = e_0(:,i)'*reshape(Ling5(p,:,:),length(Btemp),length(Btemp))*e_0(:,i);
            % ratio test
            Tr_q5(p) = T_q5(p)/crivalue_new;
        end
        i_Tr_q5 = find(Tr_q5 == max(Tr_q5));
        %%%%%%2015-09-30%%%%%%%%%%%%%%%%%%% 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % v + eta + delta
        crivalue_new = K_a(3);
        for w = 1:length(Btemp)-1
            T_q2(w) = e_0(:,i)'*reshape(Ling3(w,:,:),length(Btemp),length(Btemp))*e_0(:,i);
            % ratio test
            Tr_q2(w) = T_q2(w)/crivalue_new;
        end
        i_Tr_q2 = find(Tr_q2 == max(Tr_q2));
        
        %%%%%%%%%%%%%%%%%%%%%% 2014-09-10
        % Tr_a = [Tr_q1(i) Tr_q2(i_Tr_q2) Tr_q3(i_Tr_q3)];
        % add breakpoint model
         Tr_a = [Tr_q1(i) Tr_q2(i_Tr_q2) Tr_q3(i_Tr_q3) Tr_q5(i_Tr_q5)];
        Tr_q5(i_Tr_q5)
         imax = find(Tr_a == max(Tr_a));
        if Tr_a(imax) > 1
            
        % %%%%%%%%%%%%%%%%%%% 2014-09-10 %%%%%%%%%%%%
        doplot= 'n';
            switch imax
                case 1
                    [x1 e1] =  defovtemper(defo(i,:),vuw,doplot,temp,Bperp,Btemp,dates,Q,lambda);
                    fid = fopen([project_id '_defo_model_param.csv'],'a+');
                    x1 = [i; 2; v(i); x1(1);x1(2); 0; x1(3); 0; 0];
                    fprintf(fid,'%12.0f, %12.0f, %12.8f, %12.8f, %12.8f, %12.8f, %12.8f, %12.8f, %12.0f\n',x1);
                    fclose(fid);
                case 2
                    [x2 e2] = defovtemperdelt(defo(i,:),vuw,doplot,temp,Bperp,Btemp,dates,i_Tr_q2,Q,lambda);
                    x2 = [i; 3; v(i); x2(1);x2(3);x2(2); x2(4); 0; i_Tr_q2];
                    fid = fopen([project_id '_defo_model_param.csv'],'a+');
                    fprintf(fid,'%12.0f, %12.0f, %12.8f, %12.8f, %12.8f, %12.8f, %12.8f, %12.8f, %12.0f\n',x2);
                    fclose(fid);
                case 3
                    [x3 e3] = defovdelt(defo(i,:),vuw,doplot,Bperp,Btemp,dates,i_Tr_q3,Q,lambda);
                    x3 = [i; 4; v(i); x3(1); 0; x3(2); x3(3); 0; i_Tr_q3];
                    fid = fopen([project_id '_defo_model_param.csv'],'a+');
                    fprintf(fid,'%12.0f, %12.0f, %12.8f, %12.8f, %12.8f, %12.8f, %12.8f, %12.8f, %12.0f\n',x3);
                    fclose(fid);
                case 4  %%%%%%breakpoint model
                    [x4 e4] =defobreakpoint(defo(i,:),vuw,doplot,Bperp,Btemp,dates,i_Tr_q5,Q,lambda);
                    x4 = [i; 5; v(i); x4; 0; 0; 0; i_Tr_q5];
                    fid = fopen([project_id '_defo_model_param.csv'],'a+');
                    fprintf(fid,'%12.0f, %12.0f, %12.8f, %12.8f, %12.8f, %12.8f, %12.8f, %12.8f, %12.0f\n',x4);
                    fclose(fid);
            end
        elseif Tr_a(imax) <= 1
            fid = fopen([project_id '_defo_model_param.csv'],'a+');
            para = [i; 1; v(i); 0; 0; 0; 0; 0; 0];
            fprintf(fid,'%12.0f, %12.0f, %12.8f, %12.8f, %12.8f, %12.8f, %12.8f, %12.8f, %12.0f\n',para);
            fclose(fid);
        end
    end
end
    toc
    
    
 