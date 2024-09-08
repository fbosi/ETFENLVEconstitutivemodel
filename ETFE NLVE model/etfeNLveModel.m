%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script simulates ETFE behavior according to the isotropic or the 
% orthotropic model produced in the LIGHTEN campaign - www.lighten-itn.eu
%
% Written by: A. Comitti, 2024 ~  a.comitti@ucl.ac.uk
%
% Copyright F. Bosi, 2024; A. Comitti, 2024; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear;
% path = 'C:\Users\Alessandro Comitti\OneDrive - University College London\Documenti personali\PhD\tests\DATA\MATLAB OVERALL SCRIPTS';
path = pwd;
%%%%%%%%%%%% CHOOSE THE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = "ISO"; %choose between "ISO" and "ORTHO"

%%%%%%%%%%%% CHOOSE THE FORMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

formulation = "strain"; %choose between "stress" to obtain stresses,
% "strains" in order to get strains

%%%%%%%%%%%% loading the model constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if model == "ORTHO"
    opt1 = readmatrix(strcat(path,'\D11_Opt_1e5Mpa.csv'));
    opt2 = readmatrix(strcat(path,'\D22_Opt_1e5Mpa.csv'));
    opt6 = readmatrix(strcat(path,'\D66_Opt_1e5Mpa.csv'));
    ea11 = readmatrix(strcat(path,'\EA11_1e5Mpa.csv'));    
    ea22 = readmatrix(strcat(path,'\EA22_1e5Mpa.csv'));
    ea66 = readmatrix(strcat(path,'\EA66_1e5Mpa.csv'));
    activationvolume = 2.74;
else
    opt1 = readmatrix(strcat(path,'\ISO.csv'));
    opt2 = readmatrix(strcat(path,'\ISO.csv'));
    opt6 = readmatrix(strcat(path,'\ISO.csv'));
    ea11 = readmatrix(strcat(path,'\EAISO.csv'));    
    ea22 = readmatrix(strcat(path,'\EAISO.csv'));
    ea66 = readmatrix(strcat(path,'\EAISO.csv'));
    activationvolume = readmatrix(strcat(path,'\240304_EyringISO2.csv'));
end
cte = readmatrix(strcat(path,'\CTE.csv'));
poisson = 0.43;
TRef = 20 + 273.15;
ea11 = ea11/3 + ea22/3 + ea66/3; %I average the three values
factor = 2.303*8.31446261815324;



%% %%%%%%% TEST DATA LOADING; comment out the test not needed %%%%%%%%%%%%%
%%%% the raw data are all expressed in engineering quantities %%%%%%%%%%%%

%%% in order to simulate the example test data, make sure to install the 
%%% Curve Fitting Toolbox: 
%%% https://www.mathworks.com/products/curvefitting.html %%%%%


%1.
% %%%% uniaxial constant strain rate tests data from Instron, MD, T = 313 K
% rawData = rmmissing(readmatrix(strcat(path,'\OUTRaw_MD_CSR_01_T40.csv')));
% tempData = rmmissing(readmatrix( ...
%     strcat(path,'\OUTTemp_MD_CSR_01_T40.csv')));
% % temperature vector synchronisation
% for i = 1: size(tempData,2)-1
%     T{i} = rmmissing(interp1(tempData(:,1)-1,tempData(:,i+1), ...
%         rawData(:,1+15*(i-1)),'linear','extrap'))+273.15; 
% end
% 
% % interpolation on a strain vector
% e11 = transpose(0:0.001:0.02);
% for i = 1: size(tempData,2)-1
%     [e11C{i}, index,~] = unique(rawData(:,2+15*(i-1)),'stable');
%     trefC(:,i) = interp1(e11C{i},rawData(index,1+15*(i-1)),e11,'linear');
%     e22C(:,i) =  interp1(e11C{i},rawData(index,5+15*(i-1)),e11,'linear');
%     e12C(:,i) = zeros(length(e11),1);
%     s11C(:,i) = interp1(e11C{i},rawData(index,4+15*(i-1)),e11,'linear');
%     s22C(:,i) = zeros(length(e11),1);
%     s12C(:,i) = zeros(length(e11),1);
%     tempC(:,i) = interp1(e11C{i},T{i}(index),e11,'linear');
% end
% tref = mean(trefC,2,'omitnan'); 
% e22 = mean(e22C,2,'omitnan'); 
% e12 = zeros(length(e11),1);
% s11 = mean(s11C,2,'omitnan'); 
% s22 = zeros(length(e11),1);
% s12 = zeros(length(e11),1);
% temp = mean(tempC,2,'omitnan'); 
% s11Sd =  std(s11C,0,2,'omitnan'); 
% e11Sd =  zeros(length(e11),1);
% seqSd = s11Sd;
% eeqSd = e11Sd;
% e33 = -poisson/(1-poisson).*(e11+e22);
% eeq = 1/(1+poisson).*(1/2.*((e11-e22).^2+(e11-e33).^2+ ...
%     (e33-e22).^2+6*(e12).^2)).^0.5;
% seq = (1/2*((s11-s22).^2+(s22).^2+(s11).^2+6*(s12).^2)).^0.5;



%2.
%%%% uniaxial relaxation tests data from Instron, MD, T = 298 K %%%%%%%%%%%
% rawData = rmmissing(readmatrix(strcat( ...
%     path,'\OUTRaw_MD_Rel_9MPa_T20.csv')));
% tempData = rmmissing(readmatrix(strcat( ...
%     path,'\OUTTemp_MD_Rel_9MPa_T20.csv')));
% 
% % temperature vector synchronisation
% for i = 1: size(tempData,2)-1
%     T{i} = rmmissing(interp1(tempData(:,1),tempData(:,i+1), ...
%         rawData(:,1+9*(i-1)),'linear',tempData(end,i+1)))+273.15; 
% end
% 
% % interpolation on a time vector
% tref = transpose(logspace(-1,4,40000));
% for i = 1: size(tempData,2)-1
%     [t{i},index,~] = unique(rawData(:,1+9*(i-1)),'stable');
%     e11C(:,i)  =  interp1(t{i}, ...
%         smooth(rawData(:,2+9*(i-1)),10),tref,'nearest');
%     e22C(:,i) = interp1(t{i},rawData(:,4+9*(i-1)),tref,'nearest');
%     e12C(:,i) = zeros(length(tref),1);
%     s11C(:,i) = interp1(t{i}, ...
%         smooth(rawData(:,3+9*(i-1)),10),tref,'nearest');
%     s22C(:,i) = zeros(length(tref),1);
%     e12C(:,i) = zeros(length(tref),1);
%     tempC(:,i)  = interp1(t{i},T{i}(index),tref,'nearest');
% end    
% e11 = mean(e11C,2,'omitnan'); 
% e22 = mean(e22C,2,'omitnan'); 
% s11 = mean(s11C,2,'omitnan'); 
% temp = mean(tempC,2,'omitnan'); 
% s11Sd =  std(s11C,0,2,'omitnan');
% e11Sd =  std(e11C,0,2,'omitnan');
% cut = 17;
% e11 = e11(17:end);
% e22 = e22(17:end);
% s11 = s11(17:end);
% temp = temp(17:end);
% s11Sd = s11Sd(17:end);
% e11Sd = e11Sd(17:end);
% tref = tref(17:end)-tref(17);
% e12 = zeros(length(e11),1);
% s22 = zeros(length(e11),1);
% s12 = zeros(length(e11),1);
% seqSd = s11Sd;
% eeqSd = e11Sd;
% e33 = -poisson/(1-poisson).*(e11+e22);
% eeq = 1/(1+poisson).*(1/2.*((e11-e22).^2+(e11-e33).^2+ ...
%     (e33-e22).^2+6*(e12).^2)).^0.5;
% seq = (1/2*((s11-s22).^2+(s22).^2+(s11).^2+6*(s12).^2)).^0.5;



%3.
%%%%%%%%%%%%% uniaxial cyclic data from Instron, MD, T = 295 K %%%%%%%%%%%% 
% rawData = rmmissing(readmatrix(strcat( ...
%     path,'\OUTRaw_MD_CYC_CSR_01_T20.csv')));
% tempData = rmmissing(readmatrix(strcat( ...
%     path,'\OUTTemp_MD_CYC_CSR_01_T20.csv')));
% for i = 1: size(tempData,2)-1
%     T{i} = rmmissing(interp1(tempData(:,1)-1,tempData(:,i+1), ...
%             rawData(:,1+8*(i-1)),'linear','extrap'))+273.15; 
% end
% 
% %%%%%%%%choose the specimen, between 1 and 5
% z = 1;
% for i = z
%     e11 = rawData(:,3+8*(i-1));
%     tref = rawData(:,1+8*(i-1));
%     e22 =  rawData(:,6+8*(i-1));
%     e12 = rawData(:,8+8*(i-1));
%     s11 = rawData(:,4+8*(i-1));
%     s22 = zeros(length(e11),1);
%     s12 = zeros(length(e11),1);
%     temp = T{i};   
% end
% s11Sd = zeros(length(e11),1);
% e11Sd = zeros(length(e11),1);
% seqSd = s11Sd;
% eeqSd = e11Sd;
% e33 = -poisson/(1-poisson).*(e11+e22);
% eeq = 1/(1+poisson).*(1/2.*((e11-e22).^2+(e11-e33).^2+ ...
%     (e33-e22).^2+6*(e12).^2)).^0.5;
% seq = (1/2*((s11-s22).^2+(s22).^2+(s11).^2+6*(s12).^2)).^0.5;




%4.
%%%%%%%%%%%%% biaxial inflation on elliptical specimen, T = 298K %%%%%%%%%% 
% rawData = readmatrix(strcat(path,"\OUTRaw_ELLI_TD_T25_SR001.csv"));
% e11fit = 0:0.0001:0.015;
% e22MAT = [];
% s11MAT = [];
% s22MAT = [];
% seqMAT = [];
% e22RMAT = [];
% e11RMAT = [];
% TMAT=[];
% tMAT = [];
% 
% for i = 1:size(rawData,2)/13
%     eunique=[];
%     eequnique = [];
%     t{i} = rmmissing(rawData(:,1+(i-1)*13));
%     TC{i} = rmmissing(rawData(:,13+(i-1)*13));
%     e11C{i}=rmmissing(rawData(:,5+(i-1)*13));
%     e22C{i}=rmmissing(rawData(:,6+(i-1)*13));
%     s11C{i}=rmmissing(rawData(:,7+(i-1)*13));
%     s22C{i} =rmmissing(rawData(:,8+(i-1)*13));
% if length(s11C{i})<length(e11C{i})
%     t{i} = rmmissing(rawData(1:length(s11C{i}),1+(i-1)*13));
%     TC{i} = rmmissing(rawData(1:length(s11C{i}),13+(i-1)*13));
%     e11C{i}=rmmissing(rawData(1:length(s11C{i}),5+(i-1)*13));
%     e22C{i}=rmmissing(rawData(1:length(s11C{i}),6+(i-1)*13));
% end
%     seqC{i} = ((s11C{i}).^2+(s22C{i}).^2-(s11C{i}.*s22C{i})).^0.5;
%     [~,eunique] = unique(e11C{i});
%     tRi = interp1(e11C{i}(eunique),t{i}(eunique),e11fit,'linear');
%     TRi = interp1(e11C{i}(eunique),TC{i}(eunique),e11fit,'linear');
%     e22i = interp1(e11C{i}(eunique),e22C{i}(eunique),e11fit,'linear');
%     s11i = interp1(e11C{i}(eunique),s11C{i}(eunique),e11fit,'linear');
%     s22i = interp1(e11C{i}(eunique),s22C{i}(eunique),e11fit,'linear');
%     tMAT = [tMAT tRi'];
%     TMAT = [TMAT TRi'];
%     e22MAT = [e22MAT e22i'];
%     s11MAT = [s11MAT s11i'];
%     s22MAT = [s22MAT s22i'];
%     seqi = interp1(e11C{i}(eunique),seqC{i}(eunique),e11fit,'linear');
%     seqMAT = [seqMAT seqi'];
% end
% e11 = e11fit';
% tref = mean(tMAT,2);
% temp = mean(TMAT,2)+273.15; 
% e22 = mean(e22MAT,2);
% s11 = mean(s11MAT,2);
% s22 = mean(s22MAT,2);
% e12 = zeros(length(e11),1);
% s12 = zeros(length(e11),1);
% s11Std = std(s11MAT,[],2);
% s22Std = std(s22MAT,[],2);
% e33 = -poisson/(1-poisson).*(e11+e22);
% eeq = 1/(1+poisson).*(1/2.*((e11-e22).^2+(e11-e33).^2+ ...
%     (e33-e22).^2+6*(e12).^2)).^0.5;
% seq = mean(seqMAT,2);
% seqSd= std(seqMAT,[],2);
% eeqSd= zeros(length(eeq),1);



%5 - 6 - 7. 
%%%%%%%%%%%%%% uniaxial DMA data along MD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%5.
%%%%%%%%%%%%%% creep test at 9 MPa, 273.15K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rawData = readmatrix(strcat(path,'\OUTRaw_MD_Creep_9MPa_T00.csv'));
% tlim = 1.825750000000000e+04;
% [~,indd] = min(abs(tlim-rawData(:,1)));
% rawData = rawData(1:indd,:);
% rawData(:,2:6:18) = rawData(:,2:6:18)*1e-2;

%6.
%%%%%%%%%%%%%% relaxation test at 6 MPa, 333.15K %%%%%%%%%%%%%%%%%%%%%%%%%%
% rawData = readmatrix(strcat(path,'\OUTRaw_MD_Relax_6MPa_T60.csv'));
% tlim = 1.825750000000000e+04;
% [~,indd] = min(abs(tlim-rawData(:,1)));
% rawData = rawData(1:indd,:);

%7.
%%%%%%%%%%%%%% creep test at 7.5 MPa, 293.15K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawData1 = readmatrix(strcat(path,'\T_20_S_1200_Sample_01_CreepRecovery.csv'));
rawData2 = readmatrix(strcat(path,'\T_20_S_1200_Sample_02_CreepRecovery.csv'));
rawData3 = readmatrix(strcat(path,'\T_20_S_1200_Sample_03_CreepRecovery.csv'));
LLL = max([size(rawData1,1),size(rawData2,1),size(rawData3,1)]);
rawData = [[rawData1(:,1:6);nan(LLL-size(rawData1,1),6)], [rawData2(:,1:6);nan(LLL-size(rawData2,1),6)],[rawData3(:,1:6);nan(LLL-size(rawData3,1),6)]];
for i = 1:length(rawData)
    for k = 1:3
        if rawData(i,2+6*(k-1))>0.1
            rawData(i,2+6*(k-1))=nan(1,1);
        end
    end
end

rawData = rmmissing(rawData);
TempData = [ones(size(rawData,1),1) rawData(:,6:6:18)];
% temperature vector DMA
for  i = 1: size(TempData,2)-1
    T{i} = TempData(:,i+1)+273.15; 
end
tref = rawData(:,1);
for i = 1: size(TempData,2)-1
    [t{i},index,~] = unique(rawData(:,1+6*(i-1)),'stable');
    eC(:,i)  =  interp1(t{i},rawData(index,2+6*(i-1)),tref,'nearest');%*1e-2;
    e22C(:,i) = -0.43*eC(:,i);
    s11C(:,i) = interp1(t{i},rawData(index,3+6*(i-1)),tref,'nearest');
    TempC(:,i)  = interp1(t{i},T{i}(index),tref,'nearest');
end    
e11 = rmmissing((mean(eC,2,'omitnan'))); 
e22 = rmmissing(mean(e22C,2,'omitnan')); 
e12 = zeros(length(e11),1);
s11 = rmmissing(mean(s11C,2,'omitnan')); 
s22 = zeros(length(e11),1);
s12 = zeros(length(e11),1);
temp = rmmissing(mean(TempC,2,'omitnan')); 
e11Sd =  rmmissing(std(eC,0,2,'omitnan')); 
s11Sd =  rmmissing(std(s11C,0,2,'omitnan')); 
seqSd = s11Sd;
eeqSd = e11Sd;
e33 = -poisson/(1-poisson).*(e11+e22);
eeq = 1/(1+poisson).*(1/2.*((e11-e22).^2+(e11-e33).^2+ ...
    (e33-e22).^2+6*(e12).^2)).^0.5;
seq = (1/2*((s11-s22).^2+(s22).^2+(s11).^2+6*(s12).^2)).^0.5;





%% %%%%%%%%%%% Recursive algorithm integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if model == "ORTHO"
    for q = 1:length(opt1)-1
            tau(q)=opt1(q+1,2);
            D110=opt1(1,1);
            D11j(q)=opt1(q+1,1);
            D220=opt2(1,1); 
            D22j(q)=opt2(q+1,1);
            D120=-poisson*D110;
            D12j(q)=-poisson*D11j(q);
            D210= D120; 
            D21j(q)=D12j(q);
            D130 =  D120;
            D13j(q) = D12j(q);
            D230 =  D120;
            D23j(q) = D12j(q);
            D660 =  opt6(1,1);
            D66j(q) = opt6(q+1,1);
    end
else
    for q = 1:length(opt1)-1
        tau(q)=opt1(q+1,2);
        D110=opt1(1,1);
        D11j(q)=opt1(q+1,1);
        D220=opt2(1,1); 
        D22j(q)=opt2(q+1,1);
        D120=-poisson*(D110+D220)/2;
        D12j(q)=-poisson*(D11j(q)+D22j(q))/2;
        D210= D120; 
        D21j(q)=D12j(q);
        D130 =  D120;
        D13j(q) = D12j(q);
        D230 =  D120;
        D23j(q) = D12j(q);
        D660 =  opt6(1,1)*(1+poisson);
        D66j(q) = opt6(q+1,1)*(1+poisson);
    end
end

%%%%% initialisation: I initialise the values for both the formulations
%simulated stresses and related quantities
s11_R = zeros(length(e11),1) ;% engineering stress MD  
s22_R = zeros(length(e11),1) ;% engineering stress TD
s12_R = zeros(length(e11),1) ;% engineering stress ID (in-plane shear)
s_R = zeros(length(e11),1) ;% engineering equivalent stress
%simulated stresses and related quantities
e11_R = zeros(length(e11),1) ;% engineering strain MD  
e22_R = zeros(length(e11),1) ;% engineering strain TD
e12_R = zeros(length(e11),1) ;% engineering strain ID (in-plane shear)
e_R = zeros(length(e11),1) ;% engineering equivalent stress
%common quantities
e33_R = zeros(length(e11),1); % out of plane engineering strain
%thermal strains
etherm1 = 0;
etherm2 = 0;
etherm3 = 0;
%mechanical strains
e11M = zeros(length(e11),1);
e22M = zeros(length(e11),1);
e33M_R = zeros(length(e11),1);
%hereditary integrals
q11 = zeros(length(tau),1) ; 
q11old = zeros(length(tau),1) ;
q12 = zeros(length(tau),1) ;
q12old = zeros(length(tau),1) ;
q21 = zeros(length(tau),1) ;
q21old = zeros(length(tau),1) ;
q22 = zeros(length(tau),1) ;
q22old = zeros(length(tau),1) ;
q13 = zeros(length(tau),1) ;
q13old = zeros(length(tau),1) ;
q23 = zeros(length(tau),1) ;
q23old = zeros(length(tau),1) ;
q66 = zeros(length(tau),1) ;
q66old = zeros(length(tau),1) ;
%maximum of the equivalent stress and first value for Eyring stress
s_Rmax = 0;
sEyring = 0;
%threshold for the in plane strain rate
DeDt_thresh = 0;
%strain rate equivalent and in plane
vE = zeros(length(e11),3);
eD = 0;
eD2 = 0;
DeDt = 0;
De2Dt = 0;


if formulation == "stress"
    for n = 2:1:length(e11)
            %temperature shift factor 
            aT = 10^(-ea11/factor*(-1/temp(n) + 1/TRef));
            %stress shift factor
            %for the UMAT, I introduced a stress rate needed to detect the
            %unloading. 
            if n == 2
                eD = e_R(n-1);
                eD2 = e_R(n-1);
            else
                eD = (e_R(n-1) - e_R(n-2));
                vE(n,:) = [e11(n), e22(n),e12(n)];
                [maxEps,index] = max(vE(n),[],2);
                eD2 = vE(n,index) - vE(n-1,index);
            end
            DeDt = eD/((tref(n)-tref(n-1)));
            De2Dt = eD2/((tref(n)-tref(n-1)));
            %starting condition for the experiments: 
            % the sense of this is to avoid unstable values of the shift 
            % factor at the first step with zero stresses
            if n == 2  || n == 3 %I consider the first 2 steps
                aSigma = 1;
            else

            %%% Freezing of the value used to calculate the shift factor in
            %%% case of unloading detected

            if De2Dt<DeDt_thresh
                sEyring = s_Rmax;
            else
                sEyring = s_R(n-1);                          
            end 

            if sEyring == 0
                aSigma=1;
            end
            %here I compute the argument of the hyperbolic sine of the
            %Eyring function, and if it is too small, I approximate it with
            %its taylor expansion 
            x = (activationvolume*(sEyring)/(8.31446262*1e-3*temp(n)));
            arg = sinh(x);
            if arg<1e-10
                arg = x + 1/6*x^3 + 1/120*x^5;
            end
                aSigma = (activationvolume*sEyring/(8.31446262*1e-3* ...
                    temp(n)*arg));
            end
            %reduced time
            dtr = (tref(n)-tref(n-1))/(aT*aSigma);
            storage(n) = dtr;
            %Initialization of parameters
            F11 = 0;        
            D11 = D110;
            F12 = 0;        
            D12 = D120;
            F21 = 0;        
            D21 = D210;
            F22 = 0;        
            D22 = D220;
            F13 = 0;
            D13 = D130;
            F23 = 0;
            D23 = D230; 
            F66 = 0;
            D66 = D660; 
            %hereditary storage 
            q11old = q11;
            q12old = q12;
            q21old = q21;
            q22old = q22;
            q13old = q13;
            q23old = q23;
            q66old = q66;
            %looping through the Prony series
            for c = 1:length(tau)
                dtrexp = dtr/tau(c);
                q = dtrexp;
                if dtrexp<1e-10
                    g = 1-dtrexp/2+(dtrexp)^2/6-(dtrexp)^3/24;
                else
                    g = (1-exp(-dtrexp))/dtrexp;
                end
                D11 = D11 + D11j(c)*(1-g);
                D12 = D12 + D12j(c)*(1-g);
                D22 = D22 + D22j(c)*(1-g);
                D21 = D21 + D21j(c)*(1-g);
                D13 = D13 + D13j(c)*(1-g);
                D23 = D23 + D23j(c)*(1-g);
                D66 = D66+ D66j(c)*(1-g);
                F11 = F11 + D11j(c)*(exp(-dtrexp)*q11old(c)-g*s11_R(n-1));
                F12 = F12 + D12j(c)*(exp(-dtrexp)*q12old(c)-g*s22_R(n-1));
                F21 = F21 + D21j(c)*(exp(-dtrexp)*q21old(c)-g*s11_R(n-1));
                F22 = F22 + D22j(c)*(exp(-dtrexp)*q22old(c)-g*s22_R(n-1));
                F13 = F13 + D13j(c)*(exp(-dtrexp)*q13old(c)-g*s11_R(n-1));
                F23 = F23 + D23j(c)*(exp(-dtrexp)*q23old(c)-g*s22_R(n-1));
                F66 = F66 + D66j(c)*(exp(-dtrexp)*q66old(c)-g*s12_R(n-1)); 
            end
    %obtaining the mechanical strains
            etherm1 = polyval(cte(1,:),temp(n)) - polyval(cte(1,:),temp(1));
            etherm2 = polyval(cte(2,:),temp(n)) - polyval(cte(2,:),temp(1));
            etherm3 = polyval(cte(3,:),temp(n)) - polyval(cte(3,:),temp(1));
            e11M(n) = e11(n)- etherm1;
            e22M(n) = e22(n)- etherm2;
    
    %calculation of stresses 
            detD = D11*D22-D21*D12;
            s11_R(n) = 1/detD *(D22*(e11M(n)+F11+F12)-D12*(e22M(n)+(F21+F22)));
            s22_R(n) = 1/detD *(D11*(e22M(n)+F22+F21)-D21*(e11M(n)+(F12+F11)));
            s12_R(n) = (e12(n)+F66)/D66;
            s_R(n) = (1/2*((s11_R(n)-s22_R(n)).^2+(s22_R(n)).^2+(s11_R(n)).^2+6*(s12_R(n)).^2)).^0.5;
    % calculation of mechanical and total out of plane strain
            e33M_R(n) = D13*s11_R(n) +  D23*s22_R(n)-F13-F23; 
            e33_R(n) = e33M_R(n) + etherm3;
            e_R(n) = 1/(1+poisson).*(1/2.*((e11M(n)-e22M(n)).^2+(e11M(n)-e33M_R(n)).^2+(e33M_R(n)-e22M(n)).^2+6*(e12(n)).^2)).^0.5;
            DeDt = (e_R(n)-e_R(n-1))/(tref(n)-tref(n-1));
            %update of the Eyring variable that contains the max VM stress
            % for the loading step.
            if De2Dt<DeDt_thresh
                s_Rmax = s_Rmax;
            else
                s_Rmax = s_R(n);                          
            end      

    %calculation of the hereditary integral at the end of the step
            for c = 1:length(tau)
                dtrexp = dtr/tau(c);
                if dtrexp<1e-10
                    g = 1-dtrexp/2+(dtrexp)^2/6-(dtrexp)^3/24;
                else
                    g = (1-exp(-dtrexp))/dtrexp;
                end
                q11(c) = exp(-dtrexp)*q11old(c)+(s11_R(n)-s11_R(n-1))*g;
                q12(c) = exp(-dtrexp)*q12old(c)+(s22_R(n)-s22_R(n-1))*g;        
                q21(c) = exp(-dtrexp)*q21old(c)+(s11_R(n)-s11_R(n-1))*g;
                q22(c) = exp(-dtrexp)*q22old(c)+(s22_R(n)-s22_R(n-1))*g;
                q13(c) = exp(-dtrexp)*q13old(c)+(s11_R(n)-s11_R(n-1))*g;        
                q23(c) = exp(-dtrexp)*q23old(c)+(s22_R(n)-s22_R(n-1))*g;   
                q66(c) = exp(-dtrexp)*q66old(c)+(s12_R(n)-s12_R(n-1))*g; 
            end
    end

else
    for n = 2:1:length(e11)
            %temperature shift factor 
            aT = 10^(-ea11/factor*(-1/temp(n) + 1/TRef));
            %stress shift factor       
            if n == 2
                eD = e_R(n-1);
                eD2 = e_R(n-1);
            else
                eD = (e_R(n-1) - e_R(n-2));
                vE(n,:) = [e11_R(n-1), e22_R(n-1),e12_R(n-1)];
                [maxEps,index] = max(vE(n),[],2);
                eD2 = vE(n,index) - vE(n-1,index);
            end
            DeDt = eD/((tref(n)-tref(n-1)));
            De2Dt = eD2/((tref(n)-tref(n-1)));
            %starting condition
            if n == 2  || n == 3 
                aSigma = 1;
            else
            %%% switch in loading made to acount for loading and
            %%% unloading. 
            if De2Dt<DeDt_thresh
                sEyring = s_Rmax;
            else
                sEyring = seq(n-1);                          
            end          
            %here I compute the argument of the hyperbolic sine of the
            %eyring function, and if it is too small, I approximate it with
            %its taylor expansion 
            x = (activationvolume*(sEyring)/(8.31446262*1e-3*temp(n)));
            arg = sinh(x);
            if arg<1e-10
                arg = x + 1/6*x^3 + 1/120*x^5;
            end
                aSigma = (activationvolume*sEyring/(8.31446262*1e-3 ...
                    *temp(n)*arg));
            end
            if sEyring == 0
                aSigma = 1;
            end
            %reduced time
            dtr = (tref(n)-tref(n-1))/(aT*aSigma);
            %Initialization of parameters
            F11 = 0;        
            D11 = D110;
            F12 = 0;        
            D12 = D120;
            F21 = 0;        
            D21 = D210;
            F22 = 0;        
            D22 = D220;
            F13 = 0;
            D13 = D130;
            F23 = 0;
            D23 = D230; 
            F66 = 0;
            D66 = D660; 
            %hereditary storage 
            q11old = q11;
            q12old = q12;
            q21old = q21;
            q22old = q22;
            q13old = q13;
            q23old = q23;
            q66old = q66;
            %looping through the Prony series
            for c = 1:length(tau)
                dtrexp = dtr/tau(c);
                q = dtrexp;
                if dtrexp<1e-10
                    g = 1-dtrexp/2+(dtrexp)^2/6-(dtrexp)^3/24;
                else
                    g = (1-exp(-dtrexp))/dtrexp;
                end
                D11 = D11 + D11j(c)*(1-g);
                D12 = D12 + D12j(c)*(1-g);
                D22 = D22 + D22j(c)*(1-g);
                D21 = D21 + D21j(c)*(1-g);
                D13 = D13 + D13j(c)*(1-g);
                D23 = D23 + D23j(c)*(1-g);
                D66 = D66+ D66j(c)*(1-g);
                F11 = F11 + D11j(c)*(exp(-dtrexp)*q11old(c)-g*s11(n-1));
                F12 = F12 + D12j(c)*(exp(-dtrexp)*q12old(c)-g*s22(n-1));
                F21 = F21 + D21j(c)*(exp(-dtrexp)*q21old(c)-g*s11(n-1));
                F22 = F22 + D22j(c)*(exp(-dtrexp)*q22old(c)-g*s22(n-1));
                F13 = F13 + D13j(c)*(exp(-dtrexp)*q13old(c)-g*s11(n-1));
                F23 = F23 + D23j(c)*(exp(-dtrexp)*q23old(c)-g*s22(n-1));
                F66 = F66 + D66j(c)*(exp(-dtrexp)*q66old(c)-g*s12(n-1)); 
            end
    %thermal component of strains
            etherm1 = polyval(cte(1,:),temp(n)) - polyval(cte(1,:),temp(1));
            etherm2 = polyval(cte(2,:),temp(n)) - polyval(cte(2,:),temp(1));
            etherm3 = polyval(cte(3,:),temp(n)) - polyval(cte(3,:),temp(1));
    
    %calculation of mechanical strains 
            detD = D11*D22-D21*D12;
            e11M_R(n) = D11*s11(n)+D12*s22(n)-F11-F12;
            e22M_R(n) = D21*s11(n)+D22*s22(n)-F22-F21;
            e33M_R(n) = D13*s11(n) +  D23*s22(n)-F13-F23;
            e12_R(n) = D66*s12(n) -  F66;
    %calculation of total strains  
            e11_R(n) =  e11M_R(n) + etherm1;
            e22_R(n) =  e22M_R(n) + etherm2;
            e33_R(n) =  e33M_R(n) + etherm3;
            e_R(n) = 1/(1+poisson).*(1/2.*((e11M_R(n)-e22M_R(n)).^2+(e11M_R(n)-e33M_R(n)).^2+(e33M_R(n)-e22M_R(n)).^2+6*(e12_R(n)).^2)).^0.5;
            DeDt = (e_R(n)-e_R(n-1))/(tref(n)-tref(n-1));
            vE(n,:) = [e11_R(n), e22_R(n),e33_R(n)];
            [maxEps,index] = max(vE(n),[],2);
            eD2 = vE(n,index) - vE(n-1,index);
            De2Dt = eD2/((tref(n)-tref(n-1)));
            %update of the Eyring variable that contains the max VM stress for
            %the loading step.
            if De2Dt<DeDt_thresh
                s_Rmax = s_Rmax;
            else
                s_Rmax = seq(n);                          
            end      
    
    %calculation of the hereditary integral at the end of the step
            for c = 1:length(tau)
                dtrexp = dtr/tau(c);
                if dtrexp<1e-10
                    g = 1-dtrexp/2+(dtrexp)^2/6-(dtrexp)^3/24;
                else
                    g = (1-exp(-dtrexp))/dtrexp;
                end
                q11(c) = exp(-dtrexp)*q11old(c)+(s11(n)-s11(n-1))*g;
                q12(c) = exp(-dtrexp)*q12old(c)+(s22(n)-s22(n-1))*g;        
                q21(c) = exp(-dtrexp)*q21old(c)+(s11(n)-s11(n-1))*g;
                q22(c) = exp(-dtrexp)*q22old(c)+(s22(n)-s22(n-1))*g;
                q13(c) = exp(-dtrexp)*q13old(c)+(s11(n)-s11(n-1))*g;        
                q23(c) = exp(-dtrexp)*q23old(c)+(s22(n)-s22(n-1))*g;   
                q66(c) = exp(-dtrexp)*q66old(c)+(s12(n)-s12(n-1))*g; 
            end
    end
end




%% %%% Stress - strain plot %%%%%%
fig = figure();  
Col = ['k','r','b','m','g'];
Blue = 1/255*[150,150,245];
Orange = 1/255*[245,150,150];
plot(eeq,seq,'LineStyle', '-','Color',Blue,'LineWidth',0.5, ...
    'DisplayName', strcat('Experimental'))
hold on
erev = [eeq' flip(eeq')];
inbetween = [seq-seqSd;flip(seqSd+seq)];
h1=  fill(erev, inbetween, Blue,'FaceAlpha',0.5, 'EdgeAlpha',0.2, ...
    "DisplayName", "Experimental");
if formulation == "stress"
    h2 =   plot(eeq,s_R,'LineStyle','-.','Color',Col(4),'LineWidth',2, ...
        'DisplayName',strcat('Model'));
else
    h2 =   plot(e_R,seq,'LineStyle','-.','Color',Col(4),'LineWidth',2, ...
        'DisplayName',strcat('Model'));
end
hold off
xlim([0 max(eeq)*1.2])
ylim([0 max(seq)*1.2])
set(gca,'fontname','Times New Roman') 
set(gcf, 'color', [1 1 1])
fontsize(gca,16,"points")
fig.Units = 'centimeters';
set(gcf,'Position',[3 3 19 15.2])
xlabel("$\epsilon$",'FontName','Times New Roman','FontSize',18, ...
    'Interpreter','latex')
ylabel("$\sigma$ [MPa]",'FontName','Times New Roman','FontSize',18, ...
    'Interpreter','latex')
hold off
    
lgd = legend([h1 h2]);
lgd.NumColumns = 1;
lgd.Location = "northwest";



%% %%%%%% Time plot %%%%%%%%%%
fig = figure();
if formulation == "stress"
    yyaxis right
    h1 =  semilogx(tref,eeq,'LineStyle', ':','Marker','none','LineWidth',2, ...
        'DisplayName',strcat('\epsilon_{input - experimental}'))  ;
    hold all
    trev = [tref+1e-10 ;flip(tref+1e-10)]';
    inbetween = [(eeq-eeqSd);flip(eeqSd +eeq)];
    h4 = fill(trev, inbetween, Orange,'FaceAlpha',0.5, 'EdgeAlpha',0.2, ...
        'Marker','none',"DisplayName", "S. Deviation");
    ylabel("$\epsilon$",'FontName','Times New Roman','FontSize', ...
        18,'Interpreter','latex')
    xlim([1e-2, max(tref)*1.2])
    ylim([0 max(eeq)*2])
    
    yyaxis left
    h2 = semilogx(tref,seq,'LineStyle', '-','Marker','none','LineWidth', ...
        1,'Color',Col(1),'DisplayName',strcat('\sigma_{experimental}'))  ;  
    h3 =  semilogx(tref,s_R,'LineStyle', ':','Marker','none', ...
        'LineWidth',2,'Color',Col(3),'DisplayName',strcat('\sigma_{model}'));
    inbetween = [(seq-seqSd);flip(seqSd+seq)];
    h4 = fill(trev, inbetween, Blue,'FaceAlpha',0.5, 'EdgeAlpha',0.2, ...
        'Marker','none',"DisplayName", "S. Deviation");
    hold off
    ylim([0 max(seq)*1.2])
    set(gca,'fontname','Times New Roman')
    fontsize(gca,16,"points")
    set(fig,'defaultAxesColorOrder',[[0 0 1]; [1 1 1]]);
    set(gcf, 'color', [1 1 1])
    fig.Units = 'centimeters';
    set(gcf,'Position',[3 3 19 15.2])
    xlabel("t [s]",'FontName','Times New Roman','FontSize',18, ...
        'Interpreter','latex')
    ylabel("$\sigma$ [MPa]",'FontName','Times New Roman','FontSize', ...
        18,'Interpreter','latex')
    hold off
else
    yyaxis right
    h1 =  semilogx(tref,eeq,'LineStyle', ':','Marker','none','LineWidth',2, ...
        'DisplayName',strcat('\epsilon_{experimental}'))  ;
    hold all
    trev = [tref+1e-10 ;flip(tref+1e-10)]';
    inbetween = [(eeq-eeqSd);flip(eeqSd +eeq)];
    h4 = fill(trev, inbetween, Orange,'FaceAlpha',0.5, 'EdgeAlpha',0.2, ...
        'Marker','none',"DisplayName", "S. Deviation");
    h3 =  semilogx(tref,e_R,'LineStyle', ':','Marker','none', ...
        'LineWidth',2,'Color',Col(3),'DisplayName',strcat('\epsilon_{model}'));
    ylabel("$\epsilon$",'FontName','Times New Roman','FontSize', ...
        18,'Interpreter','latex')
    xlim([1.9, max(tref)*1.2])
    ylim([0 max(eeq)*2])
    yyaxis left
    h2 = semilogx(tref,seq,'LineStyle', '-','Marker','none','LineWidth', ...
        1,'Color',Col(1),'DisplayName',strcat('\sigma_{input - experimental}'))  ;  
    inbetween = [(seq-seqSd);flip(seqSd+seq)];
    h4 = fill(trev, inbetween, Blue,'FaceAlpha',0.5, 'EdgeAlpha',0.2, ...
        'Marker','none',"DisplayName", "S. Deviation");
    hold off
    ylim([0 max(s11)*1.2])
    set(gca,'fontname','Times New Roman')
    fontsize(gca,16,"points")
    set(fig,'defaultAxesColorOrder',[[0 0 1]; [1 1 1]]);
    set(gcf, 'color', [1 1 1])
    fig.Units = 'centimeters';
    set(gcf,'Position',[3 3 19 15.2])
    xlabel("t [s]",'FontName','Times New Roman','FontSize',18, ...
        'Interpreter','latex')
    ylabel("$\sigma$ [MPa]",'FontName','Times New Roman','FontSize', ...
        18,'Interpreter','latex')
    hold off
end
    
lgd = legend([h1 h2 h3]);
lgd.NumColumns = 1;
lgd.Location = "southeast";
