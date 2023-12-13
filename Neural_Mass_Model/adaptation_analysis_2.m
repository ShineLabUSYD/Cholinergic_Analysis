%% Neural Mass Model
% All corresponding code requires Brain Dynamics Toolbox to run the
% following script

%Running test case for parameter fitting
sigma = 0.06; %noise term
d = 0.25;
a=1;
sys = StochasticWilsonWTA_Tash_Sigmoid_adapt(sigma,d,a); %change c params. as adaptation
%gui open
gui = bdGUI(sys)
 %solver caller to generate timeseries wihtout using gui
sol=bdSolve(sys);
%extract time-series
i=[1,3]; %whichever population of interest
%we want E_L and E_R
ts=sol.y(i,:);
ts_orig = ts;
%need to downsample the timeseries
ts_downsample = downsample(ts_orig(1,:),1000); %downsample to 10ms bin
ts_downsample(2,:) = downsample(ts_orig(2,:),1000);
% check params. attractor
    %% Determine params:
autocorr_ts= autocorr(ts_downsample(2,:)); %check both E pops.
figure
plot(autocorr_ts') %find point of inflection in plot (0 point in y-axis)
TRprepost = 20;
cortSig= ts_downsample(1,:)';
for dt = -1*TRprepost+1:TRprepost
    MSD = mean((cortSig(dt+TRprepost:end-TRprepost,:) - cortSig(TRprepost:end-dt-TRprepost,:)).^2,2) ;
end
%     figure
%     histogram(MSD)

%% Sdaptation param. changing between 0.25-0.4
sigma = 0.05; %noise param
d = 0.25; % adaptation
a = 1; %bias the e-population
downsam = 1000;
nTR = 20;
[nrgSig_adapt_model,ts_adapt,sys,f] = adaptation_analysis(sigma,d,a,downsam,nTR);

%save model variables and associated params. 
savefilename = sprintf('%s%d%s%d%s','model_adapt_d_',d,'_a_',a,'.mat');
save([savefilename],'sys','ts_adapt','nrgSig_adapt_model');

figname = sprintf('%s%d%s%d%s','adapt_d_',d,'_sigma_',sigma,'.fig');
savefig(f,figname)

%create figure plot for the model

figure
set(gcf,'Color','w');
%load('model_adapt_d_2.000000e-01_a_1.mat')
avg_nrgSig_1=mean(nrgSig_adapt_model);
plot(avg_nrgSig_1,'Color',[0.4,0.7,0.5],'LineWidth',3) %green colour is adapt 0.25
%load('model_adapt_d_4.000000e-01_a_1.mat')
%avg_nrgSig_3=mean(nrgSig_adapt_model);
hold on
plot(avg_nrgSig_3,'Color',[1.0,0.7,0.0],'LineWidth',3) %orange is adapt 0.4


%% Change in excitability
sigma = 0.05; %noise param
d = 0.4; % adaptation
a = 1.3; %bias the e-population
downsam = 1000;
nTR = 20;
[nrgSig_adapt_model,ts_adapt,sys,f] = adaptation_analysis(sigma,d,a,downsam,nTR);

%save model variables and associated params. 
savefilename = sprintf('%s%d%s%d%s','model_adapt_d_',d,'_a_',a,'.mat');
save([savefilename],'sys','ts_adapt','nrgSig_adapt_model');

figname = sprintf('%s%d%s%d%s%d%s','adapt_d_',d,'_sigma_',sigma,'_excit_',a,'.fig');
savefig(f,figname)

%create figure plot for the model

figure
set(gcf,'Color','w');
load('model_adapt_d_4.000000e-01_a_1.300000e+00.mat')
avg_nrgSig_1=mean(nrgSig_adapt_model);
plot(avg_nrgSig_1,'Color',[0.4,0.7,0.5],'LineWidth',3) %green colour is decreased adaptation
load('model_adapt_d_4.000000e-01_a_1.100000e+00.mat')
avg_nrgSig_2=mean(nrgSig_adapt_model);
hold on
plot(avg_nrgSig_2,'Color',[0.8,0.7,0.2],'LineWidth',3) %inbetween is 1.
load('model_adapt_d_4.000000e-01_a_9.000000e-01.mat')
avg_nrgSig_3=mean(nrgSig_adapt_model);
hold on
plot(avg_nrgSig_3,'Color',[1.0,0.7,0.0],'LineWidth',3) %orange is 0.9 excit
% log these values
figure
set(gcf,'Color','w');
%load('model_adapt_d_4.000000e-01_a_1.300000e+00.mat')
log_nrgSig_1 = log(avg_nrgSig_1);
plot(log_nrgSig_1,'Color',[0.4,0.7,0.5],'LineWidth',3) %green colour is decreased adaptation
%load('model_adapt_d_4.000000e-01_a_1.100000e+00.mat')
log_nrgSig_2 = log(avg_nrgSig_2);
hold on
plot(log_nrgSig_2,'Color',[0.8,0.7,0.2],'LineWidth',3) %inbetween is 1.
%load('model_adapt_d_4.000000e-01_a_9.000000e-01.mat')
log_nrgSig_3 = log(avg_nrgSig_3);
hold on
plot(log_nrgSig_3,'Color',[1.0,0.7,0.0],'LineWidth',3) %orange is 0.9 excit




%% Change both adaptibility & excitation
sigma = 0.05; %noise param
d = 0.25; % adaptation
a = 1.3; %bias the e-population
downsam = 1000;
nTR = 20;
[nrgSig_adapt_model,ts_adapt,sys,f] = adaptation_analysis(sigma,d,a,downsam,nTR);

%save model variables and associated params. 
savefilename = sprintf('%s%d%s%d%s','model_adapt_d_',d,'_a_',a,'.mat');
save([savefilename],'sys','ts_adapt','nrgSig_adapt_model');

figname = sprintf('%s%d%s%d%s%d%s','adapt_d_',d,'_sigma_',sigma,'_excit_',a,'.fig');
savefig(f,figname)

%create figure plot for the model
figure
set(gcf,'Color','w');
avg_nrgSig_1=mean(nrgSig_adapt_model);
plot(avg_nrgSig_1,'Color',[0.4,0.7,0.5],'LineWidth',3)
hold on
plot(avg_nrgSig_2,'Color',[1.0,0.7,0.0],'LineWidth',3) %is 1

%% Attractor landscape
nTR = 20; % number displacement into time future
nMSD = 0.35; %msd range calculated across
rMSD = 0.01;

[nrgSig_combine,MSD_combine]=msd_calc_model(ts_downsample',nMSD,rMSD,nTR);
[nrgSig_combine,MSD_combine]=msd_calc_model(ts_orig',nMSD,rMSD,nTR);
[nrgSig_combine_EL,MSD_combine_EL]=msd_calc_model(ts_downsample(1,:)',nMSD,rMSD,nTR);
[nrgSig_combine_ER,MSD_combine_ER]=msd_calc_model(ts_downsample(2,:)',nMSD,rMSD,nTR);

[nrgSig_orig,MSD_orig]=msd_calc_model(ts_orig',nMSD,rMSD,nTR);

%plot
load('colormap.mat')
figure
set(gcf,'Color','w');
x = 1:size(nrgSig_combine,1); %both populations
y = 0:rMSD:nMSD;
[X,Y] = meshgrid(x,y);
mesh(X,Y,nrgSig_combine','FaceAlpha',0.5,'FaceColor','flat')
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
colormap(grad)
title('both populations adaptation 1')


figure
subplot(1,2,1)
set(gcf,'Color','w');
x = 1:size(nrgSig_combine_EL,1);
y = 0:rMSD:nMSD;
[X,Y] = meshgrid(x,y);
mesh(X,Y,nrgSig_combine_EL','FaceAlpha',0.5,'FaceColor','flat')
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
colormap(grad)

subplot(1,2,2)
set(gcf,'Color','w');
x = 1:size(nrgSig_combine_ER,1);
y = 0:rMSD:nMSD;
[X,Y] = meshgrid(x,y);
mesh(X,Y,nrgSig_combine_ER','FaceAlpha',0.5,'FaceColor','flat')
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
colormap(grad)
title('Adaptaion Er = 1.5 and El = 0.2')




