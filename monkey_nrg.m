%% Attractor Landscape Analysis 

%Load data
load('ts_monkeyFb');
load('ts_monkeyZb');
% Data order file
load('data_order_combo');
% 1=left
% 2 = right
% 3 = sham

% Calculate Functional Connectivity - correlation matrices
for ii = 1:29
ts_corr_monkeyF(:,:,ii) = corr(ts_monkeyFb(:,:,ii)');
end

for ii = 1:45
ts_corr_monkeyZ(:,:,ii) = corr(ts_monkeyZb(:,:,ii)');
end

ts_corr_combo = ts_corr_monkeyF;
ts_corr_combo(:,:,30:74) = ts_corr_monkeyZ;


[hh_corr,pp_corr] = ttest2(ts_corr_combo(:,:,data_order_combo~=3),ts_corr_combo(:,:,data_order_combo==3),'dim',3);

diff_ts_corr = mean(ts_corr_combo(:,:,data_order_combo~=3),3) - mean(ts_corr_combo(:,:,data_order_combo==3),3);
figure
set(gcf,'color','w');
hex = ['#68ad84';'#9fb645';'#ffa600'];

imagesc(diff_ts_corr.*hh_corr)
xlabel('ROI')
ylabel('ROI')

%figure for functional connectivity
limits = [-1 1];
clims= limits;
figure
set(gcf,'color','w');
subplot(2,4,1)
imagesc(mean(ts_corr_monkeyF(:,:,data_order_combo(1:29)~=3),3),clims)
title('all stim.')
caxis(limits)
colormap(grad)
subplot(2,4,2)
imagesc((mean(ts_corr_monkeyF(:,:,data_order_combo(1:29)~=3),3)-(mean(ts_corr_monkeyF(:,:,data_order_combo(1:29)==3),3))),clims)
title('all stim. - sham')
caxis(limits)
colormap(grad)
subplot(2,4,3)
imagesc(mean(ts_corr_monkeyF(:,:,data_order_combo(1:29)==1),3),clims)
title('stim. left')
caxis(limits)
colormap(grad)
subplot(2,4,4)
imagesc(mean(ts_corr_monkeyZ(:,:,data_order_combo(30:74)==2),3),clims)
title('stim. right')
caxis(limits)
colormap(grad)

subplot(2,4,5)
imagesc(mean(ts_corr_monkeyZ(:,:,data_order_combo(30:74)~=3),3),clims)
title('all stim.')
caxis(limits)
colormap(grad)
subplot(2,4,6)
imagesc((mean(ts_corr_monkeyZ(:,:,data_order_combo(30:74)~=3),3)-(mean(ts_corr_monkeyZ(:,:,data_order_combo(1:29)==3),3))),clims)
title('all stim. - sham')
caxis(limits)
colormap(grad)
subplot(2,4,7)
imagesc(mean(ts_corr_monkeyZ(:,:,data_order_combo(30:74)==1),3),clims)
title('stim. left')
caxis(limits)
colormap(grad)
subplot(2,4,8)
imagesc(mean(ts_corr_monkeyZ(:,:,data_order_combo(30:74)==2),3),clims)
title('stim. right')
caxis(limits)
colormap(grad)
%% Attractor landscape
%determine params
for i=1:266
    for kk=1:29
    autocorr_ts_Fb(i,:,kk) = autocorr(ts_monkeyFb(i,:,kk));
    end
end
%nTR - should be 8, this is the inflection point of autocorr


%number of tr into future calc msd 
%ndt - number displacement into time future
ndt = 25;
ds = 0:1:50; % the msd range calculated across

rand_time = 1:20:650;
nTR = 8; % number displacement into time future
nMSD = 5; %msd range calculated across
% Run Attractor Landscape Over respective monkey time-series
for ii = 1:29
    nrg_monkeyFc(:,:,ii) = nrg_calc(ts_monkeyFb(:,:,ii)',rand_time',nMSD,nTR); %output 8x6x29 (sessions)
end

for ii = 1:45
    nrg_monkeyZc(:,:,ii) = nrg_calc(ts_monkeyZb(:,:,ii)',rand_time',nMSD,nTR);
end


% nrg landscape per session
%combine attractor landscape for both monkeys
nrg_combo = nrg_monkeyFc;
nrg_combo(:,:,30:74) = nrg_monkeyZc;
%average attractor landscape based on sham vs stim
mean_nrg_sham = mean(nrg_combo(:,:,(data_order_combo==3)),3);
mean_nrg_stim = mean(nrg_combo(:,:,(data_order_combo~=3)),3);

% MSD calc. for each ROI
nrg_monkeyFc_roi = zeros(8,6,29,266);
nrg_monkeyZc_roi = zeros(266,8,6,45);

for rr = 1:266
    for ii = 1:29
    nrg_monkeyFc_roi(:,:,ii,rr) = nrg_calc(squeeze(ts_monkeyFb(rr,:,ii))',rand_time',nMSD,nTR); 
    end
end
mean_nrg_monkeyFC_roi_stim = squeeze(mean(nrg_monkeyFc_roi(:,:,data_order_combo(1:29,:)~=3,:),3));
mean_nrg_monkeyFC_roi_sham = squeeze(mean(nrg_monkeyFc_roi(:,:,data_order_combo(1:29,:)==3,:),3));


for rr=1:266
     for ii = 1:45
    nrg_monkeyZc_roi(:,:,ii,rr) = nrg_calc(squeeze(ts_monkeyZb(rr,:,ii))',rand_time',nMSD,nTR);
    end
end

nrg_combo_roi = nrg_monkeyFc_roi;
nrg_combo_roi(:,:,30:74,:) = nrg_monkeyZc_roi;

%stat test on ROI MSD
[hh_nrg_roi,pp_nrg_roi] = ttest2(nrg_combo_roi(:,:,data_order_combo~=3,:),nrg_combo_roi(:,:,data_order_combo==3,:),'dim',3);

%figures of MSD for each ROI
x = 1:nTR;
y = 0:1:nMSD;
[X,Y] = meshgrid(x,y);

for ii=1:266
figure(ii)
set(gcf,'Color','w');
subplot(1,3,1)
mesh(X,Y,squeeze(mean(nrg_combo_roi(:,:,(data_order_combo~=3),ii),3))','EdgeColor',[255,166,0]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
% xlim([1 xmax])
% zlim([2 20])
% ylim([1 50])
view(-15,30)   % XZ
title('Average stim. both monkeys')
subplot(1,3,2)
mesh(X,Y,squeeze(mean(nrg_combo_roi(:,:,(data_order_combo==3),ii),3))','EdgeColor', [159,182,69]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
% xlim([1 xmax])
% zlim([2 20])
% ylim([1 50])
view(-15,30)   % XZ
title('Average sham. both monkeys')
subplot(1,3,3)
mesh(X,Y,squeeze(mean(nrg_combo_roi(:,:,(data_order_combo~=3),ii),3)-mean(nrg_combo_roi(:,:,(data_order_combo==3),ii),3))','EdgeColor',  [105,173,132]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
% xlim([1 xmax])
% zlim([2 20])
% ylim([1 50])
view(-15,30)   % XZ
title('Average stim - sham both monkeys')
ax = gcf; exportgraphics(ax, 'avg_msd_roi.pdf','Append',true)
end
%% Rugosity calculation
% calculate mean MSD for each monkey fo
nrg_stim_Fc = mean(nrg_monkeyFc(:,:,data_order_combo(1:29,:)~=3),3);
nrg_stim_Zc = mean(nrg_monkeyZc(:,:,data_order_combo(30:74,:)~=3),3);
nrg_sham_Fc = mean(nrg_monkeyFc(:,:,data_order_combo(1:29,:)==3),3);
nrg_sham_Zc = mean(nrg_monkeyZc(:,:,data_order_combo(30:74,:)==3),3);

avg_nrg_stim_both = mean(nrg_combo(:,:,data_order_combo~=3),3);
avg_nrg_sham_both = mean(nrg_combo(:,:,data_order_combo==3),3);
%difference in the gradient calculation plus figures saved
[px1,py1,px2,py2,D1,D2,] = diff_attract(avg_nrg_stim_both,avg_nrg_sham_both,'stimulation','sham',5,1,8,0,1);

max_diff_grad_x = max(diff_grad_px,[],'all')
max_diff_grad_y = max(diff_grad_py,[],'all')


x = 1:nTR;
y = 0:1:nMSD;
[X,Y] = meshgrid(x,y);

figure
set(gcf,'Color','w'); %set background to white
subplot(2,2,1)
mesh(X,Y,nrg_stim_Fc','EdgeColor', [105,173,132]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
view(-15,30)   % XZ
title('Avg. stim. Monkey F')
subplot(2,2,2)
mesh(X,Y,nrg_stim_Zc','EdgeColor', [105,173,132]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
view(-15,30)   % XZ
title('Avg. stim. Monkey Z')
subplot(2,2,3)
mesh(X,Y,nrg_sham_Fc','EdgeColor', [255,166,0]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
view(-15,30)   % XZ
title('Avg. sham. Monkey F')
subplot(2,2,4)
mesh(X,Y,nrg_sham_Zc','EdgeColor', [255,166,0]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
view(-15,30)   % XZ
title('Avg. sham. Monkey Z')



%% stats

% stim vs. sham
[hh_corr,pp_corr] = ttest2(ts_corr_combo(:,:,data_order_combo~=3),ts_corr_combo(:,:,data_order_combo==3),'dim',3);
[hh_nrg,pp_nrg] = ttest2(nrg_combo(:,:,data_order_combo~=3),nrg_combo(:,:,data_order_combo==3),'dim',3);



%% Hemisphere Differences
% 1=left
% 2 = right
% 3 = sham

%time-series separated by Left vs Right
ts_monkeyF_L = ts_monkeyFb(:,:,data_order_combo(1:29,:)==1);
ts_monkeyF_R = ts_monkeyFb(:,:,data_order_combo(1:29,:)==2);
ts_monkeyZ_L = ts_monkeyZb(:,:,data_order_combo(29:74,:)==1);
ts_monkeyZ_R = ts_monkeyZb(:,:,data_order_combo(29:74,:)==2);
% correlation matrices
for ii = 1:size(ts_monkeyF_L,3)
ts_corr_monkeyF_L(:,:,ii) = corr(ts_monkeyF_L(:,:,ii)');
end

for ii = 1:size(ts_monkeyF_R,3)
ts_corr_monkeyF_R(:,:,ii) = corr(ts_monkeyF_R(:,:,ii)');
end

for ii = 1:size(ts_monkeyZ_L,3)
ts_corr_monkeyZ_L(:,:,ii) = corr(ts_monkeyZ_L(:,:,ii)');
end

for ii = 1:size(ts_monkeyZ_R,3)
ts_corr_monkeyZ_R(:,:,ii) = corr(ts_monkeyZ_R(:,:,ii)');
end

%attractor landscape
%determine params
for i=1:266
    for kk=1:size(ts_monkeyF_L,3) %for left side
    autocorr_ts_Fb(i,:,kk) = autocorr(ts_monkeyF_L(i,:,kk));
    end
end
for i=1:266
    for kk=1:size(ts_monkeyZ_R,3) %for left side
    autocorr_ts_Zb(i,:,kk) = autocorr(ts_monkeyZ_R(i,:,kk));
    end
end
%figure
%plot(autocorr_ts_Zb(:,:,1)') %go through each timing session - determine
%inflection point
%nTR= 10 for Left/Right Side both monkey F/Z


%number of tr into future calc msd 
%ndt - number displacement into time future
ndt = 25;
ds = 0:1:50; % the msd range calculated across

rand_time = 1:20:650;
nTR = 8; % number displacement into time future
nMSD = 5; %msd range calculated across

for ii = 1:size(ts_monkeyF_L,3)
    nrg_monkeyF_L(:,:,ii) = nrg_calc(ts_monkeyF_L(:,:,ii)',rand_time',nMSD,nTR); %output 8x6x29 (sessions)
end

for ii = 1:size(ts_monkeyF_R,3)
    nrg_monkeyF_R(:,:,ii) = nrg_calc(ts_monkeyF_R(:,:,ii)',rand_time',nMSD,nTR); %output 8x6x29 (sessions)
end

for ii = 1:size(ts_monkeyZ_L,3)
    nrg_monkeyZ_L(:,:,ii) = nrg_calc(ts_monkeyZ_L(:,:,ii)',rand_time',nMSD,nTR); %output 8x6x29 (sessions)
end

for ii = 1:size(ts_monkeyZ_R,3)
    nrg_monkeyZ_R(:,:,ii) = nrg_calc(ts_monkeyZ_R(:,:,ii)',rand_time',nMSD,nTR); %output 8x6x29 (sessions)
end

% attractor landscape combine left & right for both monkeys
%combine attractor landscape for both monkeys
nrg_combo_L(:,:,1:11) = nrg_monkeyF_L;
nrg_combo_L(:,:,12:26) = nrg_monkeyZ_L;
nrg_combo_R(:,:,1:10) = nrg_monkeyF_R;
nrg_combo_R(:,:,11:28) = nrg_monkeyZ_R;
%average attractor landscape based on sham vs stim
mean_nrg_L = mean(nrg_combo_L,3);
mean_nrg_R = mean(nrg_combo_R,3);




% plots
%figures of MSD for each ROI
x = 1:nTR;
y = 0:1:nMSD;
[X,Y] = meshgrid(x,y);

figure
set(gcf,'Color','w');
subplot(1,2,1)
mesh(X,Y,mean_nrg_L','EdgeColor',[255,166,0]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
% xlim([1 xmax])
% zlim([2 20])
% ylim([1 50])
view(-15,30)   % XZ
title('Average stim. to Left')
subplot(1,2,2)
mesh(X,Y,mean_nrg_R','EdgeColor', [159,182,69]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
% xlim([1 xmax])
% zlim([2 20])
% ylim([1 50])
view(-15,30)   % XZ
title('Average sham. to Right')


subplot(1,3,3)
mesh(X,Y,squeeze(mean(nrg_combo_roi(:,:,(data_order_combo~=3),ii),3)-mean(nrg_combo_roi(:,:,(data_order_combo==3),ii),3))','EdgeColor',  [105,173,132]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')


%% Attractor Landscape across all pulled sessions


%number of tr into future calc msd 
%ndt - number displacement into time future
ndt = 25;
ds = 0:1:50; % the msd range calculated across

rand_time = 1:20:650;
nTR = 8; % number displacement into time future
nMSD = 5; %msd range calculated across

%concatenate the time-series into one
concate_tsFb = ts_monkeyFb(:,:,1:21,1);
concate_tsFb = reshape(concate_tsFb,[266,713*21]);
concate_tsZb = ts_monkeyZb(:,:,data_order_combo(30:62,1));
concate_tsZb = reshape(concate_tsZb,[266,713*33]);
% Run Attractor Landscape Over respective monkey time-series

nrg_monkeyFc = nrg_calc(concate_tsFb',rand_time',nMSD,nTR); %output 8x6xtime(sessions)

nrg_monkeyZc = nrg_calc(concate_tsZb',rand_time',nMSD,nTR);




%averages
nrg_stim_Fc = mean(nrg_monkeyFc(:,:,:,3);
nrg_stim_Zc = mean(nrg_monkeyZc(:,:,:,3);
nrg_sham_Fc = mean(nrg_monkeyFc(:,:,:,3);
nrg_sham_Zc = mean(nrg_monkeyZc(:,:,:,3);

%% Figure

x = 1:nTR;
y = 0:1:nMSD;
[X,Y] = meshgrid(x,y);

figure
subplot(1,2,1)
mesh(X,Y,(mean(nrg_combo(:,:,(data_order_combo~=3)),3))','EdgeColor', [105,173,132]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
% xlim([1 xmax])
% zlim([2 20])
% ylim([1 50])
view(-15,30)   % XZ
title('Average stim. both monkeys')
subplot(1,2,2)
mesh(X,Y,(mean(nrg_combo(:,:,(data_order_combo==3)),3))','EdgeColor', [105,173,132]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
% xlim([1 xmax])
% zlim([2 20])
% ylim([1 50])
view(-15,30)   % XZ
title('Average sham. both monkeys')

%average the landscape and take difference
diff_nrg = mean(nrg_combo(:,:,(data_order_combo~=3)),3) - mean(nrg_combo(:,:,(data_order_combo==3)),3);
diff_nrg_flip = diff_nrg';
diff_nrg_reorder = zeros(6,8);
diff_nrg_reorder(1,:) = diff_nrg_flip(6,:);
diff_nrg_reorder(2,:) = diff_nrg_flip(5,:);
diff_nrg_reorder(3,:) = diff_nrg_flip(4,:);
diff_nrg_reorder(4,:) = diff_nrg_flip(3,:);
diff_nrg_reorder(5,:) = diff_nrg_flip(2,:);
diff_nrg_reorder(6,:) = diff_nrg_flip(1,:);


figure

mesh(X,Y,diff_nrg','EdgeColor', [105,173,132]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
% xlim([1 xmax])
% zlim([2 20])
% ylim([1 50])
view(-15,30)   % XZ
f.CurrentAxes.YDir = 'Reverse'; %reverse the y-axis
title('Average difference between stim. vs sham')
hold on
mesh(X,Y,flat_zero','facealpha',0.2,'edgealpha',0.2) %adds lines at 0 with transparency
set(gcf,'Color','w'); %set background to white

% figure left vs. right
[hh_nrg_L,pp_nrg_L] = ttest2(nrg_combo(:,:,data_order_combo==1),nrg_combo(:,:,data_order_combo==3),'dim',3);
diff_nrg_L =  mean(nrg_combo(:,:,(data_order_combo==1)),3) - mean(nrg_combo(:,:,(data_order_combo==3)),3);
[hh_nrg_R,pp_nrg_R] = ttest2(nrg_combo(:,:,data_order_combo==2),nrg_combo(:,:,data_order_combo==3),'dim',3);
diff_nrg_R =  mean(nrg_combo(:,:,(data_order_combo==2)),3) - mean(nrg_combo(:,:,(data_order_combo==3)),3);


figure
subplot(1,2,1)
set(gcf,'Color','w');
mesh(X,Y,diff_nrg_L','EdgeColor', [105,173,132]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
% xlim([1 xmax])
% zlim([2 20])
% ylim([1 50])
view(-15,30)   % XZ
title('Difference between stim. L vs sham')
subplot(1,2,2)
mesh(X,Y,diff_nrg_R','EdgeColor', [105,173,132]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
% xlim([1 xmax])
% zlim([2 20])
% ylim([1 50])
view(-15,30)   % XZ
title('Difference between stim. R vs sham')


%% Final Figure

figure
set(gcf,'color','w')
mesh(X,Y,(diff_nrg'+0.035),'EdgeColor', [105,173,132]./255)
xlabel('TR')
ylabel('MSD')
zlabel('MSD  energy')
zticklabels([-0.035,-0.025,-0.015,0,0.005,0.01])
view(-15,30)   % XZ
hold all
c= diff_nrg.*hh_nrg;
imagesc(c')
axis('xy')
xlim([0 8])
ylim([0 6])
%xlabel('TR')
%ylabel('MSD')
 %flips the axis so they are in alignment
%color selection - lower green , mid green [0.4,0.7,0.5],
%upper orange [1,0.7,0.0]
