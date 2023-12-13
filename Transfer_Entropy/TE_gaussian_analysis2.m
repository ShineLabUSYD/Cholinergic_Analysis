%% TE Gaussian Analysis - Number 2 - 1000 permutation significance iterations

%Analysis of monkey F

inhibit_monkeyF = load('TE_gaussian_1000_inhibit_monkeyF.mat');
sham_monkeyF =load('TE_gaussian_1000_sham_monkeyF.mat');

%significance from 10 point permutation in JIDT toolbox
sig_pval_inhibit = double(inhibit_monkeyF.pval<0.05);
sig_inhibit_monkeyF = (inhibit_monkeyF.TE_results).*sig_pval_inhibit;
sig_pval_sham = double(sham_monkeyF.pval<0.05);
sig_sham_monkeyF = (sham_monkeyF.TE_results).*sig_pval_sham;

c = [0.4 0.7 0.5;0 0 0;0.7 0.7 0.7];%grey color 0.7,0.7,0.7
%create grad. colour

%standard
figure
set(gcf,'Color','w');
imagesc(sig_inhibit_monkeyF)
xlabel('ROIs')
ylabel('ROIs')
title('sig. inhibited TE monkey F')

figure
set(gcf,'Color','w');
imagesc(sig_sham_monkeyF)
xlabel('ROIs')
ylabel('ROIs')
title('sig. sham TE monkey F')

%flatten matrices for ease of plotting
flat_sig_inhibit_monkeyF = reshape(sig_inhibit_monkeyF,[266*266,1]);
flat_sig_sham_monkeyF = reshape(sig_sham_monkeyF,[266*266,1]);

diff_monkeyF = flat_sig_inhibit_monkeyF - flat_sig_sham_monkeyF;

%plot of sig. difference all ROIs inhibited - sham monkeyF
figure
set(gcf,'Color','w');
imagesc(reshape(diff_monkeyF,[266,266]))
title('Sig. Difference TE inhibited - sham monkeyF')

% direct vs. indirect results
load('ChAM_AL_id2.mat');

%directed inhibit
sig_direct_inhibit_monkeyF = sig_inhibit_monkeyF(ChAM_AL_id2(:,1)==1,ChAM_AL_id2(:,1)==1);
%directed sham
sig_direct_sham_monkeyF = sig_sham_monkeyF(ChAM_AL_id2(:,1)==1,ChAM_AL_id2(:,1)==1);
%indirect inhibit
sig_indirect_inhibit_monkeyF = sig_inhibit_monkeyF(ChAM_AL_id2(:,1)==0,ChAM_AL_id2(:,1)==0);
%indirect sham
sig_indirect_sham_monkeyF = sig_sham_monkeyF(ChAM_AL_id2(:,1)==0,ChAM_AL_id2(:,1)==0);
% plots
figure
set(gcf,'Color','w');
a = reshape(sig_direct_inhibit_monkeyF,[30*30,1]);
a(a==0)=nan;
var_plot = [1:size(reshape(sig_direct_inhibit_monkeyF,[30*30,1]),1)];
set(gcf,'color','w');
c= [0.4,0.7,0.5]; %green from other plots
scatter(var_plot,a,[],c,'filled')
hold on
c= [0.7,0.7,0.7]; %grey scale color
b= reshape(sig_direct_sham_monkeyF,[30*30,1]);
b(b==0)=nan; %stops plotting zeros
scatter(var_plot,b,[],c,'filled')
ylabel('TE')
xlabel('Each ROI-ROI edge')
hold on
c= [0.8,0.8,0.8]; %light grey
for ii=1:(30*30)
    plot([var_plot(1,ii),var_plot(1,ii)],[a(ii,1),b(ii,1)],'Color',c)
end
title('Sig. direct edges of TE inhibition (green) vs sham (grey)')


% indirect
figure
set(gcf,'Color','w');
a = reshape(sig_indirect_inhibit_monkeyF,[236*236,1]);
a(a==0)=nan;
var_plot = [1:size(reshape(sig_indirect_inhibit_monkeyF,[236*236,1]),1)];
set(gcf,'color','w');
c= [0.4,0.7,0.5]; %green from other plots
scatter(var_plot,a,[],c,'filled')
hold on
c= [0.7,0.7,0.7]; %grey scale color
b= reshape(sig_indirect_sham_monkeyF,[236*236,1]);
b(b==0)=nan; %stops plotting zeros
scatter(var_plot,b,[],c,'filled')
ylabel('TE')
xlabel('Each ROI-ROI edge')
hold on
c= [0.8,0.8,0.8]; %light grey
for ii=1:size(var_plot,2)
    plot([var_plot(1,ii),var_plot(1,ii)],[a(ii,1),b(ii,1)],'Color',c)
end
title('Sig. indirect edges of TE inhibition (green) vs sham (grey)')

%% Plot results
%inhibition direct vs indirect TE


a = reshape(sig_direct_inhibit_monkeyF,[30*30,1]);
a(a==0)=nan;
b= reshape(sig_indirect_inhibit_monkeyF,[236*236,1]);
b(b==0)=nan;

figure
set(gcf,'color','w');
c= [0.4,0.7,0.5]; %green from other plots
%swarmchart(ones(size(a,1),1),a,[],c,'filled')
scatter(ones(size(a,1),1),a,[],c)
hold on
c= [0.7,0.7,0.7]; %grey scale color
scatter((ones(size(b,1),1)+1),b,[],c)
%% Figure Final
figure
set(gcf,'color','w');
c= [0.4,0.7,0.5]; %green from other plots
x1 =ones(size(a,1),1);
scatter(x1,a,[],c,'XJitter','randn') %a = reshape(sig_direct_inhibit_monkeyF,[30*30,1]); a(a==0)=nan;
hold on
c= [0.7,0.7,0.7]; %grey scale color
x2 = ones(size(b,1),1)+1;
%make it all twos
scatter(x2,b,[],c,'XJitter','randn') %b is sig indirect inhibit
%xlim([0.95,1.15])
xlabel('Direct vs Non-direct')
ylabel('Transfer Entropy')
title('MonkeyF inhibition direct vs indirect')

%% Figure Final - with sham injection too 

a1 = reshape(sig_direct_sham_monkeyF,[30*30,1]);
a1(a1==0)=nan;
b2= reshape(sig_indirect_sham_monkeyF,[236*236,1]);
b2(b2==0)=nan;

figure
set(gcf,'color','w');
c= [0.4,0.7,0.5]; %green from other plots
x1 =ones(size(a1,1),1);
scatter(x1,a1,[],c,'XJitter','randn') %a = reshape(sig_direct_inhibit_monkeyF,[30*30,1]); a(a==0)=nan;
hold on
c= [0.7,0.7,0.7]; %grey scale color
x2 = ones(size(b2,1),1)+1;
%make it all twos
scatter(x2,b2,[],c,'XJitter','randn') %b is sig indirect inhibit
%xlim([0.95,1.15])
xlabel('Direct vs Non-direct')
ylabel('Transfer Entropy')
title('MonkeyF sham direct vs indirect')

%% Figure Direct Inhibit vs Sham; Indirect inhibit vs Sham
%% Figure Final
figure
set(gcf,'color','w');
c= [0.4,0.7,0.5]; %green from other plots
x1 =ones(size(a,1),1);
scatter(x1,a,[],c,'XJitter','randn') %a = reshape(sig_direct_inhibit_monkeyF,[30*30,1]); a(a==0)=nan;
hold on
%c= [0.7,0.7,0.7]; %grey scale color
x2 = ones(size(a1,1),1)+1;
%make it all twos
scatter(x2,a1,[],c,'XJitter','randn') %a1 is direct sham
%xlim([0.95,1.15])
xlabel('Inhibit vs Sham')
ylabel('Transfer Entropy')
title('MonkeyF inhibition direct inhibit vs sham')
% significant difference
%run stats - paired t-test on the already significant te's
flat_sig_inhibit_direct_monkeyF = a;
flat_sig_sham_direct_monkeyF = a1;
[hh_direct_monkeyF,pp_direct_monkeyF,~,stats_direct_monkeyF]=ttest2(flat_sig_inhibit_direct_monkeyF,flat_sig_sham_direct_monkeyF);

figure
set(gcf,'color','w');
%c= [0.4,0.7,0.5]; %green from other plots
c= [0.7,0.7,0.7]; %grey scale color
x1 =ones(size(b,1),1); %looking indirect inhibit
scatter(x1,b,[],c,'XJitter','randn') 
hold on
c= [0.7,0.7,0.7]; %grey scale color
x2 = ones(size(b2,1),1)+1;
%make it all twos
scatter(x2,b2,[],c,'XJitter','randn') %b2 is sig indirect sham
%xlim([0.95,1.15])
xlabel('Direct vs Non-direct')
ylabel('Transfer Entropy')
title('MonkeyF indirect inhibit vs sham')
flat_sig_inhibit_indirect_monkeyF = b;
flat_sig_sham_indirect_monkeyF = b2;
[hh_indirect_monkeyF,pp_indirect_monkeyF,~,stats_indirect_monkeyF]=ttest2(flat_sig_inhibit_indirect_monkeyF,flat_sig_sham_indirect_monkeyF);



%scatter plots of comparison
figure
set(gcf,'Color','w');
var_plot = [1:size(flat_sig_inhibit_monkeyF,1)];
set(gcf,'color','w');
c= [0.4,0.7,0.5]; %green from other plots
scatter(var_plot,flat_sig_inhibit_monkeyF,[],c,'filled')
hold on
c= [0.7,0.7,0.7]; %grey scale color
scatter(var_plot,flat_sig_sham_monkeyF,[],c,'filled')
ylabel('TE')
xlabel('ROIs')
hold on
c= [0.8,0.8,0.8]; %light grey
for ii=1:size(flat_sig_inhibit_monkeyF,1)
    plot([var_plot(1,ii),var_plot(1,ii)],[flat_sig_inhibit_monkeyF(ii,1),flat_sig_sham_monkeyF(ii,1)],'Color',c)
end





%Analysis of monkey Z

inhibit_monkeyZ = load('TE_gaussian_1000_inhibit_monkeyZ.mat');
sham_monkeyZ =load('TE_gaussian_1000_sham_monkeyZ.mat');

%significance from 10 point permutation in JIDT toolbox
sig_pval_inhibit = double(inhibit_monkeyZ.pval<0.05);
sig_inhibit_monkeyZ = (inhibit_monkeyZ.TE_results).*sig_pval_inhibit;
sig_pval_sham = double(sham_monkeyZ.pval<0.05);
sig_sham_monkeyZ = (sham_monkeyZ.TE_results).*sig_pval_sham;

c = [0.4 0.7 0.5;0 0 0;0.7 0.7 0.7];%grey color 0.7,0.7,0.7
%create grad. colour

%standard
figure
set(gcf,'Color','w');
imagesc(sig_inhibit_monkeyZ)
xlabel('ROIs')
ylabel('ROIs')
title('sig. inhibited TE monkey Z')

figure
set(gcf,'Color','w');
imagesc(sig_sham_monkeyZ)
xlabel('ROIs')
ylabel('ROIs')
title('sig. sham TE monkey Z')

%flatten matrices for ease of plotting
flat_sig_inhibit_monkeyZ = reshape(sig_inhibit_monkeyZ,[266*266,1]);
flat_sig_sham_monkeyZ = reshape(sig_sham_monkeyZ,[266*266,1]);

diff_monkeyZ = flat_sig_inhibit_monkeyZ - flat_sig_sham_monkeyZ;

%plot of sig. difference all ROIs inhibited - sham monkeyF
figure
set(gcf,'Color','w');
imagesc(reshape(diff_monkeyZ,[266,266]))
title('Sig. Difference TE inhibited - sham monkeyZ')

% direct vs. indirect results
load('ChAM_AL_id2.mat');

%directed inhibit
sig_direct_inhibit_monkeyZ = sig_inhibit_monkeyF(ChAM_AL_id2(:,2)==1,ChAM_AL_id2(:,2)==1);
%directed sham
sig_direct_sham_monkeyZ = sig_sham_monkeyF(ChAM_AL_id2(:,2)==1,ChAM_AL_id2(:,2)==1);
%indirect inhibit
sig_indirect_inhibit_monkeyZ = sig_inhibit_monkeyF(ChAM_AL_id2(:,2)==0,ChAM_AL_id2(:,2)==0);
%indirect sham
sig_indirect_sham_monkeyZ = sig_sham_monkeyF(ChAM_AL_id2(:,2)==0,ChAM_AL_id2(:,2)==0);


%inhibition direct vs indirect TE


c = reshape(sig_direct_inhibit_monkeyZ,[20*20,1]);
c(c==0)=nan;
d= reshape(sig_indirect_inhibit_monkeyZ,[246*246,1]);
d(d==0)=nan;


figure
set(gcf,'color','w');
color= [0.4,0.7,0.5]; %green from other plots
x1 =ones(size(c,1),1);
scatter(x1,c,[],color,'XJitter','randn')
hold on
color= [0.7,0.7,0.7]; %grey scale color
x2 = ones(size(d,1),1)+1;
%make it all twos
scatter(x2,d,[],color,'XJitter','randn')
%xlim([0.95,1.15])
xlabel('Direct vs Non-direct')
ylabel('Transfer Entropy')
title('MonkeyZ inhibition direct vs indirect')

%% Figure Final

c1 = reshape(sig_direct_sham_monkeyZ,[20*20,1]);
c1(c1==0)=nan;
d1= reshape(sig_indirect_sham_monkeyZ,[246*246,1]);
d1(d1==0)=nan;


figure
set(gcf,'color','w');
color= [0.7,0.7,0.7]; %grey scale color
%color= [0.4,0.7,0.5]; %green from other plots
x1 =ones(size(d,1),1); %indirect inhibit
scatter(x1,d,[],color,'XJitter','randn') %indirect
hold on
%color= [0.7,0.7,0.7]; %grey scale color
x2 = ones(size(d1,1),1)+1;
%make it all twos
scatter(x2,d1,[],color,'XJitter','randn') % indirect sham
%xlim([0.95,1.15])
xlabel('Inhibit vs Sham')
ylabel('Transfer Entropy')
title('MonkeyZ inhibition indirect inhibit vs sham')
% significant difference
%run stats - paired t-test on the already significant te's
flat_sig_inhibit_indirect_monkeyZ = d;
flat_sig_sham_indirect_monkeyZ = d1;
[hh_indirect_monkeyZ,pp_indirect_monkeyZ,~,stats_indirect_monkeyZ]=ttest2(flat_sig_inhibit_indirect_monkeyZ,flat_sig_sham_indirect_monkeyZ);

figure
set(gcf,'color','w');
%c= [0.4,0.7,0.5]; %green from other plots
c= [0.7,0.7,0.7]; %grey scale color
x1 =ones(size(b,1),1); %looking indirect inhibit
scatter(x1,b,[],c,'XJitter','randn') 
hold on
c= [0.7,0.7,0.7]; %grey scale color
x2 = ones(size(b2,1),1)+1;
%make it all twos
scatter(x2,b2,[],c,'XJitter','randn') %b2 is sig indirect sham
%xlim([0.95,1.15])
xlabel('Direct vs Non-direct')
ylabel('Transfer Entropy')
title('MonkeyF indirect inhibit vs sham')
flat_sig_inhibit_indirect_monkeyF = b;
flat_sig_sham_indirect_monkeyF = b2;
[hh_indirect_monkeyF,pp_indirect_monkeyF,~,stats_indirect_monkeyF]=ttest2(flat_sig_inhibit_indirect_monkeyF,flat_sig_sham_indirect_monkeyF);















%% Combine analysis
 combine_direct_inhibit = vertcat(a,c);
 combine_indirect_inhibit = vertcat(b,d);

 
combine_direct_sham = vertcat((reshape(sig_direct_sham_monkeyF,[30*30,1])),reshape(sig_direct_sham_monkeyZ,[20*20,1]));
combine_indirect_sham = vertcat((reshape(sig_indirect_sham_monkeyF,[236*236,1])),reshape(sig_indirect_sham_monkeyZ,[246*246,1]));

%run stats - paired t-test on the already significant te's
[hh_combine_direct_indirect_inhibit,pp_combine_direct_indirect_inhibit,~,stats_combine_inhibit]=ttest2(combine_direct_inhibit,combine_indirect_inhibit);

 %% Final Figure
figure
set(gcf,'color','w');
color= [0.4,0.7,0.5]; %green from other plots
x1 =ones(size(combine_direct_inhibit,1),1);
scatter(x1,combine_direct_inhibit,[],color,'XJitter','randn')
hold on
color= [0.7,0.7,0.7]; %grey scale color
x2 = ones(size(combine_indirect_inhibit,1),1)+1;
%make it all twos
scatter(x2,combine_indirect_inhibit,[],color,'XJitter','randn')
ylim([-0.01,0.8])
%xlim([0.95,1.15])
xlabel('Direct vs Non-direct')
ylabel('Transfer Entropy')
title('Inhibition direct vs indirect both monkeys')


figure
set(gcf,'color','w');
color= [0.4,0.7,0.5]; %green from other plots
x1 =ones(size(combine_direct_inhibit,1),1);
scatter(x1,combine_direct_inhibit,[],color,'XJitter','randn')
hold on
color= [0.4,0.7,0.5]; %grey scale color
x2 =ones(size(combine_direct_sham,1),1)+1;
scatter(x2,combine_direct_sham,[],color,'XJitter','randn')
hold on
color= [0.7,0.7,0.7];
x3 =ones(size(combine_indirect_inhibit,1),1)+2;
scatter(x3,combine_indirect_inhibit,[],color,'XJitter','randn')
hold on
color= [0.7,0.7,0.7]; %grey scale color
x4 =ones(size(combine_indirect_sham,1),1)+3;
scatter(x4,combine_indirect_sham,[],color,'XJitter','randn')
ylim([-0.02,1.2])
%xlim([0.95,1.15])
%xlabel('Direct vs Non-direct')
ylabel('Transfer Entropy')
title('Inhibition direct vs indirect both monkeys')

%stats

[hh_combine_direct,pp_combine_direct,~,stats_combine_direct]=ttest2(combine_direct_inhibit,combine_direct_sham);

[hh_combine_indirect,pp_combine_indirect,~,stats_combine_indirect]=ttest2(combine_indirect_inhibit,combine_indirect_sham);

[hh_combine_sham,pp_combine_isham,~,stats_combine_sham]=ttest2(combine_direct_sham,combine_indirect_sham);

% sham injection TE plot
figure
set(gcf,'color','w');
color= [0.4,0.7,0.5]; %green from other plots
x1 =ones(size(combine_direct_sham,1),1);
scatter(x1,combine_direct_sham,[],color,'XJitter','randn')
hold on
color= [0.7,0.7,0.7]; %grey scale color
x2 = ones(size(combine_indirect_sham,1),1)+1;
%make it all twos
scatter(x2,combine_indirect_sham,[],color,'XJitter','randn')
ylim([-0.01,1.2])
%xlim([0.95,1.15])
xlabel('Direct vs Non-direct')
ylabel('Transfer Entropy')
title('Sham direct vs indirect both monkeys')

[hh_combine_direct_indirect_sham,pp_combine_direct_indirect_sham,~,stats_combine_sham]=ttest2(combine_direct_sham,combine_indirect_sham);

