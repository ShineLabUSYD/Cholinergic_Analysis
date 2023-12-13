function [nrgSig_adapt_model,ts_adapt,sys,f] = adaptation_analysis(sigma,d,a,downsam,nTR)
        %sigma = noise param
        % d = adaptation param in equation of model
        % downsam = play around with downsampling data
    sys = StochasticWilsonWTA_Tash_Sigmoid_adapt(sigma,d,a); %change d params. as adaptation
    sol=bdSolve(sys);
    %extract time-series
    i=[1,3]; %whichever population of interest
    %we want E_L and E_R
    ts=sol.y(i,:);
    ts_orig = ts;
    %need to downsample the timeseries
    ts_downsample = downsample(ts_orig(1,:),downsam); %downsample to 10ms bin
    ts_downsample(2,:) = downsample(ts_orig(2,:),downsam);
    ts_adapt =ts_orig;
    %figure timeseries
   f(1)= figure;
        plot(ts_downsample')



    %% Attractor landscape
    %nTR = 20; % number displacement into time future
    nMSD = 0.4; %msd range calculated across
    rMSD = 0.01;
    
    [nrgSig_adapt_model,MSD_adapt_model]=msd_calc_model(ts_downsample',nMSD,rMSD,nTR);
    

    load('colormap.mat')
    f(2)=figure;
        set(gcf,'Color','w');
        x = 1:size(nrgSig_adapt_model,1); %both populations
        y = 0:rMSD:nMSD;
        [X,Y] = meshgrid(x,y);
        mesh(X,Y,nrgSig_adapt_model','FaceAlpha',0.5,'FaceColor','flat')
        xlabel('TR')
        ylabel('MSD')
        zlabel('MSD  energy')
        colormap(grad)
        titlename = sprintf('%s%d%s%d','Adaptation Param =',d,' Sigma =',sigma);
        title(titlename)

end
