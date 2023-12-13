function [nrgSig,MSD] = msd_calc_model(data,nMSD,rMSD,nTR)

    % Inputs
    % data = data organised at TIME x NODES
    % dsmtx = onsets for different events - list of TRs associated with each event type
    % nMSD = maximum MSD range
    % nTR = number of TRs to track
    %rMSD = range of msd
    % Outputs
    % nrg = MSD x TR x Energy

    %% Determine params:
%    for i=1:size(ts,2) ts = timeseries
%         autocorr_ts(i,:) = autocorr(doublets(:,i));
%     end
%     figure
%     plot(autocorr_ts') %find point of inflection in plot (0 point in y-axis)
%     % MSD calculation
%     TRprepost = 10;
%     ts_cort = ts(:,1:400);
%     cortSig = ts_cort;
%     
%     for dt = -1*TRprepost+1:TRprepost
%         MSD = mean((cortSig(dt+TRprepost:end-TRprepost,:) - cortSig(TRprepost:end-dt-TRprepost,:)).^2,2) ;
%     end
%     figure
%     histogram(MSD)

    %nReg = size(dsmtx,2);
    ds = 0:rMSD:nMSD;
    nrgSig = nan(nTR,numel(ds));
    %nrgBase = nan(nTR,numel(ds));
    %nT = size(data,1);

    for dt = 1:nTR
        %% MSD calculation
        MSD = mean((data(1+dt:end,:) - data(1:end-dt,:)).^2,2);
        %MSDSig = MSD(dsmtx);
        %% Calculate probability distribution  and energy for each dt
        nrgSigdt = -1.*log(pdf(fitdist(MSD,'Kernel','BandWidth',4),ds));
    
        % Pool results across time
        nrgSig(dt,:) = nrgSigdt;
        
    end

end