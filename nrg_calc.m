function nrgSig = nrg_calc(data,dsmtx,nMSD,nTR)

    % Inputs
    % data = data organised at TIME x NODES
    % dsmtx = onsets for different events - list of TRs associated with each event type
    % nMSD = maximum MSD range
    % nTR = number of TRs to track

    % Outputs
    % nrg = MSD x TR x Energy

    nReg = size(dsmtx,2);
    ds = 0:1:nMSD;
    nrgSig = nan(nTR,numel(ds));
    nrgBase = nan(nTR,numel(ds));

    for dt = 1:nTR

        %% MSD calculation
        MSD = mean((data(1+dt:end,:) - data(1:end-dt,:)).^2,2);
        MSDSig = MSD(dsmtx);
        
        %% Calculate probability distribution  and energy for each dt and each neuromod
        
        nrgSigdt = -1.*log(pdf(fitdist(MSDSig,'Kernel','BandWidth',4),ds));
    
        % Pool results across time
        nrgSig(dt,:) = nrgSigdt;
        
    end

end