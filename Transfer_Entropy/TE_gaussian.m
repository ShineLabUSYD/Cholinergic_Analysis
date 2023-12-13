%% Transfer Entropy Calculation
%Set-up - need installation of JIDT toolbox https://github.com/jlizier/jidt
%Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
%javaaddpath('C:\\Users\\natas\\Documents\\PhD\\Code\\infodynamics-dist-1.6.1\\infodynamics.jar');
javaaddpath('C:\Users\natas\Documents\MATLAB\infodynamics-dist-dev\infodynamics.jar');
% Add utilities to the path
addpath('C:\Users\natas\Documents\MATLAB\infodynamics-dist-dev\demos\octave');
% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian');
% 0. Load/prepare the data:
load('C:\Users\natas\Documents\PhD\Turchi_nrg_data_and_code\data_order_combo1.mat') %order for trials, 1 & 2 correspond to inhibition for monkeyF
%load('C:\Users\natas\Documents\PhD\Turchi_nrg_data_and_code\ts_monkeyZb.mat');
%other monkey
load('C:\Users\natas\Documents\PhD\Turchi_nrg_data_and_code\ts_monkeyFb.mat');

% 2. Load variables

% define the data of interest
%data_order=1&&2; %load in the session for inhibited nbM
ts = ts_monkeyFb(:,:,data_order_combo1==1 | data_order_combo1==2);
%ts = ts_monkeyFb(:,:,data_order_combo1==3); sham injection
%new_data_order = data_order_combo1(30:end,1); order for monkeyZ
%ts = ts_monkeyZb(:,:,new_data_order==6); %L & R sham conditions


%predefine variable sizes
ks=zeros(size(ts,1),1);
pval=zeros(size(ts,1));
TE_results = zeros(size(ts,1));
acf_roi=zeros(size(ts,1),51);


calc.setProperty('NORMALISE', 'true'); % Normalise the individual variables
for d = 1:size(ts,1)
    [acf, lags] = autocorr(ts(d,:,1), 'NumLags', 50); %up to 50 lags
    acfDecayTime(d)=50;
    for t = 1 : 50
		if (acf(t) < exp(-1))
			acfDecayTime(d) = t;
			break;
		end
	end
    calc.setProperty('DYN_CORR_EXCL', string(acfDecayTime(d))); %dynamic measure to 'overcome' autocorrelation, for the target, we set the autocorrelation
    % 2. Set any properties to non-default values:
    calc.setProperty('AUTO_EMBED_METHOD', 'MAX_CORR_AIS_DEST_ONLY');
    calc.setProperty('AUTO_EMBED_K_SEARCH_MAX', '10'); %optimising up to 10 different params.
    calc.setProperty('AUTO_EMBED_TAU_SEARCH_MAX', '2');
    for s = 1:size(ts,1)
        
		% For each source-dest pair:
	    if (s == d)
			continue;
        end
       if d ==1 %this only needs to be done once (for first run, rather than looping)
            calc.setProperty('AUTO_EMBED_METHOD', 'MAX_CORR_AIS_DEST_ONLY');
            calc.setProperty('AUTO_EMBED_K_SEARCH_MAX', '10'); %optimising up to 10 different params.
            calc.setProperty('AUTO_EMBED_TAU_SEARCH_MAX', '2');
       else
           calc.setProperty('AUTO_EMBED_METHOD', 'NONE');    %turns autoembedding off after setting parameter      
           calc.setProperty('k_HISTORY', string(ks(d)));    %turns autoembedding off        
        end
		% 3. Initialise the calculator for (re-)use:
		calc.initialise();
        calc.startAddObservations(); %be ready for multiple observations
		% 4. Supply the sample data:
		for i=1:size(ts,3)
            source = octaveToJavaDoubleArray(ts(s,:,i));
		    destination = octaveToJavaDoubleArray(ts(d,:,i));
            calc.addObservations(source, destination);
        end
        calc.finaliseAddObservations();
		% 5. Compute the estimate:
		result = calc.computeAverageLocalOfObservations();
        ks(d)=str2num(calc.getProperty('k_HISTORY'));
		% 6. Compute the (statistical significance via) null distribution empirically (e.g. with 100 permutations):
		measDist = calc.computeSignificance(1000);
        pval(s,d)=measDist.pValue; %pval comparison from the above line of shuffling source
        TE_results(s,d)=result;
        acf_roi(d,:)=acf;
		%fprintf('TE_Gaussian (KSG)(col_%d -> col_%d) = %.4f nats\n', ...
            %s, d, result);
	end
end
%inhibition monkeyF
save('TE_gaussian_1000_inhibition_monkeyF.mat',"acf_roi","TE_results","pval","ks")
% sham monkeyF
%save('TE_gaussian_1000_sham_monkeyF.mat',"acf_roi","TE_results","pval","ks")
%inhibition monkeyZ
%save('TE_gaussian_1000_inhibition_monkeyZ.mat',"acf_roi","TE_results","pval","ks")
%sham monkeyZ
%save('TE_gaussian_1000_sham_monkeyZ.mat',"acf_roi","TE_results","pval","ks")