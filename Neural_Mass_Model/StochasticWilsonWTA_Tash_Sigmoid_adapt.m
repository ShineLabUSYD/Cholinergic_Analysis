% Model of Cholinergic system influence on 2 excitatory and one strong
% inhibitory node

% Example:
%sys = StochasticWilsonWTA_Tash_Sigmoid_bd(tau_H_L,tau_H_R,c);
%   gui = bdGUI(sys)

function sys = StochasticWilsonWTA_Tash_Sigmoid_adapt(sigma,d,a)
    % Handles to our SDE functions
    sys.sdeF   = @sdeF;                 % deterministic coefficients
    sys.sdeG   = @sdeG;                 % stochastic coefficients

    
    % ODE params 
    sys.pardef = [
        struct('name','a','value', a) %equivalent, as we don't want un
        struct('name','b','value', 1)
        struct('name','d','value', d)
        struct('name','cL','value', 1)
        struct('name','cR','value', 1)
        struct('name','tau','value', 1)
        struct('name','tau_H_L','value', 1000)
        struct('name','tau_H_R','value', 1000)
        struct('name','alpha','value', 1.5)
        struct('name','theta','value', 1)
        struct('name','delta','value', 6)
        struct('name','sigma','value', sigma) %std noise term
        ];
    
    % state variables
    sys.vardef = [
        struct('name', 'E_L', 'value', 0)
        struct('name', 'H_L', 'value', 0)
        struct('name', 'E_R', 'value', 0)
        struct('name', 'H_R', 'value', 0)
        ];

    % Default time span
    sys.tspan = [0 10000];

    % Specify SDE solvers and default options
    sys.sdesolver = {@sdeEM,@sdeSH};    % Relevant SDE solvers
    sys.sdeoption.InitialStep = 0.01;   % SDE solver step size (optional)
    sys.sdeoption.NoiseSources = 4;     % Number of driving Wiener processes, number of noise terms (4 populations of neurons)

    % Display panels
    sys.panels = struct();
    sys.panels.bdTimePortrait = [];
    
    % Default ODE options
%     sys.odeoption.RelTol = 1e-5;
%     sys.odeoption.InitialStep = 0.1;
    
        % Latex Panel
    sys.panels.bdLatexPanel.latex = {
        '$\textbf{Noise driven oscillator based on Wilson 2007}$'
        ''
        ''
        '$\tau_{L} \dot E_{L} = -E_{L} + f[a E_{L} - b E_{R} - c_{L} H_{L}]$'
        '$\tau_{H_L} \dot H_{L} = -H_{L} + d E_{L}$' %new params. that varies excitation with inhibition
        '$\tau_{R} \dot E_{R} = -E_{R} + f[a E_{R} - b E_{L} - c_{R} H_{R}]$'
        '$\tau_{H_R} \dot H_{R} = -H_{R} + d E_{R}$'
        };
    
    % Other Panels
    sys.panels.bdTimePortrait = [];
    sys.panels.bdPhasePortrait = [];
    sys.panels.bdSolverPanel = [];

end 

%needs to change to sdeF (deterministic coefficients)
function dY = sdeF(~,Y,a,b,d,cR,cL,tau,tau_H_L,tau_H_R,alpha,theta,delta,~)

    % extract state variables
    E_L = Y(1);
    H_L = Y(2);
    E_R = Y(3);
    H_R = Y(4);
    

    % system of equations
    dE_L = (-E_L + (1-E_L)*F(a*E_L - b*E_R - cL*H_L,alpha,theta,delta))*(1/tau); %cL adaptation for Left Pop.
    dH_L = (-H_L + d*E_L)*(1/tau_H_L);
    dE_R = (-E_R + (1-E_R)*F(a*E_R - b*E_L - cR*H_R,alpha,theta,delta))*(1/tau); %cR adaptation for Right Pop.
    dH_R = (-H_R + d*E_R)*(1/tau_H_R);
    
    dY = [dE_L; dH_L; dE_R; dH_R];

end 

% Sigmoidal firing-rate function
function y = F(x,alpha, theta, delta)
    y = 1./(1+exp(-alpha*(delta*x-theta)));
end

% include diffusion function sdeG --> adding noise solver (stochastic
% coefficient)
% The stochastic part of the equation.
%identity matrix


function G = sdeG(~,~,~,~,~,~,~,~,~,~,~,~,~,sigma)
noise_mat = [1 0 0 0;
             0 0 0 0;
             0 0 1 0;
             0 0 0 0]; %only impacting the noise for firing rates, not impacting adaptation
    %G = sigma .* eye(numel(Y)); %each population is getting same amount of noise (across the diagonal of identity matrix for all terms)
    G = sigma .* noise_mat;
end


