
clc
clear all
close all
%rng(921111)
for n1 = 1:3
    for n2 = 1:1
    % Population of N cells migrating on a square domain Omega
    dt    = 1;         % one time unit is one hour; one time step is a second
    freq  = 180;              % frequency of observations in seconds
    h     = 24;           % total number of hours the simulation will be run
    Kobs  = h*3600/freq;      % number of of observations
    K     = round(h*3600);    % total number of discrete time steps [seconds]
    scale = 50;             % one length unit is one average cell radius. size
    % domain
    N     = 32*2^n1;             % size of initial population
    % particle parameters
    Vmin  = 0.000075;           % depth of potential
    alpha = 4;                % steepness of potential
    sigma = 0.01;
    % population growth parameters
    lambda0 = 0.00/3600; % 0.05 / 0.05
    lambda1 = 0.00/3600; % 0.20 / 0.50
    omega   = 0.00/3600; % 0.03 / 0.12
    % plot_potentials(Vmin,alpha) % if you want to take a look at Morse and
    % LJ given the three variables above.
    
    population = cell(N,1);
    observed_cells = cell(N,1);
    Omega = scale*[0 0;1 0; 1 1; 0 1]';
    for i = 1:N
        population{i} = create_cell(0,[],[],15+rand(2,1)*20,0,0);
        observed_cells{i} = create_cell(0,[],[],zeros(2,Kobs),[],[]);
    end
    
    population = calculate_neighbours(population,1);
    N_old      = N;
    tagged_particle = population{1};
    tagged_particle.location = zeros(2,K);
    tagged_particle.local_density = zeros(1,K);
    %
    U_param = [Vmin alpha];
    
    % ALTERNATIVE PARAMETRISATION
    %a1 = 0.0016/dt;
    %a2 = 0.0001/dt;
    %k1 = 7;
    %k2 = 0.75;
    %l1 = 0.3;
    %l2 = 1.25;
    %U_param = [k1 a1 l1 k2 a2 l2];
    tic
    N_old = N;
    for k = 1:K
        % photograph the situation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if mod(k,freq) == 0
            obs = k/(freq);
            for i = 1:N
                observed_cells{i}.b_time            = ceil(population{i}.b_time/(freq));
                observed_cells{i}.d_time            = ceil(population{i}.d_time/(freq));
                 observed_cells{i}.parent            = population{i}.parent;
                 observed_cells{i}.children          = population{i}.children;
                observed_cells{i}.location(:,obs)   = population{i}.location;
%                 observed_cells{i}.neighbours(obs,:) = population{i}.neighbours;
            end
        end

        % add cells loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        newcomers   = cell(0);
        newcomers_o = cell(0);
        exiters     = cell(0);
        newindex    = N + 1;
        any_birth   = 0;
        for i = 1:min(N,N_old)
            tmp      = population{i}.birth_evolution;
            tmp2     = population{i}.birth_barrier;
            tmpD     = population{i}.death_evolution;
            tmpD2    = population{i}.death_barrier;
            % compute local density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            loc_dens = local_density(population,i,N,N_old,3.068891070360678);
            population{i}.local_density = loc_dens;
            % update internal states
            tmp  = tmp + dt*(1-loc_dens)*(lambda0+lambda1*loc_dens)*(1-tmp);
            tmpD = tmpD + dt*omega*(1-tmpD);
            % check if birth will occur %%%%%%%%%%%%%%%%
            if tmp2 < tmp
                any_birth = any_birth + 1;
                %disp(['cell number ' num2str(i) ' gave birth to cell number ' num2str(newindex) ' at time ' num2str(k)])
                crand     = 2*pi*rand;
                placetmp  = [cos(crand);sin(crand)]*(1 + 0.2*abs(randn));
                newcomers = [newcomers;...
                    create_cell(k,i,[],population{i}.location + placetmp,0,0)];
                tmp3      = population{i}.children;
                tmp3      = [tmp3 newindex];
                newindex  = newindex + 1;
                population{i}.children        = tmp3;
                population{i}.birth_evolution = 0;
                population{i}.birth_barrier   = rand;
                
                %%% newcomer copy for observation loop
                newcomers_o = [newcomers_o;...
                    create_cell(0,[],[],(population{i}.location + placetmp).*ones(2,Kobs),[],[])];
                %%% newcomer copy for observation loop
            else
                population{i}.birth_evolution = tmp;
            end
            % check if death will occur %%%%%%%%%%%%%%%%
            if tmpD2 < tmpD
                population{i}.d_time = k;
                population{i}.location = 10^30*[1;1];
                population{i}.birth_evolution = -10^20;
                population{i}.birth_barrier = 10^20;
                population{i}.death_evolution = -10^20;
            else
                population{i}.death_evolution = tmpD;
            end
        end
        population = [population; newcomers];
        observed_cells = [observed_cells; newcomers_o];
        if any_birth > 0 % recalculate neighbours if anyone got born and shit
            population = calculate_neighbours(population,k);
            N_old = N;        
        end
        
        % propagate loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N = length(population); % to account for newcomers
        if Vmin > 0 % if no interaction takes place, let's skip this part.
            F = interactions(population,U_param);
        else
            F = zeros(2,N);
        end% quite slow but of ok structure
        for i = 1:N
                if population{i}.d_time > k
                    population{i}.location(1) = max(min(population{i}.location(1)...
                        + dt*F(1,i) + sqrt(dt)*sigma*randn(1,1),scale),0);
                    population{i}.location(2) = max(min(population{i}.location(2)...
                        + dt*F(2,i) + sqrt(dt)*sigma*randn(1,1),scale),0);
                    % addition of the randomness at this stage is ok since the
                    % interactions have already been calculated.
                end
        end
        if mod(k,freq) == 0
            disp(['Hour: ' num2str(k/3600)])
            %%% RECALCULATE NEIGHBOURS HERE
            population = calculate_neighbours(population,k);
            N_old      = N;
            %%% RECALCULATE NEIGHBOURS HERE
            toc
        end
        
        
        tagged_particle.location(:,k) = population{1}.location;
        tagged_particle.local_density(k) = population{1}.local_density;
    end
    %%% %%% %%% SAVE OBSERVATIONS %%% %%% %%%
    for i = 1:N
        observed_cells{i} = rmfield(observed_cells{i},{'local_density','birth_evolution','birth_barrier'});
    end
    observed_cells{i+1} = freq;
    
    c = clock;
    c1 = num2str(c(1)); % år
    c2 = num2str(c(2)); % månad
    tmp = length(c2);
    if tmp == 1
        c2 = strcat('0',c2);
    end
    c3 = num2str(c(3)); % dag
    tmp = length(c3);
    if tmp == 1
        c3 = strcat('0',c3);
    end
    c4 = num2str(c(4)); % timme
    tmp = length(c4);
    if tmp == 1
        c4 = strcat('0',c4);
    end
    c5 = num2str(c(5)); % sekund
    tmp = length(c5);
    if tmp == 1
        c5 = strcat('0',c5);
    end
    filename = strcat('dataset_',c1,c2,c3,'_',c4,c5,'.mat');
    
    save(strcat('C:\Users\guslindw\Documents\MATLAB\Forskning, höst 2021\Datasets\',filename),'observed_cells')
    % observed_cells is the object we shall save
    % row 108 down to here is literally book keeping
    end
end




