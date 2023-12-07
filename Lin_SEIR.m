%%% Author: Jolyne Lin (based on Prof de Pillis's SIR model code)
%%% Class: MathBio 119
%%% Date: Dec 8, 2023
%%% Description: SEIR Model for Trypanosomiasis
%%% Reference: 
%%% Meisner, Julianne et al. "A mathematical model for evaluating the role
%%% of trypanocide treatment of cattle in the epidemiology and control of
%%% Trypanosoma brucei rhodesiense and T. b. gambiense sleeping sickness in
%%% Uganda." Parasite epidemiology and control vol. 5 e00106. 16 Apr. 2019,
%%% doi:10.1016/j.parepi.2019.e00106

%%%Time range
t0 = 0;          % initial time in days
tf = 350;        % final time in days
timerange = t0:tf;

%%%Parameters
%HUMAN
B_h = (43/1000)/365; % birth rate
mu_h = (1/58.45)/365; % mortality rate
i_h = 1/12; % 1/incubation period
beta_h = 0.18; % probability of bite is of a human (0.18 or 0.3)
beta_vh = 0.0076; % probability of transmission from fly to human
gamma_ha = 1/20; % duration of stage I with care seeking
gamma_h1 = 1/526; % stage I duration
gamma_h2 = 1/252; % stage II duration
delta_h = 1/50; % duration of immunity
A_h = 0.3; % probability a given stage I HAT case seeks care (0.1, 0.3, or 0.5)
s = 0.87; % sensitivity of diagnostics
epsilon_1 = 0.94; % probability of treatment success, stage I gHAT
epsilon_2 = 0.965; % probability of treatment success, stage II gHAT
mu_h2 = 1-(0.005*s*epsilon_2); % probability of death from stage II

%CATTLE
B_c = 0.29/365; % birth rate
mu_c = 0.29/365; % mortality rate
i_c = 1/12; % incubation
beta_c = 0.38; % probability of bite is of an ox
beta_vc = 0.1345; % probability of transmission from fly to ox
gamma_c = 1/100; % duration of infection w/o treatment
gamma_ct = 12/365.25; % duration of infection w treatment
delta_c = 1/50; % duration of immunity
T_c = 0.5; % probability of treatment (0, 0.5, or 0.75)

%FLIES
B_f = 0.05; % birth rate
mu_f = 1/30; % adult mortality rate
p = 1/20; % rate of expupation
i_f = 1/25; % incubation
alpha = 0.33; % bite rate
beta_f = 0.2; % probability of transmission from human/ox to fly

parlist = [B_h, mu_h, i_h, beta_h, beta_vh, gamma_ha, gamma_h1, gamma_h2, delta_h, ...
           A_h, s, epsilon_1, mu_h2, mu_h2, B_c, mu_c, i_c, beta_c, beta_vc, gamma_c, ...
           gamma_ct, delta_c, T_c, B_f, mu_f, p, i_f, alpha, beta_f];

%Initial Conditions
S_h_0 = 5148882-604;
E_h_0 = 0;
I1_h_0 = 604;
I2_h_0 = 0;
R_h_0 = 0;
N_h_0 = S_h_0+E_h_0+I1_h_0+I2_h_0+R_h_0;
S_c_0 = 1338912*(1-0.125);
E_c_0 = 0;
I_c_0 = 1338912*0.125;
R_c_0 = 0;
N_c_0 = S_c_0+E_c_0+I_c_0+R_c_0;
pupae_0 = 0;
teneral_0 = 5148882*2;
NS_0 = 0;
E_f_0 = 0;
I_f_0 = 0;
N_f_0 = teneral_0+NS_0+E_f_0+I_f_0;

y0 = [S_h_0; E_h_0; I1_h_0; I2_h_0; R_h_0;...
      S_c_0; E_c_0; I_c_0; R_c_0;...
      pupae_0; teneral_0; NS_0; E_f_0; I_f_0];

[tout,yout] = ode45(@seir_model, timerange, y0, [], parlist);

%Plot
figure;
plot(tout,yout(:,1),'r',tout,yout(:,2),'b', tout,(yout(:,3)+yout(:,4)),'g',tout,yout(:,5),'c','LineWidth',3);
xlabel('Time in Days','FontSize',14);
ylabel('Number in Population','FontSize',14);
legend('Susceptible','Exposed','Infected','Recovered','FontSize',14);
xlim([0 tf]);
title('Human Population at A_H = 0.3 and T_C = 0.5');

function ydot = seir_model(t,y,params)
    %seir_model: Local function that constructs the SEIR model from parameters
    %and initial conditions

    %%%Extract population values
    %HUMAN
    S_h = y(1);
    E_h = y(2);
    I1_h = y(3);
    I2_h = y(4);
    R_h = y(5);
    N_h = S_h+E_h+I1_h+I2_h+R_h;

    %CATTLE
    S_c = y(6);
    E_c = y(7);
    I_c = y(8);
    R_c = y(9);
    N_c = S_c+E_c+I_c+R_c;

    %FLIES
    pupae = y(10);
    teneral = y(11);
    NS = y(12);
    E_f = y(13);
    I_f = y(14);
    N_f = teneral+NS+E_f+I_f;
        
    %%%Extract parameters
    %HUMAN
    B_h = params(1);
    mu_h = params(2);
    i_h = params(3);
    beta_h = params(4);
    beta_vh = params(5);
    gamma_ha = params(6);
    gamma_h1 = params(7);
    gamma_h2 = params(8);
    delta_h = params(9);
    A_h = params(10);
    s = params(11);
    epsilon_1 = params(12);
    mu_h2 = params(13);
    
    %CATTLE
    B_c = params(14);
    mu_c = params(15);
    i_c = params(16);
    beta_c = params(17);
    beta_vc = params(18);
    gamma_c = params(19);
    gamma_ct = params(20);
    delta_c = params(21);
    T_c = params(22);
    
    %FLIES
    B_f = params(23);
    mu_f = params(24);
    mu_p = 0.3; % probability of pupal mortality
    p = params(25);
    i_f = params(26);
    alpha = params(27);
    beta_f = params(28);
    
    %%%Useful equations
    cattle_fly = N_f/N_c;
    fly_human = N_f/N_h;
    O_h = max([(0.18/(0.18+0.38))*(N_h/N_c) 1]); % first bite human
    lambda_h = O_h*alpha*beta_h*beta_vh*fly_human*(I_f/N_f); % force of infection for human
    lambda_c = alpha*beta_c*beta_vc*cattle_fly*(I_f/N_f); % force of infection for ox
    under_1_day_old = 1-exp(-(mu_f+alpha*(beta_h+beta_c))); % proportion of unfed flies which are less than one day old
    phi = (1-under_1_day_old)+alpha*under_1_day_old; % rate at which teneral flies become non-susceptible
    lambda_f = alpha*beta_f*(beta_h*(I1_h/N_h)+beta_c*(I_c/N_c)); % force of infection for fly

    %%%Set up differential equations
    %HUMAN
    S_h_dot = B_h*N_h-lambda_h*S_h-mu_h*S_h+delta_h*R_h;
    E_h_dot = lambda_h*S_h-i_h*E_h-mu_h*E_h;
    I1_h_dot = i_h*E_h-gamma_h1*I1_h*(1-A_h*s*epsilon_1)-A_h*s*epsilon_1*gamma_ha*I1_h-mu_h*I1_h;
    I2_h_dot = gamma_h1*I1_h-gamma_h2*I2_h;
    R_h_dot = (1-mu_h2)*gamma_h2*I2_h+A_h*s*epsilon_1*gamma_ha*I1_h-delta_h*R_h-mu_h*R_h;

    %CATTLE
    S_c_dot = B_c*N_c-lambda_c*S_c-mu_c*S_c+delta_c*R_c;
    E_c_dot = lambda_c*S_c-i_c*E_c-mu_c*E_c;
    I_c_dot = i_c*E_c-gamma_c*I_c*(1-T_c)-gamma_ct*I_c*T_c;
    R_c_dot = T_c*gamma_ct*I_c-delta_c*R_c-mu_c*R_c;

    %FLIES
    pupae_dot = B_f*N_f-p*pupae;
    teneral_dot = (1-mu_p)*p*pupae-phi*teneral-lambda_f*teneral-mu_f*teneral;
    NS_dot = phi*teneral-mu_f*NS;
    E_f_dot = lambda_f*teneral-i_f*E_f-mu_f*E_f;
    I_f_dot = i_f*E_f-mu_f*I_f;

    
    ydot = [S_h_dot;
            E_h_dot;
            I1_h_dot;
            I2_h_dot;
            R_h_dot;
            S_c_dot;
            E_c_dot;
            I_c_dot;
            R_c_dot;
            pupae_dot;
            teneral_dot;
            NS_dot;
            E_f_dot;
            I_f_dot];
end

