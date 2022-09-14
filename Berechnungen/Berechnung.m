%% TODO:
% 
%
%delta zwischen q's in sim von franz und diesem file
%
%
%
%
%
%
%
%
%%

%%
clc
clear
close all

%% Initalbedingungen

filepath.dir='V:\DDE\Exchange\Kundendienst_Bandende\Praktikant\PCAA\Bachelorarbeit_MatlabSimulink\Test_Files\';
filepath.filename='V350243_02061-02087_2022_KW8_RESET';
filepath.fullpath=append(filepath.dir,filepath.filename,'.mat');


%% Messsignale laden

ASMod_dmEGFld_3=struct2array(load(filepath.fullpath,'ASMod_dmEGFld_3'));%[kg/h]
VehV_v=struct2array(load(filepath.fullpath,'VehV_v')); %[km/h]
%EGSSig_ratLamEngineDs=struct2array(load(filepath.fullpath,'EGSSig_ratLamEngineDs')); %[-]
EGSSig_ratLamEngineDs=struct2array(load(filepath.fullpath,'UEGO_rLamS1B1'));    %[-]
%UCatDsT_t=struct2array(load(filepath.fullpath,'UCatDsT_t'));   %[°C]
UCatDsT_t=struct2array(load(filepath.fullpath,'T_nSCR'));   %[°C]
SCR_tUCatUsT=struct2array(load(filepath.fullpath,'SCR_tUCatUsT'));  %[°C]
EnvP_p=struct2array(load(filepath.fullpath,'EnvP_p'));  %[hpa]
EnvT_t=struct2array(load(filepath.fullpath,'EnvT_t'));  %[°C]
reset(:,2)=struct2array(load(filepath.fullpath,'reset'));
reset(:,1)=EnvT_t(:,1);

%% Messignale Zeitlichsynchronisieren und überschreiben

dt=0.02;
sim("Messwerte_Sync.slx");
clearvars -except filepath ans dt 

dm_abgas=ans.ASMod_dmEGFld_3(1/dt:end)/3600;%[kg/s]
v_vehicle=ans.VehV_v(1/dt:end); %[km/h]
lambda=ans.EGSSig_ratLamEngineDs(1/dt:end); %[-]
T_abgas_nscr(:,1)=ans.UCatDsT_t(1/dt:end);   %[°C]
T_abgas_vscr=ans.SCR_tUCatUsT(1/dt:end);  %[°C]
p_umg=ans.EnvP_p(1/dt:end)*100; %[N/m2]
T_umg=ans.EnvT_t(1/dt:end); %[°C]
t=ans.t(1/dt:end); %[s]
reset=ans.reset(1/dt:end); %[s]

clear ans
for i=1:length(t)
T_umg(i,1)=10;
end


%% Materialparameter

mat_cons.c_edelstahl = 500;                          % [J/kgK]       Wärmekapa. Edelstahl
mat_cons.c_Isomatte = 1000;                          % [J/kgK]       Wärmekapa. Isomatte der Bricks im Canning
mat_cons.deltaH_verd = 2.4e6;                        % [J/kg]        Wasser-Verdampfungsenthalpie
mat_cons.R_sLuft = 287;                              % [J/kgK]       Spez. Gaskonstante Luft
mat_cons.R_sH2Og = 461;                              % [J/kgK]       Spez. Gaskonstante Wasserdampf
mat_cons.waermekapazitaet_wasser_gasfoermig = 2020;  % [J/kgK]
mat_cons.waermekapazitaet_wasser_fluessig = 4200;    % [J/kgK]
mat_cons.waermekapazitaet_luft = 1005;               % [J/kgK]

lambda_st=14.7;
WuezRohr_Dynamikfaktor=2;
m_KonWasser=[0 0 0 0];
%T_abgas_nscr(1,1:5)=T_umg(2,1);


% Rohr nach SCR
Rohr.AnzahlZonen = 4;                                          % [-]       Anzahl der Rohrsegmente für 1D-Diskretisierung
Rohr.Laenge_Abschnitt = 1.5 / Rohr.AnzahlZonen;               % [m]
Rohr.Mantelflaeche_innen = 0.075 * 3.14 * Rohr.Laenge_Abschnitt;          % [m2]
Rohr.Mantelflaeche_aussen = 0.080 * 3.14 * Rohr.Laenge_Abschnitt;         % [m2]
Rohr.Masse_Wasser_start = 0;                                   % [kg]
Rohr.StroemungsdurchmesserEff = 0.075;                         % [m]
Rohr.LambdaRohr = 15;                                          % [W/mK]
Rohr.MasseRohr = 1.8 / Rohr.AnzahlZonen;                      % [kg]
Rohr.Waermekapazitaet = mat_cons.c_edelstahl;                           % [J/kgK]
Rohr.Waermeleitquerschnitt = Rohr.StroemungsdurchmesserEff^2*pi/4;                           % [m2]
Rohr.thermWidRohr = 0.0015 / Rohr.LambdaRohr;                 % [m2K/W]   thermischer Wid. der Rohrwand (Wandstärke / Lambda)
Rohr.thermWidIso = 0 / 0.15;                                   % [m2K/W]   darf 0 sein; thermischer Wid. der Isolierung (Wandstärke / Lambda)
Rohr.Querschnitt= Rohr.StroemungsdurchmesserEff^2*pi;      %[m2]
Rohr.Dicke=0.005;                                           %[m]


%% Vorabberechnung

% Simulationslänge
simlength=length(t);

% Krafstoffmassenstrom, Wasserentsthehung und spez.Feucht
dm_diesel=dm_abgas./lambda/lambda_st; %[kg/s] Gesamteinspritzmenge pro Sekunde
% dm_diesel(isnan(dm_diesel))=0.000001;
% dm_diesel(isinf(dm_diesel))=0.000001;
dm_wasser=dm_diesel.*1.05; %[kg/s] Wasserentstehung durch Verbrennung
for i=1:simlength
    if dm_abgas(i)~=0
        spez_feuchte(i)=dm_wasser(i)/dm_abgas(i); %[kg/kg]
    else
        spez_feuchte(i)=0;
    end
end
spez_feuchte=spez_feuchte';

%Abgasmassestrom in Abgasvolumenstrom

dv_abgas=dm_abgas.*(T_abgas_nscr+273)*287/100000;

% Erstellung Dampfsättigungskurve
Kennlinien.temp_gas=[-20 -15 -10 -5 0 5 10 15 20 25 30 35 40 50 60 70 80 90 100]; %Temperatur der Luft/Abgases in °C
Kennlinien.wasser_gehalt=[0.9 1.4 2.1 3.3 4.8 6.8 9.4 12.8 17.3 23 30.3 39.6 51.1 80 125 195 290 420 590]/1000;% Wassergehalt, in kg/m3
Kennlinien.saettigungskurve=polyfit(Kennlinien.temp_gas,Kennlinien.wasser_gehalt,4);% Erstellung der Wasserkurve durch Polynome 4tn Grades

% Polyonom Wärmeübertraung Wand zur Umgebung
Kennlinien.v_alpha_w_u=[8, 18, 32, 66, 116, 200, 280];%W/m2K
Kennlinien.v_w_u=[0, 10, 20, 50, 100, 200, 300];%/km/h
Kennlinien.alpha_wand_umg=polyfit(Kennlinien.v_w_u,Kennlinien.v_alpha_w_u,3);%W/m2K

% Abgasgeschwindigkeit
v_abgas=(dm_abgas.*(T_abgas_nscr+273.15)./p_umg)*mat_cons.R_sLuft/Rohr.Querschnitt; %[m/s]

% Thermische Größen
alpha_gas_wand=(4.13 + 0.23*T_abgas_nscr./100 - 0.0077*(T_abgas_nscr./100).^2).*v_abgas.^0.75/Rohr.StroemungsdurchmesserEff^0.25;%[W/m2K]
alpha_gas_wand(isnan(alpha_gas_wand))=0;
alpha_wand_umg=polyval(Kennlinien.alpha_wand_umg,v_vehicle); %[W/m2K]
rth_wand_umg=1./alpha_wand_umg+Rohr.LambdaRohr; %thermischer Gesamtwiderstand

dm_wasser_max=polyval(Kennlinien.saettigungskurve,T_abgas_nscr).*dv_abgas;%[kg/s]
for i=1:simlength
    if dm_wasser_max(i)<dm_wasser(i)
        dm_wasser(i)=dm_wasser_max(i);
    end
end

T_Wand=[T_umg(2,1)-0.001 T_umg(2,1)-0.001 T_umg(2,1)-0.001 T_umg(2,1)-0.001];
Q_kondensation=0;

%T_Wand=[10 10 10 10; 10 10 10 10];
%% Iterative Kondeswassermengenbestimmung
tic
t_sim=28000; %Dauer der Sim in s
Q_waermeleitung(t_sim/dt,4)=0;
Q_gas_wand(t_sim/dt,4)=0;
Q_wand_umg(t_sim/dt,4)=0;
m_KonWasser(t_sim/dt,4)=0;
m_KonWasserGes(1,t_sim/dt)=0;
T_Wand(t_sim/dt,4)=0;
T_Wand(1,1:4)=[T_umg(2,1)-0.001 T_umg(2,1)-0.001 T_umg(2,1)-0.001 T_umg(2,1)-0.001];
for x=2:t_sim/dt
    %1400/dt
    if  x==364
        disp('test')
    end
        if reset(x)~=1

    for i=1:Rohr.AnzahlZonen
        if i==1
            Q_waermeleitung(x,i)=Rohr.LambdaRohr*Rohr.Querschnitt/Rohr.Laenge_Abschnitt*(T_abgas_nscr(x-1,i)-T_Wand(x-1,i)); %[W]
        else
            Q_waermeleitung(x,i)=Rohr.LambdaRohr*Rohr.Querschnitt/Rohr.Laenge_Abschnitt*(T_Wand(x-1,i-1)-T_Wand(x-1,i)); %[W]
        end
        Q_gas_wand(x,i)=(T_abgas_nscr(x-1,i)-T_Wand(x-1,i))*alpha_gas_wand(x,1)*Rohr.Mantelflaeche_innen*WuezRohr_Dynamikfaktor; %[W]
        Q_wand_umg(x,i)= Rohr.Mantelflaeche_aussen*(T_Wand(x-1,i)-T_umg(x-1))/(Rohr.thermWidRohr+1/alpha_wand_umg(x));


        if dm_abgas(x-1)~=0
            T_abgas_nscr(x-1,i+1)=T_abgas_nscr(x-1,i)-(Q_gas_wand(x,i)/((spez_feuchte(x-1)*mat_cons.waermekapazitaet_wasser_gasfoermig+mat_cons.waermekapazitaet_luft)*dm_abgas(x-1))*dt);
            dm_KonWasser=(dm_wasser(x,i)-polyval(Kennlinien.saettigungskurve,T_Wand(x-1,i))*dv_abgas(x-1))*dt;
            if dm_KonWasser<0
                dm_KonWasser=0;
            end

            m_KonWasser(x,i)=m_KonWasser(x-1,i)+dm_KonWasser;

            if dm_KonWasser>0
                Q_kondensation=2088*dm_KonWasser*1000;
            else
                Q_kondensation=0;
            end


            dm_wasser(x,i+1)=dm_wasser(x,i)-m_KonWasser(x,i);
        else
            T_abgas_nscr(x-1,i+1)=T_abgas_nscr(x-1,i);
            m_KonWasser(x,i)=m_KonWasser(x-1,i);
            dm_wasser(x,i+1)=dm_wasser(x,i);
          
        end




        T_Wand(x,i)=T_Wand(x-1,i)+(Q_kondensation+Q_waermeleitung(x,i)+Q_gas_wand(x,i)-Q_wand_umg(x,i))*dt/mat_cons.c_edelstahl/Rohr.MasseRohr;
    end
    m_KonWasserGes(x)=sum((m_KonWasser(x,1:4)));
    else
        T_Wand(x,1:4)=T_umg(x,1);
        m_KonWasserGes(x)=m_KonWasserGes(x-1);
        m_KonWasser(x,1)=m_KonWasser(x-1,1);
        m_KonWasser(x,2)=m_KonWasser(x-1,2);
        m_KonWasser(x,3)=m_KonWasser(x-1,3);
        m_KonWasser(x,4)=m_KonWasser(x-1,4);
        T_abgas_nscr(x-1,2:5)=T_abgas_nscr(x-1,1);
        
    end
end
toc
%Q_waermeleitung(x,i)+
%Q_kond(x,i)=
%       Q_verdunst


%% Testereien

% figure
% plot(t(1:13900),Q_waermeleitung(1:13900,2))
% grid on
% legend('Q_Wärmeleitung')
% figure
% plot(t(1:x),T_abgas_nscr(1:x,:))
% grid on
% legend('T_abgas')
figure
plot(t(1:x),dm_wasser(1:x))
grid on
legend('dm_wasser')
figure
plot(t(1:x),Q_gas_wand(1:x,:))
grid on
legend('Q_gas_wand')
% figure
% plot(t(1:x),spez_feuchte(1:x))
% grid on
% legend('spez_feuchte')
figure
plot(t(1:x),T_Wand(1:x,:))
hold on
legend('T_Wand')
figure
plot(t(1:x),T_abgas_nscr(1:x,:))
grid on
legend('T_abgas')
% figure
% plot(t(1:x),m_KonWasser(1:x,1:4))
% grid on
% legend('m_KonWasser')
figure
plot(t(1:x),m_KonWasserGes(1:x))
grid on
legend('m_KonWasserGes')
figure
plot(t(1:x),Q_wand_umg(1:x,:))
grid on
legend('Q_wand_umg')
figure
plot(0:100,polyval(Kennlinien.saettigungskurve,(0:100)))


