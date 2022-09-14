%%
clc
clear
close all

%% Initalbedingungen

filepath.dir='V:\DDE\Exchange\Kundendienst_Bandende\Praktikant\PCAA\Bachelorarbeit_MatlabSimulink\Test_Files\';
filepath.filename='V350243_02061-02087_2022_KW8';
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

%% Messignale Zeitlichsynchronisieren und überschreiben
dt=0.0125;
sim("Messwerte_Sync.slx");
clearvars -except filepath ans dt

dm_abgas=ans.ASMod_dmEGFld_3(1/dt:end)/3600;%[kg/s]
VehV_v=ans.VehV_v(1/dt:end); %[km/h]
lambda=ans.EGSSig_ratLamEngineDs(1/dt:end); %[-]
T_abgas_nscr(:,1)=ans.UCatDsT_t(1/dt:end);   %[°C]
T_abgas_vscr=ans.SCR_tUCatUsT(1/dt:end);  %[°C]
p_umg=ans.EnvP_p(1/dt:end)*100; %[N/m2]
T_umg=ans.EnvT_t(1/dt:end); %[°C]
t=ans.t(1/dt:end); %[s]
clear ans


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
KondVerd_KLKF_mf_StSt_kgps = [0, 25, 100, 1000]/3600; 
Kond_KL_Faktor_ = [1, 0.3, 0.2, 0.1];  
deltaH_verd = 2.4e6; 
VerdRohr_offset = 0;  
Verd_KL_Faktor = [1,1; 0.3,0; 0.2,0; 0.1,0;] + VerdRohr_offset;
Verd_KF_Dampfdruck_StSt = [0, 1];  

%T_abgas_nscr(1,1:5)=T_umg(2,1);


% Rohr nach SCR
Rohr.mGas_imRohrabschnitt_kg = 0.005;  
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
% spez_feuchte(isnan(spez_feuchte))=0.0001;
% spez_feuchte(isinf(spez_feuchte))=0.0001;

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
alpha_wand_umg=polyval(Kennlinien.alpha_wand_umg,VehV_v); %[W/m2K]
rth_wand_umg=1./alpha_wand_umg+Rohr.LambdaRohr; %thermischer Gesamtwiderstand

dm_wasser_max=polyval(Kennlinien.saettigungskurve,T_abgas_nscr).*dm_abgas;%[kg/s]
for i=1:simlength
    if dm_wasser_max(i)<dm_wasser(i)
        dm_wasser(i)=dm_wasser_max(i);
    end
end
%%

T_Wand=[T_umg(2,1)-0.001 T_umg(2,1)-0.001 T_umg(2,1)-0.001 T_umg(2,1)-0.001];

time=t.*dt;
time=time';

T_abgas_nscr(:,2)=T_abgas_nscr;
spez_feuchte(:,2)=spez_feuchte;
T_umg(:,2)=T_umg;
p_umg(:,2)=p_umg;
dm_abgas(:,2)=dm_abgas;
VehV_v(:,2)=VehV_v;


T_abgas_nscr(:,1)=time;
spez_feuchte(:,1)=time;
T_umg(:,1)=time;
p_umg(:,1)=time;
dm_abgas(:,1)=time;
VehV_v(:,1)=time;



T_abgas_nscr=double(T_abgas_nscr);
spez_feuchte=double(spez_feuchte);
T_umg=double(T_umg);
p_umg=double(p_umg);
dm_abgas=double(dm_abgas);
VehV_v=double(VehV_v);

switch_TInit_Tstart_bei_true_sonst_TInitAusMessung=1;
T_start=T_umg(1,1);
