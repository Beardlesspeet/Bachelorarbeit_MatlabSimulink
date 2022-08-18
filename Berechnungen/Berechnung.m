%% TODO:
%%
clc
clear
close all

%% Initalbedingungen

filepath.dir='V:\DDE\Exchange\Kundendienst_Bandende\Praktikant\PCAA\Bachelorarbeit_MatlabSimulink\Test_Files\';
filepath.filename='V350243_02061-02087_2022_KW8';
filepath.fullpath=append(filepath.dir,filepath.filename,'.mat');
stepsize=1;

%% Messsignale laden

ASMod_dmEGFld_3=struct2array(load(filepath.fullpath,'ASMod_dmEGFld_3'));%[kg/h] 
VehV_v=struct2array(load(filepath.fullpath,'VehV_v')); %[km/h] 
%EGSSig_ratLamEngineDs=struct2array(load(filepath.fullpath,'EGSSig_ratLamEngineDs')); %[-]  
EGSSig_ratLamEngineDs=struct2array(load(filepath.fullpath,'UEGO_rLamS1B1'));    %[-]  
%UCatDsT_t=struct2array(load(filepath.fullpath,'UCatDsT_t'));   %[°C]  
UCatDsT_t=struct2array(load(filepath.fullpath,'T_nSCR'));   %[°C]  
SCR_tUCatUsT=struct2array(load(filepath.fullpath,'SCR_tUCatUsT'));  %[°C] 
EnvP_p=struct2array(load(filepath.fullpath,'EnvP_p'));  %[hpa] 

%% Messignale Zeitlichsynchronisieren und überschreiben

sim("Messwerte_Sync.slx");
clearvars -except filepath ans

dm_abgas=ans.ASMod_dmEGFld_3/3600;%[kg/s] 
v_vehicle=ans.VehV_v; %[km/h] 
lambda=ans.EGSSig_ratLamEngineDs; %[-] 
t_abgas_nscr=ans.UCatDsT_t;   %[°C]
t_abgas_vscr=ans.SCR_tUCatUsT;  %[°C] 
p_umg=ans.EnvP_p*100; %[N/m2]
t=ans.t;
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
Rohr.Waermeleitquerschnitt = 0.0003;                           % [m2]
Rohr.thermWidRohr = 0.0015 / Rohr.LambdaRohr;                 % [m2K/W]   thermischer Wid. der Rohrwand (Wandstärke / Lambda)
Rohr.thermWidIso = 0 / 0.15;                                   % [m2K/W]   darf 0 sein; thermischer Wid. der Isolierung (Wandstärke / Lambda)
Rohr.Querschnitt= Rohr.StroemungsdurchmesserEff^2*pi; %[m2]


%% Vorabberechnung
%Krafstoffmassenstrom
dm_diesel=dm_abgas./lambda./lambda_st; %[kg/s] Gesamteinspritzmenge pro Sekunde

% Erstellung Dampfsättigungskurve
Kennlinien.temp_gas=[-20 -15 -10 -5 0 5 10 15 20 25 30 35 40 50 60 70 80 90 100]; %Temperatur der Luft/Abgases in °C
Kennlinien.wasser_gehalt=[0.9 1.4 2.1 3.3 4.8 6.8 9.4 12.8 17.3 23 30.3 39.6 51.1 80 125 195 290 420 590];% Wassergehalt, in g/m3
Kennlinien.saettigungskurve=polyfit(Kennlinien.temp_gas,Kennlinien.wasser_gehalt,4);% Erstellung der Wasserkurve durch Polynome 4tn Grades

%Polyonom Wärmeübertraung Wand zur Umgebung
Kennlinien.v_alpha_w_u=[8, 18, 32, 66, 116, 200, 280];%W/m2K
Kennlinien.v_w_u=[0, 10, 20, 50, 100, 200, 300];%/km/h
Kennlinien.alpha_wand_umg=polyfit(Kennlinien.v_w_u,Kennlinien.v_alpha_w_u,3);%W/m2K

%Abgasgeschwindigkeit
v_abgas=(dm_abgas.*(t_abgas_nscr+273.15)./p_umg)*mat_cons.R_sLuft/Rohr.Querschnitt; %[m/s]


