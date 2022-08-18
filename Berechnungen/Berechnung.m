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

%% Messignale Zeitlichsynchronisieren und überschreiben

sim("Messwerte_Sync.slx");

ASMod_dmEGFld_3=ans.ASMod_dmEGFld_3/3600;%[kg/s] 
VehV_v=ans.VehV_v; %[km/h] 
EGSSig_ratLamEngineDs=ans.EGSSig_ratLamEngineDs; %[-] 
UCatDsT_t=ans.UCatDsT_t;   %[°C]
SCR_tUCatUsT=ans.SCR_tUCatUsT;  %[°C]  
t=ans.t;
clear ans
%% Materialparameter

c_edelstahl = 500;                          % [J/kgK]       Wärmekapa. Edelstahl
c_Isomatte = 1000;                          % [J/kgK]       Wärmekapa. Isomatte der Bricks im Canning
deltaH_verd = 2.4e6;                        % [J/kg]        Wasser-Verdampfungsenthalpie
R_sLuft = 287;                              % [J/kgK]       Spez. Gaskonstante Luft 
R_sH2Og = 461;                              % [J/kgK]       Spez. Gaskonstante Wasserdampf
waermekapazitaet_wasser_gasfoermig = 2020;  % [J/kgK]
waermekapazitaet_wasser_fluessig = 4200;    % [J/kgK]
waermekapazitaet_luft = 1005;               % [J/kgK]    

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
Rohr.Waermekapazitaet = c_edelstahl;                           % [J/kgK]
Rohr.Waermeleitquerschnitt = 0.0003;                           % [m2]
Rohr.thermWidRohr = 0.0015 / Rohr.LambdaRohr;                 % [m2K/W]   thermischer Wid. der Rohrwand (Wandstärke / Lambda)
Rohr.thermWidIso = 0 / 0.15;                                   % [m2K/W]   darf 0 sein; thermischer Wid. der Isolierung (Wandstärke / Lambda)



%% Vorabberechnung
%Krafstoffmassenstrom
dm_diesel=ASMod_dmEGFld_3./EGSSig_ratLamEngineDs./lambda_st; %[kg/s] Gesamteinspritzmenge pro Sekunde
% Erstellung Dampfsättigungskurve
Kennlinien.temp_gas=[-20 -15 -10 -5 0 5 10 15 20 25 30 35 40 50 60 70 80 90 100]; %Temperatur der Luft/Abgases in °C
Kennlinien.wasser_gehalt=[0.9 1.4 2.1 3.3 4.8 6.8 9.4 12.8 17.3 23 30.3 39.6 51.1 80 125 195 290 420 590];% wassergehalt, in g/m3
Kennlinien.saettigungskurve=polyfit(Kennlinien.temp_gas,Kennlinien.wasser_gehalt,4);% Erstellung der Wasserkurve durch Polynome 4tn Grades
