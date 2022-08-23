%% TODO:
% Wärmeenergiebilanzen aufstellen
%
%
%
%
%
%
%
%
%
%

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
dt=0.001;
sim("Messwerte_Sync.slx");
clearvars -except filepath ans dt

dm_abgas=ans.ASMod_dmEGFld_3/3600;%[kg/s]
v_vehicle=ans.VehV_v; %[km/h]
lambda=ans.EGSSig_ratLamEngineDs; %[-]
T_abgas_nscr(:,1)=ans.UCatDsT_t;   %[°C]
T_abgas_vscr=ans.SCR_tUCatUsT;  %[°C]
p_umg=ans.EnvP_p*100; %[N/m2]
T_umg=ans.EnvT_t; %[°C]
t=ans.t; %[s]
clear ans
i=1;
while dm_abgas(i)==0
    dm_abgas(i)=[];
    lambda(i)=[];
    p_umg(i)=[];
    t(i)=[];
    T_abgas_nscr(i)=[];
    T_abgas_vscr(i)=[];
    T_umg(i)=[];
    v_vehicle(i)=[];
end
dm_abgas=dm_abgas+0.000001;
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
spez_feuchte=dm_wasser./dm_abgas; %[kg/kg]
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
% alpha_gas_wand(isnan(alpha_gas_wand))=0;
alpha_wand_umg=polyval(Kennlinien.alpha_wand_umg,v_vehicle); %[W/m2K]
rth_wand_umg=1./alpha_wand_umg+Rohr.LambdaRohr; %thermischer Gesamtwiderstand

dm_wasser_max=polyval(Kennlinien.saettigungskurve,T_abgas_nscr).*dm_abgas;%[kg/s]
for i=1:simlength
    if dm_wasser_max(i)<dm_wasser(i)
        dm_wasser(i)=dm_wasser_max(i);
    end
end
T_Wand=[T_umg(2,1)+0.001 T_umg(2,1)+0.001 T_umg(2,1)+0.001 T_umg(2,1)+0.001];

%% Iterative Kondeswassermengenbestimmung
tic
for x=2:13984
    if x==8500
    disp('test')
    end

    for i=1:Rohr.AnzahlZonen
        if i==1
            Q_waermeleitung(x,i)=Rohr.LambdaRohr*Rohr.Querschnitt/Rohr.Laenge_Abschnitt*(T_abgas_nscr(x-1,i)-T_Wand(x-1,i)); %[W]
        else
            Q_waermeleitung(x,i)=Rohr.LambdaRohr*Rohr.Querschnitt/Rohr.Laenge_Abschnitt*(T_Wand(x-1,i-1)-T_Wand(x-1,i)); %[W]
        end
        Q_gas_wand(x,i)=(T_abgas_nscr(x-1,i)-T_Wand(x-1,i))*alpha_gas_wand(x-1,1)*Rohr.Mantelflaeche_innen*WuezRohr_Dynamikfaktor; %[W]

        Q_wand_umg(x,i)= Rohr.Mantelflaeche_aussen*(T_Wand(x-1,i)-T_umg(x-1)/(Rohr.thermWidRohr+alpha_wand_umg(x-1)));
  
        T_abgas_nscr(x-1,i+1)=T_abgas_nscr(x-1,i)-(Q_gas_wand(x-1,i)/(spez_feuchte(x-1)*mat_cons.waermekapazitaet_wasser_gasfoermig*dm_wasser(x-1)*dm_abgas(x-1)*dt+mat_cons.waermekapazitaet_luft*dm_abgas(x-1)*dt));
        
        m_KonWasser(x,i)=m_KonWasser(x-1,i)+(dm_wasser(x,i)-polyval(Kennlinien.saettigungskurve,T_abgas_nscr(x-1,i+1))*dm_abgas(x-1)*dt);
%         m_KonWasser(isnan(m_KonWasser))=0;
        dm_wasser(x,i+1)=dm_wasser(x,i)-m_KonWasser(x,i)/(dm_abgas(x-1)*dt);
        T_Wand(x,i)=T_Wand(x-1,i)+(Q_waermeleitung(x,i)+Q_gas_wand(x,i)-Q_wand_umg(x,i))*dm_abgas(x-1)*dt/mat_cons.c_edelstahl;
    end
    
end
toc

%Q_kond(x,i)=
%       Q_verdunst


%% Testereien
plot(t(1:13900),alpha_wand_umg(1:13900,:))
