% Allgemein, WUEZ
SensWuez_v_mps = [0, 20, 40, 80, 200, 600, 1500]/3.6;                   % m/s
SensWuez_Wuez_Wpm2pK = [10, 207, 278, 383, 600, 1060, 1700];            % W/m2K  mit GT-Power Thermocouplemodell gerechnet

%% SENSOR Thermoelement NickelChromNickel 1.5mm Fühler 
Thermoelement_NiCrNi_1mm5.dsensor = 0.0015;                             % m     Sensordurchmesser in m
Thermoelement_NiCrNi_1mm5.lsensor = 0.02;                               % m     Sensorlänge(Teil für Simulation, nicht gesamtlänge) in m
Thermoelement_NiCrNi_1mm5.cpsensor = 444 ;                              % J/kgK c von Mantelmaterial(Inconel 600)
Thermoelement_NiCrNi_1mm5.msensor = 1.1e-4;                             % kg    Masse des Thermoelements (zum Kalibrieren verwendet)in kg
Thermoelement_NiCrNi_1mm5.ASensor = Thermoelement_NiCrNi_1mm5.dsensor * pi * Thermoelement_NiCrNi_1mm5.lsensor;     %Wärmeübergangsfläche aus Sensorgrößen in m2
%Kennfeld Stromungsgeschwindigkeitsabh. Wärmeübergangszahl
%Thermoelement_NiCrNi_1mm5.wuez_v = [0.08373016, 0.6925936, 2.0955987, 4.2046504, 8.419387, 16.83488, 25.250525, 29.458845, 42.08447];
%Die Zahlen hier können nicht stimmen, durch Standard ersetzt; Thermoelement_NiCrNi_1mm5.wuez_alpha = [63.48279,102.79919, 157.99336, 210.85403, 288.03427, 401.5474, 420, 450, 500];
Thermoelement_NiCrNi_1mm5.wuez_v = SensWuez_v_mps;
Thermoelement_NiCrNi_1mm5.wuez_alpha = SensWuez_Wuez_Wpm2pK;

%% Denso NTC Abgastemperatursensor
%Sensordurchmesser in m
Denso_NTC.dsensor = 0.002;
%Sensorlänge(Teil für Simulation, nicht gesamtlänge) in m
Denso_NTC.lsensor = 0.02;
%cp von Mantelmaterial(Inconel 600)
Denso_NTC.cpsensor = 444 ;
%Masse des Thermoelements (zum Kalibrieren verwendet)in kg
Denso_NTC.msensor = 5e-4;
%Wärmeübergangsfläche aus Sensorgrößen in m2
Denso_NTC.ASensor = Denso_NTC.dsensor*pi*Denso_NTC.lsensor;
%Kennfeld Stromungsgeschwindigkeitsabh. Wärmeübergangszahl
Denso_NTC.wuez_v = SensWuez_v_mps;
Denso_NTC.wuez_alpha = SensWuez_Wuez_Wpm2pK;

%% Denso NTC Abgastemperatursensor
%Sensordurchmesser in m
Sensata_PT1000.dsensor = 0.002;
%Sensorlänge(Teil für Simulation, nicht gesamtlänge) in m
Sensata_PT1000.lsensor = 0.02;
%cp von Mantelmaterial(Inconel 600)
Sensata_PT1000.cpsensor = 444 ;
%Masse des Thermoelements (zum Kalibrieren verwendet)in kg
Sensata_PT1000.msensor = 7.3e-4;
%Wärmeübergangsfläche aus Sensorgrößen in m2
Sensata_PT1000.ASensor = Sensata_PT1000.dsensor*pi*Sensata_PT1000.lsensor;
%Kennfeld Stromungsgeschwindigkeitsabh. Wärmeübergangszahl
Sensata_PT1000.wuez_v = SensWuez_v_mps;
Sensata_PT1000.wuez_alpha = SensWuez_Wuez_Wpm2pK;