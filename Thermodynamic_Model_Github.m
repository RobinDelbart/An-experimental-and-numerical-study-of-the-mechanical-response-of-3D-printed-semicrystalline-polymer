clc
close all    
clear all

disp(' ')
disp('=================================================================')
disp(' Le Code des champions                                           ')
disp('=================================================================')
%% 

global Ti To a Profil Temp_mean_top Temp_mean_bot

%% Parametres du PB -- La ou tu peux jouer

% >>>>>> ATTENTION -- Changer Temperature piece et buse dans section ``plot_figure'' <<<<<<<

Data2Plot          = [130, 110, 115, 120]; % [deg]: lignes verticales pour les delta
Printing_Speed=30;
Flow_Rate_LH02=((0.2*1)+(0.1*0.1*pi))*Printing_Speed;
Flow_Rate_LH01=((0.1*1)+(0.05*0.1*pi))*Printing_Speed;
Flow_Rate_LH005=((0.05*1)+(0.025*0.1*pi))*Printing_Speed;
Ti  = 20  +273.15;  % [deg]:...Temperature piece
To  = 225 +273.15;  % [deg]:...Temperature buse 
a   = 0.105 *1e-6;  % [mm2/s]:.Thermal diffusivity
R1  = 1.75/2;       % [mm]:....Rayon cylindre 1
R2  = 1.375/2;      % [mm]:....Rayon cylindre 2
R3  = 1/2;          % [mm]:....Rayon cylindre 3
NB_Boudin_LH02=10;
NB_Boudin_LH01=20;
NB_Boudin_LH005=40;
L1=4.125;
L2=0.375;
L3=1.5;
% Rapide LH0.2 mm 
    tf1 = L1/(Flow_Rate_LH02/(pi*R1^2));      % [s]:.....Duree de passage dans cylindre 1
    tf2 = L2/(Flow_Rate_LH02/(pi*R2^2));       % [s]:.....Duree de passage dans cylindre 2
    tf3 = L3/(Flow_Rate_LH02/(pi*R3^2));      % [s]:.....Duree de passage dans cylindre 3
    h_top = 0.2;       % [mm] hauteur boudin top
    h_bot = h_top*NB_Boudin_LH02;  % [mm] hauteur boudin bottom
    tf_attente = 10;  % [s]  duree attente pour simu boudin sur boudin
% Moyen LH 0.1 mm
    tf1__ = L1/(Flow_Rate_LH01/(pi*R1^2));      % [s]:.....Duree de passage dans cylindre 1
    tf2__ = L2/(Flow_Rate_LH01/(pi*R2^2));       % [s]:.....Duree de passage dans cylindre 2
    tf3__ = L3/(Flow_Rate_LH01/(pi*R3^2));      % [s]:.....Duree de passage dans cylindre 3
    h_top__ = 0.1;        % [mm] hauteur boudin top
    h_bot__ = h_top__*NB_Boudin_LH01; % [mm] hauteur boudin bottom
    tf_attente__ = tf_attente;   % [s]  duree attente pour simu boudin sur boudin
% Lent LH 0.05 mm
    tf1_ = L1/(Flow_Rate_LH005/(pi*R1^2));      % [s]:.....Duree de passage dans cylindre 1
    tf2_ = L2/(Flow_Rate_LH005/(pi*R2^2));       % [s]:.....Duree de passage dans cylindre 2
    tf3_ = L3/(Flow_Rate_LH005/(pi*R3^2)) ;     % [s]:.....Duree de passage dans cylindre 3
    h_top_ = 0.05;       % [mm] hauteur boudin top
    h_bot_ = h_top__*NB_Boudin_LH005; % [mm] hauteur boudin bottom
    tf_attente_ = tf_attente;   % [s]  duree attente pour simu boudin sur boudin
    
 
for hide_rapide = []  
    %% Solveur -- La ou tu ne peux plus jouer

    tf = tf1 + tf2+ tf3;

    Temperature_centre = [];
    Dates_centre       = [];
    Position_centre    = [];

    %% I.- Cylindre 1
        R1 = R1*1e-3; %[m]
        NBptsR =50; NBptsT = 50; % Pour affiner le maillage de l'image
        x1 = 0:R1/(NBptsR-1):R1; % linspace(0,L,20);
        p1 = 0:4.1250/(NBptsT-1):4.125;
        t1 = 0:tf1/(NBptsT-1):tf1; % [linspace(0,0.02,20), linspace(0.5,1.5,10)];
        m=1; % Geometrie cylindrique
        options = odeset('RelTol',1e-5);
        sol1 = pdepe(m,@heatpde,@heatic,@heatbc,x1,t1,options);
    % Save data
    Temperature_centre = [Temperature_centre; sol1(:,1)];
    Dates_centre       = [Dates_centre        ,t1];
    Position_centre    = [Position_centre     ,p1];

    %% II.- Cylindres 2
    Profil.xp = x1;
    Profil.Tp = sol1(end,:);

        R2 = R2*1e-3; %[m]
        NBptsR =50; NBptsT = 50; % Pour affiner le maillage de l'image
        x2 = 0:R2/(NBptsR-1):R2; % linspace(0,L,20);
        p2 = 4.125:0.375/(NBptsT-1):4.125+0.375;
        t2 = 0:tf2/(NBptsT-1):tf2; % [linspace(0,0.02,20), linspace(0.5,1.5,10)];
        m=1; % Geometrie cylindrique
        options = odeset('RelTol',1e-5);
        sol2 = pdepe(m,@heatpde,@heaticProf,@heatbc,x2,t2,options);
    % Save data
    Temperature_centre = [Temperature_centre; sol2(:,1)];
    Dates_centre       = [Dates_centre        ,t2+tf1];
    Position_centre    = [Position_centre     ,p2];


    %% III.- Cylindres 3
    Profil.xp = x2;
    Profil.Tp = sol2(end,:);

        R3 = R3*1e-3; %[m]
        NBptsR =50; NBptsT = 50; % Pour affiner le maillage de l'image
        x3 = 0:R3/(NBptsR-1):R3; % linspace(0,L,20);
        p3 = 4.5:1.5/(NBptsT-1): 6;
        t3 = 0:tf3/(NBptsT-1):tf3; % [linspace(0,0.02,20), linspace(0.5,1.5,10)];
        m=1; % Geometrie cylindrique
        options = odeset('RelTol',1e-5);
        sol3 = pdepe(m,@heatpde,@heaticProf,@heatbc,x3,t3,options);
    % Save data
    Temperature_centre = [Temperature_centre;sol3(:,1)];
    Dates_centre       = [Dates_centre, t3+tf1+tf2];
    Position_centre    = [Position_centre     ,p3];
    
    Temperature_centre_rapide = Temperature_centre;
    Dates_centre_rapide       = Dates_centre;
    Position_centre_rapide    = Position_centre;
    
    %% IV.- Empilement de boudins
    sol3_rapide      = sol3(end,:);
    sol3_rapide_mean = mean(sol3(end,:));
    %sol3_rapide_mean = 167+273.15;
    Temp_mean_top = sol3_rapide_mean; % Initialisation temperature layer top a moyenne du profil de temp sorant de la buse
    Temp_mean_bot =   60+273.15;            %Ti;   %   "              "          "     bottom a temperature de la piece
    NBptsh = 800;   
    h4 = -h_bot:(h_bot+h_top)/(NBptsh-1):(h_top); % Discretisation du PB
    h4 = h4*1e-3;
    t4 = 0:0.01:tf_attente;
    t4 = [t4(1),0.001 ,t4(2:end)];
    m=0; % Geometrie euclidienne
    options = odeset('RelTol',1e-5);
    sol4 = pdepe(m,@heatpde,@heatic_boudin,@heatbc_boudin,h4,t4,options);
    
    % extraire donnee au temps donne
    t4_rapide   = t4;
    sol4_rapide = sol4-273.15;
    h4_rapide   = h4; 
    
    %% Analyse cristalisation
    % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    t4 = 0:0.01:tf_attente;  %             \
    t4 = [t4(1),t4(2:end)];  %  <__________/
    m=0; % Geometrie euclidienne
    options = odeset('RelTol',1e-5);
    % Temperature moyenne [K] du boudin a chaque increment de temps: 
    sol4 = pdepe(m,@heatpde,@heatic_boudin,@heatbc_boudin,h4,t4,options);
    sol4_deg = sol4 - 273.15;
    T_mean_rapide = mean(sol4_deg(:,ceil((NBptsh-NBptsh/(NB_Boudin_LH02+1))):end)');  
    % Interval de temps entre deux mesure de moyenne: 
    dt = 0.01; % [s]
    crist(1) = 0;
    %%
    cellules_r = [];
    t4_tot     = [];
    vol        = []; 
    Tc         = []; 
    for i =1:length(t4) 
        t4_tot = [t4_tot, t4(i)];
        Tc = [Tc, T_mean_rapide(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 2
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_rapide = mean(sol4_deg(:,800- 1*73-72:800- 1*73)'); 
    t4_last = t4_tot(end);  
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_rapide(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 3
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_rapide = mean(sol4_deg(:,800- 2*73-72:800- 2*73)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_rapide(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 4
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_rapide = mean(sol4_deg(:,800- 3*73-72:800- 3*73)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_rapide(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 5
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_rapide = mean(sol4_deg(:,800- 4*73-72:800- 4*73)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_rapide(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    
    %% Chauffage 6
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_rapide = mean(sol4_deg(:,800- 5*73-72:800- 5*73)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_rapide(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    
    %% Chauffage 7
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_rapide = mean(sol4_deg(:,800- 6*73-72:800- 6*73)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_rapide(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 8
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_rapide = mean(sol4_deg(:,800- 7*73-72:800- 7*73)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_rapide(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
      %%
    save Data_Rapide.mat        ...
         Position_centre_rapide ...
         Dates_centre_rapide    ...
         Temperature_centre_rapide ...
         tf ... 
         t4_rapide  sol4_rapide h4_rapide ...
          Data2Plot ...
         sol3_rapide  sol3_rapide_mean 
    %%
    t4_rapide    = t4;
    % Imprimer figure
    figure
    csvwrite('Time_LH02.csv',transpose(t4_tot))
    csvwrite('Temperature_LH02.csv',transpose(Tc))
    plot(t4_tot, vol*100)
        hXLabel = xlabel('Time [s]');
        hYLabel = ylabel('Cristalisation [\%]');
        set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca,'FontSize',10);
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
        set(gca,'Layer','Top');
%         set(gca, 'Position', [0.3, 0.3, 0.65, 0.65])
        set(gcf, 'PaperUnits', 'centimeters');
        x_width=4.45*2; y_width=4*2;  
        set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
        print('-painters','-dpdf','Cristalisation_LH_0_2')
    figure
    plot(t4_tot, Tc)
        hXLabel = xlabel('Time [s]');
        hYLabel = ylabel('Temprature [C]');
        set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca,'FontSize',10);
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
        set(gca,'Layer','Top');
%         set(gca, 'Position', [0.3, 0.3, 0.65, 0.65])
        set(gcf, 'PaperUnits', 'centimeters');
        x_width=4.45*2; y_width=4*2;  
        set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
        print('-painters','-dpdf','temp_LH_0_2')


end  
for hide_moyen  = [] 
    %% Solveur -- La ou tu ne peux plus jouer

    tf_moyen = tf1__ + tf2__+ tf3__;

    Temperature_centre__ = [];
    Dates_centre__       = [];
    Position_centre__    = [];

    %% I.- Cylindre 1
        R1 = R1*1e-3; %[m]
        NBptsR =50; NBptsT = 50; % Pour affiner le maillage de l'image
        x1 = 0:R1/(NBptsR-1):R1; % linspace(0,L,20);
        p1__ = 0:4.125/(NBptsT-1):4.125;
        t1__ = 0:tf1__/(NBptsT-1):tf1__; % [linspace(0,0.02,20), linspace(0.5,1.5,10)];
        m=1; % Geometrie cylindrique
        options = odeset('RelTol',1e-5);
        sol1__ = pdepe(m,@heatpde,@heatic,@heatbc,x1,t1__,options);
    % Save data
    Temperature_centre__ = [Temperature_centre__; sol1__(:,1)];
    Dates_centre__       = [Dates_centre__        ,t1__];
    Position_centre__    = [Position_centre__     ,p1__];

    %% II.- Cylindres 2
    Profil.xp = x1;
    Profil.Tp = sol1__(end,:);

        R2 = R2*1e-3; %[m]
        NBptsR =50; NBptsT = 50; % Pour affiner le maillage de l'image
        x2 = 0:R2/(NBptsR-1):R2; % linspace(0,L,20);
        p2__ = 4.125:0.375/(NBptsT-1): 4.5;
        t2__ = 0:tf2__/(NBptsT-1):tf2__; % [linspace(0,0.02,20), linspace(0.5,1.5,10)];
        m=1; % Geometrie cylindrique
        options = odeset('RelTol',1e-5);
        sol2__ = pdepe(m,@heatpde,@heaticProf,@heatbc,x2,t2__,options);
    % Save data
    Temperature_centre__ = [Temperature_centre__; sol2__(:,1)];
    Dates_centre__       = [Dates_centre__        ,t2__+tf1__];
    Position_centre__    = [Position_centre__     ,p2__];


    %% III.- Cylindres 3
    Profil.xp = x2;
    Profil.Tp = sol2__(end,:);

        R3 = R3*1e-3; %[m]
        NBptsR =50; NBptsT = 50; % Pour affiner le maillage de l'image
        x3 = 0:R3/(NBptsR-1):R3; % linspace(0,L,20);
        p3__ = 4.5:1.5/(NBptsT-1): 6;
        t3__ = 0:tf3__/(NBptsT-1):tf3__; % [linspace(0,0.02,20), linspace(0.5,1.5,10)];
        m=1; % Geometrie cylindrique
        options = odeset('RelTol',1e-5);
        sol3__ = pdepe(m,@heatpde,@heaticProf,@heatbc,x3,t3__,options);
    % Save data
    Temperature_centre__ = [Temperature_centre__;sol3__(:,1)];
    Dates_centre__       = [Dates_centre__, t3__+tf1__+tf2__];
    Position_centre__    = [Position_centre__     ,p3__];
    
    Temperature_centre_moyen = Temperature_centre__;
    Dates_centre_moyen       = Dates_centre__;
    Position_centre_moyen    = Position_centre__;
    
    %% IV.- Empilement de boudins
    sol3_moyen         = sol3__(end,:);
    x3                 = x3;
    sol3_moyen_mean    = mean(sol3__(end,:));
    %sol3_moyen_mean    = 225+273.15;
    Temp_mean_top = sol3_moyen_mean; % Initialisation temperature layer top a moyenne du profil de temp sorant de la buse
    Temp_mean_bot =   60+273.15;            %Ti;   %   "              "          "     bottom a temperature de la piece
    NBptsh = 800;   
    h4 = -h_bot__:(h_bot__+h_top__)/(NBptsh-1):(h_top__); % Discretisation du PB
    h4 = h4*1e-3;
    t4 = 0:0.01:tf_attente;
    m=0; % Geometrie euclidienne
    options = odeset('RelTol',1e-5);
    sol4__ = pdepe(m,@heatpde,@heatic_boudin,@heatbc_boudin,h4,t4,options);
    
    % extraire donnee au temps donne
    
    % extraire donnee au temps donne
    t4_moyen   = t4;
    sol4_moyen = sol4__-273.15;
    h4_moyen   = h4; 
    
    
 %% Analyse cristalisation
    % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    t4 = 0:0.01:tf_attente;  %             \
    t4 = [t4(1),t4(2:end)];  %  <__________/
    m=0; % Geometrie euclidienne
    options = odeset('RelTol',1e-5);
    % Temperature moyenne [K] du boudin a chaque increment de temps: 
    sol4 = pdepe(m,@heatpde,@heatic_boudin,@heatbc_boudin,h4,t4,options);
    sol4_deg = sol4 - 273.15;
    T_mean_moyen = mean(sol4_deg(:,762:end)');  
    % Interval de temps entre deux mesure de moyenne: 
    dt = 0.01; % [s]
    crist(1) = 0;
    %%
    cellules_r = [];
    t4_tot     = [];
    vol        = []; 
    Tc         = []; 
    for i =1:length(t4) 
        t4_tot = [t4_tot, t4(i)];
        Tc = [Tc, T_mean_moyen(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 2
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_moyen = mean(sol4_deg(:,800- 1*39-38:800- 1*39)'); 
    t4_last = t4_tot(end);  
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_moyen(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 3
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_moyen = mean(sol4_deg(:,800- 2*39-38:800- 2*39)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_moyen(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 4
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_moyen = mean(sol4_deg(:,800- 3*39-38:800- 3*39)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_moyen(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 5
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_moyen = mean(sol4_deg(:,800- 4*39-38:800- 4*39)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_moyen(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    
    %% Chauffage 6
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_moyen = mean(sol4_deg(:,800- 5*39-38:800- 5*39)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_moyen(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    
    %% Chauffage 7
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_moyen = mean(sol4_deg(:,800- 6*39-38:800- 6*39)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_moyen(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 8
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_moyen = mean(sol4_deg(:,800- 7*39-38:800- 7*39)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_moyen(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %%
    save Data_Moyen.mat        ...
         Position_centre_moyen ...
         Dates_centre_moyen    ...
         Temperature_centre_moyen ...
         tf_moyen ... 
         t4_moyen  sol4_moyen h4_moyen ...
         sol3_moyen  sol3_moyen_mean x3
    %%
    t4_moyen    = t4;
    % Imprimer figure
    csvwrite('Time_LH01.csv',transpose(t4_tot))
    csvwrite('Temperature_LH01.csv',transpose(Tc))
    figure
    plot(t4_tot, vol*100)
        hXLabel = xlabel('Time [s]');
        hYLabel = ylabel('Cristalisation [\%]');
        set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca,'FontSize',10);
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
        set(gca,'Layer','Top');
%         set(gca, 'Position', [0.3, 0.3, 0.65, 0.65])
        set(gcf, 'PaperUnits', 'centimeters');
        x_width=4.45*2; y_width=4*2;  
        set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
        print('-painters','-dpdf','Cristalisation_LH_01')
    figure
    plot(t4_tot, Tc)
        hXLabel = xlabel('Time [s]');
        hYLabel = ylabel('Temprature [C]');
        set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca,'FontSize',10);
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
        set(gca,'Layer','Top');
%         set(gca, 'Position', [0.3, 0.3, 0.65, 0.65])
        set(gcf, 'PaperUnits', 'centimeters');
        x_width=4.45*2; y_width=4*2;  
        set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
        print('-painters','-dpdf','temp_LH_01')
%     
   
end  
for hide_lent   = [] 

    %% Solveur -- La ou tu ne peux plus jouer

    tf_ = tf1_ + tf2_+ tf3_;

    Temperature_centre_ = [];
    Dates_centre_       = [];
    Position_centre_    = [];

    %% I.- Cylindre 1
        R1 = R1*1e-3; %[m]
        NBptsR =50; NBptsT = 50; % Pour affiner le maillage de l'image
        x1_ = 0:R1/(NBptsR-1):R1; % linspace(0,L,20);
        p1_ = 0:4.125/(NBptsT-1):4.125;
        t1_ = 0:tf1_/(NBptsT-1):tf1_; % [linspace(0,0.02,20), linspace(0.5,1.5,10)];
        m=1; % Geometrie cylindrique
        options = odeset('RelTol',1e-4);
        sol1_ = pdepe(m,@heatpde,@heatic,@heatbc,x1_,t1_,options);
    % Save data
    Temperature_centre_ = [Temperature_centre_; sol1_(:,1)];
    Dates_centre_       = [Dates_centre_        ,t1_];
    Position_centre_    = [Position_centre_     ,p1_];

    %% II.- Cylindres 2
    Profil.xp = x1_;
    Profil.Tp = sol1_(end,:);

        R2 = R2*1e-3; %[m]
        NBptsR =50; NBptsT = 50; % Pour affiner le maillage de l'image
        x2_ = 0:R2/(NBptsR-1):R2; % linspace(0,L,20);
        p2_ = 4.125:0.375/(NBptsT-1): 4.5;
        t2_ = 0:tf2_/(NBptsT-1):tf2_; % [linspace(0,0.02,20), linspace(0.5,1.5,10)];
        m=1; % Geometrie cylindrique
        options = odeset('RelTol',1e-4);
        sol2_ = pdepe(m,@heatpde,@heaticProf,@heatbc,x2_,t2_,options);
    % Save data
    Temperature_centre_ = [Temperature_centre_; sol2_(:,1)];
    Dates_centre_       = [Dates_centre_        ,t2_+tf1_];
    Position_centre_    = [Position_centre_     ,p2_];


    %% III.- Cylindres 3
    Profil.xp = x2_;
    Profil.Tp = sol2_(end,:);

        R3 = R3*1e-3; %[m]
        NBptsR =50; NBptsT = 50; % Pour affiner le maillage de l'image
        x3_ = 0:R3/(NBptsR-1):R3; % linspace(0,L,20);
        p3_ = 4.5:1.5/(NBptsT-1+147): 6; % 147
        t3_ = 0:tf3/(NBptsT-1):tf3_; % [linspace(0,0.02,20), linspace(0.5,1.5,10)];
        m=1; % Geometrie cylindrique
        options = odeset('RelTol',1e-7);
        sol3_ = pdepe(m,@heatpde,@heaticProf,@heatbc,x3_,t3_,options);
    % Save data
    Temperature_centre_ = [Temperature_centre_;sol3_(:,1)];
    Dates_centre_       = [Dates_centre_, t3_+tf1_+tf2_];
    Position_centre_    = [Position_centre_     ,p3_];


    Temperature_centre_lent = Temperature_centre_;
    Dates_centre_lent       = Dates_centre_;
    Position_centre_lent    = Position_centre_;
    
    
    
    %% IV.- Empilement de boudins
    sol3_lent      = sol3_(end,:);
    sol3_lent_mean = mean(sol3_(end,:));
%     sol3_lent_mean = 225+273.15;
    Temp_mean_top = sol3_lent_mean; % Initialisation temperature layer top a moyenne du profil de temp sorant de la buse
    Temp_mean_bot =   60+273.15;            %Ti;   %   "              "          "     bottom a temperature de la piece
    NBptsh = 800;   
    h4 = -h_bot_:(h_bot_+h_top_)/(NBptsh-1):(h_top_); % Discretisation du PB
    h4 = h4*1e-3;
    t4 = 0:0.01:tf_attente;
    m=0; % Geometrie euclidienne
    options = odeset('RelTol',1e-5);
    sol4_ = pdepe(m,@heatpde,@heatic_boudin,@heatbc_boudin,h4,t4,options);
    
    % extraire donnee au temps donne
    t4_lent   = t4;
    sol4_lent = sol4_-273.15;
    h4_lent   = h4; 
    
     %% Analyse cristalisation
    % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    t4 = 0:0.01:tf_attente;  %             \
    t4 = [t4(1),t4(2:end)];  %  <__________/
    m=0; % Geometrie euclidienne
    options = odeset('RelTol',1e-5);
    % Temperature moyenne [K] du boudin a chaque increment de temps: 
    sol4 = pdepe(m,@heatpde,@heatic_boudin,@heatbc_boudin,h4,t4,options);
    sol4_deg = sol4 - 273.15;
    T_mean_lent = mean(sol4_deg(:,791:end)');  
    % Interval de temps entre deux mesure de moyenne: 
    dt = 0.01; % [s]
    crist(1) = 0;
    %%
    cellules_r = [];
    t4_tot     = [];
    vol        = []; 
    Tc         = []; 
    for i =1:length(t4) 
        t4_tot = [t4_tot, t4(i)];
        Tc = [Tc, T_mean_lent(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 2
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_lent = mean(sol4_deg(:,800- 1*10-9:800- 1*10)'); 
    t4_last = t4_tot(end);  
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_lent(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 3
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_lent = mean(sol4_deg(:,800- 2*10-9:800- 2*10)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_lent(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 4
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_lent = mean(sol4_deg(:,800- 3*10-9:800- 3*10)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_lent(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 5
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_lent = mean(sol4_deg(:,800- 4*10-9:800- 4*10)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_lent(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    
    %% Chauffage 6
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_lent = mean(sol4_deg(:,800- 5*10-9:800- 5*10)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_lent(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    
    %% Chauffage 7
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_lent = mean(sol4_deg(:,800- 6*10-9:800- 6*10)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_lent(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %% Chauffage 8
     % On re resout sans le trick __________  (i.e. le petit step de 0.001)
    T_mean_lent = mean(sol4_deg(:,800- 7*10-9:800- 7*10)');
    t4_last = t4_tot(end);    
    %                                    ^            ^
    for i =1:length(t4) 
        t4_tot = [t4_tot , t4(i)+t4_last]; 
        Tc = [Tc, T_mean_lent(i)];
        % nb de cristaux entrain de se dvp par micro m:
        N_current = T2N(Tc(end)); 
        G         = T2G(Tc(end));
        for j = 1:round(N_current)
            try 
                isempty(cellules_r(j));
            catch
                cellules_r(j) = 0;
            end
            cellules_r(j) = cellules_r(j) + G*dt;
        end
        vol = [vol, 0];
        for j = 1:length(cellules_r)
            vol(end) = vol(end) + cellules_r(j)^2*pi;
        end
    end
    %%
    
    save Data_Lent.mat        ...
         Position_centre_lent ...
         Dates_centre_lent    ...
         Temperature_centre_lent ...
         tf_ ... 
         t4_lent  sol4_lent  h4_lent...
         sol3_lent  sol3_lent_mean
    %%
    t4_lent    = t4;
%     % Imprimer figure
    figure
    csvwrite('Time_LH005.csv',transpose(t4_tot))
    csvwrite('Temperature_LH005.csv',transpose(Tc))
    plot(t4_tot, vol*100)
        hXLabel = xlabel('Time [s]');
        hYLabel = ylabel('Cristalisation [\%]');
        set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca,'FontSize',10);
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
        set(gca,'Layer','Top');
%         set(gca, 'Position', [0.3, 0.3, 0.65, 0.65])
        set(gcf, 'PaperUnits', 'centimeters');
        x_width=4.45*2; y_width=4*2;  
        set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
        print('-painters','-dpdf','Cristalisation_LH_005')
    figure
    plot(t4_tot, Tc)
        hXLabel = xlabel('Time [s]');
        hYLabel = ylabel('Temprature [C]');
        set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca,'FontSize',10);
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
        set(gca,'Layer','Top');
%         set(gca, 'Position', [0.3, 0.3, 0.65, 0.65])
        set(gcf, 'PaperUnits', 'centimeters');
        x_width=4.45*2; y_width=4*2;  
        set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
        print('-painters','-dpdf','temp_LH_005')
    
   
end         
for plot_figure = [] 
    %%
    clear all
    clc
    close all
    
    Ti  = 20  +273.15;  % [deg]:...Temperature piece  <<<<<<<<<<<<<<
    To  = 225 +273.15;  % [deg]:...Temperature buse   <<<<<<<<<<<<<<
    
    load Data_Rapide.mat
    load Data_Moyen.mat
    load Data_Lent.mat
    
%     %%
    figure(1)
    hold on
    plot([0 tf_], [1 1]*(To-273.15), 'k--','Linewidth',1)
    plot([0 tf_], [1 1]*(60), 'c-','Linewidth',1)
    plot([0 tf_], [1 1]*(150),'m-','Linewidth',1)
    
    plot(Dates_centre_lent,Temperature_centre_lent-273.15, 'r','Linewidth',2)
    plot(Dates_centre_moyen,Temperature_centre_moyen-273.15, 'b','Linewidth',2)
    plot(Dates_centre_rapide,Temperature_centre_rapide-273.15, 'k','Linewidth',2)
    ylim([250 550]-273.15)
    xlim([0 tf_])
        hXLabel = xlabel('Time [s]');
        hYLabel = ylabel('Temperature [${}^{\circ}$C]');
        set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca,'FontSize',10);
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
        set(gca,'Layer','Top');
%         set(gca, 'Position', [0.3, 0.3, 0.65, 0.65])
        set(gcf, 'PaperUnits', 'centimeters');
        x_width=4.45*2; y_width=4*2;  
        set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
        print('-painters','-dpdf','Time_vs_Temperature_Dans_Buse')


    %%
    figure(2)
    hold on
    plot(Position_centre_lent,Temperature_centre_lent-273.15,'r-o','MarkerIndices',1:10:length(Temperature_centre_lent),'Linewidth',2)
    plot(Position_centre_moyen,Temperature_centre_moyen-273.15, 'b-d','MarkerIndices',1:10:length(Temperature_centre_moyen),'Linewidth',2)
    plot(Position_centre_rapide,Temperature_centre_rapide-273.15, 'k-p','MarkerIndices',1:10:length(Temperature_centre_rapide),'Linewidth',2)
    plot([0 6],[1 1]*225,'k--')
    ylim([10 250])
    xlim([0 6])
        hXLabel = xlabel('Position Z [mm]');
        hYLabel = ylabel('Temperature [${}^{\circ}$C]');
        %legend('LH 0.05 mm','LH 0.1 mm', 'LH 0.2 mm','Location','NorthWest','FontSize',10,'Interpreter','latex')
        set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca,'FontSize',10);
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
        set(gca,'Layer','Top');
%         set(gca, 'Position', [0.3, 0.3, 0.65, 0.65])
        set(gcf, 'PaperUnits', 'centimeters');
        x_width=4.45*2; y_width=4*2;  
        set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
        print('-painters','-dpdf','Position_vs_Temperature_Dans_Buse')
     %%  
     
     % Petit recadrage des donnees (convertion et rearangement de matrice)
     
        h4_moyen  = h4_moyen*1e3;
        h4_moyen  = h4_moyen';
        h4_rapide = h4_rapide*1e3;
        h4_rapide = h4_rapide';
        h4_lent   = h4_lent*1e3;
        h4_lent   = h4_lent';

     %%  
        figure(3)
    subplot(1,6,[1 2])
    pcolor(t4_rapide,h4_rapide,sol4_rapide')
%     colorbar
    caxis([70 110]); % Controller le degrade de couleur caxis([Tmin,Tmax])
    hold on
    shading interp; 
    H4   = repmat(h4_rapide,1,size(sol4_rapide',2));
    T4   = repmat(t4_rapide,size(sol4_rapide',1),1);
    [C,h] =  contour(T4,H4,sol4_rapide',[70:10:130],'k');
    clabel(C,h,'FontSize',9,'Color','k','LabelSpacing',1000)
    plot([0, t4_lent(end)], [0 0], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.2*[1 1], 'k--','LineWidth',2)
    ylim([-0.2 0.2])
    xlim([0 4])
        xticks([0 1 2 3 4])
        hXLabel = xlabel('Time [s]');
        hYLabel = ylabel('Height [mm]');
        
        set(gca,'Layer','Top');
        set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca,'FontSize',10);
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
        
    subplot(1,6,[3 4])
    pcolor(t4_moyen,h4_moyen,sol4_moyen')
    caxis([70 110]); % Controller le degrade de couleur caxis([Tmin,Tmax])
    hold on
    shading interp; 
    H4   = repmat(h4_moyen,1,size(sol4_moyen',2));
    T4   = repmat(t4_moyen,size(sol4_moyen',1),1);
    [C,h] =  contour(T4,H4,sol4_moyen',[70:10:130],'k');
    clabel(C,h,'FontSize',9,'Color','k','LabelSpacing',1000)
    ylim([-0.2 0.2])
    xlim([0 4])
    plot([0, t4_lent(end)], [0 0], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.1*[1 1], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.1*2*[1 1], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.1*3*[1 1], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.1*4*[1 1], 'k--','LineWidth',2)
        xticks([0 1 2 3 4])
        hXLabel = xlabel('Time [s]');
        set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca ,'yticklabel',[])
        set(gca,'Layer','Top');
        set(gca,'FontSize',10);
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
    
    subplot(1,6,[5 6])
    pcolor(t4_lent,h4_lent,sol4_lent')
    caxis([70 110]); % Controller le degrade de couleur caxis([Tmin,Tmax])
    hold on
    shading interp; 
%     H4   = repmat(h4_lent,1,size(sol4_lent',2));
%     T4   = repmat(t4_lent,size(sol4_lent',1),1);
%     [C,h] =  contour(T4,H4,sol4_lent','k');
%     clabel(C,h,'FontSize',9,'Color','k','LabelSpacing',1000)
    ylim([-0.2 0.2])
    xlim([0 4])
    plot([0, t4_lent(end)], [0 0], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.05*[1 1], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.05*2*[1 1], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.05*3*[1 1], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.05*4*[1 1], 'k--','LineWidth',2)
        xticks([0 1 2 3 4])
        hXLabel = xlabel('Time [s]');
        set([hXLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca,'FontSize',10);
        set(gca ,'yticklabel',[])
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
        set(gca,'Layer','Top');
    
        
        hp4 = get(subplot(1,6,[5 6]),'Position');
        h = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.025   hp4(4)]);  
        set(gcf, 'PaperUnits', 'centimeters');
        x_width=13*3; y_width=4*3;  
        set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
        print('-painters','-r620','-djpeg','Position_vs_Temperature_Dans_Layers_Buse_1')
  
    %%     
    figure
    
    pcolor(t4_lent,h4_lent,sol4_lent')
    caxis([70 110]); % Controller le degrade de couleur caxis([Tmin,Tmax])
    hold on
    shading interp; 
    H4   = repmat(h4_lent,1,size(sol4_lent',2));
    T4   = repmat(t4_lent,size(sol4_lent',1),1);
    [C,h] =  contour(T4,H4,sol4_lent',[70:10:130],'k');
    clabel(C,h,'FontSize',9,'Color','k','LabelSpacing',1000)
    ylim([-0.2 0.2])
    plot([0, t4_lent(end)], [0 0], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.05*[1 1], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.05*2*[1 1], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.05*3*[1 1], 'k--','LineWidth',2)
    plot([0, t4_lent(end)], -0.05*4*[1 1], 'k--','LineWidth',2)
    xlim([0 0.75])
        xticks([0 0.25 0.5 0.75])
        hXLabel = xlabel('Time [s]');
        hYLabel = ylabel('Height [mm]');
        set(gca,'Layer','Top');
        set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca,'FontSize',10);
%         set(gca ,'yticklabel',[])
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
        set(gca,'Layer','Top');
    
        
%         hp4 = get(subplot(1,3,3),'Position');
%         h = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.025   hp4(4)]);  
        set(gcf, 'PaperUnits', 'centimeters');
        x_width=13; y_width=4*3;  
        set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
        print('-painters','-r620','-djpeg','Position_vs_Temperature_Dans_Layer_0_05_Buse_1_ZOOM')
     
     
     
     
     
     
     
     
%      
%      %%%
%     figure(3)
%     pcolor(t4_rapide,h4_rapide,sol4_rapide')
%     colorbar
%     caxis([70 110]); % Controller le degrade de couleur caxis([Tmin,Tmax])
%     hold on
%     shading interp; 
%     H4   = repmat(h4_rapide,1,size(sol4_rapide',2));
%     T4   = repmat(t4_rapide,size(sol4_rapide',1),1);
%     [C,h] =  contour(T4,H4,sol4_rapide','k');
%     clabel(C,h,'FontSize',9,'Color','k','LabelSpacing',1000)
%     ylim([-0.2 0.2])
%         hXLabel = xlabel('Time [s]');
%         hYLabel = ylabel('Height [mm]');
%         set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
%         set(gca,'FontSize',10);
%         set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
%         set(gca,'Layer','Top');
%         set(gcf, 'PaperUnits', 'centimeters');
%         x_width=13; y_width=4*2;  
%         set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
%         print('-painters','-dpdf','Position_vs_Temperature_Dans_Layer_0_2_Buse_1')
%      %%  
%     figure(4)
%     pcolor(t4_moyen,h4_moyen,sol4_moyen')
%     colorbar
%     caxis([70 110]); % Controller le degrade de couleur caxis([Tmin,Tmax])
%     hold on
%     shading interp; 
%     H4   = repmat(h4_moyen,1,size(sol4_moyen',2));
%     T4   = repmat(t4_moyen,size(sol4_moyen',1),1);
%     [C,h] =  contour(T4,H4,sol4_moyen','k');
%     clabel(C,h,'FontSize',9,'Color','k','LabelSpacing',1000)
%     ylim([-0.1 0.1])
%         hXLabel = xlabel('height [mm]');
%         hYLabel = ylabel('Temperature [${}^{\circ}$C]');
%         set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
%         set(gca,'FontSize',10);
%         set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
%         set(gca,'Layer','Top');
%         set(gcf, 'PaperUnits', 'centimeters');
%         x_width=13; y_width=4*2;  
%         set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
%         print('-painters','-dpdf','Position_vs_Temperature_Dans_Layer_0_1_Buse_1')
%     
%      %%  
%     figure(5)
%     pcolor(t4_lent,h4_lent,sol4_lent')
%     colorbar
%     caxis([70 110]); % Controller le degrade de couleur caxis([Tmin,Tmax])
%     hold on
%     shading interp; 
%     H4   = repmat(h4_lent,1,size(sol4_lent',2));
%     T4   = repmat(t4_lent,size(sol4_lent',1),1);
%     [C,h] =  contour(T4,H4,sol4_lent','k');
%     clabel(C,h,'FontSize',9,'Color','k','LabelSpacing',1000)
%     ylim([-0.05 0.05])
%         hXLabel = xlabel('height [mm]');
%         hYLabel = ylabel('Temperature [${}^{\circ}$C]');
%         set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
%         set(gca,'FontSize',10);
%         set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
%         set(gca,'Layer','Top');
%         set(gcf, 'PaperUnits', 'centimeters');
%         x_width=13; y_width=4*2;  
%         set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
%         print('-painters','-dpdf','Position_vs_Temperature_Dans_Layer_0_05_Buse_1')
%     
        %%
    figure
    hold on
%          sol3_lent  sol3_lent_mean
    plot(x3*1e3, sol3_lent-273.15,'r-o','MarkerIndices',1:10:length(sol3_lent),'Linewidth',2); 
    plot(x3*1e3, sol3_moyen-273.15,'b-d','MarkerIndices',1:10:length(sol3_moyen),'Linewidth',2); 
    plot(x3*1e3, sol3_rapide-273.15,'k-p','MarkerIndices',1:10:length(sol3_rapide),'Linewidth',2); 
    plot(-x3*1e3, sol3_lent-273.15,'r-o','MarkerIndices',1:10:length(sol3_lent),'Linewidth',2); 
    plot(-x3*1e3, sol3_moyen-273.15,'b-d','MarkerIndices',1:10:length(sol3_moyen),'Linewidth',2); 
    plot(-x3*1e3, sol3_rapide-273.15,'k-p','MarkerIndices',1:10:length(sol3_rapide),'Linewidth',2); 
        %legend('LH 0.05 mm','LH 0.1 mm', 'LH 0.2 mm','Location','NorthWest','FontSize',10,'Interpreter','latex')
        hXLabel = xlabel('Radius [mm]');
        ylim([140 230])
        hYLabel = ylabel('Temperature [${}^{\circ}$C]');
        set([hXLabel, hYLabel],'FontSize', 10 ,'Interpreter','LaTex');
        set(gca,'FontSize',10);
        set(gca,'FontName','Helvetica','TickLength',[.02 .02],'LineWidth',1);
        set(gca,'Layer','Top');
%         set(gca, 'Position', [0.3, 0.3, 0.65, 0.65])
        set(gcf, 'PaperUnits', 'centimeters');
        x_width=4.45*2; y_width=4*2;  
        set(gcf, 'PaperPosition', [0 0 x_width y_width],'papersize',[x_width y_width]);
        print('-painters','-dpdf','Position_vs_Temperature_Fin_de_buse')
        
end   


%%
function u0 = heatic(x)
    global Ti 
    u0 = Ti;  % a t0: Temperature de 20deg uniforme de -R a R
end
function u0 = heatic_boudin(h)
    global Temp_mean_top Temp_mean_bot
    if      h < 0
        u0 = Temp_mean_bot;
    elseif  h == 0
        u0 = (Temp_mean_bot+Temp_mean_top)/2;
    else
        u0 = Temp_mean_top;
    end
end
function u0 = heaticProf(x)
    global Profil
    xp = Profil.xp;
    Tp = Profil.Tp;
    u0 = interp1(xp,Tp,x);
end
function [c,f,s] = heatpde(x,t,u,dudx)
    global a
    c = 1/a;  % coefs de l'equation differentielle
    f = dudx; %
    s = 0;    % 
end
function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t)
    global To 
    pl = ul;
    ql = 0;       % Pas d'energie qui arrive de la gauche 
    pr = ur - To;
    qr = 0;       % Pas d'energie qui arrive de la droite
end
function [pl,ql,pr,qr] = heatbc_boudin(xl,ul,xr,ur,t)
    global  Ti 
    pl = 0;
    ql = 1;       % Pas d'energie qui arrive de la gauche 
    pr = 0;
    qr = 1;       % Pas d'energie qui arrive de la droite
end
