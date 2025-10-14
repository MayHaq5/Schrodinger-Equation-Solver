clc;
clear all;
close all;


% Defining necessary values
nm = 10^-9 ;
J = 1;
sec = 1;
eV = 1.60218*(10^-19) ;
mass = 9.10938*(10^-31) ;
h_cut = (6.62607*(10^-34))/(2*pi) *(J*sec);
eV = 1.60218E-19 *(J);
dx = 0.01 *(nm);
V_low = 0.00 *(eV);


% Taking inputs

prompt = 'Potential Well Type(Infinite/ Finite):';
str = input(prompt, 's');
well_width = input('Width of the Potential Well (nm)(between 1-10): ')*(nm);
Range = well_width *4;
no_of_waves = input('How many wavefunctions want to see in well: ') ;


% For Infinite Well

if strcmp(str,'Infinite')

    % Creating matrix for position

    x = ( 0:dx:well_width );
    a = ((-Range/2)+(well_width):dx:(Range/2));


    % creating the Hamiltonian matrix

    H = eye(length(x),length(x))*(h_cut^2/(mass*dx^2));

    for i = 1: length(x)-1
        H(i,i+1) = 1*(-h_cut^2/(2*mass*dx^2));
        H(i+1,i) = 1*(-h_cut^2/(2*mass*dx^2));
    end


    % Determining the wavefunction and Energy levels

    [shi,En] = eig(H);
    E = diag(En);


    % Plotting wave function and energy
    xt = [];
    yt = [];

    for plot_no = 1:no_of_waves
        shi_plot_scale = ((E(1,1)-V_low))/(max(shi(:,plot_no))-min(shi(:,plot_no)));
        shi_plot = shi(:,plot_no)*shi_plot_scale + E(plot_no);
        plot(x/nm, shi_plot/eV,'r','Linewidth', 1.5); hold on;
        plot(x/nm, E(plot_no)*ones(size(x))/eV,'b','Linewidth', 1.5) ; hold on;

        % printing Energy state no in plot
        xt(plot_no) = well_width/nm;
        yt(plot_no) = E(plot_no)/eV;
        Energy_state_no = cell(1,plot_no);
        string1 = '\leftarrow n = ';
        string2 = num2str(plot_no);
        string3 = append(string1,string2);
        Energy_state_no (1,plot_no)  = {string3};
        text(xt,yt,Energy_state_no);



    end
    legend('Wavefunction','Energy level')

    % printing Energy Values

        fprintf('\n<strong>Output:\n corresponding Energy Values (eV):</strong>\n');
        for plot_no = 1:no_of_waves
            fprintf('E%d = %d\n',plot_no,E(plot_no));
        end


    % Creating the matrix for potential

    V_high = E(no_of_waves)*1.5;
    V = zeros(length(a),length(a));
    for i = 1:length(a)-1
        V(i,i)= V_high;
    end

    % plotting the well

    V_plot = diag(V)';
    V_plot(a>=0  & a<=well_width) = V_low;
    plot(a/nm, V_plot/eV, 'k','Linewidth', 2,'DisplayName','Potential Well') ; hold on;
    xlim([min(a/nm),max(a/nm)]) , ylim([min(V_plot/eV),max(V_plot/eV)])
    xlabel('x [nm]','FontWeight','bold')
    ylabel('Energy  (eV)','FontWeight','bold')
    title('Distance Vs Energy and Wavefunction Plot' )
    legend



    % For Finite Well

elseif strcmp(str,'Finite')

    V_high = input('Value of potential outside the well (between 0.10-0.50): ') *(eV);

    % Creating matrix for position
    x = ((-Range/2):dx:(Range/2));

    % Creating the matrix for potential
    V = zeros(length(x),length(x));
    for i = 1:length(x)
        V(i,i)= V_high;
    end

    % plotting the well
    V_plot = diag(V)';
    V_plot(x>=0  & x<=well_width/2) = V_low;
    V_plot(x<=0 & x>=-well_width/2) = V_low;
    plot(x/nm, V_plot/eV, 'k','Linewidth', 2) ; hold on;
    xlim([min(x/nm),max(x/nm)]), ylim([min(V_plot/eV),max(1.1*V_plot/eV)])
    xlabel('x [nm]','FontWeight','bold')
    ylabel('Energy  (eV)','FontWeight','bold')
    title('Distance Vs Energy and Wavefunction Plot' )
    hold on;

    % creating the Hamiltonian matrix
    H = zeros(length(x),length(x));
    H = diag((h_cut^2/(mass*dx^2)) + V_plot);

    for i = 1: length(x)-1
        H(i,i+1) = 1*(-h_cut^2/(2*mass*dx^2));
        H(i+1,i) = 1*(-h_cut^2/(2*mass*dx^2));
    end

    % Determining the wavefunction and Energy levels
    [shi,En] = eig(H);
    E = diag(En);

    % Plotting wave function and energy
    xt = [];
    yt = [];
    total_energy_level = sum(E > V_low & E <V_high);
    ytmin = (min(E)/eV)/4;
    if (no_of_waves <= total_energy_level)

        for plot_no = 1:no_of_waves

            shi_plot_scale = ((E(1,1)-V_low))/(max(shi(:,plot_no))-min(shi(:,plot_no)));
            shi_plot = shi(:,plot_no)*shi_plot_scale + E(plot_no);
            plot(x/nm, shi_plot/eV,'r','Linewidth', 1.5) ; hold on;
            plot(x/nm, E(plot_no)*ones(size(x))/eV,'b','Linewidth', 1.5) ; hold on;

            % printing Energy state no in plot
            xt(plot_no) = (well_width/2)/nm;
            yt(plot_no) = (E(plot_no)/eV)+ ytmin;
            Energy_state_no = cell(1,plot_no);
            string1 = '\leftarrow n = ';
            string2 = num2str(plot_no);
            string3 = append(string1,string2);
            Energy_state_no (1,plot_no)  = {string3};
            text(xt,yt, Energy_state_no);


        end

         % printing Energy Values

        fprintf('\n<strong>Output:\n Corresponding Energy Values (eV):</strong>\n');
        for plot_no = 1:no_of_waves
            fprintf('E%d = %d\n',plot_no,E(plot_no));
        end


    elseif (no_of_waves > total_energy_level)

        fprintf('<strong>Comment: Only %d wavefunctions exist in this well</strong>\n', total_energy_level);
        no_of_waves = total_energy_level;

        for plot_no = 1:no_of_waves
            shi_plot_scale = ((E(1,1)-V_low))/(max(shi(:,plot_no))-min(shi(:,plot_no)));
            shi_plot = shi(:,plot_no)*shi_plot_scale + E(plot_no);
            plot(x/nm, shi_plot/eV,'r','Linewidth', 1.5) ; hold on;
            plot(x/nm, E(plot_no)*ones(size(x))/eV ,'b','Linewidth', 1.5 ) ; hold on;

            % printing Energy state no in plot
            xt(plot_no) = (well_width/2)/nm;
            yt(plot_no) = (E(plot_no)/eV) + ytmin;
            Energy_state_no = cell(1,plot_no);
            string1 = '\leftarrow n = ';
            string2 = num2str(plot_no);
            string3 = append(string1,string2);
            Energy_state_no (1,plot_no)  = {string3};
            text(xt,yt, Energy_state_no);

        end

         % printing Energy Values

        fprintf('\n<strong>Output:\n Corresponding Energy Values (eV):</strong>\n');
        for plot_no = 1:no_of_waves
            fprintf('E%d = %d\n',plot_no,E(plot_no));
        end

    end

    legend('Potential Well','Wavefunction','Energy level' )

end

