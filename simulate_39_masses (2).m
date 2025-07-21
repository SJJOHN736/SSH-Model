% By simulating 39 masses (or any odd number), we end with a mass that is grounded with a K2
% spring.  It is ends with a K2 spring, where K2>K1, that show edge modes.
% The harmonic forcing has been changed to this location as well.
function simulate_39_masses
% constants/intial conditons and vectors
    N = 39;
    m = 0.1;
    K1 = 3750;
    K2 = 8750;
    F0 = 1.0;
    omega_hz = 55;  % freq. in Hz
    omega_rad = 2 * pi * omega_hz; %freq. in rad/s
    u0 = zeros(N, 1);
    v0 = zeros(N, 1);
    y0 = zeros(2*N, 1);      
    y0(1:2:end) = u0;        
    y0(2:2:end) = v0;        
   
    tspan = [0 5];
   
    %output matrix
    ode_func = @(t, y) mass_spring_ode(t, y, N, m, K1, K2, F0, omega_rad);
    [t, Y] = ode45(ode_func, tspan, y0);
   
    max_displacement = max(abs(Y(:, 1:2:end)), [], 'all');
    plot_limit = max_displacement * 1.2;
    if plot_limit < 1e-4
        plot_limit = 0.02;
    end
%graphing code
    figure(1);
    h = plot(1:N, Y(1, 1:2:end), 'o-', 'LineWidth', 2);
    axis([1 N -plot_limit plot_limit]);
    grid on;
    xlabel('Mass index');
    ylabel('Displacement (m)');
    title_str = sprintf('Forced Response at %.1f Hz', omega_hz);
    title(title_str);
   
    %animation code
    for frame = 1:5:length(t)
        set(h, 'YData', Y(frame, 1:2:end));
        drawnow;
        pause(0.01);
    end
   
    displacements = Y(:, 1:2:end);
    SS_start_index = find(t > tspan(2)/2, 1);
    max_amplitudes = max(abs(displacements(SS_start_index:end, :)));
    figure(2);
    semilogy(1:N, max_amplitudes, 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('Mass Index');
    ylabel('Maximum Displacement (m)');
    title('Decay of Amplitude Along the Chain');
end

%Equations of motion
function dy = mass_spring_ode(t, y, N, m, K1, K2, F0, omega)
% intializing varibale and dy vector
    dy = zeros(2*N, 1);
%dv of first mass 
    dy(2) = (K2*(y(3) - y(1)) - K1*y(1)) / m;
%dv for middle masses
    for i = 2:N-1
    %check even/odd
        if mod(i, 2) ~= 0
            k_L = K1;
            k_R = K2;
        else
            k_L = K2;
            k_R = K1;
        end
        dy(2*i) = ((-(k_L + k_R)*y((2*i)-1)) + (k_L*y((2*i)-3)) + (k_R*y((2*i)+1))) / m;
    end
%check if last mass is even/odd
    if mod(N-1, 2) ~= 0
        k_L_of_N = K2;
        k_R_of_N = K1;
    else
        k_L_of_N = K1;
        k_R_of_N = K2;
    end
   
    %last mass
    F_shaker = F0 * sin(omega * t);
    dy(2*N) = (-(k_L_of_N + k_R_of_N)*y((2*N)-1) + (k_L_of_N*y((2*N)-3)) + F_shaker) / m;
   
    % odd of dy vector is the evens of y vector
    dy(1:2:end) = y(2:2:end);
end
