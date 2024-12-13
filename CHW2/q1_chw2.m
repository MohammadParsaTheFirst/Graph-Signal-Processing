%% Part 1: Creating graph
% Define graph
W = [0 1 1 0 0 0 0 1;
      0 0 1 1 1 0 0 1;
      0 0 0 1 0 0 0 0;
      0 0 0 0 1 1 0 1;
      0 0 0 0 0 1 1 1;
      0 0 0 0 0 0 1 0;
      0 0 0 0 0 0 0 1;
      0 0 0 0 0 0 0 0];
W = W + W';
figure('Position', [100, 100, 800, 400]);
G = gsp_graph(W);

N = G.N; 
theta = linspace(0, 2*pi, N+1); 
theta(end) = [];
radius = 1;
G.coords = [radius * cos(theta)', radius * sin(theta)'];

G = gsp_create_laplacian(G);
G.plotting.vertex_size = 20;
% Plot graph
gsp_plot_graph(G);
title('Graph');
hold on;
for i = 1:G.N
    text(G.coords(i, 1)-0.03, G.coords(i, 2)-0.03, ...
        sprintf('%d', i), 'FontSize', 12, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Color', 'red');
end
hold off;

title('Graph $G$', Interpreter='latex');
%% part 2
G = gsp_compute_fourier_basis(G);
U12 = G.U(:, 1:2);

x = 2*U12(:,1) + U12(:,2);
x_n = awgn(x, 10);
noise = x_n - x;
signal_power = sum(x.^2) / length(x);
noise_power = sum(noise.^2) / length(noise);
snr_value = 10 * log10(signal_power / noise_power);
disp(['SNR: ', num2str(snr_value), ' dB']);

figure('Position', [100, 100, 600, 500]);

subplot(3,1,1);
gsp_plot_signal(G, x);
title('signal $\mathbf{x} = 2\mathbf{u}_1 + \mathbf{u}_2$ on graph $G$', Interpreter='latex');
hold on;
for i = 1:G.N 
    text(G.coords(i, 1)-0.03, G.coords(i, 2)-0.03, sprintf('%d', i), 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red');
end 
hold off;

subplot(3,1,2);
gsp_plot_signal(G, x_n);
title('the noisy signal $x_n = x + n$ with SNR $= 10dB$', Interpreter='latex');
hold on;
for i = 1:G.N 
    text(G.coords(i, 1)-0.03, G.coords(i, 2)-0.03, sprintf('%d', i), 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red');
end 
hold off;

subplot(3,1,3);
gsp_plot_signal(G, x-x_n);
title('the noisse signal $n = x- x_n$ on graph $G$', Interpreter='latex');
hold on;
for i = 1:G.N 
    text(G.coords(i, 1)-0.03, G.coords(i, 2)-0.03, sprintf('%d', i), 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red');
end 
hold off;
%% part 2

figure;
subplot(3, 1, 1);
stem(x, 'o-', 'LineWidth', 2);
title('Original Graph Signal');
xlabel('Node Index');
ylabel('Signal Value');
ylim([-.5, 1.5]);
grid('minor');

subplot(3, 1, 2);
stem(x_n, 'o-', 'LineWidth', 2);
title('Noisy Graph Signal');
xlabel('Node Index');
ylabel('Signal Value');
ylim([-.5, 1.5]);
grid('minor');

subplot(3, 1, 3);
stem(x -x_n, 'o-', 'LineWidth', 2);
title('Noisy Signal');
xlabel('Node Index');
ylabel('Signal Value');
ylim([-.5, 1.5]);
grid('minor');
%% part 3
for i= 1:G.N
    d(i) = sum(W(i,:));
end

D = diag(d);
L_N = D^(-1/2) * (D - W) * D^(-1/2);

Gn = G;
Gn = gsp_create_laplacian(Gn, 'normalized');

[U_lap, e_lap] = eig(L_N);
Lambda1 = diag(e_lap);
Gn.L = L_N;
Gn = gsp_compute_fourier_basis(Gn);

W_N = D^(-1/2) * W * D^(-1/2);
H = gsp_graph(W_N);
H = gsp_compute_fourier_basis(H);
[U_wn, e_wn] = eig(W_N);
Lambda2 = diag(e_wn);

figure();
subplot(2,1,1);
stem(Lambda1);
title('fourier of graph $G$ with shift operator $L_N$', Interpreter='latex');
grid("minor");
subplot(2,1,2);
stem(Lambda2);
title('fourier of graph $G$ with shift operator $W_N$', Interpreter='latex');
grid("minor");
%% Part 4, 5
Filter_H_L = [1 1 0 0 0 0 0 0];
x_n_hat_L = (U_lap)' *x_n;
y_hat_L = diag(Filter_H_L)* x_n_hat_L;
y_L = (U_lap)* y_hat_L;

% plottings 
subplot(3,1,1);
gsp_plot_signal(G, x);
title('signal $\mathbf{x} = 2\mathbf{u}_1 + \mathbf{u}_2$ on graph $G$', Interpreter='latex');
hold on;
for i = 1:G.N 
    text(G.coords(i, 1)-0.03, G.coords(i, 2)-0.03, sprintf('%d', i), 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red');
end 
hold off;

subplot(3,1,2);
gsp_plot_signal(Gn, x_n);
title('the noisy signal $x_n = x + n$ on graph $G$', Interpreter='latex');
hold on;
for i = 1:G.N 
    text(G.coords(i, 1)-0.03, G.coords(i, 2)-0.03, sprintf('%d', i), 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red');
end 
hold off;

subplot(3,1,3);
gsp_plot_signal(Gn, y_L);
title('the filtered signal $y$ with filter $H$ on signal $x_n$ with shift operator $L_N$', Interpreter='latex');
hold on;
for i = 1:G.N 
    text(G.coords(i, 1)-0.03, G.coords(i, 2)-0.03, sprintf('%d', i), 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red');
end 
hold off;

noise_component = y_L - x; 
filtered_noise_power = sum(noise_component.^2) / length(noise_component); 
snr_filtered_L = 10 * log10(signal_power / filtered_noise_power);
%% part 5

Filter_H_W = [0 0 0 0 0 0 1 1];
x_n_hat_W = (U_wn)' *x_n;
y_hat_W = diag(Filter_H_W)* x_n_hat_W;
y_W = (U_wn)* y_hat_W;

% plottings 
subplot(3,1,1);
gsp_plot_signal(G, x);
title('signal $\mathbf{x} = 2\mathbf{u}_1 + \mathbf{u}_2$ on graph $G$', Interpreter='latex');
hold on;
for i = 1:G.N 
    text(G.coords(i, 1)-0.03, G.coords(i, 2)-0.03, sprintf('%d', i), 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red');
end 
hold off;

subplot(3,1,2);
gsp_plot_signal(G, x_n);
title('the noisy signal $x_n = x + n$ on graph $G$', Interpreter='latex');
hold on;
for i = 1:G.N 
    text(G.coords(i, 1)-0.03, G.coords(i, 2)-0.03, sprintf('%d', i), 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red');
end 
hold off;

subplot(3,1,3);
gsp_plot_signal(G, y_W);
title('the filtered signal $y$ with filter $H$ on signal $x_n$ with shift operator $W_N$', Interpreter='latex');
hold on;
for i = 1:G.N 
    text(G.coords(i, 1)-0.03, G.coords(i, 2)-0.03, sprintf('%d', i), 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red');
end 
hold off;

noise_component = y_W - x; 
filtered_noise_power = sum(noise_component.^2) / length(noise_component); 
snr_filtered_W = 10 * log10(signal_power / filtered_noise_power);
%% part 6

subplot(2,1,1);
stem(Filter_H_L);
title('filter for normalized laplacian shift operator');
subplot(2,1,2);
stem(Filter_H_W);
title('filter for normalized weight shift operator');
%% part 2 cont'd

x_hat_L = (U_lap)' * x;
x_hat_W = (U_wn)' * x;

figure();
subplot(2,2,1);
stem(x_hat_L);
title('fourier of signal graph $x$ on graph $G$ with shift operator $L_N$', Interpreter='latex');
grid("minor");
subplot(2,2,2);
stem(x_hat_W);
title('fourier of signal graph $x$ on graph $G$ with shift operator $W_N$', Interpreter='latex');
grid("minor");
subplot(2,2,3);
stem(x_n_hat_L);
title('fourier of signal graph $x_n$ on graph $G$ with shift operator $L_N$', Interpreter='latex');
grid("minor");
subplot(2,2,4);
stem(x_n_hat_W);
title('fourier of signal graph $x_n$ on graph $G$ with shift operator $W_N$', Interpreter='latex');
grid("minor");

%% part 7

if snr_filtered_W > snr_filtered_L
    disp('The $W_N$ matrix outperforms with respect to $L_N$');
elseif snr_filtered_W == snr_filtered_L
    disp('Both matrices have the same performance');
else
    disp('The $L_N$ matrix outperforms with respect to $W_N$');
end
disp(snr_filtered_L);
disp(snr_filtered_W);

%%
n = 3; 
S123_L = Lambda1 .^ (0:n-1);
S123_Wn =  Lambda2 .^ (0:n-1);
h_L = pinv(S123_L) * Filter_H_L';
h_L
FIR_L = S123_L * h_L;
h_Wn = pinv(S123_Wn) * Filter_H_W';
h_Wn
FIR_Wn = S123_Wn * h_Wn;


% Plotting
figure('Position', [100, 100, 1000, 400]);

% Subplot 1: FIR_L vs Ideal Filter (Laplacian)
subplot(1, 2, 1);
stem(Lambda1, FIR_L, 'LineWidth', 1.5, 'DisplayName', 'FIR_L');
hold on;
stem(Lambda2, Filter_H_L', 'LineWidth', 1.5, 'DisplayName', 'Ideal Filter');
hold off;
xlabel('Eigenvalue');
ylabel('Amplitude');
title('FIR Filter vs Ideal Filter (Laplacian)');
legend('Location', 'northeast');
grid on;

% Subplot 2: FIR_Wn vs Ideal Filter (Norm-Weight)
subplot(1, 2, 2);
stem(Lambda1, FIR_Wn, 'LineWidth', 1.5, 'DisplayName', 'FIR_{Wn}');
hold on;
stem(Lambda2, Filter_H_W', 'LineWidth', 1.5, 'DisplayName', 'Ideal Filter');
hold off;
xlabel('Eigenvalue');
ylabel('Amplitude');
title('FIR Filter vs Ideal Filter (Norm-Weight)');
legend('Location', 'northwest');
grid on;


%% Part 9: Testing FIR filters
x_FIR_L = U_lap * (FIR_L .* x_n_hat_L);
x_FIR_Wn = U_wn * (FIR_Wn .* x_n_hat_W);

% Calculate SNR
snr_FIR_L = db(snr(x_FIR_L, x));
snr_FIR_Wn = db(snr(x_FIR_Wn, x));

% Plotting
figure('Position', [100, 100, 1400, 600]);

subplot(2, 4, 1);
gsp_plot_signal(G, x);
title('Original Signal on Graph')

subplot(2, 4, 2);
gsp_plot_signal(G, x_n);
title('Noisy Signal on Graph', ['SNR: ', num2str(10, 4)])

subplot(2, 4, 3);
gsp_plot_signal(G, x_FIR_L);
title('FIR-Filtered Signal(Laplacian) on Graph', ['SNR: ' num2str(snr_FIR_L, 4)])

subplot(2, 4, 4);
gsp_plot_signal(G, x_FIR_Wn);
title('FIR-Filtered Signal(Norm-Weight) on Graph', ['SNR: ' num2str(snr_FIR_Wn, 4)])

% Plot the signals as time series
y_max = max([x; x_n; x_FIR_Wn; x_FIR_L]) + 0.1;
y_min = min(0, min([x; x_n; x_FIR_Wn; x_FIR_L]) - 0.1);

stem_subplot(x, 2, 4, 'Original Signal', 9, y_min, y_max, 5);
stem_subplot(x_n, 2, 4, 'Noisy Signal', 9, y_min, y_max, 6);
stem_subplot(x_FIR_L, 2, 4, 'FIR-Filtered Signal(Laplacian)', 9, y_min, y_max, 7);
stem_subplot(x_FIR_Wn, 2, 4, 'FIR-Filtered Signal(Norm-Weight)', 9, y_min, y_max, 8);

% Add an overall title for the subplots
sgtitle('Signal, Noisy Signal and FIR-Filtered Signals', 'FontSize', 16);
%% Helper functions
% Function to create a stem plot in a subplot
function stem_subplot(x, n, m, title_str, xlim_val, y_min, y_max, subplot_pos)
    subplot(n, m, subplot_pos);
    stem(x);
    title(title_str);
    xlim([0, xlim_val]);
    ylim([y_min, y_max]);
end

function stem_filter(eigenvalues, filter, title_str)
    stem(eigenvalues, filter);
    hold on
%     plot(eigenvalues, filter);
    title(title_str);
end





