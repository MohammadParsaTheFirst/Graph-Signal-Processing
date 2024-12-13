%% Part 1: Creating the original graph
N = 150;
p = 8 * log(N) / N;
q = log(N) / N;
k = 3;
dim = k;

seed = 20;
rng(seed)
sig = randi([0, k - 1], [N, 1]);
sig = 2 * rescale(sig) - 1;
P = (sig == sig') * p + (sig ~= sig') * q;
A = ones(N) - (rand(N) > P);
A = triu(A);
A = abs(A - A');

G = gsp_graph(A);
G = gsp_create_laplacian(G);
G = gsp_compute_fourier_basis(G);

U23 = G.U;
G.coords = U23(:, 2:max(3, dim));

param = struct;
param.colorbar = 0;
G.plotting.edge_width = 0.1;
figure;
gsp_plot_signal(G, sig, param);
title('SBM Graph Clustering');
%% part 2: generating the filtered smooth signals

T = 10;  
r = 5;   
d_max = max(sum(A, 2)); 
alpha = 1 / (2 * d_max);
X = randn(n, T);

U = G.U;        
lambda = G.e;   
H_lambda = (1 - alpha * lambda).^(r-1);

X_filtered = zeros(size(X));
for t = 1:T
    x_hat = U' * X(:, t);
    y_hat = H_lambda .* x_hat;
    X_filtered(:, t) = U * y_hat;
end

figure();
subplot(1, 2, 1);
plot(X(:, 1), 'b');
title('Original Signal (First Sample)');
xlabel('Node Indices');
ylabel('Signal Values');
grid("minor");

subplot(1, 2, 2);
plot(X_filtered(:, 1), 'r');
title('Filtered Signal (First Sample)');
xlabel('Node Indcies');
ylabel('Signal Values');
grid("minor");

%% Part 2: Generating the filtered smooth signals

T = 1000;
r = 31;
% r = 200
%r = 800;
X_smooth = gen_smooth_filtered(G, T, r);

figure();
gsp_plot_signal(G, X_smooth(:,1), param);
title('Smooth signal on graph');
%% Part 3: Clustering, approximating Laplacian
T_values = round(logspace(1, 3, 10));
r_values = [11, 21, 26, 31];
num_r_values = numel(r_values);

figure('Position', [100, 100, 1400, 600]);

hold on;
dim = k;
for i = 1:num_r_values
    r = r_values(i);
    accuracy = zeros(numel(T_values), 1);
    for j = 1:numel(T_values)
        t = T_values(j);
        [clusters, ~] = run_single_experiment(G, t, r, k, dim);
        predicted_clusters = clusters';
        true_labels = (sig - min(sig) + 1)';
        accuracy(j) = measure_accuracy(true_labels, predicted_clusters);
        fprintf('t: %d, r: %d, Accuracy: %.2f%%\n', t, r, accuracy(j) * 100);
    end
    
    plot(T_values, accuracy, 'LineWidth', 1.5, 'DisplayName', sprintf('r: %d', r));
end

coords = G.L(:, 2:max(dim, 3));

% Perform k-means clustering
[clusters, ~] = kmeans(coords, k);
predicted_clusters = clusters';
true_labels = (sig - min(sig) + 1)';
known_accuracy = measure_accuracy(true_labels, predicted_clusters);

fprintf('Known graph Accuracy: %.2f%%\n', known_accuracy * 100);

yline(known_accuracy, '--', 'Known graph', 'LineWidth', 1.5,...
    'DisplayName', 'Known graph accuracy', 'FontSize', 14);
ylim([0, 1.1]);

xlabel('t');
ylabel('Accuracy');
title('Clustering Accuracy vs. t for different r values', ['dim: ', num2str(dim-1)]);
legend('Location', 'best');
grid on;
hold off;
%% usedd functions
function x_smooth = gen_smooth_filtered(G, t, r)
    % initial signal
    x = randn(G.N, t);
    % gaussian filter
    d_max = max(G.d);
    alpha = 1 / (2 * d_max);
    H = (eye(G.N) - alpha * G.L) ^ (r-1);
    x_smooth = H * x;
end

function [U, E] = eigen(M)
    % Use covarience to estimate eigenvectors
    [eigenvectors,eigenvalues] = eig(M);
    [E,inds] = sort(diag(eigenvalues),'descend');
    eigenvectors=eigenvectors(:,inds);
    
    % Set first component of each eigenvector to be nonnegative
    signs=sign(eigenvectors(1,:));
    signs(signs==0)=1;
    U = eigenvectors*diag(signs);
end

function [clusters, U] = run_single_experiment(G, t, r, k, dim)
    % Generate smooth signal
    X_smooth = gen_smooth_filtered(G, t, r);

    % Part 3.1: Approximating Laplacian
    % here we forget we actually have the graph
    covarience = cov(X_smooth', 1);
    [U, ~] = eigen(covarience);
    
    % Part 3.2: Finding clusters
    % set coordinates as eigenvectors
    coords = U(:, 2:max(dim, 3));

    % Perform k-means clustering
    [clusters, ~] = kmeans(coords, k);
end

function accuracy = measure_accuracy(true_labels, predicted_clusters)
    % Compute confusion matrix
    num_labels = max([max(true_labels), max(predicted_clusters)]);
    conf_matrix = zeros(num_labels);

    for i = 1:length(true_labels)
        conf_matrix(true_labels(i), predicted_clusters(i)) =...
            conf_matrix(true_labels(i), predicted_clusters(i)) + 1;
    end

    % Find the best matching between labels
    [~, match_col] = max(conf_matrix', [], 1);
    matched_labels = match_col(true_labels);

    % Compute accuracy
    accuracy = sum(matched_labels == predicted_clusters) / length(predicted_clusters);
end