
filePath = 'C:\Users\USER\Downloads\HW3_comp\Q1\Data_City.csv';
city_data = readtable(filePath);
disp(city_data(1:5, :));
%% Q10 

thr = 400;
tav = 10000;
capitals_data = readtable('C:\Users\USER\Downloads\HW3_comp\Q1\Data_City.csv');

num_capitals = height(capitals_data);
weight = zeros(num_capitals, num_capitals);
coord_matrix = [capitals_data.lat, capitals_data.lng];

weight_generator = @(distance, thr, tav) max(0, tav * (1 - distance / thr));

% Compute the weight, laplacian and degree matrix
for i = 1:num_capitals
    for j = (i + 1):num_capitals
        tmp_ = getDistance(coord_matrix(i, 1), coord_matrix(i, 2), coord_matrix(j, 1), coord_matrix(j, 2));
        weight(i, j) = weight_generator(tmp_, thr, tav);
        weight(j, i) = weight(i, j); 
    end
end
D_ = sum(weight, 2);
laplacian = diag(D_) - weight;

% Create the graph and assigning the coords
city_graph = gsp_graph(laplacian);
city_graph.coords = coord_matrix;

figure;
gsp_plot_graph(city_graph);
title('City Graphs of Iran');
%% Q11
signal = capitals_data.Temp;
signal = double(signal);
disp(class(signal(1)));

figure;
gsp_plot_signal(city_graph, signal);
title('Temperature Signal on City Graph');
colorbar;
%% Q12

N = 3;
flag = 0;
flag_0_interp_error = [];

output = MyPyramidAnalysis(signal, city_graph, N, flag);


%%
% ----------------------------------------------------------------------------------

function output_x = MyPyramidSynthesis(An_pyra_outputs, flag)
    N = length(An_pyra_outputs);
    iters = N - 2:-1:0;
    output_x = {};

    for i = iters
        if i == 0
            current_input = An_pyra_outputs{1};
            downsampled_x = An_pyra_outputs{2}.downsampled_x;
            V = An_pyra_outputs{2}.V;
            interpolation_error = An_pyra_outputs{2}.interpolation_error;
        else
            current_input = An_pyra_outputs{i}.reduced_graph;
            downsampled_x = An_pyra_outputs{i + 1}.downsampled_x;
            V = An_pyra_outputs{i + 1}.V;
            interpolation_error = An_pyra_outputs{i + 1}.interpolation_error;
        end

        U = current_input.U;
        Lambda = current_input.e;

        interpolated_x = MySynthesis(downsampled_x, current_input, V, interpolation_error, flag, U, Lambda);

        output_x{i + 1} = interpolated_x;
    end
end

function interpolated_x = MySynthesis(downsampled_x, G, V, interpolation_error, flag, U, Lambda)
    % Perform interpolation
    temp_x = MyInterpolate(downsampled_x, V, G, flag, U, Lambda);
    
    % Include interpolation error if needed
    interpolated_x = temp_x; % + interpolation_error
end

function outputs = MyPyramidAnalysis(x, G, N, flag)
    x_i = x;
    G_i = G;
    outputs = {};

    % Store initial graph
    outputs{1} = G;

    for i = 1:N
        fprintf('%d level in progress ...\n', i);

        % Perform analysis at the current level
        G_i = gsp_compute_fourier_basis(G_i);  % Ensure Fourier basis is computed
        output_i_next = MyAnalysis(G_i, x_i, flag, G_i.U, G_i.e);

        % Update x and G for the next iteration
        x_i = output_i_next.downsampled_x;
        G_i = output_i_next.reduced_graph;

        % Append current output to the results
        outputs{i + 1} = output_i_next;
    end
end

function outputs = MyAnalysis(G, x, flag, U, Lambda)
    outputs = struct();
    
    V = MyVertexSelection(G, U, Lambda);
    V_not = setdiff(1:G.N, V);

    % Reduced graph using the Kron reduction
    outputs.reduced_graph = Graph(MySKReduction(G, V));

    % Display nodes to keep
    disp('Nodes to keep: ');
    disp(V);

    % Downsample the signal
    downsampled_x = MyDS(x, V);

    % Interpolation error
    interpolated_x = MyInterpolate(downsampled_x, V, G, flag, U, Lambda);
    outputs.interpolation_error = x(:) - interpolated_x(:);

    % Store additional outputs
    outputs.V = V;
    outputs.downsampled_x = downsampled_x;
end

function [x_interp] = MyInterpolate(x_down, V, G, flag, U, Lambda)
    if flag==1
        k = floor(G.N/2);
        K = 1:k;
        %[U, ~] = eigs(G.L, k);
        U_k_tilde = U(V,K);
        U_k = U(:,K);
        x_interp = U_k * pinv(U_k_tilde) * x_down';
    elseif flag==0
        x_interp = article_interpolate(x_down', V, G, U, Lambda);
    end
end


function [x_interp] = article_interpolate(x_down, V, G, U, Lambda)
    e = 1;
    n = G.N;
    L_bar = G.L + e * eye(n);
    V_not = setdiff(1:n, V);
    % Ensure x_down is a column vector
    x_down = x_down(:);
    alpha_star = MySKReduction(L_bar, V, V_not) * x_down;

    x_interp = phi_generator(e, V, U, Lambda) * alpha_star;
end

function [V] = MyVertexSelection(G, U, Lambda)
    %L = G.L;
    %[U, Lambda] = eigs(L, 1, 'largestabs');
    [~, max_eig_index] = max(Lambda);
    u_max = U(:, max_eig_index);
    x = u_max > 0;
    V = find(x == 1);
end

function [G_new] = ReduceGraphByNodes(G, V)
    all_nodes = 1:G.N;
    V_not = setdiff(all_nodes, V);
    W = G.W;
    W(V_not,:) = [];
    W(:,V_not) = [];
    G_new = gsp_graph(W);
end

function [KronReducted] = MySKReduction(laplacian, v_1, v_1_c)
    KronReducted_1 = Laplacian_A_B(laplacian, v_1, v_1);

    negative_term = Laplacian_A_B(laplacian, v_1, v_1_c);
    negative_term = negative_term * pinv(Laplacian_A_B(laplacian, v_1_c, v_1_c));
    negative_term = negative_term * Laplacian_A_B(laplacian, v_1_c, v_1);

    KronReducted = KronReducted_1 - negative_term;
end

function [Reducted_L] = Laplacian_A_B(laplacian, A, B)
    Reducted_L = laplacian(A, B);
end

function [phi] = phi_generator(e, v_1, u, lambd)
    delta = zeros(size(u, 1), length(v_1));

    for i = 1:length(v_1)
        delta(:, i) = delta_j(size(u, 1), v_1(i));
    end

    phi = u * green_kernel(lambd, e) * inv(u) * delta;
end

function [d] = delta_j(n, m)
    d = zeros(n, 1);
    d(m) = 1;
end

function [g] = green_kernel(lambd, e)
    g_ = 1 ./ (lambd + e);
    g = diag(g_);
end

function [H] = Myfilter(G, U, Lambda)
    %[U,Lambda] = eig(G.L);
    H_hat = diag(1 + 2* Lambda) ^-1;
    H = U * H_hat * U';
end

function  [yy] = MyDS(x, V)
    yy = x(V);
end


function D = getDistance(lat1, lon1, lat2, lon2)
    rad = pi / 180;
    radius = 6371; %earth radius in kilometers
    D = abs(acos(sin(lat2 * rad) * sin(lat1 * rad) + cos(lat2 * rad) * cos(lat1 * rad) * cos(lon2 * rad - lon1 * rad)) * radius); %result in Kilometers
end


% G = gsp_sensor(n);
% % n = 58;
% % W = zeros(n,n);
% % W(1,2:n) = 1;
% % W(2,n) = 1;
% % for i=2:n-1
% %     W(i+1,i)=1;
% % end
% % W = W + W';
% % G = gsp_graph(W);
% G.coords = rand(n, 2);
% G = gsp_compute_fourier_basis(G);
% U = G.U;
% Lambda = G.e;
% 
% %%
% 
% %[U, Lambda] = eig(G.L);
% [V] = MyVertexSelection(G, U, Lambda);
% disp('Selected nodes for ''laplacian'' option:');
% disp(V);
% G_new = ReduceGraphByNodes(G,V);
% 
% subplot(2,1,1);
% gsp_plot_graph(G);
% hold on;
% plot(G.coords(:, 1), G.coords(:, 2), 'bo');  
% plot(G.coords(V, 1), G.coords(V, 2), 'ro', 'MarkerSize', 10);  
% title('Selected Nodes (Laplacian)');
% hold off;
% 
% subplot(2,1,2);
% G_new.coords = G. coords;
% gsp_plot_graph(G_new);
% title('The Graph after Downsampling')
% hold off;
% %%
% x = 1:n;
% [x_down] = MyDS(x, V);
% % Interpolate the signal using method from article (flag = 0)
% x_interp1 = MyInterpolate(x_down, V, G, 0, U, Lambda);
% disp('Interpolated Signal (Article Method):');
% disp(x_interp1);
% 
% % Interpolate the signal using method from class (flag = 1)
% x_interp2 = MyInterpolate(x_down, V, G, 1, U, Lambda);
% disp('Interpolated Signal (Class Method):');
% disp(x_interp2);