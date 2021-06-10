function [Full_U, all_infer_output, all_Error, all_G] = meshsize_compare(sysInfo, infer, all_dir, all_M, big_M)
% Compare the inference result for different space grid sparsity
%
% 1. Genereate data with the LCM of all_M
% 2. Extract the corresponding sparse solution
% 3. Get uniform optimal dimension from the finest solution.
% 4. For each solution, apply the inference algorithm and compare the performance
%% set all M

N = length(all_M);
%% Generate solution with big_M

sysInfo.M = big_M;
sysInfo = get_sysInfo_details(sysInfo);
fprintf(['Generate Data U with large M = ', num2str(big_M), ' for the meshsize compare (slow)...\n'])
Full_U = generate_data(sysInfo, all_dir, 1);


%% extract the sparse solution for all M
all_U = cell(N,1);

for kM = 1:1:N
    M = all_M(kM);
    stepsize = 3000/M;
    all_U{kM} = Full_U(1:stepsize:end, :);
end

%% Normalization
do_normalization = 1;
if do_normalization
    sysInfo.sparse_normalization = 1;
    for kM = 1:N
        int_U = sum(all_U{kM}*(2*sysInfo.L/all_M(kM)));
        all_U{kM} = all_U{kM}./int_U;
    end
end


%% Inference
all_infer_output = cell(N, 1);
all_Error = cell(N, 1);
all_G = cell(N, 1);

% Use the finest data U to find optimal dimension
U = all_U{end};
sysInfo.M = all_M(end);
sysInfo = get_sysInfo_details(sysInfo);
fprintf(['Inference for different meshsize: \nSparse M = ', num2str(all_M(end)),'...'])
G = generate_integration_kernel(sysInfo, U, all_dir, big_M);
all_G{end} = G;



candidate_N = 30;
[all_infer_output{end}, infer] = inference_kernel_optimal_dimension(U, sysInfo, infer, G, candidate_N, all_dir, 0, big_M);
all_Error{end} = evaluation_error(sysInfo, all_infer_output{end}, G);

% Use the updated parameters for the rest inference
for ind = 1:N-1
    M = all_M(end-ind);
    sysInfo.M = M;
    sysInfo = get_sysInfo_details(sysInfo);
    fprintf(['\nSparse M = ', num2str(M),'...'])
    
    U = all_U{end-ind};
    G = generate_integration_kernel(sysInfo, U, all_dir, big_M);
    
    all_G{end-ind} = G;
    all_infer_output{end-ind} = generate_inference_result(infer, sysInfo, U, G, 0, all_dir, big_M);
    all_Error{end-ind} = evaluation_error(sysInfo, all_infer_output{end-ind}, G);
end

fprintf('\n')
end




