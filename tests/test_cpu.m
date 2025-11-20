% test_cpu.m
% Simple CPU benchmark for MATLAB
% - Matrix multiplication benchmark (approx GFLOPS)
% - Simple integer-loop throughput test
% Usage:
%   matlab -batch "run('tests/test_cpu.m')"
% Optional override of matrix size:
%   export MATLAB_CPU_BENCH_N=800
%   matlab -batch "run('tests/test_cpu.m')"

% Configuration: default matrix size
n = 800; % default matrix size (800x800 is safe on most machines)
envN = getenv('MATLAB_CPU_BENCH_N');
if ~isempty(envN)
    v = str2double(envN);
    if ~isnan(v) && v>0
        n = round(v);
    end
end

fprintf('\nMATLAB CPU benchmark (lightweight) - %s\n', datestr(now));

% Number of cores (best-effort)
numCores = NaN;
try
    numCores = feature('numcores');
catch
    % feature may not be available in some builds; ignore
end
if ~isnan(numCores)
    fprintf('Detected cores (feature(''numcores'')): %d\n', numCores);
else
    fprintf('Detected cores: unknown\n');
end

% Matrix multiplication benchmark
fprintf('\nMatrix multiplication benchmark: %d x %d\n', n, n);
try
    A = rand(n, n, 'double');
    B = rand(n, n, 'double');
    fprintf('Allocations done. Running A*B... ');
    t = tic;
    C = A * B; %#ok<NASGU>
    elapsed = toc(t);
    flops_est = 2 * (double(n)^3); % approx FLOPs for dense multiply
    gflops = flops_est / (elapsed * 1e9);
    fprintf('done.\nTime: %.4f s | Estimated GFLOPS: %.3f\n', elapsed, gflops);
catch ME
    fprintf('Matrix benchmark failed: %s\n', ME.message);
end

% Integer loop throughput (simple scalar loop)
% Keep iterations modest so this finishes quickly on CI/low-resource boxes
iters = 1e6; % one million iterations
fprintf('\nInteger add loop: %d iterations\n', iters);
try
    s = uint64(0);
    t2 = tic;
    for k = uint64(1):uint64(iters)
        s = s + k; %#ok<NASGU>
    end
    elapsed2 = toc(t2);
    ops_per_sec = double(iters) / elapsed2;
    fprintf('Time: %.4f s | Throughput: %.0f ops/s\n', elapsed2, ops_per_sec);
catch ME
    fprintf('Integer loop failed: %s\n', ME.message);
end

fprintf('\nSummary:\n');
fprintf('  Matrix size: %d x %d\n', n, n);
if ~isnan(numCores)
    fprintf('  Detected cores: %d\n', numCores);
end
fprintf('  Matrix GFLOPS (approx): %s\n', num2str(gflops));
fprintf('  Integer ops/sec (approx): %s\n', num2str(ops_per_sec));
fprintf('\nNote: This is a lightweight, best-effort test. For accurate FLOPS measurements, use a dedicated benchmark on a machine with MATLAB configured.\n\n');
