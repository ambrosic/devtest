% buildfile.m
% Simple build runner for CI: executes the CPU test script and exits with
% status 0 on success or 1 on failure.
%
% Usage (from shell):
%   matlab -batch "run('buildfile.m')"
% or:
%   matlab -batch "buildfile"

% Path to the test script
testScript = fullfile(pwd, 'tests', 'test_cpu.m');

fprintf('\nRunning CPU test via %s\n', testScript);

try
    run(testScript);
    fprintf('\nCPU test completed successfully.\n');
    exit(0);
catch ME
    fprintf('\nCPU test failed: %s\n', ME.message);
    if isfield(ME, 'stack') && ~isempty(ME.stack)
        fprintf('Stack (most recent call first):\n');
        for k = 1:numel(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(k).file, ME.stack(k).line);
        end
    end
    % Ensure non-zero exit so CI can detect failure
    exit(1);
end
