function stim = rawStim(label)

stim = @(trial) trial.(label);