function [out, posterior] = singleChannelFrameNoiseReduce(in, history, prior)
    amp = abs(in);
    phase = angle(in);
    noise = min(abs(history), [], 2);
    amp = max(min(amp, 0.01), amp - 16 * noise);
    out = amp .* exp(1i .* phase);
    posterior = prior;
end
