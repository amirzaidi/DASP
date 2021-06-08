function [out, posterior] = singleChannelFrameNoiseReduce(in, prior)
    out = in;
    cut = 100;
    out(1:cut) = 0;
    out((end-cut):end) = 0;
    posterior = prior;
end
