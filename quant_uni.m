function encoded = quant_uni(v,minv,maxv,bits)
% quant_uni performs uniform quantization for x
%   rounds down if exact middle
%   v:    value to quantize
%   minv: min value for quantizer
%   maxv: max value for quantizer
%   bits: bitcount for quantizer

%fprintf('quantizer called with:\n\tv: %f\n\tmin: %f\n\tmax: %f\n\tbits: %d\n',v,minv,maxv,bits)
%% trivial quantizer:
% clamp to min:
if v <= minv, encoded = minv; return; end
% clamp to max:
if v >= maxv, encoded = maxv; return; end
%% nontrivial:
% uniform stepsize
%d = (maxv - minv) / 2^bits; 
b = 1; % loop variable
maxt = maxv; % temp value for max limit
mint = minv; % temp value for min limit
while b <= bits
    d = (maxv - minv)/(b*2); % calculate distance
    if v > d + mint % over halfway between mint -> maxt
        temp = maxt;
        mint = mint + d; % raise minimum limit by half
    else % under or at halfway between mint -> maxt
        temp = mint; 
        maxt = maxt - d; % lower maximum limit by half
    end
    b = b + 1; 
end
% move the value to the middle of the final uniform distribution
if temp == maxt
    temp = maxt - d;
else
    temp = mint + d;
end
%fprintf('quantized to: %f, did %d rounds\n',temp,b - 1)
encoded = temp; 
return;
