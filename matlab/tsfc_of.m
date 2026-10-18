function S = tsfc_of(in, alpha)
%TSFC_OF  TSFC of the separate-flow turbofan as a function of bypass ratio
%(helper for fminbnd in hp1_carpet.m).
in.alpha = alpha;
r = separate_flow_turbofan(in);
S = r.S;
end
