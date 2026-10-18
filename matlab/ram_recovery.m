function eta_r = ram_recovery(M0)
%RAM_RECOVERY  MIL-E-5008B inlet shock recovery (Mattingly eq. 6.6).
%   eta_r = 1                          for M0 <= 1
%   eta_r = 1 - 0.075 (M0 - 1)^1.35    for 1 < M0 < 5
if M0 <= 1
    eta_r = 1;
else
    eta_r = 1 - 0.075*(M0 - 1)^1.35;
end
end
