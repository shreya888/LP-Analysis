[s,fs_org]  =   audioread('haaladeebaabbas.wav'); s=s(:,1);
sig =   s(5600:6800);
sig   =   resample(sig,8000,fs_org);
fs    =   8000;
p=12;
[ lpc, lpres, pitch, formants ] = lpAnalysis( sig, p, fs );
