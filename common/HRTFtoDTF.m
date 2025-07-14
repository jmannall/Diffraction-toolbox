function [dtf,ctf]=HRTFtoDTF(hrtf, fs)

%% Check Input
f1 = 0;
f2 = fs/2;

%% Settings
N = size(hrtf,1);
nfft = N; %2^nextpow2(N);

%% Frequency bounds
df = fs/nfft;
f = 0:df:fs-df;
idx = f >= f1 & f <= f2;
idx(nfft/2+2:end) = fliplr(idx(2:nfft/2));

%% CTF calculation

% get magnitude response
hrtff=fft(hrtf,nfft);
ctffavg=mean(log(abs(hrtff)+eps),2);


% Force minimum phase 
ctfflog=mean(log(abs(hrtff)+eps),2);
ctfcep = ifft(ctfflog,nfft);
ctfcep(nfft/2+2:nfft,:) = 0;    % flip acausal part to causal part or simply multiply
ctfcep(2:nfft/2,:) = 2*ctfcep(2:nfft/2,:);    % causal part by 2 (due to symmetry)
ctfflog = fft(ctfcep,nfft);
ctfp = exp(ctfflog);

ctff = exp(ctffavg) .*exp(1j*angle(ctfp));

% get IR:
ctfmtx = real(ifft(ctff, nfft));

%% DTF calculation
dtff = hrtff;
dtff(idx,:) = hrtff(idx,:)./ctff(idx,:);
dtfmtx = ifft(dtff,nfft);

%% Attenuate to avoid clipping
ctfmtx = ctfmtx / 10;
dtfmtx = dtfmtx / 10;

%% Output Objects
dtf = dtfmtx;
ctf = ctfmtx;

function f=myifftreal(c,N) % thanks goto the LTFAT <http://ltfat.sf.net>
if rem(N,2)==0
  f=[c; flipud(conj(c(2:end-1,:)))];
else
  f=[c; flipud(conj(c(2:end,:)))];
end;
f=real(ifft(f,N,1));
