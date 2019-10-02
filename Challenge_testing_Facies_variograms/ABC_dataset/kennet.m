function [wz,pz,tf]=kennet(lyr,wvlt,dt,mopt,fs,disp)
% KENNET - Synthetic seismograms for plane wave, normal incidence propagation
% in 1-D layered media using Kennet's invariant imbedding method.
%
%[WZ,PZ,TF]=KENNET(LYR,WVLT,DT,MOPT,FS,DISP)
%
%	PZ:   seismogram at bottom of layers (particle displacement)
%	WZ:   seismogram at the top (particle displacement)
%	LYR:  [velocity, density, thickness] of layers
%	      LYR is a matrix of 3 columns with
%             number of rows = number of layers
%	WVLT: input source wavelet vector
%	      use [] to specify default wavelet (see also SOURCEWVLT)
%	DT:   time sampling interval of wavelet
%       MOPT: =0 for primaries only, =1 for 1st order multiples
%             anything else calculates all reverberations
%	FS:   option for free surface reverberations.
%             FS=0 for no free surface multiples, FS=1 for free surface effect
%	TF:   [frequency, reflectivity, transmissivity]
%	      Transfer function (complex) of the layered medium
%	DISP: =n (n < number of layers) (optional parameter)
%	      display option to show the reflectivity and transmissivity
%	      as it is being calculated recursively from the bottom of
%	      the layer stack to the top. Not shown when DISP is unspecified
%	      = -1 prevents display of any seismogram.
%
%Without any output arguments, KENNET plots the transmitted and
%reflected seismograms

%Written by T. Mukerji

i=sqrt(-1);
if (length(wvlt)==0)
    wvlt=sourcewvlt;
end;
wvlt=[wvlt(:)]';

[nlr,nc]=size(lyr);
ro=lyr(:,2);
v=lyr(:,1);
d=lyr(:,3);
n=length(wvlt);
om=(2*pi/(n*dt))*[0:1:(n/2 -1)];

freq=om/(2*pi);
p0=ifft(wvlt);
w0=(1/(ro(1)*v(1)))*p0;
rdhat=zeros(size(om));
tdhat=ones(size(om));
j=2:nlr;
deno=ro(j).*v(j)+ro(j-1).*v(j-1);
rd=(ro(j-1).*v(j-1)-ro(j).*v(j))./deno;
td=2*sqrt(ro(j).*v(j).*ro(j-1).*v(j-1))./deno;
if fs==1, rd=[-1;rd]; td=[1;td]; else rd=[0;rd]; td=[1;td]; end
ru=-rd; tu=td;
for j=nlr:-1:1
    ed=exp(i*(d(j)/v(j))*om);
    
    if mopt==0
        reverb=ones(size(om));
    elseif mopt==1
        reverb=ones(size(om))+ru(j).*ed.*rdhat.*ed;
    else
        reverb=1./(1-ru(j).*ed.*rdhat.*ed);
    end
    
    rdhat=rd(j)+tu(j).*ed.*rdhat.*ed.*reverb.*td(j);
    tdhat=tdhat.*ed.*reverb.*td(j);
    
    if mopt==1
        %  dx=find(abs(rdhat)>1); rdhat(dx)=ones(size(dx)).*exp(i*phase(rdhat(dx)));
        %  dx=find(abs(tdhat)>1); tdhat(dx)=ones(size(dx)).*exp(i*phase(tdhat(dx)));
        dx=find(real(rdhat)>1); rdhat(dx)=ones(size(dx))+i*imag(rdhat(dx));
        dx=find(imag(rdhat)>1); rdhat(dx)=real(rdhat(dx))+i*ones(size(dx));
        dx=find(real(tdhat)>1); tdhat(dx)=ones(size(dx))+i*imag(tdhat(dx));
        dx=find(imag(tdhat)>1); tdhat(dx)=real(tdhat(dx))+i*ones(size(dx));
    end
    
    if   (nargin>5)
        if disp~=-1
            if (rem(j,disp)==0)
                figure(1),plot(freq, abs(tdhat));
                title(['Transmissivity: layer ',num2str(j)]);
                xlabel('freq. Hz'), ylabel('|T(f)|'); drawnow
                figure(2),plot(freq, abs(rdhat));
                title(['Reflectivity: layer ',num2str(j)]);
                xlabel('freq. Hz'), ylabel('|R(f)|'); drawnow
            end, end, end
    
end

% Changed to just output reflectivity transmission function
% tf=[freq(:), rdhat(:), tdhat(:)];
tf=[rdhat.*p0(1:length(om))];

pz=tdhat.*p0(1:length(om));
wz=rdhat.*p0(1:length(om));


%
%   fltr=hanning(n/2); fltr=[ones(n/4,1);fltr(n/4+1:n/2)]';
%   pz=pz.*fltr; wz=wz.*fltr;
pz=[pz(1),pz(2:length(pz)),0,conj(fliplr(pz(2:length(pz))))];
wz=[0,wz(2:length(wz)),0,conj(fliplr(wz(2:length(wz))))];


pz=real(fft(pz)); wz=real(fft(wz));
%if ( (nargin>5) & (disp~=-1) )
if nargin==5
    figure(1),plot([0:dt:dt*(length(pz)-1)],(pz)),title('transmitted seismogram');
    figure(2),plot([0:dt:dt*(length(wz)-1)],wz),title('reflected seismogram');
end;
