function wrightOmegaq_test
%WRIGHTOMEGAQ_TEST  Error tests for wrightOmegaq function.
%   WRIGHTOMEGAQ_TEST generates a series of plots that show the numerical error
%   of the WRIGHTOMEGAQ function over various ranges along the real axis and the
%   complex plane. Results are compared to Matlab's built-in WRIGHTOMEGA
%   function where possible.
%   
%   See also: WRIGHTOMEGAQ, WRIGHTOMEGAQINV, WRIGHTOMEGA

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-14-12
%   Revision: 1.0, 3-12-13


% Plot function over complex plane
v = linspace(-6*pi,6*pi,1e2);
[x,y] = meshgrid(v);
z = x+y*1i;
w = wrightOmegaq(z);

figure('Renderer','zbuffer')
h=surf(v,v,abs(w));
set(h,'EdgeColor','none','FaceColor','interp','FaceLighting','phong');
light('Position',[20 20 4]);
light('Position',[-20 -20 4]);
material dull
hold on
contour(v,v,abs(w),12)
colormap(hot(256))
axis square
xlabel('X')
ylabel('Y')
zlabel('\omega(x+xi)')
title('Wright \omega function')

% set random seed for noise used to jitter inputs
s = RandStream('mt19937ar','Seed',0);

% Real axis

% -6*pi to 6*pi
x = linspace(-6*pi,6*pi,1e3);
d = v(2)-v(1);
x = sort([x x+d*randn(s,size(x))]);
wq = wrightOmegaq(x);
eaq = x-wrightOmegaqInv(wq);
w = wrightOmega(x);
ea = x-wrightOmegaqInv(w);

figure
subplot(211)
plot(x,wq,'b',x,w,'r:')
ylabel('\omega(X)')
title('\omega(X)')
legend('wrightOmegaq','wrightOmega','Location','North');
legend(gca,'boxoff')

subplot(212)
h = plot(x,abs(ea),'r:',x,abs(eaq),'b',[x(1) x(end)],[0 0]+eps,'k--',...
    x,eps(x),'k');
axis([-20 20 -eps 1.1*max(max(abs(eaq)),max(abs(ea)))+eps])
xlabel('X')
ylabel('|X-\omega-log(\omega)|')
title('Absolute Error')
legend(h(3:4),'EPS','EPS(X)','Location','North');
legend(gca,'boxoff')


% -715 to 715
x = linspace(-715,715,1e3);
d = v(2)-v(1);
x = sort([x x+d*randn(s,size(x))]);
wq = wrightOmegaq(x);
eaq = x-wrightOmegaqInv(wq);
w = wrightOmega(x);
ea = x-wrightOmegaqInv(w);

figure
subplot(211)
plot(x,wq,'b',x,w,'r:')
ylabel('\omega(X)')
title('\omega(X)')
legend('wrightOmegaq','wrightOmega','Location','North');
legend(gca,'boxoff')

subplot(212)
h = plot(x,abs(ea),'r:',x,abs(eaq),'b',x,eps(ones(size(x))),'k--',x,eps(x),'k');
axis([x(1) x(end) 0 1.1*max(abs(ea))])
xlabel('X')
ylabel('|X-\omega-log(\omega)|')
title('Absolute Error')
legend(h(3:4),'EPS','EPS(X)','Location','North');
legend(gca,'boxoff')


% -750 to 700, absolute error blows up due to exponential of small values
x = linspace(-750,-700,1e3);
d = v(2)-v(1);
x = sort([x x+d*randn(s,size(x))]);
wq = wrightOmegaq(x);
eaq = x-wrightOmegaqInv(wq);
w = wrightOmega(x);
ea = x-wrightOmegaqInv(w);

figure
subplot(211)
plot(x,wq,'b',x,w,'r:')
axis([x(1) x(end) 0 1.1*max(w)])
ylabel('\omega(X)')
title('\omega(X)')
legend('wrightOmegaq','wrightOmega','Location','North');
legend(gca,'boxoff')

subplot(212)
h = plot(x,abs(eaq),'b',x,abs(ea),'r:',x,eps(ones(size(x))),'k--',x,eps(x),'k');
axis([x(1) x(end) 0 1.1*max(abs(ea))])
xlabel('X')
ylabel('|X-\omega-log(\omega)|')
title('Absolute Error')
legend(h(3:4),'EPS','EPS(X)','Location','North');
legend(gca,'boxoff')


% eps(realmin) to realmax/2 logarithmically spaced, large error around 2^30 due
% to log(), wrightOmega() broken for X > about 2^28, returns bad values
x = 2.^(-1074:1023);
xx = log2(x);
w = wrightOmegaq(x);
ea = x-wrightOmegaqInv(w);

figure
subplot(211)
plot(xx,log10(w))
axis([xx(1) xx(end) 0 max(log10(w))])
ylabel('log_1_0(\omega(X))')
title('\omega(X)')

subplot(212)
plot(xx,abs(ea),'b')
axis([xx(1) xx(end) 0 max(abs(ea))])
xlabel('log_2(X)')
ylabel('|X-\omega-log(\omega)|')
title('Absolute Error')


% Lines of Discontinuity

% -20 to -1 (minus a bit because wrightOmega() only has resolution of a bit
% better than 1e-13; returns -1, value at the singularity, over a large region).
% The ramps from wrightOmega() are due to floating point error.
x = linspace(-20,-1-2^-44,1e3);
z = [x+pi*1i;x-pi*1i];
wq = wrightOmegaq(z);
eaq = z-wrightOmegaqInv(wq);
w = wrightOmega(z);
ea = z-wrightOmegaqInv(w);

figure
subplot(211)
h = plot(x,abs(wq),'b',x,abs(w),'r:');
ylabel('|\omega(X\pmi\pi)|')
title('|\omega(Z)| Lines of Discontinuity')
legend(h([1 3]),'wrightOmegaq','wrightOmega','Location','NorthEast');
legend(gca,'boxoff')

subplot(212)
h = plot(x,abs(ea),'r:',x,abs(eaq),'b',[x(1) x(end)],[0 0]+eps,'k--',...
    x,eps(x),'k');
axis([-20 -1 -eps 1.1*max(max(abs(eaq(:))),max(abs(ea(:))))+eps])
xlabel('X')
ylabel('|Z-\omega-log(-\omega)|')
title('Absolute Error')
legend(h(5:6),'EPS','EPS(X)','Location','NorthWest');
legend(gca,'boxoff')


% -2.5 to -1 (minus a bit because wrightOmega() only has resolution of a bit
% better than 1e-13; returns -1, value at the singularity, over a large region).
% The ramps from wrightOmega() are due to floating point error.
x = linspace(-2.5,-1-2^-44,1e3);
z = [x+pi*1i;x-pi*1i];
wq = wrightOmegaq(z);
eaq = z-wrightOmegaqInv(wq);
w = wrightOmega(z);
ea = z-wrightOmegaqInv(w);

figure
subplot(211)
h = plot(x,abs(wq),'b',x,abs(w),'r:');
ylabel('|\omega(X\pmi\pi)|')
title('|\omega(Z)| Lines of Discontinuity')
legend(h([1 3]),'wrightOmegaq','wrightOmega','Location','NorthEast');
legend(gca,'boxoff')

subplot(212)
h = plot(x,abs(ea),'r:',x,abs(eaq),'b',[x(1) x(end)],[0 0]+eps,'k--',...
    x,eps(x),'k');
axis([-2.5 -1 -eps 1.1*max(max(abs(eaq(:))),max(abs(ea(:))))+eps])
xlabel('X')
ylabel('|Z-\omega-log(-\omega)|')
title('Absolute Error')
legend(h(5:6),'EPS','EPS(X)','Location','NorthWest');
legend(gca,'boxoff')



% -715 to 715
x = linspace(-715,715,1e3);
d = v(2)-v(1);
x = sort([x x+d*randn(s,size(x))]);
wq = wrightOmegaq(x);
eaq = x-wrightOmegaqInv(wq);
w = wrightOmega(x);
ea = x-wrightOmegaqInv(w);

figure
subplot(211)
plot(x,wq,'b',x,w,'r:')
ylabel('\omega(X)')
title('\omega(X)')
legend('wrightOmegaq','wrightOmega','Location','North');
legend(gca,'boxoff')

subplot(212)
h = plot(x,abs(ea),'r:',x,abs(eaq),'b',x,eps(ones(size(x))),'k--',x,eps(x),'k');
axis([x(1) x(end) 0 1.1*max(abs(ea))])
xlabel('X')
ylabel('|X-\omega-log(\omega)|')
title('Absolute Error')
legend(h(3:4),'EPS','EPS(X)','Location','North');
legend(gca,'boxoff')


% -750 to 700, absolute error blows up due to exponential of small values
x = linspace(-750,-700,1e3);
d = v(2)-v(1);
x = sort([x x+d*randn(s,size(x))]);
wq = wrightOmegaq(x);
eaq = x-wrightOmegaqInv(wq);
w = wrightOmega(x);
ea = x-wrightOmegaqInv(w);

figure
subplot(211)
plot(x,wq,'b',x,w,'r:')
axis([x(1) x(end) 0 1.1*max(w)])
ylabel('\omega(X)')
title('\omega(X)')
legend('wrightOmegaq','wrightOmega','Location','North');
legend(gca,'boxoff')

subplot(212)
h = plot(x,abs(eaq),'b',x,abs(ea),'r:',x,eps(ones(size(x))),'k--',x,eps(x),'k');
axis([x(1) x(end) 0 1.1*max(abs(ea))])
xlabel('X')
ylabel('|X-\omega-log(\omega)|')
title('Absolute Error')
legend(h(3:4),'EPS','EPS(X)','Location','North');
legend(gca,'boxoff')


% eps(realmin) to realmax/2 logarithmically spaced, large error around 2^30 due
% to log(), wrightOmega() broken for X > about 2^28, returns bad values
x = 2.^(-1074:1023);
xx = log2(x);
w = wrightOmegaq(x);
ea = x-wrightOmegaqInv(w);

figure
subplot(211)
plot(xx,log10(w))
axis([xx(1) xx(end) 0 max(log10(w))])
ylabel('log_1_0(\omega(X))')
title('\omega(X)')

subplot(212)
plot(xx,abs(ea),'b')
axis([xx(1) xx(end) 0 max(abs(ea))])
xlabel('log_2(X)')
ylabel('|X-\omega-log(\omega)|')
title('Absolute Error')



% Complex plane, wrightOmega plotted at lower resolution because so much slower

% -6*pi to 6*pi in real and imaginarary
v = linspace(-6*pi,6*pi,1e3);
[x,y] = meshgrid(v);
d = v(2)-v(1);
z = x+d*randn(s,size(x))+(y+d*randn(s,size(x)))*1i;
w = wrightOmegaq(z);
ea = abs(z-wrightOmegaqInv(w));
em = abs(ea-eps(abs(z)));
em(em < eps) = eps;

ear = abs(real(z)-real(wrightOmegaqInv(w)));
emr = abs(ear-eps(abs(real(z))));
emr(emr < eps) = eps;
eai = abs(imag(z)-imag(wrightOmegaqInv(w)));
emi = abs(ea-eps(abs(imag(z))));
emi(emi < eps) = eps;

v2 = linspace(-6*pi,6*pi,5e1);
[x2,y2] = meshgrid(v2);
d2 = v2(2)-v2(1);
z2 = x2+d2*randn(s,size(x2))+(y2+d2*randn(s,size(x2)))*1i;
w2 = wrightOmega(z2);
ea2 = abs(z2-(wrightOmegaqInv(w2)));

ear2 = abs(real(z2)-real(wrightOmegaqInv(w2)));
eai2 = abs(imag(z2)-imag(wrightOmegaqInv(w2)));

% Magnitude
figure
subplot(221)
imagesc(v,v,abs(w))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','|\omega(Z)|')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Magnitude')

subplot(222)
imagesc(v,v,log10(em))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','Log_1_0(|Z-\omega-log(\omega)|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Absolute Error > EPS(|Z|)')

subplot(223)
imagesc(v,v,log10(ea))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','log_1_0(|Z-\omega-log(\omega)|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Log of Absolute Error')

subplot(224)
imagesc(v2,v2,log10(ea2))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','Log_1_0(|Z-\omega-log(\omega)|)')
xlabel('X')
ylabel('Y')
title('wrightOmega: Log of Absolute Error')

% Magnitude of real part
figure
subplot(221)
imagesc(v,v,abs(real(w)))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','|Re(\omega(Z))|')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Magnitude of Real Part')

subplot(222)
imagesc(v,v,log10(emr))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','Log_1_0(|Re(Z-\omega-log(\omega))|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Absolute Error > EPS(|Z|)')

subplot(223)
imagesc(v,v,log10(ear))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','log_1_0(|Re(Z-\omega-log(\omega))|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Log of Absolute Error')

subplot(224)
imagesc(v2,v2,log10(ear2))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','Log_1_0(|Re(Z-\omega-log(\omega))|)')
xlabel('X')
ylabel('Y')
title('wrightOmega: Log of Absolute Error')

% Magnitude of imaginary part
figure
subplot(221)
imagesc(v,v,abs(imag(w)))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','|Im(\omega(Z))|')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Magnitude of Imaginary Part')

subplot(222)
imagesc(v,v,log10(emi))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','Log_1_0(|Im(Z-\omega-log(\omega))|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Absolute Error > EPS(|Z|)')

subplot(223)
imagesc(v,v,log10(eai))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','log_1_0(|Im(Z-\omega-log(\omega))|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Log of Absolute Error')

subplot(224)
imagesc(v2,v2,log10(eai2))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','Log_1_0(|Im(Z-\omega-log(\omega))|)')
xlabel('X')
ylabel('Y')
title('wrightOmega: Log of Absolute Error')


% -800 to 800 in real and imaginarary
v = linspace(-800,800,1e3);
[x,y] = meshgrid(v);
d = v(2)-v(1);
z = x+d*randn(s,size(x))+(y+d*randn(s,size(x)))*1i;
w = wrightOmegaq(z);
ea = abs(z-wrightOmegaqInv(w));
em = abs(ea-eps(abs(z)));
em(em < eps) = eps;

ear = abs(real(z)-real(wrightOmegaqInv(w)));
emr = abs(ear-eps(abs(real(z))));
emr(emr < eps) = eps;
eai = abs(imag(z)-imag(wrightOmegaqInv(w)));
emi = abs(ea-eps(abs(imag(z))));
emi(emi < eps) = eps;

v2 = linspace(-800,800,5e1);
[x2,y2] = meshgrid(v2);
d2 = v2(2)-v2(1);
z2 = x2+d2*randn(s,size(x2))+(y2+d2*randn(s,size(x2)))*1i;
w2 = wrightOmega(z2);
ea2 = abs(z2-wrightOmegaqInv(w2));

ear2 = abs(real(z2)-real(wrightOmegaqInv(w2)));
eai2 = abs(imag(z2)-imag(wrightOmegaqInv(w2)));

% Magnitiude
figure
subplot(221)
imagesc(v,v,abs(w))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','|\omega(Z)|')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Magnitude')

subplot(222)
imagesc(v,v,log10(em))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','Log_1_0(|Z-\omega-log(\omega)|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Absolute Error > EPS(|Z|)')

subplot(223)
imagesc(v,v,log10(ea))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','log_1_0(|Z-\omega-log(\omega)|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Log of Absolute Error')

subplot(224)
imagesc(v2,v2,log10(ea2))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','Log_1_0(|Z-\omega-log(\omega)|)')
xlabel('X')
ylabel('Y')
title('wrightOmega: Log of Absolute Error')

% Magnitiude of real part
figure
subplot(221)
imagesc(v,v,abs(real(w)))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','|Re(\omega(Z))|')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Magnitude of Real Part')

subplot(222)
imagesc(v,v,log10(emr))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','Log_1_0(|Re(Z-\omega-log(\omega))|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Absolute Error > EPS(|Z|)')

subplot(223)
imagesc(v,v,log10(ear))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','log_1_0(|Re(Z-\omega-log(\omega))|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Log of Absolute Error')

subplot(224)
imagesc(v2,v2,log10(ear2))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','Log_1_0(|Re(Z-\omega-log(\omega))|)')
xlabel('X')
ylabel('Y')
title('wrightOmega: Log of Absolute Error')

% Magnitiude of imaginary part
figure
subplot(221)
imagesc(v,v,abs(imag(w)))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','|Im(\omega(Z))|')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Magnitude of Imaginary Part')

subplot(222)
imagesc(v,v,log10(emi))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','Log_1_0(|Im(Z-\omega-log(\omega))|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Absolute Error > EPS(|Z|)')

subplot(223)
imagesc(v,v,log10(eai))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','log_1_0(|Im(Z-\omega-log(\omega))|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Log of Absolute Error')

subplot(224)
imagesc(v2,v2,log10(eai2))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','Log_1_0(|Im(Z-\omega-log(\omega))|)')
xlabel('X')
ylabel('Y')
title('wrightOmega: Log of Absolute Error')



% Regularization test around A+/-pi*i, A <= -1
vx = linspace(-2.5,-0.5,1e3);
vy = linspace(pi-0.1,pi+0.1,1e3);
[x,y] = meshgrid(vx,vy);
dx = vx(2)-vx(1);
dy = vy(2)-vy(1);
z = x+dx*randn(s,size(x))+(y+dy*randn(s,size(y)))*1i;
w = wrightOmegaq(z);
ea = abs(z-wrightOmegaqInv(w));

% Magnitude
figure
subplot(221)
imagesc(vx,vy,abs(w))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','|\omega(Z)|')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Magnitude - Regularization Test, +Pi*i')

subplot(222)
imagesc(vx,vy,log10(ea))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','log_1_0(|Z-\omega-log(\omega)|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Log of Absolute Error, Pi*i')

vx = linspace(-2.5,-0.5,1e3);
vy = linspace(-pi-0.1,-pi+0.1,1e3);
[x,y] = meshgrid(vx,vy);
dx = vx(2)-vx(1);
dy = vy(2)-vy(1);
z = x+dx*randn(s,size(x))+(y+dy*randn(s,size(y)))*1i;
w = wrightOmegaq(z);
ea = abs(z-wrightOmegaqInv(w));

subplot(223)
imagesc(vx,vy,abs(w))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','|\omega(Z)|')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Magnitude - Regularization Test, -Pi*i')

subplot(224)
imagesc(vx,vy,log10(ea))
axis square
hc = colorbar;
set(get(hc,'YLabel'),'String','log_1_0(|Z-\omega-log(\omega)|)')
xlabel('X')
ylabel('Y')
title('wrightOmegaq: Log of Absolute Error, -Pi*i')