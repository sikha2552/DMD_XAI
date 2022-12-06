function [XLow XSparse rl iml TD Phi lamda ModeAmplitudes ModeFrequencies GrowthRates orl oimg]= DMDsvd(Y)
nF = size(Y,2)
vF=Y;
X1 = vF(:,1:nF-1);
X2 = vF(:,2:nF);
[U, E, V] = svd(X1, 'econ');
figure;plot(diag(E)/sum(diag(E)),'ro');
disp('1 reached')
r=nF-1;
%r=250;
 Ur = U(:,1:r);
 Sr = E(1:r,1:r);
 Vr = V(:,1:r);
 size(Ur)
 size(Vr)
 size(Sr)
 %Stild = U' * X2 * V * pinv(E);
 Stild =Ur'*X2*Vr/Sr;
[EigVec, EigVal] = eig(Stild);
Phi = X2*Vr/Sr*EigVec; % DMD modes
%%  Plot DMD spectrum -Discrete Lamda
%%%%%%%%%%%%Signal 1
%neeth code
dt=1;
eigs=EigVal;
lamda = diag(eigs);

L = diag(eigs);
omega=log(lamda)/dt;
% figure();plot(L, 'k*');rectangle('Position', [-1 -1 2 2], 'Curvature', 1, ...
%     'EdgeColor', 'r', 'LineStyle', '-.');axis(1.2*[-1 1 -1 1]);axis square;axis tight;title('Visualisation of \lambda')

figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(lamda),imag(lamda),'o');title('DMD Spectrum-Dicrete Eigen value');
axis([-6.1 6.1 -6.1 6.1]);
rl=real(lamda);
iml=imag(lamda);
%figure;eigplot(rl,iml)
% Plot Continuos eigen values (PHI)
figure; scatter(real(omega)/(2*pi),imag(omega)/(2*pi)); grid on; axis equal;title('DMD Spectrum-Continuos Eigen value');
figure;plot(omega./(2*3.14),'r*','Markersize',12); grid on; axis equal;
orl=real(omega);
oimg=imag(omega);

omega = log(EigVal);
Omega = exp(omega);
Psi = Ur * EigVec;
%b = pinv(Phi)*X1(:, 1);
b = Phi\X1(:, 1);
%XDMDt = Psi * Omega * b;
omegaD = abs(diag(omega));
LowRank = exp(min(omegaD))
disp LowRank
%XDMD = zeros(height * width, nF - 1);
for t = 1:nF 
    XDMD(:, t) = Psi * Omega.^t * b;%change Phi to Psi
end
disp('2 reached')
%XLow = zeros(height * width, nF - 1);
for t = 1 :nF 
        XLow(:, t) = Psi * LowRank.^t * b;%change Phi to Psi
end
disp('3 reached')
XLow = abs(XLow);
XSparse = abs(XDMD - XLow);
T=size(X1,2);
%%%%%%%Signal 1
D = X1(:, 1);       % time = 0
B= Phi\D;   % initial value; \: pseudo inverse
TD = zeros(r, T);
size(B)
time=(0:T-1)*dt;
omega=log(lamda)/dt;
for i = 1:T 
   TD(:, i) = (B.*exp(omega*time(i))); 
end
% Gets the mode amplitudes
    ModeAmplitudes=Phi\X1(:,1);
    
    %Gets the frequencies associated with the modes
    fNY=1/(2*dt);
    ModeFrequencies=(angle(L)/pi)*fNY;
        %Gets the growth rates    
    GrowthRates=log(abs(L))/dt;
end