close all;
clear all;
clc;
path1 = 'D:\RESEARCH\POST_PHD\UPC\000_DMD_DL\Poorly_trained\plane_68\';
images=dir(strcat(path1,'*.png'));
file_name=images.name;
I=strcat(path1,file_name);
img=imread(I);
img=imresize(img,[128,128]);
[m n c]=size(img);
img1=zeros(128,128,3);
img1=uint8(img1);
x=zeros(m*n*c,numel(images));
% x=zeros(m*n*c,1000);
% % % 
 i=1;
 j=1;
 for ks=1:numel(images)
 file_name=images(ks).name;
 I=strcat(path1,file_name);
 img=imread(I);
 img=imresize(img,[128,128]);
% figure;imshow(img);
 x(i:m*n*c,j)=img(:);
 j=j+1;
 %imshow(I)
 img1=img1+img;
 end
% x1=x(:,[2 1 4 3 32 14 6 8 9 10 12]);
 x1=x;
[sp1 low1 rl iml Psi Phi lamda ModeAmplitudes ModeFrequencies GrowthRates orl oimg]=dmdcompute_XAI(x1,m,n,c);
%  figure;imshow(sp1,[]);title('sparse');
%  figure;imshow(low1,[]);title('lowrank');
 sm1=imfilter(sp1, fspecial('gaussian', [3,3], .5));
 sm1=(sm1-min(sm1(:)))./(max(sm1(:))-min(sm1(:)));
%  figure;imshow(sm1);
sm1 = morphSmooth(sm1,5);
 final1 = enhanceContrast(sm1, 10);
 figure;imshow(final1);title('Sparse');
 sm2=imfilter(low1, fspecial('gaussian', [3,3], .5));
 sm2=(sm2-min(sm2(:)))./(max(sm2(:))-min(sm2(:)));
%  figure;imshow(sm1);
%  %sm1 = morphSmooth(sm1,5);
 final2= enhanceContrast(sm2, 25);
 figure;imshow(final2,jet(256));title('Low-rank');
 final1=uint8(1 * mat2gray(final1));
 %%
% nbins=10;
% eigs=rl.^2+iml.^2;
% figure;histogram(eigs,nbins);
 %% REconstruction
%Finds dominant oscillatory mode
dominantModeNo=find(ModeAmplitudes==max(ModeAmplitudes));
r1=1,r2=26;
s1= Phi(:,r1:r2)*Psi(r1:r2,:);
re=zeros(m,n,c); rows=m; cols=n;  
%s3=sum(s1,2);
s3=s1(:,1);
re=zeros(m,n,c); rows=m; cols=n;           
R = abs(reshape(s3(1:rows*cols,1),rows,cols));
G = abs(reshape(s3(rows*cols+1:2*rows*cols,1),rows,cols));
B = abs(reshape(s3(2*rows*cols+1:end,:),rows,cols));
Reconstruct= cat(3, mat2gray(R),mat2gray(G),mat2gray(B));
Reconstruct = enhanceContrast(Reconstruct, 15);
figure;imshow(Reconstruct,jet); %title('Reconstructed using' num2str( r1 ) 'to'num2str(r2)); 
str = sprintf('Reconstructed using %d to %d', r1,r2)
title(str)
% % % ds=1,de=ds+1;
% % % r=1;
% % % s1= Phi(:,ds:de)*Psi(ds:de,:); %nth mode reconstruction
% % % s2= Phi(:,1:r)*Psi(1:r,:);
% % % s1=real(s1);
% % % s2=real(s2);
% % % re=zeros(m,n,c);
% % % %figure;
% % % for i=1:size(s2,2)
% % %   Reconstruct(:,:,:,i)=  reshape(s2(:,i), [m, n,c]);
% % %   %figure;
% % %   %imshow(Reconstruct(:,:,:,i),[]);title(strcat(num2str(i),'Sparse'));
% % %   re=re+ Reconstruct(:,:,:,i);
% % % end
% % %  figure;imshow(re,[]);colormap(gca, jet); % Ignore pink map and use jet instead.
% % % colorbar(gca);
%% Entropy
 dis=(rl.*rl)+(iml.*iml);
 entropy(dis)
 dis2=(orl.^2+oimg.^2);
 entropy= rdf(dis2)
 mean(dis)

%%
%==========Mode Amplitudes===========
figure('color', 'w');
stem(abs(ModeFrequencies),abs(ModeAmplitudes), 'k^', 'filled'); hold on
plot(abs(ModeFrequencies(dominantModeNo)),abs(ModeAmplitudes(dominantModeNo)),'ro','MarkerSize',10);
set(gca, 'YScale', 'log')
title('Mode amplitude versus Frequency');
xlabel('f [Hz]');ylabel('Amplitude');

%Radial Distribution
% radialDistribution = zeros(1, length(dis));
% for k = 1 : length(rl)
%     thisDistance = round(dis(k)); % Convert distance into an index.
%     radialDistribution(thisDistance) = radialDistribution(thisDistance) + 1;
% end
% bar(radialDistribution, 'BarWidth', 1);
% grid on;

%  figure;plot(real(s1(:,1)),'blue')
%  final2=uint8(1 * mat2gray(final2));
%  imwrite(final1,'Sparse.png');
%  [xLow xSparse]= DMDsvd(x);
% low=zeros(m,n,c);
% sp=zeros(m,n,c);
% %   for i=1:size(xSparse,2)
% %   lo=  reshape(xLow(:,i), [m, n,c]);
% %   low=low+lo;
% %   imshow(lo,[]);title(strcat(num2str(i),'low'));
% %   end
% for i=1:size(xSparse,2)
%   lo(:,:,:,i)=  reshape(xLow(:,i), [m, n,c]);
%   imshow(lo(:,:,:,i),[]);title(strcat(num2str(i),'low'));
%   low=low+lo(:,:,:,i);
% 
% end
%  figure;imshow(low,[]); 
%  
% for i=1:size(xSparse,2)
%   o(:,:,:,i)=  reshape(xSparse(:,i), [m, n,c]);
%   imshow(o(:,:,:,i),[]);title(strcat(num2str(i),'low'));
%   sp=sp+o(:,:,:,i);
% end
%  figure;imshow(low,[]); 
%  figure;imshow(sp,[]); 