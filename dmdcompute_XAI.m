function[sp low rl iml Psi Phi lamda ModeAmplitudes ModeFrequencies GrowthRates orl oimg ]=dmdcompute_XAI(x,m,n,c)
[xLow xSparse rl iml Psi Phi lamda ModeAmplitudes ModeFrequencies GrowthRates orl oimg]= DMDsvd(x);
size(xLow)
low=zeros(m,n,c);
sp=zeros(m,n,c);
% % %  %-----------------------RESHAPE-----------------------------------------%
for i=1:size(xSparse,2)
  lo(:,:,:,i)=  reshape(xLow(:,i), [m, n,c]);
  %imshow(lo(:,:,:,i),[]);title(strcat(num2str(i),'low'));
  low=low+lo(:,:,:,i);

end
%figure;imshow(low,[]);title('lo alone');
for i=1:size(xSparse,2)
  o(:,:,:,i)=  reshape(xSparse(:,i), [m, n,c]);
  %figure;imshow(o(:,:,:,i),[]);title(strcat(num2str(i),'Sparse'));
  sp=sp+o(:,:,:,i);
end
%figure;imshow(sp,[]);title('sparse alone');



end