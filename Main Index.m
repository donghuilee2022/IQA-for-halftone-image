function IMADVal = IMAD(ref,dst)
%Improved MAD value caculation
% d(SDSIM) is HiIndex (SDSIM map)
% d(appear) is LoIndex (log-Gabor statistics comparsion)
% Input should be RGB color image

%The flowchat of the algorithm, part1, part2, part3.
%----Part1: Color space transformation----
%image to double
OriImg = double(ref);
DisImg = double(dst);
% RGB to Lab
ori_imgLAB = rgb2lab(OriImg);
dis_imgLAB = rgb2lab(DisImg);
%reference image channel seperation
ori_imgL = ori_imgLAB(:,:,1);
ori_imgA = ori_imgLAB(:,:,2);
ori_imgB = ori_imgLAB(:,:,3);
%distored image channel seperation
dis_imgL = dis_imgLAB(:,:,1);
dis_imgA = dis_imgLAB(:,:,2);
dis_imgB = dis_imgLAB(:,:,3);

%----Part2: Human Visual Frequency Response----
a = 131.6;
b = 0.3188;
c = 0.525;
d = 3.91;
%Average luminance value of L channel image
oriL_AverLumi = mean2(ori_imgL);
disL_AverLumi = mean2(dis_imgL);
% alpha_Lu
oriL_alpha_Lu =1/(c*log(oriL_AverLumi) + d);
disL_alpha_Lu =1/(c*log(disL_AverLumi) + d);
% K_Lu
oriL_K_Lu = a*(oriL_AverLumi^b);
disL_K_Lu = a*(disL_AverLumi^b);

[h, w] = size(ori_imgL);

%images processing with HVS filters
refLf = real( ifft2( ifftshift( fftshift( fft2( ori_imgL ) ).* hvs_filter1(h,w,32,oriL_K_Lu,oriL_alpha_Lu) ) ) );
refAf = real( ifft2( ifftshift( fftshift( fft2( ori_imgA ) ).* hvs_filter2(h,w,32) ) ) );
refBf = real( ifft2( ifftshift( fftshift( fft2( ori_imgB ) ).* hvs_filter2(h,w,32) ) ) );
dstLf = real( ifft2( ifftshift( fftshift( fft2( dis_imgL ) ).* hvs_filter1(h,w,32,disL_K_Lu,disL_alpha_Lu) ) ) );
dstAf = real( ifft2( ifftshift( fftshift( fft2( dis_imgA ) ).* hvs_filter2(h,w,32) ) ) );
dstBf = real( ifft2( ifftshift( fftshift( fft2( dis_imgB ) ).* hvs_filter2(h,w,32) ) ) );
%channel images combination
refRGB = lab2rgb(cat(3,refLf,refAf,refBf));
dstRGB = lab2rgb(cat(3,dstLf,dstAf,dstBf));

%IMAD index values caculation for R, G and B
IMAD_R = IMADComb(refRGB(:,:,1),dstRGB(:,:,1));
IMAD_G = IMADComb(refRGB(:,:,2),dstRGB(:,:,2));
IMAD_B = IMADComb(refRGB(:,:,3),dstRGB(:,:,3));

%average value caculation using index ssim_WG

IMADVal =1/3*(IMAD_R+IMAD_G+IMAD_B);

%----Part3: Functions impletation----
%--------------------------------------------------------------------------
% -----HVS filter for L channel image-----
function [res1] = hvs_filter1(x, y, nfreq, K_L, alpha_Lu)
[xplane,yplane]=meshgrid(-x/2+0.5:x/2-0.5, -y/2+0.5:y/2-0.5);	
plane=(xplane+1i*yplane)/y*2*nfreq;
radfreq=abs(plane);	
omega=0.7;
s=(1-omega)/2*cos(4*angle(plane))+(1+omega)/2;
radfreq=radfreq./s;

% Now generate the hvs_filter
% csf = 2.6*(0.0192+0.114*radfreq).*exp(-(0.114*radfreq).^1.1);
hvs_f1 = K_L*exp(-(alpha_Lu*radfreq));

res1 = hvs_f1;
end

% -----HVS filter for A and B channel images-----
function [res2] = hvs_filter2(x, y, nfreq)
[xplane,yplane]=meshgrid(-x/2+0.5:x/2-0.5, -y/2+0.5:y/2-0.5);	
plane=(xplane+1i*yplane)/y*2*nfreq;
radfreq=abs(plane);	
alpha_Lu=0.419;
A=100;

% Now generate the hvs_filter
% csf = 2.6*(0.0192+0.114*radfreq).*exp(-(0.114*radfreq).^1.1);
hvs_f2 = A*exp(-(alpha_Lu*radfreq));

res2 = hvs_f2;
end

%combination part of function IMAD
function [Val] = IMADComb(OrgImg,DstImg)
HiIndex = hi_index(OrgImg, DstImg);
LoIndex = lo_index(OrgImg, DstImg); 

b1        = exp(-2.55/3.35);
b2        = 1/(log(10)*3.35);     
alpha       = 1 ./ ( 1 + b1*(HiIndex).^b2 ) ;

% IMAD values
Val = HiIndex.^(alpha) * LoIndex.^(1-alpha);
% disp(HiIndex);
% disp(LoIndex);
% disp('------');
end
end
%--------------------------------------------------------------------------