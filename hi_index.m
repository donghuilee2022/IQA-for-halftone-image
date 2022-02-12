function [I, mp, imgs] = hi_index( ref_img, dst_img )
% Calculate the high quality index by calculating the masking map and
% then approximating the local statistics, using filtering.
% Author: Eric Larson
% Department of Electrical and Computer Engineering
% Oklahoma State University, 2008
% Image Coding and Analysis Lab

% Masking/luminance parameters
k = 0.02874;
G = 0.5 ;       % luminance threshold
C_slope = 1;    % slope of detection threshold
Ci_thrsh= -5;   % contrast to start slope, rather than const threshold
Cd_thrsh= -5;   % saturated threshold
ms_scale= 1;    % scaling constant
if(nargin==0)% for debug only,
  I={'collapsing (mask) and raw lmse',...
    'using c code',...
    'using log scale',...
    'two norm',...
    sprintf('Ci = %.2f',Ci_thrsh),...
    sprintf('Cslope = %.2f',C_slope),...
    sprintf('Cd = %.2f',Cd_thrsh),...
    sprintf('G = %.2f',G),...
    'LO QUALITY',...
    };
  return;
end
% log(Contrast of  ref-dst)   vs.   log(Contrast of reference)
%              /
%             /
%            /  _| <- Cslope
%           /
%----------+ <- Cdthresh (y axis height)
%          /\
%          ||
%       Ci thresh (x axis value)

% TAKE TO LUMINANCE DOMAIN USING LUT
if isinteger(ref_img)
  LUT = 0:1:255; %generate LUT
  LUT = k .* LUT .^ (2.2/3);
  ref = LUT( ref_img + 1 );
  dst = LUT( dst_img + 1 );
else % don't use the speed up
  ref = k .* ref_img .^ (2.2/3);
  dst = k .* dst_img .^ (2.2/3);
end

[M N] = size( ref );

% ACCOUNT FOR CONTRAST SENSITIVITY
% csf = make_csf( M, N, 32 )';
% ref = real( ifft2( ifftshift( fftshift( fft2( ref ) ).* csf ) ) );
% dst = real( ifft2( ifftshift( fftshift( fft2( dst ) ).* csf ) ) );
refS = ref;
dstS = dst;

% Use c code to get fast local stats
[std_2 std_1 m1_1] = ical_std( dst-ref, ref );

BSIZE = 16;

Ci_ref = log(std_1./m1_1); % contrast of reference (also a measure of entropy)
Ci_dst = log(std_2./m1_1); % contrast of distortions (also a measure of entropy)
Ci_dst( find( m1_1 < G ) ) = -inf;

msk       = zeros( size(Ci_dst) );
idx1      = find( (Ci_ref > Ci_thrsh) ...
  & (Ci_dst > (C_slope * (Ci_ref - Ci_thrsh) + Cd_thrsh) ) );
idx2      = find( (Ci_ref <= Ci_thrsh) & (Ci_dst > Cd_thrsh) );

msk(idx1) = ( Ci_dst( idx1 ) - (C_slope * (Ci_ref( idx1 )-Ci_thrsh) + Cd_thrsh) ) ./ ms_scale;
msk(idx2) = ( Ci_dst( idx2 ) - Cd_thrsh ) ./ ms_scale;
%= ( Contrast of heavy Dst - 0.75 * Contrast Ref ) / normalize
%= ( Contrast of low Dst  - Threshold ) / normalize

% Use lnonssim(local non-SSIM) and weight by distortion mask
win = ones( BSIZE ) ./ BSIZE^2;
% [~,ssimmap] = ssim(ref_img,dst_img);
[~,ssimmap] = ssim(ref_img,dst_img);

% % Structural dissimilarity is filtered by a 16 x 16 window
% lsdsim  = ( imfilter( ((1.00-double(ssimmap))*dst_img).^2, ...
%   win,  'symmetric', 'same', 'conv' ) );

lsdsim  = ( imfilter( (double(ssimmap)).^2, ...
  win,  'symmetric', 'same', 'conv' ) );

% Visibility caculation of structural dissimilarity
mp1    = msk .* lsdsim;

% kill the edges
mp2 = mp1( BSIZE+1:end-BSIZE-1, BSIZE+1:end-BSIZE-1);

I = norm( mp2(:) , 2 ) ./ sqrt( length( mp2(:) ) ) .* 10;

if( nargout > 2)
  imgs.ref = refS;
  imgs.dst = dstS;
end
end
