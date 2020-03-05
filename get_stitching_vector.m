function [Stitching_transformation,Intensity_peak] = get_stitching_vector(Ref_Image,Moving_Image,Hanning_window)
%Compute best translation vector for the stitching of two neighbours
%pictures
%Based on phase correlation method and fast-fourrier transform

if nargin < 3
    Hanning_window = 'Inversed';
end

image_width = size(Ref_Image,1);
image_height = size(Ref_Image,2);

%%Performing 0 padding to have the exact value of the correct translation
%%vector...

Mov_X = size(Moving_Image,1);
Mov_Y = size(Moving_Image,2);
Ref_X = size(Ref_Image,1); 
Ref_Y = size(Ref_Image,2);

if Hanning_window=="Centered"
    %Computing the hanning window :
    Hanning_mov_x = 0.5*(1-cos(2*pi*(1:Mov_X)/Mov_X)) ;
    Hanning_mov_y = 0.5*(1-cos(2*pi*(1:Mov_Y)/Mov_Y)) ;
    Hanning_mov = repmat(Hanning_mov_x,Mov_X,1).*repmat(Hanning_mov_y.',1,Mov_Y);
    
    Hanning_ref_x = 0.5*(1-cos(2*pi*(1:Ref_X)/Ref_X)) ;
    Hanning_ref_y = 0.5*(1-cos(2*pi*(1:Ref_Y)/Ref_Y)) ;
    Hanning_ref = repmat(Hanning_ref_x,Ref_X,1).*repmat(Hanning_ref_y.',1,Ref_Y);
    
    Moving_Image = Moving_Image.* Hanning_mov;
    Ref_Image = Ref_Image.*Hanning_ref;

end


if Hanning_window=="Inversed"
    %Computing the hanning window :
    Hanning_mov_x = 0.5*(1-cos(2*pi*(1:Mov_X)/Mov_X)) ;
    Hanning_mov_y = 0.5*(1-cos(2*pi*(1:Mov_Y)/Mov_Y)) ;
    Hanning_mov = repmat(Hanning_mov_x,Mov_X,1).*repmat(Hanning_mov_y.',1,Mov_Y);
    
    Hanning_ref_x = 0.5*(1-cos(2*pi*(1:Ref_X)/Ref_X)) ;
    Hanning_ref_y = 0.5*(1-cos(2*pi*(1:Ref_Y)/Ref_Y)) ;
    Hanning_ref = repmat(Hanning_ref_x,Ref_X,1).*repmat(Hanning_ref_y.',1,Ref_Y);
    
    Moving_Image = Moving_Image.* (1-Hanning_mov);
    Ref_Image = Ref_Image.*(1-Hanning_ref);

end


Moving_Image_pad = padarray(Moving_Image, [Mov_X,Mov_Y],'post');
Ref_Image_pad = padarray(Ref_Image, [Ref_X,Ref_Y],'post');


%Computing cross-correlation using FFT transform
FFT_Ref = fft2(Ref_Image_pad);
FFT_Moving = fft2(Moving_Image_pad);
IMF = FFT_Ref.*conj(FFT_Moving);
CPS = IMF./(abs(FFT_Moving).*abs(FFT_Ref));

Magnitude = (ifft2(CPS,'symmetric'));
Magnitude = Magnitude - min(Magnitude(:));
Magnitude = fftshift(Magnitude);

%Detecting the peak and quantifying its intensity

[X Y] = find((Magnitude == (max(max(gather(Magnitude))))));
X= X-Mov_X;
X = gather(X);
Y = Y-Mov_Y;
Y = gather(Y);

%For the peak intensity : normalised cross correlation is computed
U = Ref_Image;
U = U -mean(U(:));

V = circshift(Moving_Image_pad, [X Y]);
V = V(1:Ref_X,1:Ref_Y);
V = V - mean(V(:));


Normalised_cross_product = U.*V;
Normalised_cross_product = sum(Normalised_cross_product(:)) / sqrt((sum(V.^2,'all')*sum(U.^2,'all')));
Normalised_cross_product = gather(Normalised_cross_product);

Intensity_peak = gather(Normalised_cross_product);

%Creating the affine2D object that we can apply on 
Stitching_transformation = affine2d([1 0 0 ; 0 1 0 ; Y X 1]);


end

