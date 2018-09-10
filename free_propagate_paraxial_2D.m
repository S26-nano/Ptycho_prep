function [field_paraxial, x_there, y_there, delta_x_there, delta_y_there] = free_propagate_paraxial_2D(field_in, x_here, y_here,distance, k, lambda, fft_operations) %propagation paraxiale

[Nx,Ny] = size(field_in) ;
x0 = x_here(Nx/2+1,Ny/2+1) ;
y0 = y_here(Nx/2+1,Ny/2+1) ;
x_here = x_here-x0 ;
y_here = y_here-y0 ;
delta_x_here = x_here(2,1)-x_here(1,1) ;
delta_y_here = y_here(1,2)-y_here(1,1) ;



%figure(10); %subplot(2,2,1); %imagesc(abs(B.*field_in)); colorbar; %subplot(2,2,2); %imagesc(abs(fourier)); colorbar; clear B field_in

delta_qx = 2*pi/(x_here(Nx, Ny)-x_here(1,1)) ; delta_qy = 2*pi/(y_here(Nx, Ny)-y_here(1,1)) ;

fourier_space_x = x_here/delta_x_here*delta_qx ;
fourier_space_y = y_here/delta_y_here*delta_qy ;

if distance>0
    x_there = fourier_space_x*distance/k ;                                    
    y_there = fourier_space_y*distance/k ;
else
    x_there = fourier_space_x*abs(distance)/k ;                                    
    y_there = fourier_space_y*abs(distance)/k ;

end;
%display([x0,Nx,delta_x_here])
%display([fourier_space_x(1,1),fourier_space_x(512,512),k,lambda])

%one FT method, good for intermediate distance propa
if(fft_operations==1)
    B = exp(1i*(k/2)*(x_here.^2+y_here.^2)/distance) ;

    if distance>0
        fourier = fftshift(fft2(fftshift(B.*field_in)))/(Nx*Ny)^(1/2) ;
    else
        fourier = fftshift(ifft2(fftshift(B.*field_in)))*(Nx*Ny)^(1/2) ;
    end;

    A = -1i*exp(1i*k*distance).*exp(1i*k*(x_there.^2+y_there.^2)/(2*distance));
    field_paraxial = A.*fourier ;

    x_there = x_there+x0 ;
    y_there = y_there+y0 ;
    delta_x_there = x_there(2,1)-x_there(1,1) ; delta_y_there = y_there(1,2)-y_there(1,1) ;
end

%two FT method, good for very short distance propagations
if(fft_operations==2)
    
    [Xf Yf] = meshgrid([-Nx/2:Nx/2-1]/(Nx*delta_x_here));
    B2 = exp(-1i*pi*lambda*distance * (Xf.^2 + Yf.^2));
    Psi_2ft = fftshift(fftn(fftshift(field_in))) .* B2;
    Psi_2ft = exp(i*k*distance) * fftshift(ifftn(fftshift(Psi_2ft)));
    field_paraxial = Psi_2ft; 

    x_there = x_here + x0;
    y_there = x_here + y0;
    delta_x_there= delta_x_here;
    delta_y_there= delta_y_here;
end



end

