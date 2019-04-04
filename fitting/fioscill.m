function fp = fioscill(x, y, ctrl)
% function fp = fioscill(x, y, ctrl)
% ctrl = 1: offset, cos coef, sin coef, freq, shift (always 0), decay prefac = 1/mean(abs(x));
%   fp: 
% ctrl = 2: offset, amplitude, phase, freq, shift, decay prefac = 1/mean(abs(x));
%   fp: 

switch ctrl 
    case 1
        dx = mean(diff(x));
        window=1;
        ft = fft((y-mean(y)).*window) .* exp(1i * x(1) * (0:length(x)-1) * 2*pi /(dx*length(x)))/length(x);
        % hack to get phase shift, factors determined experimentally
        ft = ft(1:round(end/2));
        [~, mi] = max(abs(ft(2:end)));
        xi=max(1,mi-2):min(length(ft)-1,mi+2);           
        mia=sum(abs(ft(xi+1)).*xi)/sum(abs(ft(xi+1)));
        fp(4) = 2* pi* (mia)/((length(x)+1) * dx);        
        f=ft(mi+1);
        fp(3) = -3*real(f);
        fp(2) = -3*imag(f);        
        fp(1) = mean(y);
        fp(5) = 0;
        fp(6) = .5/mean(abs(x));        
    case 2
        dx = mean(diff(x));
        window=1;
        ft = fft((y-mean(y)).*window) .* exp(1i * x(1) * (0:length(x)-1) * 2*pi /(dx*length(x)))/length(x);
        ft = ft(1:round(end/2));
        [~, mi] = max(abs(ft(2:end)));
        xi=max(1,mi-2):min(length(ft)-1,mi+2);           
        mia=sum(abs(ft(xi+1)).*xi)/sum(abs(ft(xi+1)));
        fp(4) = 2* pi* (mia)/((length(x)+1) * dx);        
        f=ft(mi+1);
        fp(3) = angle(f);
        fp(2) = -3*abs(f);        
        fp(1) = mean(y);
        fp(5) = 0;
        fp(6) = .5/mean(abs(x));  
end