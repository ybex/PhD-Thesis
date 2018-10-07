function h = Rabi_fit(params,f)
%Parameters = [freq Position Peakheight Zerolevel]

fres = params(1); 
tp = params(2);
frabi = params(3);
offset = params(4);

delta =  fres-f;

h=frabi^2./(frabi.^2+delta.^2).*sin(2*pi*sqrt(frabi.^2+delta.^2)/2.*tp).^2+offset;
 
end