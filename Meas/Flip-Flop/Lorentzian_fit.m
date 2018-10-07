function h = Lorentzian_fit(params,t)
%Parameters = [Width Position Peakheight Zerolevel]

a = params(1); 
b = params(2); 
c = params(3).*a;
d = params(4);
h = c.*(a./((t-b).^2+a^2))+d;
end
