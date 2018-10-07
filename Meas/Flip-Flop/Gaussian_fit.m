function h = Gaussian_fit(params,t)
%Parameters = [Width Position Peakheight Zerolevel]

a = params(1); 
b = params(2); 
c = params(3);
d = params(4);

h=c.*exp(-((t-b)./a).^2)+d;
 
end