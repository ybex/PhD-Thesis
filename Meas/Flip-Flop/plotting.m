%% Mateuzs data

hfig = openfig('combined_spectrum.fig');
d = findall(hfig, '-property', 'xdata');
xydatas = arrayfun(@(h) get(h, {'xdata','ydata', 'type'}), d, 'Uniform', 0);

%% Arne data
