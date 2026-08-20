function invertspec(filename)
% INVERTSPEC    Invert a JASON spectrum point-by-point
%   INVERTSPEC(filename) inverts a spectrum in the JASON .jjh5 
%       file given by filename
%
%   Can be called from the JASON external command processing item as:
%     Cmd:       matlab
%     Arguments: -nojvm -batch "invertspec('$TMPFILE')"
% 
%   (Double quotes should be omited on macOS)

dataset = '/JasonDocument/DataPoints/';

disp(['Opening dataset: ', filename])
dsetre = h5read(filename, strcat(dataset, '/0'));
dsetim = h5read(filename, strcat(dataset, '/1'));

disp(['dataset size: ', num2str(length(dsetre))])

h5write(filename, strcat(dataset, '/0'), -dsetre);
h5write(filename, strcat(dataset, '/1'), -dsetim);
disp('dataset changed')

end
