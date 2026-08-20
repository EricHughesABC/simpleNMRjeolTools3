function doublespec(filename)
% DOUBLESPEC    Double a JASON spectrum by appending an inverted copy
%   DOUBLESPEC(filename) adds an inverted copy to the end of a 
%       spectrum in the JASON .jjh5 file given by filename
%
%   Can be called from the JASON external command processing item as:
%     Cmd:       matlab
%     Arguments: -nojvm -batch "doublespec('$TMPFILE')"
% 
%   (Double quotes should be omited on macOS)

dataset = '/JasonDocument/DataPoints/';

disp(['Opening dataset: ', filename])
dsetre = h5read(filename, strcat(dataset, '0'));
dsetim = h5read(filename, strcat(dataset, '1'));
dlength = h5readatt(filename, '/JasonDocument', 'Length');

disp(['dataset size: ', num2str(dlength(1))])

dsetre = [dsetre, -dsetre];
dsetim = [dsetim, -dsetim];
dlength(1) = 2.0 * dlength(1);
disp(['dataset size: ', num2str(dlength(1))])

fileid = H5F.open(filename, "H5F_ACC_RDWR", "H5P_DEFAULT");
H5L.delete(fileid, strcat(dataset, '0'), "H5P_DEFAULT")
H5L.delete(fileid, strcat(dataset, '1'), "H5P_DEFAULT")
H5F.close(fileid);

h5create(filename, strcat(dataset, '0'), double(dlength(1)));
h5create(filename, strcat(dataset, '1'), double(dlength(1)));
h5write(filename, strcat(dataset, '0'), dsetre)
h5write(filename, strcat(dataset, '1'), dsetim)
h5writeatt(filename, '/JasonDocument', 'Length', dlength);
disp('dataset changed')

end
