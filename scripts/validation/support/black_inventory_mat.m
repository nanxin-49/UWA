function inventory = black_inventory_mat(filename)
%BLACK_INVENTORY_MAT Inspect HDF5 classes before MATLAB deserialization.
info = h5info(filename);
classes = local_classes(info);
allowed = {'struct','cell','char','double','single','logical', ...
    'uint64','uint32','uint16','uint8','int64','int32','int16','int8'};
% MATLAB's HDF5 reference pool uses "canonical empty" for empty numeric data.
assert(all(ismember(classes,[allowed,{'canonical empty'}])), 'black:MatClass', ...
    'MAT file contains unexpected classes; do not load it.');
variables = whos('-file',filename);
assert(all(ismember({variables.class},allowed)), 'black:MatVariables','Unexpected top-level MAT variable.');
inventory = struct('filename',filename,'hdf5_classes',{unique(classes)},'variables',variables);
end

function classes = local_classes(node)
classes = {};
for ii = 1:numel(node.Attributes)
    attribute = node.Attributes(ii);
    assert(~contains(lower(attribute.Name),'function') && ...
        ~contains(lower(attribute.Name),'object'), 'black:MatObject', 'Unexpected object attribute.');
    if strcmp(attribute.Name,'MATLAB_class')
        classes{end+1} = char(attribute.Value(:).'); %#ok<AGROW>
    end
end
if isfield(node,'Datasets')
    for ii = 1:numel(node.Datasets)
        classes = [classes,local_classes(node.Datasets(ii))]; %#ok<AGROW>
    end
end
if isfield(node,'Groups')
    for ii = 1:numel(node.Groups)
        classes = [classes,local_classes(node.Groups(ii))]; %#ok<AGROW>
    end
end
end
