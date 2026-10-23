function b=closeout_read_bie_snapshot(file,collection,index)
% Read only numeric datasets of a selected HDF5 struct reference.
fid=H5F.open(file,'H5F_ACC_RDONLY','H5P_DEFAULT'); cl=onCleanup(@()H5F.close(fid));
did=H5D.open(fid,['/validation/' collection]); dl=onCleanup(@()H5D.close(did)); refs=H5D.read(did);
assert(index<=size(refs,2)); obj=H5R.dereference(did,'H5R_OBJECT',refs(:,index));
group=H5I.get_name(obj); H5O.close(obj); info=h5info(file,group); b=struct;
for j=1:numel(info.Datasets)
    name=info.Datasets(j).Name; path=[group '/' name]; inf=h5info(file,path);
    if any(strcmp(inf.Datatype.Class,{'H5T_FLOAT','H5T_INTEGER','H5T_COMPOUND'}))
        value=h5read(file,path);
        if isstruct(value) && all(isfield(value,{'real','imag'})), value=complex(value.real,value.imag); end
        if isnumeric(value), b.(name)=value; end
    end
end
b.hdf5_group=group; b.collection=collection; b.index=index; b.source_file=file;
end
