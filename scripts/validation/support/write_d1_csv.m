function write_d1_csv(T,file)
%WRITE_D1_CSV D1 tables with 17-digit double text; strict audits remain unchanged.
arguments
    T table
    file (1,:) char
end
assert(~isfile(file),'d1:CSVOverwrite','New CSV path required');
fid=fopen(file,'w','n','UTF-8'); assert(fid>=0); guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',strjoin(T.Properties.VariableNames,','));
for row=1:height(T)
    for col=1:width(T)
        if col>1, fprintf(fid,','); end
        v=T{row,col};
        if iscell(v), assert(isscalar(v)); v=v{1}; end
        if isnumeric(v)
            assert(isa(v,'double') && isscalar(v) && isreal(v),'D1 CSV requires real scalar doubles');
            fprintf(fid,'%.17g',v);
        elseif islogical(v)
            assert(isscalar(v)); fprintf(fid,'%d',v);
        else
            s=char(string(v)); assert(isrow(s) || isempty(s));
            fprintf(fid,'"%s"',strrep(s,'"','""'));
        end
    end
    fprintf(fid,'\n');
end
end
