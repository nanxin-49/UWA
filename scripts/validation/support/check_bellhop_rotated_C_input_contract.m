function [checks,receipt]=check_bellhop_rotated_C_input_contract(root,reference,cfg)
%CHECK_BELLHOP_ROTATED_C_INPUT_CONTRACT Physical SSP identical; derived guard exact.
checks=table(strings(0,1),false(0,1),'VariableNames',{'check','pass'});
receipt=struct('status','INPUT_CONTRACT_PASS','reference_B',reference,'current_C',root,'step_B_m',2*cfg.step_m,'step_C_m',cfg.step_m,'numeric_guard_change',false);
for ext={'.ssp','.iwcov','.rotq','.r25d'}
    old=[reference ext{1}]; current=[root ext{1}]; assert(isfile(old)==isfile(current),'B/C input presence mismatch');
    same=true; if isfile(old), same=strcmpi(black_file_hash(old),black_file_hash(current)); end
    if ~same && strcmp(ext{1},'.ssp')
        [nrB,valuesB]=read_ssp(old); [nrC,valuesC]=read_ssp(current);
        assert(nrB==nrC && numel(valuesB)==3*nrB && numel(valuesC)==3*nrC,'Vertical SSP shape mismatch');
        assert(isequal(valuesB(1:nrB-1),valuesC(1:nrC-1)) && isequal(valuesB(nrB+1:end),valuesC(nrC+1:end)),'Physical SSP coordinates or c changed');
        expectedB=(cfg.range_limit_m+max(8*cfg.step_m,1))/1000;
        expectedC=(cfg.range_limit_m+max(4*cfg.step_m,1))/1000;
        assert(valuesB(nrB)==expectedB && valuesC(nrC)==expectedC,'Original step-derived guard formula mismatch');
        assert(valuesB(nrB-1)==cfg.range_limit_m/1000 && valuesC(nrC-1)==cfg.range_limit_m/1000);
        assert(all(diff(valuesB(1:nrB))>0) && all(diff(valuesC(1:nrC))>0));
        assert(expectedB>cfg.range_limit_m/1000 && expectedC>cfg.range_limit_m/1000);
        receipt.numeric_guard_change=true; receipt.guard_B_m=1000*valuesB(nrB); receipt.guard_C_m=1000*valuesC(nrC); receipt.Beam_Box_r_m=cfg.range_limit_m;
        receipt.physical_SSP_and_all_sound_speeds_bitwise_equal=true; receipt.changed_SSP_indices=find(valuesB~=valuesC).';
        checks(end+1,:)={"SSP_physical_identical_and_original_guard_formula",true}; %#ok<AGROW>
    else
        assert(same,'B/C non-step physical input changed'); checks(end+1,:)={string(['byte_identical' ext{1}]),true}; %#ok<AGROW>
    end
end
oldlines=splitlines(string(fileread([reference '.env']))); newlines=splitlines(string(fileread([root '.env'])));
nonempty=find(strlength(strtrim(newlines))>0); assert(isequal(find(oldlines~=newlines),nonempty(end)),'Only ENV step row may differ');
oldlast=sscanf(char(oldlines(nonempty(end))),'%f'); newlast=sscanf(char(newlines(nonempty(end))),'%f');
assert(newlast(1)==oldlast(1)/2 && newlast(1)==cfg.step_m && isequal(newlast(2:end),oldlast(2:end)),'ENV step-only contract mismatch');
assert(newlast(3)==cfg.range_limit_m/1000,'Beam.Box.r changed');
checks(end+1,:)={'ENV_only_step_halved_Box_unchanged',true};
end
function [nr,values]=read_ssp(path)
fid=fopen(path,'r'); assert(fid>=0); c=onCleanup(@()fclose(fid)); nr=fscanf(fid,'%d',1); values=fscanf(fid,'%f');
end