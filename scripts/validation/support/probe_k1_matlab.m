function receipt=probe_k1_matlab(output_file)
%PROBE_K1_MATLAB Bounded environment/license query; no physics, pool or BLAS job.
arguments
    output_file (1,:) char
end
v=ver; names={v.Name};
receipt=struct('status','ENVIRONMENT_QUERY_COMPLETE','release',version('-release'),'version',version, ...
    'architecture',computer,'products',{names},'PCT_installed',any(strcmp(names,'Parallel Computing Toolbox')), ...
    'MATLAB_license_test',logical(license('test','MATLAB')), ...
    'PCT_license_test',logical(license('test','Distrib_Computing_Toolbox')), ...
    'default_max_computational_threads',maxNumCompThreads,'BLAS',version('-blas'), ...
    'pool_created',false,'acoustic_calls',0,'parallel_process_license_limit','unknown; no total-seat inference from license(test)');
receipt.PCT_checkout_attempted=false; receipt.PCT_checkout_succeeded=false;
if receipt.PCT_installed && receipt.PCT_license_test
    receipt.PCT_checkout_attempted=true;
    try
        receipt.PCT_checkout_succeeded=logical(license('checkout','Distrib_Computing_Toolbox'));
    catch err
        receipt.PCT_checkout_error_identifier=err.identifier;
    end
end
try
    inuse=license('inuse'); receipt.inuse_features=unique({inuse.feature});
catch err
    receipt.inuse_error_identifier=err.identifier;
end
fid=fopen(output_file,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(receipt,PrettyPrint=true));
fprintf('K1_NONACOUSTIC_DIAGNOSTIC=%s\n',jsonencode(receipt));
end
