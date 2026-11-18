function receipt = validate_d0_preparation(package_root,out)
%VALIDATE_D0_PREPARATION No physical source evaluations, synthetic fixtures only.
arguments
    package_root (1,:) char
    out (1,:) char
end
assert(~isfolder(out)); mkdir(out); t0=tic;
oldpath=path; guard=onCleanup(@()path(oldpath));
diary(fullfile(out,'matlab_diary.log')); dg=onCleanup(@()diary('off'));
addpath(fullfile(package_root,'code'));
sentinels=fullfile(out,'sentinels'); mkdir(sentinels);
blocked={'baffled_strip_piston_field_vertical','piston_independent_angular_spectrum_vertical', ...
    'd0_incident_spectrum_reference','run_piston_pm_r1_job'};
for j=1:numel(blocked)
    fid=fopen(fullfile(sentinels,[blocked{j} '.m']),'w'); assert(fid>=0);
    fprintf(fid,'function varargout=%s(varargin)\nerror(''D0:ForbiddenFieldCall'',''Formal field call in preparation'');\nend\n',blocked{j}); fclose(fid);
end
addpath(sentinels,'-begin');
for j=1:numel(blocked)
    assert(startsWith(strrep(which(blocked{j}),'\','/'),strrep(sentinels,'\','/'),'IgnoreCase',true));
end
preflight=run_d0_incident_depth(package_root,'preflight');
assert(preflight.field_calls==0 && ~isfolder(fullfile(package_root,'output')));
rejected=false;
try
    run_d0_incident_depth(package_root,'execute');
catch err
    rejected=contains(err.message,'Formal D0 requires') || contains(err.message,'owned remote supervisor');
end
assert(rejected && ~isfolder(fullfile(package_root,'output')));
cfg=jsondecode(fileread(fullfile(package_root,'config.json')));
testcfg=cfg; testcfg.kind='TEST_ONLY_SYNTHETIC'; testcfg.depths_m=[2 4 6 8 10];
sample_root=fullfile(out,'TEST_ONLY_SYNTHETIC'); mkdir(sample_root);
assertions=7; parseval_max=0; inverse_max=0;
for j=1:5
    source=struct('frequency_hz',25,'c_mps',1000,'rho_kgpm3',1000, ...
        'v0_mps',.0002,'total_width_m',.02,'z_source_m',testcfg.depths_m(j));
    x=(-220:.05:220)'; k=2*pi*source.frequency_hz/source.c_mps; q=k/3;
    amplitude=1+j/10; p=amplitude*exp(1i*(q*x+.3*j)); px=1i*q*p; pz=-1i*sqrt(k^2-q^2)*p;
    data=d0_incident_diagnostics(x,p,px,pz,source,testcfg);
    data.fixture='TEST_ONLY analytic plane wave, no piston field or formal result';
    assert(max(abs(abs(data.center_normalized_f)-1))<1e-12);
    assert(max(abs(abs(data.rms_normalized_f)-1))<1e-12);
    assert(max(abs(data.relative_phase_rad-q*x))<1e-12);
    assert(max(abs(data.local_angle_deg-asind(1/3)))<1e-10);
    assert(max(abs(data.phase_curvature_rad_per_m2(2:end-1)))<1e-10);
    % RMS normalization retains the given complex center phase.
    assert(abs(angle(data.rms_normalized_f(4401))-angle(p(4401)))<1e-12);
    for a=1:numel(data.window_spectra)
        w=data.window_spectra{a};
        parseval_max=max(parseval_max,w.parseval_relative); inverse_max=max(inverse_max,w.fft_inverse_relative);
        assert(w.parseval_relative<cfg.gates.fft_relative && w.fft_inverse_relative<cfg.gates.fft_relative);
    end
    folder=fullfile(sample_root,sprintf('depth_%04dm',source.z_source_m)); mkdir(folder);
    save(fullfile(folder,'field.mat'),'data','-v7'); writetable(data.points,fullfile(folder,'points.csv'));
    assertions=assertions+10;
end
figures=render_d0_incident_depth(sample_root,testcfg,'TEST_ONLY SYNTHETIC - NOT D0 RESULTS');
paths={'d0_incident_diagnostics.m','run_d0_incident_depth.m','d0_incident_spectrum_reference.m', ...
    'audit_d0_incident_case.m','render_d0_incident_depth.m'};
issues=cell(size(paths));
for j=1:numel(paths), issues{j}=checkcode(fullfile(package_root,'code',paths{j}),'-id'); end
pm_sensitivity_io('write',fullfile(out,'code_analyzer.json'),struct('paths',{paths},'issues',{issues}));
receipt=struct('status','PREPARATION_TEST_PASS','release',version('-release'),'PID',feature('getpid'), ...
    'elapsed_s',toc(t0),'assertions',assertions,'field_calls',0,'SSH_calls',0, ...
    'formal_mode_rejected',rejected,'FFT_Parseval_max',parseval_max,'FFT_inverse_max',inverse_max, ...
    'synthetic_figures',figures,'unrun',{{'formal field','deep AS numerical convergence','rough field gradients', ...
    'dx .025 physical checks','Ubuntu native launch/locks/signals','final PNG visual review after Fetch'}});
pm_sensitivity_io('write',fullfile(out,'receipt.json'),receipt); disp(jsonencode(receipt));
end
