function receipt=migrate_r1_bie_checkpoints(root)
%MIGRATE_R1_BIE_CHECKPOINTS No density solve; C3 saved-density audit only.
arguments
    root (1,:) char
end
previous_threads=maxNumCompThreads(1); thread_guard=onCleanup(@()maxNumCompThreads(previous_threads));
r1_bie_continuation_io('verify_package',root);
origin_check=r1_bie_continuation_io('verify_origin',root);
if isfile(fullfile(root,'migration_complete.json'))
    r1_bie_continuation_io('verify_migration',root);
    receipt=struct('status','SKIP_MIGRATION_AUDITED','new_acoustic_solves',0); return;
end
assert(~isfile(fullfile(root,'migration_started.json')) && ~isfolder(fullfile(root,'jobs')), ...
    'r1v2:PartialMigration','Partial migration requires review; do not overwrite or recompute');
r1_bie_continuation_io('write',fullfile(root,'migration_started.json'),struct('status','STARTED','new_acoustic_solves',0));
origin=fullfile(root,'origin_evidence'); q=load(fullfile(root,'common.mat'));
s=q.cfg.source;
assert(s.total_width_m==.1 && s.frequency_hz==4000 && s.z_source_m==100 && s.c_mps==1500 && s.rho_kgpm3==1025 && s.v0_mps==.001);
assert(numel(q.x)==411 && all(q.z==3) && height(q.profile)==12 && strcmp(q.cfg.phase_contract_version,'2.0.0'));
assert(abs(max(q.profile.k_rad_per_m)-.471238898038469)<1e-14);
ids={'bie_coarse','bie_main','bie_self','bie_space'}; recipes=[60 10 96;60 12 96;60 12 64;60 14 96];
fields=cell(4,1); midpoint_counts=[340 408 408 476]; fp=r1_bie_continuation_io('sha',fullfile(root,'manifest.json'));
ph=r1_bie_continuation_io('sha',fullfile(root,'policy.json')); policy=r1_bie_continuation_io('read',fullfile(root,'policy.json'));
for j=1:4
    src=fullfile(origin,'jobs',ids{j}); dest=fullfile(root,'jobs',ids{j}); mkdir(dest);
    for name={'data.mat','points.csv','config.json','timings.json','audit.json','audit_reconstruction.mat'}
        copyfile(fullfile(src,name{1}),fullfile(dest,name{1}));
    end
    b=load(fullfile(dest,'data.mat')); fields{j}=b.field; bn=b.boundary;
    assert(bn.finite && all(bn.linear_residual<1e-10) && all(bn.offgrid_boundary_residual<1e-4));
    assert(isequal(bn.receiver_x_m,q.x) && isequal(bn.receiver_z_m,q.z));
    assert(bn.half_width_m==recipes(j,1) && bn.points_per_wavelength_requested==recipes(j,2) && bn.self_quadrature_order==recipes(j,3));
    assert(numel(bn.boundary_check_x_m)==midpoint_counts(j));
end
work=fullfile(root,'migration_audits','C3'); mkdir(work);
for name={'data.mat','points.csv'}, copyfile(fullfile(root,'jobs','bie_space',name{1}),fullfile(work,name{1})); end
timer=tic; diary(fullfile(work,'diary.log')); cleanup=onCleanup(@()diary('off'));
audit=audit_r1_bie_manual_stage(work,fullfile(root,'common.mat'),recipes(4,:));
elapsed=toc(timer); diary('off'); clear cleanup;
assert(strcmp(audit.status,'PASS') && audit.failed_checks==0);
for name={'audit.json','audit_reconstruction.mat'}
    copyfile(fullfile(work,name{1}),fullfile(root,'jobs','bie_space',name{1}));
end
gates=r1_bie_continuation_gates(fields,q.x,q.z,policy); gates.policy_sha256=ph;
assert(strcmp(gates.status,'PASS') && gates.grid_diagnostic.grid_nonmonotonic);
legacy=r1_bie_continuation_io('read',fullfile(origin,'jobs','bie_space','gates.json'));
assert(max(abs([gates.pairs.relative]-[legacy.pairs([1 3 2]).relative]))<1e-12);
pred=fp;
for j=1:4
    dest=fullfile(root,'jobs',ids{j});
    r1_bie_continuation_io('write',fullfile(dest,'gates_v2.json'),gates);
    r1_bie_continuation_io('write',fullfile(dest,'effective_config.json'), ...
        struct('stage',j-1,'manifest_sha256',fp,'policy_sha256',ph,'source_root','origin_evidence', ...
        'original_config_unchanged',true,'new_acoustic_solves',0,'audit_reexecuted',j==4,'legacy_policy_status',legacy.status));
    status='IMPORTED_AUDITED_PASS'; if j==4, status='IMPORTED_AUDITED_PASS_WITH_NONMONOTONIC_DIAGNOSTIC'; end
    sig=[pred ':C' num2str(j-1) ':' fp ':' ph];
    r1_bie_continuation_io('seal',dest,'migration_commit.json',sig, ...
        {'data.mat','points.csv','config.json','timings.json','audit.json','audit_reconstruction.mat','gates_v2.json','effective_config.json'},status);
    pred=r1_bie_continuation_io('sha',fullfile(dest,'migration_commit.json'));
end
receipt=struct('status','MIGRATION_AUDITED_PASS','origin',origin_check,'manifest_sha256',fp,'policy_sha256',ph, ...
    'completed_steps',{{'C0 imported','C1 imported','C2 imported','C3 saved-density independent audit PASS'}}, ...
    'C3_audit_s',elapsed,'max_threads',maxNumCompThreads,'historical_completed_solves',3,'migration_new_acoustic_solves',0, ...
    'legacy_policy_status','FAIL','grid_diagnostic',gates.grid_diagnostic, ...
    'next_allowed_step','Local acceptance tests then user-started C4 only; not R2');
r1_bie_continuation_io('write',fullfile(root,'migration_checkpoint.json'),receipt);
r1_bie_continuation_io('seal',root,'migration_complete.json',[pred ':MIGRATED:' fp ':' ph], ...
    {'migration_checkpoint.json','migration_audits/C3/audit.json','migration_audits/C3/audit_reconstruction.mat'},'MIGRATION_AUDITED_PASS');
r1_bie_continuation_io('verify_migration',root);
end
