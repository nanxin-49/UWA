function out=k1_pm_io(action,varargin)
%K1_PM_IO K1-only immutable inputs and audit-bound restart state machine.
switch action
    case 'sha'
        out=pm_sensitivity_io('sha',varargin{1});
    case {'read','write','path'}
        out=r1_bie_manual_io(action,varargin{:});
    case 'seal'
        root=varargin{1}; name=varargin{2}; signature=varargin{3}; paths=varargin{4}; status=varargin{5};
        assert(~isfile(fullfile(root,name)),'k1:ReplaceSeal','Never replace a seal');
        files=struct('path',{},'sha256',{});
        for j=1:numel(paths)
            p=k1_pm_io('path',root,paths{j});
            files(j)=struct('path',strrep(paths{j},'\','/'),'sha256',k1_pm_io('sha',p));
        end
        out=struct('schema_version',1,'status',status,'signature',signature,'files',files);
        k1_pm_io('write',fullfile(root,name),out);
    case 'verify_seal'
        root=varargin{1}; out=k1_pm_io('read',fullfile(root,varargin{2}));
        assert(out.schema_version==1 && strcmp(out.signature,varargin{3}) && ~isempty(out.files),'k1:StaleSeal','Inputs/code/predecessor changed');
        local_files(root,out.files);
    case 'resume_state'
        root=varargin{1}; signature=varargin{2};
        if isfile(fullfile(root,'committed.json'))
            seal=k1_pm_io('verify_seal',root,'committed.json',signature); assert(strcmp(seal.status,'AUDITED_PASS'));
            a=k1_pm_io('read',fullfile(root,'audit.json')); g=k1_pm_io('read',fullfile(root,'gates.json'));
            assert(strcmp(a.status,'PASS') && a.failed_checks==0 && strcmp(g.status,'PASS')); out='SKIP_AUDITED';
        elseif isfile(fullfile(root,'solved.json'))
            seal=k1_pm_io('verify_seal',root,'solved.json',signature); assert(strcmp(seal.status,'SOLVED_UNAUDITED')); out='AUDIT_ONLY';
        elseif isfile(fullfile(root,'solve_invocation.json')) || isfile(fullfile(root,'data.mat')) || isfile(fullfile(root,'points.csv'))
            error('k1:UnsealedInvocation','Unsealed invocation/data: inspect; do not repeat calculation');
        else
            out='SOLVE';
        end
    case 'verify_manifest'
        root=varargin{1}; out=k1_pm_io('read',fullfile(root,'manifest.json'));
        assert(strcmp(out.scope,'K1_PM_CASE') && out.seed==260001 && out.max_case_processes==1 && out.threads==1);
        assert(~out.parpool && ~out.automatic_retry && out.RAM_reserve_fraction==.3);
        local_files(root,out.files);
        frozen=k1_pm_io('read',fullfile(root,'input','frozen_config.json'));
        expected=struct('frequency_hz',4000,'c_mps',1500,'rho_kgpm3',1025,'v0_mps',.001,'total_width_m',.1,'z_source_m',100);
        assert(isequal(orderfields(frozen.source),orderfields(expected)) && frozen.seed==260001);
        matrix=[6 7;6 40;5 12;8 12;5 40;8 40];
        assert(ismember([frozen.U_PM_mps frozen.modes],matrix,'rows') && strcmp(frozen.case_id,out.case_id));
        assert(frozen.period_m==160 && frozen.C2_plateau_halfwidth_m==42 && frozen.C2_support_halfwidth_m==50);
        assert(strcmp(frozen.phase_contract_version,'2.0.0') && frozen.phase_mask_fraction_of_BIE_RMS==.01);
        assert(strcmp(frozen.family,'legacy12_gaussian_extension_v1') && strcmp(frozen.wind_label,'U_PM'));
        policy=k1_pm_io('read',fullfile(root,'policy.json')); r1_bie_continuation_io('policy',policy);
        assert(strcmp(k1_pm_io('sha',fullfile(root,'policy.json')),out.policy_sha256));
        p=readtable(fullfile(root,'input','coefficients.csv')); original=readtable(fullfile(root,'input','k0_coefficients.csv'));
        assert(height(p)==frozen.modes && isequal(p.k_rad_per_m,original.k_rad_per_m) && ...
            isequal(p.cos_coeff_m,original.cos_depth_coeff_m) && isequal(p.sin_coeff_m,original.sin_depth_coeff_m) && isequal(p.Sk_m3,original.Sk_m3));
        assert(all(isfinite(p{:,:}),'all') && max(abs(p.k_rad_per_m-(1:height(p))'*2*pi/160))<1e-14);
        sk=.0081./(2*p.k_rad_per_m.^3).*exp(-.74*9.81^2./(frozen.U_PM_mps^4*p.k_rad_per_m.^2));
        assert(max(abs(p.Sk_m3./sk-1))<1e-11);
        km=k1_pm_io('read',fullfile(root,'input','k0_metrics.json')); spectrum=readtable(fullfile(root,'input','k0_spectra.csv'));
        assert(km.N==frozen.modes && km.U_PM_mps==frozen.U_PM_mps && km.sample_count==8192 && km.dx_m==.01953125);
        assert(height(spectrum)==4097 && abs(sum(spectrum.raw_bin_power_m2)-km.raw_rms_height_m^2)<1e-12);
        assert(abs(sum(spectrum.embedded_bin_power_m2)-km.embedded_rms_height_m^2)<1e-12);
        assert(abs(sum(spectrum.embedded_bin_power_m2(frozen.modes+2:end))-km.embedded_power_above_Kmax_m2)<1e-12);
        assert(abs(sum(abs(spectrum.embedded_bin_power_m2-spectrum.raw_bin_power_m2))-km.embedding_spectral_L1_power_change_m2)<1e-12);
        r=readtable(fullfile(root,'input','receivers.csv'));
        assert(isequal(r.x_m,(-205:205)'*.1953125) && isequal(r.z_m,3*ones(411,1)));
        [z,slope,second]=closeout_s3a_pm_surface(p,linspace(-70,70,14001)');
        assert(all(isfinite([z;slope;second])) && min(z)>-2 && max(z)<3);
        proof=k1_pm_io('read',fullfile(root,'input','r1_anchor_proof.json'));
        assert(strcmp(proof.status,'STRICT_R1_K0_ANCHOR_MATCH_PASS') && proof.failed_checks==0);
    otherwise
        error('Unknown K1 IO action');
end
end

function local_files(root,files)
for j=1:numel(files)
    f=files(j); assert(strcmp(k1_pm_io('sha',k1_pm_io('path',root,f.path)),f.sha256),'k1:SHADrift','SHA drift: %s',f.path);
end
end
