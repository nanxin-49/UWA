function state=validate_deepsea3d_wp0_wp1(stage)
%VALIDATE_DEEPSEA3D_WP0_WP1 Checkpointed minimal 3D physics campaign.
% setup_vertical_project; stages wp0 / pilot / extend / reference / audit.
arguments
    stage (1,1) string = "wp0"
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_wp0_wp1');
if ~isfolder(out), mkdir(out); end
file=fullfile(out,'state.mat');
if isfile(file)
    a=load(file,'state'); state=a.state;
else
    state=struct('schema','1.0','matlab_version',version,'wp0_pass',false, ...
        'pilot_pass',false,'extended',false,'reference_done',false,'runs',struct([]), ...
        'wp0',struct([]),'guards',struct([]),'references',struct([]));
end
switch stage
    case "wp0"
        rows=struct([]);
        for q=[0 .1]
            for f=[4000 6000 8000 6037.5]
                c=local_cfg(q,f,"flat",.1,16); k=2*pi*f/1500;
                for h0=[-.2 0 .2]
                    b=deepsea3d_flat(c,h0); image=c.tx_m; image(3)=2*h0-c.tx_m(3);
                    direct=deepsea3d_point_field(c.rx_m,c.tx_m,k);
                    reflected=-deepsea3d_point_field(c.rx_m,image,k);
                    rot=[cos(.713) -sin(.713);sin(.713) cos(.713)];
                    cr=c; cr.rx_m(1:2)=c.rx_m(1:2)*rot.'; br=deepsea3d_flat(cr,h0);
                    e=max(abs([direct-b.direct_native reflected-b.reflect_native br.total_native-b.total_native]))/abs(reflected);
                    hd=-exp(-1i*2*pi*f*b.excess_delay_s)/(4*pi*b.Lr_m);
                    ep=abs(b.reflect_dsp-hd)/abs(reflected);
                    row=struct('q',q,'frequency_hz',f,'height_m',h0,'algebra_error',e,'phase_error',ep, ...
                        'Ld_m',b.Ld_m,'Lr_m',b.Lr_m,'delay_s',b.excess_delay_s, ...
                        'direct_real',real(direct),'direct_imag',imag(direct), ...
                        'reflect_real',real(reflected),'reflect_imag',imag(reflected), ...
                        'total_real',real(b.total_native),'total_imag',imag(b.total_native), ...
                        'total_condition',(abs(direct)+abs(reflected))/abs(b.total_native));
                    rows=local_append(rows,row);
                end
            end
        end
        state.wp0=rows; state.wp0_pass=all([rows.algebra_error]<1e-10 & [rows.phase_error]<1e-10);
        writetable(struct2table(rows),fullfile(out,'wp0.csv'));
    case "pilot"
        assert(state.wp0_pass,'deepsea3d:G0','WP0 must pass first.');
        state.runs=struct([]); state.guards=struct([]);
        state.extended=false; state.reference_done=false; state.references=struct([]);
        for field=["spectral" "cancellation" "completion_scope" "finite_roughness_G1_pass"]
            if isfield(state,field), state=rmfield(state,field); end
        end
        for q=[0 .1]
            for name=["flat" "weak" "dual"]
                vals=struct([]);
                for h=[.1 .05 .025]
                    c=local_cfg(q,6000,name,h,16); u=solve_deepsea3d_fullfk(c);
                    state.runs=local_run(state.runs,c,u,"grid"); vals=local_append(vals,u);
                    local_save(state,file,out);
                end
                denom=abs(vals(3).reflect_native);
                correction_scale=max(abs(vals(3).rough_correction),1e-12*denom);
                g=struct('q',q,'case_name',char(name),'frequency_hz',6000, ...
                    'grid_reflect_change',abs(vals(2).reflect_native-vals(3).reflect_native)/denom, ...
                    'grid_correction_change',abs(vals(2).rough_correction-vals(3).rough_correction)/correction_scale, ...
                    'grid_total_change',abs(vals(2).total_native-vals(3).total_native)/abs(vals(3).total_native), ...
                    'aperture_change',0,'rotation_change',0,'gauss_change',0);
                for L=[8 12 20]
                    c=local_cfg(q,6000,name,.025,L); u=solve_deepsea3d_fullfk(c);
                    state.runs=local_run(state.runs,c,u,"aperture");
                    if L>=12, g.aperture_change=max(g.aperture_change,abs(u.reflect_native-vals(3).reflect_native)/denom); end
                end
                c=local_cfg(q,6000,name,.05,16); c.quadrature='gauss4'; u=solve_deepsea3d_fullfk(c);
                state.runs=local_run(state.runs,c,u,"gauss");
                g.gauss_change=abs(u.reflect_native-vals(3).reflect_native)/denom;
                c=local_cfg(q,6000,name,.025,16); c=local_rotate(c);
                u=solve_deepsea3d_fullfk(c); state.runs=local_run(state.runs,c,u,"rotation90");
                g.rotation_change=abs(u.reflect_native-vals(3).reflect_native)/denom;
                state.guards=local_append(state.guards,g);
                local_save(state,file,out);
            end
        end
        % Embedding sensitivity changes the physical surface, never a convergence test.
        for support=[18 24]
            c=local_cfg(.1,6000,"dual",.025,support+4);
            c.profile.support_m=support; c.profile.plateau_m=2*support/3;
            u=solve_deepsea3d_fullfk(c); state.runs=local_run(state.runs,c,u,"embedding_sensitivity");
        end
        gs=state.guards;
        state.pilot_pass=all([gs.grid_reflect_change]<1e-3 & [gs.grid_correction_change]<1e-3 ...
            & [gs.aperture_change]<1e-3 & [gs.gauss_change]<1e-3 & [gs.rotation_change]<1e-10);
    case "extend"
        assert(state.pilot_pass,'deepsea3d:G1Numeric','6 kHz numerical pilot must pass.');
        state.runs=state.runs(~ismember(string({state.runs.purpose}),["extend_coarse" "extend_fine"]));
        state.guards=state.guards([state.guards.frequency_hz]==6000);
        for f=[4000 8000]
            for q=[0 .1]
                for name=["flat" "weak" "dual"]
                    c=local_cfg(q,f,name,.05,16); u=solve_deepsea3d_fullfk(c);
                    state.runs=local_run(state.runs,c,u,"extend_coarse");
                    c.spacing_m=.025; v=solve_deepsea3d_fullfk(c);
                    state.runs=local_run(state.runs,c,v,"extend_fine");
                    g=struct('q',q,'case_name',char(name),'frequency_hz',f, ...
                        'grid_reflect_change',abs(u.reflect_native-v.reflect_native)/abs(v.reflect_native), ...
                        'grid_correction_change',abs(u.rough_correction-v.rough_correction)/max(abs(v.rough_correction),1e-12*abs(v.reflect_native)), ...
                        'grid_total_change',abs(u.total_native-v.total_native)/abs(v.total_native), ...
                        'aperture_change',NaN,'rotation_change',NaN,'gauss_change',NaN);
                    state.guards=local_append(state.guards,g);
                    local_save(state,file,out);
                end
            end
        end
        state.extended=true;
    case "reference"
        assert(state.pilot_pass,'deepsea3d:G1Numeric','6 kHz numerical pilot must pass.');
        state.runs=state.runs(~ismember(string({state.runs.purpose}),["reference_fk" "reference_negative_fk"]));
        rows=struct([]);
        for q=[0 .1]
            for name=["weak" "dual"]
                for scale=[1 .5 .25]
                    c=local_cfg(q,6000,name,.025,16); c.profile.amplitudes_m=c.profile.amplitudes_m*scale;
                    fk=solve_deepsea3d_fullfk(c); state.runs=local_run(state.runs,c,fk,"reference_fk");
                    r=solve_deepsea3d_shape_reference(c); c.spacing_m=.05; rc=solve_deepsea3d_shape_reference(c);
                    c.spacing_m=.025; c.half_width_m=20; rw=solve_deepsea3d_shape_reference(c);
                    cp=c; cp.half_width_m=16; cp.profile.amplitudes_m=-cp.profile.amplitudes_m;
                    fm=solve_deepsea3d_fullfk(cp); state.runs=local_run(state.runs,cp,fm,"reference_negative_fk");
                    tangent=(fk.rough_correction-fm.rough_correction)/2;
                    row=struct('q',q,'case_name',char(name),'scale',scale, ...
                        'reference_real',real(r.rough_correction),'reference_imag',imag(r.rough_correction), ...
                        'fk_real',real(fk.rough_correction),'fk_imag',imag(fk.rough_correction), ...
                        'ref_grid_relative',abs(r.rough_correction-rc.rough_correction)/abs(r.rough_correction), ...
                        'ref_support_relative',abs(r.rough_correction-rw.rough_correction)/abs(r.rough_correction), ...
                        'tangent_relative',abs(tangent-r.rough_correction)/abs(r.rough_correction), ...
                        'perturbation_relative',abs(fk.rough_correction-r.rough_correction)/abs(r.rough_correction), ...
                        'reflected_relative',abs(fk.reflect_native-r.reflect_native)/abs(r.reflect_native), ...
                        'ref_wall_s',r.wall_s+rc.wall_s+rw.wall_s,'ref_cpu_s',r.cpu_s+rc.cpu_s+rw.cpu_s, ...
                        'ref_nodes',r.node_count);
                    rows=local_append(rows,row); state.references=rows;
                    local_save(state,file,out);
                end
            end
        end
        state.reference_done=true;
    case "rotation"
        for ii=1:numel(state.guards)
            g=state.guards(ii); if g.frequency_hz~=6000, continue; end
            c=local_cfg(g.q,6000,string(g.case_name),.025,16); c=local_rotate(c);
            u=solve_deepsea3d_fullfk(c);
            ix=find([state.runs.q]==g.q & [state.runs.frequency_hz]==6000 & ...
                strcmp({state.runs.case_name},g.case_name) & strcmp({state.runs.purpose},'grid') & ...
                [state.runs.spacing_m]==.025,1);
            baseline=complex(state.runs(ix).reflect_real,state.runs(ix).reflect_imag);
            state.guards(ii).rotation_change=abs(u.reflect_native-baseline)/abs(baseline);
            ir=find([state.runs.q]==g.q & strcmp({state.runs.case_name},g.case_name) & strcmp({state.runs.purpose},'rotation90'),1);
            rr=local_run(struct([]),c,u,"rotation90"); state.runs(ir)=rr;
        end
    case "cancellation"
        assert(state.pilot_pass,'deepsea3d:G1Numeric','Numerical pilot first.');
        state.runs=state.runs(~ismember(string({state.runs.purpose}),["cancellation_coarse" "cancellation_fine"]));
        rows=struct([]);
        for q=[0 .1]
            for name=["weak" "dual"]
                c=local_cfg(q,6000,name,.05,16); bg=deepsea3d_flat(c,0);
                f0=round(6000*bg.excess_delay_s)/bg.excess_delay_s;
                timer=tic; cpu0=cputime;
                [fmin,~,~,opt]=fminbnd(@(f)local_null_objective(c,f),f0-5,f0+5, ...
                    optimset('TolX',.001,'Display','off'));
                search_wall=toc(timer); search_cpu=cputime-cpu0;
                c.frequency_hz=fmin; a=solve_deepsea3d_fullfk(c);
                state.runs=local_run(state.runs,c,a,"cancellation_coarse");
                c.spacing_m=.025; b=solve_deepsea3d_fullfk(c);
                state.runs=local_run(state.runs,c,b,"cancellation_fine");
                row=struct('q',q,'case_name',char(name),'frequency_hz',fmin, ...
                    'total_abs',abs(b.total_native),'reflect_abs',abs(b.reflect_native), ...
                    'total_condition',(abs(b.direct_native)+abs(b.reflect_native))/abs(b.total_native), ...
                    'total_relative_change',abs(a.total_native-b.total_native)/abs(b.total_native), ...
                    'reflect_relative_change',abs(a.reflect_native-b.reflect_native)/abs(b.reflect_native), ...
                    'search_wall_s',search_wall,'search_cpu_s',search_cpu,'search_calls',opt.funcCount);
                rows=local_append(rows,row);
            end
        end
        state.cancellation=rows;
        writetable(struct2table(rows),fullfile(out,'cancellation.csv'));
    case "spectral"
        assert(state.reference_done,'deepsea3d:Reference','Real-space reference first.');
        c=local_cfg(0,6000,"weak",.025,16); r=solve_deepsea3d_shape_reference(c);
        rows=struct([]);
        for L=[32 48 64]
            a=deepsea3d_shape_spectral_check(c,L,.0625);
            row=struct('half_width_m',L,'spacing_m',a.spacing_m,'nodes',a.nodes, ...
                'spectral_real',real(a.rough_correction),'spectral_imag',imag(a.rough_correction), ...
                'relative_to_rayleigh',abs(a.rough_correction-r.rough_correction)/abs(r.rough_correction), ...
                'wall_s',a.wall_s,'cpu_s',a.cpu_s,'fft_array_MiB',a.named_fft_array_MiB);
            rows=local_append(rows,row);
        end
        state.spectral=rows;
        writetable(struct2table(rows),fullfile(out,'spectral_reference_check.csv'));
    case "audit"
        assert(state.wp0_pass && state.pilot_pass && state.extended && state.reference_done, ...
            'deepsea3d:Incomplete','Finish all scoped stages.');
        assert(all(isfinite([state.runs.reflect_real state.runs.reflect_imag])));
        assert(all([state.guards.grid_reflect_change]<1e-3 & [state.guards.grid_correction_change]<1e-3));
        assert(isfield(state,'spectral') && max([state.spectral.relative_to_rayleigh])<1e-5);
        assert(isfield(state,'cancellation') && all([state.cancellation.total_relative_change]<1e-3));
        assert(all([state.references.ref_grid_relative]<1e-5 & [state.references.ref_support_relative]<1e-5));
        assert(all([state.runs.identity_error]<1e-10));
        assert(numel(state.wp0)==24 && numel(state.guards)==18 && numel(state.references)==12 ...
            && numel(state.cancellation)==4 && numel(state.spectral)==3);
        for g=state.guards
            candidates=[state.runs.q]==g.q & [state.runs.frequency_hz]==g.frequency_hz & ...
                strcmp({state.runs.case_name},g.case_name);
            if g.frequency_hz==6000
                candidates=candidates & strcmp({state.runs.purpose},'grid');
            else
                candidates=candidates & ismember(string({state.runs.purpose}),["extend_coarse" "extend_fine"]);
            end
            coarse=state.runs(candidates & [state.runs.spacing_m]==.05);
            fine=state.runs(candidates & [state.runs.spacing_m]==.025);
            assert(isscalar(coarse) && isscalar(fine));
            pa=complex(coarse.reflect_real,coarse.reflect_imag); pb=complex(fine.reflect_real,fine.reflect_imag);
            ca=complex(coarse.correction_real,coarse.correction_imag); cb=complex(fine.correction_real,fine.correction_imag);
            assert(abs(abs(pa-pb)/abs(pb)-g.grid_reflect_change)<1e-13);
            assert(abs(abs(ca-cb)/max(abs(cb),1e-12*abs(pb))-g.grid_correction_change)<1e-13);
        end
        test_table=readtable(fullfile(out,'tests.csv'));
        assert(height(test_table)==8 && all(test_table.Passed) && ~any(test_table.Failed));
        state.completion_scope='WP0 plus numerical compact-surface 3D FK and first-order physical reference';
        state.finite_roughness_G1_pass=false;
        checks=struct('G0_pass',true,'numeric_6k_pass',true,'numeric_4_8k_pass',true, ...
            'shape_reference_numerics_pass',true,'spectral_crosscheck_pass',true, ...
            'cancellation_numerics_pass',true,'finite_roughness_G1_pass',false, ...
            'WP2_allowed',false,'run_count',numel(state.runs),'guard_count',numel(state.guards), ...
            'reference_comparisons',numel(state.references),'scope',state.completion_scope);
        fid=fopen(fullfile(out,'artifact_assertions.json'),'w');
        fprintf(fid,'%s',jsonencode(checks,PrettyPrint=true)); fclose(fid);
    otherwise
        error('deepsea3d:Stage','Unknown stage.');
end
local_save(state,file,out);
fprintf('DEEPSEA3D stage=%s G0=%d numeric6k=%d extended=%d reference=%d runs=%d\n', ...
    stage,state.wp0_pass,state.pilot_pass,state.extended,state.reference_done,numel(state.runs));
end
function v=local_null_objective(c,f)
c.frequency_hz=f; u=solve_deepsea3d_fullfk(c);
v=abs(u.total_native)/abs(u.direct_native);
end
function c=local_rotate(c)
c.rx_m(1:2)=[-c.rx_m(2) c.rx_m(1)];
c.center_m=[-c.center_m(2) c.center_m(1)];
c.profile.wavevectors_radpm=[-c.profile.wavevectors_radpm(:,2) c.profile.wavevectors_radpm(:,1)];
end
function c=local_cfg(q,f,name,h,L)
rho=q*1505; center=[rho*1500/1505 0];
p=struct('height_m',0,'plateau_m',8,'support_m',12,'amplitudes_m',0, ...
    'wavevectors_radpm',[.4 0],'phases_rad',.2);
if name=="weak", p.amplitudes_m=.001; end
if name=="dual", p.amplitudes_m=[.002 .001]; p.wavevectors_radpm=[.4 0;.2 .6]; p.phases_rad=[.2 .7]; end
c=struct('frequency_hz',f,'c0_mps',1500,'tx_m',[0 0 1500],'rx_m',[rho 0 5], ...
    'q',q,'case_name',char(name),'center_m',center,'half_width_m',L,'spacing_m',h, ...
    'quadrature','trap','profile',p);
end
function rows=local_run(rows,c,u,purpose)
row=struct('q',c.q,'frequency_hz',c.frequency_hz,'case_name',c.case_name,'purpose',char(purpose), ...
    'spacing_m',c.spacing_m,'half_width_m',c.half_width_m,'support_m',c.profile.support_m, ...
    'quadrature',c.quadrature,'nodes',u.node_count,'wall_s',u.wall_s,'cpu_s',u.cpu_s, ...
    'block_MiB_estimate',u.estimated_block_payload_MiB, ...
    'reflect_real',real(u.reflect_native),'reflect_imag',imag(u.reflect_native), ...
    'total_real',real(u.total_native),'total_imag',imag(u.total_native), ...
    'correction_real',real(u.rough_correction),'correction_imag',imag(u.rough_correction), ...
    'raw_flat_error',abs(u.raw_truncated_flat-u.flat_exact)/abs(u.flat_exact), ...
    'identity_error',u.quadrature_identity_residual/abs(u.flat_exact), ...
    'integral_cancellation_condition',u.cancellation_condition, ...
    'total_condition',(abs(u.direct_native)+abs(u.reflect_native))/abs(u.total_native));
rows=local_append(rows,row);
end
function rows=local_append(rows,row)
if isempty(rows), rows=row; else, rows(end+1)=row; end
end
function local_save(state,file,out)
save(file,'state','-v7');
if ~isempty(state.runs), writetable(struct2table(state.runs),fullfile(out,'runs.csv')); end
if ~isempty(state.guards), writetable(struct2table(state.guards),fullfile(out,'guards.csv')); end
if ~isempty(state.references), writetable(struct2table(state.references),fullfile(out,'shape_reference.csv')); end
end
