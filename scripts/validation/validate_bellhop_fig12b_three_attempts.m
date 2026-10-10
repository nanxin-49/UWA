function output = validate_bellhop_fig12b_three_attempts(stage)
%VALIDATE_BELLHOP_FIG12B_THREE_ATTEMPTS Bounded three native beam/source tests.
% Each attempt includes 50/500/2000 beams, a step check and a matched control.
% No additional search is initiated by this workflow after its fixed plan.
arguments
    stage (1,:) char = 'status'
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();
outdir=fullfile(root,'results','validation','bellhop_fig12b_three_attempts');
if ~isfolder(outdir), mkdir(outdir); end
template=struct('output_dir',outdir,'tag','','beam_count',2000, ...
    'angles_deg',[-179.91 179.91],'run_type','CB R','step_m',1.25, ...
    'rbox_m',100,'matched_control',false,'depth_m',linspace(0,5000,1001), ...
    'toolbox_root','E:/MISC/BELLHOP/AcousticsToolbox_2020', ...
    'range_m',linspace(-5,5,401),'water_depth_m',5000,'source_depth_m',4950, ...
    'frequency_hz',10000,'munk_c0_mps',1500,'munk_axis_m',1300, ...
    'munk_scale_m',1300,'munk_epsilon',.00737,'ssp_spacing_m',10, ...
    'bottom_cp_mps',1600,'bottom_rho_gcm3',1.8,'bottom_atten_db_wavelength',.8, ...
    'uniform_control',false,'cerveny_type','MS','epsilon_multiplier',1, ...
    'loop_range_m',5000,'image_count',1,'beam_window',5);
plan=repmat(template,1,15); k=0;
types={'CB R','CR R','CB X'};
for attempt=1:3
    for beams=[50 500 2000]
        k=k+1; p=template; p.run_type=types{attempt}; p.beam_count=beams;
        delta=360/beams;
        if mod(beams,4)==0, p.angles_deg=[-180+delta/2 180-delta/2];
        else, p.angles_deg=[-180 180-delta]; end
        p.tag=sprintf('a%d_n%d_h1.25',attempt,beams); plan(k)=p;
    end
    k=k+1; p=template; p.run_type=types{attempt}; p.step_m=.625;
    p.tag=sprintf('a%d_n2000_h0.625',attempt); plan(k)=p;
    k=k+1; p=template; p.run_type=types{attempt}; p.matched_control=true;
    p.depth_m=[1000 3000 4900]; p.step_m=.625;
    p.tag=sprintf('a%d_control',attempt); plan(k)=p;
end
save(fullfile(outdir,'attempt_plan.mat'),'plan','types');
switch lower(stage)
    case {'next','full'}
        output=[];
        for ii=1:numel(plan)
            target=fullfile(outdir,[plan(ii).tag '.mat']);
            if isfile(target)
                s=load(target,'output'); old=s.output.cfg; fields=fieldnames(plan(ii));
                for jj=1:numel(fields)
                    assert(isequaln(old.(fields{jj}),plan(ii).(fields{jj})), ...
                        'Saved case differs: %s',target);
                end
                continue
            end
            fprintf('Attempt %d of 3: %s\n',ceil(ii/5),plan(ii).tag);
            output=validate_bellhop_fig12b('case',plan(ii));
            if strcmpi(stage,'next'), return; end
        end
        if strcmpi(stage,'full'), output=generate_bellhop_fig12b_three_attempt_visuals(outdir); end
    case 'render'
        output=generate_bellhop_fig12b_three_attempt_visuals(outdir);
    case 'status'
        output=plan;
        for ii=1:numel(plan)
            fprintf('%-25s %d\n',plan(ii).tag,isfile(fullfile(outdir,[plan(ii).tag '.mat'])));
        end
    otherwise
        error('stage must be status, next, full or render.');
end
end
