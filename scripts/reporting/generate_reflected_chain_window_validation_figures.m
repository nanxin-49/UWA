function files = generate_reflected_chain_window_validation_figures(input_value)
%GENERATE_REFLECTED_CHAIN_WINDOW_VALIDATION_FIGURES Plot reflected-chain audit.
if ischar(input_value) || isstring(input_value)
    data=load(input_value,'validation'); v=data.validation;
else
    v=input_value;
end
out=fullfile(v.config.output_dir,'figures');
if ~exist(out,'dir'), mkdir(out); end
files=strings(0,1);

f=figure('Visible','off','Color','w','Position',[100 100 1200 800]);
t=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ref=v.cases(end);
nexttile; imagesc(ref.x_m,ref.y_m,ref.surface_elevation_xy); axis image xy; colorbar;
title(sprintf('Fixed surface, W=%.1f m, H_s=%.3f m',ref.width_actual_m,ref.crop_hs_m)); xlabel('x (m)'); ylabel('y (m)');
nexttile; imagesc(ref.x_m,ref.y_m,20*log10(abs(ref.surface_incident_xy)/max(abs(ref.surface_incident_xy(:)))+eps)); axis image xy; colorbar; clim([-80 0]);
title('Surface incident field (dB normalized)'); xlabel('x (m)'); ylabel('y (m)');
nexttile; imagesc(ref.x_m,ref.y_m,20*log10(abs(ref.surface_reflected_xy)/max(abs(ref.surface_reflected_xy(:)))+eps)); axis image xy; colorbar; clim([-80 0]);
title('Surface reflected field (dB normalized)'); xlabel('x (m)'); ylabel('y (m)');
nexttile; imagesc(ref.x_m,ref.y_m,20*log10(abs(ref.receiver_reflected_xy)/max(abs(ref.receiver_reflected_xy(:)))+eps)); axis image xy; colorbar; clim([-80 0]);
title('Downward field at receiver depth (dB normalized)'); xlabel('x (m)'); ylabel('y (m)');
title(t,'Full reflected-chain fields');
files(end+1)=local_save(f,out,'01_reflected_chain_fields.png'); close(f);

f=figure('Visible','off','Color','w','Position',[100 100 1150 720]);
t=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
stages={'surface_incident_xy','surface_reflected_xy','receiver_reflected_xy'};
labels={'incident at surface','reflected at surface','reflected at receiver'};
for ss=1:3
    nexttile; hold on; grid on;
    for ii=1:numel(v.cases)
        c=v.cases(ii); iy=round(numel(c.y_m)/2)+1;
        plot(c.x_m,20*log10(abs(c.(stages{ss})(iy,:))/max(abs(c.(stages{ss})(:)))+eps), ...
            'DisplayName',sprintf('W=%.1f m',c.width_actual_m));
    end
    xlim([-40 40]); ylim([-100 5]); xlabel('x (m)'); ylabel('normalized magnitude (dB)'); title(labels{ss});
end
nexttile; hold on; grid on;
plot(v.case_table.width_actual_m,v.case_table.reflect_tl_db,'o-','DisplayName','reflected TL');
ylabel('TL (dB)'); yyaxis right; plot(v.case_table.width_actual_m,v.case_table.reflect_phase_rad,'s-','DisplayName','phase'); ylabel('phase (rad)'); xlabel('actual W (m)');
title('Receiver scalar convergence');
title(t,'Window convergence of reflected chain');
files(end+1)=local_save(f,out,'02_window_convergence_fields_and_receiver.png'); close(f);

f=figure('Visible','off','Color','w','Position',[100 100 1100 720]);
t=tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on;
for ii=1:numel(v.cases)
    tr=v.cases(ii).incident_energy_trace;
    semilogy(v.config.z_tx_m-tr.z_m,max(tr.edge5_fraction,realmin),'DisplayName',sprintf('up W=%.1f',v.cases(ii).width_actual_m));
end
xlabel('upward propagation distance (m)'); ylabel('outer 5% energy fraction'); title('Upward PE edge energy'); legend('Location','best');
nexttile; hold on; grid on;
for ii=1:numel(v.cases)
    tr=v.cases(ii).reflected_energy_trace;
    semilogy(tr.z_m,max(tr.edge5_fraction,realmin),'DisplayName',sprintf('down W=%.1f',v.cases(ii).width_actual_m));
end
xlabel('depth below surface (m)'); ylabel('outer 5% energy fraction'); title('Downward reflected PE edge energy'); legend('Location','best');
title(t,'Stage-localized boundary diagnostics');
files(end+1)=local_save(f,out,'03_stage_edge_energy.png'); close(f);

if v.converged && ~isempty(v.metric_table)
    rows=v.metric_table(v.metric_table.radius_m==2,:);
    f=figure('Visible','off','Color','w','Position',[100 100 1100 650]);
    t=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    nexttile; plot(rows.width_actual_m,abs(rows.reflect_tl_error_db),'o-'); yline(v.config.tl_limit_db,'r--'); grid on; xlabel('W (m)'); ylabel('|Delta TL| (dB)');
    nexttile; plot(rows.width_actual_m,abs(rows.reflect_phase_error_rad),'o-'); yline(v.config.phase_limit_rad,'r--'); grid on; xlabel('W (m)'); ylabel('|Delta phase| (rad)');
    nexttile; semilogy(rows.width_actual_m,rows.surface_incident_complex_l2,'o-',rows.width_actual_m,rows.surface_reflected_complex_l2,'s-',rows.width_actual_m,rows.receiver_reflected_complex_l2,'^-'); yline(v.config.center_l2_limit,'r--'); grid on; xlabel('W (m)'); ylabel('center complex L2'); legend('incident','surface reflected','receiver','limit');
    nexttile; semilogy(v.case_table.width_actual_m,max(v.case_table.max_edge5_fraction,realmin),'o-'); yline(v.config.edge5_fraction_limit,'r--'); grid on; xlabel('W (m)'); ylabel('max outer-5% fraction');
    title(t,'Preregistered convergence gates (rho <= 2 m)');
    files(end+1)=local_save(f,out,'04_convergence_gates.png'); close(f);
end
end

function path=local_save(f,out,name)
path=fullfile(out,name); exportgraphics(f,path,'Resolution',180);
end
