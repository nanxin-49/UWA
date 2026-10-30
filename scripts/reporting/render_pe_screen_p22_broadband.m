function audit = render_pe_screen_p22_broadband(run_id)
%RENDER_PE_SCREEN_P22_BROADBAND Reopen saved PE evidence and plot diagnostics.
arguments
    run_id (1,:) char
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outdir=fullfile(root,'results','validation','pe_screen_p22_comm',run_id);
assert(isfolder(outdir));
s=load(fullfile(outdir,'pe_screen_p22_broadband.mat'));
t=readtable(fullfile(outdir,'frequency_response.csv'));
td=readtable(fullfile(outdir,'cir_pdp_129.csv'));
assert(height(t)==258 && height(td)==256);
Hcsv=complex(t.H_total_real,t.H_total_imag);
Hdcsv=complex(t.H_direct_real,t.H_direct_imag);
Hrcsv=complex(t.H_reflect_real,t.H_reflect_imag);
hcsv=complex(td.cir_real,td.cir_imag);
csv_field_error=max(abs(Hcsv-s.H_total(:)));
csv_cir_error=max(abs(hcsv-s.h129(:)));
csv_pdp_error=max(abs(td.pdp-s.pdp129(:)));
component_error=max(abs(Hcsv-Hdcsv-Hrcsv));
fft_error=max([norm(fft(s.h65(:,1))-s.H_total(1:2:127,1))/norm(s.H_total(1:2:127,1)), ...
    norm(fft(s.h65(:,2))-s.H_total(1:2:127,2))/norm(s.H_total(1:2:127,2)), ...
    norm(fft(s.h129(:,1))-s.H_total(1:128,1))/norm(s.H_total(1:128,1)), ...
    norm(fft(s.h129(:,2))-s.H_total(1:128,2))/norm(s.H_total(1:128,2))]);
assert(csv_field_error<1e-13 && csv_cir_error<1e-13 && ...
    csv_pdp_error<1e-13 && component_error<1e-13 && fft_error<1e-10);
figdir=fullfile(outdir,'figures');
assert(~isfolder(figdir) && ~isfile(figdir)); mkdir(figdir);
f=figure('Visible','off','Color','w');
tiledlayout(2,1);
nexttile;
plot(s.f129,abs(s.H_total(:,1)),s.f129,abs(s.H_total(:,2)));
grid on; xlabel('Frequency (Hz)'); ylabel('|H_{total}|');
legend('Flat','Weak sine','Location','best');
title('P22 PE phase-screen complex frequency response');
nexttile;
plot(s.f129,unwrap(angle(s.H_total(:,1))), ...
    s.f129,unwrap(angle(s.H_total(:,2))));
grid on; xlabel('Frequency (Hz)'); ylabel('Unwrapped phase (rad)');
legend('Flat','Weak sine','Location','best');
savefig(f,fullfile(figdir,'frequency_response.fig'));
exportgraphics(f,fullfile(figdir,'frequency_response.png'),'Resolution',160);
close(f);
f=figure('Visible','off','Color','w');
tiledlayout(2,1);
for isurf=1:2
    nexttile;
    plot((0:63)*.25,s.pdp65(:,isurf),'o-', ...
        (0:127)*.25,s.pdp129(:,isurf),'-');
    grid on; xlabel('Delay after 97 m/c (ms)'); ylabel('PDP (relative)');
    if isurf==1, title('Flat'); else, title('Weak sine'); end
    legend('65 point / 16 ms','129 point / 32 ms','Location','best');
end
savefig(f,fullfile(figdir,'cir_pdp.fig'));
exportgraphics(f,fullfile(figdir,'cir_pdp.png'),'Resolution',160);
close(f);
audit=struct('status','SAVED_MAT_CSV_FIGURE_AUDIT_PASS', ...
    'csv_field_max_abs',csv_field_error,'csv_cir_max_abs',csv_cir_error, ...
    'csv_pdp_max_abs',csv_pdp_error,'component_max_abs',component_error, ...
    'fft_roundtrip_max_relative',fft_error,'figure_count',2, ...
    'cir_gate_status',s.receipt.status);
fid=fopen(fullfile(outdir,'audit.json'),'w'); assert(fid>0);
cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));
disp(audit);
end
