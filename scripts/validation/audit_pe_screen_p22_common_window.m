function audit = audit_pe_screen_p22_common_window(run_id,parent_run_id,refined_run_id)
%AUDIT_PE_SCREEN_P22_COMMON_WINDOW Independently reopen and recompute gate.
arguments
    run_id (1,:) char
    parent_run_id (1,:) char
    refined_run_id (1,:) char
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outdir=fullfile(root,'results','validation','pe_screen_p22_common_window',run_id);
p=load(fullfile(root,'results','validation','pe_screen_p22_comm', ...
    parent_run_id,'pe_screen_p22_broadband.mat'));
q=load(fullfile(root,'results','validation','pe_screen_p22_aliasing', ...
    refined_run_id,'pe_screen_p22_257.mat'));
s=load(fullfile(outdir,'common_window_cir.mat'));
t=readtable(fullfile(outdir,'common_window_cir.csv'));
assert(height(t)==512 && all(s.signed_t_ms==(-128:127)'*.25));
h65csv=complex(t.h65_real,t.h65_imag);
h129csv=complex(t.h129_real,t.h129_imag);
h257csv=complex(t.h257_real,t.h257_imag);
h65s=fftshift(s.h65,1); h129s=fftshift(s.h129,1);
h257s=fftshift(s.h257,1);
csv_error=max([max(abs(h65csv-h65s(:))), ...
    max(abs(h129csv-h129s(:))),max(abs(h257csv-h257s(:)))]);
assert(csv_error<1e-13);
tau=6/1500;
error65_129=zeros(1,2); error65_257=error65_129;
error129_257=error65_129; energy_outside=error65_129;
mask=s.signed_t_ms>=-4 & s.signed_t_ms<12;
for j=1:2
    f65=p.f65; f129=p.f129; f257=q.f257;
    d65=interp1(f65,p.H_direct(1:2:end,j),f257,'linear');
    r65=interp1(f65,p.H_reflect(1:2:end,j).* ...
        exp(1i*2*pi*f65*tau),f257,'linear').*exp(-1i*2*pi*f257*tau);
    d129=interp1(f129,p.H_direct(:,j),f257,'linear');
    r129=interp1(f129,p.H_reflect(:,j).* ...
        exp(1i*2*pi*f129*tau),f257,'linear').*exp(-1i*2*pi*f257*tau);
    ha=ifft(d65(1:256)+r65(1:256));
    hb=ifft(d129(1:256)+r129(1:256));
    hc=ifft(q.H_direct(1:256,j)+q.H_reflect(1:256,j));
    assert(norm(ha-s.h65(:,j))/norm(ha)<1e-12);
    assert(norm(hb-s.h129(:,j))/norm(hb)<1e-12);
    assert(norm(hc-s.h257(:,j))/norm(hc)<1e-12);
    error65_129(j)=norm(ha-hb)/norm(hb);
    error65_257(j)=norm(ha-hc)/norm(hc);
    error129_257(j)=norm(hb-hc)/norm(hc);
    signed=fftshift(hc);
    energy_outside(j)=sum(abs(signed(~mask)).^2)/sum(abs(signed).^2);
end
assert(max([error65_129,error65_257,error129_257])<=1e-3);
assert(max(energy_outside)<=1e-3);
assert(~s.receipt.original_raw_gate_pass);
audit=struct('status','PE_SCREEN_COMMON_WINDOW_CIR_INDEPENDENT_AUDIT_PASS', ...
    'csv_cir_max_abs',csv_error,'cir65_129_relative',error65_129, ...
    'cir65_257_relative',error65_257,'cir129_257_relative',error129_257, ...
    'outside_signed_window_energy',energy_outside, ...
    'unchanged_threshold',1e-3, ...
    'original_raw_gate_pass',false,'common_window_gate_pass',true);
figdir=fullfile(outdir,'figures');
assert(~isfolder(figdir)); mkdir(figdir);
f=figure('Visible','off','Color','w');
tiledlayout(2,1);
for j=1:2
    nexttile;
    semilogy(s.signed_t_ms,max(abs(h65s(:,j)).^2,1e-13),'o', ...
        s.signed_t_ms,max(abs(h129s(:,j)).^2,1e-13),'-', ...
        s.signed_t_ms,max(abs(h257s(:,j)).^2,1e-13),'--');
    xlim([-8,12]); grid on;
    xlabel('Delay relative to 97 m/c (ms)'); ylabel('PDP (relative)');
    if j==1, title('Flat, common 64 ms window');
    else, title('Weak sine, common 64 ms window'); end
    legend('65 predicted','129 predicted','257 PE','Location','best');
end
savefig(f,fullfile(figdir,'common_window_pdp.fig'));
exportgraphics(f,fullfile(figdir,'common_window_pdp.png'),'Resolution',160);
close(f);
fid=fopen(fullfile(outdir,'audit.json'),'w'); assert(fid>0);
cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));
disp(audit);
end
