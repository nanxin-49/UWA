function audit = audit_pe_screen_p22_aliasing(run_id,parent_run_id)
%AUDIT_PE_SCREEN_P22_ALIASING Recompute saved alias and quadrature evidence.
arguments
    run_id (1,:) char
    parent_run_id (1,:) char
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outdir=fullfile(root,'results','validation','pe_screen_p22_aliasing',run_id);
parent=fullfile(root,'results','validation','pe_screen_p22_comm',parent_run_id);
q=load(fullfile(outdir,'pe_screen_p22_257.mat'));
p=load(fullfile(parent,'pe_screen_p22_broadband.mat'));
t=readtable(fullfile(outdir,'frequency_response_257.csv'));
td=readtable(fullfile(outdir,'cir_pdp_257.csv'));
assert(height(t)==514 && height(td)==512);
Hcsv=complex(t.H_total_real,t.H_total_imag);
hcsv=complex(td.cir_real,td.cir_imag);
csv_H=max(abs(Hcsv-q.H_total(:)));
csv_h=max(abs(hcsv-q.h257(:)));
csv_pdp=max(abs(td.pdp-abs(q.h257(:)).^2));
assert(csv_H<1e-13 && csv_h<1e-13 && csv_pdp<1e-13);
assert(isequal(q.f257(1:2:end),p.f129));
assert(max(abs(q.H_total(1:2:end,:)-p.H_total),[],'all')==0);
window_ms=[-4,12];
time_s=(window_ms(1):.25:window_ms(2))'*.001;
fields={p.H_direct,p.H_reflect,p.H_total};
fields257={q.H_direct,q.H_reflect,q.H_total};
names=["direct";"reflect";"total"];
rows=zeros(6,7);
for k=1:3
    H129=fields{k}; H257=fields257{k};
    for j=1:2
        H65=H129(1:2:end,j);
        h65=ifft(H65(1:64));
        h129=ifft(H129(1:128,j));
        h257=ifft(H257(1:256,j));
        alias1=norm(h65-h129(1:64)-h129(65:128))/norm(h129);
        alias2=norm(h129-h257(1:128)-h257(129:256))/norm(h257);
        c65=local_band_integral(p.f65,H65,time_s);
        c129=local_band_integral(p.f129,H129(:,j),time_s);
        c257=local_band_integral(q.f257,H257(:,j),time_s);
        r=(k-1)*2+j;
        rows(r,:)=[alias1,alias2, ...
            norm(c65-c129)/norm(c129),norm(c129-c257)/norm(c257), ...
            sum(abs(h257(65:128)).^2)/sum(abs(h257).^2), ...
            sum(abs(h257(129:256)).^2)/sum(abs(h257).^2), ...
            abs(H129(end,j)-H129(1,j))];
    end
end
surface=repmat(["flat";"weak_sine"],3,1);
component=repelem(names,2);
metrics=table(component,surface,rows(:,1),rows(:,2),rows(:,3), ...
    rows(:,4),rows(:,5),rows(:,6),rows(:,7), ...
    'VariableNames',{'component','surface','alias65_129_relative', ...
    'alias129_257_relative','band_integral65_129_relative', ...
    'band_integral129_257_relative','tail16_32_energy_fraction', ...
    'tail32_64_energy_fraction','endpoint_complex_absolute'});
writetable(metrics,fullfile(outdir,'audit_metrics_v2.csv'));
total=rows(5:6,:);
tau_reflect_s=6/1500;
frequency_prediction=zeros(2,2);
for j=1:2
    d65=p.H_direct(1:2:end,j);
    r65=p.H_reflect(1:2:end,j).*exp(1i*2*pi*p.f65*tau_reflect_s);
    predicted129=interp1(p.f65,d65,p.f129,'linear')+ ...
        interp1(p.f65,r65,p.f129,'linear').* ...
        exp(-1i*2*pi*p.f129*tau_reflect_s);
    frequency_prediction(j,1)=norm(predicted129-p.H_total(:,j))/norm(p.H_total(:,j));
    r129=p.H_reflect(:,j).*exp(1i*2*pi*p.f129*tau_reflect_s);
    predicted257=interp1(p.f129,p.H_direct(:,j),q.f257,'linear')+ ...
        interp1(p.f129,r129,q.f257,'linear').* ...
        exp(-1i*2*pi*q.f257*tau_reflect_s);
    frequency_prediction(j,2)=norm(predicted257-q.H_total(:,j))/norm(q.H_total(:,j));
end
last4=zeros(1,2);
for j=1:2
    last4(j)=sum(abs(q.h257(241:256,j)).^2)/sum(abs(q.h257(:,j)).^2);
end
audit=struct('status','ALIASING_CAUSE_AUDITED_CIR_GATE_FAILED', ...
    'csv_H_max_abs',csv_H,'csv_cir_max_abs',csv_h, ...
    'csv_pdp_max_abs',csv_pdp,'delay_window_ms',window_ms, ...
    'delay_sample_spacing_ms',.25, ...
    'max_alias_identity_relative',max(rows(:,1:2),[],'all'), ...
    'total_band_integral65_129_relative',total(:,3).', ...
    'total_band_integral129_257_relative',total(:,4).', ...
    'total_tail16_32_energy_fraction',total(:,5).', ...
    'total_tail32_64_energy_fraction',total(:,6).', ...
    'total_last4ms_energy_fraction',last4, ...
    'delay_aware_frequency_prediction_relative',frequency_prediction, ...
    'delay_aware_frequency_gate_pass',max(frequency_prediction,[],'all')<=1e-3, ...
    'original_65_129_gate_pass',max(p.receipt.cir_refinement_relative)<=1e-3, ...
    'nonperiodic_65_129_gate_pass',max(total(:,3))<=1e-3, ...
    'nonperiodic_129_257_gate_pass',max(total(:,4))<=1e-3, ...
    'bie_calls',0,'fullfk_calls',0,'three_dimensional_calls',0);
figdir=fullfile(outdir,'figures_v2');
assert(~isfolder(figdir)); mkdir(figdir);
f=figure('Visible','off','Color','w');
tiledlayout(2,1);
for j=1:2
    nexttile;
    tms=(-128:127)*.25;
    h_signed=[q.h257(129:256,j);q.h257(1:128,j)];
    semilogy(tms,max(abs(h_signed).^2,1e-13),'.-');
    grid on; xlim([-16,16]); xlabel('Delay relative to 97 m/c (ms)');
    ylabel('PDP (relative)');
    if j==1, title('Flat, 257 point signed delay');
    else, title('Weak sine, 257 point signed delay'); end
end
savefig(f,fullfile(figdir,'signed_delay_257.fig'));
exportgraphics(f,fullfile(figdir,'signed_delay_257.png'),'Resolution',160);
close(f);
fid=fopen(fullfile(outdir,'audit_v2.json'),'w'); assert(fid>0);
cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));
disp(audit);
end

function h=local_band_integral(f,H,t)
weights=ones(numel(f),1);
weights([1,end])=.5;
band=f(end)-f(1);
h=exp(1i*2*pi*t*(f(:).'-f(1)))*(weights.*H(:))*(f(2)-f(1))/band;
end
