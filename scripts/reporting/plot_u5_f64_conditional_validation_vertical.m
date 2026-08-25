function plot_u5_f64_conditional_validation_vertical()
%PLOT_U5_F64_CONDITIONAL_VALIDATION_VERTICAL Receiver-generator diagnostics.
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root); setup_vertical_project();
folder=fullfile(root,'results','validation','u5_conditional_channel_f64');
M=load(fullfile(folder,'u5_conditional_channel_model_f64_full.mat'),'model'); model=M.model;
E=load(fullfile(folder,'f64_full_ensembles.mat'),'ensemble'); e=E.ensemble;
[draw,~]=sample_conditional_channel_vertical(model,10000,2900001,struct('path','proper','rank','full'));
Hg=draw.H_ref_sca_f; Hk=e.Hkd_test; Hj=e.Hj_test; F=size(Hg,1); ii=ceil(F/2);

fig=figure('Visible','off','Color','w','Position',[50 50 1200 380]);
labels={'Real part','Imaginary part','Magnitude'};
for k=1:3
    subplot(1,3,k);
    if k==1, a=real(Hk(ii,:)); b=real(Hg(ii,:));
    elseif k==2, a=imag(Hk(ii,:)); b=imag(Hg(ii,:));
    else, a=abs(Hk(ii,:)); b=abs(Hg(ii,:)); end
    n=400; p=((1:n)-0.5)/n; qa=local_empirical_quantile(a,p); qb=local_empirical_quantile(b,p);
    plot(qa,qb,'LineWidth',1.2); hold on; lo=min([qa,qb]); hi=max([qa,qb]); plot([lo hi],[lo hi],'k--');
    axis equal; grid on; xlabel('held-out kdomain quantile'); ylabel('generator quantile'); title(labels{k});
end
exportgraphics(fig,fullfile(folder,'f64_generator_qq.png'),'Resolution',180); close(fig);

Sk=local_stats(Hk); Sj=local_stats(Hj); Sg=local_stats(Hg);
fig=figure('Visible','off','Color','w','Position',[50 50 1200 380]);
names={'kdomain','joint','generator'}; values={Sk.R,Sj.R,Sg.R};
for k=1:3, subplot(1,3,k); imagesc(abs(values{k}),[0,1]); axis image; colorbar; title([names{k} ' |R|']); xlabel('f_j'); ylabel('f_i'); end
exportgraphics(fig,fullfile(folder,'f64_generator_correlation_matrices.png'),'Resolution',180); close(fig);

fig=figure('Visible','off','Color','w'); semilogy(Sk.eig/max(Sk.eig),'-o'); hold on;
semilogy(Sj.eig/max(Sj.eig),'-s'); semilogy(Sg.eig/max(Sg.eig),'-^'); grid on;
xlabel('mode'); ylabel('normalized covariance eigenvalue'); legend(names,'Location','best');
exportgraphics(fig,fullfile(folder,'f64_generator_eigenvalues.png'),'Resolution',180); close(fig);

f=model.frequency_axis; ek=local_lfm(Hk,f); ej=local_lfm(Hj,f); eg=local_lfm(Hg,f);
fig=figure('Visible','off','Color','w'); plot(ek,'LineWidth',1.2); hold on; plot(ej,'LineWidth',1.2); plot(eg,'LineWidth',1.2);
grid on; xlabel('matched-filter sample'); ylabel('normalized envelope'); legend(names);
exportgraphics(fig,fullfile(folder,'f64_generator_lfm.png'),'Resolution',180); close(fig);
end

function q=local_empirical_quantile(x,p)
x=sort(x(:)); t=linspace(0,1,numel(x)); q=interp1(t,x,p,'linear','extrap');
end
function s=local_stats(H)
X=H-mean(H,2); C=(X*X')/(size(H,2)-1); d=sqrt(max(real(diag(C)),0));
s=struct('R',C./max(d*d.',eps),'eig',sort(max(real(eig(0.5*(C+C'))),0),'descend'));
end
function env=local_lfm(H,f)
fs=12000; N=512; t=(0:N-1).'/fs; D=.02; active=t<D; tx=zeros(N,1);
tx(active)=exp(1i*pi*((f(end)-f(1))/D)*(t(active)-D/2).^2);
fb=(-N/2:N/2-1).'*fs/N; Hb=complex(zeros(N,size(H,2)));
for m=1:size(H,2), Hb(:,m)=interp1(f-mean(f),H(:,m),fb,'linear',0); end
TX=fftshift(fft(tx)); rx=ifft(ifftshift(TX.*Hb),[],1); mf=ifft(fft(rx).*conj(fft(tx)),[],1);
env=mean(abs(mf),2); env=env/max(env);
end
