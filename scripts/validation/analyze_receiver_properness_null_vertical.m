function properness=analyze_receiver_properness_null_vertical(mat_path,n_mc,seed)
%ANALYZE_RECEIVER_PROPERNESS_NULL_VERTICAL Calibrate finite-sample P/C bias.
if nargin<1 || isempty(mat_path)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    mat_path=fullfile(root,'results','validation','cached_joint_kstat_pe_receiver', ...
        'cached_joint_kstat_pe_receiver_validation.mat');
end
if nargin<2, n_mc=2000; end
if nargin<3, seed=900001; end
S=load(mat_path,'validation'); v=S.validation; modes={'kdomain','joint','independent'};
rng(seed,'twister'); properness=struct('n_mc',n_mc,'seed',seed, ...
    'null_hypothesis','proper complex Gaussian with train-set covariance and test-set sample count');
null_by_mode=cell(size(modes));
for kk=1:numel(modes)
    mode=modes{kk}; C=v.stats.([mode '_train']).C; L=size(v.H_ref_sca.([mode '_test']),2);
    [Q,D]=eig(0.5*(C+C'),'vector'); keep=real(D)>max(real(D))*1e-12;
    A=Q(:,keep).*sqrt(max(real(D(keep)).',0)); null_ratio=zeros(n_mc,1);
    for bb=1:n_mc
        Z=(randn(sum(keep),L)+1i*randn(sum(keep),L))/sqrt(2);
        X=A*Z; X=X-mean(X,2); Cs=(X*X')/(L-1); Ps=(X*X.')/(L-1);
        null_ratio(bb)=norm(Ps,'fro')/max(norm(Cs,'fro'),eps);
    end
    q=sort(null_ratio); null_by_mode{kk}=q;
    ci90=local_quantiles(q,[0.05,0.95]);
    ci95=local_quantiles(q,[0.025,0.975]);
    ci99=local_quantiles(q,[0.005,0.995]);
    observed=v.stats.([mode '_test']).pseudo_ratio;
    properness.(mode)=struct('observed_ratio',observed,'proper_null_mean',mean(q), ...
        'proper_null_median',median(q),'proper_null_ci90',ci90, ...
        'proper_null_ci95',ci95,'proper_null_ci99',ci99, ...
        'upper_tail_p',(1+sum(q>=observed))/(n_mc+1), ...
        'observed_inside_ci95',observed>=ci95(1) && observed<=ci95(2), ...
        'reject_proper_at_5pct',observed>ci95(2), ...
        'reject_proper_at_1pct',observed>ci99(2), ...
        'decision',local_decision(observed,ci95,ci99), ...
        'test_sample_count',L,'null_ratio_samples',q);
end
[folder,~,~]=fileparts(mat_path);
save(fullfile(folder,'receiver_properness_null.mat'),'properness');
local_plot_null(null_by_mode,modes,properness,folder);
fprintf('P/C observed | null mean median | CI90 | CI95 | CI99 | p_upper | decision\n');
for kk=1:numel(modes)
    x=properness.(modes{kk}); fprintf(['%s %.6f %.6f %.6f ', ...
        '[%.6f %.6f] [%.6f %.6f] [%.6f %.6f] %.6f %s\n'], ...
        modes{kk},x.observed_ratio,x.proper_null_mean,x.proper_null_median, ...
        x.proper_null_ci90,x.proper_null_ci95,x.proper_null_ci99, ...
        x.upper_tail_p,x.decision);
end
end

function values=local_quantiles(sorted_values,p)
n=numel(sorted_values); index=max(1,min(n,round(1+(n-1).*p)));
values=sorted_values(index).';
end

function value=local_decision(observed,ci95,ci99)
if observed>ci99(2)
    value='proper hypothesis rejected above 99% upper bound';
elseif observed>ci95(2)
    value='proper hypothesis rejected at 5% upper-tail level';
else
    value='proper hypothesis not rejected';
end
end

function local_plot_null(null_by_mode,modes,properness,folder)
colors=lines(numel(modes));
fig=figure('Visible','off','Color','w','Position',[80 80 1100 420]);
subplot(1,2,1); hold on;
h_pdf=gobjects(numel(modes),1); h_obs=gobjects(numel(modes),1);
for kk=1:numel(modes)
    h_pdf(kk)=histogram(null_by_mode{kk},40,'Normalization','pdf','DisplayStyle','stairs', ...
        'LineWidth',1.4,'EdgeColor',colors(kk,:));
    h_obs(kk)=xline(properness.(modes{kk}).observed_ratio,'--','Color',colors(kk,:), ...
        'LineWidth',1.2);
end
grid on; xlabel('||P||_F / ||C||_F'); ylabel('null PDF'); title('Proper-null histogram');
subplot(1,2,2); hold on;
for kk=1:numel(modes)
    q=null_by_mode{kk}; plot(q,(1:numel(q)).'/numel(q),'LineWidth',1.4,'Color',colors(kk,:));
    xline(properness.(modes{kk}).observed_ratio,'--','Color',colors(kk,:),'LineWidth',1.2);
end
grid on; xlabel('||P||_F / ||C||_F'); ylabel('empirical CDF'); title('Proper-null empirical CDF');
legend([h_pdf;h_obs],[strcat(modes,' null'),strcat(modes,' observed')], ...
    'Location','bestoutside');
exportgraphics(fig,fullfile(folder,'receiver_properness_null.png'),'Resolution',180); close(fig);
end
