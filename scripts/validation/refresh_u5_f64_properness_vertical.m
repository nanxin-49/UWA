function properness=refresh_u5_f64_properness_vertical()
%REFRESH_U5_F64_PROPERNESS_VERTICAL Refresh prescribed-rule F=64 null test.
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root);
folder=fullfile(root,'results','validation','u5_conditional_channel_f64');
rf=fullfile(folder,'u5_conditional_channel_validation_f64_full.mat');
ef=fullfile(folder,'f64_full_ensembles.mat'); mf=fullfile(folder,'u5_conditional_channel_model_f64_full.mat');
R=load(rf,'result'); E=load(ef,'ensemble'); result=R.result; e=E.ensemble;
properness=struct();
properness.kdomain=properness_null_test_vertical(result.stats.kdomain_train.C,e.Hkd_test,2000,result.seeds.properness_base+1);
properness.joint=properness_null_test_vertical(result.stats.joint_train.C,e.Hj_test,2000,result.seeds.properness_base+2);
properness.independent=properness_null_test_vertical(result.stats.independent_train.C,e.Hi_test,2000,result.seeds.properness_base+3);
result.properness=properness; result.pass.proper_joint_not_rejected=~properness.joint.reject_proper_at_5pct;
save(rf,'result','-v7.3');
M=load(mf,'model'); model=M.model; model.stats.properness_test=properness.joint;
model.validation.properness=properness; save(mf,'model','-v7.3');
local_plot(properness,folder);
for name={'kdomain','joint','independent'}
    q=properness.(name{1}); fprintf('%s observed %.6f CI95 [%.6f %.6f] p %.6f %s; one-sided: %s\n', ...
        name{1},q.observed_ratio,q.null_interval_95,q.upper_tail_p_value,q.decision,q.one_sided_decision);
end
end

function local_plot(p,folder)
names={'kdomain','joint','independent'}; colors=lines(3);
fig=figure('Visible','off','Color','w','Position',[60 60 1100 420]);
subplot(1,2,1); hold on; hp=gobjects(3,1); ho=gobjects(3,1);
for k=1:3, q=p.(names{k}); hp(k)=histogram(q.null_ratio_samples,40,'Normalization','pdf', ...
    'DisplayStyle','stairs','LineWidth',1.3,'EdgeColor',colors(k,:));
    ho(k)=xline(q.observed_ratio,'--','Color',colors(k,:),'LineWidth',1.2); end
grid on; xlabel('||P||_F/||C||_F'); ylabel('null PDF'); title('F=64 properness null');
subplot(1,2,2); hold on;
for k=1:3, q=sort(p.(names{k}).null_ratio_samples); plot(q,(1:numel(q)).'/numel(q),'Color',colors(k,:),'LineWidth',1.3); xline(p.(names{k}).observed_ratio,'--','Color',colors(k,:)); end
grid on; xlabel('||P||_F/||C||_F'); ylabel('empirical CDF');
legend([hp;ho],[strcat(names,' null'),strcat(names,' observed')],'Location','bestoutside');
exportgraphics(fig,fullfile(folder,'f64_properness_null.png'),'Resolution',180); close(fig);
end
