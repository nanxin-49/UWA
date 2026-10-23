function export_two_dimensional_closeout_s2b_sensitivity(out)
% Export every computed sensitivity complex value without reevaluation.
assert(~isfile(fullfile(out,'data/sensitivity_fields.csv'))); rows={};
p=load(fullfile(out,'data/p0_data.mat'),'fields'); f=p.fields;
for id={'bie1','bie2','bie3','pe1','pe2','pe3','pe4','fk1','fk2','fk3','fk4','flat'}
    q=load(fullfile(out,'checkpoints',[id{1} '.mat']),'payload'); d=q.payload;
    variants=struct('primary',d.u); xs=d.x_m; zs=d.z_m;
    if startsWith(id{1},'bie')
        names=fieldnames(d.level_fields);
        for j=1:numel(names)
            uv=d.level_fields.(names{j});
            for k=1:size(uv,2), variants.(sprintf('%s_%d',names{j},k))=uv(:,k); end
        end
    elseif startsWith(id{1},'fk')
        for k=1:4, variants.(sprintf('L%d_N%d',d.runs{k}.surface_x_m(end),numel(d.runs{k}.surface_x_m)))=d.runs{k}.receiver_field; end
    elseif startsWith(id{1},'pe')
        variants.N1968_same_source=d.grid_double_u;
        for i=1:411
            rows{end+1,1}=table(string(id{1}),"half_step_line411",i,f.x_m(i),3,real(d.step_half_line(i)),imag(d.step_half_line(i)), ...
                'VariableNames',{'job_id','variant','probe_index','x_m','z_m','u_real','u_imag'}); %#ok<AGROW>
        end
    end
    names=fieldnames(variants);
    for j=1:numel(names)
        v=variants.(names{j});
        for i=1:numel(v)
            rows{end+1,1}=table(string(id{1}),string(names{j}),i,xs(i),zs(i),real(v(i)),imag(v(i)), ...
                'VariableNames',{'job_id','variant','probe_index','x_m','z_m','u_real','u_imag'}); %#ok<AGROW>
        end
    end
end
writetable(vertcat(rows{:}),fullfile(out,'data/sensitivity_fields.csv'));
end
