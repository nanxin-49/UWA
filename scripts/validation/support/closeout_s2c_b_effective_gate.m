function gate = closeout_s2c_b_effective_gate(run_dir)
%CLOSEOUT_S2C_B_EFFECTIVE_GATE Apply only the explicit K_HIGH budget amendment.
arguments
    run_dir (1,:) char
end
assert(isfolder(run_dir) && contains(run_dir,'20261023_s2c_b_khigh_'));
initial=fullfile(run_dir,'protocol','P1_initial_gate.json');
gate=jsondecode(fileread(initial));
assert(gate.subbudget_s==1800 && gate.core_remaining_s==5791.339895);
gate.original_subbudget_s=gate.subbudget_s;
gate.budget_amendment_sha256='';
amendment=fullfile(run_dir,'protocol','P1_budget_amendment_2400.json');
if ~isfile(amendment), return; end
a=jsondecode(fileread(amendment));
auth=fullfile(run_dir,'protocol','authorization.json');
assert(strcmp(a.authorized_case,'K_HIGH') && a.explicit_user_authorization);
assert(strcmp(a.scope,'S2c-B K_HIGH only') && ~a.other_candidates_authorized);
assert(a.original_subbudget_s==1800 && a.new_subbudget_s==2400);
assert(strcmp(a.initial_gate_sha256,closeout_s2b_hash(initial)));
assert(strcmp(a.original_authorization_sha256,closeout_s2b_hash(auth)));
assert(a.shared_core_cap_s==7200 && a.S3A_reserve_s==3000);
assert(a.new_subbudget_s<=gate.core_remaining_s-gate.S3A_reserve_s);
gate.subbudget_s=a.new_subbudget_s;
gate.budget_amendment_sha256=closeout_s2b_hash(amendment);
end
