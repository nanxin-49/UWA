function value = fast_fk_height_interpolate(nodes, samples, query)
%FAST_FK_HEIGHT_INTERPOLATE Barycentric Chebyshev-Lobatto interpolation.
% Validation-only, O(q*N*r) with a fixed, small collocation count q.
nodes=nodes(:); query=query(:); q=numel(nodes);
assert(size(samples,1)==q && q>=2);
assert(all(query>=min(nodes)-1e-12 & query<=max(nodes)+1e-12));
w=(-1).^(0:q-1)'; w([1 end])=w([1 end])/2;
value=complex(zeros(numel(query),size(samples,2)));
for first=1:256:numel(query)
 rows=first:min(first+255,numel(query)); diff=query(rows)-nodes.';
 [distance,nearest]=min(abs(diff),[],2);
 exact=distance<=16*eps(max(1,max(abs(nodes))));
 ratio=w.'./diff; ratio(exact,:)=0;
 denominator=sum(ratio,2); denominator(exact)=1;
 value(rows,:)=(ratio*samples)./denominator;
 ii=find(exact);
 if ~isempty(ii), value(rows(ii),:)=samples(nearest(ii),:); end
end
end
