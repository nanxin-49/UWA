function arrivals = read_bellhop_arrivals_ascii_vertical(file)
%READ_BELLHOP_ARRIVALS_ASCII_VERTICAL Read standard ASCII Bellhop arrivals.
fid=fopen(file,'r');
if fid<0, error('Cannot open Bellhop arrivals: %s',file); end
cleanup=onCleanup(@()fclose(fid));
format_line=strtrim(fgetl(fid));
if ~contains(upper(format_line),'2D'), error('Expected a 2-D ASCII arrivals file.'); end
f=fscanf(fid,'%f',1); ns=fscanf(fid,'%d',1); sd=fscanf(fid,'%f',ns);
nd=fscanf(fid,'%d',1); rd=fscanf(fid,'%f',nd);
nr=fscanf(fid,'%d',1); rr=fscanf(fid,'%f',nr);
template=struct('frequency_hz',NaN,'source_depth_m',NaN,'receiver_depth_m',NaN, ...
    'receiver_range_m',NaN,'amplitude_complex',complex(NaN),'amplitude',NaN,'phase_deg',NaN, ...
    'delay_s',NaN,'delay_imag_s',NaN,'source_angle_deg',NaN,'receiver_angle_deg',NaN, ...
    'top_bounce_count',NaN,'bottom_bounce_count',NaN);
max_arrivals=fscanf(fid,'%d',1);
rows=repmat(template,max(1,ns*nd*nr*max_arrivals),1);
row_count=0;
for is=1:ns
    for id=1:nd
        for ir=1:nr
            n=fscanf(fid,'%d',1); values=fscanf(fid,'%f',[8,n]);
            if size(values,2)~=n, error('Incomplete arrivals block.'); end
            for aa=1:n
                q=values(:,aa); z=q(1)*exp(1i*deg2rad(q(2)));
                row_count=row_count+1;
                rows(row_count)=struct('frequency_hz',f,'source_depth_m',sd(is), ...
                    'receiver_depth_m',rd(id),'receiver_range_m',rr(ir), ...
                    'amplitude_complex',z,'amplitude',abs(z),'phase_deg',q(2), ...
                    'delay_s',q(3),'delay_imag_s',q(4),'source_angle_deg',q(5), ...
                    'receiver_angle_deg',q(6),'top_bounce_count',q(7), ...
                    'bottom_bounce_count',q(8));
            end
        end
    end
end
rows=rows(1:row_count);
arrivals=struct2table(rows);
clear cleanup
end
