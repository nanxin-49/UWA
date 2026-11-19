function [value,owner] = d1b_tiles(count, tile_rows, workers, kernel)
%D1B_TILES Execute disjoint row tiles in fixed order on the client or pool.
blocks = ceil(count/tile_rows);
if workers==0
    slots = cell(blocks,1);
    for b = 1:blocks
        slots{b} = kernel(local_rows(b,tile_rows,count));
    end
    value = vertcat(slots{:}); owner = zeros(count,1);
    return
end
assert(blocks>=workers && gcp('nocreate').NumWorkers==workers);
spmd
    block_ids = spmdIndex:spmdSize:blocks;
    pieces = cell(numel(block_ids),1);
    for j = 1:numel(block_ids)
        pieces{j} = kernel(local_rows(block_ids(j),tile_rows,count));
    end
    payload = struct('pieces',{pieces},'pid',feature('getpid'));
end
slots = cell(blocks,1); owner = zeros(count,1);
for w = 1:workers
    part = payload{w};
    ids = w:workers:blocks;
    slots(ids) = part.pieces;
    for b = ids
        owner(local_rows(b,tile_rows,count)) = part.pid;
    end
end
value = vertcat(slots{:});
end

function rows = local_rows(block,width,count)
rows = (block-1)*width+1:min(block*width,count);
end
