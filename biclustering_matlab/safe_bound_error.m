function [ub0, ub] = safe_bound_error(blk, At, C, b, y, Z2, Bt, y2, l, ub_max_eig)
  
  % Post-processing through error bounds
  
  Aty = sdpnalAtyfun(blk, At, y);
  Znew = ops(C, '-', Aty);
  if ~isempty(Bt)
    Bty = sdpnalAtyfun(blk, Bt, y2);
    Znew = ops(Znew, '-', Bty);
  end
  if ~isempty(Z2)
     Znew = ops(Znew, '-', Z2); 
  end

  eigtmp = eig(full(Znew{1}));
  idx = eigtmp < -1e-8; 
  ub0 = -(b'*y + l'*y2);
  pert = ub_max_eig * sum(eigtmp(idx));
  ub = ub0 - pert;
  
end