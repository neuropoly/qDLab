function acq = scd_scheme2struct(scheme)
    acq.G=scheme(:,4);
    acq.Delta=scheme(:,5);
    acq.delta=scheme(:,6);
    acq.TE=scheme(:,7);
    acq.scheme=scheme;
    acq.q=scheme(:,8);
end