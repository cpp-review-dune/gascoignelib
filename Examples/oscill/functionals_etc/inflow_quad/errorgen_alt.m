function errorgen_alt(file_exa, file_approx, div)
    A = dlmread(file_approx);
    E = dlmread(file_exa);
    Time = A(:,1);
    A = A(:,2:end);
    E = E(div:div:end,2:end);
    %sigma = std(E)
    %mu = mean(E)
    maxv = max(E);
    rel_err = abs((A-E))./maxv;
    out = [Time rel_err];
    outstr = strcat('ALTfct_error_seq_p', file_approx(16:end));
    dlmwrite(outstr, out, ' ');
end

