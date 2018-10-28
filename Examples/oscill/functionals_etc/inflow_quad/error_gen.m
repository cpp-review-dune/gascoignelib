function error_gen(file_exa, file_approx, div)
    A = dlmread(file_approx);
    E = dlmread(file_exa);
    Time = A(:,1);
    A = A(:,2:end);
    E = E(div:div:end,2:end);
    rel_err = abs((A-E)./(E+1));
    out = [Time rel_err];
    outstr = strcat('fct_error_seq_p', file_approx(16:end));
    dlmwrite(outstr, out, ' ');
end