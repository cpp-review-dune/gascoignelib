function L2error_gen()
    % exact results
    u1 = dlmread('u.1.bup', ' ', [1 0 1072 5]);
    u2 = dlmread('u.2.bup', ' ', [1 0 1072 5]);
    u3 = dlmread('u.3.bup', ' ', [1 0 1072 5]);
    u4 = dlmread('u.4.bup', ' ', [1 0 1072 5]);
    u5 = dlmread('u.5.bup', ' ', [1 0 1072 5]);
    u6 = dlmread('u.6.bup', ' ', [1 0 1072 5]);
    u7 = dlmread('u.7.bup', ' ', [1 0 1072 5]);
    u8 = dlmread('u.8.bup', ' ', [1 0 1072 5]);
    % Results from iteration no. 1
    p1I1 = dlmread('C02tc52I1/p.0.bup', ' ', [1 0 1072 5]);
    p2I1 = dlmread('C02tc52I1/p.1.bup', ' ', [1 0 1072 5]);
    p3I1 = dlmread('C02tc52I1/p.2.bup', ' ', [1 0 1072 5]);
    p4I1 = dlmread('C02tc52I1/p.3.bup', ' ', [1 0 1072 5]);
    p5I1 = dlmread('C02tc52I1/p.4.bup', ' ', [1 0 1072 5]);
    p6I1 = dlmread('C02tc52I1/p.5.bup', ' ', [1 0 1072 5]);
    p7I1 = dlmread('C02tc52I1/p.6.bup', ' ', [1 0 1072 5]);
    p8I1 = dlmread('C02tc52I1/p.7.bup', ' ', [1 0 1072 5]);
    % Results from iteration no. 2
    p1I2 = dlmread('C02tc52I2/p.0.bup', ' ', [1 0 1072 5]);
    p2I2 = dlmread('C02tc52I2/p.1.bup', ' ', [1 0 1072 5]);
    p3I2 = dlmread('C02tc52I2/p.2.bup', ' ', [1 0 1072 5]);
    p4I2 = dlmread('C02tc52I2/p.3.bup', ' ', [1 0 1072 5]);
    p5I2 = dlmread('C02tc52I2/p.4.bup', ' ', [1 0 1072 5]);
    p6I2 = dlmread('C02tc52I2/p.5.bup', ' ', [1 0 1072 5]);
    p7I2 = dlmread('C02tc52I2/p.6.bup', ' ', [1 0 1072 5]);
    p8I2 = dlmread('C02tc52I2/p.7.bup', ' ', [1 0 1072 5]);
    % Results from iteration no. 3
    p1I3 = dlmread('C02tc52I3/p.0.bup', ' ', [1 0 1072 5]);
    p2I3 = dlmread('C02tc52I3/p.1.bup', ' ', [1 0 1072 5]);
    p3I3 = dlmread('C02tc52I3/p.2.bup', ' ', [1 0 1072 5]);
    p4I3 = dlmread('C02tc52I3/p.3.bup', ' ', [1 0 1072 5]);
    p5I3 = dlmread('C02tc52I3/p.4.bup', ' ', [1 0 1072 5]);
    p6I3 = dlmread('C02tc52I3/p.5.bup', ' ', [1 0 1072 5]);
    p7I3 = dlmread('C02tc52I3/p.6.bup', ' ', [1 0 1072 5]);
    p8I3 = dlmread('C02tc52I3/p.7.bup', ' ', [1 0 1072 5]);

    e1I1 = (u1 - p1I1);
    e1I2 = (u1 - p1I2);
    e1I3 = (u1 - p1I3);

    e2I1 = (u2 - p2I1);
    e2I2 = (u2 - p2I2);
    e2I3 = (u2 - p2I3);

    e3I1 = (u3 - p3I1);
    e3I2 = (u3 - p3I2);
    e3I3 = (u3 - p3I3);

    e4I1 = (u4 - p4I1);
    e4I2 = (u4 - p4I2);
    e4I3 = (u4 - p4I3);

    e5I1 = (u5 - p5I1);
    e5I2 = (u5 - p5I2);
    e5I3 = (u5 - p5I3);

    e6I1 = (u6 - p6I1);
    e6I2 = (u6 - p6I2);
    e6I3 = (u6 - p6I3);

    e7I1 = (u7 - p7I1);
    e7I2 = (u7 - p7I2);
    e7I3 = (u7 - p7I3);

    e8I1 = (u8 - p8I1);
    e8I2 = (u8 - p8I2);
    e8I3 = (u8 - p8I3);

    err = zeros(8,6);

    v_norm = [norm(u1(:,2:3), 'fro'), norm(u2(:,2:3), 'fro'), norm(u3(:,2:3), 'fro'), norm(u4(:,2:3), 'fro'), norm(u5(:,2:3), 'fro'), norm(u6(:,2:3), 'fro'), norm(u7(:,2:3), 'fro'), norm(u8(:,2:3), 'fro')]

    p_norm = [norm(u1(:,1), 'fro'), norm(u2(:,1), 'fro'), norm(u3(:,1), 'fro'), norm(u4(:,1), 'fro'), norm(u5(:,1), 'fro'), norm(u6(:,1), 'fro'), norm(u7(:,1), 'fro'), norm(u8(:,1), 'fro')]

    err(1,1) = norm(e1I1(:,2:3), 'fro');
    err(1,2) = norm(e1I2(:,2:3), 'fro');
    err(1,3) = norm(e1I3(:,2:3), 'fro');

    err(2,1) = norm(e2I1(:,2:3), 'fro');
    err(2,2) = norm(e2I2(:,2:3), 'fro');
    err(2,3) = norm(e2I3(:,2:3), 'fro');

    err(3,1) = norm(e3I1(:,2:3), 'fro');
    err(3,2) = norm(e3I2(:,2:3), 'fro');
    err(3,3) = norm(e3I3(:,2:3), 'fro');

    err(4,1) = norm(e4I1(:,2:3), 'fro');
    err(4,2) = norm(e4I2(:,2:3), 'fro');
    err(4,3) = norm(e4I3(:,2:3), 'fro');

    err(5,1) = norm(e5I1(:,2:3), 'fro');
    err(5,2) = norm(e5I2(:,2:3), 'fro');
    err(5,3) = norm(e5I3(:,2:3), 'fro');

    err(6,1) = norm(e6I1(:,2:3), 'fro');
    err(6,2) = norm(e6I2(:,2:3), 'fro');
    err(6,3) = norm(e6I3(:,2:3), 'fro');

    err(7,1) = norm(e7I1(:,2:3), 'fro');
    err(7,2) = norm(e7I2(:,2:3), 'fro');
    err(7,3) = norm(e7I3(:,2:3), 'fro');

    err(8,1) = norm(e8I1(:,2:3), 'fro');
    err(8,2) = norm(e8I2(:,2:3), 'fro');
    err(8,3) = norm(e8I3(:,2:3), 'fro');
    %%Pressure
    err(1,4) = norm(e1I1(:,1), 'fro');
    err(1,5) = norm(e1I2(:,1), 'fro');
    err(1,6) = norm(e1I3(:,1), 'fro');

    err(2,4) = norm(e2I1(:,1), 'fro');
    err(2,5) = norm(e2I2(:,1), 'fro');
    err(2,6) = norm(e2I3(:,1), 'fro');

    err(3,4) = norm(e3I1(:,1), 'fro');
    err(3,5) = norm(e3I2(:,1), 'fro');
    err(3,6) = norm(e3I3(:,1), 'fro');

    err(4,4) = norm(e4I1(:,1), 'fro');
    err(4,5) = norm(e4I2(:,1), 'fro');
    err(4,6) = norm(e4I3(:,1), 'fro');

    err(5,4) = norm(e5I1(:,1), 'fro');
    err(5,5) = norm(e5I2(:,1), 'fro');
    err(5,6) = norm(e5I3(:,1), 'fro');

    err(6,4) = norm(e6I1(:,1), 'fro');
    err(6,5) = norm(e6I2(:,1), 'fro');
    err(6,6) = norm(e6I3(:,1), 'fro');

    err(7,4) = norm(e7I1(:,1), 'fro');
    err(7,5) = norm(e7I2(:,1), 'fro');
    err(7,6) = norm(e7I3(:,1), 'fro');

    err(8,4) = norm(e8I1(:,1), 'fro');
    err(8,5) = norm(e8I2(:,1), 'fro');
    err(8,6) = norm(e8I3(:,1), 'fro');

    err_tbl = array2table(err, 'VariableNames', {'V_L2_I1', 'V_L2_I2','V_L2_I3','P_L2_I1', 'P_L2_I2', 'P_L2_I3'})
    writetable(err_tbl, 'abs_errData.csv','Delimiter',',')


    err(1,1) = err(1,1)/v_norm(1);
    err(1,2) = err(1,2)/v_norm(1);
    err(1,3) = err(1,3)/v_norm(1);

    err(2,1) = err(2,1)/v_norm(2);
    err(2,2) = err(2,2)/v_norm(2);
    err(2,3) = err(2,3)/v_norm(2);

    err(3,1) = err(3,1)/v_norm(3);
    err(3,2) = err(3,2)/v_norm(3);
    err(3,3) = err(3,3)/v_norm(3);

    err(4,1) = err(4,1)/v_norm(4);
    err(4,2) = err(4,2)/v_norm(4);
    err(4,3) = err(4,3)/v_norm(4);

    err(5,1) = err(5,1)/v_norm(5);
    err(5,2) = err(5,2)/v_norm(5);
    err(5,3) = err(5,3)/v_norm(5);

    err(6,1) = err(6,1)/v_norm(6);
    err(6,2) = err(6,2)/v_norm(6);
    err(6,3) = err(6,3)/v_norm(6);

    err(7,1) = err(7,1)/v_norm(7);
    err(7,2) = err(7,2)/v_norm(7);
    err(7,3) = err(7,3)/v_norm(7);

    err(8,1) = err(8,1)/v_norm(8);
    err(8,2) = err(8,2)/v_norm(8);
    err(8,3) = err(8,3)/v_norm(8);
    %%Pressure
    err(1,4) = err(1,4)/p_norm(1);
    err(1,5) = err(1,5)/p_norm(1);
    err(1,6) = err(1,6)/p_norm(1);

    err(2,4) = err(2,4)/p_norm(2);
    err(2,5) = err(2,5)/p_norm(2);
    err(2,6) = err(2,6)/p_norm(2);

    err(3,4) = err(3,4)/p_norm(3);
    err(3,5) = err(3,5)/p_norm(3);
    err(3,6) = err(3,6)/p_norm(3);

    err(4,4) = err(4,4)/p_norm(4);
    err(4,5) = err(4,5)/p_norm(4);
    err(4,6) = err(4,6)/p_norm(4);

    err(5,4) = err(5,4)/p_norm(5);
    err(5,5) = err(5,5)/p_norm(5);
    err(5,6) = err(5,6)/p_norm(5);

    err(6,4) = err(6,4)/p_norm(6);
    err(6,5) = err(6,5)/p_norm(6);
    err(6,6) = err(6,6)/p_norm(6);

    err(7,4) = err(7,4)/p_norm(7);
    err(7,5) = err(7,5)/p_norm(7);
    err(7,6) = err(7,6)/p_norm(7);

    err(8,4) = err(8,4)/p_norm(8);
    err(8,5) = err(8,5)/p_norm(8);
    err(8,6) = err(8,6)/p_norm(8);

    err_tbl = array2table(err, 'VariableNames', {'V_L2_I1', 'V_L2_I2','V_L2_I3','P_L2_I1', 'P_L2_I2', 'P_L2_I3'})
    writetable(err_tbl, 'rel_errData.csv','Delimiter',',')
end
