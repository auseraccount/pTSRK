    %% Verify the order conditions of pTSRK schemes developed in H. Zhang, X. Qian, J. Xia, S. Song, Unconditionally maximum-principle-preserving parametric integrating factor two-step Runge-Kutta schemes for parabolic equations, 2022
    
    phi = cell(9, 1);
    phi{1} = @(x) 1;
    for i = 1:8
        phi{i+1} = @(x) phi{i}(x) + x^i/factorial(i);
    end
    TSRK_flag = 324; % 3 stages, 2 steps, 4-th order
    stage = floor(TSRK_flag/100);
    step  = mod(floor(TSRK_flag/10), 10);
    order = mod(TSRK_flag, 10);
    if ~exist('SSPIF-TSRK-methods-master')
        fprintf('Download TSRK file from https://github.com/SSPmethods/SSPIF-TSRK-methods');
        urlwrite('https://codeload.github.com/SSPmethods/SSPIF-TSRK-methods/zip/refs/heads/master', 'SSPIF-TSRK-methods-master.zip');
        unzip('SSPIF-TSRK-methods-master.zip', '.');
    end
    tsrkfilename = ['./SSPIF-TSRK-methods-master/eSSPTSRKplus methods/' ...
        num2str(stage) 's' num2str(step) 'k' num2str(order) 'pSSPTSRK+.mat'];
    load(tsrkfilename);
    fprintf('TSRK method loaded: step = %d, stage = %d, order = %d\n', step, stage, order);
    % Transform the TSRK coefficients into Butcher tableau
    Dtheta = [D;theta(2:-1:1)]; ABhat = [Ahat; Bhat]; AB = [A; B];
    tildeD = [1 0; 0 1; D]; tildeA = [0 zeros(1,size(A,2)); Ahat, A]; tildeB = [Bhat B];
    e = ones(size(A,2)+1,1); L = [1 0]';
    tildeD = [tildeD; theta]; tildeAB = [tildeA; tildeB];
    c = tildeAB*e - tildeD*L;
    fprintf('Butcher tableau of TSRK(%d, %d):\n', stage, order);
    [tildeD tildeAB]
    syms x;
    modify_flag = 2; %1: pTSRK1; 2: pTSRK2; 0: original IFRK
    hatD(1, 1) = sym(1); hatD(2, 2) = sym(1); orderp = order;
    if modify_flag == 1 % pTSRK1
        hatc(1) = sym(-1); hatc(2) = sym(0); psi(1) = sym(1); psi(2) = exp((1+c(2))*x);
        for i = 3:length(c)
            % To accelerate computations, we apply Taylor expansions to psi
            psi(i) = taylor((tildeD(i,1) * exp((1+c(1))*x) + tildeD(i,2)*exp((1+c(2))*x)), x, 0, 'order', orderp+2);
            for j = 1:i-1
                psi(i) = taylor(psi(i) + tildeAB(i,j).* exp((1+c(j))*x) .* x, x, 0, 'order', orderp+2);
            end
            for j = 1:2
                hatD(i,j) = taylor(1/psi(i) * (tildeD(i,j) .* exp((1+c(j))*x)), x, 0, 'order', orderp+2);
                for k = 1:i-1
                    hatD(i,j) = taylor(hatD(i,j) + 1/psi(i) * x * tildeAB(i,k)*exp((1+c(k))*x) * hatD(k,j), x, 0, 'order', orderp+2);
                end
            end
            hatc(i) = -hatD(i,1);
            for j = 1:i-1
                hatAB(i,j) = taylor(1/psi(i) * (tildeAB(i,j)*exp((1+c(j))*x)), x, 0, 'order', orderp+2);
                for k = 3: i-1
                    hatAB(i,j) = taylor(hatAB(i,j) + 1/psi(i) * x * tildeAB(i,k) * exp((1+c(k))*x) * hatAB(k,j), x, 0, 'order', orderp+2);
                end
                hatc(i) = hatc(i) + hatAB(i,j);
            end
        end
    elseif modify_flag == 2
        hatc(1) = sym(-1); hatc(2) = sym(0);
        psi(1) = sym(phi{orderp}((1+c(1))*x)); psi(2) = phi{orderp}((1+c(2))*x);
        for i = 3:length(c)
            psi(i) =  taylor(tildeD(i,1)*psi(1) + tildeD(i,2)*psi(2), x, 0, 'order', orderp+2);
            for j = 1:i-1
                psi(i) = taylor(psi(i) + tildeAB(i,j).* psi(j) .* x, x, 0, 'order', orderp+2);
            end
            for j = 1:2
                hatD(i,j) = taylor(1/psi(i) * tildeD(i,j) * psi(j), x, 0, 'order', orderp+2);
                for k = 1:i-1
                    hatD(i,j) = taylor(hatD(i,j) + 1/psi(i) * x * tildeAB(i,k) * psi(k) * hatD(k,j), x, 0, 'order', orderp+2);
                end
            end
            hatc(i) = -hatD(i,1);
            for j = 1:i-1
                hatAB(i,j) = taylor(1/psi(i) * tildeAB(i,j) * psi(j), x, 0, 'order', orderp+2);
                for k = 3: i-1
                    hatAB(i,j) = taylor(hatAB(i,j) + 1/psi(i) * x * tildeAB(i,k) * psi(k) * hatAB(k,j), x, 0, 'order', orderp+2);
                end
                hatc(i) = hatc(i) + hatAB(i,j);
            end
        end
    end
    
    hatA = hatAB(1:end-1, :); hatB = hatAB(end, :);
    hattheta = hatD(end, :); hatc = hatc(:);
    hatc1 = hatc(1:end-1);  c1 = c(1:end-1);
    %% Verification of order conditions
    for k = 1:order - 1
        for m = 0:k-1
            tmp = taylor(1/factorial(k) * (c.^k - hatD*(-L).^k) - 1/factorial(k-1) * hatAB*(c1.^m .* hatc1.^(k-1-m)), x, 0, 'order', orderp);
            tau(1:length(c),k,m+1) = tmp;
        end
    end
    disp('Order 1:')
    vpa(taylor(hatB * e - (1 + hattheta*L), x, 0, 'order', order+1), 6)
    if order == 1
        return;
    end
    disp('Order 2:')
    for m = 0:1
        vpa(taylor(hatB * (c1.^m .* hatc1.^(1-m)) - (1 - hattheta * L.^2)/2, x, 0, 'order', order), 6)
    end
    if order == 2
        return;
    end
    disp('Order 3:')
    for m = 0:2
        vpa(taylor(hatB * (c1.^m .* hatc1.^(2-m)) - (1 + hattheta*L.^3)/3, x, 0, 'order', order-1), 6)
    end
    for m = 0:1
        vpa(taylor(hatB * tau(1:end-1, 2, 1+m), x, 0, 'order', order-1), 6)
    end
    if order == 3
        return;
    end
    disp('Order 4:')
    for m = 0:3
        vpa(taylor(hatB * (c1.^m .* hatc1.^(3-m)) - (1 - hattheta*L.^4)/4, x, 0, 'order', order-2), 6)
    end
    for m = 0:1
        vpa(taylor( hatB * hatA * tau(1:end-1, 2,1+m), x, 0, 'order', order-2), 6)
    end
    for m1 = 0:1
        for m2 = 0:1
            vpa(taylor(( hatB .* (c1.^m .* hatc1.^(1-m))') * tau(1:end-1, 2, 1+m2), x, 0, 'order', order-2), 6)
        end
    end
    for m = 1:2
        vpa(taylor( hatB * tau(1:end-1, 3, 1+m), x, 0, 'order', order-2), 6)
    end
    if order == 4
        return;
    end
    disp('Order 5:')
    for m = 0:4
        vpa(taylor(hatB * (c1.^m .* hatc1.^(4-m)) - (1 + hattheta * L.^5)/5, x, 0, 'order', order-3), 6)
    end
    for m = 0:2
        vpa(taylor(hatB * hatA * tau(1:end-1, 3, 1+m), x, 0, 'order', order-3), 6)
    end
    for m1 = 0:1
        for m2 = 0:2
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(1-m1))') * tau(1:end-1, 3, 1+m2), x, 0, 'order', order-3), 6)
        end
    end
    for m = 0:3
        vpa(taylor(hatB * tau(1:end-1, 4, 1+m), x, 0, 'order', order-3), 6)
    end
    for m = 0:1
        vpa(taylor(tau(1:end-1, 2, 1+m), x, 0, 'order', order-3), 6)
    end
    if order == 5
        return;
    end
    disp('Order 6:')
    for m = 0:5
        vpa(taylor(hatB * (c1.^m .* hatc1.^(5-m)) - (1 - hattheta * L.^6) / 6, x, 0, 'order', order-4), 6)
    end
    for m = 0:3
        vpa(taylor(hatB * hatA * tau(1:end-1, 4, 1+m), x, 0, 'order', order-4), 6)
    end
    for m1 = 0:1
        for m2 = 0:3
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(1-m1))') * tau(1:end-1, 4, 1+m2), x, 0, 'order', order-4), 6)
        end
    end
    for m = 0:4
        vpa(taylor(hatB * tau(1:end-1, 5, 1+m), x, 0, 'order', order-4), 6)
    end
    for m = 0:2
        vpa(taylor(hatB * hatA^2 * tau(1:end-1, 3, 1+m), x, 0, 'order', order-4), 6)
    end
    for m1 = 0:1
        for m2 = 0:2
            vpa(taylor((hatB * hatA .* (c1.^m1 .* hatc1.^(1-m1))') * tau(1:end-1, 3, 1+m2), x, 0, 'order', order-4), 6)
        end
    end
    for m1 = 0:1
        for m2 = 0:2
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(1-m1))') * hatA * tau(1:end-1, 3, 1+m2), x, 0, 'order', order-4), 6)
        end
    end
    for m1 = 0:2
        for m2 = 0:2
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(2-m1))') * tau(1:end-1, 3, 1+m2), x, 0, 'order', order-4), 6)
        end
    end
    if order == 6
        return;
    end
    disp('Order 7:')
    for m = 0:6
        vpa(taylor(hatB * (c1.^m .* hatc1.^(6-m)) - (1 + hattheta * L.^6)/7, x, 0, 'order', order-5), 6)
    end
    for m = 0:4
        vpa(taylor(hatB * hatA * tau(1:end-1, 5, 1+m), x, 0, 'order', order-5), 6)
    end
    for m1 = 0:1
        for m2 = 0:4
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(1-m1))') * tau(1:end-1, 5, 1+m2), x, 0, 'order', order-5), 6)
        end
    end
    for m = 0:5
        vpa(taylor(hatB * tau(1:end-1, 6, 1+m), x, 0, 'order', order-5), 6)
    end
    for m = 0:3
        vpa(taylor(hatB * hatA^2 * tau(1:end-1, 4, 1+m), x, 0, 'order', order-5), 6)
    end
    for m1 = 0:1
        for m2 = 0:3
            vpa(taylor((hatB * hatA .* (c1.^m1 .* hatc1.^(1-m1))') * tau(1:end-1, 4, 1+m2), x, 0, 'order', order-5), 6)
        end
    end
    for m1 = 0:1
        for m2 = 0:3
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(1-m1))') * hatA * tau(1:end-1, 4, 1+m2), x, 0, 'order', order-5), 6)
        end
    end
    for m1 = 0:2
        for m2 = 0:3
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(2-m1))') * tau(1:end-1, 4, 1+m2), x, 0, 'order', order-5), 6)
        end
    end
    for m = 0:2
        vpa(taylor(tau(1:end-1, 3, 1+m), x, 0, 'order', order-5), 6)
    end
    if order == 7
        return; end
    disp('Order 8:')
    for m = 0:7
        vpa(taylor(hatB * (c1.^m .* hatc1.^(7-m)) - (1 - hattheta * L.^8)/8, x, 0, 'order', order-6), 6)
    end
    for m = 0:5
        vpa(taylor(hatB * hatA * tau(1:end-1, 6, 1+m), x, 0, 'order', order-6), 6)
    end
    for m1 = 0:1
        for m2 = 0:5
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(1-m1))') * tau(1:end-1, 6, 1+m2), x, 0, 'order', order-6), 6)
        end
    end
    for m = 0:6
        vpa(taylor(hatB * tau(1:end-1, 7, 1+m), x, 0, 'order', order-6), 6)
    end
    
    for m = 0:3
        vpa(taylor(hatB * hatA^3 * tau(1:end-1, 4, 1+m), x, 0, 'order', order-6), 6)
    end
    for m = 0:4
        vpa(taylor(hatB * hatA^2 * tau(1:end-1, 5, 1+m), x, 0, 'order', order-6), 6)
    end
    for m1 = 0:1
        for m2 = 0:3
            vpa(taylor((hatB * hatA^2 .* (c1.^m1 .* hatc1.^(1-m1))') * tau(1:end-1, 4, 1+m2), x, 0, 'order', order-6), 6)
        end
    end
    for m1 = 0:1
        for m2 = 0:3
            vpa(taylor((hatB * hatA .* (c1.^m1 .* hatc1.^(1-m1))') * hatA * tau(1:end-1, 4, 1+m2), x, 0, 'order', order-6), 6)
        end
    end
    
    for m1 = 0:1
        for m2 = 0:4
            vpa(taylor((hatB * hatA .* (c1.^m1 .* hatc1.^(1-m1))') * tau(1:end-1, 5, 1+m2), x, 0, 'order', order-6), 6)
        end
    end
    for m1 = 0:2
        for m2 = 0:3
            vpa(taylor((hatB * hatA .* (c1.^m1 .* hatc1.^(2-m1))') * tau(1:end-1, 4, 1+m2), x, 0, 'order', order-6), 6)
        end
    end
    for m1 = 0:1
        for m2 = 0:3
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(1-m1))') * hatA^2 * tau(1:end-1, 4, 1+m2), x, 0, 'order', order-6), 6)
        end
    end
    for m1 = 0:1
        for m2 = 0:4
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(1-m1))') * hatA * tau(1:end-1, 5, 1+m2), x, 0, 'order', order-6), 6)
        end
    end
    
    for m1 = 0:1
        for m2 = 0:1
            for m3 = 0:3
                vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(1-m1))') * hatA .* (c1.^m2 .* hatc1.^(1-m2))' * tau(1:end-1, 4, 1+m3), x, 0, 'order', order-6), 6)
            end
        end
    end
    for m1 = 0:2
        for m2 = 0:3
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(2-m1))') * hatA *  tau(1:end-1, 4, 1+m2), x, 0, 'order', order-6), 6)
        end
    end
    for m1 = 0:2
        for m2 = 0:4
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(2-m1))') * tau(1:end-1, 5, 1+m2), x, 0, 'order', order-6), 6)
        end
    end
    for m1 = 0:3
        for m2 = 0:3
            vpa(taylor((hatB .* (c1.^m1 .* hatc1.^(3-m1))') * tau(1:end-1, 4, 1+m2), x, 0, 'order', order-6), 6)
        end
    end
