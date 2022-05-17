clear;
%% Plot amplification functions of stabilization IFTSRK schemes developed in H. Zhang, X. Qian, J. Xia, S. Song, Unconditionally maximum-principle-preserving parametric integrating factor two-step Runge-Kutta schemes for parabolic equations, 2022
plot_amplification_function = 1;
red = [228, 26, 28]/255;
blue = [55, 126, 184]/255;
green = [77, 175, 74]/255;
purple = [152, 78, 163]/255;
orange = [255, 127, 0]/255;
brown = [166, 86, 40]/255;
c0 = [ 255, 237, 111]*.85/255;
c1 = [ 247, 129, 191]/255;
c2 = [ 153, 153, 153]/255;
colors(1,:) = red;
colors(3,:) = blue;
colors(2,:) = orange;
colors(4,:) = green;
colors(5,:) = brown;%orange;%brown;
colors(6,:) = purple;
colors(7,:) = c0;
colors(8,:) = c1;
colors(9,:) = c2;
colors(10,:) = red;
colors(11,:) = blue;
colors(12,:) = orange;
colors(13,:) = green;
colors(14,:) = brown;
colors(15,:) = purple;
colors(16,:) = c0;
colors(17,:) = c1;
colors(18,:) = c2;

fs = 20;
 
markers = ['o', '>', 's', '.', '*', 's', 'x', 'o', '<', 'o', 's', '.', '*', '>', 'x', 'o'];

if plot_amplification_function == 1 
    TSRK_flagv = [222 323 424 525 625 726 827 927 1027 1128];
    for kk = 1:length(TSRK_flagv)
        TSRK_flag = TSRK_flagv(kk);
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
        Dtheta = [D;theta(2:-1:1)]; ABhat = [Ahat; Bhat]; AB = [A; B];
        tildeD = [1 0; 0 1; D]; tildeA = [0 zeros(1,size(A,2)); Ahat, A]; tildeB = [Bhat B];
        e = ones(size(A,2)+1,1); l = [1 0]';
        tildeD = [tildeD; theta]; tildeAB = [tildeA; tildeB];
        c = tildeAB*e - tildeD*l;
        xv = 0:0.5:30;
        syms x; g_exp = cell(length(c)); g_exp0 = g_exp;
        g_fun = zeros(length(c), length(xv)); g_fundu = g_fun;
        g_fun1 = g_fun;phi{1} = @(x) 1; phi{2} = @(x) 1 + x;phi{3} = @(x) 1 + x + 1/2*x.^2;phi{4} = @(x) 1 + x + 1/2*x.^2 + 1/6*x.^3;phi{5} = @(x) 1 + x + 1/2*x.^2 + 1/6*x.^3 + 1/24*x.^4;phi{6} = @(x) 1 + x + 1/2*x.^2 + 1/6*x.^3 + 1/24*x.^4+1/120*x.^5;
        phi{7} = @(x) 1 + x + 1/2*x.^2 + 1/6*x.^3 + 1/24*x.^4+1/120*x.^5+1/720*x.^6;
        for i = 3:length(c)
            % Compute amplificaiton functions
            g_fun(i,:) = exp(-(1+c(i))*xv) .* ( tildeD(i,1)*exp((1-l(1))*xv) + tildeD(i,2)*exp((1-l(2))*xv));
            for j = 1:i-1
                g_fun(i,:) = g_fun(i,:) + exp(-(1+c(i))*xv) .* tildeAB(i,j).* exp((1+c(j))*xv) .* xv;
            end
        end
        figure;
        for i = 3:length(c)
            plot(xv, g_fun(i,:), 'marker', markers(i), 'markersize', 3, 'color', colors(i,:), 'linewidth', 1.5);%colorv(i,:) ); hold on;
            hold on;
        end
        xlabel('z','fontsize', fs-4, 'interpreter', 'latex'); ylabel('g(z)','fontsize', fs-4, 'interpreter', 'latex');
        switch stage
            case 2
                hl = legend('g_1', 'g_2', 'location', 'east');
            case 3
                hl = legend('g_1', 'g_2', 'g_3',  'location', 'east');
            case 4
                hl = legend('g_1', 'g_2', 'g_3', 'g_4', 'location', 'east');
            case 5
                hl = legend('g_1', 'g_2', 'g_3', 'g_4', 'g_5', 'location', 'east');
            case 6
                hl = legend('g_1', 'g_2', 'g_3', 'g_4', 'g_5', 'g_6', 'location', 'east');
            case 7
                hl = legend('g_1', 'g_2', 'g_3', 'g_4', 'g_5', 'g_6', 'g_7', 'location', 'east');
            case 8
                hl = legend('g_1', 'g_2', 'g_3', 'g_4', 'g_5', 'g_6', 'g_7', 'g_8', 'location', 'east');
            case 9
                hl = legend('g_1', 'g_2', 'g_3', 'g_4', 'g_5', 'g_6', 'g_7', 'g_8', 'g_9', 'location', 'east');
            case 10
                hl = legend('g_1', 'g_2', 'g_3', 'g_4', 'g_5', 'g_6', 'g_7', 'g_8', 'g_9', 'g_{10}', 'location', 'east');
            case 11
                hl = legend('g_1', 'g_2', 'g_3', 'g_4', 'g_5', 'g_6', 'g_7', 'g_8', 'g_9', 'g_{10}', 'g_{11}', 'location', 'east');
        end
        set(hl, 'fontsize', fs-8);
        fig_init([0 0 6 6]);
    end
end

