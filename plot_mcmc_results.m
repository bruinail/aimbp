function [] = plot_mcmc_results(samples,model_true,type)

if(strcmp(type,'Phi_rate'))
    Phi = zeros([size(model_true.Phi) length(samples)]);
    for i = 1:length(samples)
        Phi(:,:,i) = samples(i).model.Phi;
    end
    figure;
    
    Phi_rs = squeeze(Phi(1,3,:));
    Phi_rs_true = model_true.Phi(1,3);
    subplot(3,2,1);
    hold on;
    histogram(Phi_rs,-0.05:0.01:1.05);
    line([Phi_rs_true,Phi_rs_true],[0,length(samples)],'Color','r','LineWidth',1);
    subplot(3,2,2);
    plot(1:length(samples),Phi_rs);
    
    Phi_rd = squeeze(Phi(2,4,:));
    Phi_rd_true = model_true.Phi(2,4);
    subplot(3,2,3);
    hold on;
    histogram(Phi_rd,-0.05:0.01:1.05);
    line([Phi_rd_true,Phi_rd_true],[0,length(samples)],'Color','r','LineWidth',1);
    subplot(3,2,4);
    plot(1:length(samples),Phi_rd);
    
    Phi_rh = squeeze(Phi(5,6,:));
    Phi_rh_true = model_true.Phi(5,6);
    subplot(3,2,5);
    hold on;
    histogram(Phi_rh,-0.05:0.01:1.05);
    line([Phi_rh_true,Phi_rh_true],[0,length(samples)],'Color','r','LineWidth',1);
    subplot(3,2,6);
    plot(1:length(samples),Phi_rh);
elseif(strcmp(type,'Upsilon'))
    Upsilon = zeros([size(model_true.Upsilon) length(samples)]);
    for i = 1:length(samples)
        Upsilon(:,:,i) = samples(i).model.Upsilon;
    end
    figure;
    
    Emax_nd = squeeze(Upsilon(1,1,:));
    Emax_nd_true = model_true.Upsilon(1,1);
    subplot(2,2,1);
    hold on;
    histogram(Emax_nd,11);
    line([Emax_nd_true,Emax_nd_true],[0,length(samples)],'Color','r','LineWidth',1);
    subplot(2,2,2);
    plot(1:length(samples),Emax_nd);
    
    EC50_nd = squeeze(Upsilon(1,2,:));
    EC50_nd_true = model_true.Upsilon(1,2);
    subplot(2,2,3);
    hold on;
    histogram(EC50_nd,11);
    line([EC50_nd_true,EC50_nd_true],[0,length(samples)],'Color','r','LineWidth',1);
    subplot(2,2,4);
    plot(1:length(samples),EC50_nd);
elseif(strcmp(type,'mu_0_baseline'))
    mu_0 = zeros([size(model_true.mu_0) length(samples)]);
    for i = 1:length(samples)
        mu_0(:,:,i) = samples(i).model.mu_0;
    end
    figure;
    
    mu_0_baseline_sbp = squeeze(mu_0(3,1,:));
    mu_0_baseline_sbp_true = model_true.mu_0(3,1);
    subplot(1,2,1);
    hold on;
    histogram(mu_0_baseline_sbp,11);
    line([mu_0_baseline_sbp_true,mu_0_baseline_sbp_true],[0,length(samples)],'Color','r','LineWidth',1);
    subplot(1,2,2);
    plot(1:length(samples),mu_0_baseline_sbp);
elseif(strcmp(type,'custom_paper'))
    Bmax = zeros([size(model_true.U.Bmax) length(samples)]);
    r_B = zeros([size(model_true.U.r_B) length(samples)]);
    Emax = zeros([size(model_true.U.Emax) length(samples)]);
    EC50 = zeros([size(model_true.U.EC50) length(samples)]);
    mu_0 = zeros([size(model_true.mu_0) length(samples)]);
    for i = 1:length(samples)
        Bmax(:,:,i) = samples(i).model.U.Bmax;
        r_B(:,:,i) = samples(i).model.U.r_B;
        Emax(:,:,i) = samples(i).model.U.Emax;
        EC50(:,:,i) = samples(i).model.U.EC50;
        mu_0(:,:,i) = samples(i).model.mu_0;
    end
    
    n_drugs = size(model_true.U.Emax,1);
    n_plots = 3+2*n_drugs;
    figure;
    
    Bmax = squeeze(Bmax(1,1,:));
    Bmax_true = model_true.U.Bmax(1,1);
    subplot(n_plots,3,1);
    hold on;
    histogram(Bmax,11);
    title('B_{max} Histogram');
    line([Bmax_true,Bmax_true],[0,length(samples)],'Color','r','LineWidth',1);
    subplot(n_plots,3,2);
    hold on;
    plot(1:length(samples),Bmax);
    chunk = floor(length(samples)/5);
    for i = 1:5
        chunk_mean = mean(Bmax((i-1)*chunk+1:i*chunk));
        line([(i-1)*chunk+1 i*chunk],[chunk_mean chunk_mean],'Color','r','LineWidth',1);
    end
    title('B_{max} Chain Trace');
    subplot(n_plots,3,3);
    autocorr(Bmax,'NumLags',100);
    title('B_{max} ACF');
    
    r_B = squeeze(r_B(1,1,:));
    r_B_true = model_true.U.r_B(1,1);
    subplot(n_plots,3,4);
    hold on;
    histogram(r_B,-0.05:0.01:1.05);
    title('r_B Histogram');
    line([r_B_true,r_B_true],[0,length(samples)],'Color','r','LineWidth',1);
    subplot(n_plots,3,5);
    hold on;
    plot(1:length(samples),r_B);
    chunk = floor(length(samples)/5);
    for i = 1:5
        chunk_mean = mean(r_B((i-1)*chunk+1:i*chunk));
        line([(i-1)*chunk+1 i*chunk],[chunk_mean chunk_mean],'Color','r','LineWidth',1);
    end
    title('r_B Chain Trace');
    subplot(n_plots,3,6);
    autocorr(r_B,'NumLags',100);
    title('r_B ACF');
    
    mu_0_baseline_sbp = squeeze(mu_0(3,1,:));
    mu_0_baseline_sbp_true = model_true.mu_0(3,1);
    subplot(n_plots,3,7);
    hold on;
    histogram(mu_0_baseline_sbp,11);
    title('mu_0^{(3)} Histogram');
    line([mu_0_baseline_sbp_true,mu_0_baseline_sbp_true],[0,length(samples)],'Color','r','LineWidth',1);
    subplot(n_plots,3,8);
    hold on;
    plot(1:length(samples),mu_0_baseline_sbp);
    chunk = floor(length(samples)/5);
    for i = 1:5
        chunk_mean = mean(mu_0_baseline_sbp((i-1)*chunk+1:i*chunk));
        line([(i-1)*chunk+1 i*chunk],[chunk_mean chunk_mean],'Color','r','LineWidth',1);
    end
    title('mu_0^{(3)} Chain Trace');
    subplot(n_plots,3,9);
    autocorr(mu_0_baseline_sbp,'NumLags',100);
    title('mu_0^{(3)} ACF');
    
    for n_d = 1:n_drugs
        Emax_nd = squeeze(Emax(n_d,1,:));
        Emax_nd_true = model_true.U.Emax(n_d,1);
        subplot(n_plots,3,9+6*(n_d-1)+1);
        hold on;
        histogram(Emax_nd,11);
        title(sprintf('E_{max}(%d) Histogram',n_d));
        line([Emax_nd_true,Emax_nd_true],[0,length(samples)],'Color','r','LineWidth',1);
        subplot(n_plots,3,9+6*(n_d-1)+2);
        hold on;
        plot(1:length(samples),Emax_nd);
        chunk = floor(length(samples)/5);
        for i = 1:5
            chunk_mean = mean(Emax_nd((i-1)*chunk+1:i*chunk));
            line([(i-1)*chunk+1 i*chunk],[chunk_mean chunk_mean],'Color','r','LineWidth',1);
        end
        title(sprintf('E_{max}(%d) Chain Trace',n_d));
        subplot(n_plots,3,9+6*(n_d-1)+3);
        autocorr(Emax_nd,'NumLags',100);
        title(sprintf('E_{max}(%d) ACF',n_d));

        EC50_nd = squeeze(EC50(n_d,1,:));
        EC50_nd_true = model_true.U.EC50(n_d,1);
        subplot(n_plots,3,9+6*(n_d-1)+4);
        hold on;
        histogram(EC50_nd,11);
        title(sprintf('EC_{50}(%d) Histogram',n_d));
        line([EC50_nd_true,EC50_nd_true],[0,length(samples)],'Color','r','LineWidth',1);
        subplot(n_plots,3,9+6*(n_d-1)+5);
        hold on;
        plot(1:length(samples),EC50_nd);
        chunk = floor(length(samples)/5);
        for i = 1:5
            chunk_mean = mean(EC50_nd((i-1)*chunk+1:i*chunk));
            line([(i-1)*chunk+1 i*chunk],[chunk_mean chunk_mean],'Color','r','LineWidth',1);
        end
        title(sprintf('EC_{50}(%d) Chain Trace',n_d));
        subplot(n_plots,3,9+6*(n_d-1)+6);
        autocorr(EC50_nd,'NumLags',100);
        title(sprintf('EC_{50}(%d) ACF',n_d));
    end
end