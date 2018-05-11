%% Final Project: Thickness
% close all;
% clear;
%% Quantitative Analysis on Ellipses
% Calculate thicknesses at known interval for each case
num_iter = 2000;
r1_minor = 100;
r0 = 50;

%% 1:1 Sample rate
for i = 1
    r1_major = i*100;
    n = r1_major+ceil(r1_major/10);
    x = -n:n;
    y = -n:n;
    [xx,yy] = meshgrid(x,y);
    
    u = ones(size(xx));
    u((xx.^2/r1_major.^2+yy.^2/r1_minor.^2)<1) = .5;
    u((xx.^2+yy.^2)<r0^2) = 0;
%     [u, L0, L1] = thickness(u, 3000);
    w = L0+L1;
    
    major_axis = w(n+1,:);
    major_axis(major_axis==0) = [];
    major_axis_l0 = L0(n+1,:);
    major_axis_l1 = L1(n+1,:);
    major_axis_l0(major_axis_l0 == 0) = [];
    major_axis_l1(major_axis_l1 == 0) = [];
    
    moving = movstd(w(n+1,:),r1_major-r0);
    std_dev_l0_a = movstd(major_axis_l0, r1_major-r0);
    std_dev_l1_a = movstd(major_axis_l1, r1_major-r0);
    
%     figure;
%     plot(major_axis_l0, std_dev_l0_a);
%     xlabel('Thickness from d_0R to d_1R');
%     ylabel('Standard Deviation from True');
%     title(['Moving Standard Deviation from (R_a = ' num2str(r1_major) ') - (R_b = ' num2str(r0) ') for L_0'])
%     figure;
%     plot(major_axis_l1, std_dev_l1_a);
%     xlabel('Thickness from d_1R to d_0R');
%     ylabel('Standard Deviation from True');
%     title(['Moving Standard Deviation from (R_a = ' num2str(r1_major) ') - (R_b = ' num2str(r0) ') for L_1'])

%     major_b(i,:) = major_axis;
    std_dev_a(i) = std(moving);
    major_means_a(i) = mean(major_axis-1);
    eccentricity_a(i) = sqrt(r1_major.^2+r1_minor.^2);
end

%% 2:1 Sample rate
% std_dev_l0_b = zeros(10,
for i = 1
    r1_major = i*100;
    n = r1_major+ceil(r1_major/10);
    x = -n:.5:n;
    y = -n:.5:n;
    [xx,yy] = meshgrid(x,y);
    
    u = ones(size(xx));
    u((xx.^2/r1_major.^2+yy.^2/r1_minor.^2)<1) = .5;
    u((xx.^2+yy.^2)<r0^2) = 0;
    [u, L0, L1] = thickness(u, 3000);
    w = L0+L1;
    
    major_axis = w(2*n+1,:);
    major_axis(major_axis==0) = [];
    major_axis_l0 = L0(2*n+1,:);
    major_axis_l1 = L1(2*n+1,:);
    major_axis_l0(major_axis_l0 == 0) = [];
    major_axis_l1(major_axis_l1 == 0) = [];
    
    moving = movstd(w(2*n+1,:),r1_major-r0);
    std_dev_l0_b = movstd(major_axis_l0, r1_major-r0);
%     figure;
%     plot(major_axis_l0, std_dev_l0_b);
%     xlabel('Points from d_0R to d_1R');
%     ylabel('Standard Deviation from True');
%     title(['Moving Standard Deviation from (R_a = ' num2str(r1_major) ') - (R_b = ' num2str(r0) ') for L_0'])

%     major_b(i,:) = major_axis;
    std_dev_b(i) = std(moving);
    major_means_b(i) = mean(major_axis-1)/2;
    eccentricity_b(i) = sqrt(r1_major.^2+r1_minor.^2);
end


