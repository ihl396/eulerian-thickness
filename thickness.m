% u must be a square matrix
% Set region of interest to 0.5 and other region potential from 0 to 1
% max_win is for quicker computation window (usually largest foci)
% floor(n/2)-max_win+2:ceil(n/2)+max_win

function [u, L0, L1] = thickness(u, max_iter) 
u_prev = u;
input_u = u;
n = length(u);
% field_e = [];
n_iter = 0;
delta_t = .25;
% if max_win <= 0
%     max_win = n;
% end

if max_iter <=0
    max_iter = 5000;
end

% Compute laplace(U) = 0

while n_iter < max_iter
    n_iter = n_iter+1;
    for i = 2:n-1
        for j = 2:n-1
            if input_u(i,j) == .5
                u_prev(i,j) = u(i,j)+ delta_t*(u(i+1,j) + u(i-1,j) + u(i,j+1) +u(i,j-1) - 4*u(i,j));
            end
        end
    end

    u = u_prev;
%     grad_phi_x = (phi(2:end,:)-phi(1:end-1,:))./2;
%     grad_phi_y = (phi(:,2:end)-phi(:,1:end-1))./2;
%     field_e(n_iter) = sqrt(sum(sum(grad_phi_x.^2))+sum(sum(grad_phi_y.^2)));
%     if n_iter == 1
%         e_ratio = 1;
%     else
%         e_ratio = (field_e(n_iter-1)-field_e(n_iter))/2;
%     end
%     fprintf('%8.6f | %8.6f\n', e_ratio, field_e(n_iter))
end

% Compute Tangent
grad_phi_x = zeros(size(u));
grad_phi_y = zeros(size(u));
T_x = zeros(size(u));
T_y = zeros(size(u));
for i = 2:n-1
    for j = 2:n-1
        if input_u(i,j) == .5
            grad_phi_x(i,j) = .5*(u_prev(i+1,j)-u_prev(i-1,j));
            grad_phi_y(i,j) = .5*(u_prev(i,j+1)-u_prev(i,j-1));
            norm_u = sqrt(grad_phi_x(i,j).^2+grad_phi_y(i,j).^2);
            T_x(i,j) = grad_phi_x(i,j)./norm_u;
            T_y(i,j) = grad_phi_y(i,j)./norm_u;
        end
    end
end

L0 = zeros(size(u));
L1 = zeros(size(u));
dL0 = L0;
dL1 = L1;
iter = 0;
while iter < max_iter
    iter = iter + 1;
    for i = 2:n-1
        for j = 2:n-1
            if input_u(i,j) == .5
                % upwind differencing
                if -T_x(i,j) < 0
                    L0_x = L0(i-1,j);
                    L1_x = L1(i+1,j);
                else
                    L0_x = L0(i+1,j);
                    L1_x = L1(i-1,j);
                end
                if -T_y(i,j) < 0
                    L0_y = L0(i,j-1);
                    L1_y = L1(i,j+1);
                else
                    L0_y = L0(i,j+1);
                    L1_y = L1(i,j-1);
                end
                dL0(i,j) = (1 + abs(T_x(i,j))*L0_x + abs(T_y(i,j))*L0_y)/(abs(T_x(i,j))+abs(T_y(i,j)));
                dL1(i,j) = (1 + abs(T_x(i,j))*L1_x + abs(T_y(i,j))*L1_y)/(abs(T_x(i,j))+abs(T_y(i,j)));
            end
        end
    end
    L0 = dL0;
    L1 = dL1;
end

end