figure(1)
plot(epsilon(t,1),stress(t,1),epsilon(t,3),stress(t,3),'--','linewidth',2), grid on;
set(gca,'fontname', 'Times New Roman');
xlabel('$\epsilon_1/\epsilon_3$','Interpreter','latex','fontsize',15);
ylabel('$\sigma_1/\sigma_3$/Pa','Interpreter','latex','fontsize',15);
title('Stress-strain relations','fontsize',15);
legend('\fontsize{15} \sigma_1','\fontsize{15} \sigma_3')

figure(2)
plot(epsilon(t,1),sqrt_J2(t),'-','color','b','linewidth',2), grid on;
set(gca,'fontname', 'Times New Roman');
xlabel('\epsilon_1','fontsize',15);
ylabel('$\sqrt{J_2}$','Interpreter','latex','fontsize',15);
title('$J_2-\epsilon_1$ relation', 'Interpreter','latex','fontsize',15');

figure(3)
plot(p(t),sqrt_J2(t),'linewidth',2), grid on;
set(gca,'fontname', 'Times New Roman');
xlabel('$p$','Interpreter','latex','fontsize',15);
ylabel('$\sqrt{J_2}$','Interpreter','latex','fontsize',15);
title('Loading path','fontsize',15);
