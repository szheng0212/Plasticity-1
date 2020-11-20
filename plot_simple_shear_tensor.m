figure(1)
plot(squeeze(epsilon(1,3,t)),squeeze(stress(1,3,t)),'linewidth',2), grid on;
set(gca,'fontname', 'Times New Roman');
xlabel('$\epsilon_s$','Interpreter','latex','fontsize',15);
ylabel('$\sigma_s$/Pa','Interpreter','latex','fontsize',15);
title('Stress train relation','fontsize',15);

figure(2)
plot(squeeze(epsilon(1,3,t)),sqrt_J2(t),'-','color','b','linewidth',2), grid on;
set(gca,'fontname', 'Times New Roman');
xlabel('\epsilon_s','fontsize',15);
ylabel('$\sqrt{J_2}$','Interpreter','latex','fontsize',15);
title('$J_2-\epsilon_s$ relation', 'Interpreter','latex','fontsize',15);
%xticks(-0.5:0.1:0);

figure(3)
plot(p(t),sqrt_J2(t),'linewidth',2), grid on;
set(gca,'fontname', 'Times New Roman');
xlabel('$p$','Interpreter','latex','fontsize',15);
ylabel('$\sqrt{J_2}$','Interpreter','latex','fontsize',15);
title('Loading path','fontsize',15);
