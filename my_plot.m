function  my_plot(data)
%my_plot plot stress strain relations, etc.
    for i=1:data.fignum
        figure(i)
        plot(data.x(i,:),data.y(i,:),'linewidth',2), grid on;
        set(gca,'fontname', 'Times New Roman');
        xlabel(data.xlabel(i,:),'Interpreter','latex','fontsize',15);
        ylabel(data.ylabel(i,:),'Interpreter','latex','fontsize',15);
        title(data.title(i,:),'Interpreter','latex','fontsize',15);
    end
end

