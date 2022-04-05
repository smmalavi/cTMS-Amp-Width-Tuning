%==========================================================================
% Plot tilde(rp) at n=n_max
%==========================================================================
%
% Seyed Mohammad Mahdi Alavi+, Stellantis (Chrysler), Canada 
% Fidel Vila-Rodriguez, Unitverisyt of British Columbia, Canada 
% Adam Mahdi, University of Oxford, UK
% Stefan M. Goetz, University of Cambridge (UK), Duke University (USA)
% +: code written by
% e-mail: mahdi.alavi.work@gmail.com
%
% April 2022
%==========================================================================


fig=figure
hold on

dt=1e-6;
dVc=0.01;



tilda_rp_curve_true=k1/(mu*true_taum^2-2*sigma*true_taum+1).*(-w*exp(-Tp_plot/true_taum)+...
      ((mu*true_taum-sigma).*sin(w*Tp_plot)+w*cos(w*Tp_plot)).*exp(-sigma*Tp_plot));


taum_=taum_est_f(end);
tilda_rp_curve_nmax=k1/(mu*taum_.^2-2*sigma*taum_+1).*(-w*exp(-Tp_plot/taum_)+...
      ((mu*taum_-sigma).*sin(w*Tp_plot)+w*cos(w*Tp_plot)).*exp(-sigma*Tp_plot));

tilda_rp_samples_nmax=k1./(mu*taum_est_f(1:n_max).^2-2*sigma*taum_est_f(1:n_max)+1).*...
        (((mu*taum_est_f(1:n_max)-sigma).*sin(w*Tp_f(1:n_max))+...
        w.*cos(w*Tp_f(1:n_max))).*exp(-sigma*Tp_f(1:n_max))-...
        w.*exp(-Tp_f(1:n_max)./taum_est_f(1:n_max)));


plot(Tp_plot,tilda_rp_curve_true,'k','LineWidth',1)
        
plot(Tp_plot,tilda_rp_curve_nmax,'--r','LineWidth',1) 
plot(Tp_f,tilda_rp_samples_nmax,'or',...
        'MarkerEdgeColor','r','MarkerSize',8);

plot(Tp_plot(index_max_tilda_rp),max_tilda_rp_curve,'o',...
        'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',8);


xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1e-6) 
xlabel('$t_p~ (\mu s)$','interpreter','latex')
ylabel('$\tilde{r}_p$','interpreter','latex')
ax=gca;
ax.FontName = 'Times New Roman';
ax.FontSize = fig_font_size;
box on

strmax = ['b'];
text(185e-6,.112,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);

saveas(fig,sprintf('fig-tilda_rp_Tp_nmax.fig'))
 saveas(fig,sprintf('fig-tilda_rp_Tp_nmax.pdf'))
 saveas(fig,sprintf('fig-tilda_rp_Tp_nmax.png'))
 saveas(gcf,'fig-tilda_rp_Tp_nmax','epsc')
