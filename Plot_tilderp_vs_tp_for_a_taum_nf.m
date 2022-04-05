%==========================================================================
% Plot tilde(rp) at n=nf
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

Tp_plot=min_Tp:dt:max_Tp;

if ~isempty(n_conv_f)
    
    tilda_rp_curve_true=k1/(mu*true_taum^2-2*sigma*true_taum+1).*(-w*exp(-Tp_plot/true_taum)+...
      ((mu*true_taum-sigma).*sin(w*Tp_plot)+w*cos(w*Tp_plot)).*exp(-sigma*Tp_plot));

    [max_tilda_rp_curve,index_max_tilda_rp]=max(tilda_rp_curve_true);
    
    taum_=taum_est_f(n_conv_f(1));
    tilda_rp_curve_nf=k1/(mu*taum_^2-2*sigma*taum_+1).*(-w*exp(-Tp_plot/taum_)+...
      ((mu*taum_-sigma).*sin(w*Tp_plot)+w*cos(w*Tp_plot)).*exp(-sigma*Tp_plot));
    
   
    
    tilda_rp_samples_nf=k1./(mu*taum_est_f(1:n_conv_f(1)).^2-2*sigma*taum_est_f(1:n_conv_f(1))+1).*...
        (((mu*taum_est_f(1:n_conv_f(1))-sigma).*sin(w*Tp_f(1:n_conv_f(1)))+...
        w.*cos(w*Tp_f(1:n_conv_f(1)))).*exp(-sigma*Tp_f(1:n_conv_f(1)))-...
        w.*exp(-Tp_f(1:n_conv_f(1))./taum_est_f(1:n_conv_f(1))));

    plot(Tp_plot,tilda_rp_curve_true,'k','LineWidth',1)
    
    plot(Tp_plot,tilda_rp_curve_nf,'--r','LineWidth',1) 
    plot(Tp_f(1:n_conv_f(1)),tilda_rp_samples_nf,'or',...
        'MarkerEdgeColor','r','MarkerSize',8);
    plot(Tp_plot(index_max_tilda_rp),max_tilda_rp_curve,'o',...
        'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',8);
end

% taum_=taum_est_f(end);
% tilda_rp_plot_nmax=k1/(mu*taum_.^2-2*sigma*taum_+1).*(-w*exp(-Tp_plot/taum_)+...
%       ((mu*taum_-sigma).*sin(w*Tp_plot)+w*cos(w*Tp_plot)).*exp(-sigma*Tp_plot));
% 
% plot(Tp_plot,tilda_rp_plot_nmax) 
% plot(Tp_f,rp_surf_tuning_nmax,'^b')


xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1e-6) 
xlabel('$t_p~ (\mu s)$','interpreter','latex')
ylabel('$\tilde{r}_p$','interpreter','latex')
ax=gca;
ax.FontName = 'Times New Roman';
ax.FontSize = fig_font_size;

legend('True', 'Estimation','Samples','Peak','Location','northwest','Box','off')

strmax = ['a'];
text(185e-6,.112,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', fig_font_size);

box on
 saveas(fig,sprintf('fig-tilda_rp_Tp_nf.fig'))
 saveas(fig,sprintf('fig-tilda_rp_Tp_nf.pdf'))
 saveas(fig,sprintf('fig-tilda_rp_Tp_nf.png'))
 saveas(gcf,'fig-tilda_rp_Tp_nf','epsc')
