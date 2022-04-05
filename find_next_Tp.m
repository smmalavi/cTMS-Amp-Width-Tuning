%==========================================================================
% Find next pulse width
% depolorization factor maximization
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



function q = find_next_Tp(z,tp_f_1tonm1,taum_nm1,g_normal, k1, mu, sigma, w)
    

%     % FIM maximization 
%     Sinm1=k1/(mu*taum_nm1^2-2*sigma*taum_nm1+1)*(w/taum_nm1.*exp(-tp_f_1tonm1/taum_nm1)+...
%         ((mu*taum_nm1-2*sigma)*w.*cos(w*tp_f_1tonm1)+...
%         (-w^2+sigma^2-sigma*mu*taum_nm1).*sin(w*tp_f_1tonm1)).*exp(-sigma*tp_f_1tonm1));
%     FIMnm1=dot(Sinm1,Sinm1);
%     
%     Sn=k1/(mu*taum_nm1^2-2*sigma*taum_nm1+1)*(w/taum_nm1*exp(-z/taum_nm1)+...
%         ((mu*taum_nm1-2*sigma)*w*cos(w*z)+(-w^2+sigma^2-sigma*mu*taum_nm1)*sin(w*z))*exp(-sigma*z));
%     FIMn=dot(Sn,Sn);
%     
%     FIM=FIMnm1+FIMn;
%     q=-real(FIM);

    
%     % rp maximization
    q=-k1/(mu*taum_nm1^2-2*sigma*taum_nm1+1)*(-w*exp(-z/taum_nm1)+...
        ((mu*taum_nm1-sigma)*sin(w*z)+w*cos(w*z))*exp(-sigma*z));
   
    
end
    
        
 
