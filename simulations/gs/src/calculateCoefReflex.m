function R_para = calculateCoefReflex(angle_rad,n1,n2)

cos_theta_i = cos(angle_rad);
cos_theta_t = sqrt(1 - (n1/n2)^2 .* (sin(angle_rad).^2));

R_para      = (n2.*cos_theta_i - n1.*cos_theta_t) ./ ...
              (n2.*cos_theta_i + n1.*cos_theta_t);
    
end