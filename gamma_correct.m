function [sonuc] = gamma_correct(u)
    if u<=0.0031308
       sonuc=323*u/25;
    else
       sonuc=(211*u^(5/12)-11)/200;
    end
end

