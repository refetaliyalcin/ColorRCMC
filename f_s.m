function [sonuc] = f_s(s)
    if s>((24/116)^3)
       sonuc=s^(1/3);
    else
       sonuc=841*s/108+16/116;
    end
end

