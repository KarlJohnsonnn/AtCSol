function f=Rate(Reak,tBegin)
switch Reak.Type
    
    case 'PHOTO'
        suntime=mod(tBegin/3600,24);% REAL(dp), PARAMETER :: SunRise=4.50_dp, SunSet=19.50_dp
        if suntime>=4.50 && suntime<=19.50
            sun=updatesun(suntime);
            Meff=1;
            k=Reak.con(1)*sun;
            f=k*Meff;
        else
            f=0;
        end
    case 'PHOTO2'
        suntime=mod(tBegin/3600,24);
        if suntime>=4.50 && suntime<=19.50
            sun=updatesun(suntime);
            Meff=1;
            k=Reak.con(1)*sun*sun;
            f=k*Meff;
        else
            f=0;
        end
    case 'PHOTO3'
        suntime=mod(tBegin/3600,24);
        if suntime>=4.50 && suntime<=19.50
            sun=updatesun(suntime);
            Meff=1;
            k=Reak.con(1)*sun*sun*sun;
            f=k*Meff;
        else
            f=0;
        end
    case 'PHOTO1'
        disp('type unknown')
    case 'CONST'
        f=Reak.con(1);
end
end

function sun=updatesun(Tlocal)
%calculation sun for SmallStratoKPP
SunRise=4.5;
SunSet=19.5;
Ttmp = (2*Tlocal-SunRise-SunSet) / (SunSet-SunRise);
if Ttmp>0
    Ttmp=Ttmp*Ttmp;
else
    Ttmp=-Ttmp*Ttmp;
end
sun = (1+cos(pi*Ttmp)) * 0.5;
end