function [CRBr,CRBmu] = computeCRB(B,Bgradx,Bgrady,theta,phi,gamma,s,b)

[mxx,myy,mzz,mxy,mxz,myz] = angle2M(theta,phi,gamma);

I = (B(:,:,1)*mxx + B(:,:,2)*myy + B(:,:,3)*mzz + ...
    B(:,:,4)*mxy + B(:,:,5)*mxz + B(:,:,6)*myz)*s+b;
I(I<1e-15) = nan;
Idx = (Bgradx(:,:,1)*mxx + Bgradx(:,:,2)*myy + Bgradx(:,:,3)*mzz + ...
    Bgradx(:,:,4)*mxy + Bgradx(:,:,5)*mxz + Bgradx(:,:,6)*myz)*s;
Idy = (Bgrady(:,:,1)*mxx + Bgrady(:,:,2)*myy + Bgrady(:,:,3)*mzz + ...
    Bgrady(:,:,4)*mxy + Bgrady(:,:,5)*mxz + Bgrady(:,:,6)*myz)*s;
FIr(1,1) = nansum(nansum(Idx.^2./I));
FIr(2,2) = nansum(nansum(Idy.^2./I));
FIr(1,2) = nansum(nansum(Idx.*Idy./I));
FIr(2,1) = FIr(1,2);
CRBr = inv(FIr);

[mxx,myy,mzz,mxy,mxz,myz] = angle2M(theta,phi,gamma+1e-6);
Ig = (B(:,:,1)*mxx + B(:,:,2)*myy + B(:,:,3)*mzz + ...
    B(:,:,4)*mxy + B(:,:,5)*mxz + B(:,:,6)*myz)*s+b;

[mxx,myy,mzz,mxy,mxz,myz] = angle2M(theta+1,phi,gamma);
It = (B(:,:,1)*mxx + B(:,:,2)*myy + B(:,:,3)*mzz + ...
    B(:,:,4)*mxy + B(:,:,5)*mxz + B(:,:,6)*myz)*s+b;

[mxx,myy,mzz,mxy,mxz,myz] = angle2M(theta,phi+1,gamma);
Ip = (B(:,:,1)*mxx + B(:,:,2)*myy + B(:,:,3)*mzz + ...
    B(:,:,4)*mxy + B(:,:,5)*mxz + B(:,:,6)*myz)*s+b;

FImu(1,1) = nansum(nansum((It-I).^2./I));
FImu(2,2) = nansum(nansum((Ip-I).^2./I));
FImu(3,3) = nansum(nansum(((Ig-I)*1e6).^2./I));
FImu(1,2) = nansum(nansum((It-I).*(Ip-I)./I));
FImu(1,3) = nansum(nansum((It-I).*((Ig-I)*1e6)./I));
FImu(2,3) = nansum(nansum((Ip-I).*((Ig-I)*1e6)./I));
FImu(2,1) = FImu(1,2);
FImu(3,1) = FImu(1,3);
FImu(3,2) = FImu(2,3);
CRBmu = inv(FImu);

end


function [mxx,myy,mzz,mxy,mxz,myz] = angle2M(theta,phi,gamma)

mux = sind(theta)*cosd(phi);
muy = sind(theta)*sind(phi);
muz = cosd(theta);
mxx = mux^2*gamma+(1-gamma)/3; myy = muy^2*gamma+(1-gamma)/3; mzz = muz^2*gamma+(1-gamma)/3; 
mxy = mux*muy*gamma; mxz = mux*muz*gamma; myz = muy*muz*gamma;

end

