function qqt = QQ1d_se(u1,v1,w,a,b)

% u1=[0:0.1:1]
% v1=[0:0.1:1]
% w=0.2
% a=0

if (size(u1,1)<size(u1,2))
    u1 = u1';
end

if (size(v1,1)<size(v1,2))
    v1 = v1';
end

u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);


% qqt = 0.5.*(u-a).*erf((u-a)/(2.*w)) ...
%     + w/sqrt(sym(pi)).*exp(-(u-a).^2/(4.*w.^2))...
%     - 0.5.*(v-u).*erf((v-u)/(2.*w)) ...
%     - w/sqrt(sym(pi)).*exp(-(v-u).^2/(4.*w.^2))...
%     + 0.5.*(v-a).*erf((v-a)/(2.*w)) ...
%     + w/sqrt(sym(pi)).*exp(-(v-a).^2/(4.*w.^2))...
%     - w/sqrt(sym(pi));

qqt = pi.*w.^2.*(u-a).*erf((u-a)/(2.*w)) ...
    + 2.*sqrt(pi).*w.^3.*exp(-(u-a).^2/(4.*w.^2))...
    - pi.*w.^2.*(v-u).*erf((v-u)/(2.*w)) ...
    - 2.*sqrt(pi).*w.^3.*exp(-(v-u).^2/(4.*w.^2))...
    + pi.*w.^2.*(v-a).*erf((v-a)/(2.*w)) ...
    + 2.*sqrt(pi).*w.^3.*exp(-(v-a).^2/(4.*w.^2))...
    - 2.*sqrt(pi).*w.^3;


end
