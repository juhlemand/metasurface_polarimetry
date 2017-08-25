clear all
close all
figure
sphere
%surf(x,y,z,'HandleVisibility','off')
hold on
axis off
colormap white

%plot3(1,0,0,'marker','o','MarkerSize',10,'MarkerFaceColor','b')
%plot3(-1,0,0,'marker','x','MarkerSize',10,'LineWidth',4)
alpha(.5)
lims = axis;
lgd=1.3;
plot3(lgd*lims(1:2),[0 0],[0 0],'k','LineWidth',1.5,'HandleVisibility','off')
plot3([0 0],lgd*lims(1:2),[0 0],'k','LineWidth',1.5,'HandleVisibility','off')
plot3([0 0],[0 0],lgd*lims(1:2),'k','LineWidth',1.5,'HandleVisibility','off')
text(1.1,0,0.1,'s1','FontSize',12)
text(0,1.1,0.1,'s2','FontSize',12)
text(0.1,0,1.1,'s3','FontSize',12)
arrow3d([0 0],[0 0],[1 lgd],.5,0.005,0.03,'k');
arrow3d([0 0],[1 lgd],[0 0],.5,0.005,0.03,'k');
arrow3d([1 lgd],[0 0],[0 0],.5,0.005,0.03,'k');

m=1;
theta = linspace(0,2*pi,64);
for i=1:length(theta)
for j=1:length(theta)
    Si(:,m)=[1;cos(theta(i))*cos(theta(j));sin(theta(i))*cos(theta(j));sin(theta(j))];
    m=m+1;
end
end

%scatter3(Si(2,:),Si(3,:),Si(4,:),'marker','o','MarkerFaceColor','b')
r=0;
q=1;
v=pi/4;
LinDiattenuator=0.5*[q+r (q-r)*cos(2*v) (q-r)*sin(2*v) 0;...
    (q-r)*cos(2*v) (q+r)*cos(2*v)^2+2*sqrt(q*r)*sin(2*v)^2 (q+r-2*sqrt(q*r))*sin(2*v)*cos(2*v) 0;...
    (q-r)*sin(2*v) (q+r-2*sqrt(q*r))*sin(2*v)*cos(2*v) (q+r)*sin(2*v)^2+2*sqrt(q*r)*cos(2*v)^2 0;...
    0 0 0 2*sqrt(q*r)];
r=0.5;
q=1;
v=0;
LinDiattenuator2=0.5*[q+r (q-r)*cos(2*v) (q-r)*sin(2*v) 0;...
    (q-r)*cos(2*v) (q+r)*cos(2*v)^2+2*sqrt(q*r)*sin(2*v)^2 (q+r-2*sqrt(q*r))*sin(2*v)*cos(2*v) 0;...
    (q-r)*sin(2*v) (q+r-2*sqrt(q*r))*sin(2*v)*cos(2*v) (q+r)*sin(2*v)^2+2*sqrt(q*r)*cos(2*v)^2 0;...
    0 0 0 2*sqrt(q*r)];
farve=['g' 'm' 'y' 'k'];
delta=[pi/4 pi/6 pi/3 2*pi/3];

for j=1:4%4
phi=0;
LinRetard=[1 0 0 0;...
    0 cos(2*phi)^2+sin(2*phi)^2*cos(delta(j)) sin(2*phi)*cos(2*phi)*(1-cos(delta(j))) -sin(2*phi)*sin(delta(j));...
    0 sin(2*phi)*cos(2*phi)*(1-cos(delta(j))) sin(2*phi)^2+cos(2*phi)^2*cos(delta(j))  cos(2*phi)*sin(delta(j));...
    0 sin(2*phi)*sin(delta(j)) -cos(2*phi)*sin(delta(j)) cos(delta(j))];


Sout1=LinDiattenuator*LinDiattenuator2*LinRetard*Si;
Sout=LinDiattenuator*LinRetard*Si;

%scatter3(Si(2,:).*Sout(1,:),Si(3,:).*Sout(1,:),Si(4,:).*Sout(1,:),'marker','o','MarkerFaceColor',farve(j))
plot3(Si(2,:).*Sout1(1,:),Si(3,:).*Sout1(1,:),Si(4,:).*Sout1(1,:),'-')
end
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off