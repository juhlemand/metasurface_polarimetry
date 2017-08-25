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

theta = linspace(0,2*pi,64);
Sinc=[ones(1,length(theta));cos(theta);sin(theta);zeros(1,length(theta))];

for i=1:length(theta)
     plot3(Sinc(2,i),Sinc(3,i),Sinc(4,i),'marker','o','MarkerSize',8,'MarkerFaceColor','b')
end
VLP=0.5*[1 -1 0 0;-1 1 0 0;0 0 0 0;0 0 0 0];
LP45=0.5*[1 0 1 0;0 0 0 0;1 0 1 0;0 0 0 0];
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

delta=[pi/4 pi/6 pi/3 2*pi/3];
phi=0;
farve=['g' 'm' 'y' 'k'];

for j=1:4
LinRetard=[1 0 0 0;...
    0 cos(2*phi)^2+sin(2*phi)^2*cos(delta(j)) sin(2*phi)*cos(2*phi)*(1-cos(delta(j))) -sin(2*phi)*sin(delta(j));...
    0 sin(2*phi)*cos(2*phi)*(1-cos(delta(j))) sin(2*phi)^2+cos(2*phi)^2*cos(delta(j))  cos(2*phi)*sin(delta(j));...
    0 sin(2*phi)*sin(delta(j)) -cos(2*phi)*sin(delta(j)) cos(delta(j))];

    Sout=LinDiattenuator2*LinRetard*Sinc;
     plot3(Sout(2,:)./Sout(1,:),Sout(3,:)./Sout(1,:),Sout(4,:)./Sout(1,:),'marker','o','MarkerSize',8,'MarkerFaceColor','r')
     Sout2=LinDiattenuator*LinDiattenuator2*LinRetard*Sinc;
     %Sout2(:,i)=LinDiattenuator*LinRetard*Sinc(:,i);
     plot3(Sout2(2,:)./Sout2(1,:),Sout2(3,:)./Sout2(1,:),Sout2(4,:)./Sout2(1,:),'marker','o','MarkerSize',8,'MarkerFaceColor',farve(j))

intafSout2(j,:)=Sout2(1,:);
end
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
figure
for j=1:4
    plot(theta*90/pi,intafSout2(j,:))
    if j==1 
        hold on 
    end
end
% for i=1:length(theta)
%      plot3(Sinc(2,i),Sinc(3,i),Sinc(4,i),'marker','o','MarkerSize',8,'MarkerFaceColor','b')
% end