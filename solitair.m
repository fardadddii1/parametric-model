clear all

R=6;                                          %stent original diameter
SL=47;                                        %stent total lenth 
SCN=4;                                        %Number of cells to be added to the original design  


b=2;                                         %each cell's width in original config (circumfrential direction)
c=5;                                         %each cell's length in original config (longitudinal)

a=b/20;                                      %mesh size horizontally--the number at denominator is the the numberof mesh used for each

u=0.05;                                      %a value for distance between ovelapping 
ua=-2/180*pi;                                %the angel by which overlapping sheets go into each other-- (+ goes in and - they don'nt over lap)

%%

marx1=[-2 -3 -4 -4 -3 -2 -1   1 2 3 4 4 3 2 1]';    %mesh of the stent, definig the design of cells -- circumferential direction--first coil
marz1=[ 4  4  4  4  4  4  4   4 4 4 4 4 4 4 4]';    %mesh of the stent, definig the design of cells -- longitudinal direction--first coil

marx2=[2 3 4 4 3 2 1   -1 -2 -3 -4 -4 -3 -2 -1]';       % mesh of the stent, definig the design of cells -- circumferential direction--second coil 
marz2=[4 4 4 4 4 4 4    4  4  4  4  4  4  4  4]';       %mesh of the stent, definig the design of cells -- longitudinal direction--second coil 

marxT1=[-1 -2 -2 -2 -1 -1 -1 -1  -1 0 0 1 1 1 1 1 1 1 1 1 1 1 1]';              %mesh of the stent, definig the design of cells -- longitudinal direction--tail pattern left
marzT1=[ 4  4  4  4  4  4  4  4   4 5 5 4 4 4 4 4 4 4 4 4 4 4 4]';              %mesh of the stent, definig the design of cells -- longitudinal direction--tail pattern left

marxT2=-[-1 -2 -2 -2 -1 -1 -1 -1  -1 0 0 1 1 1 1 1 1 1 1 1 1 1 1]';              %mesh of the stent, definig the design of cells -- longitudinal direction--tail pattern right
marzT2= [ 4  4  4  4  4  4  4  4   4 5 5 4 4 4 4 4 4 4 4 4 4 4 4]';              %mesh of the stent, definig the design of cells -- longitudinal direction--tail pattern right

numz=10;                                    %number of loops in longitudinal direction--leave it at 10 or so, bucause afterward an specific number will be chosen-let it be
numx=4;                                     %half number of wires

NC=numx*2;                                  %number of coils

LP=length(marz1)+1;

NOL1=[(2.5+SCN)*LP+1.5 (3.5+SCN)*LP+1 (3.5+SCN)*LP (3+SCN)*LP+1]+1;   %Number of cells of each coil of first type (the multiplication factor), each 16 nodes is one half sin wave (called cells here), each number of the matrix is for one coil (totally four coils here)

NOL2=[(3+SCN)*LP+1 (3+SCN)*LP+1 (3.5+SCN)*LP+1 (2.5+SCN)*LP+1]+1;   %Number of Loops of each coil of second type, each 16 nodes is one half sin wave (called loop here), each number is for one coil (totally four coils here)

SC=(sum(marz1)+4)*(3.5+SCN);                %nmber of nodes in working area
ST=sum(marzT2)+4;                           %number of nodes in nonworking area (tail)
STC=SC+ST+2*NC/4;                           %this is tne number of nodes from the tip of stent to its end--number 8 and second is because the coil the end of working elngth is at the middle coil and the tip of the stent is at the end of last coil, and coils get lower from middle to the end

%%
% az=c/64;                                  %mesh size vertically

ISL=SL+10;                                  % a large number for the length of sheet
az=SL/STC;
WL=SC*az;                                   %working length
TL=STC*az;                                  %total length
Z=0:az:ISL;                                 %Length coordinates

[kjl,zs]=size(Z);

%%
%other fomulations for for the profile section

%for example raduis at different profile sections, tapered

% r=10+t.^2;    %

%%
%%%%%%circular cross section

% for i=1:zs

r=R/2;
theta=2*asin(a/(2*r));        %the angle needed for each segment of dimentions axa

xc=0;
yc=0;                         %center of the cylinder

%profile generation

n=0;

for th = 0:theta:2*pi+ua
    
    %     theta=2*asin(a/(2*r));
    theta=a/(2*pi*r);
    r=r-u*theta;
    
    % th1=asin((R/r)*sin(th));
    n=n+1;
    xunit(1,n) = r * cos(th) + xc;
    yunit(1,n) = r * sin(th) + yc;
    
    if imag(xunit)+imag(yunit)~=0
        displaye('your design is not in reasanable range, please change the design paramters: R, L, b , c, u, or ua')
        return
    end
end
%%
%new meshing for a regular me
curl=[[xunit(1,:)]',[yunit(1,:)']];

d=size(curl);

for klo=1:d-1
    
    lement_length_X(klo)=norm(curl(klo+1,:)-curl(klo,:));
    
end

L=sum(lement_length_X);             %length of the coil

a=(L)/(NC*20);                      %refinign based on the length we need, to have the lengtth of new elements
a=(L+5*a)/(NC*20);

nn=round(L/a)+3;                    %Calculating number of segments needed
e=a/100;                            %an error for the segment's length
ny=1;

P(1,:)=curl(1,:);

for po=1:nn-1
    
    if ny>d(1,1)
        break
    end
    %a loop for choosing the segment which satisfies our desired shape
    while norm(P(po,:)-curl(ny,:))<a
        ny=ny+1;
        if ny>d(1,1)
            break
        end
    end
    if ny>d(1,1)
        break
    end
    for t=0:0.0001:1
        %Curve\line fitting using bycentric coordinates
        P(po+1,:)=(1-t)*P(po,:)+t*curl(ny,:);
        %testing element lenght
        element_length=norm(P(po+1,:)-P(po,:));
        P;
        %checkpoint for element length
        if element_length>a
            break
        end
    end
    
end

for klo=1:size(P)-1
    lement_length_P(klo)=norm(P(klo+1,:)-P(klo,:));
end

%orthogonal length
L_X=sum(lement_length_X);
%new segmentation length
L_P=sum(lement_length_P);
relative_error=(L_X-L_P)/L_X;           %relative error based on the length of coil

curl=P;
xunit1(1,:)=curl(:,1)';
yunit1(1,:)=curl(:,2)';

% end

[kjh,s]=size(xunit1);

X=repmat(xunit1,zs,1);
Y=repmat(yunit1,zs,1);
Z=repmat(Z,s,1)';

%flat plane--for testing purposes

% xf=1:a:L_P*1.1;
% sf=length(xf);
% X=repmat(xf,zs,1);
% Y=zeros(zs,sf);
% Z=repmat(Z,sf,1)';

% figure (1)
% surf(X,Y,Z,'LineStyle',':','EdgeColor',([0,0,0]),'FaceColor',([1,1,1]),'FaceAlpha',0.3)
% hold on

%%
%first coil

nx1(1)=20;                                          %initiation node on circumferential direction
nz1(1)=50;                                          %initiation node on longitudinal direction

smx=size(marx1);
smz=size(marz1);

for k=2:smx+1
    
    nx1(k)=nx1(k-1)+marx1(k-1);
    
end

for m=2:smz+1
    
    nz1(m)=nz1(m-1)+marz1(m-1);
    
end

nx1=repmat(nx1,1,numz);
nz11=nz1;
stp1=max(nz11)-min(nz11);

for ty=1:numz-1
    
    nz11=[[nz11(:)];[nz1(:)+stp1*ty+4*ty]]';
    
end
nx1=[nx1(1)+1,[nx1]];      %the last element to joint the wires
nz11=[nz11(1)-4,[nz11]];      %the last element to joint the wires

% for h=1:size(nz11')
%
%     if nz11(h)>zs
%         break
%     end
%     C1(h,:)=[X(nz11(h),nx1(h)),Y(nz11(h),nx1(h)),Z(nz11(h),nx1(h))];
%     figure (1)
%     scatter3(C1(h,1),C1(h,2),C1(h,3),'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[1 1 1])
%     hold on
%
% end

%%
%second coil

nx2(1)=22;          %nx2(1)-nx1(1)=2                                    %initiation node on circumferential direction
nz2(1)=50;          %nz2(1)=nx1(1)                                               %initiation node on longitudinal direction

smx=size(marx2);
smz=size(marz2);

for k=2:smx+1
    
    nx2(k)=nx2(k-1)+marx2(k-1);
    
end

for m=2:smz+1
    
    nz2(m)=nz2(m-1)+marz2(m-1);
    
end

nx2=repmat(nx2,1,numz);
nz22=nz2;
stp2=max(nz22)-min(nz22);

for ty=1:numz-1
    
    nz22=[[nz22(:)];[nz2(:)+stp2*ty+4*ty]]';
    
end

nx2=[nx2(1)-1,[nx2]];         %the last element to joint the wires
nz22=[nz22(1)-4,[nz22]];      %the last element to joint the wires

% for h=1:size(nz11')
%
%     if nz22(h)>zs
%         break
%     end
%     C2(h,:)=[X(nz22(h),nx2(h)),Y(nz22(h),nx2(h)),Z(nz22(h),nx2(h))];
%     figure (1)
%     scatter3(C2(h,1),C2(h,2),C2(h,3),'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[1 1 1])
%     hold on
%
% end
%%
%other coils in circumferential direction

NCX1(1,:)=nx1;

% NCX21(1,:)=nx21;

NCX2(1,:)=nx2;

% NCX1(2,:)=NCX1(1,:)+max(NCX2(1,:));
%
% NCX2(2,:)=NCX2(1,:)+max(NCX1(1,:));
%
% NCX1(3,:)=NCX1(2,:)+max(NCX2(2,:));
%
% NCX2(3,:)=NCX2(2,:)+max(NCX1(2,:));

for r=2:numx
    
    NCX1(r,:)=NCX1(r-1,:)+(2*(max(NCX2(r-1,:))-min(NCX2(r-1,:)))-2);    %the first coil reconstruction
    
    NCX2(r,:)=NCX2(r-1,:)+(2*(max(NCX1(r,:))-min(NCX1(r,:)))-2);        %the second coil reconstruction
    
end

NCZ1(1,:)=nz11;

NCZ2(1,:)=nz22;

for r=2:numx
    
    NCZ1(r,:)=NCZ1(r-1,:)+8;                                           %the first coil reconstruction
    
    NCZ2(r,:)=NCZ2(r-1,:)+8;                                           %the second coil reconstruction
    
end

mt=0;

for q=1:numx
    mt=mt+1;
    for h=1:NOL1(q)
        
        if NCZ1(q,h)>zs
            
            break
            
        end
        
        if NCX1(q,h)>s
            
            continue
            
        end
        
        CH1(h,:,q)=[X(NCZ1(q,h),NCX1(q,h)),Y(NCZ1(q,h),NCX1(q,h)),Z(NCZ1(q,h),NCX1(q,h))];
        
        
%         figure (1)
%         
%         scatter3(CH1(h,1,q),CH1(h,2,q),CH1(h,3,q),'filled','MarkerEdgeColor','k',...
%             'MarkerFaceColor',[0 0 0])
%         
%         hold on
        
        
    end
    
    EPX(mt)=NCX1(q,h);
    EPZ(mt)=NCZ1(q,h);
    
    EPX1(q)=NCX1(q,h);         %X index for end point of first coils set
    EPZ1(q)=NCZ1(q,h);         %Z index for end point of first coils set
    EPCH1(q,:)=[X(NCZ1(q,h),NCX1(q,h)),Y(NCZ1(q,h),NCX1(q,h)),Z(NCZ1(q,h),NCX1(q,h))];     %Endpoint for the
    
    %eEnd points depiction
    
%     figure (1)
%     
%     scatter3(EPCH1(q,1),EPCH1(q,2),EPCH1(q,3),'filled','MarkerEdgeColor','k',...
%             'MarkerFaceColor',[1 0 0])
%     hold on

    CH11=CH1(:,:,q);
    % writing first coils' nodes in excel files
    CH11( ~any(CH11,2), : ) = [];   %nonzero elements
    
    figure (1)
    
    plot3(CH11(:,1),CH11(:,2),CH11(:,3),'Color',uint8([17 17 17]),'linewidth',1)
    
    hold on
    filename = sprintf('coil1_%d.txt',q);
    save(filename,'CH11','-ascii')
    
    for h=1:NOL2(q)
        
        if NCZ2(q,h)>zs
            
            break
            
        end
        
        if NCX2(q,h)>s
            continue
        end
        
        CH2(h,:,q)=[X(NCZ2(q,h),NCX2(q,h)),Y(NCZ2(q,h),NCX2(q,h)),Z(NCZ2(q,h),NCX2(q,h))];
        
        
%         figure (1)
%         
%         scatter3(CH2(h,1,q),CH2(h,2,q),CH2(h,3,q),'filled','MarkerEdgeColor','k',...
%                     'MarkerFaceColor',[0 0 0])
%         
%         hold on
        
    end
    
    
    
    mt=mt+1;
    EPX(mt)=NCX2(q,h);
    EPZ(mt)=NCZ2(q,h);
    
    EPX2(q)=NCX2(q,h);       %X index for end point of the second coils set
    EPZ2(q)=NCZ2(q,h);       %Z index for end point of the second coils set
    EPCH2(q,:)=[X(NCZ2(q,h),NCX2(q,h)),Y(NCZ2(q,h),NCX2(q,h)),Z(NCZ2(q,h),NCX2(q,h))]; %Endpoint for the
    
    %Endpoints depiction
    
%     figure (1)
%     
%     scatter3(EPCH2(q,1),EPCH2(q,2),EPCH2(q,3),'filled','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[1 1 1])
%     
%     hold on
    
    %writing first coils' nodes in excel files
    
    CH22=CH2(:,:,q);
    CH22( ~any(CH22,2), : ) = [];   %nonzero elements
    
    figure (1)
    
    plot3(CH22(:,1),CH22(:,2),CH22(:,3),'Color',uint8([17 17 17]),'linewidth',1)
    
    hold on
    
    filename = sprintf('coil2_%d.txt',q);
    save(filename,'CH22','-ascii')
    
end


%%
%tail first part
nxT1(1)=EPX(length(EPX)/2+1);                                                    %initiation node on circumferential direction
nzT1(1)=EPZ(length(EPX)/2+1);                                                   %initiation node on longitudinal direction

smx=size(marxT1);
smz=size(marzT1);

for k=2:smx+1
    
    nxT1(k)=nxT1(k-1)+marxT1(k-1);
    
end

for m=2:smz+1
    
    nzT1(m)=nzT1(m-1)+marzT1(m-1);
    
end

for qw=1:size(nxT1')
    
    CT1(qw,:)=[X(nzT1(qw),nxT1(qw)),Y(nzT1(qw),nxT1(qw)),Z(nzT1(qw),nxT1(qw))];
%     
%         figure (1)
%     
%         scatter3(CT1(qw,1),CT1(qw,2),CT1(qw,3),'filled','MarkerEdgeColor','k',...
%             'MarkerFaceColor',[0 0 0])
%     
%         hold on
    
end



CT1( ~any(CT1,2), : ) = [];   %nonzero elements

figure (1)

plot3(CT1(:,1),CT1(:,2),CT1(:,3),'Color',uint8([17 17 17]),'linewidth',1)

hold on
filename = sprintf('coilT_1.txt');
save(filename,'CT1','-ascii')

%%
%tail second part
nxT2(1)=EPX(length(EPX)/2+1);                                                     %initiation node on circumferential direction
nzT2(1)=EPZ(length(EPX)/2+1);                                                     %initiation node on longitudinal direction

EPX_h=EPX(length(EPX)/2);                                                     %this point will be replaced afterward, just saving it for making the coil1_3 end
EPX_hn=EPX(length(EPX)/2+1);                                                  %this point will be replaced afterward, just saving it for making the coil1_3 end

EPZ_h=EPZ(length(EPZ)/2);                                                     %half end point--this point will be replaced afterward, just saving it for making the coil1_3 end
EPZ_hn=EPZ(length(EPZ)/2+1);                                                 %next half end point--this point will be replaced afterward, just saving it for making the coil1_3 end


smx=size(marxT2);
smz=size(marzT2);

for k=2:smx+1
    
    nxT2(k)=nxT2(k-1)+marxT2(k-1);
    
end

for m=2:smz+1
    
    nzT2(m)=nzT2(m-1)+marzT2(m-1);
    
end

for qw=1:length(nxT2)
    
    CT2(qw,:)=[X(nzT2(qw),nxT2(qw)),Y(nzT2(qw),nxT2(qw)),Z(nzT2(qw),nxT2(qw))];
    
%         figure (1)
%     
%         scatter3(CT2(qw,1),CT2(qw,2),CT2(qw,3),'filled','MarkerEdgeColor','k',...
%             'MarkerFaceColor',[0 0 0])
%         
%         hold on
    
end

CT2( ~any(CT2,2), : ) = [];   %nonzero elements

figure (1)

plot3(CT2(:,1),CT2(:,2),CT2(:,3),'Color',uint8([17 17 17]),'linewidth',1)

hold on

filename = sprintf('coilT_2.txt');

save(filename,'CT2','-ascii')

%%
%fifth coil end (coil1_3)
lep=length(EPZ)/2;

SOCl=(EPZ(lep+1)-EPZ(lep))/(EPX(lep+1)-EPX(lep));              %Slope of connecting wires: (Z_2-Z_1)/(X_2-X_1)
ZWMARl=round(4*SOCl);                            %Z direction march: X direction march multiplied by slope
%SOC(lep)=round(SOC(lep),1)                              %Rounding the slope to the first decimal

NOMl=fix((EPX(lep+1)-EPX(lep))/4);

NWZ1l(lep,1)=EPZ(lep);
NWX1l(lep,1)=EPX(lep);

wr=1;
CC11l(wr,:)=[X(NWZ1l(lep,1),NWX1l(lep,1)),Y(NWZ1l(lep,1),NWX1l(lep,1)),Z(NWZ1l(lep,1),NWX1l(lep,1))]; %Connecting coils for the first one

for tr=2:NOMl
    
    wr=wr+1;
    NWX1l(lep,tr)=NWX1l(lep,tr-1)+4;
    NWZ1l(lep,tr)=NWZ1l(lep,tr-1)+ZWMARl;
    %     NWX1(lep,tr)=ceil(NWX1(lep,tr-1)+4*SOC(lep));
    CC11l(wr,:)=[X(NWZ1l(lep,tr),NWX1l(lep,tr)),Y(NWZ1l(lep,tr),NWX1l(lep,tr)),Z(NWZ1l(lep,tr),NWX1l(lep,tr))]; %Connecting coils for the first one
    
end

lep=lep+1;
NWZ1l(lep,1)=EPZ(lep);
NWX1l(lep,1)=EPX(lep);

wr=wr+1;
CC11l(wr,:)=[X(NWZ1l(lep,1),NWX1l(lep,1)),Y(NWZ1l(lep,1),NWX1l(lep,1)),Z(NWZ1l(lep,1),NWX1l(lep,1))]; %Connecting coils for the first one

CC11l( ~any(CC11l,2), : ) = [];   %nonzero elements

% figure (1)
% scatter3(CC11l(:,1),CC11l(:,2),CC11l(:,3),'filled','MarkerEdgeColor','k',...
%     'MarkerFaceColor',[0 0 0])
% hold on

figure (1)

plot3(CC11l(:,1),CC11l(:,2),CC11l(:,3),'Color',uint8([17 17 17]),'linewidth',1)

hold on
filename = sprintf('coil1_3_e.txt');
save(filename,'CC11l','-ascii')

%%
%side coils

TJPX1=nxT1(16);    %these are the indices where the side coils meet the tail
TJPX2=nxT2(16);
TJPZ1=nzT1(16);
TJPZ2=nzT2(16);

EPX(length(EPX)/2)=TJPX1;
EPX(length(EPX)/2+1)=TJPX2;

EPZ(length(EPZ)/2)=TJPZ1;
EPZ(length(EPZ)/2+1)=TJPZ2;
% EPX=[EPX(1:3) nxT1(16) nxT2(16) EPX(6:8)];
% EPZ=[EPZ(1:3) nzT1(16) nzT2(16) EPZ(6:8)];

wr=1;
for q=1:length(EPX)-1
    
    SOC(q)=(EPZ(q+1)-EPZ(q))/(EPX(q+1)-EPX(q));          %Slope of connecting wires: (Z_2-Z_1)/(X_2-X_1)
    ZWMAR(q)=round(4*SOC(q));                            %Z direction march: X direction march multiplied by slope
    %SOC(q)=round(SOC(q),1)                              %Rounding the slope to the first decimal
    
    NOM=fix((EPX(q+1)-EPX(q))/4);
    
    NWZ1(q,1)=EPZ(q);
    NWX1(q,1)=EPX(q);
    
    wr=wr+1;
    CC11(wr,:)=[X(NWZ1(q,1),NWX1(q,1)),Y(NWZ1(q,1),NWX1(q,1)),Z(NWZ1(q,1),NWX1(q,1))]; %Connecting coils for the first one
    
    if q==length(EPX)/2
        brkp=wr;
        %         wr=wr+1;
        %         CC11(wr,:)=[X(NWZ1(q,1),NWX1(q,1)),Y(NWZ1(q,1),NWX1(q,1)),Z(NWZ1(q,1),NWX1(q,1))]; %Connecting coils for the first one
        continue
        
    end
    
    for tr=2:NOM
        
        wr=wr+1;
        NWX1(q,tr)=NWX1(q,tr-1)+4;
        NWZ1(q,tr)=NWZ1(q,tr-1)+ZWMAR(q);
        %     NWX1(q,tr)=ceil(NWX1(q,tr-1)+4*SOC(q));
        CC11(wr,:)=[X(NWZ1(q,tr),NWX1(q,tr)),Y(NWZ1(q,tr),NWX1(q,tr)),Z(NWZ1(q,tr),NWX1(q,tr))]; %Connecting coils for the first one
        
    end
    
end

%last node
q=q+1;

NWZ1(q,1)=EPZ(q);
NWX1(q,1)=EPX(q);

wr=wr+1;
CC11(wr,:)=[X(NWZ1(q,1),NWX1(q,1)),Y(NWZ1(q,1),NWX1(q,1)),Z(NWZ1(q,1),NWX1(q,1))]; %Connecting coils for the first one
%
%CC11=[CC1(:,:,1);CC1(:,:,2);CC1(:,:,3)];
CC22=CC11(brkp+1:length(CC11),:);
CC11(brkp+1:length(CC11),:)=[];

CC11( ~any(CC11,2), : ) = [];   %nonzero elements
CC22( ~any(CC22,2), : ) = [];   %nonzero elements

% figure (1)
% scatter3(CC11(:,1),CC11(:,2),CC11(:,3),'filled','MarkerEdgeColor','k',...
%     'MarkerFaceColor',[0 0 0])
% hold on

% figure (1)
% scatter3(CC22(:,1),CC22(:,2),CC22(:,3),'filled','MarkerEdgeColor','k',...
%     'MarkerFaceColor',[0 0 0])
% hold on

figure (1)

plot3(CC11(:,1),CC11(:,2),CC11(:,3),'Color',uint8([17 17 17]),'linewidth',1)

hold on

figure (1)

plot3(CC22(:,1),CC22(:,2),CC22(:,3),'Color',uint8([17 17 17]),'linewidth',1)

hold on

filename = sprintf('coilSi_1.txt');
save(filename,'CC11','-ascii')

filename = sprintf('coilSi_2.txt');
save(filename,'CC22','-ascii')
