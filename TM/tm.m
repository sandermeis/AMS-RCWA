clear
close all
clc

materials=["GaAs","Substrate","Air","Glass","InGaP","Ag","Au","Al03GaAs","MgF2","ZnS"];


for i=1:length(materials)
    FileName=strcat(["refractive indices.xlsx - "],materials{i},".csv");
  x{i} = csvread(FileName(1,1));
end

n=["Air","MgF2","ZnS","Al03GaAs","GaAs","InGaP","Al03GaAs","Ag","Air"];

d=[1E10,94,44,25,300,100,100,100,1E10];

d_total=sum(d(2):d(end-1));
x_grid=0:0.1:d_total;
find(x_grid==d(2))

for theta=(0:2)/6*pi

%n={"Air","MgF2","ZnS","GaAs","Ag","Air"};
%
%d=[1E10,94,44,300,100,1E10];

%n={"Air","MgF2","GaAs","Ag","Air"};
%
%d=[1E10,94,300,100,1E10];

len_n=length(n);

labda_1=300;
labda_2=1000;

%Absorption=zeros(labda_2-labda_1+1,len_n);
E_L=zeros(2,len_n,labda_2-labda_1+1);
E_R=zeros(2,len_n,labda_2-labda_1+1);

for labda=labda_1:labda_2
  index=labda-labda_1+1;
  
  if testtrue(materials,n,x,labda)
     
    n_0=zeros(1,len_n);
    for i=1:len_n
      n_0(i)=getdata(x,materials,n(i),labda);
      
      sin_th(i)=n_0(1)*sin(theta)/n_0(i);
      cos_th(i)=sqrt(n_0(i)^2-(n_0(1)*sin(theta))^2)/n_0(i);     
    end
      
      n_star=n_0.*cos_th;
      n_ast=n_0.*sin_th;
      
      k=2.*pi./labda.*n_star;

      delta=d.*k;
      
      M_total=I(n_star,1,2);
      for j=2:(len_n-1)
      M_total=M_total*L(delta,j)*I(n_star,j,j+1);
      end
    
      t_=1./M_total(1,1);
      r_=M_total(2,1)./M_total(1,1);

%       Absorption(index,len_n)=Poynting(n_star, n_ast,len_n,[t_;0]);
%       Absorption(index,len_n-1)=r_.*conj(r_);
%       
%       E_L=[t_;0];
%       for j=(len_n-1):-1:2
%         E_R=I(n_star,j,j+1)*E_L;
%         E_L=L(delta,j)*E_R;
%         
%         S=Poynting(n_star,n_ast,j,E_L)-Poynting(n_star, n_ast,j,E_R);
%         Absorption(index,j-1)=S;
%       end

        E_L(:,len_n,index)=[t_;0];
        for j=(len_n-1):-1:2
            E_R(:,j,index)=I(n_star,j,j+1)*E_L(:,j+1,index);
            E_L(:,j,index)=L(delta,j)*E_R(:,j,index);
        end
        E_R(:,1,index)=[1;r_];
  end
end
S1=Poynting(n_star,n_ast,E_L);
S2=Poynting(n_star, n_ast,E_R);
S=S1-S2;
S(1,:,:)=E_R(2,1,:).*conj(E_R(2,1,:));
S(len_n,:,:)=E_L(1,len_n,:).*conj(E_L(1,len_n,:));
plot_TM(n,labda_1,labda_2,S)

end

function plot_TM(n,labda_1,labda_2,Absorption)
len_n=length(n);
x_axis=[labda_1:labda_2]'*ones(1,len_n);

%Swap columns
Absorption=Absorption';

Absorption=column_to_back(Absorption,1);
n_legend=column_to_back(n,1);

lg=find(strcmp(n, 'GaAs'));
Absorption_disp=column_to_front(Absorption,lg);
n_legend=column_to_front(n_legend,lg);

figure
plot(x_axis,Absorption_disp)
legend(n_legend,'Location','eastoutside');
ylim([-0.1 1.1])
xlim([1.1*labda_1-0.1*(labda_2+1) 1.1*labda_2+0.1*(-labda_1+1)])
xlabel("Wavelength (nm)")
ylabel("Fractional Intensity (a.u.)")
xticks(0:100:labda_2)
yticks(-0.1:0.1:1.1)
end


function c=column_to_front(input,column)
    len_n=size(input,2);
    c=input(:,[column,1:(column-1),(column+1):len_n]);
end

function c=column_to_back(input,column)
    len_n=size(input,2);
    c=input(:,[1:(column-1),(column+1):len_n,column]);
end

function n=getdata(x,materials,name,lab)

index = find(strcmp(materials, name));
ind=find(x{1,index}(:,1)==lab);

n=x{1,index}(ind,2)+1i*x{1,index}(ind,3);

end


function [refl]=r(n,i,j)
  refl=(n(i)-n(j))/(n(i)+n(j));
end


function [trans]=t(n,i,j)
  trans=2*n(i)/(n(i)+n(j));
end


function [l]=L(delta,n)
  l=[exp(-1i*delta(n)),0;0,exp(1i*delta(n))];
end


function [i]=I(n_0,i,j)
  i=[1,r(n_0,i,j);r(n_0,i,j),1]*1/t(n_0,i,j);
end


function bool_found=testtrue(materials,n0,x,lab)
  bool_found=1;
  for i=1:length(n0)
    index = find(strcmp(materials, n0{i}));
    temp=x{1,index}(:,1);
    bool_found=bool_found*~isempty(find(temp==lab));
  end
end


function Sz=Poynting(n_0,n_ast,E)
  Ef=zeros(size(E));
  Ef(3,:,:)=0;
  Ef(1,:,:)=E(1,:,:);
  
  Eb=zeros(size(E));
  Eb(3,:,:)=0;
  Eb(1,:,:)=E(2,:,:);
  
  uf=zeros(3,size(n_0,2));
  uf(2,:)=n_ast;
  uf(3,:)=n_0;
  
  ub=zeros(3,size(n_0,2));
  ub(2,:)=n_ast;
  ub(3,:)=-n_0;
  B=multicross(uf,Ef)+multicross(ub,Eb);
  
  S=real(multicross(Ef+Eb,conj(B)))/n_0(1);
  Sz=squeeze(S(3,:,:));
%   S1=real(conj(E(1)+E(2))*n_0(i)*(E(1)-E(2)))/real(n_0(1));
%   S2=real((E(1)+E(2))*conj(n_0(i))*conj(E(1)-E(2)))/real(conj(n_0(1)));
%   S=(S1+S2)/2;

end

function C=multicross(A,B)
    C=zeros(size(B));
    C(1,:,:)=A(2,:,:).*B(3,:,:)-A(3,:,:).*B(2,:,:);
    C(2,:,:)=A(3,:,:).*B(1,:,:)-A(1,:,:).*B(3,:,:);
    C(3,:,:)=A(1,:,:).*B(2,:,:)-A(2,:,:).*B(1,:,:);
end