clear all; close all; clc
%%
                                                    %%%%%%%%%%%%%������������%%%%%%%%%%%%
N=10;          %��ֵ�źŵ���Ŀ
fs=2000;        %������
f=60;           %Ƶ��

A=100;          %���

sparsity=0.8;   %ϡ����
SNR=5;         %����ȣ����������ȡֵ1-100���˴���������sigma(��������)
buffer=16;    %������ȷ������˲ʱ�źŲ���˴��ص�
%%
                                                            %ģ����ͨ��ֵ�ź�
Lf=fs/(2*f)+1;
n=0:Lf-1;n=[n,Lf-1];%n��ȡֵ0,1,2...Lf-1  
Pf=A*sin(pi*n/(Lf-1));%sin(wx)����w=(2*pi/(Lf-1)),x=n/2   ����Ϊ16.666666ms����������Ϊ8.333   ���Զ�Ӧʱ��Ϊn/2
%%
                                                    %���˹��������׼�����    �õ���˹������
M=Lf;             %�������źŵĳ��ȣ������ڵ�����������
s1=0;
s2=Lf-1; 
s=Pf;      %%%ȡΪ���ź�%%%
S=sum(s.^2);%����������Ԫ��ƽ����
sigma=sqrt(S/(SNR*(s2-s1+1)));%����sigma�����ڲ���������
L=fix((M/(1-sparsity)*N)+Lf-1);            %����x����������
w=sqrt(sigma.^2)*randn(1,L);%��������Ϊsigmaƽ������ΪL�ĸ�˹������


%%%%%%%%%%%%��˹�����������ֵ���Լ�������������ݵĴ�С
max_w=max(abs(w));%%�����ķ�ֵΪ100��������ȡ��������Ϊ
fault=max_w+0.618*A;
%%%%%%%%%%%%

    w1=w(2:2:end);%��Ӧʱ��Ϊ�������ź�ȡֵ

t_wPf=0.5:0.5:L/2;%������ʱ�䷶Χ                           ��Ӧʱ��Ϊn/2
    T_wPf=t_wPf;
    t_wPf1=1:L/2;%��Ӧ������ʱ��
%%
                                            %���ź���ѡ����λ�ã�buffer��ʹ���Ƿ�ֹ���ڵ���������ֵ���
i=1:N;
a=1+(i-1)*M/(1-sparsity)+buffer;b=i*M/(1-sparsity);%%%%%%��[a,b]�������ѡȡ���źŵ�λ��
%(b-a)*rand(1)+a
tau_i=(b-a).*rand(1,N)+a;    %���������ѡ��õ�ÿ�����źŵ�����tau_i;
tau_i1=repmat(tau_i',1,length(n));%�õ�����ͬ�ľ���
tau_i=round(tau_i);%��λ�ý�����������ȡ��������Ϊ�±�
%%
                                                            %%%%%%%%%%%%�ϳ��ź� 
l=zeros(N,length(n));%��֪l����Ĵ�С�����Դ��ڴ���ֱ��ȡ�ô�С����ֹ��loop��ÿ�ζ��ڴ��ȡ��
k=0:Lf-1;
k=[k,Lf];
for j=1:N
    k1=repmat(k,N,1);
    l(j,:)=k1(j,:)+tau_i1(j,:);%�ܹ�10���壬ÿ������Lf=18�������ʾ��ÿ�����λ���ɵ�һ��λ��ȷ��,�õ�10*18 ���������
end
                                                        %%%%%%%%%%%%�õ����õ������ź�
PM=zeros(N,L);%��֪PM�Ĵ�С
for d=1:N;
    Fz=zeros(1,tau_i(d)-2);%���ڱ����0��ʼ����Ҫ�ȷ�λ��С1�����Դ˴���2����λ��ǰȫ����0
    Bz=zeros(1,L-tau_i(d)+2-length(Pf));%��֮��λ��ȫ��0
    P=[Fz,Pf,Bz];
    PM(d,:)=P;     %�õ���ľ����ʾ
end
%%
                                                                           %%%%%%%%%�Եõ��ķ��ź��Լ������źŽ��е���
sum_P=sum(PM,1);%���к͵õ�һ��������
sum_P1=sum_P(2:2:end);
[pks1,locs1]=findpeaks(sum_P);%ÿ���߶���Ѱ�����м���ֵ������
%%%%%%%%%%%%%%%%%%%�ҵ��������źŵ�λ�ã��������ѡȡ�����źŵ�λ��
not_s=find(sum_P==0);
ran_int=randperm(length(not_s),1);%��1:length(not_s)��Χ�ڲ���һ���������
fault_loc=not_s(ran_int);%��ù����źŵ�λ��
fault_loc1=1/2*fault_loc;
%%%%%%%%%%%%%%%%%%%%
locs2=t_wPf(locs1);%ԭʼ�����з��Ӧ��ʱ��

Pin=zeros(1,length(sum_P));%�����������ź�pin
Pin(fault_loc)=fault;
x=sum_P+Pin+w;%%%%%%%%%%%%%%%%%%%%%%%%ʱ��Ϊn/2
Xx=x;
ES=Xx(locs1);%ԭʼ�������źŷ��Ӧ����Ч�źŵķ�ֵ
x1=x(2:2:end);%��Ӧʱ��Ϊ�������ź�ȡֵ
 %%
figure
plot(t_wPf,sum_P,'--','LineWidth',1.5);
hold on;
plot(t_wPf,x,'color',[.7 .7 .7]);
hold on;
plot(fault_loc1,0,'gx');
hold on;
plot(locs2,0,'rx');
hold on;
xlabel('time(ms)');
ylabel('amplitude(uV)');
legend('original signal','noisy data','pulse interference','detected peaks');
hold off;



%%
%%%%%%%%%С���任
Pj=zeros(5,length(x));
time=zeros(5,length(x));
len=zeros(1,length(x));
a_1=[sqrt(1/2),sqrt(1/2)];b_1=[-sqrt(1/2),sqrt(1/2)];
c=0.382;
for im=1:5
    if mod(L,2)==1
        x=[x,0];%�źų���Ϊ����ʱ��ĩβ�����Ϊż��λ
        t_wPf=[t_wPf,max(t_wPf)+0.5*im];
    end
    X=x;
    T=t_wPf(2:2:end);
%%
%�ֽ����
    A_1=a_1;B_1=b_1;
    for ki=1:(length(x)*1/2-1)
        A_1=blkdiag(A_1,a_1);
        B_1=blkdiag(B_1,b_1);
    end
    AB=[A_1',B_1']';%�õ��ֽ����
%%
    %�õ�ƽ����ϸ���ź�
    X_sd=(AB*X')';

    X_s=X_sd(1:length(X_sd)*1/2);%ƽ���ź�
    X_d=X_sd(length(X_sd)*1/2+1:length(X_sd));%ϸ���ź�
    if (im==1)
        X_s1=X_s;%�ֽ�һ�ε�ƽ���ź�
        X_d1=X_d;%�ֽ�һ�ε�ϸ���ź�
        T1=T;%�ֽ�һ��ƽ���źŶ�Ӧ��ʱ��
    end
    if (im==2)
        X_s2=X_s;%�ֽ�2�ε�ƽ���ź�
        X_d2=X_d;%�ֽ�2�ε�ϸ���ź�
        T2=T;%�ֽ�2��ƽ���źŶ�Ӧ��ʱ��   
    end
    if (im==3)
        X_s3=X_s;
        X_d3=X_d;
        T3=T;
    end
    if(im==4)
        X_s4=X_s;
        X_d4=X_d;
        T4=T;
    end
    if(im==5)
        X_s5=X_s;
        X_d5=X_d;
        T5=T;
    end
%% 
%���źŽ����ع�
    X_1=[X_s,zeros(size(X_s))]*AB;%ƽ���ź��ع�һ��
    t_wPf_1=t_wPf;%��Ӧ��ʱ��
    if im==1
        X_r10=X_1;
        T_r10=t_wPf_1;
        %%%%����ֵ0.618max(X_r)%%%%%%%%%%%��һ����ֵ%%%%%
        gamma=c*max(X_r10);
        %%%%
        X_r10(X_r10<gamma)=0;
%         [pks10,locs10]=findpeaks(X_r10);%��ֵ
%         locs10=find(X_r10>0);%�ҵ�����0�ĵ������
%         T_r100=T_r10(locs10);%�ҵ���Ӧ��ʱ��
        T_r100=X_r10&T_r10;
    end
    
    if im==2%���ηֽ�֮��ƽ���źŵ��ع�
        X_r21=X_1;
%         T_r21=t_wPf_1;
        %%%���ع�����
        A_20=a_1;B_20=b_1;
        for ki=1:(length(X_r21)-1)
            A_20=blkdiag(A_20,a_1);
            B_20=blkdiag(B_20,b_1);
        end
        AB_20=[A_20',B_20']';%�õ��ع�����
        X_r20=[X_r21,zeros(size(X_r21))]*AB_20;
        T_r20=linspace(0.5,length(X_r21),length(X_r20));
        %%%%����ֵ0.618max(X_r)
        gamma=c*max(X_r20);
        %%%%%%%%%%%%%%%%%%%%%%
        xx=X_r20;
        xx(xx<0)=0;
        %%%%%%%%%%%%%%%%%%%%%%
        X_r20(X_r20<gamma)=0;
%         [pks20,locs20]=findpeaks(X_r20);%��ֵ
%         locs20=find(X_r20>0);%�ҵ�����0�ĵ������
%         T_r200=T_r20(locs20);%�ҵ���Ӧ��ʱ��
        T_r200=X_r20&T_r20;
    end
    
    if im==3%���ηֽ�֮��ƽ���źŵ��ع�
        X_r32=X_1;
%         T_r32=t_wPf_1;
        %%%���ع�����
        A_31=a_1;B_31=b_1;
        for ki=1:(length(X_r32)-1)
            A_31=blkdiag(A_31,a_1);
            B_31=blkdiag(B_31,b_1);
        end
        AB_31=[A_31',B_31']';%�õ��ع�����
        X_r31=[X_r32,zeros(size(X_r32))]*AB_31;
%         T_r31=linspace(0.5,length(X_r32),length(X_r31));

        A_30=a_1;B_30=b_1;
        for ki=1:(length(X_r31)-1)
            A_30=blkdiag(A_30,a_1);
            B_30=blkdiag(B_30,b_1);
        end
        AB_30=[A_30',B_30']';%�õ��ع�����
        X_r30=[X_r31,zeros(size(X_r31))]*AB_30;
        T_r30=linspace(0.5,length(X_r31),length(X_r30));
        %%%%����ֵ0.618max(X_r)
        gamma=c*max(X_r30);
        %%%%
        X_r30(X_r30<gamma)=0;
%         [pks30,locs30]=findpeaks(X_r30);%��ֵ
%         locs30=find(X_r30>0);%�ҵ�����0�ĵ������
%         T_r300=T_r30(locs30);%�ҵ���Ӧ��ʱ��
        T_r300=X_r30&T_r30;
    end
    
    if im==4%4�ηֽ�֮��ƽ���źŵ��ع�
        X_r43=X_1;
%         T_r43=t_wPf_1;
        %%%���ع�����
        A_42=a_1;B_42=b_1;
        for ki=1:(length(X_r43)-1)
            A_42=blkdiag(A_42,a_1);
            B_42=blkdiag(B_42,b_1);
        end
        AB_42=[A_42',B_42']';%�õ��ع�����
        X_r42=[X_r43,zeros(size(X_r43))]*AB_42;
%         T_r31=linspace(0.5,length(X_r32),length(X_r31));

        A_41=a_1;B_41=b_1;
        for ki=1:(length(X_r42)-1)
            A_41=blkdiag(A_41,a_1);
            B_41=blkdiag(B_41,b_1);
        end
        AB_41=[A_41',B_41']';%�õ��ع�����
        X_r41=[X_r42,zeros(size(X_r42))]*AB_41;
%         T_r41=linspace(0.5,length(X_r42),length(X_r41));

        A_40=a_1;B_40=b_1;
        for ki=1:(length(X_r41)-1)
            A_40=blkdiag(A_40,a_1);
            B_40=blkdiag(B_40,b_1);
        end
        AB_40=[A_40',B_40']';%�õ��ع�����
        X_r40=[X_r41,zeros(size(X_r41))]*AB_40;
        T_r40=linspace(0.5,length(X_r41),length(X_r40));
        %%%%����ֵ0.618max(X_r)
        gamma=c*max(X_r40);
        %%%%
        X_r40(X_r40<gamma)=0;
%         [pks40,locs40]=findpeaks(X_r40);%��ֵ
%         locs40=find(X_r40>0);%�ҵ�����0�ĵ������
%         T_r400=T_r40(locs40);%�ҵ���Ӧ��ʱ��
        T_r400=X_r40&T_r40;
    end
    
    if im==5
        X_r54=X_1;
%         T_r54=t_wPf_1;
        %%%���ع�����
        A_53=a_1;B_53=b_1;
        for ki=1:(length(X_r54)-1)
            A_53=blkdiag(A_53,a_1);
            B_53=blkdiag(B_53,b_1);
        end
        AB_53=[A_53',B_53']';%�õ��ع�����
        X_r53=[X_r54,zeros(size(X_r54))]*AB_53;
%         T_r31=linspace(0.5,length(X_r32),length(X_r31));

        A_52=a_1;B_52=b_1;
        for ki=1:(length(X_r53)-1)
            A_52=blkdiag(A_52,a_1);
            B_52=blkdiag(B_52,b_1);
        end
        AB_52=[A_52',B_52']';%�õ��ع�����
        X_r52=[X_r53,zeros(size(X_r53))]*AB_52;
%         T_r52=linspace(0.5,length(X_r53),length(X_r52));

        A_51=a_1;B_51=b_1;
        for ki=1:(length(X_r52)-1)
            A_51=blkdiag(A_51,a_1);
            B_51=blkdiag(B_51,b_1);
        end
        AB_51=[A_51',B_51']';%�õ��ع�����
        X_r51=[X_r52,zeros(size(X_r52))]*AB_51;
%         T_r51=linspace(0.5,length(X_r52),length(X_r51));

        A_50=a_1;B_50=b_1;
        for ki=1:(length(X_r51)-1)
            A_50=blkdiag(A_50,a_1);
            B_50=blkdiag(B_50,b_1);
        end
        AB_50=[A_50',B_50']';%�õ��ع�����
        X_r50=[X_r51,zeros(size(X_r51))]*AB_50;
        T_r50=linspace(0.5,length(X_r51),length(X_r50));
        %%%%����ֵ0.618max(X_r)
        gamma=c*max(X_r50);
        %%%%
        X_r50(X_r50<gamma)=0;
%         [pks50,locs50]=findpeaks(X_r50);%��ֵ
%         locs50=find(X_r50>0);%�ҵ�����0�ĵ������
%         T_r500=T_r50(locs50);%�ҵ���Ӧ��ʱ��
        T_r500=X_r50&T_r50;
    end
%     if im>=2
%         A1_1=a_1;B1_1=b_1;
%         for ki=1:(length(X_1)*1/2-1)
%             A1_1=blkdiag(A1_1,a_1);
%             B1_1=blkdiag(B1_1,b_1);
%         end
%         AB=[A1_1',B1_1']';%�õ��ֽ����
%         X_2=[X_1,zeros(size(X_1))]*AB;
%         t_wPf_2=linspace(0.5:230:460);
%     end
    % X_2=[zeros(size(X_s)),X_d]*AB;%�ֽ�һ��ϸ���ع�
    
    x=X_s;%��Ϊ�´ηֽ���ź�
    t_wPf=T;%��Ϊ�´ηֽ���źŶ�Ӧ��ʱ��
    L=length(x);%��ȡ�źŵĳ���
end
figure%�ع��ź�
set(gcf,'position',[0 0 600 400]);
subplot(511);
plot(T_r10,X_r10);
axis([0 500,-inf,inf]);
xlabel('(1)');
subplot(512);
plot(T_r20,X_r20);
axis([0 500,-inf,inf]);
xlabel('(2)');
subplot(513);
plot(T_r30,X_r30);
xlabel('(3)');
subplot(514);
plot(T_r40,X_r40);
xlabel('(4)');
subplot(515);
plot(T_r50,X_r50);
xlabel('(5)');
hold on;
% xlabel('δ�ӿ����ֵ');
%%
ML=length(T_r500);%�����ع��źŵ���󳤶�,�����Ĳ���
T_R=zeros(5,ML);
T_R(1,1:length(T_r100))=T_r100;
T_R(2,1:length(T_r200))=T_r200;
T_R(3,1:length(T_r300))=T_r300;
T_R(4,1:length(T_r400))=T_r400;
T_R(5,:)=T_r500;
T_R=T_R(:,(1:length(w)));
%%
%%%%%%%%%%%%%%%%%%%%%%%%��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�ڶ�����ֵ���ڶ��ηֽ��ع����ź������ÿ����ֵ��Ϊt_ridge=2*fix(8*0.382)-1����һ��ȥ�������ԡ�
t_ridge=fix(0.382*Lf);%������ֵ
T_2=T_R(2,:);%��ȡ���ηֽ��ع����ź���Ϣ
t_1loc=find(T_2);%��ȡT_2��1��λ������

t_1floc=find(diff([nan t_1loc nan])~=1);%1�������ֵ���ʼ����
t_1cL=diff(t_1floc)-1;%�����ĳ���-1
t_1floc=t_1floc(1:length(t_1floc)-1);%ȥ��β��

FP=find(t_1cL<t_ridge);%�ҵ�����С����ֵ�ļٷ��λ��
T_1floc=t_1floc(FP);%�ҵ������Գ�ʼ��λ��
T_1lloc=t_1floc(FP)+t_1cL(FP);%�ҵ��ü���������λ��

T_0=zeros(1,length(T_2));%����0����
for jj=1:length(T_1floc)%��Ӧ�����Ե�λ����-1
    T_0(t_1loc(T_1floc(jj)):t_1loc(T_1lloc(jj)))=-1;
end
T_R(2,:)=T_2+T_0;%ȥ��������
% % % % % % % % % % % % % % % 
x_r20=T_R(2,:).*X_r20(1:length(T_R(2,:)));
t_r20=T_r20(1:length(x_r20));
% % % % % % % % % % % % % % % 
% % % % % % % 
figure
subplot(2,1,1);
plot(T_r20,X_r20);
hold on;
subplot(2,1,2);
plot(t_r20,x_r20);
hold on;
% % % % % % % %
%%
% % % % % % % % % % % % % % % % % % % % % % % % ����©����%%%%%%%%%%%%%%%%%%%%%%
k_ridge=fix(0.618*5);%��������ֵ�����ü�����ֵΪfix(0.618*5)=3

T_L3=sum(T_R)-T_R(4,:)-T_R(5,:);%�õ�ǰ���еĺ�
T_r3=T_L3==k_ridge;%�����˵���3��ֵ����һ��4��5������ͬ��λ��ֻ������©�������(01)(10)(00)
%�ֱ�������е�����
%�ֱ���T_R(4,:)��T_R(5,:)���룬�������©������1�����п��ܳ���5��λ�ã�û����©ʱ���������ֲ���
t_4=T_r3-(T_r3&T_R(4,:))>0;t_5=T_r3-(T_r3&T_R(5,:))>0;%�ڵ�4�͵�5���ҵ���©�ķ�Ŀ���λ��
T_4=(t_4+T_R(4,:))>0;T_5=(t_5+T_R(5,:))>0;%©���λ�ò�1
%��T_4,T_5�ֱ��滻T_R(4,:)T_R(5,:)
T_R(4,:)=T_4;T_R(5,:)=T_5;
T_l5=sum(T_R);

T_R5=T_l5==5;
%%
X_R10=X_r10.*T_R5;%ȥ����������ź�
X_R10=X_R10(2:2:end);%ǰһ������λ���һ��ż��Ϊ��ͬ
T_R10=T_r10(2:2:end);

[P,LOC]=findpeaks(X_R10);%��ֵ������Ӧ��λ��

m=1;%�巶Χ����ʼλ��
ts=6;%ʱ����ֵ�����ã�����������������ʱ�䷶Χ���С��6������Ϊʵ������һ����
peak_num=0;%��¼�����
for ti=1:length(LOC)-1
    if LOC(ti+1)-LOC(ti)>ts         
        M_P=max(P(m:ti));%��Χ�����ֵ
        M_t=find(P(m:ti)==M_P)-1;%��Χ�ڷ����ֵ������
        M_T=LOC(m+M_t);%��Χ�ڷ����ֵ��Ӧ��ʱ��
        m=ti+1;
        peak_num=peak_num+1;
    end
    if peak_num>0
        P_t(1,peak_num)=M_P;%�ع������ֵ
        T_p(1,peak_num)=M_T;%������Ӧ��ʱ��
    end
    if LOC(ti+1)==max(LOC)
        M_P=max(P(m:ti+1));
        M_t=find(P(m:ti+1)==M_P)-1;%��Χ�ڷ����ֵ������
        M_T=LOC(m+M_t);%��Χ�ڷ����ֵ��Ӧ��ʱ��
%         m=ti+1;
        peak_num=peak_num+1;
    end
end
T_p(1,peak_num)=M_T;%������Ӧ��ʱ��
% P_t=x1(T_p);%�����ֵ
X_D1=X_d1(T_p);
x_D1=X_D1<0;
Xd=-0.5*x_D1;
%���Է����ֵ��Ӧ��ʱ������ֵ
Pt=T_p+Xd;

loc_Px=2*Pt;%���ֵ��Ӧ��λ������
Px=Xx(loc_Px);%���ֵ


figure(4)
set(gcf,'position',[0 200 600 400]);
axy1=plot(locs2,0,'bo');
hold on;
plot(fault_loc1,0,'gx');
axy2=plot(Pt,20,'r*');
% hold on;
% legend([axy1(1),axy2(1)],'original signal','noisy data');
for ii=1:length(Px);
    str=['(' num2str(Pt(ii)) ')'];
    text(Pt(ii)-12,20+10,str);
end
for ij=1:length(locs2);
    str=['(' num2str(locs2(ij)) ')'];
    text(locs2(ij)-12,20-40,str);
end
set(gca,'xtick',(0:100:max(T_wPf)));
set(gca,'ytick',(-100:50:100));
axis([0 max(T_wPf) 0 220]);
xlabel('����ֵ');
hold off;
% figure
% plot(T_r20,xx);
% hold on;
%%

%%%����ֵ�Ľ��

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_R(2,:)=T_2;
T_L3=sum(T_R)-T_R(4,:)-T_R(5,:);%�õ�ǰ���еĺ�
T_r3=T_L3==k_ridge;%�����˵���3��ֵ����һ��4��5������ͬ��λ��ֻ������©�������(01)(10)(00)
%�ֱ�������е�����
%�ֱ���T_R(4,:)��T_R(5,:)���룬�������©������1�����п��ܳ���5��λ�ã�û����©ʱ���������ֲ���
t_4=T_r3-(T_r3&T_R(4,:))>0;t_5=T_r3-(T_r3&T_R(5,:))>0;%�ڵ�4�͵�5���ҵ���©�ķ�Ŀ���λ��
T_4=(t_4+T_R(4,:))>0;T_5=(t_5+T_R(5,:))>0;%©���λ�ò�1
%��T_4,T_5�ֱ��滻T_R(4,:)T_R(5,:)
T_R(4,:)=T_4;T_R(5,:)=T_5;
T_l5=sum(T_R);

T_R5=T_l5==5;
%%
X_R10=X_r10.*T_R5;%ȥ����������ź�
X_R10=X_R10(2:2:end);%ǰһ������λ���һ��ż��Ϊ��ͬ
T_R10=T_r10(2:2:end);

[P,LOC]=findpeaks(X_R10);%��ֵ������Ӧ��λ��

m=1;%�巶Χ����ʼλ��
ts=6;%ʱ����ֵ�����ã�����������������ʱ�䷶Χ���С��6������Ϊʵ������һ����
peak_num=0;%��¼�����
for ti=1:length(LOC)-1
    if LOC(ti+1)-LOC(ti)>ts         
        M_P=max(P(m:ti));%��Χ�����ֵ
        M_t=find(P(m:ti)==M_P)-1;%��Χ�ڷ����ֵ������
        M_T=LOC(m+M_t);%��Χ�ڷ����ֵ��Ӧ��ʱ��
        m=ti+1;
        peak_num=peak_num+1;
    end
    if peak_num>0
        P_t(1,peak_num)=M_P;%�ع������ֵ
        T_p(1,peak_num)=M_T;%������Ӧ��ʱ��
    end
    if LOC(ti+1)==max(LOC)
        M_P=max(P(m:ti+1));
        M_t=find(P(m:ti+1)==M_P)-1;%��Χ�ڷ����ֵ������
        M_T=LOC(m+M_t);%��Χ�ڷ����ֵ��Ӧ��ʱ��
%         m=ti+1;
        peak_num=peak_num+1;
    end
end
T_p(1,peak_num)=M_T;%������Ӧ��ʱ��
% P_t=x1(T_p);%�����ֵ
X_D1=X_d1(T_p);
x_D1=X_D1<0;
Xd=-0.5*x_D1;
%���Է����ֵ��Ӧ��ʱ������ֵ
Pt=T_p+Xd;

loc_Px=2*Pt;%���ֵ��Ӧ��λ������
Px=Xx(loc_Px);%���ֵ


figure(5)
set(gcf,'position',[600 0 600 400]);
axy1=plot(locs2,0,'bo');
hold on;
plot(fault_loc1,0,'gx');
axy2=plot(Pt,20,'r*');
% hold on;
% legend([axy1(1),axy2(1)],'original signal','noisy data');
for ii=1:length(Px);
    str=['(' num2str(Pt(ii)) ')'];
    text(Pt(ii)-12,20+10,str);
end
for ij=1:length(locs2);
    str=['(' num2str(locs2(ij)) ')'];
    text(locs2(ij)-12,20-40,str);
end
set(gca,'xtick',(0:100:max(T_wPf)));
set(gca,'ytick',(-100:50:100));
axis([0 max(T_wPf) 0 220]);
xlabel('˫��ֵ');
hold off;
% figure
% plot(T_r20,xx);
% hold on;






