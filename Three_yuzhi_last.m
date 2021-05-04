clear all; close all; clc
%%
                                                    %%%%%%%%%%%%%产生仿真数据%%%%%%%%%%%%
N=10;          %峰值信号的数目
fs=2000;        %采样率
f=60;           %频率

A=100;          %振幅

sparsity=0.8;   %稀疏性
SNR=5;         %信噪比，典型信噪比取值1-100；此处用来计算sigma(噪声方差)
buffer=16;    %缓冲区确保两个瞬时信号不会彼此重叠
%%
                                                            %模拟普通峰值信号
Lf=fs/(2*f)+1;
n=0:Lf-1;n=[n,Lf-1];%n的取值0,1,2...Lf-1  
Pf=A*sin(pi*n/(Lf-1));%sin(wx)所以w=(2*pi/(Lf-1)),x=n/2   周期为16.666666ms　　半周期为8.333   所以对应时间为n/2
%%
                                                    %求高斯白噪声标准差，方差    得到高斯白噪声
M=Lf;             %单个峰信号的长度，半周期的正玄波长度
s1=0;
s2=Lf-1; 
s=Pf;      %%%取为峰信号%%%
S=sum(s.^2);%求行向量的元素平方和
sigma=sqrt(S/(SNR*(s2-s1+1)));%计算sigma，用于产生白噪声
L=fix((M/(1-sparsity)*N)+Lf-1);            %数据x的整个长度
w=sqrt(sigma.^2)*randn(1,L);%产生方差为sigma平方长度为L的高斯白噪声


%%%%%%%%%%%%高斯白噪声的最大值用以计算所需故障数据的大小
max_w=max(abs(w));%%正常的峰值为100所以我们取故障数据为
fault=max_w+0.618*A;
%%%%%%%%%%%%

    w1=w(2:2:end);%对应时间为整数的信号取值

t_wPf=0.5:0.5:L/2;%噪声的时间范围                           对应时间为n/2
    T_wPf=t_wPf;
    t_wPf1=1:L/2;%对应的整数时间
%%
                                            %在信号中选择峰的位置，buffer的使用是防止相邻的两个峰出现叠加
i=1:N;
a=1+(i-1)*M/(1-sparsity)+buffer;b=i*M/(1-sparsity);%%%%%%在[a,b]区间随机选取峰信号的位置
%(b-a)*rand(1)+a
tau_i=(b-a).*rand(1,N)+a;    %在区间随机选择得到每个峰信号的坐标tau_i;
tau_i1=repmat(tau_i',1,length(n));%得到列相同的矩阵
tau_i=round(tau_i);%对位置进行四舍五入取整方便作为下标
%%
                                                            %%%%%%%%%%%%合成信号 
l=zeros(N,length(n));%已知l矩阵的大小，所以从内存中直接取得大小，防止在loop中每次对内存的取出
k=0:Lf-1;
k=[k,Lf];
for j=1:N
    k1=repmat(k,N,1);
    l(j,:)=k1(j,:)+tau_i1(j,:);%总共10个峰，每个峰由Lf=18个坐标表示，每个峰的位置由第一个位置确定,得到10*18 的坐标矩阵
end
                                                        %%%%%%%%%%%%得到有用的连续信号
PM=zeros(N,L);%已知PM的大小
for d=1:N;
    Fz=zeros(1,tau_i(d)-2);%由于本身从0开始并且要比峰位置小1，所以此处减2，峰位置前全部补0
    Bz=zeros(1,L-tau_i(d)+2-length(Pf));%峰之后位置全补0
    P=[Fz,Pf,Bz];
    PM(d,:)=P;     %得到峰的矩阵表示
end
%%
                                                                           %%%%%%%%%对得到的峰信号以及噪声信号进行叠加
sum_P=sum(PM,1);%求列和得到一个行向量
sum_P1=sum_P(2:2:end);
[pks1,locs1]=findpeaks(sum_P);%每个尺度上寻找所有极大值的坐标
%%%%%%%%%%%%%%%%%%%找到非有用信号的位置，用以随机选取故障信号的位置
not_s=find(sum_P==0);
ran_int=randperm(length(not_s),1);%在1:length(not_s)范围内产生一个随机整数
fault_loc=not_s(ran_int);%求得故障信号的位置
fault_loc1=1/2*fault_loc;
%%%%%%%%%%%%%%%%%%%%
locs2=t_wPf(locs1);%原始数据中峰对应的时间

Pin=zeros(1,length(sum_P));%添加脉冲干扰信号pin
Pin(fault_loc)=fault;
x=sum_P+Pin+w;%%%%%%%%%%%%%%%%%%%%%%%%时间为n/2
Xx=x;
ES=Xx(locs1);%原始数据中信号峰对应的有效信号的峰值
x1=x(2:2:end);%对应时间为整数的信号取值
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
%%%%%%%%%小波变换
Pj=zeros(5,length(x));
time=zeros(5,length(x));
len=zeros(1,length(x));
a_1=[sqrt(1/2),sqrt(1/2)];b_1=[-sqrt(1/2),sqrt(1/2)];
c=0.382;
for im=1:5
    if mod(L,2)==1
        x=[x,0];%信号长度为奇数时在末尾补零成为偶数位
        t_wPf=[t_wPf,max(t_wPf)+0.5*im];
    end
    X=x;
    T=t_wPf(2:2:end);
%%
%分解矩阵
    A_1=a_1;B_1=b_1;
    for ki=1:(length(x)*1/2-1)
        A_1=blkdiag(A_1,a_1);
        B_1=blkdiag(B_1,b_1);
    end
    AB=[A_1',B_1']';%得到分解矩阵
%%
    %得到平滑和细节信号
    X_sd=(AB*X')';

    X_s=X_sd(1:length(X_sd)*1/2);%平滑信号
    X_d=X_sd(length(X_sd)*1/2+1:length(X_sd));%细节信号
    if (im==1)
        X_s1=X_s;%分解一次的平滑信号
        X_d1=X_d;%分解一次的细节信号
        T1=T;%分解一次平滑信号对应的时间
    end
    if (im==2)
        X_s2=X_s;%分解2次的平滑信号
        X_d2=X_d;%分解2次的细节信号
        T2=T;%分解2次平滑信号对应的时间   
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
%对信号进行重构
    X_1=[X_s,zeros(size(X_s))]*AB;%平滑信号重构一次
    t_wPf_1=t_wPf;%对应的时间
    if im==1
        X_r10=X_1;
        T_r10=t_wPf_1;
        %%%%求阈值0.618max(X_r)%%%%%%%%%%%第一个阈值%%%%%
        gamma=c*max(X_r10);
        %%%%
        X_r10(X_r10<gamma)=0;
%         [pks10,locs10]=findpeaks(X_r10);%求极值
%         locs10=find(X_r10>0);%找到大于0的点的坐标
%         T_r100=T_r10(locs10);%找到对应的时间
        T_r100=X_r10&T_r10;
    end
    
    if im==2%两次分解之后平滑信号的重构
        X_r21=X_1;
%         T_r21=t_wPf_1;
        %%%求重构矩阵
        A_20=a_1;B_20=b_1;
        for ki=1:(length(X_r21)-1)
            A_20=blkdiag(A_20,a_1);
            B_20=blkdiag(B_20,b_1);
        end
        AB_20=[A_20',B_20']';%得到重构矩阵
        X_r20=[X_r21,zeros(size(X_r21))]*AB_20;
        T_r20=linspace(0.5,length(X_r21),length(X_r20));
        %%%%求阈值0.618max(X_r)
        gamma=c*max(X_r20);
        %%%%%%%%%%%%%%%%%%%%%%
        xx=X_r20;
        xx(xx<0)=0;
        %%%%%%%%%%%%%%%%%%%%%%
        X_r20(X_r20<gamma)=0;
%         [pks20,locs20]=findpeaks(X_r20);%求极值
%         locs20=find(X_r20>0);%找到大于0的点的坐标
%         T_r200=T_r20(locs20);%找到对应的时间
        T_r200=X_r20&T_r20;
    end
    
    if im==3%三次分解之后平滑信号的重构
        X_r32=X_1;
%         T_r32=t_wPf_1;
        %%%求重构矩阵
        A_31=a_1;B_31=b_1;
        for ki=1:(length(X_r32)-1)
            A_31=blkdiag(A_31,a_1);
            B_31=blkdiag(B_31,b_1);
        end
        AB_31=[A_31',B_31']';%得到重构矩阵
        X_r31=[X_r32,zeros(size(X_r32))]*AB_31;
%         T_r31=linspace(0.5,length(X_r32),length(X_r31));

        A_30=a_1;B_30=b_1;
        for ki=1:(length(X_r31)-1)
            A_30=blkdiag(A_30,a_1);
            B_30=blkdiag(B_30,b_1);
        end
        AB_30=[A_30',B_30']';%得到重构矩阵
        X_r30=[X_r31,zeros(size(X_r31))]*AB_30;
        T_r30=linspace(0.5,length(X_r31),length(X_r30));
        %%%%求阈值0.618max(X_r)
        gamma=c*max(X_r30);
        %%%%
        X_r30(X_r30<gamma)=0;
%         [pks30,locs30]=findpeaks(X_r30);%求极值
%         locs30=find(X_r30>0);%找到大于0的点的坐标
%         T_r300=T_r30(locs30);%找到对应的时间
        T_r300=X_r30&T_r30;
    end
    
    if im==4%4次分解之后平滑信号的重构
        X_r43=X_1;
%         T_r43=t_wPf_1;
        %%%求重构矩阵
        A_42=a_1;B_42=b_1;
        for ki=1:(length(X_r43)-1)
            A_42=blkdiag(A_42,a_1);
            B_42=blkdiag(B_42,b_1);
        end
        AB_42=[A_42',B_42']';%得到重构矩阵
        X_r42=[X_r43,zeros(size(X_r43))]*AB_42;
%         T_r31=linspace(0.5,length(X_r32),length(X_r31));

        A_41=a_1;B_41=b_1;
        for ki=1:(length(X_r42)-1)
            A_41=blkdiag(A_41,a_1);
            B_41=blkdiag(B_41,b_1);
        end
        AB_41=[A_41',B_41']';%得到重构矩阵
        X_r41=[X_r42,zeros(size(X_r42))]*AB_41;
%         T_r41=linspace(0.5,length(X_r42),length(X_r41));

        A_40=a_1;B_40=b_1;
        for ki=1:(length(X_r41)-1)
            A_40=blkdiag(A_40,a_1);
            B_40=blkdiag(B_40,b_1);
        end
        AB_40=[A_40',B_40']';%得到重构矩阵
        X_r40=[X_r41,zeros(size(X_r41))]*AB_40;
        T_r40=linspace(0.5,length(X_r41),length(X_r40));
        %%%%求阈值0.618max(X_r)
        gamma=c*max(X_r40);
        %%%%
        X_r40(X_r40<gamma)=0;
%         [pks40,locs40]=findpeaks(X_r40);%求极值
%         locs40=find(X_r40>0);%找到大于0的点的坐标
%         T_r400=T_r40(locs40);%找到对应的时间
        T_r400=X_r40&T_r40;
    end
    
    if im==5
        X_r54=X_1;
%         T_r54=t_wPf_1;
        %%%求重构矩阵
        A_53=a_1;B_53=b_1;
        for ki=1:(length(X_r54)-1)
            A_53=blkdiag(A_53,a_1);
            B_53=blkdiag(B_53,b_1);
        end
        AB_53=[A_53',B_53']';%得到重构矩阵
        X_r53=[X_r54,zeros(size(X_r54))]*AB_53;
%         T_r31=linspace(0.5,length(X_r32),length(X_r31));

        A_52=a_1;B_52=b_1;
        for ki=1:(length(X_r53)-1)
            A_52=blkdiag(A_52,a_1);
            B_52=blkdiag(B_52,b_1);
        end
        AB_52=[A_52',B_52']';%得到重构矩阵
        X_r52=[X_r53,zeros(size(X_r53))]*AB_52;
%         T_r52=linspace(0.5,length(X_r53),length(X_r52));

        A_51=a_1;B_51=b_1;
        for ki=1:(length(X_r52)-1)
            A_51=blkdiag(A_51,a_1);
            B_51=blkdiag(B_51,b_1);
        end
        AB_51=[A_51',B_51']';%得到重构矩阵
        X_r51=[X_r52,zeros(size(X_r52))]*AB_51;
%         T_r51=linspace(0.5,length(X_r52),length(X_r51));

        A_50=a_1;B_50=b_1;
        for ki=1:(length(X_r51)-1)
            A_50=blkdiag(A_50,a_1);
            B_50=blkdiag(B_50,b_1);
        end
        AB_50=[A_50',B_50']';%得到重构矩阵
        X_r50=[X_r51,zeros(size(X_r51))]*AB_50;
        T_r50=linspace(0.5,length(X_r51),length(X_r50));
        %%%%求阈值0.618max(X_r)
        gamma=c*max(X_r50);
        %%%%
        X_r50(X_r50<gamma)=0;
%         [pks50,locs50]=findpeaks(X_r50);%求极值
%         locs50=find(X_r50>0);%找到大于0的点的坐标
%         T_r500=T_r50(locs50);%找到对应的时间
        T_r500=X_r50&T_r50;
    end
%     if im>=2
%         A1_1=a_1;B1_1=b_1;
%         for ki=1:(length(X_1)*1/2-1)
%             A1_1=blkdiag(A1_1,a_1);
%             B1_1=blkdiag(B1_1,b_1);
%         end
%         AB=[A1_1',B1_1']';%得到分解矩阵
%         X_2=[X_1,zeros(size(X_1))]*AB;
%         t_wPf_2=linspace(0.5:230:460);
%     end
    % X_2=[zeros(size(X_s)),X_d]*AB;%分解一次细节重构
    
    x=X_s;%作为下次分解的信号
    t_wPf=T;%作为下次分解的信号对应的时间
    L=length(x);%读取信号的长度
end
figure%重构信号
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
% xlabel('未加宽度阈值');
%%
ML=length(T_r500);%计算重构信号的最大长度,不够的补零
T_R=zeros(5,ML);
T_R(1,1:length(T_r100))=T_r100;
T_R(2,1:length(T_r200))=T_r200;
T_R(3,1:length(T_r300))=T_r300;
T_R(4,1:length(T_r400))=T_r400;
T_R(5,:)=T_r500;
T_R=T_R(:,(1:length(w)));
%%
%%%%%%%%%%%%%%%%%%%%%%%%降低误报率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%第二个阈值，在二次分解重构的信号中设置宽度阈值记为t_ridge=2*fix(8*0.382)-1，进一步去掉假阳性。
t_ridge=fix(0.382*Lf);%设置阈值
T_2=T_R(2,:);%读取二次分解重构的信号信息
t_1loc=find(T_2);%读取T_2中1的位置坐标

t_1floc=find(diff([nan t_1loc nan])~=1);%1连续出现的起始坐标
t_1cL=diff(t_1floc)-1;%连续的长度-1
t_1floc=t_1floc(1:length(t_1floc)-1);%去掉尾巴

FP=find(t_1cL<t_ridge);%找到长度小于阈值的假峰的位置
T_1floc=t_1floc(FP);%找到假阳性初始的位置
T_1lloc=t_1floc(FP)+t_1cL(FP);%找到该假阳性最后的位置

T_0=zeros(1,length(T_2));%生成0向量
for jj=1:length(T_1floc)%对应假阳性的位置置-1
    T_0(t_1loc(T_1floc(jj)):t_1loc(T_1lloc(jj)))=-1;
end
T_R(2,:)=T_2+T_0;%去掉假阳性
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
% % % % % % % % % % % % % % % % % % % % % % % % 降低漏报率%%%%%%%%%%%%%%%%%%%%%%
k_ridge=fix(0.618*5);%第三个阈值，设置脊长阈值为fix(0.618*5)=3

T_L3=sum(T_R)-T_R(4,:)-T_R(5,:);%得到前三行的和
T_r3=T_L3==k_ridge;%保留了等于3的值并置一（4、5行在相同峰位置只出现遗漏的情况）(01)(10)(00)
%分别检查第四行第五行
%分别与T_R(4,:)和T_R(5,:)相与，如果有遗漏，则补上1，最有可能出现5的位置，没有遗漏时，基本保持不变
t_4=T_r3-(T_r3&T_R(4,:))>0;t_5=T_r3-(T_r3&T_R(5,:))>0;%在第4和第5行找到遗漏的峰的可能位置
T_4=(t_4+T_R(4,:))>0;T_5=(t_5+T_R(5,:))>0;%漏检的位置补1
%用T_4,T_5分别替换T_R(4,:)T_R(5,:)
T_R(4,:)=T_4;T_R(5,:)=T_5;
T_l5=sum(T_R);

T_R5=T_l5==5;
%%
X_R10=X_r10.*T_R5;%去掉噪声后峰信号
X_R10=X_R10(2:2:end);%前一个奇数位与后一个偶数为相同
T_R10=T_r10(2:2:end);

[P,LOC]=findpeaks(X_R10);%求极值及所对应的位置

m=1;%峰范围的起始位置
ts=6;%时间阈值的设置，当相邻两个峰所在时间范围相距小于6我们认为实际属于一个峰
peak_num=0;%记录峰个数
for ti=1:length(LOC)-1
    if LOC(ti+1)-LOC(ti)>ts         
        M_P=max(P(m:ti));%范围内最大值
        M_t=find(P(m:ti)==M_P)-1;%范围内峰最大值的坐标
        M_T=LOC(m+M_t);%范围内峰最大值对应的时间
        m=ti+1;
        peak_num=peak_num+1;
    end
    if peak_num>0
        P_t(1,peak_num)=M_P;%重构峰最大值
        T_p(1,peak_num)=M_T;%峰最大对应的时间
    end
    if LOC(ti+1)==max(LOC)
        M_P=max(P(m:ti+1));
        M_t=find(P(m:ti+1)==M_P)-1;%范围内峰最大值的坐标
        M_T=LOC(m+M_t);%范围内峰最大值对应的时间
%         m=ti+1;
        peak_num=peak_num+1;
    end
end
T_p(1,peak_num)=M_T;%峰最大对应的时间
% P_t=x1(T_p);%峰最大值
X_D1=X_d1(T_p);
x_D1=X_D1<0;
Xd=-0.5*x_D1;
%所以峰最大值对应的时间和最大值
Pt=T_p+Xd;

loc_Px=2*Pt;%最大值对应的位置坐标
Px=Xx(loc_Px);%最大值


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
xlabel('多阈值');
hold off;
% figure
% plot(T_r20,xx);
% hold on;
%%

%%%二阈值的结果

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_R(2,:)=T_2;
T_L3=sum(T_R)-T_R(4,:)-T_R(5,:);%得到前三行的和
T_r3=T_L3==k_ridge;%保留了等于3的值并置一（4、5行在相同峰位置只出现遗漏的情况）(01)(10)(00)
%分别检查第四行第五行
%分别与T_R(4,:)和T_R(5,:)相与，如果有遗漏，则补上1，最有可能出现5的位置，没有遗漏时，基本保持不变
t_4=T_r3-(T_r3&T_R(4,:))>0;t_5=T_r3-(T_r3&T_R(5,:))>0;%在第4和第5行找到遗漏的峰的可能位置
T_4=(t_4+T_R(4,:))>0;T_5=(t_5+T_R(5,:))>0;%漏检的位置补1
%用T_4,T_5分别替换T_R(4,:)T_R(5,:)
T_R(4,:)=T_4;T_R(5,:)=T_5;
T_l5=sum(T_R);

T_R5=T_l5==5;
%%
X_R10=X_r10.*T_R5;%去掉噪声后峰信号
X_R10=X_R10(2:2:end);%前一个奇数位与后一个偶数为相同
T_R10=T_r10(2:2:end);

[P,LOC]=findpeaks(X_R10);%求极值及所对应的位置

m=1;%峰范围的起始位置
ts=6;%时间阈值的设置，当相邻两个峰所在时间范围相距小于6我们认为实际属于一个峰
peak_num=0;%记录峰个数
for ti=1:length(LOC)-1
    if LOC(ti+1)-LOC(ti)>ts         
        M_P=max(P(m:ti));%范围内最大值
        M_t=find(P(m:ti)==M_P)-1;%范围内峰最大值的坐标
        M_T=LOC(m+M_t);%范围内峰最大值对应的时间
        m=ti+1;
        peak_num=peak_num+1;
    end
    if peak_num>0
        P_t(1,peak_num)=M_P;%重构峰最大值
        T_p(1,peak_num)=M_T;%峰最大对应的时间
    end
    if LOC(ti+1)==max(LOC)
        M_P=max(P(m:ti+1));
        M_t=find(P(m:ti+1)==M_P)-1;%范围内峰最大值的坐标
        M_T=LOC(m+M_t);%范围内峰最大值对应的时间
%         m=ti+1;
        peak_num=peak_num+1;
    end
end
T_p(1,peak_num)=M_T;%峰最大对应的时间
% P_t=x1(T_p);%峰最大值
X_D1=X_d1(T_p);
x_D1=X_D1<0;
Xd=-0.5*x_D1;
%所以峰最大值对应的时间和最大值
Pt=T_p+Xd;

loc_Px=2*Pt;%最大值对应的位置坐标
Px=Xx(loc_Px);%最大值


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
xlabel('双阈值');
hold off;
% figure
% plot(T_r20,xx);
% hold on;






