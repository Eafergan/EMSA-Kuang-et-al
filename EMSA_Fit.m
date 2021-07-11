clear all;

%This code generate random sets of data (using the Boot_S_Fit function) to fit, based on the data sets from observation 
%The function Boot_S_Fit creates a matrix of the 3 variables: K, f and C,
%the values of the matrix is the sum of the 3 equations that are fitted so
%the minimum value of the matrix correspond have the indices of the K, f, C
%the function need its variable range to be adjusted for each run.
%otherwise it will run out of memory, so use the histogram to check that
%the 3 variable are in range and if one of them is maxed out so adjust the
%range
%
SuConc= [2.50000000000000;2.50000000000000;2.50000000000000;2.50000000000000;10;10;10;10;40;40;40;40;160;160;160;160];

SuAloneCSL=[0.545249786228123,0.546166270401319,0.517765436522022,0.483750263445663,0.246609877518646,0.250198551190367,0.226187190409220,0.250834473705899,0.122343994543396,0.140065078243236,0.108941241966686,0.137221866472578,0.0710222680509683,0.0998918757321835,0.0580550401825572,0.0870017091186934;0.410153595596077,0.412101178547835,0.428636799638404,0.448181006217949,0.506780386494257,0.498526109161292,0.510717133817494,0.496886262040477,0.483555227995975,0.477689350111376,0.497649580485387,0.466139232474806,0.438047074549976,0.430860371775534,0.451628181616717,0.436303206110425;0.0445966181758008,0.0417325510508459,0.0535977638395739,0.0680687303363884,0.246609735987097,0.251275339648341,0.263095675773287,0.252279264253624,0.394100777460630,0.382245571645389,0.393409177547926,0.396638901052616,0.490930657399056,0.469247752492283,0.490316778200726,0.476695084770882];
SuHCSL=[0.558881427637718,0.562661395910912,0.573161817486308,0.500182581343826,0.274849227839397,0.300228811140753,0.283927349532005,0.266533619348185,0.157623249519112,0.189697106066217,0.170067266463168,0.145671141623194,0.0751184595578876,0.104738222133878,0.0914918367382075,0.0769780664043888;0.421953673262292,0.414151992964862,0.404781744489432,0.460150769159152,0.576054611494268,0.546472236768704,0.555441594939554,0.567965818549282,0.549419345079403,0.537816048189416,0.534696374420060,0.544698525208513,0.467848653012233,0.455626189148864,0.470764773444766,0.461744684659595;0.0191648990999895,0.0231866111242262,0.0220564380242601,0.0396666494970222,0.149096160666335,0.153298952090543,0.160631055528441,0.165500562102533,0.292957405401486,0.272486845744367,0.295236359116773,0.309630333168293,0.457032887429880,0.439635588717258,0.437743389817026,0.461277248936016];
NCMCSL=[0.521494790976041,0.517759274467774,0.520584017310076,0.492424787697724,0.304712168332083,0.350377451454960,0.335722629079946,0.289902535300048,0.185522483692880,0.225650678857885,0.211842968693699,0.183435194357907,0.151012873073741,0.183563052613792,0.172758364883362,0.141263178589461;0.371863564387684,0.380614244626900,0.367807127164924,0.376566677907729,0.562604340095156,0.523511278982543,0.530277531419885,0.556573565585562,0.593732791577540,0.560665380085759,0.559303556613926,0.558424408827786,0.525991741227660,0.532211852317592,0.522283450955586,0.534143663464959;0.106641644636274,0.101626480905325,0.111608855525000,0.131008534394547,0.132683491572761,0.126111269562497,0.133999839500169,0.153523899114391,0.220744724729580,0.213683941056356,0.228853474692375,0.258140396814308,0.322995385698598,0.284225095068616,0.304958184161052,0.324593157945581];
SuAloneSPS=[0.771586446455311,0.793350337181586,0.785880054263850,0.716583455976604,0.465585568892672,0.472373835765532,0.475513077542496,0.452706492144550,0.236204506388552,0.269091882119542,0.256595569518366,0.234949548929454,0.0577296581813838,0.0845824166214856,0.0686619400476650,0.0696943220158690;0.216827191685224,0.199476664609649,0.216857658250378,0.268470687347826,0.468262117200589,0.453103164561889,0.450924245226194,0.462811572418345,0.542136700644378,0.525099230971803,0.526858071228679,0.521403285603923,0.492003904474484,0.483654307267884,0.489713146656066,0.489896378477929;0.0115863618594658,0.00717299820876441,-0.00273771251422826,0.0149458566755693,0.0661523139067395,0.0745229996725793,0.0735626772313102,0.0844819354371049,0.221658792967069,0.205808886908655,0.216546359252955,0.243647165466623,0.450266437344132,0.431763276110630,0.441624913296269,0.440409299506203];
SuHSPS=[0.846628734340396,0.810096738627880,0.825256019022269,0.747481378988168,0.539340673352500,0.554154343779605,0.539679849488197,0.511197271979013,0.328026458706950,0.369961248932509,0.336065675210485,0.319683058094860,0.107748793696603,0.121283353553596,0.108379099737312,0.114449547828589;0.138412429701340,0.168734277375328,0.161619531880675,0.215264857179509,0.407079523937102,0.388084272568737,0.406824744089898,0.410466776938580,0.524825212561232,0.504411910576795,0.520765167641638,0.512776144597272,0.504376755177958,0.493184895548807,0.509511419442579,0.483844579206034;0.0149588359582639,0.0211689839967918,0.0131244490970560,0.0372537638323234,0.0535798027103987,0.0577613836516583,0.0534954064219053,0.0783359510824070,0.147148328731818,0.125626840490696,0.143169157147877,0.167540797307869,0.387874451125439,0.385531750897597,0.382109480820109,0.401705872965377];
NCMSPS=[0.503653411402363,0.468001348962774,0.491545346536033,0.486193075008541,0.204155228229542,0.224334003167623,0.244452136883567,0.216956339954379,0.0848555579363370,0.102759443104920,0.110212812792381,0.102634321344148,0.0499515852463960,0.0615142535358324,0.0578425293750376,0.0595076285109280;0.170637220518472,0.189909879026973,0.192965165800442,0.188272946374700,0.292645533253886,0.297728938239934,0.289364448451435,0.301136054581915,0.330062335137255,0.330768816095351,0.330772665972390,0.332225972209887,0.323333189142093,0.337048971596252,0.350747960440846,0.348726885609464;0.325709368079164,0.342088772010253,0.315489487663525,0.325533978616759,0.503199238516572,0.477937058592443,0.466183414664998,0.481907605463706,0.585082106926408,0.566471740799729,0.559014521235230,0.565139706445965,0.626715225611511,0.601436774867916,0.591409510184116,0.591765485879608];


runs=1000;

TheList=zeros(3,runs);

SetToFit2=SuAloneCSL;

RanGenMean=zeros(3,16);
RanGenSD=zeros(3,16);

for i=1:3
RanGenMean(i,1:4)=mean(SetToFit2(i,1:4));
RanGenMean(i,5:8)=mean(SetToFit2(i,5:8));
RanGenMean(i,9:12)=mean(SetToFit2(i,9:12));
RanGenMean(i,13:16)=mean(SetToFit2(i,13:16));
RanGenSD(i,1:4)=std(SetToFit2(i,1:4));
RanGenSD(i,5:8)=std(SetToFit2(i,5:8));
RanGenSD(i,9:12)=std(SetToFit2(i,9:12));
RanGenSD(i,13:16)=std(SetToFit2(i,13:16));
end

Klist=zeros(1,runs);
C1list=zeros(1,runs);
Flist=zeros(1,runs);

parfor i=1:runs
RanMat = RanGenSD.*randn(3,16) + RanGenMean;
[Kdi,C1i,Fi] = Boot_S_Fit(RanMat);

Klist(i)=Kdi;
C1list(i)=C1i;
Flist(i)=Fi;
end


TheList(1,:)=Klist;
TheList(2,:)=C1list;
TheList(3,:)=Flist;

figure
hist(Klist);
title('K');

figure
hist(Flist)
title('F');

figure
hist(C1list);
title('C')


p50k=prctile(TheList(1,:),50);
p50C=prctile(TheList(2,:),50);
p50F=prctile(TheList(3,:),50);

plusminusC=0.5*(prctile(TheList(2,:),95)-prctile(TheList(2,:),5));
plusminusK=0.5*(prctile(TheList(1,:),95)-prctile(TheList(1,:),5));
plusminusF=0.5*(prctile(TheList(3,:),95)-prctile(TheList(3,:),5));


function [Kd,C1,F]=Boot_S_Fit(SetToFit)
 
%1-f squere
Alpha=@(Kd,x)  (x*(1./Kd));
fnc0=@(c1,Kd,f,x)(  (1-2.*f-f.^2).*(1./(1 + 2.*Alpha(Kd,x) + c1.*(Alpha(Kd,x).^2))) + (2.*f).*(1./(1 + Alpha(Kd,x))) + f.^2) ; 
fnc1=@(c1,Kd,f,x)(  (1-2.*f-f.^2).*((2.*Alpha(Kd,x))./(1 + 2.*Alpha(Kd,x) + c1.*(Alpha(Kd,x).^2))) + (2.*f).*((Alpha(Kd,x))./(1 + Alpha(Kd,x))) ) ; 
fnc2=@(c1,Kd,f,x)( (1-2.*f-f.^2).*(((c1.*Alpha(Kd,x).^2))./(1 + 2.*Alpha(Kd,x) + c1.*(Alpha(Kd,x).^2)))   ); 

 
 
SuConc=[2.5 2.5 2.5 2.5 10 10 10 10 40 40 40 40 160 160 160 160]';

%f for sps 0.77 for CSL 0.6

C1m=single(0.4:0.1:2);
fm=single(0.18:0.01:0.35);
Kdm=single(3:0.1:10);

C1m2=single(ones(1,1,length(C1m)));
fm2=single(ones(1,1,1,length(fm)));

for c=1:length(C1m)
    C1m2(1,1,c)=C1m(c);
end

for f=1:length(fm)
    fm2(1,1,1,f)=fm(f);
end

SetToFitM0=Alpha(ones(1,length(Kdm)),SetToFit(1,:)');
SetToFitM0=SetToFitM0.*ones(1,1,length(C1m));
SetToFitM0=SetToFitM0.*ones(1,1,1,length(fm));
SetToFitM1=Alpha(ones(1,length(Kdm)),SetToFit(2,:)');
SetToFitM1=SetToFitM1.*ones(1,1,length(C1m));
SetToFitM1=SetToFitM1.*ones(1,1,1,length(fm));
SetToFitM2=Alpha(ones(1,length(Kdm)),SetToFit(3,:)');
SetToFitM2=SetToFitM2.*ones(1,1,length(C1m));
SetToFitM2=SetToFitM2.*ones(1,1,1,length(fm));

leastSqM=squeeze(sum((SetToFitM0-fnc0(C1m2,Kdm,fm2,SuConc)).^2)+sum((SetToFitM1-fnc1(C1m2,Kdm,fm2,SuConc)).^2)+sum((SetToFitM2-fnc2(C1m2,Kdm,fm2,SuConc)).^2));


[mxv,idx] = max(1./leastSqM(:));
[Kdind,Cind,Find] = ind2sub(size(leastSqM),idx);
C1=C1m(Cind);
Kd=Kdm(Kdind);
F=fm(Find);
clear leastSqM;
end