function [dist,sharp]=newjudgedouble()
global double doublewidth linepixel sharpthreshold size2

profile=getprofile();

num=10;
% % Im=sxmopen("E:\STM\AU\Au_20200719_01_4.2K_006.sxm");
% % x=20;pic=scale(Im{1,1}.data)
% % profile=pic(x,:);
% % plot(1:512,profile)

[~,m]=size(profile);%微分
dif=zeros(num,m);
for line=1:num
    for i=1:m-1
    dif(line,i)=profile(line,i+1)-profile(line,i);
    end
end
% plot(1:512,dif)

dis=zeros(num,1);
sharpm=zeros(num,1);
sss=0;

%识别两个峰的距离
for i=1:num
    
    figure(5)
    plot(profile(i,:))
    figure(6)
    plot(dif(i,:))
    pause(1)
    
    ab=abs(dif(i,:));
    max1=max(ab);
    mid=find(ab==max1);
    max2=max(profile(i,:));
    min2=min(profile(i,:));
    [r,c]=size(mid);
    
     if (r*c==1)&&(mid>linepixel/5)&&(mid<4*linepixel/5)&&(max2-min2>1e-10)
        sss=sss+1;

        ju=0;
        for j=mid:linepixel-1
            if ((abs(dif(i,j+1))-abs(dif(i,j)))<0)&&(ju==0)
                ju=1;
            end
            if ((abs(dif(i,j+1))-abs(dif(i,j)))>0)&&(ju==1)
                ju=2;
            end
            if ((abs(dif(i,j+1))-abs(dif(i,j)))<0)&&(ju==2)
                if abs(dif(i,j))>max1/3
                    break
                else
                    ju=1;
                end
            end
        end
        jd=j-mid;

        ju=0;
        for jm=2:mid
            jj=mid+2-jm;
            if ((abs(dif(i,jj-1))-abs(dif(i,jj)))<0)&&(ju==0)
                ju=1;
            end
            if ((abs(dif(i,jj-1))-abs(dif(i,jj)))>0)&&(ju==1)
                ju=2;
            end
            if ((abs(dif(i,jj-1))-abs(dif(i,jj)))<0)&&(ju==2)
                if abs(dif(i,jj))>max1/3
                    break
                else
                    ju=1;
                end
            end
        end
        jjd=mid-jj;

        dis(sss,1)=min([jd,jjd]);


    %识别半宽             
            for j=mid:linepixel-1
                if abs(dif(i,j))<max1/2
                    break
                end
            end
            for jm=2:mid
                jj=mid+2-jm;
                 if abs(dif(i,jj))<max1/2
                    break
                 end
            end
            sharpm(sss,1)=j-jj;
     end 
end

if sss==0
    dist=0; 
    sharp=1;
    disp('here is not good,i will change scan place')
else
    dist=(size2{1,1}/linepixel)*sum(dis)/sss
    sharp=(size2{1,1}/linepixel)*sum(sharpm)/sss%转化成实际距离
end

if (dist<doublewidth)||(sharp>sharpthreshold)
   double=0;
else
   double=1;
end
end

