function varargout=boundsat(Y)
n=length(Y);
for j=1:n
if Y(j)>100
    Y(j)=100;
else
    if Y(j)<-100
        Y(j)=-100;
    end   
end
end

varargout{1}=Y;
end