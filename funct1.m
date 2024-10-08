%read vids
%this function reads in all the images into a vector
count = 0;

shell = zeros(v.Height,v.Width,30);

for x =1:30
    %this gave an error bacause i'm saving one image in a vector for size
    %30
    shell(x) = double(read(v,x));
end 
