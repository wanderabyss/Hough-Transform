function img = drawline(img, r1,c1, r2,c2)
 
% Draws a white line segment in the image 
% characterized by two points
% (r1,c1)  and (r2,c2).

image = zeros(size(img));
[h,w,c] = size(img);
rr = [r1,r2];
cc = [c1,c2];

% calculate number of samples
deltar   = rr(2) - rr(1);
deltac   = cc(2) - cc(1);
line_len = sqrt(deltar.^2+deltac.^2);

for t=linspace(0,1,line_len*2)
    
    % find pixel to color
    r = round(rr(1) + t*deltar);
    c = round(cc(1) + t*deltac);
    
    % clip to draw within image boundaries
    if(r < 1 || r > size(img,1) || ...
       c < 1 || c > size(img,2) )
        continue;
    end
    
    % color the pixel white
    image(r,c) = uint8(255);
    
end

img = uint8(max(double(image), double(img)));