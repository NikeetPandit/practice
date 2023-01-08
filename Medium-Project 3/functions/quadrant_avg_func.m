%--- Extracting quadrants of image (lower right, lower left, upper right, upper left)
%--- These quadrants are where there is only the background of the image
%... then averaging to take noise estimate where there is only background
%... taking average is hoped to have a more realistic noise estimate
%
%Written by Nikeet Pandit

function avg_image_noise = quadrant_avg_func(ImageIn)
cnt = 0; 
quadrants = [1:101; 1:101; 1:101; 400:500; 400:500; 1:101; 400:500; 400:500]; %quadrants where there is image background from rose
for i = 1:2:size(quadrants,1)
    cnt = cnt + 1;
    avg_image_noise(:,:,cnt) = ImageIn(quadrants(i,:),quadrants(i+1,:));

    if i+1 == size(quadrants,1)
        break
    end
end
avg_image_noise = mean(avg_image_noise,3); %taking average along 3rd matrix dimension
