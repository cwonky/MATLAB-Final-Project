% Solicit image from user
image = imread(input("Enter image file name: ",'s')); 

% Find boundaries between different color indices (indexing) - we also need to figure out how many different colors we want)
boundaries = edge(indexed_img, 'sobel', 0); 
line_art = ~boundaries; % Invert so lines are black, background is white
