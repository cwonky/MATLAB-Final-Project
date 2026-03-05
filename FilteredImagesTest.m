% Filtering from Lecture 6 Matrices
% Read a color image
C = imread('peppers.png'); % replace with any image

% Channel names and colors
channelNames = {'Red','Green','Blue'};
channelColors = [1 0 0; 0 1 0; 0 0 1]; % RGB for visualization

figure;
for i = 1:3
  channelImg = zeros(size(C), 'uint8'); % same size as original
  channelImg(:,:,i) = C(:,:,i); % put values in the correct channel
  subplot(1,3,i);
  imshow(channelImg); % shows colored version
  title([channelNames{i} ' Channel']);
  xlabel('Columns');
  ylabel('Rows');
end

sgtitle('RGB Channels in Color');
