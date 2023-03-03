
% Controls the amount of compression
compression_factor = 0.98;

% Read the source image and convert to grayscale
B1 = rescale(double(imread("photos\sea_1024.png")));
B1 = (B1(:,:,1) + B1(:,:,2) + B1(:,:,3))/3;

figure(1);
imshow(B1);

N = 8; % Block size in one dimension
nvec = 0:(N-1);
mvec = 0:(N-1);
[nmat, mmat] = meshgrid(nvec, mvec);

% We store the 8x8 blocks of the source image here
source_blocks = cell((size(B1,1)/8).^2,1);

% We store the compressed 8x8 blocks here
result_blocks = source_blocks;

% We store the DCT-II building blocks here 
building_blocks = cell(64,1);

% Split the source image into 8x8 blocks
bs = 8;
bc = size(B1, 1)/8; % block count
for i = 1:bc
    for j = 1:bc
        lin_idx = (i - 1) * bc + j;
        start_i = 1 + (i - 1) * bs; % i * bs
        start_j = 1 + (j - 1) * bs;
        end_i = start_i + bs - 1;
        end_j = start_j + bs - 1;
        source_blocks{lin_idx} = B1(start_i : end_i, start_j : end_j);
    end
end

for i = 1:(size(B1,1)/8)
end

idx = 1;

figure(2);

% Computing the DCT-II building blocks and
% storing these in the building_blocks cell array
for i = 0:(N-1)
	for j = 0:(N-1)
        % Sampling the cosine
		c1 = cos(pi/N*(nmat+1/2) * i);
		c2 = cos(pi/N*(mmat+1/2) * j);
		block = (c1 .* c2);

        % Normalize block
        block_length = sum(sum(block .* block)).^(1/2);
        block = block ./ block_length;

        % Plotting each building block
		subplot(8, 8, j * N + i + 1);
		imagesc(block);
		colormap("gray");

        % Store building block for later 
        building_blocks{idx} = block;
        idx = idx + 1;
	end
end

% Test that the building blocks are orthonormal
epsilon = 1 * 10.^(-10);
tests1_passed = 0; % should equal to 4032 when every block is orthonormal

for i = 1:size(building_blocks,1)
    for j = 1:size(building_blocks,1)
        if j ~= i % when not diagonal
            if sum(sum(building_blocks{i} .* building_blocks{j})) < epsilon
                tests1_passed = tests1_passed + 1;
            end
        end
    end
end

if tests1_passed == 64 * 64 - 64 % subtract diagonal elements
    disp("All blocks are orthonormal");
else
    disp("Some blocks are not orthonormal");
end

tests2_passed = 0; % should equal to 64 when 

% Test that the Frobenius inner product for each block
% with itself is close to 1 i.e each block is normalized
for i = 1:size(building_blocks, 1)
    s = sum(sum((building_blocks{1}) .* (building_blocks{1})));
    if s > s - epsilon && s < s + epsilon
        tests2_passed = tests2_passed + 1;
    end
end

if tests2_passed == 8 * 8
    disp("All block are normalized");
else
    disp("Some block are not normalized");
end

% Test that we can perform DCT on our test matrix and
% convert it back to the original form

% create random 8x8 matrix
rng(1,'twister'); % use fixed seed (becomes deterministic)
Br = randi([0,9], 8, 8);

Br_reconstructed = zeros(8,8);
for i = 1:size(building_blocks, 1)
    Br_reconstructed = Br_reconstructed + sum(sum(Br .* building_blocks{i})) .* building_blocks{i};
end

% We can just eyeball if these look the same
Br
Br_reconstructed

% Perform DCT-II on all of the source_blocks

% We will get 64 coefficient per one image block
IP = zeros(size(source_blocks,1),64);

% Do the forward DCT-II and store all coefficients in IP
for i = 1:size(source_blocks, 1)
    for j = 1:size(building_blocks, 1)
        IP(i, j) = sum(sum(source_blocks{i} .* building_blocks{j}));
    end
end

% Perform compression
tmp = sort(abs(IP(:)));
th = tmp(round(compression_factor * length(tmp)));

% Reconstruct the image from the coefficients and store them in
% result_blocks
for i = 1:size(source_blocks,1)
    reconst_block = zeros(8,8);
    for j = 1:size(building_blocks, 1)
        coeff = IP(i,j);
        if abs(coeff) < th
            coeff = 0;
        end
        reconst_block = reconst_block + coeff .* building_blocks{j};
    end
    result_blocks{i} = reconst_block;
end

% Stich the source image back together from result_blocks
result = zeros(1,size(B1,1) + 1);% add padding so concatenation doesnt fail
for i = 1:bc
    A = zeros(8,1); % also padding here for the same reason as ^
    for j = 1:bc
        lin_idx = (i - 1) * bc + j;
        A = cat(2, A, result_blocks{lin_idx});
    end
    result = cat(1, result, A);
end
result = result(2:size(result,1), 2:size(result,1)); % remove the padding

figure(3);
imshow(result);
