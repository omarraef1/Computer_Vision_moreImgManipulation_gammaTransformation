function [num_questions] = hw2(infile);

close all;
%part A
ourFile = fopen(infile);
format long g;
data = textscan(ourFile, '%f%f%f');
fclose(ourFile);
%data{:}
%whos('data')

%plot $
x = data{1};
y = data{2};
z = data{3};
whos('data')
figure, plot3(x,y,z,'.-');
sensorMat = [x y z];
%sensorMat

%random light spectra
rng(477);
randySpectra = rand(1600,101);
%randy
rgbCompute = randySpectra*sensorMat;
whos('rgbCompute')
whos('sensorMat')
whos('randySpectra')
%sensorMat
%randySpectra
%rgbCompute;
totalmax = max(rgbCompute(:));
scaleFactor = 255/totalmax;
scaledMat = rgbCompute./totalmax;
factoredMat = scaledMat.*255;
%factoredMat
whos('factoredMat')
whos('rgbCompute')
r=factoredMat(:,1);
rmax = max(r(:));
g=factoredMat(:,2);
gmax = max(g(:));
b=factoredMat(:,3);
bmax = max(b(:));

rmin = min(r(:));
gmin = min(g(:));
bmin = min(b(:));

rmax
gmax
bmax
rmin
gmin
bmin

whos('r')
whos('g')
whos('b')
%factoredMat
factoredMat(1,1)
factoredMat(1,2)
factoredMat(1,3)
%factoredMat

newImg = zeros(400,400,3, 'uint8');
lineCount = 1;
for vertical = 1:39
    for ro = 1:40
        for co = 1:3
            if vertical == 1
                if ro == 1
                    newImg(1:10,ro:ro*10,co) = factoredMat(ro,co);
                else
                    newImg(1:10,ro*10-10:ro*10+10,co) = factoredMat(ro,co);
                end
            else
                if ro == 1
                    newImg(vertical*10-10:vertical*10+10,ro:ro*10,co) = factoredMat(ro*vertical+40,co);
                else
                    newImg(vertical*10-10:vertical*10+10,ro*10-10:ro*10+10,co) = factoredMat(ro*vertical+40,co);
                end
            end
        end
    end
    lineCount = lineCount + 1;
end
whos('newImg')
figure, imshow(newImg);

% 2. least square method,
% estimate sensitivities
% plot real and estimated $
%factoredMat

%unsure
whos('sensorMat')
whos('factoredMat')
whos('rgbCompute')
psuedoInv = pinv(randySpectra);
leastSquare = psuedoInv*rgbCompute;
%leastSquare
whos('leastSquare')
%leastSquare


%THIS COMMENT STUB STATES THAT 
%THIS CODE IS THE PROPERTY OF OMAR R.G. (UofA Student)


figure, plot3(x,y,z, '.-');
hold on
plot3(leastSquare(:,1),leastSquare(:,2),leastSquare(:,3));
hold off
%whos('leastSquare')

% compute RMS error between each three estimated sensors and actual ones, 
xxx=sqrt(((leastSquare.^2)+(sensorMat).^2))./101;
xxx % report value

% compute RMS error between rgb calculated using actual sensors and
% estimated ones
rgbCompute1 = randySpectra*leastSquare;
xxx2=sqrt(((rgbCompute.^2)+(rgbCompute1).^2))./101;
%rgbCompute = randySpectra*sensorMat;
%---------_RMS ERRORS_end_---------------%



% 3. noise, min max, rms errors(2),
noise = 10*randn(1600,1);
%noise
noisedrgb = rgbCompute + noise;
noisedReal = psuedoInv*noisedrgb;
noisedReal
figure, plot3(x, y, z, '.-');
hold on
plot3(noisedReal(:,1), noisedReal(:,2), noisedReal(:,3));
hold off

%mins max
noisedReal(noisedReal>255) = 255;
noisedReal(noisedReal<0) = 0;
noisedReal
figure, plot3(x, y, z, '.-');
hold on
plot3(noisedReal(:,1), noisedReal(:,2), noisedReal(:,3));
hold off

xxx3=sqrt(((sensorMat.^2)+(noisedReal).^2))./101;
rgbCompute3 = randySpectra*noisedReal;
xxx4=sqrt(((rgbCompute.^2)+(rgbCompute3).^2))./101;
%rms for rgb xxx4
%---------_RMS ERRORS_end_---------------%

% 4. noise plotting?
%i = 0;
noise1 = 5*10;
noisedrgb1 = rgbCompute + noise1;
noisedReal1 = psuedoInv*noisedrgb1;
noisedReal1
figure, plot3(x, y, z, '.-');
hold on
plot3(noisedReal1(:,1), noisedReal1(:,2), noisedReal1(:,3));
hold off
noise2 = 10*10;
noisedrgb2 = rgbCompute + noise2;
noisedReal2 = psuedoInv*noisedrgb2;
noisedReal2
figure, plot3(x, y, z, '.-');
hold on
plot3(noisedReal2(:,1), noisedReal2(:,2), noisedReal2(:,3));
hold off



% part b. 5. gamma calibration

%whos('sensorMat')
%whos('noise')
%whos('leastSquare')
%whos('factoredMat')
%whos('randySpectra')
%whos('noisedReal')
%whos('rgbCompute')
%whos('xxx4')

gammaCorr = 255*(noisedReal/255).^(1/2.2);
gammaCorr

whos('randySpectra')

end

