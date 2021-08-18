% Task 1
%
% Elias Ferreira - up201502818
%
% Nuno Minhoto - up201604509
%
%

close all;
clear all;
load('ground_truth.mat');
 a=40;
 TruePositives=0;
 FalsePositives=0;
TotalRealFaces =0;
%scale=3;
for i = 1:20
     flag=0;
    fprintf(" img nº:%d \n",i);
   
imgOriginal = imread(sprintf('%d.png',i));      %%resize da imagem
[lin,col,~] = size(imgOriginal);
img = imresize(imgOriginal,[1024 1559]);

while(flag==0)
    img2=img;
    img2(:,:,1) = medfilt2(img(:,:,1), [3 3]);   %%filtro de média
    img2(:,:,2) = medfilt2(img(:,:,2), [3 3]);
    img2(:,:,3) = medfilt2(img(:,:,3), [3 3]);
    
Igray = rgb2gray(img2);                          %%aplicação do metodo de canny e operador de prewitt
Iedge = im2uint8(edge(Igray,'canny',0.2));
[Gmag,Gdir] = imgradient(Iedge,'prewitt');
%  figure,imshow(Gmag);title("bordas");
Iedge = repmat(Iedge,[1 1 1]);
Ifinal =  Iedge;  

    
ycbcr = rgb2ycbcr(img);
ycbcrdouble = im2double(ycbcr);


y = ycbcrdouble(:,:,1);                         %%obtenção das componentes ycbcr
cb = ycbcrdouble(:,:,2);
cr = ycbcrdouble(:,:,3);


mask = (cr >= 133/255) & (cr <= 173/255);
mask = mask & (cb >= 77/255) & (cb <= 127/255);


erstrel = strel('line',a,90);                   %%elemento estruturante em cruz
erstrel22 = strel('line',a,0);
disc = strel('disk', a, 0);                     %%elemento estruturante em disco
dilstrel = strel('disk', 20, 0);

mask = imerode(mask, erstrel);                  %%aplicação das erosões e dilatações
mask = imerode(mask, erstrel22);
mask = imdilate(mask, dilstrel);
st = stdfilt(y);
stmask = (st >= 1.82/255);
%figure,imshow(y);
%figure,imshow(st);
%figure,imshow(stmask);

m = mask & stmask ;
op = bwareaopen (m, 4000);
dilstrel2 = strel('disk', 4, 0);
cl = imclose(op, dilstrel2);
lel = bwareaopen (~cl, 40000);
lel = ~lel;
finalmask = lel;
% figure,
% %
% imshow(finalmask); title("final mask");


%applied_mask = regionprops(finalmask, 'ConvexImage').ConvexImage;

%applied_mask = finalmask & convex;

%figure, imshow(applied_mask);


 Ifinal = imresize(Ifinal ,[lin col]);
 
 finalmask = imresize(finalmask,[lin col]);   %resize de volta para o tamanho original
 
    

 finalmask2 = finalmask .* (~Ifinal);
  finalmask2 =imfill(finalmask2 ,'holes');
%  figure,
% imshow(finalmask2); title("perimetros");

finalmask2 = logical(finalmask2);

 %finalmask2 = imopen(finalmask2, se);

L4 = bwlabel(finalmask2,4);
%  figure,
% imshow(L4); title("L4");
stats = regionprops(L4,'centroid','area','perimeter','boundingbox','orientation','MajorAxisLength','MinorAxisLength');
%%obtensão das caracteristicas gerais da região

if isempty(stats)% teste para ver se nenhuma regiao foi detetada e reduzir erosao inicial 
    a=a-5;
    
 
else
 flag=1;
  break;
end    

end

%filtragem das regioes identificadas pelas caracteristicas

thblobarea = 0.05;
areasortedstats = sortStruct(stats, 'Area', -1); %sort blobs by area, descending order
maxarea = areasortedstats(1).Area; % blob with biggest area
% new stats containing only blobs with area greater than thblobarea * 100 % of the biggest
% blob
filteredblobs = stats([stats.Area] > maxarea * thblobarea);
nblobs = length(filteredblobs);
filteredimage = bwpropfilt(finalmask2, 'Area', nblobs);
% figure,
% imshow(filteredimage); title("filteredimage");
orientationThreshold = 40; % areas mais inclinadas que 40º são eliminadas
ratioThreshold = 1.8;
filteredimageorientationpositive = bwpropfilt(filteredimage, 'Orientation', [orientationThreshold 90]);
filteredimageorientationnegative = bwpropfilt(filteredimage, 'Orientation', [-90 -1 * orientationThreshold]);
filteredorientation = filteredimageorientationpositive + filteredimageorientationnegative;
orientationstruct = filteredblobs(([filteredblobs.Orientation] > orientationThreshold ) | ([filteredblobs.Orientation] < -1 * orientationThreshold));
blobratio = orientationstruct([orientationstruct.MajorAxisLength] ./ [orientationstruct.MinorAxisLength] < ratioThreshold);

% ncaras = length(orientationstruct);

    filteredorientation = logical(filteredorientation);
    
    delete_pos = [];

    for h = 1:length(orientationstruct)

        bounding_box = orientationstruct(h).BoundingBox;
        ratio = bounding_box(3) / bounding_box(4);

        if((bounding_box(3) / bounding_box(4) < 0.35) || (bounding_box(3) / bounding_box(4) > 1.3)) 

                ratio = bounding_box(4) / bounding_box(3);

                ratio_filter = bwpropfilt(filteredorientation, 'Area', [orientationstruct(h).Area-10 orientationstruct(h).Area+10]);

                filteredorientation = filteredorientation & ~ratio_filter;    
                


                delete_pos(end+1) = h;

        end

    end
    
    del = length(delete_pos);
    
    while del > 0
       
         orientationstruct(delete_pos(del))=[];
         
         del = del - 1;
        
    end

    ncaras = length(orientationstruct);
    
%     figure; imshow(filteredorientation);
    
    
figure,
imshow(imgOriginal );%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on;
%criação das bounding boxes
for k = 1:length(orientationstruct)
    bbox = orientationstruct(k).BoundingBox;
    rectangle('Position', bbox, 'EdgeColor', 'r', 'LineWidth', 3);
end
hold off

filteredlabel = bwlabel(filteredorientation, 4);


rgbfiltered = label2rgb(filteredlabel);

RGB4 = label2rgb(L4);
% [B,L] = bwboundaries(finalmask,'noholes');

%figure, imshow(label2rgb(L, @jet, [.5 .5 .5]))
%hold on
%for k = 1:length(B)
%   boundary = B{k};
%   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
%end

centroids = cat(1,orientationstruct.Centroid);
%figure,
%imshow(rgbfiltered)
%hold on
%plot(centroids(:,1),centroids(:,2),'b*')
%hold off


%figure,
%imshow(finalmask)////////////////////////////////////////////////////////////////

%comparação dos resultados com o ground_truth

bb = ground_truth_store(i).ground_truth;
n_faces = size(bb,1);

%figure(),imshow(img)
 %   hold on
    jaccard_index = zeros(ncaras,n_faces);%info n regioes, 
    %[lin,col] = size(finalmask);
    for k = 1 : ncaras
        BB = orientationstruct(k).BoundingBox;
   %     rectangle('Position', [BB(1),BB(2),BB(3),BB(4)],'EdgeColor','r','LineWidth',2);
        
        for j = 1 : n_faces
            A = zeros(lin,col);
            B = zeros(lin,col);
            A(ceil(BB(2)):floor(BB(2)+BB(4)),ceil(BB(1)):floor(BB(1)+BB(3))) = 1;
            B(bb(j,1):bb(j,2),bb(j,3):bb(j,4)) = 1;
            jaccard_index(j, k) = jaccard(A,B);
        end
    end
  %  hold off

    %saveas(gcf,sprintf('task1_results/%d.png',i));
    TP = length( find( ( jaccard_index(:,:) >= 0.5 ) ) );
    TruePositives = TruePositives + TP;
    FP = ncaras-TP;
    FalsePositives = FalsePositives + FP;

    %fprintf('TruePositve: %d --- FalsePositive: %d --- Missed: %d --- Percentage: %.2f %% ... ',TP, FP, n_faces-TP, TP/n_faces*100);

    TotalRealFaces = TotalRealFaces + n_faces;

    %figure; 
   % subplot(2, 1, 1);imshow(filteredimage, []);title('Resultado final th');
  % subplot(2, 1, 2);imshow(yBinary2, []);title('Y binarizado');
end

fprintf("Finished\n");
fprintf('TruePositve: %d --- FalsePositive: %d --- Missed: %d --- Percentage: %.2f %% ... ',TP, FP, n_faces-TP, TP/n_faces*100);
fprintf('Hit Percentage: %.2f %% --- FalsePositives: %d\n',TruePositives/TotalRealFaces*100, FalsePositives);