%% Recortar as caras de todas as imagens
% [3º, 1º, 4º-3º , 2º-1º ]

%Variaveis inicias

close all
clear all
clc
GT = load('ground_truth.mat');

TP = 0;
FP = 0;
TSigns = 0;
Badmask_count_1 = 0;
Badmask = 0;
Gmask_count_1 = 0;
Gmask_count = 0;
Gmask = 0;
Nmask_count = 0;
Nmask_count_T = 0;
Nmask = 0;
Badmask_count_FP = 0;
Badmask_count_M = 0;
Badmask_count = 0;
CMatrix=zeros(1,5);


for i = 1:30
    
    Img_N = GT.ground_truth_store(i).file;
    Img = imread (Img_N);
    
    CaraAnaliz = GT.ground_truth_store(i).ground_truth;
    
    for j = 1:size(CaraAnaliz,1)

        x = CaraAnaliz(j,3);
        y = CaraAnaliz(j,1);
        w = CaraAnaliz(j,4)-CaraAnaliz(j,3);
        l = CaraAnaliz(j,2)-CaraAnaliz(j,1);
        
        Img_fin = imresize(imcrop(Img,[x,y,w,l]),[400 400]);
        
        Resposta = GT.ground_truth_store(i).mask(j);
        
        YCbCr = rgb2ycbcr(Img_fin);

        Y=YCbCr( :,:,1);
        Cb=YCbCr( :,:,2);
        Cr=YCbCr( :,:,3);

        Ppele = (77<Cb) & (173>Cr) & (Cb<127) & (Cr>133); %imagem que indica o pixeis q sao pele
 
        %Faz uma possivel correçao de buracos no meio da pele
        tam=2;
        se=strel('disk',tam);
        imgFilled=imfill(Ppele,'holes');
        imgOpen = imopen(imgFilled, se);
        
      
        imgFiltered = medfilt2(imgOpen);

        Imj_Avl = imgFiltered;

        Calc=Imj_Avl((floor(size(Imj_Avl,1)/2):size(Imj_Avl,1)),:); 
        PeleNUM=sum(Calc(:)==1);
        MaskNUM=sum(Calc(:)==0);

        if PeleNUM <= MaskNUM
                    Ret=1; %passa para descobrir se esta bem posta ou nao
                    teste = 'mascara';

        else
                    Ret=0;  % e detetado sem mascara
                    teste = 'sem mascara';
                    if(strcmp(Resposta,'without_mask'))
                        Nmask_count_T=Nmask_count_T + 1;
                        TP=TP+1;
                    elseif(strcmp(Resposta,'mask_weared_incorrect')) 
                        Badmask_count_FP = Badmask_count_FP + 1;
                    elseif(strcmp(Resposta,'with_mask'))
                        Badmask_count_M = Badmask_count_M + 1; 
                    end
        end
        if Ret == 1
                LinhaMas=(edge(Cr./Cb,'sobel','horizontal'))+(edge(histeq(Y),'sobel','horizontal')); %Usamos sobel horizontal para detetar a linnha da mascara, em Cr/Cb a mascara esta bem delimitada, entao fazemos um OU com a imaguem Y para detetarmos melhor a mascara
                LinhaMas=bwmorph(LinhaMas,'bridge'); %Usamos bridge para ligar pontos proximos
                LinhaMas=bwmorph(LinhaMas,'skel');%Usamos skel para tornar a linha da mascara o mais fina possivel depois de todas as opracoes
                Aux = bwconncomp(LinhaMas,8);
                Meds = regionprops(Aux, 'MajorAxisLength');
                CaixaMed = regionprops(Aux, 'BoundingBox');
                
                
                max_h=0;
                max_i=0;
                for k = 1 : length(Meds)
                        Caixa_fin = Meds(k).MajorAxisLength;
                        if Caixa_fin > max_h
                                max_h=Caixa_fin;
                                max_i=k;
                        end
                end
                Caixa_fin = CaixaMed(max_i).BoundingBox;
                
                

                %calculo do ratio, sendo M2T a distancia da linha da mascara ate a
                %testa e Caixa_fin(3) e a largura da boundingbox
                %detetada pelo traco da mask
                linhafinal=round(Caixa_fin(2));
                colunameio=round((Caixa_fin(1)+Caixa_fin(3))/2);
                M2T=sum(Imj_Avl((1:linhafinal),colunameio));
                Rat = M2T/Caixa_fin(3);
                %determinacao do ratio atravez de testes
                if(strcmp(Resposta,'mask_weared_incorrect'))
                    Badmask_count_1 = Badmask_count_1 + 1;
                    Badmask(Badmask_count_1)=Rat;
                elseif(strcmp(Resposta,'without_mask'))
                    Nmask_count = Nmask_count + 1;
                    FP = FP +1;
                    Nmask(Nmask_count) = Rat;
                elseif(strcmp(Resposta,'with_mask'))
                    Gmask_count_1 = Gmask_count_1 + 1;
                    Gmask(Gmask_count_1) = Rat;
                end



                %Diferenca entre mal posta e bem
                if Rat<0.1
                    Ret=2;
                    if(strcmp(Resposta,'with_mask'))
                        TP=TP+1;
                        Gmask_count = Gmask_count + 1;
                        teste = 'bem posta';
                    end
                else
                    Ret=3;
                    if(strcmp(Resposta,'mask_weared_incorrect'))
                        TP=TP+1;
                        Badmask_count = Badmask_count + 1;
                        teste = 'mal posta';
                    end
                end
        end 
    end    
end


%Bem colocada/Sem máscara/Mal Colocada/VerdadeirosPositivos/FalsosPositivos

CMatrix(1,1) = Gmask_count_1;
CMatrix(1,2) = Nmask_count_T;
CMatrix(1,3) = Nmask_count;
CMatrix(1,4) = TP;
CMatrix(1,5) = FP;

CMatrix
