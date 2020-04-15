%% Audio input

clear; clc;

%Exemplo de audio com duracao de 60s e 3 locutores
audio_file = "plenario-01_01-10-14-58_1031.MP3";

[audioIn, fs] = audioread(audio_file);
time = (1:length(audioIn))/fs;
%% Feature Extraction

win = round(0.02*fs); % fs = 44100 amostras/s. Janela de 200ms, o que equivale a 882 amostras
step = round(0.01*fs); % Step de 100ms (441 amostras)


%Fun��es de extra��o de features retiradas do livro AUDIO ANALYSIS: A MATLAB Approach (Giannakopoulos T e Pikrakis A)
stEnergy = ShortTimeEnergy(audioIn, win, step); %Energia
stZcr = ShortTimeZCR(audioIn, win, step); %ZCR

% figure
% subplot(2,1,1);
% plot(stEnergy);
% title("Energy");
% 
% subplot(2,1,2);
% plot(stZcr);
% title("ZCR");
%% Voiced Segments

% M�todo de detec��o de sil�ncio baseado em Rabiner (1975).
% Apesar de muito simples, funciona bem em �udios com boa SNR

% speechDetection(): funcao retirada - por�m modificada - de Theodoros Giannakopoulos (2020)
% Silence removal in speech signals (https://www.mathworks.com/matlabcentral/fileexchange/28826-silence-removal-in-speech-signals), 
% MATLAB Central File Exchange. Retrieved April 15, 2020.



[Limits] = speechDetection(audioIn, stEnergy, stZcr, win, step, fs); % [Limits]: matriz em que cada linha indica as amostras
                                                                     % de inicio e fim do segmento de voz
                                                                     % Exemplo: 3o segmento encontra-se entre as amostras 337366 e 394695
                                                                    
                                                                     
audioOutSamples = sum(Limits(:,2)-Limits(:,1)); 
audioOut = []; % Arquivo de audio com todos os per�odos de sil�ncio retirados


%Evitar gerar gr�ficos com arquivos de �udio maiores (10 min, por exemplo).
figure
hax=axes; 
plot(audioIn)
hold on 

% audioOut � iterativamente criado
for i=1:length(Limits)
    range = Limits(i,1):Limits(i,2); range = transpose(range); % amostras que correspondem ao periodo de fala
    voicedSegment = audioIn(range);                            % audio de entrada correspondente �s amostras de fala
    audioOut = [audioOut; voicedSegment, range];               % iterativamente audioOut vai sendo criado. 
                                                               % A primeira coluna corresponde �s amplitures das amostras de audio
                                                               % e a segunda coluna, ao �ndice da amostra no �udio original

    
    line([Limits(i,1) Limits(i,1)],get(hax,'YLim'),'Color',[1 0 0])
    line([Limits(i,2) Limits(i,2)],get(hax,'YLim'),'Color',[1 0 0])
end
%% MFCC Feature Extraction
% 
% 
% Features do �udio que ser�o efetivamente utilizadas no processo de segmenta��o 
% e clusteriza��o.

%Fun��o j� dispon�vel no matlab. Por padr�o, essa fun�ao usa uma janela de
%30ms com passo de 10ms, mas � poss�vel  modificar.
coeffs = mfcc(audioOut(:,1),fs,"LogEnergy","Ignore"); %corresponde a uma matriz com 13 colunas, cada uma representando uma feature
                                                      %o numero de linhas � o n�mero de amostras dispon�veis no dom�nio
                                                    


%% Segmentation/Clusterization
% Abordagem inicial baseada no artigo de Imseng D e Friedland G 
%(Robust Speaker Diarization for Short Speech Recordings)
%
% Proximo passo, se necess�rio: abordagem do The ICSI RT-09 Speaker Diarization System


% Initial segmentation generated by uniformly partitioning the audio into k 
% segments of the same length.
seg_window = 200; % 200 amostras de mfcc (gera 21 segmentos)
L = length(coeffs);
n_segs = floor(L/seg_window);
initial_segments = cell(n_segs, 1);

curPos = 1;

options = statset('Display','final');
g = 5;

for i=1:n_segs
    frame = coeffs(curPos:curPos+seg_window,:); %Iteracao em cada frame
    frame_gmm = fitgmdist(frame,g,'CovarianceType','diagonal','Options',options); % Cada cluster � modelado com uma GMM contendo g gaussianas
    initial_segments{i,1} = frame_gmm; %segmentos iniciais
    
    curPos = curPos + seg_window;
end

% Pr�ximas etapas:

% 1) Re-Segmentation: Run Viterbi alignment to find the opti-mal path of frames and models. 
%    The classifications based on 10 ms frames are very noisy; a minimum duration of 2:5 seconds is assumed 
%    for each speech segment.

% 2) Re-Training: Given the new segmentation of the audio track, compute new GMMs for each of the clusters.

% 3) Cluster Merging: Given the new GMMs, try to find the two clusters that most likely represent the same speaker. 
%    This is done by computing a score based on the Bayesian Information Criterion (BIC) of each of the clusters 
%    and the BIC score of a new GMM trained on the merged segments for two clusters. 
%    If the BIC score of the merged GMM is larger than or equal to the sum of the individual BIC scores, 
%    the two models are merged and the algorithm continues at the re-segmentation step using the merged GMM. 
%    If no pair is found, the algorithm stops