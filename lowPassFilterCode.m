function [filtData] = FftRft(trace,FPS,sigmaf)
%% trace = input signal; FPS = frames per second; sigmaf = width of your filter 

    X = fft(trace); 
    X2sided = abs(X); %this gives you the 2 sided FT 
    X1sided = X2sided(1:length(X2sided)/2); % the one sided FT 
    N = length(X1sided);
    
    %determine frequency limit we can resolve
    Lsec = N/FPS;
    f = FPS*((Lsec/2))/Lsec;

    % plot the frequency amplitudes
    freqAxis = linspace(0,f,N);
    figure; plot(freqAxis,X1sided); 
    title('Amplitudes as a function of frequency');
    
    %% 
    N2 = length(X2sided);
    figure;
    %create the mask/filter
    mask = zeros(1, N2);
    f = 0:N2/2;
    mask(1:N2/2+1) = exp(-(f/(2*sigmaf)).^2);
    mask(N2:-1:N2/2+2) = mask(2:N2/2);
    subplot(2,2,1); plot(freqAxis,mask(1:length(mask)/2)); title('Low-pass mask'); 
    % Multiply the FT by the mask
    Xfilt = X .* mask;
    subplot(2,2,2); plot(freqAxis,abs(Xfilt(1:length(Xfilt)/2))); title('Weighted transform');
    subplot(2,2,3); plot(trace); title('Original Signal');
    % Transform back
    filtData = ifft(Xfilt);
    subplot(2,2,4); plot(filtData); title('Low-pass filtered signal');
    
end 
