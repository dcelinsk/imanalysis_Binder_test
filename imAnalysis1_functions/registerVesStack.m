function [J, registeredTransformations] = registerVesStack(regStack,regTemp)

% Register Stacks

%red = importedImages(:,:,1:550);
%green = importedImages(:,:,551:1100);

%redMean = meanproj(red(:,:,250:300));

% Adapted Registration Code from Chris' GUI
%regStack = red;
%regTemp = redMean;
subpixelFactor=50;
registeredTransformations = zeros(4,size(regStack,3));

parfor n=1:size(regStack,3)
            imReg = regStack(:,:,n);
            [out1,~] = dftregistration(fft2(regTemp),fft2(imReg),subpixelFactor);
            registeredTransformations(:,n) = out1;
            %assignin('base',['registeredTransformations(:,' num2str(n) ')'],out1);
            %regStackString(:,:,n) = uint16(abs(ifft2(out2)));
            %assignin('base',[regStackString '(:,:,' num2str(n) ')'],uint16(abs(ifft2(out2))));

end
        %clear regTempC
%Conventional transform
parfor ind = 1:size(regStack,3);
dxdy = [registeredTransformations(4,:)',registeredTransformations(3,:)'];
 J(:,:,ind) = imtranslate(regStack(:,:,ind),dxdy(ind,:),'nearest','FillValues',0,'OutputView','same');
end

end