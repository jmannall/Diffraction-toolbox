function [ir, tfmag, tfcomplex, tvec, fvec] = ProcessBTMResults(inFilePath, filehandlingparameters, controlparameters, cadFilePath, savePath)
    

    % Load BTM results
    load([inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_ir.mat'], 'irdirect', 'irgeom', 'irdiff');
    numReceivers = size(irdirect, 2);

    % Combine diffraction orders
    difforder = controlparameters.difforder;
    nir = length(irdirect);
    ndiff = 0;
    if difforder > 1
        load([inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_irhod.mat'], 'irhod');
        maxdifforder = length(irhod);
        ndiff = size(irhod{2}, 1);
        padding = zeros(ndiff - nir, numReceiers, 1);
        irdirect = [irdirect; padding];
        irgeom = [irgeom; padding];
        irdiff = [irdiff; padding];
%         irdiff = [irdiff; padding, zeros(ndiff, difforder - 1)];
        for i = 2:maxdifforder
            irdiff(:,:,i) = irhod{i};
        end
    end

    % Check if modelling source as a plane wave
    if isfield(controlparameters, 'Rstart')
        r = controlparameters.Rstart;
        irdirect = r * irdirect;
        irgeom = r * irgeom;
        irdiff = r * irdiff;
    end

    % Create frequency and time vectors
    nfft = controlparameters.nfft;
    fs = controlparameters.fs;
    fvec = fs/nfft*[0:nfft/2-1];
    tvec = 1/fs*[0:max(nir,ndiff)-1];

    % Create template
    template.complete = [];
    template.direct = [];
    template.geom = [];
    for i = 1:difforder
        idx = ['diff', num2str(i)];
        template.(idx) = [];
    end
    [ir, tfmag, tfcomplex] = deal(repmat(template, 1, 1));
    
    % Calculate frequency responses from the impulse responses
    [tfmag.complete, tfcomplex.complete, ir.complete] = IrToTf(irdirect + irgeom + sum(irdiff, 3), nfft);
    [tfmag.direct, tfcomplex.direct, ir.direct] = IrToTf(irdirect, nfft);
    [tfmag.geom, tfcomplex.geom, ir.geom] = IrToTf(irgeom, nfft);
    for i = 1:difforder
        idx = ['diff', num2str(i)];
        [tfmag.(idx), tfcomplex.(idx), ir.(idx)] = IrToTf(irdiff(:,:,i), nfft);
    end

    % Save results

    if controlparameters.saveFiles >= 1
        binary = fliplr(dec2bin(controlparameters.saveFiles));
        if binary(1) == '1'
            save(savePath, "ir", "tfmag", "tvec", "fvec", "tfcomplex");
            disp('Result saved')
        end
    end

    % Clear up and delete files
    delete(cadFilePath);
    path = [inFilePath,filesep,'results',filesep,filehandlingparameters.filestem];
    delete([path, '_eddata.mat']);
    delete([path, '_ir.mat']);
    delete([path, '_paths.mat']);
    delete([path, '_Rdata.mat']);
    delete([path, '_Sdata.mat']);
    if difforder > 1
        delete([path, '_irhod.mat']);
        delete([path, '_ed2data.mat']);
    end
end