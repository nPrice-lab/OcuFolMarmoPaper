function allFig(repoPath)
%%  function for plotting all figures in the marmo behav paper
%   Input: repoPath -path for the OcuFolMarmoPaper repo. Assume current
%   directory is the repo if no input.
%   
%   Required package: 
%   "panel" at https://au.mathworks.com/matlabcentral/fileexchange/20003-panel
%   "colorcet" at https://colorcet.com/
%   "shaded error bar" at https://au.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
%   "violin plot" at https://github.com/bastibe/Violinplot-Matlab
%   
%   "shaded error bar" and "violin plot" are incorporated in the /src folder in this repo 
%
%   Note: need to use image editing software to combine legend/colorbar
%   into the figures for Fig 2-4.

    %% General parameters
    if nargin<1
        repoPath = cd;
    end
    
    subj = {'butch' 'brick' 'm1899'};
    subjS = {'Bu' 'Br' 'Ni'};
    pd = {'postsac' 'sftf' 'direction'};
    ext = 'condsCB';
    filePath = fullfile(repoPath,'data');
    figPath = fullfile(repoPath,'figures/');
    if ~exist(figPath, 'dir')
        mkdir(figPath)
    end
    
    %% Method (Fig 1)
    
    clear p
    close all
    p = panel();

    methodPlot(subj,subjS,pd{1},ext,filePath,figPath,p);    
    
    %% Experiment 1 results (Fig 2)
    
    clear p
    close all
    p = panel();
    expt1Plot(subj,subjS,pd{1},ext,filePath,figPath,p);
    
    % add legend
    p(3,2).select()
    lgd = legend({'10' '30' '50' '100' '200' '300'});
    lgd.Title.String = 'ms';
    lgd.Location = 'southeast';
    p.export(fullfile(figPath,'Expt1Legend.tif'),'-rp');      
    
    %% Experiment 2 results (Fig 3)
    
    clear p
    close all
    p = panel();
    subj2 = subj(1:2);
    subjS2 = subjS(1:2);
    expt2Plot(subj2,subjS2,pd{2},ext,filePath,figPath,p);
    
    for irow = 1:2
        for icol = 1:2
            p(irow,3,icol,5).select()
            colorbar
            yticks([])
            p(irow,3,icol,4).select()
            yticks([])
        end
    end

    p.export(fullfile(figPath,'Expt2colorbar.tif'),'-rp','-a1.2');

    p(2,2,1,3).select()
    p.de.margin = 10;
    lgd = legend({'0.04','0.08','0.16','0.31','0.62','1.24','2.48'});
    lgd.Title.String = {'cpd'};
    lgd.Location = 'westoutside';
        
    p(2,4,1,3).select()
    p.de.margin = 10;
    lgd = legend({'1.56','3.13','6.13','12.25','18.75','25'});
    lgd.Title.String = {'Hz'};
    lgd.Location = 'westoutside';
    
    p.export(fullfile(figPath,'Expt2Legend.tif'),'-rp','-a1.2');
    %% Experiment 3  results (Fig. 4)

    clear p
    close all
    p = panel();
    expt3Plot(subj,subjS,pd{3},ext,filePath,figPath,p);
    
    % add legend
    p(3,2,2).select()
    lgd = legend({'Same','Diff'});
%     lgd.Title.String = {'Direction','Congruency'};
    lgd.Location = 'southeast';
    p.export(fullfile(figPath,'Expt3Legend.tif'),'-rp');    
end