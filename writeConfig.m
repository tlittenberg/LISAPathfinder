function writeConfig(filename,delimiter)
% DESCRIPTION: Writes a configuration file to text for the MCMC
% micrometeorite code.
%
% USAGE: writeConfig(filename)
%        writeConfig(filename,delimter)
%
% Ira Thorpe
% 2016.03.02

% parse inputs
if nargin < 2
  delimiter = '%';
end
if nargin < 1
  filename = 'config.dat';
end

% gather data

groupNames = {...
  'Spacecraft Geometry',...
  'Housing 1 Geometry',...
  'Housing 2 Geometry',...
  'Spacecraft Mass Properties'};

geoParams = getGeometry();
groupParams{1} = [...
  geoParams,...
  LPFParam.EOM_RB_X,...
  LPFParam.EOM_RB_Y,...
  LPFParam.EOM_RB_Z];

groupParams{2} = [...
  LPFParam.EOM_H1SC_X,...
  LPFParam.EOM_H1SC_Y,...
  LPFParam.EOM_H1SC_Z];

groupParams{3} = [...
  LPFParam.EOM_H2SC_X,...
  LPFParam.EOM_H2SC_Y,...
  LPFParam.EOM_H2SC_Z];

[IH1,IH2] = getHousingInertias();

groupParams{4} = [...
  LPFParam.EOM_SC_M,...
  IH1,...
  IH2];
  



%

% open file & write header
fid = fopen(filename,'w');
fprintf(fid,['%s Configuration file for MCMC LPF Micrometeorite code\n'....
  '%s Date written %s'], delimiter,delimiter,datestr(now,'yyyy-mm-dd HH:MM:SS'));



% write
for ii = 1:numel(groupNames)
  % write group header
  fprintf(fid,'\n%s %s',delimiter,groupNames{ii});
  % write group contents
  for jj = 1:numel(groupParams{ii})
    p = groupParams{ii}(jj);
    ustr = char(p.units);
    ustr = strrep(ustr,'[','');
    ustr = strrep(ustr,']','');
    fprintf(fid,'\n%s, %8.6e, %s, %s',p.name,p.double,ustr,p.description);
  end
end

% close
fclose(fid);


end

function geomParams = getGeometry()
% load geometry file
data = load('corners.mat');

% height
geomParams(1) = LPFParam(...
  'SC_H',data.h,'m','Height of spacecraft');

% corners
kk = 1;
for ii = 1:8
  kk = kk+1;
  geomParams(kk) = LPFParam(...
    sprintf('SC_BOT_CORNER_%i_X',ii),data.pt(ii,1),'m',...
    sprintf('x coordinate of spacecraft bottom deck corner %i',ii));
    kk = kk+1;
  geomParams(kk) = LPFParam(...
    sprintf('SC_BOT_CORNER_%i_Y',ii),data.pt(ii,2),'m',...
    sprintf('y coordinate of spacecraft bottom deck corner %i',ii));
end

end

function [IH1params, IH2params] = getHousingInertias()
% function to get the moments of inertia in the housing frames

M = double(LPFParam.EOM_SC_M);

IM = reshape(double([...
  [LPFParam.EOM_SC_IXX; LPFParam.EOM_SC_IXY; LPFParam.EOM_SC_IXZ]...
  [LPFParam.EOM_SC_IYX; LPFParam.EOM_SC_IYY; LPFParam.EOM_SC_IYZ]...
  [LPFParam.EOM_SC_IZX; LPFParam.EOM_SC_IZY; LPFParam.EOM_SC_IZZ]...
  ]'),3,3);

RBM = double([LPFParam.EOM_RB_X,LPFParam.EOM_RB_Y,LPFParam.EOM_RB_Z]);

RH1M = double([LPFParam.EOM_H1SC_X,LPFParam.EOM_H1SC_Y,LPFParam.EOM_H1SC_Z]);
RH2M = double([LPFParam.EOM_H2SC_X,LPFParam.EOM_H2SC_Y,LPFParam.EOM_H2SC_Z]);

% parallel axis theorem to body frame
IB = IM - M*(RBM*RBM'*eye(3)-RBM'*RBM);

% parallel axis theorem to housing frames
RH1B = RH1M-RBM;
RH2B = RH2M-RBM;

IH1 = IB + M*(RH1B*RH1B'*eye(3)-RH1B'*RH1B);
IH2 = IB + M*(RH2B*RH2B'*eye(3)-RH2B'*RH2B);

% build LPFParam arrays
coordNames = {'X','Y','Z'};
kk = 0;
for ii = 1:3
  for jj = 1:3
    kk = kk+1;
    % H1
    IH1params(kk) = LPFParam(...
      sprintf('EOM_SC_IH1_%s%s',coordNames{ii},coordNames{jj}),... %name
      IH1(ii,jj), ...% value
      unit('kg*m^2'), ...% units
      sprintf('%s%s component of spacecraft moment of inertia tensor about H1',...
      coordNames{ii},coordNames{jj}));
    %H2
    IH2params(kk) = LPFParam(...
      sprintf('EOM_SC_IH2_%s%s',coordNames{ii},coordNames{jj}),... %name
      IH2(ii,jj), ...% value
      unit('kg*m^2'), ...% units
      sprintf('%s%s component of spacecraft moment of inertia tensor about H2',...
      coordNames{ii},coordNames{jj}));
  end
end


end