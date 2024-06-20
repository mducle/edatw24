%% Initialise spinw4
swroot = '/mnt/ceph-training/course_materials/spinw4';
addpath(genpath(fullfile(swroot, 'swfiles')));
addpath(genpath(fullfile(swroot, 'external')));
addpath(genpath(fullfile(swroot, 'dat_files')));
addpath('/mnt/ceph-training/course_materials/spinw_files');

% Set SpinW to use mex files
swpref().setpref('usemex', 1);

%%
% An example of using the form factor to determine the contribution of
% different ions to a spectrum from Sr3NiIrO6
%
% This spin chain material has both magnetic Ni2+ and Ir4+ ions. 
% The excitations associated with each are well separated in energy so
% we're going to use the form factor calculation and substitution to 
% find out which bands belong to which ion. 
% The data and model is taken from
% 
% Toth et al., Phys. Rev. B, 93 174422 (2016).
%
% https://link.aps.org/pdf/10.1103/PhysRevB.93.174422

sni = spinw();
sni.genlattice('lat_const', [9.61 9.61 11.1658], ...
               'angled', [90 90 120], 'spgr', 'R -3 c')
sni.addatom('r', [0 0 0], 'S', 0.5, 'label', 'MIr4');
sni.addatom('r', [0 0 0.25], 'S', 1, 'label', 'MNi2');
sni.gencoupling('maxDistance', 6)

% Values from PRB 93 174422 (2016)
Jxy = 21.6; Jz = 46.6; A = 4.95;
J2b = 0; J3a = -2.83; J3b = -1.37; Jtri = 1.46;
J3c = (Jtri - J3b + 2*J3a*0.5) / (0.5^2);

sni.addmatrix('label','J1','value',diag([Jxy Jxy Jz]),'color','white');
sni.addmatrix('label','J2b','value',J2b,'color','gray');
sni.addmatrix('label','J3a','value',J3a,'color','lightgray');
sni.addmatrix('label','J3b','value',J3b,'color','purple');
sni.addmatrix('label','J3c','value',J3c,'color','red');
sni.addmatrix('label','A','value',diag([0 0 A]));

sni.addcoupling('mat', 'J1', 'bond', 1);
sni.addcoupling('mat', 'J2b', 'bond', 3);
sni.addcoupling('mat', 'J3a', 'bond', 4);
sni.addcoupling('mat', 'J3b', 'bond', 5);
sni.addcoupling('mat', 'J3c', 'bond', 6);

sni.addaniso('A', 'MNi2');

sni.genmagstr('mode', 'direct', ... 
              'S', [0  0  0  0  0  0  0  0  0  0  0  0;
                    0  0  0  0  0  0  0  0  0  0  0  0;
                   -1 -1  1 -1  1 -1  1  1 -1  1 -1  1]);
plot(sni)

% Plot a spin wave spectrum and verify that there are bands at around 
% 30meV and 90meV as expected from the paper

% Examine the unit_cell field of the spinw object (type: sni.unit_cell) 
% – information on the sub-fields are here: 
% https://spinw.org/SWproperties/#unit_cell – make a copy of the form factor field. 
% Then set the form factor coefficients corresponding to Ir4+ to zero 
% (you can run: [~,coeff] = sw_mff('MIr4') to get an idea what the 
% coefficients of the form factor of Ir4+ is.
%
% With the Ir4+ form factor set to zero, calculate a powder spectrum 
% with 'formfact', true and check that the main spectral weight is now 
% from the 30meV mode. Integrate your powder spectrum over Q to get a 
% plot similar to the red shaded area in Fig 4d. of the PRB paper (
% reproduced above).
%
% Repeat the calculation with the Ni2+ form factor set to zero and see 
% if you can reproduce Fig 4d. (Use the saved form factor coefficients 
% you made earlier to restore the Ir4+ form factor, or just run the 
% setup script in the blue box above again).


ff0 = sni.unit_cell.ff
spec = sni.spinwave({[0 0 0] [1 1 1] [1 0 0] 200})
figure; sw_plotspec(spec)

powspec1 = sni.powspec(linspace(0.5, 4, 100), 'nRand', 100, ... 
                       'Evect', linspace(0,120,500), 'formfact', true)
figure; sw_plotspec(powspec1, 'dE', 10); caxis([0 1])
title('Powder spectrum with both Ni2+ and Ir4+'); 

%%

sni.unit_cell.ff(:,:,1) = 0;
powspec2 = sni.powspec(linspace(0.5, 4, 100), 'nRand', 100, ... 
                      'Evect', linspace(0,120,500), 'formfact', true)
figure; sw_plotspec(powspec2, 'dE', 10); caxis([0 1])
title('Powder spectrum with just Ni2+');

%%

sni.unit_cell.ff = ff0;
sni.unit_cell.ff(:,:,2) = 0;
powspec3 = sni.powspec(linspace(0.5, 4, 100), 'nRand', 100, ... 
                      'Evect', linspace(0,120,500), 'formfact', true)
figure; sw_plotspec(powspec3, 'dE', 10); caxis([0 1])
title('Powder spectrum with just Ir4+');
%%
figure; hold all;
ee = (powspec1.Evect(1:end-1) + powspec1.Evect(2:end)) / 2;
plot(ee, sum(powspec1.swConv, 2));
plot(ee, sum(powspec2.swConv, 2));
plot(ee, sum(powspec3.swConv, 2));
legend({'Full spectrum', 'Ni2+ only', 'Ir4+ only'})


%%
% An example of twinning... 

% Create the Pr(Ca,Sr)Mn2O7 model
prcasrmn2o7

% You can view the model file by uncommenting this:
%edit prcasrmn2o7

% Code to add twin and plot - comment / uncomment the code to see what the
% effect of the twin is.
% What happens if you change the twin axis or rotation?
% Why does the effect you see happen?

pcsmo.addtwin('axis', [0 0 1], 'phid', 90, 'vol', 0.5);

nQ = 101; Qhv = linspace(-2,2,nQ); Qkv = linspace(-2,2,nQ); Qlv = 0;
nE = 201; Ev = linspace(0,100,nE);
Ecut = [34 36];

[Qh, Qk, Ql] = ndgrid(Qhv,Qkv,Qlv);
spec = pcsmo.spinwave([Qh(:) Qk(:) Ql(:)]', 'sortMode', false, 'hermit', false);
spec = sw_egrid(spec,'component','Sperp','Evect',Ev);
spec = sw_instrument(spec,'dE',5);
spec3D = reshape(spec.swConv,nE-1,nQ,nQ);
Eidx = find(Ev>Ecut(1) & Ev<Ecut(2));
figure;
cut1 = squeeze(sum(spec3D(Eidx,:,:),1))/numel(Eidx)/(Ev(2)-Ev(1));
imagesc(Qhv/2,Qkv/2,cut1);
xlabel('Qh'); ylabel('Qk');
