% The aim of this tutorial is to set up a SpinW model of a fairly complex
% "real-world" system - the half doped manganite Pr(Ca0.9Sr0.1)2Mn2O7,
% which has a lattice of checkerboard charge ordered Mn3+ and Mn4+ ions.
% Orbital ordering also means that the superexchange interactions have
% defined directions and so the symmetry-based shortcuts used by SpinW
% don't quite work, and more effort needs to be put in to define the model.
%
% The model is based on the paper: 
% Ground State in a Half-Doped Manganite Distinguished by Neutron Spectroscopy
% G.E. Johnstone, T.G. Perring, O. Sikora, D. Prabhakaran, and A.T. Boothroyd
% Phys. Rev. Lett. 109 237202 (2012), https://arxiv.org/abs/1210.7108

pcsmo = spinw;

% The actual spacegroup adopted by Pr(Ca0.9Sr0.1)2Mn2O7 is (a non-standard)
% 'Amam' (no. 63). If we want to put this in SpinW, we would need to use
% the standard setting, 'Cmcm', and to swap the axes [abc]->[cba].
% See PRB 80 064402 ( https://arxiv.org/abs/0907.1167 ) for the full
% crystal structure information. The most important is that:
lat = [5.408 5.4599 19.266];
alf = [90 90 90];
% And the Mn atom is in the 8g site with position:
mn_pos = [0.75 0.74812 0.09956];

% We can set this up to plot what the structure looks like
% We use the 'perm' option of genlattice to permute 'Cmcm' to 'Amam'
pcsmo.genlattice('lat_const', lat, 'angled', alf, 'spgr', 'C m c m', ...
    'perm', 'cba'); 
pcsmo.addatom('label', 'MMn4', 'r', mn_pos, 'S', 2, 'color', 'gold');
plot(pcsmo, 'range',[2 2 1])

% You should be able to see that the structure has two bilayers, one
% with layers at c~0.4 and c~0.6, and one at c~0.1 and c~0.9.

% Note that the "spgr" string can also be a list of the generators
% rather than just the group name. You can find generators for non-standard
% settings in the Bilbao Crystallographic Server: 
% http://www.cryst.ehu.es/cryst/get_gen.html (select generators not general
% positions). You need the generators in x,y,z form and note that you don't
% have to include the identity "x,y,z" generator.
% You could try to look up the generators for Amam (no. 63) and see that it
% agrees with the coordinate permutations above.

%%
% Now look at just the middle bilayer
plot(pcsmo, 'range', [0 2; 0 2; 0.3 0.7])

%%
% We can simplify this model by considering only a single bilayer instead
% the two in the original unit cell.
% However, a complication, is that the CE magnetic structure in this
% material is a 2-k structure so it can only be represented in SpinW
% using the supercell ('direct') method.
%
% In addition, the complex geometry of the exchange interactions means that
% we have two options for how to deal with defining them in SpinW:
%   1. Use the 'subidx' option in addcoupling() but this means inspecting
%      the structure to determine which bond is which.
%   2. Define the different atoms with different labels and use the 'atom'
%      option in addcoupling() to specify the bonds exactly.
%
% In the following we will use option 2. (The commented out section at the
% end of this mfile shows how it would be done using option 1).
%
% There is another problem with the nature of this model, however, in that
% not only is the structure 2-k, but the exchange interactions also do not
% obey the translational symmetry of the high temperature Amam unit cell.
% This means that to express them, we have to use a "structural" unit cell 
% which is the same size as the magnetic unit cell (since SpinW only allows
% you to define atoms within the first unit cell, and forces the bonds to 
% be translationally symmetric in the supercell). This means we have to use
% lattice parameters a and b which are twice as large as in the high 
% temperature Amam structure, which means that the hk0 indices outputted 
% by SpinW should be divided by 2 to get the equivalent positions in the 
% paper.

SM4 = 7/4;   % Spin length for Mn4+
SM3 = 7/4;   % Spin length for Mn3+

pcsmo = spinw;
pcsmo.genlattice('lat_const', lat.*[2 2 1], 'angled', alf, 'spgr', 'x,y+1/2,-z'); 

% In addition to defining the symmetry, because we are using custom labels
% (e.g. not of the form 'MXXN' where XX is the element abbreviation and N
% is the valence), we have to define the form factor and scattering length
[~,ffn3] = sw_mff('MMn3');
[~,ffn4] = sw_mff('MMn4');
myaddatom3 = @(x,y,z) pcsmo.addatom('label', x, 'r', y, 'S', SM3, 'color', z, ...
    'formfactn', ffn3, 'formfactx', 'MMn3', 'Z', 25, 'b', sw_nb('MMn3'));
myaddatom4 = @(x,y,z) pcsmo.addatom('label', x, 'r', y, 'S', SM4, 'color', z, ...
    'formfactn', ffn4, 'formfactx', 'MMn4', 'Z', 25, 'b', sw_nb('MMn4'));
myaddatom4('Mn4-up', [0 0 0.1], 'gold');
myaddatom4('Mn4-up', [0.5 0.5 0.1], 'gold');
myaddatom4('Mn4-dn', [0 0.5 0.1], 'gold');
myaddatom4('Mn4-dn', [0.5 0 0.1], 'gold');
myaddatom3('Mn3-up', [0.25 0.75 0.1], 'black');
myaddatom3('Mn3-up', [0.75 0.75 0.1], 'black');
myaddatom3('Mn3-dn', [0.25 0.25 0.1], 'black');
myaddatom3('Mn3-dn', [0.75 0.25 0.1], 'black');

% Plots the atoms
plot(pcsmo, 'range', [1 1 1])

%%
% Generate the magnetic structure
S0 = [0; 1; 0];

% Get the list of indices of which are spin up and which are down atoms:
spin_up = find(~cellfun(@isempty, strfind(pcsmo.table('matom').matom, 'up')));
spin_dn = find(~cellfun(@isempty, strfind(pcsmo.table('matom').matom, 'dn')));

SS = zeros(3, 16);
SS(:, spin_up) = repmat(S0, 1, numel(spin_up));
SS(:, spin_dn) = repmat(-S0, 1, numel(spin_dn));
pcsmo.genmagstr('mode', 'direct', 'S', SS)

% Plot only a single bilayer
plot(pcsmo, 'range', [0 1; 0 1; 0 0.2])

% Do you understand the above two sections of code?

% How would you define the atoms and magnetic structure if you didn't label
% the atoms as "up" or "down" when you define them?

%%
% Now define the exchange interactions for the Goodenough model.
% We have to force P1 symmetry because the nearest neighbour bonds (intra-
% and inter-chain) do not obey any translational symmetry of the structural
% orthorhombic cell.
pcsmo.gencoupling('forceNoSym', true)
% Print out table to determine which bond is which
pcsmo.table('bond', 1)

JF1 = -11.39;
JA = 1.5;
JF2 = -1.35;
JF3 = 1.5;
Jperp = 0.88;
D = 0.074;

pcsmo.addmatrix('label', 'JF1', 'value', JF1, 'color', 'green');
pcsmo.addmatrix('label', 'JA', 'value', JA, 'color', 'yellow');
pcsmo.addmatrix('label', 'JF2', 'value', JF2, 'color', 'white');
pcsmo.addmatrix('label', 'JF3', 'value', JF3, 'color', 'red');
pcsmo.addmatrix('label', 'Jperp', 'value', Jperp, 'color', 'blue');
pcsmo.addmatrix('label', 'D', 'value', diag([0 0 D]), 'color', 'white');

% The zig-zag chains couple Mn3-Mn4 with same spin.
pcsmo.addcoupling('mat', 'JF1', 'bond', 1, 'atom', {'Mn3-up', 'Mn4-up'})
pcsmo.addcoupling('mat', 'JF1', 'bond', 1, 'atom', {'Mn3-dn', 'Mn4-dn'})
% And opposite spins for the inter-chain interaction
pcsmo.addcoupling('mat', 'JA', 'bond', 1, 'atom', {'Mn3-up', 'Mn4-dn'})
pcsmo.addcoupling('mat', 'JA', 'bond', 1, 'atom', {'Mn3-dn', 'Mn4-up'})
% Second neighbour is the inter-layer coupling with our lattice parameters
pcsmo.addcoupling('mat', 'Jperp', 'bond', 2)
% JF3 couples Mn3 within the same zig-zag (same spin)
pcsmo.addcoupling('mat', 'JF3', 'bond', 3, 'atom', 'Mn3-up')
pcsmo.addcoupling('mat', 'JF3', 'bond', 3, 'atom', 'Mn3-dn')
% For the Mn4+ intra-zigzag next nearest neighbour we cannot use just the
% label name, because then we also get an inter-chain coupling (uncomment:
% pcsmo.gencoupling('forceNoSym', true)
% pcsmo.addcoupling('mat', 'JF2', 'bond', 8, 'atom', 'Mn4-up')
% plot(pcsmo, 'range', [0 2; 0 2; 0 0.2])
% And see what it produces).
% We actually only want interactions going in the +b direction which
% originates from the Mn4+ atom which have a=0.5.
% Find indexes of the Mn4+ atoms which have a=0.5:
idmid = find((~cellfun(@isempty, strfind(pcsmo.table('matom').matom, 'Mn4'))) ...
    .* (pcsmo.table('matom').pos(:,1)==0.5));
bond8 = pcsmo.table('bond', 8);
% Finds the bonds which start on one of these atoms and goes along +b
idstart = find(ismember(bond8.idx1, idmid) .* (bond8.dr(:,2)>0));
% Finds the bonds which ends on one of these atoms and goes along -b
idend = find(ismember(bond8.idx2, idmid) .* (bond8.dr(:,2)<0));
pcsmo.addcoupling('mat', 'JF2', 'bond', 8, 'subIdx', [idstart; idend]')

pcsmo.addaniso('D')
plot(pcsmo, 'range', [0 1; 0 1; 0 0.2])

%%
% Check that the structure we defined is optimum for the give exchanges
res = pcsmo.optmagsteep()
plot(pcsmo, 'range', [0 1; 0 1; 0 0.2])

% How many iterations did optmagsteep take?
% What does this mean?

%%
% Plots the spin wave along some directions

spec = pcsmo.spinwave({[0 0 0] [2 0 0] [2 0 2] [0 0 2] [0 0 0] 500}, 'hermit', false);
figure; sw_plotspec(spec);
specg = sw_egrid(spec, 'Evect', linspace(0,100,2000));
figure; sw_plotspec(specg, 'mode', 'color', 'dE',0.5);

% Does the dispersion agree with the paper?

%%
%{
% In this section we use the 'subIdx' option of addcoupling to specify
% the exchange interaction, rather through the atom labels.

SM4 = 7/4;   % Spin length for Mn4+
SM3 = 7/4;   % Spin length for Mn3+

pcsmo = spinw;
pcsmo.genlattice('lat_const', [[a b]*2 c], 'angled', alf, 'spgr', 'C 2/m'); 
pcsmo.addatom('label', 'MMn4', 'r', [0 0 0.1], 'S', SM4, 'color', 'gold');
pcsmo.addatom('label', 'MMn4', 'r', [0 0.5 0.1], 'S', SM4, 'color', 'gold');
pcsmo.addatom('label', 'MMn3', 'r', [0.25 0.25 0.1], 'S', SM3, 'color', 'black');
plot(pcsmo, 'range', [1 1 1])%, 'range', [0 2; 0 2; -0.2 0.2])

% Define the CE structure. There are 16 atoms in the magnetic unit cell.
% Atoms 1,2,5,6,9,10,15,16 are in one layer. Atoms 1,2,5,6 are Mn4+
S0 = [0; 1; 0];
pcsmo.genmagstr('mode', 'direct', 'nExt', [1 1 1], ...
    'S', [S0 S0 -S0 -S0  -S0 -S0 S0 S0  -S0 S0 S0 -S0  -S0 S0 S0 -S0]);
% Plot one of the bilayer only
plot(pcsmo, 'range', [0 1; 0 1; 0 0.2])

%%
% Now define the exchange interactions for the Goodenough model.
% We have to force P1 symmetry because the nearest neighbour bonds (intra-
% and inter-chain) do not obey any symmetry.
% (We only used C2/m above to save having to type in all the coordinates
% individually).
pcsmo.gencoupling('forceNoSym', true)
% Print out table to determine which bond is which
pcsmo.table('bond', 1)
% Plot the structure and click on each atom to figure out which is which
plot(pcsmo, 'range', [0 1; 0 1; 0 0.2]);   % Upper layer
plot(pcsmo, 'range', [0 1; 0 1; -0.2 0]);  % Lower layer

JF1 = -11.39;
JA = 1.5;
JF2 = -1.35;
JF3 = 1.5;
Jperp = 0.88;
D = 0.074;

pcsmo.addmatrix('label', 'JF1', 'value', JF1, 'color', 'green');
pcsmo.addmatrix('label', 'JA', 'value', JA, 'color', 'yellow');
pcsmo.addmatrix('label', 'JF2', 'value', JF2, 'color', 'white');
pcsmo.addmatrix('label', 'JF3', 'value', JF3, 'color', 'red');
pcsmo.addmatrix('label', 'Jperp', 'value', Jperp, 'color', 'yellow');
pcsmo.addmatrix('label', 'D', 'value', diag([0 0 D]), 'color', 'white');

%%
pcsmo.gencoupling('forceNoSym', true)
pcsmo.addcoupling('mat', 'JF1', 'bond', 1, 'subIdx', [1 2 8 17 18 19 29 32 ...
                                                      5 6 11 21 22 24 27 28])
pcsmo.addcoupling('mat', 'JA', 'bond', 1, 'subIdx', [3 7 9 10 15 16 30 31 ...
                                                     4 12 13 14 20 23 25 26])
plot(pcsmo, 'range', [0 1; 0 1; -0.2 0])
pcsmo.addcoupling('mat', 'Jperp', 'bond', 2)
pcsmo.addcoupling('mat', 'JF2', 'bond', 8, 'subIdx', [2 3 10 11  6 7 14 15])
pcsmo.addcoupling('mat', 'JF3', 'bond', 3, 'atom', 'MMn3')
pcsmo.addaniso('D')
plot(pcsmo, 'range', [0 1; 0 1; -0.2 0])
%}