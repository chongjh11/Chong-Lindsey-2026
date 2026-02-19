%% Get the slip rate with stress-constrained - see Lindsey 2021 slip rate
% 
% 
%   Last modified on 24-Sep-2025
%   by Jeng Hann, Chong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sliprate_back,sliprate_ESPM,Iunlocked] = get_sliprate_stressconstrain(lockstart, lockend, driven,inputsliprate,KK, KK_sig0, faultcoords_fwd)

    % Get the depths points from locked fault
    vy = faultcoords_fwd(:,2);

    % Resize the stresses (from traction from Pylith)
    KK = squeeze(KK(2,1:end-1,1:end-1));
    KK_sig0 = squeeze(KK_sig0(2,1:size(KK,1)))'; % the sig0 due to ESPM model (influence on locked patch from slabbot etc)
    Vpl = inputsliprate*ones(length(KK),1); % locked patches slip rates

    sig0_backslip   = KK*Vpl; % synthetic stress for the backslip case 
    sig0_ESPM       = -KK_sig0.*Vpl; % stress calculated from Pylith for ESPM [optional, depends on the input slip in Pylith]

    % Set up the slip rates for the backslip and slab
    sliprate_back = zeros(size(sig0_backslip)); 
    sliprate_ESPM = zeros(size(sig0_ESPM)); 
    
    % Set the deep slip rate and compute the stress caused by it 
    sliprate_back(vy<=driven) = Vpl(vy<=driven); % impose creep deep slip below the locking depth
    sliprate_ESPM(vy<=driven) = Vpl(vy<=driven); % impose creep deep slip below the locking depth
    Iunlocked = (vy >= lockstart  |  (vy <= lockend & vy >= driven)); % define the unlocked patches (above lockstart and below lockend and above drive creep)

    %%% Calculate the stress again if you have changed the locking rates in
    %%% this script (not the mesh & slip impulses in Pylith) %%%
    deepcreepstressrate = -KK * sliprate_back; % calculated stressing rate from backslip (or locked portion only)
    KKsub = KK(Iunlocked,Iunlocked); % get the stress kernel from non-locked portion
    substressrate = sig0_backslip(Iunlocked) + deepcreepstressrate(Iunlocked); % calculate stress rate from non-locked portion
    substressrate_ESPM = sig0_ESPM(Iunlocked) + deepcreepstressrate(Iunlocked); % calculate stress rate from non-locked portion due to 
    
    % Use this subset of patches to balance the stress rate
    subsliprate = KKsub\substressrate;
    subsliprate_ESPM = KKsub\substressrate_ESPM;
    
    % Compute the final slip rate (give slip where it's unlocked)
    sliprate_back(Iunlocked) = subsliprate; % note: for 3D case add strike
    sliprate_ESPM(Iunlocked) = subsliprate_ESPM; 


end