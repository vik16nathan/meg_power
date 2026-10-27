function bst_headless_init()
%BST_HEADLESS_INIT: Gets a headless MATLAB session ready to run Brainstorm
%scripting calls (bst_process, bst_report, in_bst_data, ...), which none of
%them resolve to until Brainstorm itself has been started -- a fresh
%headless MATLAB session (e.g. launched by run_pipeline.sh) does not have
%brainstorm3 on the path or a protocol loaded by default, unlike an
%interactive session where the user already has the Brainstorm GUI open.
%
%Consolidates the fixes proven necessary in batch_cat12.m during headless
%debugging on 2026-07-27: path setup, starting Brainstorm in server mode,
%selecting the OmegaSubset protocol, disabling plugin auto-updates, and
%forcing off Brainstorm's online version-check (which hangs for a long time
%even with auto-updates disabled). Does NOT include batch_cat12.m's CAT12
%toolbox-symlink dedup: that fixes a duplicate registration that only gets
%recreated by bst_plugin('Install', 'cat12', ...), i.e. only relevant to
%scripts that actually call process_segment_cat12.

ProtocolName = 'OmegaSubset';

if ~exist('brainstorm', 'file')
    addpath('/export02/data/vikramn/brainstorm3');
end

% Server mode fully disables GUI popups, so this runs fine unattended with
% no display. If Brainstorm is already open (e.g. running this
% interactively from the GUI), this is a no-op and the already-loaded
% protocol/session is reused.
if ~brainstorm('status')
    brainstorm server
end

iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    error('Protocol "%s" not found in the current Brainstorm database.', ProtocolName);
end
gui_brainstorm('SetCurrentProtocol', iProtocol);

% Several scripts in this pipeline load plugins (Iso2mesh, DUNEuro) besides
% CAT12 -- disable auto-updates and the online version-check broadly so none
% of them can hit the same "reinstall + hang" behavior CAT12 did.
bst_set('AutoUpdates', 0);
global GlobalData;
GlobalData.Program.isInternet = 0;

end
