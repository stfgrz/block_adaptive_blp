function ok = ea_fetch_v2_series(opts)
% EA_FETCH_V2_SERIES  Download (or verify) the additional monthly series of
% the v2 design from the ECB Data Portal.  Keys VERIFIED on 2026-09-13.
%
%   hicp_sa_ea.csv     ICP.M.U2.Y.000000.3.INX   HICP overall index, ECB,
%                      working-day and seasonally adjusted, 2015 = 100
%   loans_nfc_ea.csv   BSI.M.U2.Y.U.A20T.A.I.U2.2240.Z01.E   adjusted loans
%                      to NFCs, index of notional stocks, SA (from 2003-01)
%   ccb_nfc_ea.csv     MIR.M.U2.B.A2I.AM.R.A.2240.EUR.N   composite cost of
%                      borrowing, NFCs, percent (from 2003-01)
%   eonia_ea.csv       FM.M.U2.EUR.4F.MM.EONIA.HSTA   EONIA, monthly average
%   euribor1m_ea.csv   FM.M.U2.EUR.RT.MM.EURIBOR1MD_.HSTA   1M Euribor,
%                      monthly average (fallback indicator; credit premium)
%
% The 1-month OIS LEVEL (Refinitiv EUREON1M=) is not on the portal: see
% empirical/docs/DATA.md and import_ois_daily.m.
%
% opts.dry_run (default false) prints the URLs and does nothing; opts.only
% (cell of file names) restricts the set; opts.force re-downloads.  Each
% downloaded file is validated with ea_check_series over its own required
% window and gets a <file>.source.txt sidecar.  Returns true if every
% requested file is present and valid.

if nargin < 1, opts = struct(); end
if exist('ea_paths', 'file') ~= 2
    addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
end
P = ea_paths();
if ~isfield(opts, 'dry_run'), opts.dry_run = false; end
if ~isfield(opts, 'force'), opts.force = false; end
base = 'https://data-api.ecb.europa.eu/service/data/';
spec = { ...
  'hicp_sa_ea.csv',   'ICP/M.U2.Y.000000.3.INX',              [1999 1],  'HICP overall index, ECB, working-day and seasonally adjusted, 2015=100'; ...
  'loans_nfc_ea.csv', 'BSI/M.U2.Y.U.A20T.A.I.U2.2240.Z01.E',  [2003 1],  'Adjusted loans to euro area NFCs, index of notional stocks, SA'; ...
  'ccb_nfc_ea.csv',   'MIR/M.U2.B.A2I.AM.R.A.2240.EUR.N',     [2003 1],  'Composite cost of borrowing indicator, NFCs, percent p.a.'; ...
  'eonia_ea.csv',     'FM/M.U2.EUR.4F.MM.EONIA.HSTA',         [1999 1],  'EONIA, monthly average, percent p.a.'; ...
  'euribor1m_ea.csv', 'FM/M.U2.EUR.RT.MM.EURIBOR1MD_.HSTA',   [1999 1],  '1-month Euribor, monthly average, percent p.a.'};
if isfield(opts, 'only') && ~isempty(opts.only)
    spec = spec(ismember(spec(:, 1), opts.only), :);
end
ok = true;
for i = 1:size(spec, 1)
    fpath = fullfile(P.raw, spec{i, 1});
    url = sprintf('%s%s?format=csvdata&startPeriod=1999-01', base, spec{i, 2});
    need0 = 12 * spec{i, 3}(1) + spec{i, 3}(2);  need1 = 12 * 2019 + 12;
    if opts.dry_run
        fprintf('%-18s <- %s\n', spec{i, 1}, url);  continue
    end
    if exist(fpath, 'file') == 2 && ~opts.force
        try
            [ym, v] = read_sdmx_csv(fpath);  ea_check_series(spec{i, 1}, ym, v, need0, need1);
            fprintf('  %s present and valid\n', spec{i, 1});  continue
        catch
            fprintf('  %s present but invalid; re-downloading\n', spec{i, 1});
        end
    end
    try
        websave(fpath, url);
        [ym, v] = read_sdmx_csv(fpath);
        ea_check_series(spec{i, 1}, ym, v, need0, need1);
        ea_write_provenance(fpath, sprintf('ECB %s (%s)', strrep(spec{i, 2}, '/', '.'), spec{i, 4}), url, true);
        fprintf('  %s downloaded and validated (%d-%02d..%d-%02d)\n', spec{i, 1}, ...
                floor((ym(1) - 1) / 12), ym(1) - 12 * floor((ym(1) - 1) / 12), floor((ym(end) - 1) / 12), ym(end) - 12 * floor((ym(end) - 1) / 12));
    catch err
        ok = false;
        fprintf(2, '  %s FAILED: %s\n    manual route: open %s in a browser, save as %s\n', spec{i, 1}, err.message, url, fpath);
    end
end
end
