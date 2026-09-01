function ok = fetch_outcome_data(opts)
% FETCH_OUTCOME_DATA  Download (or verify) the four monthly outcome series.
%
% TARGET FILES (all in empirical/data/raw/, SDMX-CSV with TIME_PERIOD and
% OBS_VALUE columns)
% ----------------------------------------------------------------------
%   ip_ea.csv      euro-area industrial production, volume index,
%                  industry excl. construction (B-D), seasonally and
%                  calendar adjusted.        [Eurostat sts_inpr_m]
%   hicp_ea.csv    euro-area HICP all-items index (2015 = 100).
%                                             [Eurostat prc_hicp_midx]
%   rate1y_ea.csv  1-year nominal rate LEVEL, monthly average.  Baseline:
%                  12-month EURIBOR (available from 1999; caveat: contains
%                  a bank credit premium, material 2008-2012 -- documented
%                  in docs/CH7_DESIGN.md Sec. 2).  Preferred if you can
%                  locate it in the ECB portal: a 1-year OIS/EONIA-swap
%                  level series.               [ECB Data Portal, FM.*]
%   stoxx50_ea.csv EURO STOXX 50 price index, monthly close/average level.
%                                              [ECB Data Portal, FM.*]
%
% BEHAVIOUR
% ---------
% For each series: if the file exists it is parsed AND VALIDATED; if it
% validates it is kept, otherwise it is moved aside (to raw/rejected/) and
% re-downloaded.  Downloads try the candidate URLs in order (websave,
% falling back to urlwrite on older Octave).  Anything that downloads is
% validated before being accepted.  If a series cannot be obtained --
% offline machine, changed series keys, blocked network -- MANUAL
% INSTRUCTIONS are printed and ok is returned false.
%
% WHY THE VALIDATION EXISTS (read this before disabling it)
% ---------------------------------------------------------
% A file that merely PARSES is not usable data.  This package originally
% shipped four placeholder csvs -- correct filenames, correct columns,
% 1999-01..2025-12 with no gaps -- whose contents were generated, not
% observed: HICP rose by exactly 0.13 index points every month for 323
% consecutive months, the 1-year rate sat at exactly -0.5000 for years at
% a time, and neither industrial production nor the STOXX 50 showed the
% 2008 crisis or the 2020 collapse.  The old version of this function
% reported "all four outcome files ready" for them, and the whole chapter
% would have been estimated on fabricated numbers with nothing in any
% output saying so.  ea_check_series now refuses such files; its header
% documents the three fingerprints and the thresholds.  Set
% opts.skip_plausibility = true only if you have looked at a series and
% decided a fingerprint is a false alarm.
%
% INPUTS (opts, all optional)
% ---------------------------
%   opts.raw_dir           default <repo>/empirical/data/raw
%   opts.need_y0m0         first month the pipeline needs, default [2000 1]
%   opts.need_y1m1         last month,                     default [2019 12]
%                          (must match assemble_dataset's window)
%   opts.skip_plausibility default false; keeps coverage checks, drops the
%                          synthetic-data fingerprints
%   opts.force             default false; re-download even if the local
%                          file is present and valid
%
% OUTPUT
% ------
%   ok : true only when all four files are present, parse, and validate.
%
% IMPORTANT -- THE URLS ARE UNVERIFIED
% ------------------------------------
% The series keys below were written from portal documentation and could
% NOT be tested from the build environment (both ec.europa.eu and
% data-api.ecb.europa.eu are unreachable from it).  Portals do re-key
% series and Eurostat has been migrating STS datasets from the 2015 = 100
% base (I15) to 2021 = 100 (I21), so several candidates are tried per
% series.  If they all fail, follow the printed manual route, export as
% SDMX-CSV / csvdata, and save under the target filename; read_sdmx_csv
% only needs TIME_PERIOD and OBS_VALUE columns, whatever else is present.
% A different index base year is harmless here: the variables enter as
% 100*log(index) and a rescaling of the index is absorbed by the constant.

if nargin < 1, opts = struct(); end

% Self-sufficient on the path, whatever the working directory is.
if exist('ea_paths', 'file') ~= 2
    addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
end
P = ea_paths();

if ~isfield(opts, 'raw_dir'),   opts.raw_dir = P.raw; end
if ~isfield(opts, 'need_y0m0'), opts.need_y0m0 = [2000 1]; end
if ~isfield(opts, 'need_y1m1'), opts.need_y1m1 = [2019 12]; end
if ~isfield(opts, 'skip_plausibility'), opts.skip_plausibility = false; end
if ~isfield(opts, 'force'),     opts.force = false; end
if exist(opts.raw_dir, 'dir') ~= 7
    [okdir, msg] = mkdir(opts.raw_dir);
    assert(okdir == 1, 'fetch_outcome_data: cannot create %s (%s)', opts.raw_dir, msg);
end
need0 = 12 * opts.need_y0m0(1) + opts.need_y0m0(2);
need1 = 12 * opts.need_y1m1(1) + opts.need_y1m1(2);

eust = 'https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/';
ecb  = 'https://data-api.ecb.europa.eu/service/data/';

spec = struct('file', {}, 'urls', {}, 'manual', {});

spec(1).file = 'ip_ea.csv';
spec(1).urls = { ...
  [eust 'sts_inpr_m/M.PROD.B-D.SCA.I15.EA20/?format=SDMX-CSV&startPeriod=1999-01'], ...
  [eust 'sts_inpr_m/M.PROD.B-D.SCA.I15.EA19/?format=SDMX-CSV&startPeriod=1999-01'], ...
  [eust 'sts_inpr_m/M.PROD.B-D.SCA.I21.EA20/?format=SDMX-CSV&startPeriod=1999-01'], ...
  [eust 'sts_inpr_m/M.PROD.B-D.SCA.I21.EA19/?format=SDMX-CSV&startPeriod=1999-01']};
spec(1).manual = ['Eurostat Data Browser (ec.europa.eu/eurostat/databrowser), ' ...
  'dataset sts_inpr_m "Production in industry, monthly".  Filter: NACE ' ...
  'B-D (industry except construction), s_adj = SCA (seasonally and ' ...
  'calendar adjusted), unit = index (I15 if offered, else I21 -- the base ' ...
  'year does not matter, the series enters as 100*log), geo = EA20 (or ' ...
  'EA19), period from 1999-01.  Download > SDMX-CSV.'];

spec(2).file = 'hicp_ea.csv';
spec(2).urls = { ...
  [eust 'prc_hicp_midx/M.I15.CP00.EA/?format=SDMX-CSV&startPeriod=1999-01'], ...
  [eust 'prc_hicp_midx/M.I15.CP00.EA19/?format=SDMX-CSV&startPeriod=1999-01'], ...
  [eust 'prc_hicp_midx/M.I15.CP00.EA20/?format=SDMX-CSV&startPeriod=1999-01']};
spec(2).manual = ['Eurostat Data Browser, dataset prc_hicp_midx "HICP ' ...
  'monthly index".  Filter: coicop = CP00 (all-items), unit = I15 ' ...
  '(2015 = 100), geo = EA (euro area, changing composition), period from ' ...
  '1999-01.  Download > SDMX-CSV.  (This is the NSA index; see ' ...
  'CH7_DESIGN Sec. 2 for the seasonality discussion -- p = 12 absorbs it.)'];

spec(3).file = 'rate1y_ea.csv';
spec(3).urls = { ...
  [ecb 'FM/M.U2.EUR.RT.MM.EURIBOR1YD_.HSTA?format=csvdata&startPeriod=1999-01'], ...
  [ecb 'FM/M.U2.EUR.RT.MM.EURIBOR1YD.HSTA?format=csvdata&startPeriod=1999-01']};
spec(3).manual = ['ECB Data Portal (data.ecb.europa.eu), search "Euribor ' ...
  '1-year" and pick the MONTHLY series (frequency M), historical close, ' ...
  'average of observations through period, euro area.  Export > CSV.  ' ...
  'The series key is of the form FM.M.U2.EUR.RT.MM.EURIBOR1YD_.HSTA; if ' ...
  'the portal shows a different key, use whatever it gives you.  If you ' ...
  'find a 1-year OIS / EONIA swap LEVEL series with 1999+ coverage, ' ...
  'prefer it under the same filename and note the swap in the thesis ' ...
  '(it removes the bank credit premium discussed in CH7_DESIGN Sec. 2).'];

spec(4).file = 'stoxx50_ea.csv';
spec(4).urls = { ...
  [ecb 'FM/M.U2.EUR.DS.EI.DJES50I.HSTA?format=csvdata&startPeriod=1999-01']};
spec(4).manual = ['ECB Data Portal, search "Dow Jones Euro Stoxx 50" and ' ...
  'pick the MONTHLY price index (frequency M), historical close, average ' ...
  'of observations through period.  Export > CSV.  Series key of the form ' ...
  'FM.M.U2.EUR.DS.EI.DJES50I.HSTA.'];

ok = true;
for i = 1:numel(spec)
    fpath = fullfile(opts.raw_dir, spec(i).file);

    % --- keep an existing file only if it parses AND validates -----------
    if exist(fpath, 'file') == 2 && ~opts.force
        [good, chk] = try_validate(fpath, need0, need1, spec(i).file, opts);
        if good
            fprintf('fetch_outcome_data: %-14s OK, %d obs %s..%s (roughness %.2f)\n', ...
                    spec(i).file, chk.n_obs, ym_str(chk.ym0), ym_str(chk.ym1), ...
                    chk.roughness);
            continue
        end
        rej = fullfile(opts.raw_dir, 'rejected');
        if exist(rej, 'dir') ~= 7, mkdir(rej); end
        moved = fullfile(rej, spec(i).file);
        [mv, mvmsg] = movefile(fpath, moved, 'f');
        fprintf(2, '\nfetch_outcome_data: REJECTED the existing %s:\n', spec(i).file);
        for k = 1:numel(chk.msgs), fprintf(2, '    %s\n', chk.msgs{k}); end
        if mv == 1
            fprintf(2, '    moved to %s; will try to download a replacement.\n\n', moved);
        else
            fprintf(2, '    (could not move it aside: %s)\n\n', mvmsg);
        end
    end

    % --- download --------------------------------------------------------
    got = false;
    for u = 1:numel(spec(i).urls)
        try
            if exist('websave', 'file') || exist('websave', 'builtin')
                websave(fpath, spec(i).urls{u});
            else
                urlwrite(spec(i).urls{u}, fpath);  %#ok<URLWR>
            end
            [good, chk] = try_validate(fpath, need0, need1, spec(i).file, opts);
            if ~good
                fprintf(2, 'fetch_outcome_data: %s downloaded but rejected:\n', spec(i).file);
                for k = 1:numel(chk.msgs), fprintf(2, '    %s\n', chk.msgs{k}); end
                error('rejected');            % fall through to the next url
            end
            fprintf('fetch_outcome_data: downloaded %-14s %d obs %s..%s\n', ...
                    spec(i).file, chk.n_obs, ym_str(chk.ym0), ym_str(chk.ym1));
            got = true;
            break
        catch
            if exist(fpath, 'file') == 2, delete(fpath); end
        end
    end

    if ~got
        ok = false;
        fprintf(2, ['\nfetch_outcome_data: could NOT obtain %s automatically.\n' ...
                    'MANUAL ROUTE: %s\nSave the file as:\n    %s\n' ...
                    'Then rerun fetch_outcome_data.\n\n'], ...
                spec(i).file, spec(i).manual, fpath);
    end
end

if ok
    fprintf('fetch_outcome_data: all four outcome files ready and validated.\n');
else
    fprintf(2, ['fetch_outcome_data: NOT READY -- follow the manual instructions ' ...
                'above, then rerun.  assemble_dataset will refuse to run until ' ...
                'this reports all four ready.\n']);
end
end

% -------------------------------------------------------------------------
function [good, chk] = try_validate(fpath, need0, need1, name, opts)
% Parse and validate; never throws, so a bad file is a rejection, not a
% crash.
chk = struct('ok', false, 'msgs', {{sprintf('%s: could not be parsed.', name)}}, ...
             'n_obs', 0, 'ym0', NaN, 'ym1', NaN, 'roughness', NaN);
try
    [ym, val] = read_sdmx_csv(fpath);
catch err
    chk.msgs = {sprintf('%s: %s', name, err.message)};
    good = false;  return
end
chk = ea_check_series(name, ym, val, need0, need1);
if opts.skip_plausibility
    % Keep coverage/duplicate failures, drop the synthetic fingerprints.
    keep = {};
    for k = 1:numel(chk.msgs)
        if isempty(strfind(chk.msgs{k}, 'LOOKS SYNTHETIC'))  %#ok<STREMP>
            keep{end + 1} = chk.msgs{k};  %#ok<AGROW>
        end
    end
    chk.msgs = keep;
    chk.ok = isempty(keep);
end
good = chk.ok;
end

function s = ym_str(v)
if isnan(v), s = '?'; return; end
s = sprintf('%d-%02d', floor((v - 1) / 12), v - 12 * floor((v - 1) / 12));
end
